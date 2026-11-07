# Resumable night runner: for each candidate article, list its spreadsheets in the PMC open data bucket,
# download them, screen them with detect.py, store the result in SQLite and delete the files
# One committed SQLite transaction per article is the checkpoint, so it can be killed and restarted at any time
#
# Usage: run_batch.py [--config config.ini] [--ignore-window | --wait-for-window] [--workers N] [--max-mbps X] [--limit N]
#                     [--candidates FILE] [--db FILE] [--stop-after-gene-papers N --species mouse,rat]
import argparse
import configparser
import datetime
import json
import os
import re
import shutil
import sqlite3
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lowprio

BASE_DIR = Path(__file__).resolve().parent.parent
WORK_DIR = BASE_DIR / "work"
PAUSE_FILE = BASE_DIR / "PAUSE"
NIGHTS_LOG = BASE_DIR / "logs" / "nights.log"
S3_URL = "https://pmc-oa-opendata.s3.amazonaws.com/"
MAX_ATTEMPTS = 2

SCHEMA = """
CREATE TABLE IF NOT EXISTS articles (pmcid TEXT PRIMARY KEY, status TEXT, n_files INTEGER, bytes INTEGER,
	n_gene_columns INTEGER, n_flags INTEGER, species TEXT, tiers TEXT, error TEXT, seconds REAL, attempts INTEGER,
	finished_utc TEXT);
CREATE TABLE IF NOT EXISTS files (pmcid TEXT, file TEXT, bytes INTEGER, status TEXT);
CREATE TABLE IF NOT EXISTS gene_columns (pmcid TEXT, file TEXT, sheet TEXT, col TEXT, species TEXT, top20_human INTEGER,
	top20_mouse INTEGER, top20_rat INTEGER, non_empty INTEGER, symbol_cells INTEGER, flagged_cells INTEGER,
	new_name_cells INTEGER, old_name_text_cells INTEGER);
CREATE TABLE IF NOT EXISTS flags (pmcid TEXT, file TEXT, sheet TEXT, cell TEXT, species TEXT, kind TEXT, rule TEXT,
	value TEXT, guess_original TEXT, current_symbol TEXT, row_ensembl TEXT, row_ensembl_symbols TEXT, tier TEXT,
	left_value TEXT, right_value TEXT);
CREATE TABLE IF NOT EXISTS sessions (started_utc TEXT, ended_utc TEXT, reason TEXT, articles INTEGER, errors INTEGER,
	bytes INTEGER);
"""


def ParseTime(text):
	hours, minutes = text.strip().split(":")
	return datetime.time(int(hours), int(minutes))


def InWindow(now, start, end):
	# The window can cross midnight (23:30 -> 07:30)
	current = now.time()
	if start <= end:
		return start <= current < end
	return current >= start or current < end


# ---------------- worker side ----------------

SETTINGS = None
SESSION = None


def InitWorker(settings):
	global SETTINGS, SESSION
	import requests
	import detect
	lowprio.SetLowPriority()
	SETTINGS = settings
	SESSION = requests.Session()
	SESSION.headers["User-Agent"] = "PaperRerun-pilot/0.1 (volunteer research audit)"
	detect.InitWorker()


def ListSpreadsheets(pmcid):
	for attempt in range(3):
		try:
			response = SESSION.get(S3_URL, params={"list-type": "2", "prefix": pmcid + "."}, timeout=60)
			response.raise_for_status()
			break
		except Exception:
			if attempt == 2:
				raise
			time.sleep(5 * (attempt + 1))
	keys = re.findall(r"<Key>([^<]+)</Key><LastModified>[^<]*</LastModified>.*?<Size>(\d+)</Size>", response.text)
	if not keys:
		return None
	latestVersion = max(int(key.split("/")[0].split(".")[1]) for key, size in keys)
	return [(key, int(size)) for key, size in keys
		if int(key.split("/")[0].split(".")[1]) == latestVersion and re.search(r"\.xlsx?$", key, re.I)]


def DownloadCapped(key, targetPath, bytesPerSecond):
	# Stream in 1 MB chunks and sleep when ahead of the cap (cap is per worker)
	targetPath.parent.mkdir(parents=True, exist_ok=True)
	startTime = time.monotonic()
	written = 0
	with SESSION.get(S3_URL + key, timeout=300, stream=True) as response:
		response.raise_for_status()
		with open(targetPath, "wb") as handle:
			for chunk in response.iter_content(chunk_size=1024 * 1024):
				handle.write(chunk)
				written += len(chunk)
				if bytesPerSecond > 0:
					aheadSeconds = written / bytesPerSecond - (time.monotonic() - startTime)
					if aheadSeconds > 0:
						time.sleep(aheadSeconds)
	return written


def ProcessArticle(pmcid):
	import detect
	startTime = time.monotonic()
	articleDir = WORK_DIR / pmcid
	result = {"pmcid": pmcid, "status": "", "n_files": 0, "bytes": 0, "files": [], "gene_columns": [], "flags": [],
		"error": ""}
	try:
		spreadsheets = ListSpreadsheets(pmcid)
		if spreadsheets is None:
			result["status"] = "not_in_bucket"
			return result
		if not spreadsheets:
			result["status"] = "no_spreadsheet"
			return result
		totalBytes = sum(size for key, size in spreadsheets)
		if totalBytes > SETTINGS["max_article_bytes"]:
			result["status"] = "too_big"
			result["bytes"] = 0
			result["files"] = [{"file": key, "bytes": size, "status": "skipped_article_too_big"} for key, size in spreadsheets]
			return result
		for key, size in spreadsheets:
			if size > SETTINGS["max_file_bytes"]:
				result["files"].append({"file": key, "bytes": size, "status": "skipped_too_big"})
				continue
			targetPath = articleDir / key
			result["bytes"] += DownloadCapped(key, targetPath, SETTINGS["bytes_per_second_per_worker"])
			result["n_files"] += 1
			fileRow, geneColumns, flags = detect.ScreenFile(targetPath, rootDir=articleDir)
			result["files"].append({"file": key, "bytes": size, "status": fileRow["status"]})
			result["gene_columns"].extend(geneColumns)
			result["flags"].extend(flags)
			# Delete right away, the screening result is all we keep
			targetPath.unlink(missing_ok=True)
		result["status"] = "done"
	except Exception as error:
		result["status"] = "error"
		result["error"] = type(error).__name__ + ": " + str(error)[:200]
	finally:
		shutil.rmtree(articleDir, ignore_errors=True)
		result["seconds"] = time.monotonic() - startTime
	return result


# ---------------- main side ----------------

def OpenState(databasePath):
	databasePath.parent.mkdir(parents=True, exist_ok=True)
	connection = sqlite3.connect(databasePath)
	connection.executescript(SCHEMA)
	return connection


def SaveArticle(connection, result):
	# One transaction per article: either all of it is stored or none of it
	tiers = {}
	for flag in result["flags"]:
		tiers[flag["tier"]] = tiers.get(flag["tier"], 0) + 1
	species = sorted({column["species"] for column in result["gene_columns"]})
	with connection:
		previousAttempts = connection.execute("SELECT attempts FROM articles WHERE pmcid = ?", (result["pmcid"],)).fetchone()
		attempts = (previousAttempts[0] if previousAttempts else 0) + 1
		for table in ("files", "gene_columns", "flags"):
			connection.execute(f"DELETE FROM {table} WHERE pmcid = ?", (result["pmcid"],))
		connection.execute("INSERT OR REPLACE INTO articles VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", (result["pmcid"],
			result["status"], result["n_files"], result["bytes"], len(result["gene_columns"]), len(result["flags"]),
			",".join(species), json.dumps(tiers), result["error"], round(result["seconds"], 2), attempts,
			datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")))
		connection.executemany("INSERT INTO files VALUES (?,?,?,?)",
			[(result["pmcid"], row["file"], row["bytes"], row["status"]) for row in result["files"]])
		connection.executemany("INSERT INTO gene_columns VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
			[(row["pmcid"], row["file"], row["sheet"], row["column"], row["species"], row["top20_human"], row["top20_mouse"],
			row["top20_rat"], row["non_empty"], row["symbol_cells"], row["flagged_cells"], row["new_name_cells"],
			row["old_name_text_cells"]) for row in result["gene_columns"]])
		connection.executemany("INSERT INTO flags VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
			[(row["pmcid"], row["file"], row["sheet"], row["cell"], row["species"], row["kind"], row["rule"], row["value"],
			row["guess_original"], row["current_symbol"], row["row_ensembl"], row["row_ensembl_symbols"], row["tier"],
			row["left"], row["right"]) for row in result["flags"]])


def FolderBytes(folder):
	# Workers delete their files while this walks the folder, so a file can vanish between rglob and stat
	totalBytes = 0
	for root, dirs, files in os.walk(folder):
		for fileName in files:
			try:
				totalBytes += os.stat(os.path.join(root, fileName)).st_size
			except OSError:
				continue
	return totalBytes


def LoadCandidates(candidatesPath):
	import csv
	with open(candidatesPath, encoding="utf-8") as handle:
		return [row["pmcid"] for row in csv.DictReader(handle)]


def CountGenePapers(connection, speciesList):
	placeholders = ",".join("?" for unused in speciesList)
	return connection.execute(f"SELECT COUNT(DISTINCT pmcid) FROM gene_columns WHERE species IN ({placeholders})",
		speciesList).fetchone()[0]


def WriteNightLog(sessionStart, reason, doneCount, errorCount, bytesTotal, statusCounts, remaining, connection, databaseName):
	genePapers = connection.execute("SELECT COUNT(DISTINCT pmcid) FROM gene_columns").fetchone()[0]
	flaggedPapers = connection.execute("SELECT COUNT(DISTINCT pmcid) FROM flags").fetchone()[0]
	line = (f"{sessionStart:%Y-%m-%d %H:%M} -> {datetime.datetime.now():%Y-%m-%d %H:%M} | db {databaseName} | "
		f"articles {doneCount} {json.dumps(statusCounts, sort_keys=True)} | errors {errorCount} | "
		f"GB {bytesTotal / 1e9:.2f} | total gene-list papers {genePapers}, flagged {flaggedPapers} | "
		f"remaining {remaining} | stop: {reason}")
	NIGHTS_LOG.parent.mkdir(parents=True, exist_ok=True)
	with open(NIGHTS_LOG, "a", encoding="utf-8") as handle:
		handle.write(line + "\n")
	print(line, flush=True)


def Main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--config", default=str(BASE_DIR / "config.ini"))
	parser.add_argument("--ignore-window", action="store_true")
	parser.add_argument("--wait-for-window", action="store_true", help="sleep until WINDOW_START instead of exiting")
	parser.add_argument("--workers", type=int)
	parser.add_argument("--max-mbps", type=float)
	parser.add_argument("--limit", type=int, default=0, help="stop after this many articles in this session")
	parser.add_argument("--candidates")
	parser.add_argument("--db")
	parser.add_argument("--stop-after-gene-papers", type=int, default=0)
	parser.add_argument("--species", default="human,mouse,rat,mouse_or_rat")
	arguments = parser.parse_args()

	config = configparser.ConfigParser()
	config.read(arguments.config, encoding="utf-8")
	runConfig = config["run"]
	windowStart = ParseTime(runConfig.get("WINDOW_START", "23:30"))
	windowEnd = ParseTime(runConfig.get("WINDOW_END", "07:30"))
	workers = arguments.workers or runConfig.getint("WORKERS", 4)
	maxMbps = arguments.max_mbps if arguments.max_mbps is not None else runConfig.getfloat("MAX_MBPS", 0)
	maxDiskBytes = runConfig.getfloat("MAX_DISK_GB", 20) * 1e9
	minFreeBytes = runConfig.getfloat("MIN_FREE_GB", 30) * 1e9
	maxArticleBytes = runConfig.getfloat("MAX_ARTICLE_GB", 2) * 1e9
	candidatesPath = BASE_DIR / (arguments.candidates or runConfig.get("CANDIDATES"))
	databasePath = BASE_DIR / (arguments.db or runConfig.get("STATE_DB"))
	# Every article in flight can hold up to MAX_ARTICLE_GB, so in-flight work must fit under the disk cap
	if workers * maxArticleBytes > maxDiskBytes * 0.9:
		print(f"WORKERS x MAX_ARTICLE_GB must stay under MAX_DISK_GB, lower one of them", flush=True)
		return 2
	settings = {"max_article_bytes": maxArticleBytes, "max_file_bytes": runConfig.getfloat("MAX_FILE_MB", 300) * 1e6,
		"bytes_per_second_per_worker": maxMbps * 1e6 / workers if maxMbps > 0 else 0}

	lowprio.SetLowPriority()
	if arguments.wait_for_window and not arguments.ignore_window:
		# Started before bed: sleep until the window opens, checking PAUSE once a minute
		while not InWindow(datetime.datetime.now(), windowStart, windowEnd):
			if PAUSE_FILE.exists():
				print("PAUSE file present while waiting, nothing to do", flush=True)
				return 0
			time.sleep(60)
	if not arguments.ignore_window and not InWindow(datetime.datetime.now(), windowStart, windowEnd):
		print(f"outside the window {windowStart:%H:%M}-{windowEnd:%H:%M}, nothing to do", flush=True)
		return 0
	if PAUSE_FILE.exists():
		print("PAUSE file present, nothing to do", flush=True)
		return 0

	# Leftovers from a killed run: the database never saw them, so they get redone
	shutil.rmtree(WORK_DIR, ignore_errors=True)
	WORK_DIR.mkdir(parents=True, exist_ok=True)
	connection = OpenState(databasePath)
	finished = {pmcid for pmcid, status, attempts in connection.execute("SELECT pmcid, status, attempts FROM articles")
		if status != "error" or attempts >= MAX_ATTEMPTS}
	pending = [pmcid for pmcid in LoadCandidates(candidatesPath) if pmcid not in finished]
	speciesList = arguments.species.split(",")
	print(f"{len(pending)} articles left in {candidatesPath.name}, {workers} workers, cap {maxMbps or 'none'} MB/s", flush=True)

	sessionStart = datetime.datetime.now()
	sessionStartUtc = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
	doneCount = 0
	errorCount = 0
	bytesTotal = 0
	statusCounts = {}
	reason = "candidates finished"
	position = 0
	inFlight = {}
	with ProcessPoolExecutor(max_workers=workers, initializer=InitWorker, initargs=(settings,)) as pool:
		try:
			while position < len(pending) or inFlight:
				# Decide whether to start more work; a stop only blocks new articles, the ones in hand finish
				stopReason = ""
				if not arguments.ignore_window and not InWindow(datetime.datetime.now(), windowStart, windowEnd):
					stopReason = "window end"
				elif PAUSE_FILE.exists():
					stopReason = "PAUSE file"
				elif arguments.limit and doneCount + len(inFlight) >= arguments.limit:
					stopReason = f"limit {arguments.limit}"
				elif arguments.stop_after_gene_papers and CountGenePapers(connection, speciesList) >= arguments.stop_after_gene_papers:
					stopReason = f"{arguments.stop_after_gene_papers} gene-list papers ({arguments.species})"
				elif shutil.disk_usage(BASE_DIR).free < minFreeBytes:
					stopReason = "drive free space below MIN_FREE_GB"
				# Room check: what is on disk now plus the worst case of every article in flight and the next one
				usedBytes = FolderBytes(WORK_DIR) + databasePath.stat().st_size
				hasRoom = usedBytes + (len(inFlight) + 1) * maxArticleBytes <= maxDiskBytes
				if not stopReason and not hasRoom and not inFlight:
					stopReason = "disk cap"
				if stopReason:
					reason = stopReason
				canSubmit = not stopReason and hasRoom and position < len(pending) and len(inFlight) < workers
				if canSubmit:
					pmcid = pending[position]
					position += 1
					inFlight[pool.submit(ProcessArticle, pmcid)] = pmcid
					continue
				if not inFlight:
					break
				completed, unused = wait(list(inFlight), timeout=30, return_when=FIRST_COMPLETED)
				for future in completed:
					pmcid = inFlight.pop(future)
					result = future.result()
					SaveArticle(connection, result)
					doneCount += 1
					bytesTotal += result["bytes"]
					statusCounts[result["status"]] = statusCounts.get(result["status"], 0) + 1
					if result["status"] == "error":
						errorCount += 1
					if doneCount % 50 == 0:
						print(f"{datetime.datetime.now():%H:%M:%S} {doneCount} done, {bytesTotal / 1e9:.2f} GB, {statusCounts}", flush=True)
				if stopReason and not inFlight:
					break
		except KeyboardInterrupt:
			reason = "interrupted (Ctrl+C), articles in hand dropped"
			pool.shutdown(wait=False, cancel_futures=True)
	if position >= len(pending) and not inFlight and reason == "candidates finished":
		reason = "candidates finished"
	with connection:
		connection.execute("INSERT INTO sessions VALUES (?,?,?,?,?,?)", (sessionStartUtc,
			datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"), reason, doneCount, errorCount, bytesTotal))
	remaining = len(pending) - doneCount
	WriteNightLog(sessionStart, reason, doneCount, errorCount, bytesTotal, statusCounts, remaining, connection, databasePath.name)
	shutil.rmtree(WORK_DIR, ignore_errors=True)
	return 0


def SelfCheck():
	assert InWindow(datetime.datetime(2026, 11, 5, 23, 45), datetime.time(23, 30), datetime.time(7, 30))
	assert InWindow(datetime.datetime(2026, 11, 5, 3, 0), datetime.time(23, 30), datetime.time(7, 30))
	assert not InWindow(datetime.datetime(2026, 11, 5, 7, 30), datetime.time(23, 30), datetime.time(7, 30))
	assert not InWindow(datetime.datetime(2026, 11, 5, 14, 0), datetime.time(23, 30), datetime.time(7, 30))
	assert InWindow(datetime.datetime(2026, 11, 5, 14, 0), datetime.time(9, 0), datetime.time(17, 0))
	print("self check ok")


if __name__ == "__main__":
	if len(sys.argv) > 1 and sys.argv[1] == "--check":
		SelfCheck()
		sys.exit(0)
	sys.exit(Main())
