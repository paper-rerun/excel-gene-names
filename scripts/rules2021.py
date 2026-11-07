# The 2021 screen (Abeysooriya et al., GeneNameErrors2020/gene_names.sh) redone in Python, to measure what our detector misses
# Their pipeline: readxl with col_types="text" -> TSV -> for each column, head -20 (header included) ->
# grep -cxFf against species gene name lists, >4 hits = gene list -> awk regexes on the column minus its first 2 lines
# Differences we cannot copy: no R here, so readxl's text rendering is imitated (dates become the serial number,
# floats get 17 significant digits with a capital E), and their Ensembl name lists become our current HGNC/MGI/RGD symbols
import csv
import datetime
import functools
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import detect
import lowprio

AWK_RULES = [re.compile(r"^[0-9][0-9][0-9][0-9][0-9]$"), re.compile(r"[0-9][0-9]/[0-9][0-9]/[0-9][0-9]"),
	re.compile(r"[0-9][0-9]-[0-9][0-9]-[0-9][0-9]"),
	re.compile(r"[0-9]\-(JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC|Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec|jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)"),
	re.compile(r"[0-9]\.[0-9][0-9]E\+[0-9][0-9]")]
CURRENT = None


def InitWorker():
	global CURRENT
	lowprio.SetLowPriority()
	symbolSets = detect.LoadSymbolSets()
	CURRENT = {species: {symbol for symbol, current in symbolMap.items() if symbol == current}
		for species, symbolMap in symbolSets.items()}


def RenderLikeReadxl(value):
	if value is None:
		return "NA"
	if isinstance(value, (datetime.datetime, datetime.date)) and not isinstance(value, datetime.time):
		serial = detect.ExcelSerialFromDate(value)
		return str(serial)
	if isinstance(value, bool):
		return "TRUE" if value else "FALSE"
	if isinstance(value, int):
		return str(value)
	if isinstance(value, float):
		if value.is_integer() and abs(value) < 1e15:
			return str(int(value))
		# Excel shows big numbers as 2.31E+19 in a narrow General column, which is what the sci regex was written for
		if abs(value) >= 1e10:
			return format(value, ".2E")
		return repr(value)
	return str(value).replace("\t", " ").replace("\n", " ")


def ScreenFile2021(path, rootDir):
	relativePath = path.relative_to(rootDir).as_posix()
	pmcid = relativePath.split(".")[0]
	geneColumns = 0
	flags = []
	try:
		for sheetName, columns in detect.ReadSheets(path):
			for columnIndex, column in enumerate(columns):
				lines = [RenderLikeReadxl(value) for value in column]
				topLines = lines[:20]
				hits = {species: sum(1 for line in topLines if line in names) for species, names in CURRENT.items()}
				if max(hits.values()) <= 4:
					continue
				geneColumns += 1
				for rowIndex, line in enumerate(lines[2:], start=3):
					if any(rule.search(line) for rule in AWK_RULES):
						flags.append({"pmcid": pmcid, "file": relativePath, "sheet": sheetName,
							"cell": detect.ColumnLetter(columnIndex) + str(rowIndex), "value": line[:40]})
		status = "ok"
	except Exception as error:
		status = "error: " + type(error).__name__
	return pmcid, relativePath, status, geneColumns, flags


def Main():
	supportDir = Path(sys.argv[1])
	outDir = Path(sys.argv[2])
	workers = int(sys.argv[3]) if len(sys.argv) > 3 else 2
	lowprio.SetLowPriority()
	outDir.mkdir(parents=True, exist_ok=True)
	files = sorted(path for path in supportDir.rglob("*") if path.is_file() and re.search(r"\.xlsx?$", path.name, re.I))
	flagRows = []
	columnCounts = []
	with ProcessPoolExecutor(max_workers=workers, initializer=InitWorker) as pool:
		for pmcid, relativePath, status, geneColumns, flags in pool.map(functools.partial(ScreenFile2021, rootDir=supportDir), files, chunksize=2):
			columnCounts.append({"pmcid": pmcid, "file": relativePath, "status": status, "gene_columns": geneColumns})
			flagRows.extend(flags)
	for name, rows in (("files_2021.csv", columnCounts), ("flags_2021.csv", flagRows)):
		with open(outDir / name, "w", newline="", encoding="utf-8") as handle:
			writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else ["empty"])
			writer.writeheader()
			writer.writerows(rows)
	print(len(files), "files", len(flagRows), "flags")


if __name__ == "__main__":
	if len(sys.argv) > 1 and sys.argv[1] == "--check":
		assert RenderLikeReadxl(datetime.datetime(2011, 3, 1)) == "40603"
		assert RenderLikeReadxl(2.310009e19) == "2.31E+19"
		assert AWK_RULES[4].search(RenderLikeReadxl(2.310009e19))
		assert RenderLikeReadxl(40603.0) == "40603"
		print("self check ok")
		sys.exit(0)
	sys.exit(Main())
