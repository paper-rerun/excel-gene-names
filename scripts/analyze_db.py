# Rates from a run_batch.py database, by species and by evidence tier
# Usage: analyze_db.py state/lab_mouse.sqlite [more.sqlite ...] > data/<name>_summary.txt
import json
import math
import sqlite3
import sys
from collections import Counter, defaultdict


def WilsonInterval(successes, total, z=1.96):
	if total == 0:
		return 0.0, 0.0, 0.0
	proportion = successes / total
	denominator = 1 + z * z / total
	center = (proportion + z * z / (2 * total)) / denominator
	halfWidth = z * math.sqrt(proportion * (1 - proportion) / total + z * z / (4 * total * total)) / denominator
	return proportion, center - halfWidth, center + halfWidth


def FormatRate(successes, total):
	proportion, low, high = WilsonInterval(successes, total)
	return f"{successes}/{total} = {proportion * 100:.1f}% (IC95 Wilson {low * 100:.1f}-{high * 100:.1f}%)"


def PaperVerdicts(flags):
	# (pmcid, species) -> "auto" or "manual"; species "any" covers all lists of the paper
	evidence = defaultdict(lambda: {"A": 0, "B_date": 0, "B_number": 0, "other": 0})
	for pmcid, species, tier, rule, guess, value, fileName, sheet, cell in flags:
		for key in ((pmcid, species), (pmcid, "any")):
			if tier == "A_ensembl_match":
				evidence[key]["A"] += 1
			elif tier == "B_guess_no_id" and rule == "five_digit":
				evidence[key]["B_number"] += 1
			elif tier == "B_guess_no_id":
				evidence[key]["B_date"] += 1
			else:
				evidence[key]["other"] += 1
	verdicts = {}
	for key, counts in evidence.items():
		isAuto = counts["A"] > 0 or counts["B_date"] > 0 or counts["B_number"] >= 2
		verdicts[key] = "auto" if isAuto else "manual"
	return verdicts


def ReadManual(path):
	# CSV with pmcid,verdict (real / false_positive)
	import csv
	with open(path, encoding="utf-8") as handle:
		return {row["pmcid"]: row["verdict"] for row in csv.DictReader(handle)}


def Main():
	articles = []
	geneColumns = []
	flags = []
	sessions = []
	for databasePath in [argument for argument in sys.argv[1:] if argument.endswith(".sqlite")]:
		connection = sqlite3.connect(databasePath)
		articles += connection.execute("SELECT pmcid, status, n_files, bytes, seconds FROM articles").fetchall()
		geneColumns += connection.execute("SELECT pmcid, species FROM gene_columns").fetchall()
		flags += connection.execute("SELECT pmcid, species, tier, rule, guess_original, value, file, sheet, cell FROM flags").fetchall()
		sessions += connection.execute("SELECT started_utc, ended_utc, reason, articles, errors, bytes FROM sessions").fetchall()
	print("inputs:", " ".join(sys.argv[1:]))
	print("articles processed:", len(articles), dict(Counter(row[1] for row in articles)))
	print(f"spreadsheet GB downloaded: {sum(row[3] for row in articles) / 1e9:.2f}")
	print(f"seconds per article (worker time): {sum(row[4] or 0 for row in articles) / max(len(articles), 1):.2f}")
	withSheets = [row for row in articles if row[1] in ("done", "too_big")]
	print(f"articles with spreadsheets: {len(withSheets)} of {len(articles)}")
	speciesByPaper = defaultdict(set)
	for pmcid, species in geneColumns:
		speciesByPaper[pmcid].add(species)
	print("papers with a gene list:", len(speciesByPaper))
	print("flags by tier:", dict(Counter(row[2] for row in flags)))
	# Paper verdict per species: auto = an A cell, a B date cell, or 2+ B five-digit numbers (a lone 5-digit
	# number was a false positive in the rat sample); a paper with flags but no auto verdict goes to manual review
	verdicts = PaperVerdicts(flags)
	manualPaths = [argument for argument in sys.argv[1:] if argument.endswith(".csv")]
	manualVerdicts = ReadManual(manualPaths[0]) if manualPaths else {}
	for species in ("any", "human", "mouse", "rat", "mouse_or_rat"):
		papers = [pmcid for pmcid, speciesSet in speciesByPaper.items() if species == "any" or species in speciesSet]
		if not papers:
			continue
		autoPapers = [pmcid for pmcid in papers if verdicts.get((pmcid, species)) == "auto"]
		manualPapers = [pmcid for pmcid in papers if verdicts.get((pmcid, species)) == "manual"]
		anyFlag = autoPapers + manualPapers
		print(f"{species}: papers with list {len(papers)} | auto {FormatRate(len(autoPapers), len(papers))} | "
			f"any flag {FormatRate(len(anyFlag), len(papers))} | waiting for manual {len(manualPapers)}")
		if manualVerdicts:
			confirmed = [pmcid for pmcid in anyFlag if manualVerdicts.get(pmcid) == "real"]
			autoWrong = [pmcid for pmcid in autoPapers if manualVerdicts.get(pmcid) != "real"]
			manualReal = [pmcid for pmcid in manualPapers if manualVerdicts.get(pmcid) == "real"]
			print(f"   manual check: confirmed {FormatRate(len(confirmed), len(papers))} | auto yes, manual no: "
				f"{len(autoWrong)} of {len(autoPapers)} {autoWrong} | real among manual-only: {len(manualReal)} of {len(manualPapers)}")
	print("sessions:")
	for session in sessions:
		print("  ", session)


if __name__ == "__main__":
	sys.exit(Main())
