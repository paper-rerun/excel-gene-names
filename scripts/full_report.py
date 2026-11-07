# Numbers for REPORT.md from the full run, the precision sample and the manual-only review
# Usage: full_report.py  (reads state/full.sqlite, data/full_review.csv, data/candidates_full_genom.csv)
# Writes data/full_report_numbers.txt and data/per_paper_table.csv
import csv
import math
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze_db
import precision_sample

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
REFERENCE_2021 = {"all": (3436, 11117), "human": (2419, 7936), "mouse": (609, 1577), "rat": (68, 327)}
REFERENCE_2016 = (704, 3597)
MIN_JOURNAL_PAPERS = 30


def Wilson(successes, total, z=1.96):
	if total == 0:
		return 0.0, 0.0, 0.0
	proportion = successes / total
	denominator = 1 + z * z / total
	center = (proportion + z * z / (2 * total)) / denominator
	halfWidth = z * math.sqrt(proportion * (1 - proportion) / total + z * z / (4 * total * total)) / denominator
	return proportion, center - halfWidth, center + halfWidth


def Rate(successes, total):
	proportion, low, high = Wilson(successes, total)
	return f"{successes:.0f}/{total} = {proportion * 100:.1f}% (IC95 {low * 100:.1f}-{high * 100:.1f}%)"


def Family(guessOriginal):
	# Old-symbol family of the first candidate; the candidates inside one date are ambiguous (MARCH1/MARC1)
	first = guessOriginal.split("/")[0].upper()
	letters = "".join(character for character in first if character.isalpha())
	if letters == "SEP" and first.endswith("15"):
		return "SEP15 (SELENOF)"
	return {"SEPT": "SEPT* (SEPTIN*)", "SEP": "SEPT* (SEPTIN*)", "MARCH": "MARCH* (MARCHF*)", "MAR": "MARCH* (MARCHF*)",
		"MARC": "MARC1/2 (MTARC1/2)", "DEC": "DEC1 (DELEC1)", "AGO": "AGO1-4, Excel em pt/es/it (ago = agosto)",
		"MEI": "MEI1/MEI4, Excel em holandes (mei = maio)"}.get(letters, "outros: " + letters)


def Main():
	connection = sqlite3.connect(BASE_DIR / "state" / "full.sqlite")
	flags = precision_sample.LoadFlags(connection)
	geneColumns = connection.execute("SELECT pmcid, species FROM gene_columns").fetchall()
	candidates = {row["pmcid"]: row for row in csv.DictReader(open(DATA_DIR / "candidates_full_genom.csv", encoding="utf-8"))}
	reviews = {row["pmcid"]: row for row in csv.DictReader(open(DATA_DIR / "full_review.csv", encoding="utf-8"))}
	strata, manualOnly = precision_sample.PaperStrata(flags)
	verdicts = analyze_db.PaperVerdicts(flags)
	speciesByPaper = defaultdict(set)
	for pmcid, species in geneColumns:
		speciesByPaper[pmcid].add(species)
	listPapers = set(speciesByPaper)
	lines = []

	# 1. Precision by stratum
	strataPopulation = Counter(strata.values())
	sampled = [row for row in reviews.values() if row["set"] == "precision_sample"]
	lines.append("== precision of the automatic call (seeded sample, stratified by species x tier)")
	precisionByStratum = {}
	for stratum in sorted(strataPopulation):
		reviewed = [row for row in sampled if (row["species"], row["tier"]) == stratum]
		real = sum(1 for row in reviewed if row["verdict"] == "real")
		precisionByStratum[stratum] = (real, len(reviewed))
		lines.append(f"  {stratum[0]:13s} {stratum[1]}: population {strataPopulation[stratum]}, " + Rate(real, len(reviewed)))
	realAll = sum(1 for row in sampled if row["verdict"] == "real")
	lines.append("  unweighted, all 100: " + Rate(realAll, len(sampled)))
	totalAuto = sum(strataPopulation.values())
	weighted = 0.0
	variance = 0.0
	for stratum, size in strataPopulation.items():
		real, reviewedCount = precisionByStratum[stratum]
		share = size / totalAuto
		precision = real / reviewedCount
		weighted += share * precision
		correction = (size - reviewedCount) / (size - 1) if size > 1 else 0
		variance += share * share * precision * (1 - precision) / reviewedCount * correction
	lines.append(f"  weighted by stratum size: {weighted * 100:.1f}% (normal IC95 {max(0, weighted - 1.96 * math.sqrt(variance)) * 100:.1f}-{min(1, weighted + 1.96 * math.sqrt(variance)) * 100:.1f}%)")
	for tier in ("A", "B"):
		real = sum(1 for row in sampled if row["tier"] == tier and row["verdict"] == "real")
		total = sum(1 for row in sampled if row["tier"] == tier)
		lines.append(f"  tier {tier} pooled: " + Rate(real, total))
	falseRows = [row for row in sampled if row["verdict"] != "real"]
	lines.append(f"  false positives in the sample: {[(row['pmcid'], row['note']) for row in falseRows]}")

	# 2. Manual-only papers
	manualRows = [reviews[pmcid] for pmcid in manualOnly if pmcid in reviews]
	manualCounts = Counter(row["verdict"] for row in manualRows)
	lines.append("")
	lines.append(f"== manual-only papers: {len(manualOnly)}, reviewed {len(manualRows)}: {dict(manualCounts)}")
	lines.append(f"  reasons: {dict(Counter(row['category'] for row in manualRows))}")
	manualReal = {row["pmcid"] for row in manualRows if row["verdict"] == "real"}
	manualProbable = {row["pmcid"] for row in manualRows if row["verdict"] == "probable"}
	lines.append(f"  real by reason: {dict(Counter(row['category'] for row in manualRows if row['verdict'] == 'real'))}")
	lines.append(f"  if the {len(manualProbable)} probable ones were real too, add {len(manualProbable)} papers to the estimate")

	# 3. Estimated affected papers: each auto paper weighs its stratum precision, manual ones count if confirmed
	def PaperWeight(pmcid):
		if pmcid in manualReal:
			return 1.0
		if pmcid in strata:
			if pmcid in reviews and reviews[pmcid]["set"] == "precision_sample":
				return 1.0 if reviews[pmcid]["verdict"] == "real" else 0.0
			real, reviewedCount = precisionByStratum[strata[pmcid]]
			return real / reviewedCount
		return 0.0

	def GroupLine(label, papers):
		papers = [pmcid for pmcid in papers if pmcid in listPapers]
		autoCount = sum(1 for pmcid in papers if pmcid in strata)
		estimate = sum(PaperWeight(pmcid) for pmcid in papers)
		return (f"{label}: gene-list papers {len(papers)} | auto {Rate(autoCount, len(papers))} | "
			f"estimated {Rate(round(estimate), len(papers))}")

	lines.append("")
	lines.append("== rates (auto = detector call; estimated = auto x stratum precision + confirmed manual-only)")
	lines.append(GroupLine("all species", listPapers))
	for species in ("human", "mouse", "rat", "mouse_or_rat"):
		# Per species: the paper has a list of that species and the auto call is for that species (or manual confirmed)
		speciesPapers = [pmcid for pmcid in listPapers if species in speciesByPaper[pmcid]]
		autoSpecies = [pmcid for pmcid in speciesPapers if verdicts.get((pmcid, species)) == "auto"]
		estimate = sum(PaperWeight(pmcid) for pmcid in autoSpecies) + sum(1 for pmcid in speciesPapers
			if pmcid in manualReal and pmcid not in autoSpecies)
		lines.append(f"{species}: gene-list papers {len(speciesPapers)} | auto {Rate(len(autoSpecies), len(speciesPapers))} | "
			f"estimated {Rate(round(estimate), len(speciesPapers))}")
	for name, (affected, total) in REFERENCE_2021.items():
		lines.append(f"2021 reference {name}: " + Rate(affected, total))
	lines.append("2016 reference all: " + Rate(*REFERENCE_2016))
	estimateAll = sum(PaperWeight(pmcid) for pmcid in listPapers)
	refAffected, refTotal = REFERENCE_2021["all"]
	pooled = (estimateAll + refAffected) / (len(listPapers) + refTotal)
	zValue = (estimateAll / len(listPapers) - refAffected / refTotal) / math.sqrt(pooled * (1 - pooled) * (1 / len(listPapers) + 1 / refTotal))
	lines.append(f"z vs 2021 (estimated, all species): z={zValue:.1f}, p={math.erfc(abs(zValue) / math.sqrt(2)):.2g}")
	lines.append("all species, estimate counting the probable manual ones as real: " + Rate(round(estimateAll) + len(manualProbable), len(listPapers)))

	lines.append("")
	lines.append("== by publication year")
	for year in ("2024", "2025", "2026"):
		lines.append(GroupLine(year, [pmcid for pmcid in listPapers if candidates[pmcid]["pub_year"] == year]))
	lines.append("")
	lines.append(f"== by journal (>= {MIN_JOURNAL_PAPERS} gene-list papers)")
	journalCounts = Counter(candidates[pmcid]["journal"] for pmcid in listPapers)
	journalLines = []
	for journal, count in journalCounts.items():
		if count >= MIN_JOURNAL_PAPERS:
			papers = [pmcid for pmcid in listPapers if candidates[pmcid]["journal"] == journal]
			estimate = sum(PaperWeight(pmcid) for pmcid in papers)
			journalLines.append((estimate / count, GroupLine(journal, papers)))
	for unused, line in sorted(journalLines, reverse=True):
		lines.append("  " + line)

	# Same journal, then and now: Abeysooriya 2021 Table 5 (articles with Excel gene lists, affected), 2014-2020
	table5 = {"Nat Commun": (920, 345), "PLoS One": (946, 244), "Sci Rep": (767, 227), "BMC Genomics": (660, 166),
		"PLoS Genet": (448, 134), "Oncotarget": (326, 107), "Front Genet": (313, 94), "Elife": (243, 89),
		"Proc Natl Acad Sci U S A": (155, 73), "Cell Rep": (158, 71), "Genome Biol": (193, 66), "Nature": (118, 52),
		"Nat Genet": (140, 48), "Genome Med": (137, 44), "PeerJ": (137, 39), "Cell": (74, 39), "Clin Epigenetics": (109, 38),
		"Nucleic Acids Res": (120, 36), "BMC Med Genomics": (117, 31), "Front Oncol": (85, 31), "Transl Psychiatry": (73, 29),
		"BMC Cancer": (105, 28), "PLoS Pathog": (80, 27), "Commun Biol": (74, 27), "PLoS Biol": (66, 26), "Aging": (56, 26),
		"EBioMedicine": (51, 26), "Epigenetics Chromatin": (64, 25), "PLoS Comput Biol": (97, 24), "Oncogene": (53, 22),
		"iScience": (58, 20), "Sci Adv": (56, 20), "BMC Bioinformatics": (77, 19), "G3 (Bethesda)": (74, 15),
		"Hum Mol Genet": (53, 15), "BMC Plant Biol": (52, 6), "Front Plant Sci": (75, 5)}
	lines.append("")
	lines.append(f"== same journal in 2021 Table 5 (journals with >= {MIN_JOURNAL_PAPERS} gene-list papers here)")
	lowerCount = 0
	compared = 0
	sumNow = 0.0
	sumNowTotal = 0
	sumThenAffected = 0
	sumThenTotal = 0
	for journal, (thenTotal, thenAffected) in table5.items():
		papers = [pmcid for pmcid in listPapers if candidates[pmcid]["journal"] == journal]
		if len(papers) < MIN_JOURNAL_PAPERS:
			continue
		estimate = sum(PaperWeight(pmcid) for pmcid in papers)
		compared += 1
		lowerCount += int(estimate / len(papers) < thenAffected / thenTotal)
		sumNow += estimate
		sumNowTotal += len(papers)
		sumThenAffected += thenAffected
		sumThenTotal += thenTotal
		lines.append(f"  {journal}: now {estimate / len(papers) * 100:.1f}% of {len(papers)} | 2014-2020 {thenAffected / thenTotal * 100:.1f}% of {thenTotal}")
	lines.append(f"  lower now in {lowerCount} of {compared} journals; pooled now {Rate(round(sumNow), sumNowTotal)} vs 2014-2020 {Rate(sumThenAffected, sumThenTotal)}")

	# 4. Old symbols: papers per family, from auto and confirmed papers
	lines.append("")
	lines.append("== old symbol families behind the errors (papers with >= 1 cell of that family; auto + confirmed manual)")
	affectedPapers = set(strata) | manualReal
	papersByFamily = defaultdict(set)
	localeExamples = defaultdict(set)
	for pmcid, species, tier, rule, guess, value, fileName, sheet, cell in flags:
		if pmcid not in affectedPapers or not guess or tier.startswith("C_"):
			continue
		papersByFamily[Family(guess)].add(pmcid)
		if "locale" in rule:
			localeExamples[rule].add(value)
	for family, paperSet in sorted(papersByFamily.items(), key=lambda item: -len(item[1])):
		lines.append(f"  {family}: {len(paperSet)} ({len(paperSet) / len(affectedPapers) * 100:.1f}% of {len(affectedPapers)})")
	lines.append(f"  text cells in non-English month format: {{rule: examples}} {dict((rule, sorted(values)[:8]) for rule, values in localeExamples.items())}")
	# Same count without 5-digit numbers, which can be IDs that happen to read as a date
	papersByFamilyDates = defaultdict(set)
	for pmcid, species, tier, rule, guess, value, fileName, sheet, cell in flags:
		if pmcid in affectedPapers and guess and not tier.startswith("C_") and rule != "five_digit":
			papersByFamilyDates[Family(guess)].add(pmcid)
	lines.append("  same, date cells and text only (no 5-digit numbers):")
	for family, paperSet in sorted(papersByFamilyDates.items(), key=lambda item: -len(item[1])):
		lines.append(f"    {family}: {len(paperSet)} ({len(paperSet) / len(affectedPapers) * 100:.1f}%)")
	cellKinds = Counter(rule for pmcid, species, tier, rule, *unused in flags if pmcid in affectedPapers and not tier.startswith("C_"))
	lines.append(f"  cells by kind in affected papers: {dict(cellKinds)}")

	# 5. Run facts
	articleCounts = dict(connection.execute("SELECT status, COUNT(*) FROM articles GROUP BY status").fetchall())
	totalBytes = connection.execute("SELECT SUM(bytes) FROM articles").fetchone()[0]
	lines.append("")
	lines.append(f"== run: {sum(articleCounts.values())} articles {articleCounts}, {totalBytes / 1e9:.1f} GB streamed")
	errorRows = connection.execute("SELECT pmcid, error FROM articles WHERE status = 'error'").fetchall()
	lines.append(f"  errors: {errorRows}")
	fileStatus = Counter(status.split(":")[0] for status, in connection.execute("SELECT status FROM files"))
	lines.append(f"  files: {dict(fileStatus)}")

	(DATA_DIR / "full_report_numbers.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
	print("\n".join(lines))

	# 6. Per-paper table, not for publication yet
	tierRank = {"A_ensembl_match": 0, "B_guess_no_id": 1}
	examples = {}
	cellCounts = Counter()
	for pmcid, species, tier, rule, guess, value, fileName, sheet, cell in flags:
		# Manual-confirmed papers can rest on cells the guesser could not name (TAMM41, MARS1, text dates)
		isUsable = not tier.startswith("C_") and guess
		if pmcid not in affectedPapers or not (isUsable or pmcid in manualReal):
			continue
		cellCounts[pmcid] += 1
		rank = (tierRank.get(tier, 9), 0 if rule == "date_cell" else 1)
		if pmcid not in examples or rank < examples[pmcid][0]:
			examples[pmcid] = (rank, species, tier, rule, guess, value, fileName, sheet, cell)
	currentByFlag = {}
	for pmcid, current, guess in connection.execute("SELECT pmcid, current_symbol, guess_original FROM flags"):
		currentByFlag[(pmcid, guess)] = current
	with open(DATA_DIR / "per_paper_table.csv", "w", newline="", encoding="utf-8") as handle:
		writer = csv.writer(handle)
		writer.writerow(["pmcid", "doi", "journal", "pub_year", "call", "species", "evidence", "n_cells", "file", "sheet",
			"cell", "value_in_file", "old_symbol_candidates", "current_symbol_candidates"])
		for pmcid in sorted(affectedPapers):
			if pmcid not in examples:
				continue
			unused, species, tier, rule, guess, value, fileName, sheet, cell = examples[pmcid]
			call = "manual_confirmed" if pmcid in manualReal else ("reviewed_real" if pmcid in reviews and reviews[pmcid]["verdict"] == "real" else "auto")
			if pmcid in reviews and reviews[pmcid]["verdict"] != "real":
				call = "reviewed_false_positive"
			candidate = candidates[pmcid]
			writer.writerow([pmcid, candidate["doi"], candidate["journal"], candidate["pub_year"], call, species,
				"ensembl_id" if tier == "A_ensembl_match" else "date_pattern", cellCounts[pmcid], fileName.split("/")[-1],
				sheet, cell, value, guess, currentByFlag.get((pmcid, guess), "")])


if __name__ == "__main__":
	sys.exit(Main())
