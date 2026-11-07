# Draw the precision sample: 100 auto-flagged papers from state/full.sqlite, stratified by species x tier,
# plus the list of manual-only papers. Writes data/precision_sample.csv and data/manual_only.csv
# Stratum of a paper: the species of its auto-flagged list (human > mouse > rat > mouse_or_rat when several)
# and tier A if that list has an Ensembl-confirmed cell, B otherwise
import csv
import random
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze_db

BASE_DIR = Path(__file__).resolve().parent.parent
SEED = 20261107
SAMPLE_SIZE = 100
MIN_PER_STRATUM = 5
SPECIES_ORDER = ["human", "mouse", "rat", "mouse_or_rat"]


def LoadFlags(connection):
	return connection.execute("SELECT pmcid, species, tier, rule, guess_original, value, file, sheet, cell FROM flags").fetchall()


def PaperStrata(flags):
	verdicts = analyze_db.PaperVerdicts(flags)
	hasTierA = {(row[0], row[1]) for row in flags if row[2] == "A_ensembl_match"}
	strata = {}
	manualOnly = []
	papers = sorted({row[0] for row in flags})
	for pmcid in papers:
		if verdicts.get((pmcid, "any")) != "auto":
			manualOnly.append(pmcid)
			continue
		autoSpecies = [species for species in SPECIES_ORDER if verdicts.get((pmcid, species)) == "auto"]
		# Auto only when the lists are taken together (two lone 5-digit numbers in two lists): first flagged species
		flaggedSpecies = [species for species in SPECIES_ORDER if (pmcid, species) in verdicts]
		species = (autoSpecies or flaggedSpecies)[0]
		tier = "A" if (pmcid, species) in hasTierA else "B"
		strata[pmcid] = (species, tier)
	return strata, manualOnly


def Allocate(populationByStratum, total, minimum):
	# Proportional, with a floor so small strata (rat, mouse_or_rat) get looked at; capped at the stratum size
	allocation = {stratum: min(size, minimum) for stratum, size in populationByStratum.items()}
	remaining = total - sum(allocation.values())
	bigTotal = sum(size - allocation[stratum] for stratum, size in populationByStratum.items())
	shares = {stratum: (size - allocation[stratum]) / bigTotal * remaining for stratum, size in populationByStratum.items()}
	for stratum, share in shares.items():
		allocation[stratum] += int(share)
	# Hand out what rounding left to the largest remainders
	leftover = total - sum(allocation.values())
	for stratum in sorted(shares, key=lambda key: -(shares[key] - int(shares[key])))[:leftover]:
		allocation[stratum] += 1
	return allocation


def Main():
	connection = sqlite3.connect(BASE_DIR / "state" / "full.sqlite")
	flags = LoadFlags(connection)
	strata, manualOnly = PaperStrata(flags)
	population = defaultdict(list)
	for pmcid, stratum in strata.items():
		population[stratum].append(pmcid)
	populationSizes = {stratum: len(members) for stratum, members in population.items()}
	allocation = Allocate(populationSizes, SAMPLE_SIZE, MIN_PER_STRATUM)
	randomGenerator = random.Random(SEED)
	rows = []
	for stratum in sorted(population):
		members = sorted(population[stratum])
		for pmcid in randomGenerator.sample(members, allocation[stratum]):
			rows.append({"pmcid": pmcid, "species": stratum[0], "tier": stratum[1],
				"stratum_population": populationSizes[stratum], "stratum_sample": allocation[stratum]})
	with open(BASE_DIR / "data" / "precision_sample.csv", "w", newline="", encoding="utf-8") as handle:
		writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
		writer.writeheader()
		writer.writerows(rows)
	with open(BASE_DIR / "data" / "manual_only.csv", "w", newline="", encoding="utf-8") as handle:
		writer = csv.writer(handle)
		writer.writerow(["pmcid"])
		writer.writerows([pmcid] for pmcid in manualOnly)
	print("seed", SEED)
	for stratum in sorted(population):
		print(f"stratum {stratum[0]:13s} tier {stratum[1]}: population {populationSizes[stratum]}, sample {allocation[stratum]}")
	print(f"auto papers {len(strata)}, manual-only papers {len(manualOnly)}")


if __name__ == "__main__":
	sys.exit(Main())
