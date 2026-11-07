# Independent check of flagged cells: when an Ensembl gene ID sits next to the date, does it belong to the guessed gene?
import csv
import re
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
symbolByEnsembl = {}
with open(BASE_DIR / "dl" / "ref" / "hgnc_complete_set.txt", encoding="utf-8") as handle:
	for row in csv.DictReader(handle, delimiter="\t"):
		if row.get("ensembl_gene_id"):
			symbolByEnsembl[row["ensembl_gene_id"]] = row["symbol"]
checked = 0
agree = 0
disagree = []
for flag in csv.DictReader(open(BASE_DIR / "data" / "flags.csv", encoding="utf-8")):
	ensemblIds = re.findall(r"ENSG\d{11}", flag["left"] + " " + flag["right"])
	if not ensemblIds or not flag["current_symbol"]:
		continue
	checked += 1
	neighbourSymbols = {symbolByEnsembl.get(ensemblId, "?") for ensemblId in ensemblIds}
	if neighbourSymbols & set(flag["current_symbol"].split("/")):
		agree += 1
	else:
		disagree.append((flag["pmcid"], flag["cell"], flag["value"], flag["current_symbol"], sorted(neighbourSymbols)))
print(f"flagged cells with an Ensembl ID beside them: {checked}; guess matches the ID: {agree}")
for item in disagree:
	print("mismatch", item)

# Which genes the Ensembl-backed cells really are (current symbol: papers)
from collections import defaultdict
papersPerGene = defaultdict(set)
for flag in csv.DictReader(open(BASE_DIR / "data" / "flags.csv", encoding="utf-8")):
	for ensemblId in re.findall(r"ENSG\d{11}", flag["left"] + " " + flag["right"]):
		gene = symbolByEnsembl.get(ensemblId, "?")
		if gene in flag["current_symbol"].split("/"):
			papersPerGene[gene].add(flag["pmcid"])
print("genes confirmed by Ensembl ID (gene: papers):", ", ".join(f"{gene} {len(paperSet)}" for gene, paperSet in sorted(papersPerGene.items(), key=lambda item: (-len(item[1]), item[0]))))
