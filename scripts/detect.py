# Screen the downloaded spreadsheets for gene names that Excel turned into dates or numbers
# Gene column rule and base regexes copied from Abeysooriya et al. 2021 (GeneNameErrors2020/gene_names.sh):
# a column is a gene list if more than 4 of its first 20 cells are official symbols of one species
# (here: 5 distinct current symbols, see PILOTO.md for why)
import csv
import datetime
import functools
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

# -I drops the script folder from sys.path, lowprio lives next to this file
sys.path.insert(0, str(Path(__file__).resolve().parent))

BASE_DIR = Path(__file__).resolve().parent.parent
REF_DIR = BASE_DIR / "dl" / "ref"
SUPP_DIR = BASE_DIR / "dl" / "supp"
DATA_DIR = BASE_DIR / "data"
MIN_SYMBOLS_IN_TOP_20 = 5
MAX_CELLS_PER_SHEET = 5_000_000
NEW_NAME_PATTERN = re.compile(r"^(SEPTIN|MARCHF)[0-9]+$", re.I)
OLD_NAME_PATTERN = re.compile(r"^(SEPT|MARCH)[0-9]+$", re.I)

MONTHS_EN = "jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec"
# Abbreviations Excel shows in other locales (pt, es, fr, de, it, nl)
MONTHS_OTHER = "fev|abr|mai|ago|set|out|dez|ene|dic|févr|fév|avr|juil|déc|mär|mrz|okt|giu|lug|ott|mei|maa|sept"

# Base rules: the awk regexes from the 2021 screen
BASE_RULES = [
	("five_digit", re.compile(r"^[0-9]{5}$")),
	("dd/dd/dd", re.compile(r"[0-9][0-9]/[0-9][0-9]/[0-9][0-9]")),
	("dd-dd-dd", re.compile(r"[0-9][0-9]-[0-9][0-9]-[0-9][0-9]")),
	("d-Mon", re.compile(r"[0-9]-(" + MONTHS_EN + ")", re.I)),
	("sci_E", re.compile(r"[0-9]\.[0-9][0-9]E\+[0-9][0-9]", re.I)),
]
# Extra rules for formats the 2021 regex misses
EXTRA_RULES = [
	# Month token must be the bare month and the separator "-", "/" or space: with "." or a longer word these
	# matched R's make.unique names in the lab sample ("Jun.1", "Set.1", "Marcks.1")
	("Mon-d", re.compile(r"^(" + MONTHS_EN + r"|sept)[-/ ][0-9]{1,4}$", re.I)),
	("d-Mon_locale", re.compile(r"^[0-9]{1,2}[-/ ](" + MONTHS_OTHER + r")\.?([-/ ][0-9]{2,4})?$", re.I)),
	("Mon_locale-d", re.compile(r"^(" + MONTHS_OTHER + r")\.?[-/ ][0-9]{1,4}$", re.I)),
	("numeric_date", re.compile(r"^[0-9]{1,2}[./][0-9]{1,2}[./][0-9]{2,4}$")),
	("iso_date", re.compile(r"^[0-9]{4}[-/.][0-9]{1,2}[-/.][0-9]{1,2}([ T][0-9:]+)?$")),
	("sci_E_any", re.compile(r"^[0-9](\.[0-9]+)?E\+?[0-9]+$", re.I)),
]
# MARC1/MARC2 (now MTARC1/2), AGO1-4 (Portuguese "ago" = August) and MEI1/MEI4 (Dutch "mei" = May) also turn into dates
MONTH_PREFIXES = {1: ["JAN"], 2: ["FEB"], 3: ["MARCH", "MAR", "MARC"], 4: ["APR"], 5: ["MAY", "MEI"], 6: ["JUN"], 7: ["JUL"],
	8: ["AUG", "AGO"], 9: ["SEPT", "SEP"], 10: ["OCT"], 11: ["NOV"], 12: ["DEC"]}


def LoadSymbolSets():
	# Each set maps symbol -> current symbol, so old names (SEPT2) can be named after the fix (SEPTIN2)
	human = {}
	humanAliases = []
	with open(REF_DIR / "hgnc_complete_set.txt", encoding="utf-8") as handle:
		for row in csv.DictReader(handle, delimiter="\t"):
			if row["status"] != "Approved":
				continue
			human[row["symbol"]] = row["symbol"]
			for oldSymbol in row["prev_symbol"].strip('"').split("|"):
				if oldSymbol and oldSymbol not in human:
					human[oldSymbol] = row["symbol"]
			# Aliases too, only for naming the original (SEP15 is an alias of SELENOF), detection uses current symbols
			for aliasSymbol in row["alias_symbol"].strip('"').split("|"):
				if aliasSymbol:
					humanAliases.append((aliasSymbol, row["symbol"]))
	# Aliases go last so a previous symbol wins (DEC1 is DELEC1's old name and a BHLHE40 alias)
	for aliasSymbol, currentSymbol in humanAliases:
		human.setdefault(aliasSymbol, currentSymbol)
	mouse = {}
	with open(REF_DIR / "MRK_List1.rpt", encoding="utf-8") as handle:
		for row in csv.DictReader(handle, delimiter="\t"):
			if row["Marker Type"] != "Gene":
				continue
			if row["Status"] == "O":
				mouse[row["Marker Symbol"]] = row["Marker Symbol"]
			elif row["Status"] == "W" and row["Current Marker Symbol (if withdrawn)"]:
				mouse.setdefault(row["Marker Symbol"], row["Current Marker Symbol (if withdrawn)"])
	rat = {}
	with open(REF_DIR / "GENES_RAT.txt", encoding="utf-8") as handle:
		lines = (line for line in handle if not line.startswith("#"))
		for row in csv.DictReader(lines, delimiter="\t"):
			rat[row["SYMBOL"]] = row["SYMBOL"]
			for oldSymbol in (row.get("OLD_SYMBOL") or "").split(";"):
				if oldSymbol and oldSymbol not in rat:
					rat[oldSymbol] = row["SYMBOL"]
	# Drop one-letter and numeric tokens (RGD lists "C" and "1" as old symbols), they match genotype and 0/1 columns
	validSymbol = re.compile(r"^[A-Za-z][A-Za-z0-9._@-]+$")
	return {species: {symbol: current for symbol, current in symbolMap.items() if validSymbol.match(symbol)}
		for species, symbolMap in (("human", human), ("mouse", mouse), ("rat", rat))}


def LoadEnsemblMap():
	# Ensembl gene ID -> current symbol, human, mouse and rat in one map (the prefixes do not collide)
	symbolByEnsembl = {}
	with open(REF_DIR / "hgnc_complete_set.txt", encoding="utf-8") as handle:
		for row in csv.DictReader(handle, delimiter="\t"):
			if row.get("ensembl_gene_id"):
				symbolByEnsembl[row["ensembl_gene_id"]] = row["symbol"]
	with open(REF_DIR / "MRK_ENSEMBL.rpt", encoding="utf-8") as handle:
		for fields in csv.reader(handle, delimiter="\t"):
			if len(fields) > 5 and fields[5].startswith("ENSMUSG"):
				symbolByEnsembl[fields[5]] = fields[1]
	with open(REF_DIR / "GENES_RAT.txt", encoding="utf-8") as handle:
		lines = (line for line in handle if not line.startswith("#"))
		for row in csv.DictReader(lines, delimiter="\t"):
			for ensemblId in (row.get("ENSEMBL_ID") or "").split(";"):
				if ensemblId.startswith("ENSRNOG"):
					symbolByEnsembl[ensemblId] = row["SYMBOL"]
	return symbolByEnsembl


SYMBOLS = None
ENSEMBL = None
ENSEMBL_PATTERN = re.compile(r"ENS(?:MUS|RNO)?G[0-9]{11}")


NAMES = None


def LoadGeneNames():
	# Current symbol (upper case) -> every name the gene had, so an Ensembl ID can confirm "Sept6" or "Marc1"
	# even when the old name is ambiguous in the symbol maps; human, mouse and rat share keys (SEPTIN2 / Septin2)
	namesByCurrent = {}

	def Add(currentSymbol, names):
		nameSet = namesByCurrent.setdefault(currentSymbol.upper(), {currentSymbol.upper()})
		nameSet.update(name.strip().upper() for name in names if name.strip())

	with open(REF_DIR / "hgnc_complete_set.txt", encoding="utf-8") as handle:
		for row in csv.DictReader(handle, delimiter="\t"):
			if row["status"] == "Approved":
				Add(row["symbol"], row["prev_symbol"].strip('"').split("|") + row["alias_symbol"].strip('"').split("|"))
	with open(REF_DIR / "MRK_List1.rpt", encoding="utf-8") as handle:
		for row in csv.DictReader(handle, delimiter="\t"):
			if row["Marker Type"] != "Gene":
				continue
			if row["Status"] == "O":
				Add(row["Marker Symbol"], (row["Marker Synonyms (pipe-separated)"] or "").split("|"))
			elif row["Current Marker Symbol (if withdrawn)"]:
				Add(row["Current Marker Symbol (if withdrawn)"], [row["Marker Symbol"]])
	with open(REF_DIR / "GENES_RAT.txt", encoding="utf-8") as handle:
		lines = (line for line in handle if not line.startswith("#"))
		for row in csv.DictReader(lines, delimiter="\t"):
			Add(row["SYMBOL"], (row.get("OLD_SYMBOL") or "").split(";"))
	return namesByCurrent


def InitWorker():
	global SYMBOLS, ENSEMBL, NAMES
	SYMBOLS = LoadSymbolSets()
	ENSEMBL = LoadEnsemblMap()
	NAMES = LoadGeneNames()


def AssignTier(guessOriginal, currentSymbols, rowEnsemblSymbols, namesByCurrent=None):
	# A: an Ensembl ID in the same row is the guessed gene (200/200 right in the pilot), counts without review
	# B: no ID in the row but the date maps back to a date-prone symbol, counts and gets a sampled review
	# C: no guess, or the ID in the row is another gene, goes to manual review
	if rowEnsemblSymbols:
		guessed = {symbol.upper() for symbol in currentSymbols.split("/") if symbol}
		rowSymbols = {symbol.upper() for symbol in rowEnsemblSymbols}
		if guessed & rowSymbols:
			return "A_ensembl_match"
		# Or one of the guessed old names is a former name of the gene the ID points to
		guessedOriginals = {symbol.upper() for symbol in guessOriginal.split("/") if symbol}
		for rowSymbol in rowSymbols:
			if guessedOriginals & (namesByCurrent or {}).get(rowSymbol, set()):
				return "A_ensembl_match"
		return "C_ensembl_mismatch"
	if guessOriginal:
		return "B_guess_no_id"
	return "C_no_guess"


def ExcelSerialFromDate(value):
	# 1900 date system, with Excel's fake 29 Feb 1900
	baseDate = datetime.datetime(1899, 12, 30)
	if isinstance(value, datetime.time):
		return None
	if not isinstance(value, datetime.datetime):
		value = datetime.datetime(value.year, value.month, value.day)
	return (value - baseDate).days


def GuessOriginalSymbol(month, day, symbolMap, year=None):
	# 2-Sep -> SEPT2, also tried in title case for mouse and rat
	# A date on day 1 can also be month-year: "Sept5" -> Sep 2005, so the year's last two digits are candidates too
	numbers = [day]
	if year is not None and day == 1 and 1990 <= year <= 2030:
		numbers.append(year % 100)
	originals = []
	currents = []
	for number in numbers:
		for prefix in MONTH_PREFIXES.get(month, []):
			for candidate in (prefix + str(number), prefix.capitalize() + str(number), prefix + "0" + str(number)):
				if candidate in symbolMap and candidate not in originals:
					originals.append(candidate)
					currents.append(symbolMap[candidate])
	return "/".join(originals), "/".join(currents)


def ParseMonthDay(text):
	# The whole cell must be the date: "4-octyl itaconate" is not 4-Oct (it was a false hit in the mouse sample)
	match = re.match(r"^([0-9]{1,2})[-/ ]([A-Za-zéä]{3,4})\.?([-/ ][0-9]{2,4})?$", text)
	if not match:
		match2 = re.match(r"^([A-Za-zéä]{3,4})\.?[-/ ]([0-9]{1,2})$", text)
		if not match2:
			return None
		dayText, monthText = match2.group(2), match2.group(1)
	else:
		dayText, monthText = match.group(1), match.group(2)
	monthIndex = {"jan": 1, "feb": 2, "fev": 2, "fév": 2, "mar": 3, "mär": 3, "mrz": 3, "maa": 3, "apr": 4, "abr": 4,
		"avr": 4, "may": 5, "mai": 5, "mei": 5, "jun": 6, "giu": 6, "jul": 7, "jui": 7, "lug": 7, "aug": 8, "ago": 8,
		"sep": 9, "set": 9, "oct": 10, "out": 10, "okt": 10, "ott": 10, "nov": 11, "dec": 12, "dez": 12, "dic": 12,
		"déc": 12, "ene": 1}.get(monthText[:3].lower())
	if not monthIndex:
		return None
	return monthIndex, int(dayText)


def ClassifyCell(value, symbolMap):
	# Returns (kind, rule, rendered text, guessed original, current symbol) or None
	if isinstance(value, (datetime.datetime, datetime.date)) and not isinstance(value, datetime.time):
		serial = ExcelSerialFromDate(value)
		original, current = GuessOriginalSymbol(value.month, value.day, symbolMap, value.year)
		return ("date_cell", "date_cell", value.strftime("%Y-%m-%d") + " (serial " + str(serial) + ")", original, current)
	if isinstance(value, bool) or value is None:
		return None
	if isinstance(value, (int, float)):
		if isinstance(value, float) and abs(value) >= 1e10:
			return ("number", "sci_E", repr(value), "", "")
		if float(value).is_integer() and 10000 <= abs(value) <= 99999:
			serialDate = datetime.datetime(1899, 12, 30) + datetime.timedelta(days=int(value))
			original, current = GuessOriginalSymbol(serialDate.month, serialDate.day, symbolMap, serialDate.year)
			return ("number", "five_digit", str(int(value)), original, current)
		return None
	text = str(value).strip()
	if not text or len(text) > 40:
		return None
	for ruleName, pattern in BASE_RULES + EXTRA_RULES:
		if pattern.search(text):
			original, current = "", ""
			monthDay = ParseMonthDay(text)
			if monthDay:
				original, current = GuessOriginalSymbol(monthDay[0], monthDay[1], symbolMap)
			elif ruleName == "five_digit":
				serialDate = datetime.datetime(1899, 12, 30) + datetime.timedelta(days=int(text))
				original, current = GuessOriginalSymbol(serialDate.month, serialDate.day, symbolMap, serialDate.year)
			isBase = ruleName in [name for name, unused in BASE_RULES]
			return ("text", ruleName if isBase else "extra:" + ruleName, text, original, current)
	return None


def ColumnLetter(index):
	letters = ""
	index += 1
	while index:
		index, remainder = divmod(index - 1, 26)
		letters = chr(65 + remainder) + letters
	return letters


def ReadSheets(path):
	# Yields (sheet name, list of columns, each a list of raw values)
	# Go by the magic bytes, extensions lie (.xls files that are xlsx or plain TSV)
	with open(path, "rb") as handle:
		magic = handle.read(4)
	if magic == b"PK\x03\x04":
		import openpyxl
		# A file handle skips openpyxl's extension check
		# The handle is closed on the way out, Windows refuses to delete a file that is still open
		with open(path, "rb") as fileHandle:
			workbook = openpyxl.load_workbook(fileHandle, read_only=True, data_only=True)
			try:
				for sheet in workbook.worksheets:
					rows = []
					cellCount = 0
					for row in sheet.iter_rows(values_only=True):
						rows.append(row)
						cellCount += max(len(row), 1)
						if cellCount > MAX_CELLS_PER_SHEET:
							break
					# Transpose, padding short rows
					width = max((len(row) for row in rows), default=0)
					columns = [[row[columnIndex] if columnIndex < len(row) else None for row in rows] for columnIndex in range(width)]
					yield sheet.title, columns
			finally:
				workbook.close()
	elif magic != b"\xd0\xcf\x11\xe0":
		# Plain text saved with an Excel extension, read as TSV
		with open(path, encoding="utf-8", errors="replace") as handle:
			rows = [line.rstrip("\r\n").split("\t") for line in handle]
		width = max((len(row) for row in rows), default=0)
		yield "text", [[row[columnIndex] if columnIndex < len(row) else None for row in rows] for columnIndex in range(width)]
	else:
		import xlrd
		# Read from bytes, so xlrd keeps no handle on the file
		workbook = xlrd.open_workbook(file_contents=path.read_bytes())
		for sheetIndex in range(workbook.nsheets):
			sheet = workbook.sheet_by_index(sheetIndex)
			columns = []
			for columnIndex in range(sheet.ncols):
				values = []
				for rowIndex in range(min(sheet.nrows, MAX_CELLS_PER_SHEET // max(sheet.ncols, 1))):
					cell = sheet.cell(rowIndex, columnIndex)
					if cell.ctype == xlrd.XL_CELL_DATE:
						try:
							values.append(xlrd.xldate.xldate_as_datetime(cell.value, workbook.datemode))
						except Exception:
							values.append(cell.value)
					elif cell.ctype == xlrd.XL_CELL_EMPTY:
						values.append(None)
					else:
						values.append(cell.value)
				columns.append(values)
			yield sheet.name, columns


def CurrentSymbolsIn(cells, symbolMap):
	# Only current official symbols count for detection, old names are too noisy (RGD lists "TT", "PD", "K5")
	return {value.strip() for value in cells if isinstance(value, str) and symbolMap.get(value.strip()) == value.strip()}


def PickSpecies(topCells, wholeColumn):
	# Distinct symbols, so a column repeating "PC1" or "DBP" does not count as a gene list
	counts = {species: len(CurrentSymbolsIn(topCells, symbolMap)) for species, symbolMap in SYMBOLS.items()}
	bestSpecies = max(counts, key=counts.get)
	if counts[bestSpecies] < MIN_SYMBOLS_IN_TOP_20:
		return None, counts
	# Mouse and rat share most symbols, so a near tie in the top 20 (3 or less) is broken on the whole column
	if bestSpecies in ("mouse", "rat") and abs(counts["mouse"] - counts["rat"]) <= 3:
		wholeCounts = {species: len(CurrentSymbolsIn(wholeColumn, SYMBOLS[species])) for species in ("mouse", "rat")}
		if wholeCounts["mouse"] == wholeCounts["rat"]:
			return "mouse_or_rat", counts
		bestSpecies = max(wholeCounts, key=wholeCounts.get)
	return bestSpecies, counts


def ScreenFile(path, rootDir=None):
	relativePath = path.relative_to(rootDir or SUPP_DIR).as_posix()
	pmcid = relativePath.split(".")[0]
	geneColumns = []
	flags = []
	try:
		for sheetName, columns in ReadSheets(path):
			for columnIndex, column in enumerate(columns):
				topCells = column[:20]
				species, counts = PickSpecies(topCells, column)
				if not species:
					continue
				symbolMap = SYMBOLS["mouse"] if species == "mouse_or_rat" else SYMBOLS[species]
				nonEmpty = sum(1 for value in column if value not in (None, ""))
				symbolCells = sum(1 for value in column if isinstance(value, str) and value.strip() in symbolMap)
				# Did the list use the 2020 names (SEPTIN2, MARCHF1) or keep the old ones as text (SEPT2, MARCH1)?
				newNameCells = sum(1 for value in column if isinstance(value, str) and NEW_NAME_PATTERN.match(value.strip()))
				oldNameCells = sum(1 for value in column if isinstance(value, str) and OLD_NAME_PATTERN.match(value.strip()))
				columnFlags = 0
				# Skip the header row, the gene rule runs on the data
				for rowIndex, value in enumerate(column[1:], start=2):
					result = ClassifyCell(value, symbolMap)
					if not result:
						continue
					columnFlags += 1
					kind, rule, rendered, original, current = result
					# Neighbour cells help the manual check
					leftValue = columns[columnIndex - 1][rowIndex - 1] if columnIndex > 0 else ""
					rightValue = columns[columnIndex + 1][rowIndex - 1] if columnIndex + 1 < len(columns) else ""
					# Any Ensembl gene ID in the same row, version suffix dropped by the pattern
					rowText = " ".join(str(otherColumn[rowIndex - 1]) for otherColumn in columns
						if rowIndex - 1 < len(otherColumn) and isinstance(otherColumn[rowIndex - 1], str))
					rowEnsemblIds = sorted(set(ENSEMBL_PATTERN.findall(rowText)))
					rowEnsemblSymbols = {ENSEMBL[ensemblId] for ensemblId in rowEnsemblIds if ensemblId in ENSEMBL}
					flags.append({"pmcid": pmcid, "file": relativePath, "sheet": sheetName,
						"cell": ColumnLetter(columnIndex) + str(rowIndex), "species": species, "kind": kind,
						"rule": rule, "value": rendered, "guess_original": original, "current_symbol": current,
						"row_ensembl": "/".join(rowEnsemblIds[:4]), "row_ensembl_symbols": "/".join(sorted(rowEnsemblSymbols)),
						"tier": AssignTier(original, current, rowEnsemblSymbols, NAMES),
						"left": str(leftValue)[:40], "right": str(rightValue)[:40]})
				geneColumns.append({"pmcid": pmcid, "file": relativePath, "sheet": sheetName,
					"column": ColumnLetter(columnIndex), "species": species, "top20_human": counts["human"],
					"top20_mouse": counts["mouse"], "top20_rat": counts["rat"], "non_empty": nonEmpty,
					"symbol_cells": symbolCells, "flagged_cells": columnFlags,
					"new_name_cells": newNameCells, "old_name_text_cells": oldNameCells})
		status = "ok"
	except Exception as error:
		status = "error: " + type(error).__name__ + " " + str(error)[:100]
	return {"pmcid": pmcid, "file": relativePath, "status": status}, geneColumns, flags


def SelfCheck():
	global SYMBOLS
	SYMBOLS = {"human": {"SEPT2": "SEPTIN2", "MARCHF1": "MARCHF1", "MARCH1": "MARCHF1", "TP53": "TP53"},
		"mouse": {"Sept2": "Septin2"}, "rat": {}}
	assert ClassifyCell("2-Sep", SYMBOLS["human"])[3] == "SEPT2"
	assert ClassifyCell("Mar-01", SYMBOLS["human"])[1] == "extra:Mon-d"
	assert ClassifyCell(datetime.datetime(2016, 3, 1), SYMBOLS["human"])[4] == "MARCHF1"
	assert ClassifyCell(2.310009e19, SYMBOLS["human"])[1] == "sci_E"
	assert ClassifyCell("2-set", SYMBOLS["human"])[1] == "extra:d-Mon_locale"
	assert ClassifyCell("TP53", SYMBOLS["human"]) is None
	assert ClassifyCell(0.05, SYMBOLS["human"]) is None
	assert ClassifyCell("2-Sep", SYMBOLS["mouse"])[3] == "Sept2"
	SYMBOLS["human"]["SEPT5"] = "SEPTIN5"
	assert ClassifyCell(datetime.datetime(2005, 9, 1), SYMBOLS["human"])[3] == "SEPT5"
	# Tiers: ID agrees, ID disagrees, no ID with guess, nothing
	assert AssignTier("SEPT9", "SEPTIN9", {"SEPTIN9"}) == "A_ensembl_match"
	assert AssignTier("MARCH1/MARCH11", "MARCHF1/MARCHF11", {"MARCHF11"}) == "A_ensembl_match"
	assert AssignTier("Sept2", "Septin2", {"SEPTIN2"}) == "A_ensembl_match"
	assert AssignTier("MARCH1", "MARCHF1", {"TP53"}) == "C_ensembl_mismatch"
	assert AssignTier("SEPT2", "SEPTIN2", set()) == "B_guess_no_id"
	assert AssignTier("", "", set()) == "C_no_guess"
	assert AssignTier("Sept1/Sep1/Sept6", "Septin1/Selenof/Foo", {"Septin6"}, {"SEPTIN6": {"SEPTIN6", "SEPT6"}}) == "A_ensembl_match"
	# make.unique names from R are not dates
	assert ClassifyCell("Jun.1", SYMBOLS["mouse"]) is None
	assert ClassifyCell("Set.1", SYMBOLS["mouse"]) is None
	assert ClassifyCell("Marcks.1", SYMBOLS["mouse"]) is None
	assert ClassifyCell("Sep-02", SYMBOLS["human"])[1] == "extra:Mon-d"
	SYMBOLS["human"]["OCT4"] = "POU5F1"
	assert ClassifyCell("4-octyl itaconate", SYMBOLS["human"])[3] == ""
	assert ClassifyCell("4-Oct", SYMBOLS["human"])[3] == "OCT4"
	assert ENSEMBL_PATTERN.findall("ENSG00000184640.13, ENSMUSG00000001833 x") == ["ENSG00000184640", "ENSMUSG00000001833"]
	# Real reference files: the IDs the pilot leaned on must map to the right genes
	if (REF_DIR / "MRK_ENSEMBL.rpt").exists():
		ensemblMap = LoadEnsemblMap()
		assert ensemblMap["ENSG00000184640"] == "SEPTIN9"
		assert ensemblMap["ENSG00000183654"] == "MARCHF11"
		assert ensemblMap["ENSG00000117791"] == "MTARC2"
		assert ensemblMap["ENSMUSG00000026276"] == "Septin2"
		assert ensemblMap["ENSRNOG00000017952"] == "Septin2"
		geneNames = LoadGeneNames()
		assert "SEPT6" in geneNames["SEPTIN6"] and "MARC1" in geneNames["MTARC1"]
	print("self check ok")


def ParseArguments():
	# --supp DIR --out DIR --workers N --low-priority
	arguments = {"supp": SUPP_DIR, "out": DATA_DIR, "workers": 4, "low": False}
	position = 1
	while position < len(sys.argv):
		flag = sys.argv[position]
		if flag == "--low-priority":
			arguments["low"] = True
			position += 1
			continue
		value = sys.argv[position + 1]
		if flag == "--supp":
			arguments["supp"] = Path(value)
		elif flag == "--out":
			arguments["out"] = Path(value)
		elif flag == "--workers":
			arguments["workers"] = int(value)
		position += 2
	return arguments


def InitLowPriorityWorker():
	import lowprio
	lowprio.SetLowPriority()
	InitWorker()


def Main():
	if len(sys.argv) > 1 and sys.argv[1] == "--check":
		SelfCheck()
		return 0
	arguments = ParseArguments()
	supportDir = arguments["supp"]
	outDir = arguments["out"]
	outDir.mkdir(parents=True, exist_ok=True)
	if arguments["low"]:
		import lowprio
		lowprio.SetLowPriority()
	files = sorted(path for path in supportDir.rglob("*") if path.is_file() and re.search(r"\.xlsx?$", path.name, re.I))
	fileRows, columnRows, flagRows = [], [], []
	initializer = InitLowPriorityWorker if arguments["low"] else InitWorker
	with ProcessPoolExecutor(max_workers=arguments["workers"], initializer=initializer) as pool:
		for fileRow, geneColumns, flags in pool.map(functools.partial(ScreenFile, rootDir=supportDir), files, chunksize=2):
			fileRows.append(fileRow)
			columnRows.extend(geneColumns)
			flagRows.extend(flags)
			print(fileRow["status"][:60], fileRow["file"], len(geneColumns), len(flags), flush=True)
	for name, rows, fields in (("file_status.csv", fileRows, ["pmcid", "file", "status"]),
			("gene_columns.csv", columnRows, None), ("flags.csv", flagRows, None)):
		with open(outDir / name, "w", newline="", encoding="utf-8") as handle:
			fieldNames = fields or (list(rows[0]) if rows else ["empty"])
			writer = csv.DictWriter(handle, fieldnames=fieldNames)
			writer.writeheader()
			writer.writerows(rows)


if __name__ == "__main__":
	sys.exit(Main())
