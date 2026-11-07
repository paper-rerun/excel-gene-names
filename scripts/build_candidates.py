# Fetch every PMCID for a Europe PMC query and save them shuffled with a fixed seed
# Shuffled on purpose: if the run stops after N nights, what is done is still a random sample
# Usage: build_candidates.py <name> "<query>" [seed]
import csv
import datetime
import json
import random
import sys
import time
from pathlib import Path

import requests

BASE_DIR = Path(__file__).resolve().parent.parent
SEARCH_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


def FetchAll(query):
	session = requests.Session()
	session.headers["User-Agent"] = "PaperRerun-pilot/0.1 (volunteer research audit)"
	hitsByPmcid = {}
	cursor = "*"
	hitCount = None
	while True:
		params = {"query": query, "format": "json", "pageSize": 1000, "resultType": "lite", "cursorMark": cursor}
		for attempt in range(5):
			try:
				response = session.get(SEARCH_URL, params=params, timeout=120)
				response.raise_for_status()
				break
			except requests.RequestException:
				time.sleep(10 * (attempt + 1))
		else:
			raise RuntimeError("Europe PMC kept failing at cursor " + cursor)
		page = response.json()
		hitCount = page["hitCount"]
		results = page["resultList"]["result"]
		for hit in results:
			if hit.get("pmcid"):
				hitsByPmcid[hit["pmcid"]] = hit
		nextCursor = page.get("nextCursorMark")
		print(len(hitsByPmcid), "/", hitCount, flush=True)
		if not results or not nextCursor or nextCursor == cursor:
			break
		cursor = nextCursor
		time.sleep(0.3)
	return hitCount, list(hitsByPmcid.values())


def Main():
	name = sys.argv[1]
	query = sys.argv[2]
	seed = int(sys.argv[3]) if len(sys.argv) > 3 else 20261104
	runDate = datetime.datetime.now(datetime.timezone.utc).isoformat()
	hitCount, hits = FetchAll(query)
	# Same research-article filter as the pilot
	keptHits = []
	for hit in hits:
		publicationType = hit.get("pubType", "").lower()
		title = hit.get("title", "").lower()
		if "review" in publicationType or "erratum" in publicationType or "correction" in publicationType:
			continue
		if title.startswith(("correction", "erratum", "retraction", "publisher correction", "author correction")):
			continue
		keptHits.append(hit)
	keptHits.sort(key=lambda hit: hit["pmcid"])
	random.Random(seed).shuffle(keptHits)
	outPath = BASE_DIR / "data" / ("candidates_" + name + ".csv")
	with open(outPath, "w", newline="", encoding="utf-8") as handle:
		writer = csv.writer(handle)
		writer.writerow(["order", "pmcid", "journal", "doi", "pub_year"])
		for order, hit in enumerate(keptHits):
			writer.writerow([order, hit["pmcid"], hit.get("journalTitle", ""), hit.get("doi", ""), hit.get("pubYear", "")])
	logEntry = {"name": name, "query": query, "run_utc": runDate, "hit_count": hitCount, "with_pmcid": len(hits),
		"research_articles": len(keptHits), "seed": seed, "file": outPath.name}
	with open(BASE_DIR / "logs" / "candidate_queries.jsonl", "a", encoding="utf-8") as handle:
		handle.write(json.dumps(logEntry) + "\n")
	print(json.dumps(logEntry))


if __name__ == "__main__":
	sys.exit(Main())
