# Gene names turned into dates by Excel, 2024-2026

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24039868.svg)](https://doi.org/10.5281/zenodo.24039868)

**21.1% of open-access genomics papers that ship gene lists in Excel files (1,978 of 9,367, 95% CI 20.3-22.0%) have at least one gene name converted into a date or a number.** The same kind of screen found 30.9% for 2014-2020 (Abeysooriya et al. 2021, PLOS Comput Biol 17:e1008984) and 19.6% for 2005-2015 (Ziemann et al. 2016, Genome Biol 17:177). The rate dropped after HGNC renamed the SEPT and MARCH genes in 2020, but it is still 2 papers in 10.

This repository has the aggregate results and the code. It does not list individual papers. Part of the Paper Rerun project (https://paper-rerun.github.io). Built by AI agents on a home computer between 4 and 7 November 2026 and reviewed by a human before publication.

## Results

| Gene lists | Papers with lists | Affected, estimated (95% CI) | 2014-2020 |
|---|---|---|---|
| All | 9,367 | **21.1% (20.3-22.0)** | 30.9% |
| Human | 8,645 | **20.0% (19.2-20.8)** | 30.5% |
| Mouse | 952 | **29.5% (26.7-32.5)** | 38.6% |
| Rat | 159 | **12.6% (8.3-18.6)** | 20.8% |

- By year: 2024 22.9%, 2025 21.0%, 2026 19.3%.
- In the 23 journals that both screens cover with at least 30 papers, the rate fell in 22. Pooled: 23.1% now against 32.0% in 2014-2020. eLife is the exception (37.2% against 36.6%).
- Highest among the 55 journals with at least 30 gene-list papers: EMBO J 45.7%, PNAS 40.2%, eLife 37.2%. Nature Communications, the largest, has 27.8% (270 of 972).
- Old symbols behind the errors (share of affected papers, date and text cells only): MARCH1-11 69.4%, SEPT1-14 60.4%, DEC1 24.1%, SEP15 19.3%. Excel set to other languages adds new ones: AGO1-4 read as August in Portuguese, Spanish or Italian (51 papers), MEI1/MEI4 in Dutch (14), and, found only in the manual review, MARS1/MARS2 in French (6), TAMM41 in Finnish (5) and GEN1 in Italian (2).

The 2020 renaming fixed the names, not the files: most errors come from annotation made before 2020 (GTEx, GENCODE, methylation arrays, eQTL tables).

## Tables

| File | What it is |
|---|---|
| `tables/rates_by_species.csv` | rates per species, automatic call and estimate, with Wilson 95% CIs |
| `tables/rates_by_year.csv` | rates per publication year |
| `tables/rates_by_journal.csv` | rates for the 55 journals with at least 30 gene-list papers |
| `tables/journals_vs_2014_2020.csv` | the 23 journals also in Table 5 of Abeysooriya et al. 2021, then and now |
| `tables/old_symbol_families.csv` | which old symbols cause the errors |
| `tables/precision_sample.csv` | hand check of 100 automatic calls, by species and tier |
| `tables/candidates_full_genom.csv` | the 81,125 articles screened (PMCID, journal, DOI, year), in the shuffled order of the run. It does not say which ones were affected |

## Method

- **Corpus:** every open-access PMC article from 2024-2026 with supplementary material and "genom*" in the title or abstract, the same cut as the 2021 study. Europe PMC query on 2026-11-04: `(TITLE:genom* OR ABSTRACT:genom*) AND OPEN_ACCESS:y AND HAS_SUPPL:y AND PUB_YEAR:[2024 TO 2026] AND IN_PMC:y`, 83,486 hits, 81,125 after removing reviews, errata and corrections.
- **Files:** only `.xls`/`.xlsx`, from the open PMC bucket on AWS (`pmc-oa-opendata`), read and deleted right away. 158.5 GB went through the machine. 30,925 articles had spreadsheets, 9,367 had gene lists.
- **Gene column:** at least 5 distinct current official symbols (HGNC, MGI, RGD) in the first 20 cells.
- **Flagged cell:** inside a gene column, a date-typed cell, a 5-digit date serial, a date written as text (`1-Mar`, `Sep-02`, `2-set`) or scientific notation, mapped back to the symbol it came from (2-Sep to SEPT2, ago/2002 to AGO2).
- **Automatic call:** a paper counts as affected with one cell confirmed by an Ensembl ID on the same row (tier A), a date that maps back to a date-prone symbol (tier B), or two or more such 5-digit serials. Everything else goes to manual review.
- **Precision:** 100 automatic calls sampled at random (seed 20261107), stratified by species and tier, re-downloaded and checked by hand: 100 of 100 confirmed (95% CI 96.3-100%). The 103 papers sent to manual review gave 37 real, 14 probable, 7 unclear and 45 false positives.
- **Estimate:** automatic calls times the precision of their stratum, plus the confirmed manual ones. Counting the 14 probable ones too gives 21.3%.

## What the screen can't see

- Horizontal gene lists, CSV, TSV, zipped files, tables in Word or PDF.
- Gene columns with fewer than 5 distinct symbols in the first 20 rows.
- Conversions that left no trace (the row was deleted or overwritten).
- Locales the detector does not know yet: Finnish, French and Italian only showed up in the manual review.
- Species come from the symbols, not the paper: cattle, goat, pig or chicken lists in upper case count as human (4 of 73 "human" papers in the precision sample).
- Only open-access PMC papers.

So these figures are a floor.

## Rerunning it

Python 3.12 with `openpyxl`, `xlrd` and `requests`. `scripts/build_candidates.py` rebuilds the article list from Europe PMC, `scripts/run_batch.py` runs the screen (resumable, with a time window and a `PAUSE` file, see `config.ini`), `scripts/analyze_db.py` and `scripts/full_report.py` produce the numbers. `scripts/detect.py --check` runs the detector's self-test. The HGNC, MGI and RGD symbol sets are downloaded from their official sites.

## About individual papers

The per-paper list (file, sheet, cell, value, old and current symbol) exists but is not published here. We would rather tell the authors first and give them time to correct their files. If you are an author or editor and want the entries for your papers or journal, open an issue in this repository.

## Licence

CC BY 4.0 (see `LICENSE.md`). Please cite this repository and the two earlier studies above.
