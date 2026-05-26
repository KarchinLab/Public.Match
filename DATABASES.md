# Database Reference

This document describes the reference databases included in `Databases/` and how they are used by Public.Match.

> **Raw file size ≠ entries loaded for matching.**
> Each parser applies species filtering, chain-type filtering, valid amino-acid checks, and deduplication before the sequences are used for CDR3 matching. See the [Database Summary](#database-summary) at the bottom for both raw and loaded counts.

---

## IEDB — Immune Epitope Database

**File:** `Databases/IEDB/iedb.xlsx`
**Source:** https://www.iedb.org (receptor export)
**Raw rows:** ~217,000
**Loaded for matching:** ~164,240
**Format:** Excel (.xlsx), manually exported from the IEDB receptor search

### Filtering applied
- Extract CDR3β from whichever chain (1 or 2) has `Type == beta`; no pre-filter on receptor type (gammadelta/construct rows carry no beta chains and fall out naturally)
- Prefer `CDR3 Curated` over `CDR3 Calculated` where available
- Validate CDR3β as a canonical amino acid string
- Deduplicate on `(cdr3b, epitope)`
- Paired-mode CDR3α filtering is handled by `load_databases`, not the parser

### Key columns

| Column | Description |
|---|---|
| `Chain 1 - CDR3 Curated` / `Chain 1 - CDR3 Calculated` | CDR3 sequence for chain 1 (TRA or TRB) |
| `Chain 1 - Type` | Chain type (alpha / beta) |
| `Chain 2 - CDR3 Curated` / `Chain 2 - CDR3 Calculated` | CDR3 sequence for chain 2 |
| `Epitope - Name` | Epitope peptide sequence |
| `Epitope - Source Organism` | Pathogen or tissue of origin |
| `Assay - MHC Allele Names` | HLA restriction |
| `Receptor - Type` | alphabeta / gammadelta / etc. |

---

## VDJdb

**File:** `Databases/VDJdb/vdjdb_extracted/VDJdb_05302026.tsv`
**Fallback:** `Databases/VDJdb/vdjdb.slim.txt`
**Source:** https://vdjdb.cdr3.net (full export, May 2026)
**Raw rows:** ~209,000
**Loaded for matching:** ~95,256
**Format:** TSV

### Filtering applied
- Keep only `Species == HomoSapiens`
- Keep `Score >= 0` (all entries, including unverified; see note below)
- Pair TRB and TRA rows via `complex.id` where both chains are present
- Validate CDR3 sequences as canonical amino acid strings
- Deduplicate on `(cdr3b, cdr3a, epitope)`

### VDJdb confidence score

| Score | Meaning |
|---|---|
| 0 | No supporting publication; sequence is real but unverified |
| 1 | Single publication |
| 2 | Multiple independent publications |
| 3 | Structural or high-confidence functional validation |

Public.Match uses `Score >= 0` to maximise coverage. To restrict to publication-backed entries only, set `min_score=1` in `vdjdb.load()`.

### Key columns

| Column | Description |
|---|---|
| `CDR3` | CDR3 amino acid sequence |
| `Gene` | TRA or TRB |
| `V` / `J` | V and J gene calls |
| `Species` | HomoSapiens / MusMusculus |
| `MHC A` / `MHC B` | HLA alleles |
| `MHC class` | MHCI or MHCII |
| `Epitope` | Epitope peptide |
| `Epitope gene` | Antigen gene name |
| `Epitope species` | Pathogen |
| `Score` | VDJdb confidence score (0–3) |

---

## McPAS-TCR

**File:** `Databases/McPAS/McPAS-TCR.csv`
**Source:** http://friedmanlab.weizmann.ac.il/McPAS-TCR
**Raw rows:** ~40,700
**Loaded for matching:** ~29,649
**Format:** CSV (latin-1 encoded)

### Filtering applied
- Keep only `Species == Human`
- Require non-null `CDR3.beta.aa`
- Validate CDR3β as canonical amino acids
- Deduplicate on `(cdr3b, epitope)`

### Key columns

| Column | Description |
|---|---|
| `CDR3.beta.aa` | CDR3β amino acid sequence |
| `CDR3.alpha.aa` | CDR3α amino acid sequence |
| `TRBV` / `TRBJ` | V and J gene calls |
| `Epitope.peptide` | Epitope peptide |
| `MHC` | HLA restriction |
| `Pathology` | Disease/pathology category |
| `Category` | Infectious disease / Cancer / Autoimmune |
| `Species` | Human / Mouse |

### Notes
- Many entries lack an epitope peptide; use `Pathology` for disease-level annotation.

---

## 10x Genomics Dcode (Human Donors)

**Files:** `Databases/10xDcode/vdj_v1_hs_aggregated_donor{1–4}_binarized_matrix.csv`
**Long format:** `Databases/10xDcode/10xDcode_long.csv`
**Source:** 10x Genomics public dataset (4 healthy donors)
**Raw cells:** ~189,500 across 4 donors (~85,500 confirmed binder cell–epitope pairs)
**Loaded for matching:** ~18,561
**Format:** CSV

### Filtering applied
- Keep only cells where `_binder == True` for a given dextramer (~85,500 binder cell–epitope pairs)
- Cells with no CDR3β detected (~5,100) are excluded
- Deduplicate on `(cdr3b, epitope, HLA, donor)` — collapses clonally expanded cells (same CDR3β appearing across many cells of the same donor); ~85,500 pairs → ~18,800 unique clone–epitope combinations
- Parser deduplicates further on `(cdr3b, epitope)` — collapses the same TCR seen across multiple donors → ~18,561

> **Why ~190k cells → ~18k entries?** The raw files are one row per *cell*, not per *clone*. Clonally expanded T cells share the same CDR3β; the most abundant clone in the dataset appears in ~6,700 cells across one donor alone. After collapsing to unique CDR3β–epitope pairs the redundancy resolves to ~18,561 distinct entries.

### Epitope column format

```
A0201_GILGFVFTL_Flu-MP_Influenza
│     │         │      └─ Pathogen
│     │         └─ Antigen name
│     └─ Epitope peptide
└─ HLA allele
```

### Key columns (long format)

| Column | Description |
|---|---|
| `cdr3b` | CDR3β amino acid sequence |
| `cdr3a` | CDR3α amino acid sequence (first TRA chain) |
| `epitope` | Epitope peptide |
| `HLA` | HLA allele (e.g. `A0201`) |
| `antigen` | Antigen name (e.g. `Flu-MP`) |
| `pathogen` | Pathogen (e.g. `Influenza`, `CMV`, `Cancer`) |

---

## MixTCRpred

**File:** `Databases/MixTCRpred/full_training_set_146pmhc.csv`
**Source:** MixTCRpred training data (Heidelberg group)
**Raw rows:** ~17,700
**Loaded for matching:** ~6,875
**Format:** CSV

### Filtering applied
- Keep only `species == HomoSapiens`
- Require non-null `cdr3_TRB`
- Validate CDR3β as canonical amino acids
- Deduplicate on `(cdr3b, epitope)`

### Key columns

| Column | Description |
|---|---|
| `cdr3_TRB` | CDR3β amino acid sequence |
| `cdr3_TRA` | CDR3α amino acid sequence |
| `epitope` | Epitope peptide |
| `MHC` | HLA allele |
| `species` | HomoSapiens / MusMusculus |

---

## BATCAVE

**Files:** `Databases/BATCAVE/TCR_pMHCI_mutational_scan.csv`, `TCR_pMHCII_mutational_scan.csv`
**Source:** BATCAVE mutational scan dataset
**Raw rows:** ~24,875 (MHC-I) + ~5,730 (MHC-II)
**Loaded for matching:** ~34
**Format:** CSV

### Filtering applied
- Keep only `tcr_source_organism == human`
- Keep only rows where `peptide == index_peptide` (native epitope; excludes thousands of mutagenesis scan variants)
- Keep only `peptide_activity > 0` (any positive activation)
- Deduplicate on `(cdr3b, epitope)`

Most of the raw rows are mutagenesis scan variants (point mutations of the epitope tested against each TCR). After collapsing to native-epitope rows only, ~60 human TCRs with positive activation survive.

> **Why not a higher activity threshold?** BATCAVE mixes 13 assay types on incompatible scales: normalized assays (CD137, IFNg, T-Scan, NFAT-GFP) report activity on a 0–1 scale where 1.0 = full native-peptide response, while absolute assays (ELISA, NFAT luminescence) have ranges up to 34,000+. A single absolute threshold (e.g. ≥ 20) would silently drop all entries from normalized assays even when they show full binding. Since `peptide == index_peptide` already confirms these are real binders, `activity > 0` is the correct criterion.

### Key columns

| Column | Description |
|---|---|
| `cdr3b` | CDR3β amino acid sequence |
| `cdr3a` | CDR3α amino acid sequence |
| `index_peptide` | Native (wild-type) epitope peptide |
| `peptide` | Tested peptide (may be a mutant) |
| `peptide_activity` | Functional activation score |
| `mhc` | HLA allele |
| `tcr_source_organism` | `human` or `mouse` |

---

## NeoTCR

**File:** `Databases/NeoTCR/NeoTCR data-20221220.xlsx`
**Source:** NeoTCR neoantigen-reactive TCR dataset
**Raw rows:** ~1,000
**Loaded for matching:** ~916
**Format:** Excel (.xlsx)

### Filtering applied
- Require non-null `TRB_CDR3`
- Validate CDR3β as canonical amino acids
- Deduplicate on `(cdr3b, epitope)`

### Key columns

| Column | Description |
|---|---|
| `TRB_CDR3` | CDR3β amino acid sequence |
| `TRA_CDR3` | CDR3α amino acid sequence |
| `Neoepitope` | Neoantigen peptide |
| `Antigen` | Antigen gene |
| `Tumor` | Tumor type |
| `HLA Allele` | HLA restriction |

---

## CEDAR

**File:** `Databases/CEDAR/cedar.xlsx`
**Source:** CEDAR immunology database
**Raw rows:** ~76,200
**Loaded for matching:** ~41,266
**Format:** Excel (.xlsx)

### Filtering applied
- Parse paired, beta-only, and alpha-only rows by chain type columns (no receptor-type pre-filter; 5 gammadelta and 4 construct rows carry no beta chains and fall out naturally)
- Prefer `CDR3 Curated` over `CDR3 Calculated`
- Validate CDR3 sequences as canonical amino acids
- Deduplicate on `(cdr3b, cdr3a, epitope)`
- Note: alpha-only rows are excluded when running in beta or paired mode (no CDR3β present)

### Key columns

| Column | Description |
|---|---|
| `Chain 1 - CDR3 Curated` / `Chain 1 - CDR3 Calculated` | CDR3 for chain 1 |
| `Chain 1 - Type` | alpha / beta |
| `Chain 2 - CDR3 Curated` / `Chain 2 - CDR3 Calculated` | CDR3 for chain 2 |
| `Chain 2 - Type` | alpha / beta |
| `Epitope - Name` | Epitope peptide |
| `Epitope - Source Molecule` | Antigen protein |
| `Epitope - Source Organism` | Pathogen or tissue |
| `Assay - MHC Allele Names` | HLA restriction |

---

## Database Summary

Counts shown are for **beta-chain mode** (default), after all filtering and deduplication.

| Database | Raw file rows | Loaded for matching | Key filters |
|---|---|---|---|
| IEDB | ~217,000 | ~164,240 | beta chain rows only; dedup on (cdr3b, epitope) |
| VDJdb | ~209,000 | ~95,256 | HomoSapiens; score ≥ 0; dedup on (cdr3b, cdr3a, epitope) |
| McPAS | ~40,700 | ~29,649 | Human only; dedup on (cdr3b, epitope) |
| 10x Dcode | ~190k cells / ~85k binder pairs | ~18,561 | Confirmed binders; clonal dedup on (cdr3b, epitope) |
| MixTCRpred | ~17,700 | ~6,875 | HomoSapiens; dedup on (cdr3b, epitope) |
| BATCAVE | ~30,600 | ~60 | Human; native peptide only; activity > 0 |
| NeoTCR | ~1,000 | ~916 | dedup on (cdr3b, epitope) |
| CEDAR | ~76,200 | ~41,266 | beta chain rows only; dedup on (cdr3b, cdr3a, epitope) |
| **Total** | | **~356,454** | |
