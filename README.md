# Public.Match

**Hackathon Project** | TCR Repertoire Analysis

---

## The Problem

T-cell receptor (TCR) repertoires from cancer patients contain valuable signals — but identifying which TCRs are "public" (shared across individuals and described in curated databases) requires tools that are tightly coupled to a single database. **TCRMatch**, for example, only queries IEDB. Researchers working with VDJdb, McPAS-TCR, or custom epitope databases have no unified solution.

## Our Idea

**Public.Match** is a generalized TCR public-sequence matching tool that:

- Accepts patient TCR repertoires (CDR3α/β sequences) as input
- Searches across multiple curated databases — **IEDB**, **VDJdb**, **McPAS-TCR**, **10x Genomics pMHC**, and others
- Returns matched public TCRs with epitope annotations, HLA restrictions, and match scores
- Is built with **Claude Code**, using AI-assisted development to rapidly prototype and extend the tool during the hackathon

Think of it as TCRMatch — but database-agnostic.

## Why It Matters

| Today | With Public.Match |
|---|---|
| Run TCRMatch → IEDB only | Single query → IEDB + VDJdb + McPAS + 10x |
| Manual format conversion per DB | Unified input/output schema |
| No cross-database deduplication | Merged, ranked hits across all sources |

Identifying public TCRs that recognize known epitopes helps distinguish **antigen-specific** from **bystander** T cells — a key step in spatial immunology pipelines like our own [Soleil](https://github.com/Marcus-Mendes) engagement scoring framework.

## Approach

1. **Unify database schemas** — normalize all sources into a common `cdr3_alpha / cdr3_beta / epitope / mhc / v_gene / j_gene / source` format
2. **Extend or wrap TCRMatch** — reuse its edit-distance / GLIPH2-inspired scoring logic against any database
3. **Build a CLI** — `public-match --input repertoire.tsv --db iedb vdjdb mcpas 10x --score 0.97`
4. **Claude Code as co-developer** — use Claude Code to accelerate implementation, handle format parsing edge cases, and generate test cases

## Databases

| Database | Entries | CDR3α | CDR3β | Epitope | HLA | Folder |
|---|---|---|---|---|---|---|
| [IEDB](https://www.iedb.org/) | 226,280 TCR records | ✓ | ✓ | ✓ | ✓ | `Databases/IEDB/` |
| [VDJdb](https://vdjdb.cdr3.net/) | 145,408 chain records | ✓ | ✓ | ✓ | ✓ | `Databases/VDJdb/` |
| [McPAS-TCR](https://friedmanlab.weizmann.ac.il/McPAS-TCR/) | 40,779 | ✓ | ✓ | ✓ | ✓ | `Databases/McPAS/` |
| [10x Genomics pMHC](https://www.10xgenomics.com/) | 189,515 cells / 4 donors | ✓ | ✓ | ✓ (55 pMHC) | ✓ | `Databases/10xDcode/` |
| [MixTCRpred](https://github.com/GfellerLab/MixTCRpred) | 17,715 αβ pairs | ✓ | ✓ | ✓ (146 pMHC) | ✓ | `Databases/MixTCRpred/` |
| [BATCAVE](https://github.com/meyer-lab-cshl/BATMAN-paper) | 24,875 TCR–peptide measurements | ✓ | ✓ | ✓ (mutational scan) | ✓ | `Databases/BATCAVE/` |
| [OTS](https://opig.stats.ox.ac.uk/webapps/ots/) | 1.63M non-redundant paired αβ | ✓ | ✓ | — (publicness only) | — | `Databases/OTS/` (manual download) |

### Unified schema

All databases are mapped to a common record format:

```
cdr3_alpha    CDR3 amino acid sequence of the alpha chain
cdr3_beta     CDR3 amino acid sequence of the beta chain
epitope       Epitope peptide sequence
mhc           MHC/HLA restriction (e.g. HLA-A*02:01)
v_alpha       TRAV gene
j_alpha       TRAJ gene
v_beta        TRBV gene
j_beta        TRBJ gene
source        Database of origin (IEDB / VDJdb / McPAS / 10x)
```

### Databases on the roadmap

| Database | Entries | Notes |
|---|---|---|
| [TCRdb 2.0](https://guolab.wchscu.cn/TCRdb2/) | ~700M sequences | Broad clinical coverage; no epitope labels |
| [STCRDab](http://opig.stats.ox.ac.uk/webapps/stcrdab/) | ~1,000 | 3D structural data from PDB |
| [PIRD](https://db.cngb.org/pird/) | large | Pan Immune Repertoire Database; China National GeneBank |
| [ePytope-TCR datasets](https://www.cell.com/cell-genomics/fulltext/S2666-979X(25)00202-2) | 21 datasets / 762 epitopes | 2025 benchmarking collection on Zenodo |

## Installation

```bash
git clone https://github.com/Marcus-Mendes/Public.Match
cd Public.Match
conda env create -f environment.yml
conda activate public-match
pip install -e .
```

## Usage

### Beta chain (default)

```bash
python -m public_match --input sequences.fasta
```

### Alpha chain

```bash
python -m public_match --input alpha_seqs.fasta --chain alpha
```

### Paired αβ — two FASTA files (matched by sequence name)

```bash
python -m public_match --input beta.fasta --input-alpha alpha.fasta --chain paired
```

### Paired αβ — single TSV/CSV with both columns

```bash
python -m public_match --input repertoire.tsv --chain paired
# explicit column names if auto-detection fails:
python -m public_match --input repertoire.tsv --chain paired \
  --seq-col cdr3_beta --alpha-col cdr3_alpha
```

### Select specific databases

```bash
python -m public_match --input sequences.fasta --db iedb vdjdb mcpas
```

Available databases: `iedb`, `vdjdb`, `mcpas`, `10x`, `mixtcrpred`, `batcave`, `ots`

### Matching methods and thresholds

```bash
# BLOSUM62-normalised score (default, 0–1; 0.97 = near-exact match)
python -m public_match --input sequences.fasta --method blosum --threshold 0.97

# Edit distance (integer; 1 = one substitution allowed)
python -m public_match --input sequences.fasta --method edit --threshold 1

# Exact match only
python -m public_match --input sequences.fasta --method exact
```

### Custom database

```bash
python -m public_match --input sequences.fasta --custom-db my_db.csv
# specify the CDR3β column if it differs from the defaults:
python -m public_match --input sequences.fasta --custom-db my_db.tsv \
  --custom-db-cdr3-col junction_aa
```

### Output

Results are written to `public_match_results.csv` by default. Use `--output` to change the path:

```bash
python -m public_match --input sequences.fasta --output results/my_run.csv
```

Output columns: `query_name`, `query_cdr3b` (and `query_cdr3a` for paired mode), `cdr3_alpha`, `cdr3_beta`, `epitope`, `mhc`, `v_alpha`, `j_alpha`, `v_beta`, `j_beta`, `source`, `score`.

### All options

```
python -m public_match --help

  --input/-i PATH          Input file: FASTA or tabular (TSV/CSV/AIRR)
  --output/-o PATH         Output CSV (default: public_match_results.csv)
  --db DB [DB ...]         Databases to search (default: all)
  --method {blosum,edit,exact}  Matching method (default: blosum)
  --threshold FLOAT        Score threshold (default: 0.97)
  --chain {beta,alpha,paired}   Chain mode (default: beta)
  --seq-col COL            CDR3β column in tabular input
  --alpha-col COL          CDR3α column in tabular input
  --name-col COL           ID column in tabular input
  --input-alpha PATH       CDR3α FASTA for paired mode
  --custom-db PATH [PATH ...]  Custom database file(s)
  --custom-db-cdr3-col COL     CDR3β column in custom DB
```

---

## Hackathon Deliverable

A working CLI prototype that takes a CDR3 repertoire file and returns matched public TCRs from all four databases, with a unified output format and match score.

---

*Built at Hackathon · 2026 · with Claude Code*
