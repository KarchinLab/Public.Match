import pandas as pd
from pathlib import Path

# Full VDJdb export (May 2026) — larger, more up-to-date than the slim build
VDJDB_PATH = Path("Databases/VDJdb/vdjdb_extracted/VDJdb_05302026.tsv")

# Fallback to the slim build if the full export is absent
VDJDB_SLIM_PATH = Path("Databases/VDJdb/vdjdb.slim.txt")

_AA_PAT = r"^[ACDEFGHIKLMNPQRSTVWY]+$"


def _load_full(path: Path, min_score: int) -> pd.DataFrame:
    """Parse the VDJdb full export (columns: Gene, CDR3, Species, Score, Epitope, …)."""
    df = pd.read_csv(path, sep="\t", low_memory=False)

    df = df[df["Species"] == "HomoSapiens"].copy()
    df = df[df["Score"].fillna(0).astype(int) >= min_score]

    tra = df[df["Gene"] == "TRA"].copy()
    trb = df[df["Gene"] == "TRB"].copy()

    paired_ids = set(tra["complex.id"]) & set(trb["complex.id"])

    def _meta(sub: pd.DataFrame) -> dict:
        return {
            "epitope":   sub["Epitope"].astype(str).str.strip(),
            "antigen":   sub["Epitope gene"].astype(str).str.strip(),
            "pathogen":  sub["Epitope species"].astype(str).str.strip(),
            "HLA":       sub["MHC A"].astype(str).str.strip(),
            "source_db": "VDJdb",
        }

    parts = []

    # Paired: TRB rows joined with their TRA partner on complex.id
    trb_p = trb[trb["complex.id"].isin(paired_ids)].copy()
    tra_lookup = (
        tra[tra["complex.id"].isin(paired_ids)][["complex.id", "CDR3"]]
        .rename(columns={"CDR3": "cdr3a"})
        .drop_duplicates("complex.id")
    )
    merged = trb_p.merge(tra_lookup, on="complex.id", how="left")
    parts.append(pd.DataFrame({
        "cdr3b": merged["CDR3"].str.upper().str.strip(),
        "cdr3a": merged["cdr3a"].str.upper().str.strip(),
        **_meta(merged),
    }))

    # Unpaired TRB rows (no TRA partner)
    trb_u = trb[~trb["complex.id"].isin(paired_ids)]
    if len(trb_u):
        parts.append(pd.DataFrame({
            "cdr3b": trb_u["CDR3"].str.upper().str.strip(),
            "cdr3a": pd.NA,
            **_meta(trb_u),
        }))

    # Unpaired TRA rows (no TRB partner)
    tra_u = tra[~tra["complex.id"].isin(paired_ids)]
    if len(tra_u):
        parts.append(pd.DataFrame({
            "cdr3b": pd.NA,
            "cdr3a": tra_u["CDR3"].str.upper().str.strip(),
            **_meta(tra_u),
        }))

    return pd.concat(parts, ignore_index=True)


def _load_slim(path: Path, min_score: int) -> pd.DataFrame:
    """Parse the VDJdb slim build (columns: gene, cdr3, species, vdjdb.score, …)."""
    df = pd.read_csv(path, sep="\t", low_memory=False)

    df = df[df["species"] == "HomoSapiens"].copy()
    df = df[df["vdjdb.score"].fillna(0).astype(int) >= min_score]

    tra = df[df["gene"] == "TRA"].copy()
    trb = df[df["gene"] == "TRB"].copy()

    paired_ids = set(tra["complex.id"]) & set(trb["complex.id"])

    def _meta(sub: pd.DataFrame) -> dict:
        return {
            "epitope":   sub["antigen.epitope"].astype(str).str.strip(),
            "antigen":   sub["antigen.gene"].astype(str).str.strip(),
            "pathogen":  sub["antigen.species"].astype(str).str.strip(),
            "HLA":       sub["mhc.a"].astype(str).str.strip(),
            "source_db": "VDJdb",
        }

    parts = []

    trb_p = trb[trb["complex.id"].isin(paired_ids)].copy()
    tra_lookup = (
        tra[tra["complex.id"].isin(paired_ids)][["complex.id", "cdr3"]]
        .rename(columns={"cdr3": "cdr3a"})
        .drop_duplicates("complex.id")
    )
    merged = trb_p.merge(tra_lookup, on="complex.id", how="left")
    parts.append(pd.DataFrame({
        "cdr3b": merged["cdr3"].str.upper().str.strip(),
        "cdr3a": merged["cdr3a"].str.upper().str.strip(),
        **_meta(merged),
    }))

    trb_u = trb[~trb["complex.id"].isin(paired_ids)]
    if len(trb_u):
        parts.append(pd.DataFrame({
            "cdr3b": trb_u["cdr3"].str.upper().str.strip(),
            "cdr3a": pd.NA,
            **_meta(trb_u),
        }))

    tra_u = tra[~tra["complex.id"].isin(paired_ids)]
    if len(tra_u):
        parts.append(pd.DataFrame({
            "cdr3b": pd.NA,
            "cdr3a": tra_u["cdr3"].str.upper().str.strip(),
            **_meta(tra_u),
        }))

    return pd.concat(parts, ignore_index=True)


def load(
    path: Path = VDJDB_PATH,
    min_score: int = 0,
) -> pd.DataFrame:
    # Auto-detect format from column names; fall back to slim build if needed
    if path.exists():
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
        if "CDR3" in header and "Gene" in header:
            out = _load_full(path, min_score)
        else:
            out = _load_slim(path, min_score)
    elif VDJDB_SLIM_PATH.exists():
        print(f"    [VDJdb] Full file not found, falling back to {VDJDB_SLIM_PATH}", flush=True)
        out = _load_slim(VDJDB_SLIM_PATH, min_score)
    else:
        print("    [VDJdb] Warning: no VDJdb file found — returning empty table", flush=True)
        return pd.DataFrame(
            columns=["cdr3b", "cdr3a", "epitope", "antigen", "pathogen", "HLA", "source_db"]
        )

    # Validate amino-acid strings; allow NA (unpaired chain placeholder)
    cdr3b_ok = out["cdr3b"].str.match(_AA_PAT, na=True)
    cdr3a_ok = out["cdr3a"].str.match(_AA_PAT, na=True)
    out = out[cdr3b_ok & cdr3a_ok]

    return out.drop_duplicates(subset=["cdr3b", "cdr3a", "epitope"]).reset_index(drop=True)
