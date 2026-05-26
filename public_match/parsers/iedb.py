import pandas as pd
from pathlib import Path

IEDB_PATH = Path("Databases/IEDB/iedb.xlsx")

_AA_PAT = r"^[ACDEFGHIKLMNPQRSTVWY]+$"


def _coalesce_cdr3(curated: pd.Series, calculated: pd.Series) -> pd.Series:
    """Use curated CDR3 when available, fall back to calculated."""
    return curated.where(curated.notna() & (curated != ""), calculated)


def load(path: Path = IEDB_PATH) -> pd.DataFrame:
    df = pd.read_excel(path, dtype_backend="numpy_nullable")

    # Use fillna(False) to collapse nullable booleans to plain bool —
    # without this, ~c1_beta is <NA> (not False) for rows where Chain 1 Type
    # is null, causing c2_beta & ~c1_beta to evaluate as <NA> and silently
    # drop ~156k valid beta-chain rows (PyArrow three-value logic).
    #
    # We do NOT pre-filter to alphabeta: gammadelta/construct rows have no
    # beta chains and fall out naturally; restricting to alphabeta is only
    # relevant for paired-chain searches, which load_databases handles at
    # the combined level by requiring cdr3a to be present.
    c1_is_beta  = (df["Chain 1 - Type"] == "beta").fillna(False)
    c2_is_beta  = (df["Chain 2 - Type"] == "beta").fillna(False)
    c1_is_alpha = (df["Chain 1 - Type"] == "alpha").fillna(False)
    c2_is_alpha = (df["Chain 2 - Type"] == "alpha").fillna(False)

    def _cdr3(sub, curated_col, calc_col):
        return _coalesce_cdr3(
            sub[curated_col].astype(str),
            sub[calc_col].astype(str),
        ).str.upper().str.strip()

    def _meta(sub):
        return {
            "epitope":   sub["Epitope - Name"].astype(str).str.strip(),
            "antigen":   sub["Epitope - Source Molecule"].astype(str).str.strip(),
            "pathogen":  sub["Epitope - Source Organism"].astype(str).str.strip(),
            "HLA":       sub["Assay - MHC Allele Names"].astype(str).str.strip(),
            "source_db": "IEDB",
        }

    parts = []

    # Rows where Chain 2 is the beta chain (the dominant orientation in IEDB)
    if c2_is_beta.any():
        sub = df[c2_is_beta]
        parts.append(pd.DataFrame({
            "cdr3b": _cdr3(sub, "Chain 2 - CDR3 Curated", "Chain 2 - CDR3 Calculated"),
            "cdr3a": _cdr3(sub, "Chain 1 - CDR3 Curated", "Chain 1 - CDR3 Calculated").where(
                c1_is_alpha[c2_is_beta].values
            ),
            **_meta(sub),
        }))

    # Rows where Chain 1 is the beta chain and Chain 2 is not (avoids double-counting)
    mask_c1b = c1_is_beta & ~c2_is_beta
    if mask_c1b.any():
        sub = df[mask_c1b]
        parts.append(pd.DataFrame({
            "cdr3b": _cdr3(sub, "Chain 1 - CDR3 Curated", "Chain 1 - CDR3 Calculated"),
            "cdr3a": _cdr3(sub, "Chain 2 - CDR3 Curated", "Chain 2 - CDR3 Calculated").where(
                c2_is_alpha[mask_c1b].values
            ),
            **_meta(sub),
        }))

    if not parts:
        return pd.DataFrame(
            columns=["cdr3b", "cdr3a", "epitope", "antigen", "pathogen", "HLA", "source_db"]
        )

    out = pd.concat(parts, ignore_index=True)
    out = out[out["cdr3b"].str.match(_AA_PAT, na=False)]
    return out.drop_duplicates(subset=["cdr3b", "epitope"]).reset_index(drop=True)
