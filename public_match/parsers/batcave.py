import pandas as pd
from pathlib import Path

BATCAVE_MHCI_PATH = Path("Databases/BATCAVE/TCR_pMHCI_mutational_scan.csv")
BATCAVE_MHCII_PATH = Path("Databases/BATCAVE/TCR_pMHCII_mutational_scan.csv")

# BATCAVE mixes 13 assay types on incompatible scales (0–1 normalized vs 0–34,770 absolute).
# A single absolute threshold incorrectly drops all normalized-assay entries where 1.0 = full
# native-peptide response. Since peptide == index_peptide already guarantees a real binder,
# we only require strictly positive activation (> 0).
_MIN_ACTIVITY = 0.0   # strictly positive; applied as > 0


def _load_file(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)

    # filter: human TCRs, native epitope only, positive activation
    df = df[df["tcr_source_organism"].str.lower() == "human"].copy()
    df = df[df["peptide"] == df["index_peptide"]].copy()
    df = df[df["peptide_activity"] > _MIN_ACTIVITY].copy()
    df = df[df["cdr3b"].notna()].copy()

    return pd.DataFrame({
        "cdr3b":    df["cdr3b"].str.upper().str.strip(),
        "cdr3a":    df["cdr3a"].astype(str).str.upper().str.strip(),
        "epitope":  df["index_peptide"].str.strip(),
        "antigen":  df["peptide_type"].str.strip(),
        "pathogen": df["peptide_type"].str.strip(),
        "HLA":      df["mhc"].str.strip(),
        "source_db": "BATCAVE",
    })


def load(
    mhci_path: Path = BATCAVE_MHCI_PATH,
    mhcii_path: Path = BATCAVE_MHCII_PATH,
) -> pd.DataFrame:
    parts = []
    if mhci_path.exists():
        parts.append(_load_file(mhci_path))
    if mhcii_path.exists():
        parts.append(_load_file(mhcii_path))

    if not parts:
        return pd.DataFrame(columns=["cdr3b", "epitope", "antigen", "pathogen", "HLA", "source_db"])

    out = pd.concat(parts, ignore_index=True)
    out = out[out["cdr3b"].str.match(r"^[ACDEFGHIKLMNPQRSTVWY]+$", na=False)]
    return out.drop_duplicates(subset=["cdr3b", "epitope"]).reset_index(drop=True)
