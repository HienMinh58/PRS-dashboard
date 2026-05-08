import os
import subprocess

import pandas as pd


GWAS_COLUMN_ALIASES = {
    "SNP": ["SNP", "RSID", "ID", "MARKERNAME", "VARIANT_ID"],
    "CHR": ["CHR", "CHROM", "CHROMOSOME"],
    "BP": ["BP", "POS", "POSITION"],
    "A1": ["A1", "EA", "EFFECT_ALLELE", "ALT"],
    "A2": ["A2", "NEA", "OTHER_ALLELE", "REF"],
    "BETA": ["BETA", "B", "EFFECT"],
    "P": ["P", "PVAL", "P_VALUE", "PVALUE"],
    "MAF": ["MAF", "EAF", "AF", "FREQ"],
}

CORE_GWAS_COLUMNS = ["SNP", "CHR", "BP", "A1", "A2"]
OPTIONAL_GWAS_COLUMNS = ["BETA", "P", "MAF"]
BIM_COLUMNS = ["CHR", "SNP", "CM", "BP", "A1", "A2"]
VALID_ALLELES = {"A", "C", "G", "T"}

AMBIGUOUS_SNPS = {
    ("A", "T"),
    ("T", "A"),
    ("C", "G"),
    ("G", "C"),
}

def read_gwas_summary_stats(file_path):
    """Read GWAS summary statistics and standardise known PRS columns."""
    df = pd.read_csv(file_path, sep="\t")
    return standardize_gwas_columns(df)


def read_plink_bim(bim_path_or_prefix):
    """Read a PLINK .bim file and return standard BIM columns."""
    bim_path = _resolve_bim_path(bim_path_or_prefix)
    df = pd.read_csv(bim_path, sep=r"\s+", header=None, names=BIM_COLUMNS)
    return _clean_variant_columns(df)


def find_invalid_bim_variants(bim_path_or_prefix):
    """Return BIM rows with invalid, missing, or duplicate alleles."""
    bim_path = _resolve_bim_path(bim_path_or_prefix)
    df = pd.read_csv(bim_path, sep=r"\s+", header=None, names=BIM_COLUMNS)
    a1 = df["A1"].astype(str).str.upper()
    a2 = df["A2"].astype(str).str.upper()
    invalid_mask = (a1 == a2) | (~a1.isin(VALID_ALLELES)) | (~a2.isin(VALID_ALLELES))
    return df.loc[invalid_mask].copy()


def bim_has_non_numeric_chromosomes(bim_path_or_prefix):
    """Return True when BIM chromosome values include non-numeric labels like X."""
    bim_path = _resolve_bim_path(bim_path_or_prefix)
    df = pd.read_csv(bim_path, sep=r"\s+", header=None, names=BIM_COLUMNS, usecols=[0])
    return pd.to_numeric(df["CHR"], errors="coerce").isna().any()


def clean_plink_invalid_alleles(prefix, out_prefix=None, plink_cmd="plink2", keep_chrom=None):
    """
    Remove BIM variants with invalid/duplicate alleles using PLINK2.

    When ``keep_chrom`` is provided, the output is also restricted to that
    chromosome. Pass ``"autosome"`` to keep chromosomes 1-22.

    Returns ``(clean_prefix, summary)``. If no invalid variants are found,
    and no chromosome filter is requested, the original prefix is returned and
    no PLINK command is run.
    """
    invalid = find_invalid_bim_variants(prefix)
    summary = {
        "input_prefix": prefix,
        "output_prefix": out_prefix or f"{prefix}.allele_clean",
        "num_invalid_bim_variants": int(len(invalid)),
        "chrom_filter": keep_chrom,
        "exclude_file": None,
        "plink_command": None,
        "cleaned": False,
    }

    if invalid.empty and keep_chrom is None:
        return prefix, summary

    clean_prefix = summary["output_prefix"]
    os.makedirs(os.path.dirname(clean_prefix) or ".", exist_ok=True)

    cmd = [
        plink_cmd,
        "--bfile",
        prefix,
    ]
    exclude_file = None
    if not invalid.empty:
        exclude_file = f"{clean_prefix}.exclude_snps.txt"
        invalid["SNP"].to_csv(exclude_file, index=False, header=False)
        cmd.extend(["--exclude", exclude_file])
    if keep_chrom == "autosome":
        cmd.append("--autosome")
    elif keep_chrom is not None:
        cmd.extend(["--chr", str(keep_chrom)])
    cmd.extend(["--make-bed", "--out", clean_prefix])

    subprocess.run(cmd, check=True, capture_output=True, text=True)

    summary.update(
        {
            "exclude_file": exclude_file,
            "plink_command": cmd,
            "cleaned": True,
        }
    )
    return clean_prefix, summary


def match_gwas_to_bim(gwas_df, bim_df):
    """Match standardised GWAS and BIM variants by SNP ID only."""
    gwas_clean = _clean_gwas_frame(gwas_df)
    bim_clean = _clean_variant_columns(bim_df).drop_duplicates(subset=["SNP"])

    matched = pd.merge(
        gwas_clean,
        bim_clean[["SNP", "CHR", "BP", "A1", "A2"]],
        on="SNP",
        how="inner",
        suffixes=("", "_BIM"),
    )

    output_cols = CORE_GWAS_COLUMNS + [
        col for col in OPTIONAL_GWAS_COLUMNS if col in matched.columns
    ]
    audit_cols = ["CHR_BIM", "BP_BIM", "A1_BIM", "A2_BIM"]
    return matched[output_cols + audit_cols]


def run_qc_v1(gwas_path, bim_path_or_prefix, remove_ambiguous=True, maf_threshold=0.01):
    """
    Run QC v1: read GWAS, read BIM, match by SNP ID, and return matched data.

    This skeleton intentionally does not perform liftover, strand flipping,
    allele harmonisation, or LD computation.
    """
    gwas_df = read_gwas_summary_stats(gwas_path)
    bim_df = read_plink_bim(bim_path_or_prefix)
    matched_df = match_gwas_to_bim(gwas_df, bim_df)
    num_matched_snps = int(matched_df["SNP"].nunique())

    n_ambiguous_removed = 0
    if remove_ambiguous:
        allele_pairs = list(
            zip(
                matched_df["A1"].str.upper(),
                matched_df["A2"].str.upper(),
            )
        )
        # Optimizing code for big dataset
        is_ambiguous = pd.Series(
            [pair in AMBIGUOUS_SNPS for pair in allele_pairs],
            index=matched_df.index
        )
        n_ambiguous_removed = int(is_ambiguous.sum())

        matched_df = matched_df.loc[~is_ambiguous].copy()

    n_maf_filtered_snps = 0
    maf_warning = None
    if "MAF" in matched_df.columns:
        initial_len = len(matched_df)
        matched_df = matched_df.loc[
            matched_df["MAF"] >= maf_threshold
        ].copy()
        n_maf_filtered_snps = int(initial_len - len(matched_df))
    else:
        maf_warning = "MAF column not found. MAF filtering skipped."

    summary = {
        "num_gwas_snps": int(gwas_df["SNP"].nunique()),
        "num_bim_snps": int(bim_df["SNP"].nunique()),
        "num_matched_snps": num_matched_snps,
        "num_unmatched_gwas_snps": int(
            gwas_df["SNP"].nunique() - num_matched_snps
        ),
        "num_ambiguous_removed": n_ambiguous_removed,
        "n_maf_filtered_snps": n_maf_filtered_snps,
        "num_final_snps_after_qc": int(len(matched_df)),
    }
    
    if maf_warning:
        summary["warning"] = maf_warning

    return matched_df, summary


def standardize_gwas_columns(df):
    """Return a GWAS dataframe with recognised columns renamed to standard names."""
    rename_map = {}
    upper_to_original = {str(col).upper(): col for col in df.columns}

    for standard_name, aliases in GWAS_COLUMN_ALIASES.items():
        match = next((alias for alias in aliases if alias in upper_to_original), None)
        if match is not None:
            rename_map[upper_to_original[match]] = standard_name

    standardised = df.rename(columns=rename_map)
    missing = [col for col in CORE_GWAS_COLUMNS if col not in standardised.columns]
    if missing:
        raise ValueError(f"Missing required GWAS columns: {', '.join(missing)}")

    keep_cols = CORE_GWAS_COLUMNS + [
        col for col in OPTIONAL_GWAS_COLUMNS if col in standardised.columns
    ]
    return standardised[keep_cols]


def _clean_gwas_frame(df):
    clean = _clean_variant_columns(df)
    for col in ["BETA", "P", "MAF"]:
        if col in clean.columns:
            clean[col] = pd.to_numeric(clean[col], errors="coerce")
    return clean.drop_duplicates(subset=["SNP"])


def _clean_variant_columns(df):
    clean = df.copy()
    clean = clean.dropna(subset=CORE_GWAS_COLUMNS)
    clean["SNP"] = clean["SNP"].astype(str)
    clean["CHR"] = clean["CHR"].astype(str)
    clean["BP"] = pd.to_numeric(clean["BP"], errors="coerce")
    clean["A1"] = clean["A1"].astype(str).str.upper()
    clean["A2"] = clean["A2"].astype(str).str.upper()
    clean = clean.dropna(subset=["BP"])
    clean["BP"] = clean["BP"].astype(int)
    return clean

def _resolve_bim_path(bim_path_or_prefix):
    if str(bim_path_or_prefix).endswith(".bim"):
        return bim_path_or_prefix

    bim_path = f"{bim_path_or_prefix}.bim"
    if os.path.exists(bim_path):
        return bim_path

    raise FileNotFoundError(f"BIM file not found: {bim_path_or_prefix}")

def is_ambiguous_snp(a1, a2):
    if a1 is None or a2 is None:
        return False

    return (
        str(a1).upper(),
        str(a2).upper()
    ) in AMBIGUOUS_SNPS
