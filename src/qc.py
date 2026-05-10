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
    """
    Reads and standardises a GWAS summary statistics file.

    Args:
        file_path (str): Path to the GWAS file.

    Returns:
        pd.DataFrame: A standardised DataFrame with recognised columns (SNP, CHR, BP, A1, A2, etc.).
    """
    df = pd.read_csv(file_path, sep="\t")
    return standardize_gwas_columns(df)


def read_plink_bim(bim_path_or_prefix):
    """
    Reads a PLINK .bim file.

    Args:
        bim_path_or_prefix (str): Path to the .bim file or the PLINK prefix.

    Returns:
        pd.DataFrame: A DataFrame with standard BIM columns (CHR, SNP, CM, BP, A1, A2).
    """
    bim_path = _resolve_bim_path(bim_path_or_prefix)
    df = pd.read_csv(bim_path, sep=r"\s+", header=None, names=BIM_COLUMNS)
    return _clean_variant_columns(df)


def find_invalid_bim_variants(bim_path_or_prefix):
    """
    Identifies variants in a BIM file with invalid, missing, or duplicate alleles.

    Args:
        bim_path_or_prefix (str): Path to the .bim file or the PLINK prefix.

    Returns:
        pd.DataFrame: A DataFrame containing only the invalid variant rows.
    """
    bim_path = _resolve_bim_path(bim_path_or_prefix)
    df = pd.read_csv(bim_path, sep=r"\s+", header=None, names=BIM_COLUMNS)
    a1 = df["A1"].astype(str).str.upper()
    a2 = df["A2"].astype(str).str.upper()
    invalid_mask = (a1 == a2) | (~a1.isin(VALID_ALLELES)) | (~a2.isin(VALID_ALLELES))
    return df.loc[invalid_mask].copy()


def bim_has_non_numeric_chromosomes(bim_path_or_prefix):
    """
    Checks if a BIM file contains non-numeric chromosome identifiers (e.g., 'X', 'Y', 'MT').

    Args:
        bim_path_or_prefix (str): Path to the .bim file or the PLINK prefix.

    Returns:
        bool: True if non-numeric chromosomes are present, False otherwise.
    """
    bim_path = _resolve_bim_path(bim_path_or_prefix)
    df = pd.read_csv(bim_path, sep=r"\s+", header=None, names=BIM_COLUMNS, usecols=[0])
    return pd.to_numeric(df["CHR"], errors="coerce").isna().any()


def clean_plink_invalid_alleles(prefix, out_prefix=None, plink_cmd="plink2", keep_chrom=None):
    """
    Removes variants with invalid/duplicate alleles from PLINK binary files.

    Args:
        prefix (str): Input PLINK file prefix.
        out_prefix (str, optional): Output PLINK file prefix. Defaults to prefix.allele_clean.
        plink_cmd (str): Command to run PLINK (default 'plink2').
        keep_chrom (str/int, optional): If provided, filters to this chromosome (or "autosome").

    Returns:
        tuple: (clean_prefix, summary_dict)
            clean_prefix is the prefix of the resulting files.
            summary_dict contains metadata about the operation.
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
    """
    Matches GWAS summary statistics to BIM variants by SNP ID.

    Args:
        gwas_df (pd.DataFrame): Standardised GWAS DataFrame.
        bim_df (pd.DataFrame): Standardised BIM DataFrame.

    Returns:
        pd.DataFrame: Merged DataFrame containing only variants found in both.
    """
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
    Executes the standard Quality Control (QC) v1 pipeline.

    The pipeline includes:
    1. Standardising GWAS columns.
    2. Matching GWAS variants to the target BIM file.
    3. Filtering out ambiguous SNPs (optional).
    4. Filtering SNPs by Minor Allele Frequency (MAF).

    Args:
        gwas_path (str): Path to GWAS file.
        bim_path_or_prefix (str): Target genotype path/prefix.
        remove_ambiguous (bool): Whether to remove A/T and C/G SNPs.
        maf_threshold (float): MAF threshold for filtering.

    Returns:
        tuple: (matched_df, summary_dict)
            matched_df is the clean DataFrame after QC.
            summary_dict contains counts of filtered/remaining SNPs.
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
