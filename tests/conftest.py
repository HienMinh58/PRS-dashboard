"""
Shared pytest fixtures for PRS Dashboard tests.

Provides lightweight, deterministic test data for GWAS validation,
PLINK file handling, and phenotype checks. No real genetic data is
required — everything is synthesised in-memory or via tmp_path.
"""
import os
import struct
import pytest
import pandas as pd
import numpy as np


# ── GWAS fixtures ────────────────────────────────────────────────

@pytest.fixture
def gwas_valid_df():
    """Minimal valid GWAS summary-stats DataFrame."""
    return pd.DataFrame({
        "SNP":  ["rs1", "rs2", "rs3", "rs4", "rs5"],
        "A1":   ["A",   "C",   "G",   "T",   "A"],
        "A2":   ["G",   "T",   "A",   "C",   "C"],
        "BETA": [0.05, -0.12,  0.03,  0.08, -0.01],
        "P":    [0.01,  0.04,  0.50,  0.001, 0.90],
    })


@pytest.fixture
def gwas_valid_file(tmp_path, gwas_valid_df):
    """Write a valid GWAS file to disk and return its path."""
    path = tmp_path / "valid_gwas.txt"
    gwas_valid_df.to_csv(path, sep="\t", index=False)
    return str(path)


@pytest.fixture
def gwas_missing_cols_file(tmp_path):
    """GWAS file that is missing the P-value column."""
    df = pd.DataFrame({
        "SNP":  ["rs1", "rs2"],
        "A1":   ["A",   "C"],
        "A2":   ["G",   "T"],
        "BETA": [0.05, -0.12],
        # P column intentionally missing
    })
    path = tmp_path / "missing_p_gwas.txt"
    df.to_csv(path, sep="\t", index=False)
    return str(path)


@pytest.fixture
def gwas_with_ambiguous_snps_file(tmp_path):
    """GWAS containing ambiguous A/T and C/G SNPs that should be removed."""
    df = pd.DataFrame({
        "SNP":  ["rs_ok", "rs_AT", "rs_CG"],
        "A1":   ["A",     "A",     "C"],
        "A2":   ["C",     "T",     "G"],
        "BETA": [0.05,    0.10,    -0.07],
        "P":    [0.01,    0.02,     0.03],
    })
    path = tmp_path / "ambig_gwas.txt"
    df.to_csv(path, sep="\t", index=False)
    return str(path)


@pytest.fixture
def gwas_with_or_file(tmp_path):
    """GWAS using OR instead of BETA — validator should convert via log(OR)."""
    df = pd.DataFrame({
        "SNP":  ["rs1", "rs2"],
        "A1":   ["A",   "C"],
        "A2":   ["G",   "T"],
        "OR":   [1.15,  0.90],
        "P":    [0.01,  0.04],
    })
    path = tmp_path / "or_gwas.txt"
    df.to_csv(path, sep="\t", index=False)
    return str(path)


@pytest.fixture
def gwas_duplicate_snps_file(tmp_path):
    """GWAS with duplicated SNP IDs that should be deduplicated."""
    df = pd.DataFrame({
        "SNP":  ["rs1", "rs1", "rs2"],
        "A1":   ["A",   "A",   "C"],
        "A2":   ["G",   "G",   "T"],
        "BETA": [0.05,  0.05, -0.12],
        "P":    [0.01,  0.01,  0.04],
    })
    path = tmp_path / "dup_gwas.txt"
    df.to_csv(path, sep="\t", index=False)
    return str(path)


# ── PLINK fixtures ───────────────────────────────────────────────

def _write_minimal_plink(directory, prefix, snp_ids, n_samples=5):
    """
    Write a minimal but structurally valid PLINK 1 binary fileset.

    The .bed file contains the correct magic bytes and mode byte,
    followed by zero-filled genotype data (all homozygous reference).
    """
    bim_path = os.path.join(directory, f"{prefix}.bim")
    fam_path = os.path.join(directory, f"{prefix}.fam")
    bed_path = os.path.join(directory, f"{prefix}.bed")

    # .bim — one row per SNP
    bim_rows = []
    for i, snp in enumerate(snp_ids):
        bim_rows.append(f"1\t{snp}\t0\t{(i+1)*1000}\tA\tG")
    with open(bim_path, "w") as f:
        f.write("\n".join(bim_rows) + "\n")

    # .fam — one row per sample
    fam_rows = []
    for i in range(n_samples):
        fam_rows.append(f"FAM{i}\tIID_{i}\t0\t0\t1\t-9")
    with open(fam_path, "w") as f:
        f.write("\n".join(fam_rows) + "\n")

    # .bed — magic number (0x6c, 0x1b, 0x01) + zero-genotype bytes
    n_snps = len(snp_ids)
    bytes_per_snp = (n_samples + 3) // 4  # PLINK packs 4 samples per byte
    with open(bed_path, "wb") as f:
        f.write(struct.pack("BBB", 0x6C, 0x1B, 0x01))
        f.write(b"\x00" * (n_snps * bytes_per_snp))

    return os.path.join(directory, prefix)


@pytest.fixture
def plink_valid_prefix(tmp_path):
    """Valid minimal PLINK fileset with 5 SNPs / 5 samples."""
    return _write_minimal_plink(
        str(tmp_path), "test_target",
        snp_ids=["rs1", "rs2", "rs3", "rs4", "rs5"],
    )


@pytest.fixture
def plink_missing_bim(tmp_path):
    """PLINK fileset with the .bim file removed."""
    prefix = _write_minimal_plink(
        str(tmp_path), "missing_bim",
        snp_ids=["rs1"],
    )
    os.remove(f"{prefix}.bim")
    return prefix


# ── Phenotype fixtures ───────────────────────────────────────────

@pytest.fixture
def pheno_valid_file(tmp_path):
    """Valid phenotype file (FID IID PHENO)."""
    rows = ["FID\tIID\tPHENO"]
    for i in range(5):
        rows.append(f"FAM{i}\tIID_{i}\t{np.random.normal():.4f}")
    path = tmp_path / "test.pheno"
    path.write_text("\n".join(rows) + "\n")
    return str(path)


@pytest.fixture
def pheno_single_col_file(tmp_path):
    """Phenotype file with only 1 column — should fail validation."""
    path = tmp_path / "bad.pheno"
    path.write_text("value\n1.0\n2.0\n")
    return str(path)
