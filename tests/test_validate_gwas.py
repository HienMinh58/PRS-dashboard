"""
Tests for src.validation.validate_gwas.

Covers: happy-path, missing columns, ambiguous SNPs, OR→BETA conversion,
duplicate SNPs, invalid alleles, and out-of-range P-values.
"""
import os
import numpy as np
import pandas as pd
import pytest

from src.validation import validate_gwas


# ── Happy path ───────────────────────────────────────────────────

class TestValidGWAS:
    def test_returns_valid(self, gwas_valid_file):
        is_valid, msg, count, processed, snps = validate_gwas(gwas_valid_file)
        assert is_valid is True

    def test_snp_count(self, gwas_valid_file):
        _, _, count, _, snps = validate_gwas(gwas_valid_file)
        assert count == 5
        assert snps == {"rs1", "rs2", "rs3", "rs4", "rs5"}

    def test_processed_file_exists(self, gwas_valid_file):
        _, _, _, processed, _ = validate_gwas(gwas_valid_file)
        assert os.path.isfile(processed)

    def test_processed_has_standard_columns(self, gwas_valid_file):
        _, _, _, processed, _ = validate_gwas(gwas_valid_file)
        df = pd.read_csv(processed, sep="\t")
        assert list(df.columns) == ["SNP", "A1", "A2", "BETA", "P"]


# ── Missing columns ─────────────────────────────────────────────

class TestMissingColumns:
    def test_missing_p_column_fails(self, gwas_missing_cols_file):
        is_valid, msg, *_ = validate_gwas(gwas_missing_cols_file)
        assert is_valid is False
        assert "P-value" in msg


# ── Ambiguous SNPs ───────────────────────────────────────────────

class TestAmbiguousSNPs:
    def test_ambiguous_snps_removed(self, gwas_with_ambiguous_snps_file):
        is_valid, msg, count, _, snps = validate_gwas(gwas_with_ambiguous_snps_file)
        assert is_valid is True
        # Only "rs_ok" should survive; rs_AT (A/T) and rs_CG (C/G) are ambiguous
        assert count == 1
        assert "rs_ok" in snps
        assert "rs_AT" not in snps
        assert "rs_CG" not in snps


# ── OR → BETA conversion ────────────────────────────────────────

class TestORConversion:
    def test_or_converted_to_beta(self, gwas_with_or_file):
        is_valid, _, _, processed, _ = validate_gwas(gwas_with_or_file)
        assert is_valid is True
        df = pd.read_csv(processed, sep="\t")
        # BETA should be log(OR)
        assert np.isclose(df.loc[0, "BETA"], np.log(1.15), atol=1e-6)
        assert np.isclose(df.loc[1, "BETA"], np.log(0.90), atol=1e-6)


# ── Duplicate SNPs ───────────────────────────────────────────────

class TestDuplicateSNPs:
    def test_duplicates_removed(self, gwas_duplicate_snps_file):
        is_valid, msg, count, _, snps = validate_gwas(gwas_duplicate_snps_file)
        assert is_valid is True
        assert count == 2  # rs1 (deduped) + rs2
        assert "1 duplicates" in msg or "duplicate" in msg.lower()


# ── Edge cases ───────────────────────────────────────────────────

class TestEdgeCases:
    def test_invalid_p_values_filtered(self, tmp_path):
        """P-values outside [0, 1] should be dropped."""
        df = pd.DataFrame({
            "SNP":  ["rs_ok", "rs_neg", "rs_big"],
            "A1":   ["A",     "C",      "G"],
            "A2":   ["C",     "A",      "T"],
            "BETA": [0.1,     0.2,      0.3],
            "P":    [0.05,   -0.01,     1.5],
        })
        path = tmp_path / "bad_p.txt"
        df.to_csv(path, sep="\t", index=False)
        is_valid, _, count, _, _ = validate_gwas(str(path))
        assert is_valid is True
        assert count == 1  # only rs_ok survives

    def test_nonexistent_file_returns_error(self, tmp_path):
        fake = str(tmp_path / "nope.txt")
        is_valid, msg, *_ = validate_gwas(fake)
        assert is_valid is False

    def test_empty_file_returns_error(self, tmp_path):
        path = tmp_path / "empty.txt"
        path.write_text("")
        is_valid, msg, *_ = validate_gwas(str(path))
        assert is_valid is False
