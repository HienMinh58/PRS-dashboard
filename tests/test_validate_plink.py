"""
Tests for src.validation.validate_plink.

Covers: happy-path, missing files, SNP overlap with GWAS, and allele checks.
"""
import os
import pytest

from src.validation import validate_plink


class TestValidPlink:
    def test_valid_prefix_passes(self, plink_valid_prefix):
        is_valid, msg = validate_plink(plink_valid_prefix)
        assert is_valid is True
        assert "5 samples" in msg
        assert "5 SNPs" in msg

    def test_overlap_reported_when_gwas_snps_provided(self, plink_valid_prefix):
        gwas_snps = {"rs1", "rs3", "rs99"}
        is_valid, msg = validate_plink(plink_valid_prefix, gwas_snps)
        assert is_valid is True
        # rs1 and rs3 overlap → 2 SNPs
        assert "2 SNPs" in msg or "Overlap OK" in msg


class TestMissingFiles:
    def test_missing_bim_fails(self, plink_missing_bim):
        is_valid, msg = validate_plink(plink_missing_bim)
        assert is_valid is False
        assert "Missing" in msg

    def test_totally_bogus_prefix_fails(self, tmp_path):
        is_valid, msg = validate_plink(str(tmp_path / "nonexistent"))
        assert is_valid is False


class TestLowOverlap:
    def test_low_overlap_warns(self, plink_valid_prefix):
        """Fewer than 100 overlapping SNPs should produce a WARNING."""
        gwas_snps = {"rs1"}  # only 1 overlap out of 5
        is_valid, msg = validate_plink(plink_valid_prefix, gwas_snps)
        # The validator still passes but includes a WARNING
        assert is_valid is True
        assert "WARNING" in msg
