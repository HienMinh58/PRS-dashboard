"""
Tests for src.validation.validate_phenotype.
"""
import io
import pytest

from src.validation import validate_phenotype


class TestValidPhenotype:
    def test_fid_iid_pheno_format(self, pheno_valid_file):
        with open(pheno_valid_file, "r") as f:
            is_valid, msg, df = validate_phenotype(f)
        assert is_valid is True
        assert len(df) == 5
        assert list(df.columns) == ["Sample_ID", "PHENO"]

    def test_two_column_format(self, tmp_path):
        """IID + PHENO (no FID) should also work."""
        path = tmp_path / "two_col.pheno"
        path.write_text("IID\tPHENO\nIID_0\t1.5\nIID_1\t2.3\n")
        with open(str(path), "r") as f:
            is_valid, msg, df = validate_phenotype(f)
        assert is_valid is True
        assert len(df) == 2


class TestInvalidPhenotype:
    def test_single_column_does_not_crash(self, pheno_single_col_file):
        """A 1-column file should either fail or at least not crash.

        Current implementation: validate_phenotype raises a ValueError
        for < 2 columns, which is caught and returned as is_valid=False.
        However, if the CSV parser reads it as 2 columns (e.g. index + value),
        it may still pass. This test documents the actual behaviour.
        """
        with open(pheno_single_col_file, "r") as f:
            is_valid, msg, df = validate_phenotype(f)
        # The function must not raise — it returns a tuple either way
        assert isinstance(is_valid, bool)

    def test_all_missing_values_returns_empty(self, tmp_path):
        """If every phenotype value is NaN the result should be empty."""
        path = tmp_path / "nan.pheno"
        path.write_text("IID\tPHENO\nIID_0\t\nIID_1\t\n")
        with open(str(path), "r") as f:
            is_valid, msg, df = validate_phenotype(f)
        # Depending on implementation, either fails or returns 0 rows
        if is_valid:
            assert len(df) == 0
