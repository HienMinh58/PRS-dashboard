"""
Tests for src.utils helper functions.
"""
import numpy as np
import pandas as pd
import pytest
from src.utils import generate_mock_phenotype, calculate_metrics, export_to_csv, export_to_excel


class TestGenerateMockPhenotype:
    def test_continuous_shape(self):
        y = generate_mock_phenotype(n_samples=100, binary=False)
        assert len(y) == 100

    def test_binary_values(self):
        y = generate_mock_phenotype(n_samples=200, binary=True)
        assert set(np.unique(y)).issubset({0, 1})

    def test_deterministic(self):
        a = generate_mock_phenotype(50, binary=False)
        b = generate_mock_phenotype(50, binary=False)
        np.testing.assert_array_equal(a, b)


class TestCalculateMetrics:
    def test_regression_metrics(self):
        y_true = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        y_pred = np.array([1.1, 2.1, 2.9, 4.2, 4.8])
        m = calculate_metrics(y_true, y_pred, is_binary=False)
        assert "R2" in m
        assert m["R2"] > 0.9

    def test_binary_metrics(self):
        y_true = np.array([0, 0, 1, 1])
        y_pred = np.array([0.1, 0.2, 0.8, 0.9])
        m = calculate_metrics(y_true, y_pred, is_binary=True)
        assert "AUC" in m
        assert m["AUC"] == 1.0


class TestExport:
    def test_csv_roundtrip(self):
        df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
        csv_bytes = export_to_csv(df)
        assert isinstance(csv_bytes, bytes)
        assert b"a,b" in csv_bytes

    @pytest.mark.skipif(
        not pytest.importorskip("openpyxl", reason="openpyxl not installed"),
        reason="openpyxl not installed",
    )
    def test_excel_is_bytes(self):
        df = pd.DataFrame({"x": [10]})
        xlsx = export_to_excel(df)
        assert isinstance(xlsx, bytes)
        assert xlsx[:2] == b"PK"
