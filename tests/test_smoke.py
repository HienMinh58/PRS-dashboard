"""
Smoke tests — verify that key modules can be imported without error.

These run in CI without Streamlit, PLINK, or any external tool.
"""


def test_import_validation():
    from src.validation import validate_gwas, validate_plink, validate_ld_ref, validate_phenotype


def test_import_prs_methods():
    from src.prs_methods import execute_prs_pipeline, run_prs_csx


def test_import_ml_models():
    from src.ml_models import train_ml_models


def test_import_utils():
    from src.utils import (
        generate_mock_phenotype,
        export_to_csv,
        export_to_excel,
        calculate_metrics,
    )
