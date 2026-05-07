import pandas as pd

from src.qc import (
    find_invalid_bim_variants,
    match_gwas_to_bim,
    read_plink_bim,
    run_qc_v1,
    standardize_gwas_columns,
    is_ambiguous_snp,
)


def test_standardize_gwas_columns_with_aliases():
    df = pd.DataFrame(
        {
            "rsid": ["rs1"],
            "chromosome": [1],
            "pos": [1000],
            "effect_allele": ["a"],
            "other_allele": ["g"],
            "beta": [0.1],
            "pval": [0.05],
            "eaf": [0.2],
        }
    )

    standardised = standardize_gwas_columns(df)

    assert list(standardised.columns) == ["SNP", "CHR", "BP", "A1", "A2", "BETA", "P", "MAF"]


def test_read_plink_bim_accepts_prefix(plink_valid_prefix):
    bim = read_plink_bim(plink_valid_prefix)

    assert list(bim.columns) == ["CHR", "SNP", "CM", "BP", "A1", "A2"]
    assert len(bim) == 5


def test_find_invalid_bim_variants(tmp_path):
    bim_path = tmp_path / "target.bim"
    bim_path.write_text(
        "\n".join(
            [
                "1 rs_ok 0 100 A G",
                "1 rs_zero 0 200 0 0",
                "1 rs_same 0 300 C C",
                "1 rs_bad 0 400 A I",
            ]
        )
        + "\n"
    )

    invalid = find_invalid_bim_variants(str(bim_path))

    assert list(invalid["SNP"]) == ["rs_zero", "rs_same", "rs_bad"]


def test_match_gwas_to_bim_by_snp_id_only(plink_valid_prefix):
    gwas = pd.DataFrame(
        {
            "SNP": ["rs1", "rs3", "rs99"],
            "CHR": [2, 2, 2],
            "BP": [10, 30, 990],
            "A1": ["T", "C", "A"],
            "A2": ["C", "A", "G"],
            "BETA": [0.1, -0.2, 0.3],
            "P": [0.01, 0.02, 0.03],
        }
    )
    bim = read_plink_bim(plink_valid_prefix)

    matched = match_gwas_to_bim(gwas, bim)

    assert list(matched["SNP"]) == ["rs1", "rs3"]
    assert "CHR_BIM" in matched.columns
    assert "A1_BIM" in matched.columns


def test_run_qc_v1_returns_matched_dataframe_and_summary(tmp_path, plink_valid_prefix):
    gwas = pd.DataFrame(
        {
            "SNP": ["rs1", "rs2", "rs_missing"],
            "CHR": [1, 1, 1],
            "BP": [1000, 2000, 9999],
            "A1": ["A", "C", "G"],
            "A2": ["G", "T", "A"],
            "BETA": [0.1, 0.2, 0.3],
        }
    )
    gwas_path = tmp_path / "gwas.tsv"
    gwas.to_csv(gwas_path, sep="\t", index=False)

    matched, summary = run_qc_v1(str(gwas_path), plink_valid_prefix)

    assert list(matched["SNP"]) == ["rs1", "rs2"]
    assert summary == {
        "num_gwas_snps": 3,
        "num_bim_snps": 5,
        "num_matched_snps": 2,
        "num_unmatched_gwas_snps": 1,
        "num_ambiguous_removed": 0,
        "num_final_snps_after_qc": 2,
    }

def test_ambiguous_pairs():
    assert is_ambiguous_snp("A", "T") is True
    assert is_ambiguous_snp("T", "A") is True
    assert is_ambiguous_snp("C", "G") is True
    assert is_ambiguous_snp("G", "C") is True


def test_non_ambiguous_pairs():
    assert is_ambiguous_snp("A", "C") is False
    assert is_ambiguous_snp("A", "G") is False
    assert is_ambiguous_snp("T", "C") is False


def test_lowercase_input():
    assert is_ambiguous_snp("a", "t") is True
    assert is_ambiguous_snp("c", "g") is True


def test_none_input():
    assert is_ambiguous_snp(None, "A") is False
    assert is_ambiguous_snp("A", None) is False

def test_run_qc_no_ambiguous(monkeypatch):

    # mock gwas
    gwas_df = pd.DataFrame({
        "SNP": ["rs1", "rs2", "rs3"],
    })

    # mock bim
    bim_df = pd.DataFrame({
        "SNP": ["rs1", "rs2", "rs3"],
    })

    # mock matched
    matched_df = pd.DataFrame({
        "SNP": ["rs1", "rs2", "rs3"],
        "A1": ["A", "C", "G"],
        "A2": ["T", "G", "A"],
    })

    # patch functions
    monkeypatch.setattr("src.qc.read_gwas_summary_stats", lambda x: gwas_df)
    monkeypatch.setattr("src.qc.read_plink_bim", lambda x: bim_df)
    monkeypatch.setattr("src.qc.match_gwas_to_bim", lambda a, b: matched_df)

    result_df, summary = run_qc_v1("gwas.txt", "bim", remove_ambiguous=False)

    assert len(result_df) == 3
    assert summary["num_matched_snps"] == 3
    assert summary["num_final_snps_after_qc"] == 3
