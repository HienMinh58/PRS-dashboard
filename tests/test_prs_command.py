"""
Tests for PRS-CSx command construction and pipeline safety guards.

NOTE: src.prs_methods imports streamlit at module level, so we
install a lightweight stub before importing it.
"""
import sys
import types
import pytest
from unittest.mock import patch, MagicMock

# ── Stub streamlit before importing prs_methods ──────────────────
_st_stub = types.ModuleType("streamlit")
_st_stub.empty = MagicMock
_st_stub.info = MagicMock()
_st_stub.error = MagicMock()
_st_stub.progress = MagicMock(return_value=MagicMock())
_st_stub.spinner = MagicMock(return_value=MagicMock(__enter__=MagicMock(), __exit__=MagicMock()))
sys.modules.setdefault("streamlit", _st_stub)

from src.prs_methods import run_prs_csx


class TestCommandConstruction:
    """Verify the CLI command is built correctly from user inputs."""

    @patch("src.prs_methods.subprocess.Popen")
    @patch("src.prs_methods.os.makedirs")
    @patch("glob.glob", return_value=[])
    def test_single_pop_command(self, mock_glob, mock_mkdirs, mock_popen):
        proc = MagicMock()
        proc.stdout = iter([])
        proc.wait.return_value = None
        proc.returncode = 1
        mock_popen.return_value = proc

        with pytest.raises(Exception):
            run_prs_csx(
                gwas_files=["/data/eur.txt"],
                gwas_pops=["EUR"],
                gwas_ns=["100000"],
                target_plink="/data/target",
                params={"phi": "1e-2", "a": "1.0"},
            )

        cmd = mock_popen.call_args[0][0]
        assert cmd[0] == "python"
        assert "PRScsx.py" in cmd[1]
        assert "--ref_dir=/app/ld_reference" in cmd
        assert "--bim_prefix=/data/target" in cmd
        assert "--sst_file=/data/eur.txt" in cmd
        assert "--pop=EUR" in cmd

    @patch("src.prs_methods.subprocess.Popen")
    @patch("src.prs_methods.os.makedirs")
    @patch("glob.glob", return_value=[])
    def test_multi_pop_command(self, mock_glob, mock_mkdirs, mock_popen):
        proc = MagicMock()
        proc.stdout = iter([])
        proc.wait.return_value = None
        proc.returncode = 1
        mock_popen.return_value = proc

        with pytest.raises(Exception):
            run_prs_csx(
                gwas_files=["/data/eur.txt", "/data/afr.txt"],
                gwas_pops=["EUR", "AFR"],
                gwas_ns=["100000", "50000"],
                target_plink="/data/target",
                params={"phi": "1e-2", "a": "1.0"},
            )

        cmd = mock_popen.call_args[0][0]
        assert "--sst_file=/data/eur.txt,/data/afr.txt" in cmd
        assert "--pop=EUR,AFR" in cmd
        assert "--n_gwas=100000,50000" in cmd


class TestPipelineSafety:
    """Unimplemented methods must raise, never return fake data."""

    def test_tl_prs_not_implemented(self):
        from src.prs_methods import run_tl_prs
        with pytest.raises(NotImplementedError):
            run_tl_prs([], "", {})

    def test_ct_sleb_not_implemented(self):
        from src.prs_methods import run_ct_sleb
        with pytest.raises(NotImplementedError):
            run_ct_sleb([], "", {})

    def test_me_bayes_sl_not_implemented(self):
        from src.prs_methods import run_me_bayes_sl
        with pytest.raises(NotImplementedError):
            run_me_bayes_sl([], "", {})
