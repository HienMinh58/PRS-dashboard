import pytest
from unittest.mock import patch, MagicMock
from src.ld_compute import compute_ld_from_target

@patch("src.ld_compute.subprocess.run")
@patch("src.ld_compute.os.makedirs")
@patch("src.ld_compute.os.path.exists")
def test_compute_ld_command_construction(mock_exists, mock_makedirs, mock_run):
    # Setup mocks
    # We need to simulate subprocess.run for checking 'plink2 --version' and the actual command.
    # We will just let it succeed
    mock_run.return_value = MagicMock(returncode=0)
    
    # Simulate that the output file exists
    mock_exists.return_value = True 
    
    # Test single chromosome
    summary = compute_ld_from_target("mock_target", "1", window_kb=1000, r2_threshold=0.1, out_dir="/tmp/ld")
    
    assert summary["success"] is True
    
    cmd = summary["command"]
    assert "plink2" in cmd or "plink" in cmd
    assert "--bfile" in cmd
    assert "mock_target" in cmd
    assert "--ld-window-kb" in cmd
    assert "1000" in cmd
    assert "--ld-window-r2" in cmd
    assert "0.1" in cmd
    # It should have either --r2-unphased or --r2
    assert any(x in cmd for x in ["--r2-unphased", "--r2"])
    assert "--chr" in cmd
    assert "1" in cmd
    
    # Test all autosomes
    summary_all = compute_ld_from_target("mock_target", "1-22")
    cmd_all = summary_all["command"]
    assert "--autosome" in cmd_all
    assert "--chr" not in cmd_all
