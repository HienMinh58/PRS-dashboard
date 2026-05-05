"""
Tests for src.validation.validate_ld_ref.

PRS-CSx requires a specific LD reference layout::

    ld_dir/
    ├── snpinfo_mult_1kg_hm3
    └── ldblk_1kg_{pop}/
        ├── ldblk_1kg_chr1.hdf5
        └── ...

These tests verify that validate_ld_ref rejects missing, empty, or
incomplete reference directories with clear error messages.
"""
import os
import pytest

from src.validation import validate_ld_ref


# ── Fixtures ─────────────────────────────────────────────────────

@pytest.fixture
def ld_complete(tmp_path):
    """Build a minimal but structurally valid LD reference tree."""
    root = tmp_path / "ld_reference"
    root.mkdir()

    # snpinfo file (empty placeholder is fine for structure check)
    (root / "snpinfo_mult_1kg_hm3").write_text("")

    # EUR population — chromosomes 1 and 22
    eur_dir = root / "ldblk_1kg_eur"
    eur_dir.mkdir()
    (eur_dir / "ldblk_1kg_chr1.hdf5").write_bytes(b"\x00")
    (eur_dir / "ldblk_1kg_chr22.hdf5").write_bytes(b"\x00")

    return str(root)


@pytest.fixture
def ld_no_snpinfo(tmp_path):
    """LD tree missing the snpinfo_mult_1kg_hm3 file."""
    root = tmp_path / "ld_reference"
    root.mkdir()
    eur_dir = root / "ldblk_1kg_eur"
    eur_dir.mkdir()
    (eur_dir / "ldblk_1kg_chr1.hdf5").write_bytes(b"\x00")
    return str(root)


@pytest.fixture
def ld_no_pop_dir(tmp_path):
    """LD tree has snpinfo but no population subdirectory."""
    root = tmp_path / "ld_reference"
    root.mkdir()
    (root / "snpinfo_mult_1kg_hm3").write_text("")
    return str(root)


@pytest.fixture
def ld_missing_chrom(tmp_path):
    """LD tree has the population dir but is missing a specific chromosome."""
    root = tmp_path / "ld_reference"
    root.mkdir()
    (root / "snpinfo_mult_1kg_hm3").write_text("")
    eur_dir = root / "ldblk_1kg_eur"
    eur_dir.mkdir()
    # Only chr1 exists — chr2 is missing
    (eur_dir / "ldblk_1kg_chr1.hdf5").write_bytes(b"\x00")
    return str(root)


@pytest.fixture
def ld_empty(tmp_path):
    """LD directory exists but is completely empty."""
    root = tmp_path / "ld_reference"
    root.mkdir()
    return str(root)


# ── Happy path ───────────────────────────────────────────────────

class TestValidLDRef:
    def test_valid_eur_chr1(self, ld_complete):
        is_valid, msg = validate_ld_ref("EUR", "1", ld_dir=ld_complete)
        assert is_valid is True
        assert "OK" in msg

    def test_valid_eur_chr22(self, ld_complete):
        is_valid, msg = validate_ld_ref("EUR", "22", ld_dir=ld_complete)
        assert is_valid is True

    def test_pop_case_insensitive(self, ld_complete):
        """Population codes should be lowercased internally."""
        is_valid, _ = validate_ld_ref("eur", "1", ld_dir=ld_complete)
        assert is_valid is True

    def test_chrom_as_int(self, ld_complete):
        """Chromosome can be passed as an integer."""
        is_valid, _ = validate_ld_ref("EUR", 1, ld_dir=ld_complete)
        assert is_valid is True


# ── Missing root directory ───────────────────────────────────────

class TestMissingDirectory:
    def test_nonexistent_dir_fails(self, tmp_path):
        fake = str(tmp_path / "does_not_exist")
        is_valid, msg = validate_ld_ref("EUR", "1", ld_dir=fake)
        assert is_valid is False
        assert "not found" in msg

    def test_message_contains_path(self, tmp_path):
        fake = str(tmp_path / "nope")
        _, msg = validate_ld_ref("EUR", "1", ld_dir=fake)
        assert "nope" in msg


# ── Empty directory ──────────────────────────────────────────────

class TestEmptyDirectory:
    def test_empty_dir_fails(self, ld_empty):
        is_valid, msg = validate_ld_ref("EUR", "1", ld_dir=ld_empty)
        assert is_valid is False
        assert "empty" in msg.lower()


# ── Missing snpinfo file ────────────────────────────────────────

class TestMissingSNPInfo:
    def test_no_snpinfo_fails(self, ld_no_snpinfo):
        is_valid, msg = validate_ld_ref("EUR", "1", ld_dir=ld_no_snpinfo)
        assert is_valid is False
        assert "snpinfo_mult_1kg_hm3" in msg

    def test_message_suggests_download(self, ld_no_snpinfo):
        _, msg = validate_ld_ref("EUR", "1", ld_dir=ld_no_snpinfo)
        assert "Download" in msg or "download" in msg


# ── Missing population directory ─────────────────────────────────

class TestMissingPopulation:
    def test_missing_pop_dir_fails(self, ld_no_pop_dir):
        is_valid, msg = validate_ld_ref("EUR", "1", ld_dir=ld_no_pop_dir)
        assert is_valid is False
        assert "EUR" in msg or "eur" in msg

    def test_afr_not_present(self, ld_complete):
        """EUR reference exists but AFR does not."""
        is_valid, msg = validate_ld_ref("AFR", "1", ld_dir=ld_complete)
        assert is_valid is False
        assert "AFR" in msg or "afr" in msg


# ── Missing chromosome file ─────────────────────────────────────

class TestMissingChromosome:
    def test_missing_chr_file_fails(self, ld_missing_chrom):
        # chr1 exists, chr2 does not
        is_valid, msg = validate_ld_ref("EUR", "2", ld_dir=ld_missing_chrom)
        assert is_valid is False
        assert "chr2" in msg

    def test_present_chr_still_passes(self, ld_missing_chrom):
        is_valid, _ = validate_ld_ref("EUR", "1", ld_dir=ld_missing_chrom)
        assert is_valid is True


# ── Default ld_dir (production path) ─────────────────────────────

class TestDefaultPath:
    def test_default_uses_app_path(self):
        """When ld_dir is not provided, the function should check /app/ld_reference.

        On a dev machine this path won't exist, so it should return False
        with a 'not found' message — NOT silently pass.
        """
        is_valid, msg = validate_ld_ref("EUR", "1")
        # /app/ld_reference almost certainly doesn't exist on the test host
        assert is_valid is False
        assert "not found" in msg or "ld_reference" in msg
