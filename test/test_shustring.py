import shutil

import pytest

from sequana.repeats.shustring import Repeats

from . import test_dir

fasta_file = f"{test_dir}/data/fasta/measles.fa"

# shustring is an external bioconda binary; skip the whole module if unavailable.
pytestmark = pytest.mark.skipif(shutil.which("shustring") is None, reason="shustring binary not installed")


def test_repeats_basic():
    rr = Repeats(fasta_file)
    rr.threshold = 4
    assert rr.length > 0
    assert rr.longest_shustring > 0
    assert "repeat_length" in rr.df_shustring.columns
    assert isinstance(rr.begin_end_repeat_position, list)
    assert rr.list_len_repeats is not None
    assert rr.names[0] == rr.header


def test_repeats_merge():
    rr = Repeats(fasta_file, merge=True)
    rr.threshold = 4
    merged = rr.begin_end_repeat_position
    # toggling do_merge off recomputes the (un-merged) positions
    rr.do_merge = False
    assert isinstance(merged, list)


def test_repeats_threshold_negative_raises():
    rr = Repeats(fasta_file)
    with pytest.raises(ValueError):
        rr.threshold = -1


def test_repeats_bad_name_raises():
    with pytest.raises(ValueError):
        Repeats(fasta_file, name="does_not_exist")


def test_repeats_to_wig(tmp_path):
    rr = Repeats(fasta_file)
    rr.threshold = 4
    out = tmp_path / "repeats.wig"
    rr.to_wig(str(out), step=100)
    assert out.exists() and out.stat().st_size > 0


def test_repeats_plots():
    rr = Repeats(fasta_file)
    rr.threshold = 4
    rr.plot()
    rr.hist_length_repeats()
