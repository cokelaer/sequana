import pandas as pd
import pytest

import sequana.repeats.directrepeat as dr_module
from sequana.repeats.directrepeat import DirectRepeats

from . import test_dir

MEASLES = f"{test_dir}/data/fasta/measles.fa"


def _write_fasta(tmpdir, seq, name="planted.fa"):
    path = tmpdir.join(name)
    path.write(f">s\n{seq}\n")
    return str(path)


def test_dr_columns_and_invariants():
    d = DirectRepeats(MEASLES)
    d.run(progress=False)
    assert list(d.df.columns) == ["seqid", "start", "end", "repeat", "spacer", "num", "remainder", "subset", "sequence"]
    for row in d.df.itertuples(index=False):
        assert d.min_repeat <= row.repeat <= d.max_repeat
        assert 0 <= row.spacer <= d.max_spacer
        assert row.end - row.start == len(row.sequence)
        assert row.subset == int(row.spacer == 0)
        # the two arms are identical: seq[:repeat] == seq[repeat+spacer : 2*repeat+spacer]
        arm = row.sequence[: row.repeat]
        second = row.sequence[row.repeat + row.spacer : 2 * row.repeat + row.spacer]
        assert arm == second


def test_dr_planted_no_spacer(tmpdir):
    arm = "ACGTGACTGA"  # 10 bp, no internal period
    fasta = _write_fasta(tmpdir, "N" * 4 + arm + arm + "N" * 4)
    d = DirectRepeats(fasta)
    d.run(progress=False)
    assert len(d.df) == 1
    r = d.df.iloc[0]
    assert (r["start"], r["end"], r["repeat"], r["spacer"]) == (4, 24, 10, 0)
    assert (r["num"], r["remainder"], r["subset"]) == (1, 0, 1)
    assert r["sequence"] == arm + arm


def test_dr_planted_with_spacer(tmpdir):
    arm = "ACGTGACTGA"
    fasta = _write_fasta(tmpdir, "N" * 4 + arm + "TT" + arm + "N" * 4)
    d = DirectRepeats(fasta)
    d.run(progress=False)
    assert len(d.df) == 1
    r = d.df.iloc[0]
    assert (r["repeat"], r["spacer"], r["subset"]) == (10, 2, 0)
    assert r["start"] == 4 and r["end"] == 26


def test_dr_btr_expansion(tmpdir):
    # a tandem AT block expands beyond two arms (BTR), giving num > 1 or remainder
    fasta = _write_fasta(tmpdir, "N" * 4 + "AT" * 15 + "N" * 4)
    d = DirectRepeats(fasta)
    d.run(progress=False)
    assert len(d.df) == 1
    r = d.df.iloc[0]
    assert r["spacer"] == 0 and r["subset"] == 1
    # the recorded total length is repeat*num + remainder
    assert r["end"] - r["start"] >= 2 * r["repeat"]


def test_dr_too_short_no_call(tmpdir):
    # arms of 9 bp are below min_repeat (10) -> nothing
    arm = "ACGTGACTG"  # 9 bp
    fasta = _write_fasta(tmpdir, "N" * 4 + arm + arm + "N" * 4)
    d = DirectRepeats(fasta)
    d.run(progress=False)
    assert d.df.empty


def test_dr_numba_fallback_parity():
    d1 = DirectRepeats(MEASLES)
    d1.run(progress=False)
    orig = dr_module._HAS_NUMBA
    try:
        dr_module._HAS_NUMBA = False
        d2 = DirectRepeats(MEASLES)
        d2.run(progress=False)
    finally:
        dr_module._HAS_NUMBA = orig
    pd.testing.assert_frame_equal(d1.df, d2.df)


def test_dr_to_gff_and_bed(tmpdir):
    arm = "ACGTGACTGA"
    fasta = _write_fasta(tmpdir, "N" * 4 + arm + arm + "N" * 4)
    d = DirectRepeats(fasta)
    d.run(progress=False)
    gff = tmpdir.join("dr.gff")
    d.to_gff(gff)
    lines = [l for l in open(gff) if not l.startswith("#")]
    assert lines
    fields = lines[0].split("\t")
    assert fields[2] == "Direct_Repeat"
    assert int(fields[3]) == int(d.df.iloc[0].start) + 1
    assert ";repeat=" in fields[8] and ";spacer=" in fields[8] and ";subset=" in fields[8]
    # composition counts only the arm (repeat bases)
    bed = tmpdir.join("dr.bed")
    d.to_bed(bed)
    assert bed.check()


def test_dr_to_bed_empty_raises(tmpdir):
    d = DirectRepeats(_write_fasta(tmpdir, "N" * 50))
    d.run(progress=False)
    with pytest.raises(ValueError):
        d.to_bed(str(tmpdir.join("empty.bed")))
