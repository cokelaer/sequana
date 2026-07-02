import pandas as pd

import sequana.repeats.tandem as tandem_module
from sequana.repeats.tandem import ShortTandemRepeats

from . import test_dir

MEASLES = f"{test_dir}/data/fasta/measles.fa"


def _write_fasta(tmpdir, seq, name="planted.fa"):
    path = tmpdir.join(name)
    path.write(f">s\n{seq}\n")
    return str(path)


def test_str_columns_and_invariants():
    s = ShortTandemRepeats(MEASLES)
    s.run(progress=False)
    assert list(s.df.columns) == ["seqid", "start", "end", "period", "num", "remainder", "type", "sequence"]
    for row in s.df.itertuples(index=False):
        p = row.period
        assert s.min_period <= p <= s.max_period
        assert row.num >= s.min_reps
        assert row.end - row.start == row.num * p + row.remainder
        assert row.end - row.start >= s.min_span
        unit = row.sequence[:p]
        # every full copy equals the unit, and the remainder is its prefix
        for k in range(row.num):
            assert row.sequence[k * p : (k + 1) * p] == unit
        assert row.sequence[row.num * p :] == unit[: row.remainder]


def test_str_planted(tmpdir):
    fasta = _write_fasta(tmpdir, "N" * 4 + "CAG" * 5 + "N" * 4)
    s = ShortTandemRepeats(fasta)
    s.run(progress=False)
    assert len(s.df) == 1
    r = s.df.iloc[0]
    assert (r["start"], r["end"], r["period"], r["num"], r["remainder"]) == (4, 19, 3, 5, 0)


def test_str_partial_remainder(tmpdir):
    fasta = _write_fasta(tmpdir, "N" * 4 + "CAG" * 4 + "CA" + "N" * 4)
    s = ShortTandemRepeats(fasta)
    s.run(progress=False)
    assert len(s.df) == 1
    r = s.df.iloc[0]
    assert r["period"] == 3 and r["num"] == 4 and r["remainder"] == 2


def test_str_needs_min_reps_and_span(tmpdir):
    # only two copies, total 6 bp -> below min_reps (3) and min_span (10)
    fasta = _write_fasta(tmpdir, "N" * 4 + "CAG" * 2 + "N" * 4)
    s = ShortTandemRepeats(fasta)
    s.run(progress=False)
    assert s.df.empty


def test_str_type_code(tmpdir):
    # (AT)x6: even, complementary (A-T), strictly alternating R/Y, not a palindrome
    # code = 1(even) + 2(alt RY) + 0 + 8(complement) = 11
    fasta = _write_fasta(tmpdir, "G" * 4 + "AT" * 6 + "G" * 4)
    s = ShortTandemRepeats(fasta)
    s.run(progress=False)
    assert not s.df.empty
    assert s.df.iloc[0]["type"] == 11


def test_str_numba_fallback_parity():
    s1 = ShortTandemRepeats(MEASLES)
    s1.run(progress=False)
    orig = tandem_module._HAS_NUMBA
    try:
        tandem_module._HAS_NUMBA = False
        s2 = ShortTandemRepeats(MEASLES)
        s2.run(progress=False)
    finally:
        tandem_module._HAS_NUMBA = orig
    pd.testing.assert_frame_equal(s1.df, s2.df)


def test_str_to_gff_and_bed(tmpdir):
    fasta = _write_fasta(tmpdir, "N" * 4 + "CAG" * 5 + "N" * 4)
    s = ShortTandemRepeats(fasta)
    s.run(progress=False)
    gff = tmpdir.join("str.gff")
    s.to_gff(gff)
    lines = [l for l in open(gff) if not l.startswith("#")]
    assert lines
    fields = lines[0].split("\t")
    assert fields[2] == "Short_Tandem_Repeat"
    assert int(fields[3]) == int(s.df.iloc[0].start) + 1
    assert ";length=" in fields[8] and ";type=" in fields[8]
    bed = tmpdir.join("str.bed")
    s.to_bed(bed)
    assert bed.check()
