import pandas as pd
import pytest

import sequana.repeats.zdna as zdna_module
from sequana.repeats.zdna import ZDNA

from . import test_dir

MEASLES = f"{test_dir}/data/fasta/measles.fa"

# valid alternating purine/pyrimidine dinucleotides (gfa pupy, AT/TA excluded)
_VALID = {"AC", "TG", "CG", "CA", "GC", "GT"}


def _write_fasta(tmpdir, seq, name="planted.fa"):
    path = tmpdir.join(name)
    path.write(f">s\n{seq}\n")
    return str(path)


def test_zdna_columns_and_invariants():
    z = ZDNA(MEASLES)
    z.run(progress=False)
    assert list(z.df.columns) == ["seqid", "start", "end", "length", "score", "subset", "sequence"]
    for row in z.df.itertuples(index=False):
        assert row.length >= z.min_z
        assert row.end - row.start == row.length
        assert len(row.sequence) == row.length
        assert row.subset == int(row.score >= z.min_kvscore)
        # every consecutive dinucleotide of the motif is a valid Z-DNA dinuc
        s = row.sequence.upper()
        for k in range(len(s) - 1):
            assert s[k : k + 2] in _VALID


def test_zdna_planted_strong(tmpdir):
    # GC/CG only: 10 bp -> 9 dinucleotides x25 -> score = 225//2 = 112, subset=1
    fasta = _write_fasta(tmpdir, "N" * 4 + "GCGCGCGCGC" + "N" * 20)
    z = ZDNA(fasta)
    z.run(progress=False)
    assert len(z.df) == 1
    r = z.df.iloc[0]
    assert (r["start"], r["end"], r["length"]) == (4, 14, 10)
    assert r["score"] == 112
    assert r["subset"] == 1
    assert r["sequence"] == "GCGCGCGCGC"


def test_zdna_weak_low_score(tmpdir):
    # AC/CA only: 9 dinucleotides x3 -> score = 27//2 = 13 -> subset=0
    fasta = _write_fasta(tmpdir, "N" * 4 + "ACACACACAC" + "N" * 20)
    z = ZDNA(fasta)
    z.run(progress=False)
    assert len(z.df) == 1
    r = z.df.iloc[0]
    assert r["length"] == 10 and r["score"] == 13 and r["subset"] == 0


def test_zdna_excludes_at_ta(tmpdir):
    # AT/TA alternation is not Z-forming -> no motif
    fasta = _write_fasta(tmpdir, "N" * 4 + "ATATATATATAT" + "N" * 20)
    z = ZDNA(fasta)
    z.run(progress=False)
    assert z.df.empty


def test_zdna_needs_min_length(tmpdir):
    # 9 bp run -> npy = 9 < min_z (10) -> nothing
    fasta = _write_fasta(tmpdir, "N" * 4 + "GCGCGCGCG" + "N" * 20)
    z = ZDNA(fasta)
    z.run(progress=False)
    assert z.df.empty


def test_zdna_numba_fallback_parity():
    z1 = ZDNA(MEASLES)
    z1.run(progress=False)
    orig = zdna_module._HAS_NUMBA
    try:
        zdna_module._HAS_NUMBA = False
        z2 = ZDNA(MEASLES)
        z2.run(progress=False)
    finally:
        zdna_module._HAS_NUMBA = orig
    pd.testing.assert_frame_equal(z1.df, z2.df)


def test_zdna_to_gff_and_bed(tmpdir):
    fasta = _write_fasta(tmpdir, "N" * 4 + "GCGCGCGCGC" + "N" * 20)
    z = ZDNA(fasta)
    z.run(progress=False)
    gff = tmpdir.join("z.gff")
    z.to_gff(gff)
    lines = [l for l in open(gff) if not l.startswith("#")]
    assert lines
    fields = lines[0].split("\t")
    assert fields[2] == "Z_DNA_Motif"
    assert int(fields[3]) == int(z.df.iloc[0].start) + 1
    assert ";score=" in fields[8] and ";subset=" in fields[8] and ";composition=" in fields[8]
    bed = tmpdir.join("z.bed")
    z.to_bed(bed)
    assert bed.check()


def test_zdna_to_bed_empty_raises(tmpdir):
    z = ZDNA(_write_fasta(tmpdir, "N" * 50))
    z.run(progress=False)
    with pytest.raises(ValueError):
        z.to_bed(str(tmpdir.join("empty.bed")))
