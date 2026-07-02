import pandas as pd

import sequana.repeats.hdna as hdna_module
from sequana.repeats.hdna import HDNA
from sequana.tools import reverse_complement

from . import test_dir

MEASLES = f"{test_dir}/data/fasta/measles.fa"


def _write_fasta(tmpdir, seq, name="planted.fa"):
    path = tmpdir.join(name)
    path.write(f">s\n{seq}\n")
    return str(path)


def _purine(base):
    return base in "AG"


def test_hdna_columns_and_invariants():
    h = HDNA(MEASLES)
    h.run(progress=False)
    assert list(h.df.columns) == ["seqid", "start", "end", "length", "repeat", "spacer", "subset", "sequence"]
    seq = open(MEASLES).read().split("\n", 1)[1].replace("\n", "").upper()
    for row in h.df.itertuples(index=False):
        rep, spacer = row.repeat, row.spacer
        # geometry
        assert row.end - row.start == 2 * rep + spacer
        assert rep >= h.min_repeat
        assert h.min_spacer <= spacer <= h.max_spacer
        # mirror repeat: right arm is the reverse (NOT complement) of the left arm
        left = seq[row.start : row.start + rep]
        right = seq[row.end - rep : row.end]
        assert left == right[::-1]
        # arm is homopurine or homopyrimidine within 10% impurities
        pur = sum(_purine(b) for b in left)
        impure = rep - max(pur, rep - pur)
        assert impure * 10 <= rep


def test_hdna_planted(tmpdir):
    arm = "AGGAAGGAGA"  # homopurine, length 10
    motif = arm + "AA" + arm[::-1]
    fasta = _write_fasta(tmpdir, "N" * 6 + motif + "N" * 6)
    h = HDNA(fasta)
    h.run(progress=False)
    assert not h.df.empty
    hit = h.df[(h.df.start == 6) & (h.df.end == 6 + len(motif))]
    assert not hit.empty
    assert hit["repeat"].max() >= len(arm)


def test_hdna_ignores_inverted_repeat(tmpdir):
    # a (non-pure) inverted repeat is not a triplex / mirror repeat
    arm = "ACGGTTAC"
    motif = arm + "AT" + reverse_complement(arm)
    fasta = _write_fasta(tmpdir, "N" * 6 + motif + "N" * 6)
    h = HDNA(fasta)
    h.run(progress=False)
    assert h.df.empty


def test_hdna_purity_tolerance(tmpdir):
    def detect(arm):
        motif = arm + "AA" + arm[::-1]
        fasta = _write_fasta(tmpdir, "N" * 6 + motif + "N" * 6, name=f"{arm}.fa")
        h = HDNA(fasta)
        h.run(progress=False)
        return not h.df.empty

    assert detect("AGGAAGGAGA")  # 0 impurity  -> detected
    assert detect("AGGAAGGAGC")  # 1/10 = 10%  -> detected
    assert not detect("AGGAACCAGA")  # 2/10 = 20% -> rejected


def test_hdna_n_does_not_extend(tmpdir):
    # N must not act as a mirror match; a pure-N stretch yields nothing
    fasta = _write_fasta(tmpdir, "N" * 60)
    h = HDNA(fasta)
    h.run(progress=False)
    assert h.df.empty


def test_hdna_numba_fallback_parity(tmpdir):
    arm = "AGGAAGGAGA"
    fasta = _write_fasta(tmpdir, "N" * 6 + arm + "AA" + arm[::-1] + "N" * 6)
    h1 = HDNA(fasta)
    h1.run(progress=False)
    orig = hdna_module._HAS_NUMBA
    try:
        hdna_module._HAS_NUMBA = False
        h2 = HDNA(fasta)
        h2.run(progress=False)
    finally:
        hdna_module._HAS_NUMBA = orig
    pd.testing.assert_frame_equal(h1.df, h2.df)


def test_hdna_to_gff(tmpdir):
    arm = "AGGAAGGAGA"
    fasta = _write_fasta(tmpdir, "N" * 6 + arm + "AA" + arm[::-1] + "N" * 6)
    h = HDNA(fasta)
    h.run(progress=False)
    out = tmpdir.join("tpx.gff")
    h.to_gff(out)
    lines = [l for l in open(out) if not l.startswith("#")]
    assert lines
    fields = lines[0].split("\t")
    assert fields[2] == "Triplex"
    assert int(fields[3]) == int(h.df.iloc[0].start) + 1
