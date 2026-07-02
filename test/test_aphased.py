from sequana.repeats.aphased import APhasedRepeats

from . import test_dir

MEASLES = f"{test_dir}/data/fasta/measles.fa"
SPACER = "GCGCGCG"  # 7 non-AT bases -> tract centres ~10 bp apart


def _write_fasta(tmpdir, seq, name="planted.fa"):
    path = tmpdir.join(name)
    path.write(f">s\n{seq}\n")
    return str(path)


def test_aphased_columns_and_invariants():
    a = APhasedRepeats(MEASLES)
    a.run(progress=False)
    assert list(a.df.columns) == ["seqid", "start", "end", "length", "tracts", "sequence"]
    for row in a.df.itertuples(index=False):
        assert row.end > row.start
        assert row.tracts >= a.min_tracts
        assert row.length == row.end - row.start
        # reported region starts and ends on an A/T base (a tract edge)
        assert row.sequence[0] in "AT" and row.sequence[-1] in "AT"


def test_aphased_planted(tmpdir):
    # three A-tracts one helical turn apart -> one A-phased repeat
    motif = "AAA" + SPACER + "AAA" + SPACER + "AAA"
    fasta = _write_fasta(tmpdir, "C" * 5 + motif + "C" * 5)
    a = APhasedRepeats(fasta)
    a.run(progress=False)
    assert len(a.df) == 1
    assert a.df.iloc[0].tracts == 3
    assert a.df.iloc[0].start == 5
    assert a.df.iloc[0].end == 5 + len(motif)


def test_aphased_needs_min_tracts(tmpdir):
    motif = "AAA" + SPACER + "AAA"  # only two tracts
    fasta = _write_fasta(tmpdir, "C" * 5 + motif + "C" * 5)
    a = APhasedRepeats(fasta)
    a.run(progress=False)
    assert a.df.empty


def test_aphased_requires_phasing(tmpdir):
    # tracts too far apart (~19 bp) -> not in phase -> nothing
    motif = "AAA" + "GC" * 8 + "AAA" + "GC" * 8 + "AAA"
    fasta = _write_fasta(tmpdir, "C" * 5 + motif + "C" * 5)
    a = APhasedRepeats(fasta)
    a.run(progress=False)
    assert a.df.empty


def test_aphased_reverse_strand_tract(tmpdir):
    # a tract valid only on the reverse-complement strand still counts
    # ATT (forward weak) == AAT on the rc strand (valid A-tract)
    motif = "AAA" + SPACER + "AAA" + SPACER + "ATT"
    fasta = _write_fasta(tmpdir, "C" * 5 + motif + "C" * 5)
    a = APhasedRepeats(fasta)
    a.run(progress=False)
    assert len(a.df) == 1
    assert a.df.iloc[0].tracts == 3


def test_aphased_to_gff_and_bed(tmpdir):
    motif = "AAA" + SPACER + "AAA" + SPACER + "AAA"
    fasta = _write_fasta(tmpdir, "C" * 5 + motif + "C" * 5)
    a = APhasedRepeats(fasta)
    a.run(progress=False)
    gff = tmpdir.join("apr.gff")
    a.to_gff(gff)
    lines = [l for l in open(gff) if not l.startswith("#")]
    assert lines
    fields = lines[0].split("\t")
    assert fields[2] == "A_Phased_Repeat"
    assert int(fields[3]) == int(a.df.iloc[0].start) + 1
    bed = tmpdir.join("apr.bed")
    a.to_bed(bed)
    assert bed.check()
