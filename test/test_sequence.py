import shutil

import numpy as np
import pylab
import pytest

from sequana.sequence import DNA, RNA, Repeats, Sequence, TrackSpec

from . import test_dir

datafile = f"{test_dir}/data/fasta/measles.fa"


def test_dna():
    seq = "ACGTTTT"
    dna = DNA(seq)
    dna.threshold
    assert dna.sequence == seq
    assert dna.get_complement() == "TGCAAAA"
    assert dna.get_reverse() == "TTTTGCA"
    assert dna.get_reverse_complement() == "AAAACGT"
    dna.check()
    dna.stats()

    # inplace functions
    dna.reverse()
    dna.complement()
    dna.reverse_complement()
    dna.window = 0.5
    dna.window
    dna.type_window
    dna.AT_skew
    dna.GC_skew

    dna = DNA("jjjj")
    try:
        dna.check()
        assert False
    except:
        assert True

    # read a file and tests the __len__ method
    dna = DNA(datafile)
    assert len(dna) == 15894
    dna.gc_content()
    # tests for the ORF and CDS functions
    assert dna.ORF_pos.shape == (1525, 6)
    dna.type_filter = "CDS"
    assert dna.ORF_pos.shape == (328, 6)
    dna.threshold = 20
    assert dna.ORF_pos.shape == (231, 6)
    dna.type_filter = "ORF"
    assert dna.ORF_pos.shape == (964, 6)
    dna.threshold = 0
    assert dna.ORF_pos.shape == (1525, 6)
    dna.hist_ORF_CDS_linearscale()
    dna.hist_ORF_CDS_logscale()

    dna.barplot_count_ORF_CDS_by_frame()

    # test occurences
    dna._data = "ACGTGGGGGTT"
    assert dna.get_occurences("GGG", False) == [4]
    assert dna.get_occurences("GGG", True) == [4, 5, 6]

    # property
    dna.type_filter
    dna.type_filter = "CDS"
    try:
        dna.type_filter = "dummy"
        assert False
    except ValueError:
        assert True

    assert dna.entropy("AAAAAAAAA") == 0


def test_rna():
    rna = RNA("ACUG")


@pytest.mark.external_tool
@pytest.mark.skipif(shutil.which("shustring") is None, reason="shustring not installed")
def test_repeats(tmpdir):
    datafile = f"{test_dir}/data/fasta/measles.fa"
    rep = Repeats(datafile)
    rep.threshold = 11
    # export to wig
    tmpfile = tmpdir.join("test.wig")
    rep.to_wig(tmpfile, step=100)
    assert len(rep.begin_end_repeat_position) == 158
    rep.do_merge = True
    assert len(rep.begin_end_repeat_position) == 156
    rep.do_merge = False
    assert len(rep.begin_end_repeat_position) == 158
    assert rep.df_shustring.shape == (15888, 2)
    assert rep.longest_shustring == 15

    # test histogram
    rep = Repeats(datafile)
    rep.threshold = 11
    rep.hist_length_repeats()
    rep.plot()
    assert rep.length == 15894
    rep.names
    rep.header

    # test several sequence identifiers, with comments, with same starting name
    # >tig001
    # >tig001 comment # with tabs
    # >tig001    comment # with spaces
    # >tig001ABC
    datafile = f"{test_dir}/data/fasta/test_shustring.fa"
    rep = Repeats(datafile)
    rep.threshold = 10


def test_gc_skew():

    dna = DNA(datafile)
    try:
        dna.window = -1
        assert False
    except ValueError:
        assert True

    dna.window = 100
    dna.plot_all_skews()
    for x in dna:  # test the iterator
        pass


def test_dna_plot_tracks():
    dna = DNA(datafile)
    dna.window = 100

    # default tracks are the skew/content ones
    fig = dna.plot_tracks()
    assert len(fig.axes) == len(DNA.DEFAULT_TRACKS)

    # a single track may be given as a string
    assert len(dna.plot_tracks("GC_skew").axes) == 1

    # tracks from both families (precomputed skews and get_* methods)
    assert len(dna.plot_tracks(["GC_skew", "curvature", "entropy"]).axes) == 3

    # every registered track must be plottable
    for name in DNA.get_available_tracks():
        assert len(dna.plot_tracks(name).axes) == 1
        pylab.close("all")


def test_dna_plot_tracks_window():
    dna = DNA(datafile)

    # a window is required by the skew tracks only
    with pytest.raises(AttributeError):
        dna.plot_tracks("GC_skew")
    dna.plot_tracks("entropy", window=500)

    # passing a window is equivalent to setting the window attribute
    dna.plot_tracks("GC_skew", window=200)
    assert dna.window == 200


def test_dna_plot_tracks_unknown():
    dna = DNA(datafile)
    dna.window = 100
    with pytest.raises(ValueError):
        dna.plot_tracks(["GC_skew", "dummy"])
    with pytest.raises(ValueError):
        dna.plot_tracks([])


def test_dna_plot_tracks_user_defined():
    dna = DNA("ACGTACGTACGTACGT")
    try:
        DNA.TRACKS["dummy"] = TrackSpec(
            getter=lambda dna, window: np.arange(len(dna.sequence)),
            ylabel="dummy",
            color="m",
        )
        assert "dummy" in DNA.get_available_tracks()
        assert len(dna.plot_tracks("dummy").axes) == 1
    finally:
        del DNA.TRACKS["dummy"]


def test_sequence_get_statistics():
    s = Sequence("AACGGTT")
    stats = s.get_statistics()
    assert "single" in stats
    df = stats["single"]
    assert df.loc["A", "count"] == 2
    assert df.loc["G", "count"] == 2
    assert df.loc["C", "count"] == 1
    assert df.loc["T", "count"] == 2
    assert df.loc["GC", "count"] == 3
    assert abs(df.loc["A", "percentage"] - 2 / 7 * 100) < 1e-6


def test_sequence_iter():
    # test_shustring.fa has 3 sequences; __init__ consumes the first,
    # so the iterator should yield the remaining 2
    multi_fa = f"{test_dir}/data/fasta/test_shustring.fa"
    dna = DNA(multi_fa)
    chunks = list(dna)
    assert len(chunks) == 2
    assert all(isinstance(c, str) for c in chunks)


def test_dna_get_dna_flexibility():
    # Use a short known sequence; result should be a numpy array of same length
    seq = "ACGTACGTACGT"
    dna = DNA(seq)
    flex = dna.get_dna_flexibility(window=4)
    assert len(flex) == len(seq)
    assert all(f > 0 for f in flex)


def test_dna_get_twist():
    dna = DNA("ACGTACGT")
    twist = dna.get_twist()
    assert len(twist) == 7
    assert twist[0] == 34.5  # AC
    assert twist[1] == 30.0  # CG
    # unknown letters fall back to the canonical twist
    assert all(x == pytest.approx(360 / 10.5) for x in DNA("NNNN").get_twist())


def test_dna_get_curvature():
    seq = "ACGTACGTACGTACGTACGTACGT"
    curvature = DNA(seq).get_curvature(window=11)
    assert len(curvature) == len(seq) - 1
    # window//2 steps on each side are undefined
    assert all(np.isnan(curvature[:5]))
    assert all(np.isnan(curvature[-5:]))
    assert all(v >= 0 for v in curvature[5:-5])


def test_dna_get_curvature_phasing():
    # bends repeated at the helical repeat add up; at half a turn they cancel
    phased = DNA(("AAAAAA" + "CCCCC") * 20).get_curvature(window=11)
    antiphased = DNA(("AAA" + "CC") * 44).get_curvature(window=11)
    assert np.nanmean(phased) > 5 * np.nanmean(antiphased)

    # a constant roll rotating along the helix cancels out over one turn
    assert np.nanmax(DNA("A" * 100).get_curvature(window=11)) < 1


def test_dna_get_entropy():
    # Uniform sequence -> low entropy; mixed -> higher
    dna_uniform = DNA("AAAAAAAAAA")
    dna_mixed = DNA("ACGTACGTAC")
    e_uniform = dna_uniform.get_entropy(window=4)
    e_mixed = dna_mixed.get_entropy(window=4)
    assert len(e_uniform) == 10
    assert e_uniform.mean() < e_mixed.mean()


def test_dna_get_informational_entropy():
    seq = "ACGTACGTACGTACGT"
    dna = DNA(seq)
    ie = dna.get_informational_entropy(window=8)
    assert len(ie) == len(seq)
    assert all(v >= 0 for v in ie)


def test_dna_get_dinucleotide_count():
    seq = "AACCGGTT"
    dna = DNA(seq)
    result = dna.get_dinucleotide_count(window=4)
    assert len(result) == len(seq)
    assert all(isinstance(v, int) for v in result)


def test_dna_get_trinucleotide_count():
    seq = "AAACCCGGGTTT"
    dna = DNA(seq)
    result = dna.get_trinucleotide_count(window=6)
    assert len(result) == len(seq)


def test_dna_get_homopolymers():
    seq = "AAAACCCGGGTTTT"
    dna = DNA(seq)
    result = dna.get_homopolymers(N=3, window=6)
    assert len(result) == len(seq)
    assert any(v > 0 for v in result)


def test_dna_get_karlin_signature_difference():
    seq = "ACGTACGTACGTACGTACGTACGT"
    dna = DNA(seq)
    result = dna.get_karlin_signature_difference(window=10)
    assert len(result) == len(seq)
    assert all(v >= 0 for v in result)
