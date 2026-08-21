from sequana.multiqc.plots import STAR, Bowtie1Reader, Bowtie2, FeatureCounts

from . import test_dir


def test_feature_count():

    fc = FeatureCounts(f"{test_dir}/data/multiqc_featureCounts_star.txt")
    fc.plot(html_code=True)

    fc = FeatureCounts(f"{test_dir}/data/multiqc_featureCounts_bowtie2.txt")
    fc.plot(html_code=True)


def test_feature_count_assignment_plot_file():
    # this multiqc file uses column names such as 'Unassigned: No Features' instead of
    # 'Unassigned_NoFeatures' and has no multi-mapping column
    fc = FeatureCounts(f"{test_dir}/data/featureCounts_assignment_plot.txt")
    assert "Unassigned_NoFeatures" in fc.df.columns
    fig = fc.plot(html_code=True)
    assert len(fig.data) == 4


def test_bowtie1():
    b = Bowtie1Reader(f"{test_dir}/data/multiqc_bowtie1.txt")
    b.plot_bar(html_code=True)


def test_bowtie2():
    b = Bowtie2(f"{test_dir}/data/multiqc_bowtie2_unpaired.txt")
    b.plot(html_code=True)

    b = Bowtie2(f"{test_dir}/data/multiqc_bowtie2_paired.txt")
    b.plot(html_code=True)


def test_bowtie2_plot_files():
    # multiqc also exports the data of the plots themselves. Those files use the
    # legend names (e.g. 'SE mapped uniquely') rather than the keys found in
    # multiqc_bowtie2.txt. Both must be understood.
    b = Bowtie2(f"{test_dir}/data/mqc_bowtie2_se_plot_1.txt")
    fig = b.plot(html_code=True)
    assert len(fig.data) == 3
    # all bars must sum up to 100%
    assert abs(sum(trace.x[0] for trace in fig.data) - 100) < 1e-6

    b = Bowtie2(f"{test_dir}/data/mqc_bowtie2_pe_plot_1.txt")
    fig = b.plot(html_code=True)
    assert len(fig.data) == 6
    assert abs(sum(trace.x[0] for trace in fig.data) - 100) < 1e-6


def test_bowtie2_no_alignment_column():
    import pytest

    with pytest.raises(ValueError):
        Bowtie2(f"{test_dir}/data/multiqc_star.txt").plot(html_code=True)


def test_star():
    b = STAR(f"{test_dir}/data/multiqc_star.txt")
    b.plot(html_code=True)
