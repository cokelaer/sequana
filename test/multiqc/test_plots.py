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


def test_star():
    b = STAR(f"{test_dir}/data/multiqc_star.txt")
    b.plot(html_code=True)
