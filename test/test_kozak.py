from sequana import GFF3, FastA
from sequana.kozak import Kozak
from sequana.lazy import pandas as pd

from . import test_dir


def test_kozak():

    fastafile = f"{test_dir}/data/fasta/ecoli_MG1655.fa"
    gff = f"{test_dir}/data/gff/ecoli_MG1655.gff"

    k = Kozak(fastafile, gff, "gene", "ID")
    k.set_context(left_kozak=15, right_kozak=12, keep_ATG_only=True, include_start_codon=True)

    df = k.get_data()

    assert len(df["kozak_left"].iloc[0]) == 15
    assert list(df["start_codon"].unique()) == ["ATG"]
    assert k.metrics["ATG_ratio"] != 1  # original data has non-ATG start codon

    # check caching is functional (correctly reset)
    k.set_context(left_kozak=12, right_kozak=12, keep_ATG_only=False, include_start_codon=True)
    df = k.get_data()
    assert len(df["kozak_left"].iloc[0]) == 12
    assert k.metrics["ATG_ratio"] != 1
    assert list(df["start_codon"].unique()) != ["ATG"]

    dd = k.plot_logo(df)
    k.plot_logo_bits(df)
    k.plot_logo_purine_pyrimidine(df)
    k.get_gc_per_chromosome()

    assert k._get_logo_data(df).equals(k.plot_logo())

    motif = k._get_logo_data(df)
    assert k._add_purine_pyrimidine(motif).equals(k.plot_logo_purine_pyrimidine())

    k.get_information_content(motif)


def test_kozak_uniform():
    fastafile = f"{test_dir}/data/fasta/ecoli_MG1655.fa"
    gff = f"{test_dir}/data/gff/ecoli_MG1655.gff"

    k = Kozak(fastafile, gff, "gene", "ID")
    k.set_context(background_method="uniform")

    # Trigger background generation
    bg = k.get_random_contexts()
    assert len(bg) == len(k.get_data())
    assert "kozak_left" in bg.columns
    assert "start_codon" in bg.columns
    assert "kozak_right" in bg.columns

    # Check KL divergence with uniform background
    kl = k._get_KL_data()
    assert len(kl) == k.left_kozak + k.right_kozak


def test_kozak_compare_by_chromosome():
    fastafile = f"{test_dir}/data/fasta/ecoli_MG1655.fa"
    gff = f"{test_dir}/data/gff/ecoli_MG1655.gff"

    k = Kozak(fastafile, gff, "gene", "ID")
    k.set_context(left_kozak=9, right_kozak=6)

    df = k.get_data()

    # Single chromosome case: returns empty DataFrame with warning
    result = k.compare_motif_by_chromosome(df)
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 0

    summary, residuals = k.compare_motif_by_chromosome_detailed(df)
    assert isinstance(summary, pd.DataFrame)
    assert isinstance(residuals, pd.DataFrame)
    assert len(summary) == 0

    # Multi-chromosome case: simulate by splitting rows into fake chromosomes
    df_multi = df.copy()
    mid = len(df_multi) // 2
    df_multi.loc[df_multi.index[:mid], "chrom"] = "chromA"
    df_multi.loc[df_multi.index[mid:], "chrom"] = "chromB"

    result = k.compare_motif_by_chromosome(df_multi)
    assert len(result) > 0
    assert set(result.columns) == {"position", "chi2", "p_value", "dof", "significant"}
    assert (result["p_value"] >= 0).all() and (result["p_value"] <= 1).all()
    assert (result["chi2"] >= 0).all()

    summary, residuals = k.compare_motif_by_chromosome_detailed(df_multi)
    assert len(summary) > 0
    assert "most_extreme_chrom" in summary.columns
    assert "max_residual" in summary.columns
    assert set(summary["most_extreme_chrom"].unique()).issubset({"chromA", "chromB"})
    # residuals_df: rows = chromosomes, cols = positions
    assert set(residuals.index) == {"chromA", "chromB"}
    assert residuals.shape[1] == len(summary)
    assert (residuals.values >= 0).all()  # effect sizes are absolute

    # Dual-panel plot method
    import matplotlib

    matplotlib.use("Agg")
    summary2, residuals2, fig = k.plot_chromosome_divergence(df_multi)
    assert fig is not None
    assert len(fig.axes) >= 2  # top + bottom + maybe colorbar
