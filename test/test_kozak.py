import pytest

from sequana import GFF3, FastA
from sequana.kozak import ConsensusBuilder, Kozak
from sequana.lazy import pandas as pd

from . import test_dir


@pytest.fixture
def pfm():
    # position-frequency matrix used across ConsensusBuilder tests
    return pd.DataFrame(
        {
            "A": [0.90, 0.45, 0.30, 0.55, 0.10],
            "C": [0.04, 0.40, 0.30, 0.20, 0.10],
            "G": [0.03, 0.10, 0.30, 0.15, 0.70],
            "T": [0.03, 0.05, 0.10, 0.10, 0.10],
        }
    )


def test_consensus_builder_modes(pfm):
    cb = ConsensusBuilder(pfm)
    assert cb.get_consensus(mode="max_only") == "AAAAG"
    assert cb.get_consensus(mode="majority") == "AnnAG"
    assert cb.get_consensus(mode="threshold") == "AmvaG"
    assert cb.get_consensus(mode="relative") == "AmvaG"
    assert cb.get_consensus(mode="information") == "Amvag"


def test_consensus_builder_unknown_mode(pfm):
    with pytest.raises(ValueError):
        ConsensusBuilder(pfm).get_consensus(mode="does-not-exist")


def test_consensus_builder_cavener(pfm):
    cb = ConsensusBuilder(pfm)
    assert cb.get_consensus_cavener() == "AMaAG"
    assert cb.get_consensus_cavener(notation="expanded") == "AA/CaAG"


def test_consensus_builder_notation(pfm):
    cb = ConsensusBuilder(pfm)
    # iupac is the default and uses single ambiguity codes
    assert cb.get_consensus(mode="threshold", notation="iupac") == "AmvaG"
    # expanded uses slash-joined bases for ambiguous positions
    assert cb.get_consensus(mode="threshold", notation="expanded") == "Aa/ca/c/gaG"


def test_consensus_builder_iupac():
    cb = ConsensusBuilder(pd.DataFrame())
    # order independent, and the unsorted 'GC' key still resolves to S
    assert cb._iupac(["C", "G"]) == "S"
    assert cb._iupac(["G", "C"]) == "S"
    assert cb._iupac(["A"]) == "A"
    assert cb._iupac(["A", "C", "G", "T"]) == "N"
    assert cb._iupac(["C", "G"], notation="expanded") == "C/G"
    with pytest.raises(ValueError):
        cb._iupac(["A", "C"], notation="bad")


def test_consensus_builder_all_consensus(pfm, capsys):
    ConsensusBuilder(pfm).all_consensus()
    out = capsys.readouterr().out
    assert "cavener: AMaAG" in out
    assert "majority: AnnAG" in out


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
