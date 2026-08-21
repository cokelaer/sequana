import numpy as np
import pytest

from sequana.metrics import (
    brukner_flexibility,
    compute_bendability,
    compute_curvature,
    compute_helix_twist,
    compute_sidd,
    dinucleotide_roll,
    helix_twist,
    lempel_ziv_complexity,
)


def test_lempel_ziv_complexity_repetitive():
    # Highly repetitive sequence should have lower complexity
    low = lempel_ziv_complexity("AAAAAAAAAA")
    high = lempel_ziv_complexity("ACGTACGTGT")
    assert low < high


def test_lempel_ziv_complexity_returns_float():
    result = lempel_ziv_complexity("ACGT")
    assert isinstance(result, float)
    assert result > 0


def test_compute_bendability():
    seq = "ACGTA"
    scores = compute_bendability(seq, brukner_flexibility, window=3)
    assert len(scores) == len(seq) - 3 + 1
    # All trinucleotides in seq should be in the scale
    assert all(s is not None for s in scores)


def test_compute_bendability_unknown_trinuc():
    # Use a scale that doesn't cover all trinucleotides
    scores = compute_bendability("NNN", {}, window=3)
    assert scores == [None]


def test_compute_helix_twist():
    seq = "ACGT"
    scores = compute_helix_twist(seq, helix_twist)
    assert len(scores) == len(seq) - 1
    assert all(s is not None for s in scores)


def test_compute_helix_twist_unknown_dinuc():
    scores = compute_helix_twist("NN", {})
    assert scores == [None]


def _roll_profile(seq):
    return [dinucleotide_roll.get(seq[i : i + 2], 0) for i in range(len(seq) - 1)]


def test_compute_curvature_flanks_are_nan():
    seq = "ACGTACGTACGTACGTACGT"
    curvature = compute_curvature(_roll_profile(seq), window=11)
    assert len(curvature) == len(seq) - 1
    # window//2 positions at each end have no value
    assert np.isnan(curvature[:5]).all()
    assert np.isnan(curvature[-5:]).all()
    assert not np.isnan(curvature[5:-5]).any()


def test_compute_curvature_phasing():
    # identical bends add up when they are in phase (no twist between steps) and
    # cancel out when successive bends point in opposite directions
    roll = np.ones(21)
    in_phase = compute_curvature(roll, twist=np.zeros(21), window=11)
    out_of_phase = compute_curvature(roll, twist=np.full(21, 180.0), window=11)
    assert in_phase[10] == pytest.approx(11)
    assert out_of_phase[10] == pytest.approx(1)


def test_compute_curvature_empty():
    assert len(compute_curvature([])) == 0


def test_compute_sidd_at_rich_is_destabilized():
    at_rich = "ATATATATATATATATATAT"
    gc_rich = "GCGCGCGCGCGCGCGCGCGC"
    seq = at_rich + gc_rich
    energies = compute_sidd(seq)
    assert len(energies) == len(seq)
    # the AT-rich half opens at a lower cost than the GC-rich one
    assert energies[: len(at_rich)].min() < energies[len(at_rich) :].min()


def test_compute_sidd_supercoiling_helps_opening():
    seq = "ATATATATATATATATATATGCGCGCGCGC"
    relaxed = compute_sidd(seq, sigma=0)
    supercoiled = compute_sidd(seq, sigma=-0.06)
    assert np.nanmin(supercoiled) < np.nanmin(relaxed)


def test_compute_sidd_too_short():
    assert np.isnan(compute_sidd("A")).all()
