import numpy as np
import pytest

from sequana.minor_groove_width import MinorGrooveWidth, minor_groove_width
from sequana.tools import reverse_complement


@pytest.fixture
def mgw():
    return MinorGrooveWidth()


def test_pentamer_lookup(mgw):
    # exact value from the DNAshape query table
    assert mgw.pentamer_mgw("AAAAA") == 3.38
    assert mgw.pentamer_mgw("GCGCG") == 5.54


def test_reverse_complement_symmetry(mgw):
    # a pentamer and its reverse complement share the same MGW
    for pentamer in ["AAAAA", "GATTA", "ACGTA", "CCGGT"]:
        assert mgw.pentamer_mgw(pentamer) == mgw.pentamer_mgw(reverse_complement(pentamer))


def test_pentamer_invalid(mgw):
    assert np.isnan(mgw.pentamer_mgw("AAANA"))  # non-ACGT base
    assert np.isnan(mgw.pentamer_mgw("AAAA"))  # wrong length


def test_predict_flanks_are_nan(mgw):
    values = mgw.predict("GATTACAGATTACA")
    assert len(values) == 14
    assert np.isnan(values[0]) and np.isnan(values[1])
    assert np.isnan(values[-1]) and np.isnan(values[-2])
    # central positions are defined
    assert not np.isnan(values[2])


def test_predict_central_value(mgw):
    # single pentamer -> only the central base is defined
    values = mgw.predict("AAAAA")
    assert values[2] == 3.38
    assert np.isnan(values[[0, 1, 3, 4]]).all()


def test_predict_handles_N(mgw):
    values = mgw.predict("AAAAANAAAAA")
    # every window overlapping the N (5 of them) plus the 4 flanks are NaN
    assert int(np.isnan(values).sum()) == 9


def test_convenience_function():
    values = minor_groove_width("AAAAA")
    assert round(float(values[2]), 2) == 3.38
