import matplotlib

matplotlib.use("Agg")

import pytest
from matplotlib import colors
from matplotlib import pyplot as plt

from sequana.viz import upsetplot as upset

CONTENTS = {"A": ["g1", "g2", "g3"], "B": ["g2", "g3", "g4"], "C": ["g1", "g4"]}


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


def test_upset_plot():
    u = upset.UpSet(upset.from_contents(CONTENTS), sort_by="cardinality", show_counts=True)
    axes = u.plot()
    assert sorted(axes) == ["intersections", "matrix", "shading", "totals"]


def test_upset_plot_no_nan_colors():
    # UpSetPlot 0.9.0 sets its default colours with styles[col].fillna(x, inplace=True),
    # which is a no-op with the copy-on-write semantics of pandas >= 3. The colours
    # stayed NaN and matplotlib raised 'Invalid RGBA argument: nan'.
    u = upset.UpSet(upset.from_contents(CONTENTS), sort_by="cardinality")
    axes = u.plot()

    collections = [c for c in axes["matrix"].collections if len(c.get_facecolor())]
    assert collections
    for collection in collections:
        for attribute in ("get_facecolor", "get_edgecolor"):
            for rgba in getattr(collection, attribute)():
                assert colors.to_hex(rgba)  # raises on a NaN component


def test_upset_options():
    u = upset.UpSet(
        upset.from_contents(CONTENTS),
        sort_by="degree",
        sort_categories_by="cardinality",
        orientation="vertical",
        min_subset_size=1,
    )
    u.style_subsets(present="A", facecolor="red", edgecolor="black")
    u.plot()


def test_upset_reformat():
    assert len(upset.from_memberships([["A"], ["A", "B"], ["B"]], data=[1, 2, 3])) == 3
    assert upset.query(upset.from_contents(CONTENTS), sort_by="cardinality").data.shape == (4, 1)
