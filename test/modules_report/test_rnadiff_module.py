import json
import os
import re
import shutil

import matplotlib

matplotlib.use("Agg")

import pytest

from sequana.modules_report.rnadiff import RNAdiffModule
from sequana.utils import config

from . import test_dir

RNADIFF_DIR = f"{test_dir}/../data/rnadiff/rnadiff_0.15.4"


@pytest.fixture
def report(tmpdir):
    """Build the RNA-diff HTML report in a temporary directory"""
    folder = str(tmpdir.join("rnadiff"))
    shutil.copytree(RNADIFF_DIR, folder)
    os.makedirs(f"{folder}/images", exist_ok=True)

    config.output_dir = str(tmpdir)
    RNAdiffModule(folder, gff=None, output_filename="summary.html")
    return open(str(tmpdir.join("summary.html"))).read()


def test_rnadiff_module_anchors(report):
    # every entry of the table of contents must point to one and only one
    # division. The sections used to share anchors ('table', 'filters_option'),
    # so all their links landed on the same division.
    links = re.findall(r'<a href="#([^"]+)">([^<]*)</a>', report)
    ids = re.findall(r'<div id="([^"]+)"', report)

    assert links
    assert len(ids) == len(set(ids)), "duplicated division identifiers"
    for anchor, name in links:
        assert ids.count(anchor) == 1, f"the '{name}' link points to {ids.count(anchor)} divisions"


def test_rnadiff_module_sections(report):
    for name in ("Diagnostic plots", "Clusterisation", "Normalisation", "Dispersion", "DGE results"):
        assert name in report


@pytest.fixture
def report_two_comparisons(tmpdir):
    """Same report, with a second comparison so that the venn/upset are included"""
    folder = str(tmpdir.join("rnadiff"))
    shutil.copytree(RNADIFF_DIR, folder)
    shutil.copy(f"{folder}/A_vs_B_degs_DESeq2.csv", f"{folder}/A_vs_C_degs_DESeq2.csv")
    os.makedirs(f"{folder}/images", exist_ok=True)

    config.output_dir = str(tmpdir)
    RNAdiffModule(folder, gff=None, output_filename="summary.html")
    return open(str(tmpdir.join("summary.html"))).read()


def test_rnadiff_module_venn(report_two_comparisons):
    report = report_two_comparisons

    assert "Venn diagram" in report
    assert 'id="rnadiff_venn_plot"' in report

    # the data used by the interactive diagram
    data = json.loads(re.search(r"var D = (\{.*?\});\n", report, re.DOTALL).group(1))
    assert sorted(data["names"]) == ["A_vs_B", "A_vs_C"]
    assert sorted(data["layouts"]) == ["2"]

    # the two comparisons are identical, hence a gene is either in both lists
    # or in none of them
    assert all(bool(row[0]) == bool(row[1]) for row in data["rows"])
