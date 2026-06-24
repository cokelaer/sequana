import pytest

from sequana.wigtools import yield_bigwig_by_chromosome


def test_yield_bigwig_by_chromosome(tmp_path):
    pyBigWig = pytest.importorskip("pyBigWig")

    bw_file = tmp_path / "test.bw"
    bw = pyBigWig.open(str(bw_file), "w")
    bw.addHeader([("chr1", 1000), ("chr2", 1000)])
    bw.addEntries(
        ["chr1", "chr1", "chr1"],
        [0, 100, 200],
        ends=[100, 200, 300],
        values=[1.0, 2.0, 3.0],
    )
    bw.addEntries(
        ["chr2", "chr2"],
        [0, 100],
        ends=[100, 200],
        values=[4.0, 5.0],
    )
    bw.close()

    result = dict(yield_bigwig_by_chromosome(bw_file))
    assert result["chr1"] == [(0, 1.0), (100, 2.0), (200, 3.0)]
    assert result["chr2"] == [(0, 4.0), (100, 5.0)]
