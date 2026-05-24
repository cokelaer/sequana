from sequana.summary import Summary

from . import test_dir


def test_summary():
    s = Summary("test2", sample_name="chr1", data={"mean": 1})
    assert s.data == {"mean": 1}
    assert s.version
    assert s.date
    d = s.as_dict()
    assert "name" in d
    assert "version" in d
    assert "data" in d
    assert "date" in d

    # test wrong constructor
    try:
        s = Summary("test")
        assert False
    except:
        assert True

    try:
        s = Summary("test", "test")
        assert False
    except:
        assert True

    # test data_description
    s = Summary("test2", data={"mean": 1})
    s.data_description = {"mean": "mean of the data set"}
    assert s.data_description == {"mean": "mean of the data set"}
    try:
        s.data_description = {"dummy": 1}
        assert False
    except:
        assert True

    from easydev import TempFile

    with TempFile(suffix=".json") as fh:
        s.to_json(fh.name)


def test_summary_from_json():

    s = Summary(f"{test_dir}/modules_report/data/test_summary_fastq_stats.json")


def test_summary_cli_json_export():
    """Test JSON export from sequana summary CLI"""
    import json
    import tempfile
    from pathlib import Path

    from click.testing import CliRunner

    from sequana.scripts.main.summary import summary

    runner = CliRunner()

    # Test FASTQ JSON export
    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as fh:
        json_file = fh.name

    try:
        result = runner.invoke(summary, [f"{test_dir}/data/reads.fastq", "--output-json", json_file])
        assert result.exit_code == 0, f"CLI failed: {result.output}"

        # Verify JSON file created
        assert Path(json_file).exists()

        # Load and validate JSON structure
        with open(json_file) as f:
            data = json.load(f)

        assert isinstance(data, list)
        assert len(data) == 1
        assert "filename" in data[0]
        assert "format" in data[0]
        assert "timestamp" in data[0]
        assert "stats" in data[0]
        assert data[0]["format"] == "fastq"
        assert isinstance(data[0]["stats"], dict)
        assert "n_reads" in data[0]["stats"]
        assert "GC content" in data[0]["stats"]
    finally:
        Path(json_file).unlink(missing_ok=True)

    # Test FASTA JSON export
    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as fh:
        json_file = fh.name

    try:
        result = runner.invoke(
            summary,
            [f"{test_dir}/data/test_fasta_filtering.fasta", "--output-json", json_file],
        )
        assert result.exit_code == 0, f"CLI failed: {result.output}"

        with open(json_file) as f:
            data = json.load(f)

        assert len(data) == 1
        assert data[0]["format"] == "fasta"
        assert data[0]["stats"]["N"] == 4
        assert "N50" in data[0]["stats"]
    finally:
        Path(json_file).unlink(missing_ok=True)
