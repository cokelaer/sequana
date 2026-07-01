import glob

from sequana.iem import BCLConvert
from sequana.iem import SampleSheet
from sequana.iem import SampleSheet as IEM
from sequana.iem import SampleSheetFactory, get_sample_sheet_version

from . import test_dir


def test_iem():

    try:
        iem = IEM("dummy")
        assert False
    except:
        assert True

    for this in [
        "iem/wrong/test_expdesign_wrong.csv",
        "iem/good/test_expdesign_miseq_illumina_1.csv",
        "iem/good/test_expdesign_miseq_illumina2.csv",
    ]:

        filename = f"{test_dir}/data/{this}"

        iem = IEM(filename)
        iem.settings
        # iem.name
        iem.samples
        iem.index_adapters
        iem.header
        iem.to_fasta("TEST")
        try:
            iem.validate()
        except SystemExit:
            pass
        iem.checker()
        print(iem)


def test_quick_fix(tmpdir):

    fout = tmpdir.join("fix.csv")
    e = IEM(f"{test_dir}/data/iem/wrong/test_expdesign_miseq_illumina_semicommas.csv")
    e.quick_fix(fout)
    e = IEM(fout)
    e.validate()


def test_warning():

    for filename in glob.glob(f"{test_dir}/data/iem/warning/*csv"):
        print(filename)
        e = IEM(filename)
        e.validate()
        checks = e.checker()
        errors = [1 for check in checks if check["status"] == "Error"]
        assert len(errors) == 0
        warnings = [1 for check in checks if check["status"] == "Warning"]
        assert len(warnings) > 0


def test_wrong():

    for filename in glob.glob(f"{test_dir}/data/iem/wrong/*csv"):

        e = IEM(f"{filename}")

        try:
            e.validate()
            assert False, filename
        except SystemExit:
            checks = e.checker()
            errors = [1 for check in checks if check["status"] == "Error"]
            if sum(errors) == 0:
                assert False, filename


def test_iem_samplesheets():
    for filename in glob.glob(f"{test_dir}/data/iem/good/*csv"):
        print(filename)
        e = IEM(f"{filename}")
        e.validate()
        e.version
        e.instrument
        e.header
        checks = e.checker()
        for check in checks:
            if check["status"] == "Error":
                assert False


def test_bclconvert_v2():
    filename = f"{test_dir}/data/iem/test_samplesheet_v2_bclconvert.csv"

    # version detection and factory routing
    assert get_sample_sheet_version(filename) == "v2"
    ss = SampleSheetFactory(filename)
    assert isinstance(ss, BCLConvert)

    # v2 sections parsed and pointed to by the class attributes
    assert ss.data_section in ss.sections
    assert ss.settings_section in ss.sections
    assert list(ss.df.columns) == ["sample_id", "index", "index2"]
    assert "AdapterRead1" in ss.settings

    # this example is valid: no error expected
    ss.validate()
    checks = ss.checker()
    assert [c for c in checks if c["status"] == "Error"] == []


def test_bclconvert_cloud_data_malformed():
    # A missing comma in the optional [Cloud_Data] section must be flagged as a
    # Warning (not Error): BCL Convert ignores that section for demux.
    filename = f"{test_dir}/data/iem/test_samplesheet_v2_bclconvert_cloud_data_bad.csv"
    ss = SampleSheetFactory(filename)
    assert isinstance(ss, BCLConvert)

    checks = ss.checker()

    # the mandatory [BCLConvert_Data] section is still fine: no error
    assert [c for c in checks if c["status"] == "Error"] == []

    cloud = [c for c in checks if c.get("name") == "check_cloud_data_section_csv_format"]
    assert len(cloud) == 1
    assert cloud[0]["status"] == "Warning"


def test_bclconvert_quick_fix_preserves_semicolons(tmpdir):
    filename = f"{test_dir}/data/iem/test_samplesheet_v2_bclconvert.csv"
    fout = tmpdir.join("fixed.csv")
    BCLConvert(filename).quick_fix(str(fout))
    text = fout.read()
    # internal semicolons in OverrideCycles must survive quick_fix
    assert "R1:Y151;I1:I10;I2:I10;R2:Y151" in text


def test_get_sample_sheet_version_v1():
    filename = f"{test_dir}/data/iem/good/test_expdesign_miseq_illumina_1.csv"
    assert get_sample_sheet_version(filename) == "v1"
    assert isinstance(SampleSheetFactory(filename), SampleSheet)
