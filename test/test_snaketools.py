from sequana import snaketools, sequana_data
from sequana.snaketools import DOTParser
import os
from nose.plugins.attrib import attr
from sequana import Module, SequanaConfig


testing = lambda x: sequana_data(x, "testing")

def test_dot_parser():
    s = DOTParser(testing("test_dag.dot"))
    s.add_urls()
    try:os.remove("test_dag.ann.dot")
    except:pass


def test_modules():
    assert "dag" in snaketools.modules.keys()
    assert snaketools.modules['dag'].endswith("dag.rules")


def test_getcleanup_rules():
    filename =  snaketools.modules['fastq_sampling']
    try:
        snaketools.get_cleanup_rules(filename)
    except:
        pass


def test_snakemake_stats():
    # this is created using snakemake with the option "--stats stats.txt"
    s = snaketools.SnakeMakeStats(testing("test_snakemake_stats.txt"))
    s.plot()


def test_module():
    # a rule without README
    m = snaketools.Module('mark_duplicates')
    m.description
    print(m)
    m.path
    m.snakefile

    # a rule with README
    m = snaketools.Module('dag')
    m.description
    m.overview
    m.is_executable()
    m.check()

    # a pipeline
    m = snaketools.Module('quality_control')
    m.is_executable()
    m.check()
    m.snakefile
    m.name


@attr("onweb")
def test_module_onweb():
    m = snaketools.Module('quality_control')
    m.onweb()


def test_valid_config():
    config = snaketools.SequanaConfig(None)

    s = snaketools.Module("quality_control")
    config = snaketools.SequanaConfig(s.config)

    from easydev import TempFile
    with TempFile() as fh:
        config.save(fh.name)


def test_sequana_config():
    s = snaketools.Module("quality_control")
    config = snaketools.SequanaConfig(s.config)

    assert config.config.get("kraken:dummy", "test") == "test"
    assert config.config.get("kraken:dummy") == None


    # --------------------------------- tests different constructors
    config = snaketools.SequanaConfig()
    config = snaketools.SequanaConfig({"test":1}, mode="others")
    assert config.config.test == 1
    # with a dictionary
    config = snaketools.SequanaConfig(config.config, mode="others")
    # with a sequanaConfig instance
    config = snaketools.SequanaConfig(config, mode="others")
    # with a non-yaml file
    try:
        json = sequana_data('test_summary_fastq_stats.json')
        config = snaketools.SequanaConfig(json, mode="others")
        assert False
    except:
        assert True
    try:
        config = snaketools.SequanaConfig("dummy_dummy")
        assert False
    except:
        assert True


    # loop over all pipelines, read the config, save it and check the content is
    # identical. This requires to remove the templates. We want to make sure the
    # empty strings are kept and that "no value" are kept as well
    #
    #    field1: ""
    #    field2:
    #
    # is unchanged
    for pipeline in snaketools.pipeline_names:
        config_filename = Module(pipeline)._get_config()
        cfg1 = SequanaConfig(config_filename)
        cfg1.cleanup() # remove templates and strip strings
        cfg1.save("test.yaml")
        cfg2 = SequanaConfig("test.yaml")

        assert cfg2._yaml_code == cfg1._yaml_code
        cfg2._update_config()
        assert cfg1.config == cfg2.config

def test_message():
    snaketools.message("test")

def test_dummy_manager():
    ss = snaketools.DummyManager()
    ss = snaketools.DummyManager(["test1.fastq.gz", "test2.fastq.gz"])
    assert ss.paired == True
    ss = snaketools.DummyManager(["test1.fastq.gz"])
    assert ss.paired == False
    ss = snaketools.DummyManager("test1.fastq.gz")
    assert ss.paired == False

def test_pipeline_manager():

    # test missing input_directory
    cfg = SequanaConfig({}, mode="other")
    try:
        pm = snaketools.PipelineManager("custom", cfg)
        assert False
    except:
        assert True

    # normal behaviour but no input provided:
    config = Module("quality_control")._get_config()
    cfg = SequanaConfig(config)
    cfg.cleanup() # remove templates
    try:
        pm = snaketools.PipelineManager("custome", cfg)
        assert False
    except:
        assert True

    cfg = SequanaConfig(config)
    cfg.cleanup() # remove templates
    file1 = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz")
    file2 = sequana_data("Hm2_GTGAAA_L005_R2_001.fastq.gz")
    cfg.config.input_samples['file1'] = file1
    cfg.config.input_samples['file2'] = file2
    pm = snaketools.PipelineManager("custome", cfg)
    assert pm.paired == True

    cfg = SequanaConfig(config)
    cfg.cleanup() # remove templates
    file1 = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz")
    cfg.config.input_samples['file1'] = file1
    pm = snaketools.PipelineManager("custome", cfg)
    assert pm.paired == False

    pm.getlogdir("fastqc")
    pm.getwkdir("fastqc")
    pm.getrawdata()
    pm.getreportdir("test")
    pm.getname("fastqc")



def test_file_name_factory():
    import glob

    def inner_test(ff):
        len(ff)
        print(ff)
        ff.filenames
        ff.realpaths
        ff.all_extensions
        ff.pathnames
        ff.extensions

    #list
    list_files = glob.glob("*.py")
    ff = snaketools.FileFactory(list_files)
    inner_test(ff)

    # glob
    ff = snaketools.FileFactory("*py")
    inner_test(ff)


    directory = os.path.dirname(sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz"))

    ff = snaketools.FastQFactory(directory + "/*fastq.gz", verbose=True)
    assert ff.tags == ['Hm2_GTGAAA_L005']

    ff.get_file1(ff.tags[0])
    ff.get_file2(ff.tags[0])
    assert len(ff) == 1
