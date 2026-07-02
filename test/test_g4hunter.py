import glob

from sequana.repeats.G4hunter import G4Hunter, G4HunterReader

from . import test_dir

fasta_file = f"{test_dir}/data/fasta/measles.fa"


def test_base_score():
    g = G4Hunter(fasta_file)
    scores = list(g.base_score("GGGGAAAACCCC"))
    # G-runs give positive (capped at 4), C-runs negative (capped at -4), other 0
    assert scores[:4] == [4, 4, 4, 4]
    assert scores[4:8] == [0, 0, 0, 0]
    assert scores[8:] == [-4, -4, -4, -4]


def test_g4hunter_run(tmp_path):
    g = G4Hunter(fasta_file, window=25, score=1)
    g.run(str(tmp_path))
    produced = sorted(p.split("/")[-1] for p in glob.glob(f"{tmp_path}/*"))
    assert any(p.endswith("-Merged.txt") for p in produced)
    assert any(p.endswith("-W25-S1.txt") for p in produced)


def test_g4hunter_reader_and_bed(tmp_path):
    g = G4Hunter(fasta_file, window=25, score=1)
    g.run(str(tmp_path))
    merged = glob.glob(f"{tmp_path}/*-Merged.txt")[0]

    reader = G4HunterReader(merged)
    assert len(reader.data_merged) > 0
    assert sum(len(v) for v in reader.data_merged.values()) > 0

    bed = tmp_path / "g4.bed"
    reader.to_bed(str(bed))
    assert bed.exists() and bed.stat().st_size > 0


def test_g4hunter_reader_load_files(tmp_path):
    g = G4Hunter(fasta_file, window=25, score=1)
    g.run(str(tmp_path))
    reader = G4HunterReader()
    reader.load_files(f"{tmp_path}/*-Merged.txt")
    assert len(reader.data_merged) > 0
