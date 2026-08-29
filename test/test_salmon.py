from sequana.salmon import Salmon

from . import test_dir

data_dir = f"{test_dir}/data"

# an Ensembl-like annotation: identifiers are prefixed (ID=gene:GENE1), genes are
# spread over three features (gene, ncRNA_gene, pseudogene) and the sub-features
# (exon, CDS) carry the transcript_id of their transcript.
GFF = f"{data_dir}/salmon/annotation.gff"
QUANT1 = f"{data_dir}/salmon/sample1_quant.sf"
QUANT2 = f"{data_dir}/salmon/sample2_quant.sf"


def _get_data(output):
    """Return the feature counts as a dictionary of rows indexed by Geneid"""
    header, *rows = output.split("\n")
    columns = header.split("\t")
    return {row.split("\t")[0]: dict(zip(columns, row.split("\t"))) for row in rows}


def test_salmon_eukaryotes():
    salmon = Salmon(QUANT1, GFF)
    data = _get_data(salmon.get_feature_counts())

    # genes are kept whatever their feature (gene, ncRNA_gene, pseudogene)
    assert sorted(data.keys()) == ["GENE1", "GENE2", "GENE3", "GENE4"]

    # the gene:  prefix is removed but the coordinates are those of the gene
    assert data["GENE1"]["Chr"] == "1"
    assert data["GENE1"]["Start"] == "100"
    assert data["GENE1"]["End"] == "2000"
    assert data["GENE1"]["Strand"] == "+"

    # counts of the two transcripts of GENE1 are summed (10.4 + 5.2 = 15.6) and
    # rounded: DESeq2 does not accept fractional counts.
    assert data["GENE1"]["sample1_quant.sf"] == "16"
    assert data["GENE2"]["sample1_quant.sf"] == "7"
    assert data["GENE4"]["sample1_quant.sf"] == "0"
    for row in data.values():
        assert float(row["sample1_quant.sf"]) == int(row["sample1_quant.sf"])

    # TPM-weighted effective length: (100*1000 + 300*500) / (100+300) = 625
    assert float(data["GENE1"]["Length"]) == 625
    # a gene with no expression falls back on the mean effective length
    assert float(data["GENE4"]["Length"]) == 800


def test_salmon_no_subfeature_contamination():
    # exons and CDS of the annotation carry a transcript_id (as NCBI and GTF do).
    # They must not hide the transcript-to-gene relation.
    salmon = Salmon(QUANT1, GFF)
    salmon.get_feature_counts()

    assert salmon.trs2genes["TR1"] == "gene:GENE1"
    assert salmon.trs2genes["transcript:TR1"] == "gene:GENE1"
    assert "exon:TR1-1" not in salmon.trs2genes
    assert "CDS:PROT1" not in salmon.trs2genes
    assert all(x.startswith("gene:") for x in salmon.trs2genes.values())


def test_salmon_feature_filtering():
    # by default no filtering on the feature is performed
    assert len(_get_data(Salmon(QUANT1, GFF).get_feature_counts(feature=None))) == 4

    # restricting to 'gene' drops the ncRNA_gene and pseudogene entries
    assert sorted(_get_data(Salmon(QUANT1, GFF).get_feature_counts(feature="gene"))) == ["GENE1", "GENE4"]

    # several features may be provided
    data = _get_data(Salmon(QUANT1, GFF).get_feature_counts(feature="gene,ncRNA_gene"))
    assert sorted(data) == ["GENE1", "GENE2", "GENE4"]


def test_salmon_save_feature_counts(tmpdir):
    outfile = tmpdir.join("counts.out")
    Salmon(QUANT1, GFF).save_feature_counts(str(outfile))

    lines = [l for l in open(str(outfile)).read().splitlines() if l.strip()]
    assert lines[0].startswith("# Program:sequana.salmon")
    assert lines[1].startswith("Geneid\tChr\tStart\tEnd\tStrand\tLength")
    assert len(lines) == 6  # comment + header + 4 genes


def test_salmon_effective_length_is_sample_dependent():
    # the effective length of a gene depends on the isoforms expressed in the
    # sample. Merging count files on that column would drop the genes.
    data1 = _get_data(Salmon(QUANT1, GFF).get_feature_counts())
    data2 = _get_data(Salmon(QUANT2, GFF).get_feature_counts())
    assert data1["GENE1"]["Length"] != data2["GENE1"]["Length"]


def test_salmon_prokaryotes():
    # here the salmon identifiers are the genes themselves (no transcript: tag)
    salmon = Salmon(f"{data_dir}/salmon/prokaryote_quant.sf", f"{data_dir}/gff/ecoli_MG1655.gff")
    data = _get_data(salmon.get_feature_counts())

    assert sorted(data) == ["gene-b0001", "gene-b0002", "gene-b0003"]
    assert data["gene-b0001"]["Chr"] == "NC_000913.3"
    assert data["gene-b0001"]["Start"] == "190"
    assert data["gene-b0001"]["prokaryote_quant.sf"] == "120"
    for row in data.values():
        assert float(row["prokaryote_quant.sf"]) == int(row["prokaryote_quant.sf"])
