.. _userguide:

Library user guide
##################

Use this page if you intend to call **Sequana** from Python or a Jupyter
notebook — to read NGS file formats, compute metrics, or assemble report
sections from a script.

For end-users running pre-built pipelines, see :ref:`pipeline_user_guide`
instead.

.. contents::
   :local:
   :depth: 2


Test data
=========

Sequana ships a small data folder accessible from Python via
:func:`~sequana.sequana_data`::

    from sequana import sequana_data
    filename = sequana_data('JB409847.bed')

A complete list of bundled files is in :class:`sequana.datatools`.


Coverage from a BED file
========================

Read a BED file produced by ``bedtools genomecov``::

    from sequana import SequanaCoverage
    gc = SequanaCoverage(filename)

Select a chromosome, compute the running median and z-score::

    chrom = gc[0]
    chrom.running_median(n=5001, circular=True)
    chrom.compute_zscore()

Plot the coverage with its 3-sigma confidence band::

    chrom.plot_coverage()

.. plot::

    from sequana import sequana_data
    filename = sequana_data('JB409847.bed')
    from sequana import SequanaCoverage
    gc = SequanaCoverage(filename)
    chrom = gc[0]
    chrom.running_median(n=5001, circular=True)
    chrom.compute_zscore()
    chrom.plot_coverage()

A matching notebook is at `notebooks/coverage.ipynb
<https://github.com/sequana/sequana/blob/main/notebooks/coverage.ipynb>`_.


FastQ inspection
================

The :class:`~sequana.fastqc.FastQC` class exposes per-read metrics::

    from sequana import FastQC, sequana_data

    fastqc = FastQC(sequana_data("test.fastq"))
    print(fastqc.fastq)
    for x in 'ACGT':
        fastqc.get_actg_content()[x].hist(
            alpha=0.5, label=x, histtype='step', lw=3, bins=10)

.. plot::
    :include-source:

    from sequana import FastQC, sequana_data
    fastqc = FastQC(sequana_data("test.fastq"))
    for x in 'ACGT':
        fastqc.get_actg_content()[x].hist(
            alpha=0.5, label=x, histtype='step', lw=3, bins=10)
    from pylab import legend
    legend()


Reading sequence data (FASTA/FASTQ)
===================================

The :class:`~sequana.fasta.FastA` and :class:`~sequana.fastq.FastQ` classes
read and manipulate sequence files::

    from sequana import FastA, sequana_data
    fasta = FastA(sequana_data("measles.fasta"))
    for record in fasta:
        print(record.name, len(record.seq))

Access GC content and other sequence metrics via :class:`~sequana.sequence.DNA`.

Annotations (GFF3/GenBank)
===========================

Parse genome annotations with :class:`~sequana.gff3.GFF3` or
:class:`~sequana.genbank.GenBank`::

    from sequana import GFF3, sequana_data
    gff = GFF3(sequana_data("annotations.gff3"))
    genes = [r for r in gff if r.feature == "gene"]

Variants (VCF)
==============

The :class:`~sequana.vcftools.VCF` class reads and filters VCF files::

    from sequana import VCF, sequana_data
    vcf = VCF(sequana_data("variants.vcf"))
    for variant in vcf:
        print(variant.CHROM, variant.POS, variant.REF, variant.ALT)

Taxonomy (Kraken)
=================

Classify sequences and parse Kraken output with :mod:`sequana.kraken`::

    from sequana.kraken import KrakenResults
    kr = KrakenResults(sequana_data("kraken.out"))
    kr.plot()  # Krona pie chart

Building HTML report sections from Python
=========================================

Sequana's pipeline reports are assembled from reusable building blocks in
:mod:`sequana.modules_report`. You can call them on your own data.
Example for a BAM file::

    from sequana import BAM, sequana_data
    from sequana.modules_report.bamqc import BAMQCModule
    BAMQCModule(sequana_data("test.bam"), "bam.html")

The generated ``bam.html`` is a self-contained page (see
`bam.html <_static/bam.html>`_ for a rendered example).

To build a brand-new report module, see :ref:`module_reports` in the
developer guide.

Feature counting and differential expression
==============================================

Count reads per feature (gene, exon, etc.) with :class:`~sequana.featurecounts.FeatureCounts`::

    from sequana import FeatureCounts
    fc = FeatureCounts("counts.txt")
    df = fc.df  # pandas DataFrame of gene counts

For differential expression, see :mod:`sequana.rnadiff` (wraps DESeq2) and
:mod:`sequana.enrichment` for KEGG/GSEA enrichment analysis.


Where to look next
==================

- :ref:`API reference <references>` — the full module index.
- :doc:`auto_examples/index` — Sphinx-Gallery of short, runnable scripts.
- :doc:`notebooks` — Jupyter notebooks demonstrating BAM, coverage, FastQ,
  feature counts and ribodesigner workflows.
- :ref:`cli_reference` — the ``sequana`` CLI sub-commands (some of these
  wrap the Python API).
