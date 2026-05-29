.. _tutorial:

Tutorial
==========

Hands-on examples for a few representative Sequana pipelines and standalones.
Each section assumes you have followed the :ref:`installation`. For a complete
list of pipelines and the matching install commands, see :ref:`pipelines`.

.. contents::
   :local:
   :depth: 2


The ``sequana`` CLI
-------------------

The top-level ``sequana`` command groups ~30 helper sub-commands (FASTQ/FASTA
utilities, GFF/GTF fixers, enrichment tools, summary, telomark, etc.). List
them with::

    sequana --help

For shell completion (bash, zsh, fish), follow the
`Click instructions <https://click.palletsprojects.com/en/stable/shell-completion/>`_,
e.g. for bash::

    _SEQUANA_COMPLETE=bash_source sequana > ~/.sequana-complete.bash
    echo '. ~/.sequana-complete.bash' >> ~/.bashrc


fastqc pipeline
---------------

We will run the
`sequana_fastqc <https://github.com/sequana/fastqc>`_ pipeline on a pair of
FastQ files. The data is a Measles virus sequencing run (HiSeq, PCRFree
adapters, ~10% adapter content, index ``GTGAAA``). For testing, download:

- :download:`R1 <../sequana/resources/data/Hm2_GTGAAA_L005_R1_001.fastq.gz>`
- :download:`R2 <../sequana/resources/data/Hm2_GTGAAA_L005_R2_001.fastq.gz>`)

(1500 reads each.) Then::

    pip install sequana_fastqc --upgrade
    sequana_fastqc --input-directory .
    cd fastqc
    sh fastqc.sh

Open ``summary.html`` in your browser.


Taxonomy (standalone)
---------------------

Quick taxonomy classification of a FastQ file uses ``sequana_taxonomy``
(see :ref:`standalone_sequana_taxonomy`).

Download a toy Kraken database (100 FASTA files, measles + a few viruses)::

    from sequana import KrakenDownload, sequana_config_path
    kd = KrakenDownload()
    kd.download("toydb")
    database_path = sequana_config_path + "/kraken_toydb"

Then either drive it from Python::

    from sequana import KrakenPipeline
    kp = KrakenPipeline(["R1.fastq.gz", "R2.fastq.gz"],
                        database="~/.config/sequana/kraken_toydb")
    kp.run()

…or from the shell::

    sequana_taxonomy --file1 Test_R1.cutadapt.fastq.gz \
                     --file2 Test_R2.cutadapt.fastq.gz \
                     --database <database_path>

Open ``taxonomy/kraken.html`` (Krona pie chart). A reference rendering is
available `here <_static/krona.html>`_.


Variant calling pipeline
------------------------

The variant calling pipeline performs mapping, calling and annotation
(snpEff + coverage). See :ref:`pipeline_vc` for full details.

Install::

    pip install sequana_variant_calling --upgrade

Fetch a reference + GenBank annotation. Using the Measles accession
``K01711.1``::

    from sequana.snpeff import download_fasta_and_genbank
    download_fasta_and_genbank("K01711", "measles")

Then run::

    sequana_variant_calling --input-directory . \
                            --reference measles.fa \
                            --annotation measles.gbk
    cd variant_calling
    sh variant_calling.sh

When the pipeline finishes, ``index.html`` shows the multiqc summary, and each
sample has its own ``report_SAMPLENAME/summary.html``.

About the configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The reference and annotation paths are set when initiating the pipeline. Open
``config.yaml`` to tune the rest::

    annotation_file: measles.gbk
    reference_file: measles.fa

.. warning:: ``mark_duplicates`` outputs can be huge and require scratch space
   on a cluster.

.. warning:: In the ``coverage`` section, reduce the window size for short
   genomes.


De-novo assembly
----------------

::

    pip install sequana_denovo --upgrade
    sequana_denovo --input-directory . --working-directory denovo_test

Edit ``denovo_test/config.yaml`` before running. The
``digital_normalisation`` section is the main memory knob — for tests set
``max-tablesize`` to ``3e6``.


RNA-seq
-------

Full reference: :ref:`pipeline_rnaseq`. Quick recipe (single-end yeast data,
HiSeq2500)::

    pip install sequana_rnaseq --upgrade

Get the genome + GFF for `Saccer3
<http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/>`_::

    mkdir Saccer3 && cd Saccer3
    wget http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz
    tar -xvzf chromFa.tar.gz
    cat *.fa > Saccer3.fa
    wget http://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff -O Saccer3.gff
    rm -f chr*
    cd ..

.. warning:: All reference files must share the same prefix
   (``Saccer3`` here) and live in the same directory.

.. warning:: The counting step expects GFF only (GTF/SAF coming).

Run::

    sequana_rnaseq --genome-directory Saccer3 --aligner bowtie2
    cd rnaseq
    sh rnaseq.sh

On a SLURM cluster the same pipeline can be run with::

    sbatch sh rnaseq.sh --profile slurm

(see the pipeline README for the cluster profile setup.)


Apptainer containers
--------------------

Every Sequana pipeline ships an ``apptainers.yaml`` that points at containers
maintained by the `damona <https://damona.readthedocs.io>`_ project. Pipelines
pull them automatically when invoked with ``--use-apptainer``::

    sequana_fastqc --input-directory . --use-apptainer

This is the recommended way to avoid conda/system tool clashes.
