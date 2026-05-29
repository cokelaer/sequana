.. _pipelines:

Pipelines
##########

Sequana ships a family of NGS pipelines, each in its own GitHub repository and
PyPI package. They share a common installation pattern, CLI shape, and HTML
report layout.

.. contents::
   :local:
   :depth: 2


Quick start
===========

::

    pip install sequana_<name> --upgrade
    sequana_<name> --help

For example::

    pip install sequana_fastqc --upgrade
    sequana_fastqc --input-directory my_data
    cd fastqc && sh fastqc.sh

The configuration file ``config.yaml`` lives next to the snakefile in the
working directory and can be edited before launch. See
:ref:`pipeline_user_guide` for a walkthrough of the common options.

Containers
==========

Every Sequana pipeline ships an ``apptainers.yaml`` pointing at the
matching `damona <https://damona.readthedocs.io>`_ image. Pull them on demand
with::

    sequana_<name> ... --use-apptainer

This is the recommended way to avoid conda / system tool clashes.


Pipeline catalogue
==================

Quality control
---------------

.. list-table::
   :widths: 18 50 32
   :header-rows: 1

   * - Pipeline
     - Description
     - Repository
   * - ``sequana_fastqc``
     - Sequencing quality control (FastQC + MultiQC).
     - https://github.com/sequana/fastqc
   * - ``sequana_pacbio_qc``
     - PacBio long-read quality control.
     - https://github.com/sequana/pacbio_qc
   * - ``sequana_ribofinder``
     - Estimate ribosomal content of a sample.
     - https://github.com/sequana/ribofinder

Mapping and coverage
--------------------

.. list-table::
   :widths: 18 50 32
   :header-rows: 1

   * - Pipeline
     - Description
     - Repository
   * - ``sequana_mapper``
     - Map reads onto a target genome (bowtie2 / bwa / minimap2).
     - https://github.com/sequana/mapper
   * - ``sequana_multicov``
     - Multi-sample coverage analysis.
     - https://github.com/sequana/multicov

Variants and RNA-seq
--------------------

.. list-table::
   :widths: 18 50 32
   :header-rows: 1

   * - Pipeline
     - Description
     - Repository
   * - ``sequana_variant_calling``
     - Variant calling (freebayes + snpEff + coverage).
     - https://github.com/sequana/variant_calling
   * - ``sequana_rnaseq``
     - Full RNA-seq pipeline (alignment, counts, DESeq2 report).
     - https://github.com/sequana/rnaseq

Long-read pipelines
-------------------

.. list-table::
   :widths: 18 50 32
   :header-rows: 1

   * - Pipeline
     - Description
     - Repository
   * - ``sequana_lora``
     - Map long reads onto a target genome.
     - https://github.com/sequana/lora
   * - ``sequana_nanomerge``
     - Merge barcoded (or unbarcoded) Nanopore FASTQ files + report.
     - https://github.com/sequana/nanomerge
   * - ``sequana_laa``
     - Long-read amplicon analysis.
     - https://github.com/sequana/laa

Taxonomy and assembly
---------------------

.. list-table::
   :widths: 18 50 32
   :header-rows: 1

   * - Pipeline
     - Description
     - Repository
   * - ``sequana_multitax``
     - Taxonomic profiling of multiple samples.
     - https://github.com/sequana/multitax
   * - ``sequana_denovo``
     - De-novo assembly (digital normalisation + assembler + QC).
     - https://github.com/sequana/denovo

Utilities
---------

.. list-table::
   :widths: 18 50 32
   :header-rows: 1

   * - Pipeline
     - Description
     - Repository
   * - ``sequana_demultiplex``
     - Demultiplex raw Illumina runs (bcl2fastq / bcl-convert).
     - https://github.com/sequana/demultiplex
   * - ``sequana_downsampling``
     - Downsample FASTQ/BAM files.
     - https://github.com/sequana/downsampling
   * - ``sequana_depletion``
     - Remove or select reads mapping a reference.
     - https://github.com/sequana/depletion
   * - ``sequana_revcomp``
     - Reverse-complement sequence data.
     - https://github.com/sequana/revcomp
   * - ``sequana_trf``
     - Tandem repeats finder wrapper.
     - https://github.com/sequana/trf

For the complete and always-up-to-date list see the
`sequana GitHub organisation <https://github.com/sequana>`_.

Pipeline reference pages
========================

A handful of pipelines have a dedicated page in this manual. They embed the
README of the corresponding repository via the ``sequana_pipeline`` Sphinx
directive.

.. toctree::
   :maxdepth: 1

   pipeline_demultiplex
   pipeline_fastqc
   pipeline_mapper
   pipeline_pacbio_qc
   pipeline_ribofinder
   pipeline_rnaseq
   pipeline_vc


Naming convention
=================

Each pipeline is distributed as ``sequana_<name>`` (PyPI) and installs a
console script named ``sequana_<name>``. For historical pipelines also exposed
as standalone tools (``sequana_coverage``), the pipeline form is
``sequana_pipelines_<name>`` to disambiguate.
