.. _pipeline_user_guide:

Pipeline user guide
###################

This guide walks through running a Sequana pipeline end-to-end: install,
configure, run, inspect. The same skeleton applies to every pipeline.

.. contents::
   :local:
   :depth: 2


Install a pipeline
==================

Each pipeline is an independent PyPI package. Install it inside a Python
3.10+ environment alongside Sequana::

    pip install sequana_fastqc --upgrade

Verify::

    sequana_fastqc --help


Initialise the working directory
================================

Pipelines come with a config file, a snakefile, and a runner script. The
initialisation command copies them into a working directory of your choice
(default: the pipeline name)::

    sequana_fastqc --input-directory my_data --working-directory test1
    cd test1

Inside ``test1`` you will typically find:

- ``config.yaml`` — pipeline parameters (input/output, tools, resources).
- ``<name>.rules`` — the Snakefile.
- ``<name>.sh`` — convenience launcher (snakemake with sensible defaults).
- ``apptainers.yaml`` — container URIs used when ``--use-apptainer`` is set.


Common CLI options
==================

Every Sequana pipeline understands these options (with sensible defaults):

``--input-directory``
    Where to look for input FASTQ files. Default: ``.``.
``--input-pattern``
    Glob to select files. Default: ``*fastq.gz``. Use ``*/*fastq.gz`` if
    samples sit in sub-directories.
``--input-readtag``
    Pattern used to detect paired-end reads. Default: ``_R[12]_``.
``--working-directory``
    Where the pipeline files get copied. Use ``--force`` to overwrite.
``--run-mode {local,slurm}``
    Run locally or generate a SLURM-aware launcher. Auto-detected when
    ``sbatch`` is on the path.
``--use-apptainer``
    Pull and execute every rule inside the matching apptainer image.
``--deps``
    Print external dependencies and check whether they are installed.

Use ``sequana_<name> --help`` to discover pipeline-specific flags.


Edit the configuration
======================

``config.yaml`` is plain YAML. The most common fields are at the top::

    input_directory: /abs/path/to/data
    input_readtag: _R[12]_
    input_pattern: '*fastq.gz'

Every tool used by the pipeline has its own section (``cutadapt:``,
``bwa:``, ``coverage:`` …). The defaults are tuned for typical datasets; tweak
them for unusual cases (short genomes, very deep coverage, single-end data,
…).


Run the pipeline
================

Two equivalent ways::

    sh <pipeline>.sh

or directly through snakemake::

    snakemake -s <pipeline>.rules -j 4 -p

``-j N`` sets the number of parallel jobs, ``-p`` prints shell commands.
On a SLURM cluster, the generated ``<pipeline>.sh`` already includes the
cluster profile.

When the run is complete, the HTML report is at::

    ./summary.html


Clean up
========

To remove temporary files but keep the report::

    make clean


Tips
====

- Always run ``--deps`` once after installing a new pipeline.
- Re-run with ``--force`` to overwrite an existing working directory.
- For long runs on a cluster, prefer ``--use-apptainer`` — it pins the tool
  versions and eliminates conda-env clashes.
- See :ref:`pipelines` for the full pipeline catalogue and
  :ref:`tutorial` for end-to-end examples.
