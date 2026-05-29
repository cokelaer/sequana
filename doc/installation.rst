.. _installation:

Installation
##########################################

.. contents::
   :local:
   :depth: 2


Quick installation
==================

::

    pip install sequana

works out of the box on Python 3.10, 3.11 or 3.12.

For most pipelines you will also need a few non-Python tools (e.g. samtools,
bwa, fastqc, kraken2). The easiest way to get those is to combine
``pip`` with a `mamba <https://mamba.readthedocs.io>`_ / conda environment, or
to use the provided apptainer containers.


.. _installation_pip:

Recommended: virtual environment + pip
======================================

Sequana is published on `PyPI <https://pypi.org/project/sequana>`_. We strongly
recommend installing it in a dedicated environment so that dependencies do not
collide with system packages.

::

    mamba create -n sequana_env "python=3.12"
    mamba activate sequana_env
    pip install --upgrade sequana


To install a Sequana **pipeline**, install the matching PyPI package in the
same environment::

    pip install sequana_fastqc
    pip install sequana_rnaseq
    pip install sequana_variant_calling
    # ...

Each pipeline lives in its own repository under
https://github.com/sequana and ships its own non-Python dependencies.
You can let the pipeline tell you what is missing::

    sequana_fastqc --deps

A curated bundle of pipelines can be installed in one go::

    pip install "sequana[pipelines]"


.. _installation_conda:

Bioconda
========

Sequana is also published on `bioconda <https://bioconda.github.io/recipes/sequana/README.html>`_
but the bioconda recipe is usually behind the PyPI release. Use it only if you
specifically want a conda package.

Set up the channels once::

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

Then::

    mamba create -n sequana_env sequana
    mamba activate sequana_env


.. _installation_apptainer:

Apptainer / Singularity containers
==================================

All containers used by Sequana pipelines are produced by the
`damona <https://damona.readthedocs.io>`_ project. You do not need to download
them manually — Sequana pipelines pull them on demand when you pass
``--use-apptainer``.

If you want to grab a specific Sequana image yourself, browse
https://damona.readthedocs.io for the available releases.


From GitHub (development version)
=================================

If you want the unreleased code or wish to contribute, install Sequana with
`poetry <https://python-poetry.org>`_::

    git clone https://github.com/sequana/sequana
    cd sequana
    poetry install --with dev

See the :ref:`developers` section for the test, lint and docs workflow.


.. _installation_sequanix:

Sequanix (GUI)
==============

Since v0.16, **Sequanix** (the graphical front-end for Snakemake workflows)
lives in its own repository: https://github.com/sequana/sequanix. Install it
separately if you need a GUI.


History
=======

Before v0.8, all pipelines shipped inside this repository. The number of
pipelines and the variety of non-Python dependencies (BWA, Kraken, snpEff,
DESeq2…) made a single install impractical, so each pipeline now lives in its
own repository with its own release cycle. See https://github.com/sequana for
the up-to-date list.
