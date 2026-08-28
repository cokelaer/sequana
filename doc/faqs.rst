.. _faqs:

FAQ / Troubleshooting
=====================

.. contents::
   :local:
   :depth: 2


Installation
------------

What are the dependencies?
~~~~~~~~~~~~~~~~~~~~~~~~~~

Two flavours:

- **Python** libraries (numpy, pandas, matplotlib, …). Installed automatically
  by ``pip install sequana``.
- **External** tools used by pipelines (bwa, samtools, kraken2, fastqc, …).
  Install them via bioconda, your system package manager, or simply run the
  pipeline with ``--apptainer-prefix`` to skip the question entirely.

Sequana itself only needs ``kraken2``, ``cd-hit`` and ``krona`` to be on the
``$PATH`` (those are used by the ``sequana_taxonomy`` standalone).

matplotlib
~~~~~~~~~~

If you see X11 errors when matplotlib tries to open a window (e.g. on a
cluster), force the headless backend::

    mkdir -p ~/.config/matplotlib
    echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc

Then start a new shell.


Input data
----------

Expected file naming convention
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most pipelines expect gzipped FastQ files following the Illumina naming
convention::

    SAMPLE_S1_L001_R1_001.fastq.gz
    SAMPLE_S1_L001_R2_001.fastq.gz

The ``_R1_`` / ``_R2_`` tag identifies paired files. The ``input_readtag``
parameter in the pipeline config accepts custom patterns such as ``_R[12]_``
or ``_R1\.`` if your files don't follow the Illumina naming scheme.


Pipeline runs
-------------

What to do if a pipeline run fails
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Common causes, in decreasing order of frequency:

1. **Bad input pattern** — empty sample set, wrong ``--input-readtag``.
2. **Missing config value** — open ``config.yaml`` and check required fields.
3. **Cluster resources** — job killed because not enough memory was
   allocated. Bump the ``resources`` section in ``config.yaml`` or in the
   SLURM profile.
4. **Pipeline bug** — report on
   ``https://github.com/sequana/<pipeline>/issues``.

For verbose output, use snakemake's built-in flags::

    snakemake -s <pipeline>.rules --printshellcmds

The wrapper script ``<pipeline>.sh`` does not forward arbitrary flags to snakemake.

Variant Calling — snpEff "Cannot find sequence" error
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If snpEff fails with::

    java.lang.RuntimeException: Cannot find sequence for 'LN831026.gbk'

…your GenBank file is missing the embedded sequence (header only). Re-download
the file using::

    from sequana.snpeff import download_fasta_and_genbank
    download_fasta_and_genbank("LN831026", "myref")


Library (Python API)
--------------------

How do I import Sequana without slowing down my script?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sequana uses lazy imports for heavy dependencies. You can speed up imports by
being selective::

    from sequana import FastA, VCF  # Fast — only loads what you need

    from sequana import *            # Slower — loads everything

How do I plot in a headless environment (cluster, CI/CD)?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Force matplotlib to use a non-interactive backend::

    import matplotlib
    matplotlib.use('Agg')  # Before importing Sequana or matplotlib.pyplot

    from sequana import SequanaCoverage
    cov = SequanaCoverage("file.bed")
    cov[0].plot_coverage()
    # Image saves, window doesn't open

How do I handle large files (BAM, VCF)?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sequana file readers are iterator-based to keep memory usage low::

    from sequana import VCF
    vcf = VCF("large_file.vcf.gz")
    for variant in vcf:  # Processes one variant at a time
        print(variant.CHROM, variant.POS)


PacBio
------

pbindex: "read group ID not found"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If ``pbindex`` complains::

    FATAL pbindex ERROR: [pbbam] BAM header ERROR: read group ID not found: ...

…the BAM is missing the ``@RG`` line (typical after sub-sampling). Re-attach
the original header::

    samtools view -H ORIGINAL.bam | grep '@RG' > new_header.txt
    samtools reheader new_header.txt sample.bam > corrected_sample.bam


Apptainer / Singularity
-----------------------

Base home directory does not exist
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If apptainer aborts with::

    ERROR  : Base home directory does not exist within the container: /pasteur

…the container has no entry for that path. With sudo access::

    sudo apptainer shell --writable sequana.sif
    mkdir /pasteur
    exit

Without sudo, prepare the image on a machine where you have sudo and ship the
result back. The cleanest workaround is to bind the host path explicitly::

    apptainer exec -B /pasteur:/pasteur sequana.sif sequana_fastqc ...


Open an issue
-------------

If nothing above helps, please open an issue against the affected repository
(pipeline-specific issues go to ``https://github.com/sequana/<pipeline>``;
library issues to https://github.com/sequana/sequana/issues). Include:

- Sequana version (``sequana --version``)
- Pipeline name and version
- The failing command and the full traceback or snakemake log
