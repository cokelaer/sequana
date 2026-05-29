.. _cli_reference:

CLI reference
=============

The top-level ``sequana`` command groups ~30 sub-commands (FASTQ/FASTA
utilities, GFF/GTF fixers, enrichment helpers, summaries, …). The reference
below is generated directly from the source.

Each sub-command also accepts ``--help`` from the shell, e.g.::

    sequana fastq --help
    sequana enrichment-kegg --help

.. click:: sequana.scripts.main.main:main
   :prog: sequana
   :nested: full
