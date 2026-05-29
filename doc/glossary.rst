Glossary
========

.. glossary::
    :sorted:

    Apptainer
        Container runtime used by Sequana pipelines (formerly known as
        Singularity). Images are pulled on demand via ``--use-apptainer``.

    BAI
        Index file accompanying a :term:`BAM` file (non-standard extension
        ``.bai``).

    BAM
        Binary version of the :term:`SAM` alignment format.

    BED
        Tab-separated format describing genomic intervals (chrom, start, end,
        plus optional fields). Used for coverage reporting.

    CIGAR
        Compact string describing how a read aligns to a reference
        (matches, insertions, deletions, soft-clips, …).

    Conda environment
        Isolated Python + binary environment managed by
        ``conda``/``mamba``. Recommended for installing Sequana to avoid
        clashing with system packages.

    DEG
        Differentially Expressed Gene. Output of RNA-seq differential
        expression analyses (see :class:`sequana.rnadiff`).

    DSRC
        A compression tool dedicated to FastQ files.

    FASTA
        Plain-text format for nucleotide or protein sequences. Each record
        starts with a ``>`` header line, followed by one or more sequence
        lines.

    FASTQ
        Plain-text format combining a sequence with per-base quality scores.
        Typically gzipped (``.fastq.gz``).

    GFF
        General Feature Format. Tab-separated annotation file describing
        genes, exons and other features along a sequence. See also GFF3.

    GTF
        Gene Transfer Format. Tab-separated annotation similar to :term:`GFF`
        but with stricter attribute conventions.

    JSON
        Human-readable data serialisation language commonly used for
        configuration. See https://en.wikipedia.org/wiki/JSON.

    k-mer
        Substring of fixed length *k* from a sequence. Used in classification
        (Kraken), assembly (de Bruijn graphs), and quality control.

    MultiQC
        Aggregates outputs of many bioinformatics tools (FastQC, samtools,
        cutadapt …) into a single interactive HTML report. Sequana ships
        plugins for several of its pipelines.

    Module
        A directory that contains a Snakemake rule and an associated
        README. Especially relevant for the Sequana pipelines. See
        :ref:`developers`.

    NGS
        Next-generation sequencing. Catch-all term for high-throughput
        sequencing technologies (Illumina, PacBio, ONT …).

    Rule
        Smallest unit of a Snakemake workflow: declares input, output, and
        the shell/Python code that turns one into the other.

    SAM
        Sequence Alignment Map. Tab-separated format describing read
        alignments to a reference.

    Sample sheet
        Tabular description of an Illumina run (sample IDs, indices,
        adapters). Parsed by :class:`sequana.iem.SampleSheet`.

    Snakefile
        A file that contains one or several Snakemake rules.

    Snakemake
        Python-based workflow engine used by every Sequana pipeline.

    Taxon
        A taxonomic unit (species, genus, …) referenced by an NCBI taxid.
        Sequana's :mod:`sequana.taxonomy` loads the NCBI taxonomy dump.

    VCF
        Variant Call Format, used by the variant calling pipeline.

    Wrapper
        Reusable Snakemake rule shipped by
        `sequana-wrappers <https://github.com/sequana/sequana-wrappers>`_ and
        consumed by every pipeline via Snakemake's native ``wrapper:``
        directive.

    YAML
        Human-readable data serialisation language commonly used for
        configuration files. Every Sequana pipeline configuration is YAML.
