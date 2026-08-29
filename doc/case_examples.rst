.. _case_examples:

Case examples
=============

Real-world end-to-end walkthroughs combining Sequana pipelines and analysis.


Quality control + variant calling (Measles virus)
==================================================

**Goal:** Assess raw sequencing quality, then map reads and call variants on a
reference genome.

**Data:** Paired-end MiSeq run (Measles virus, ~10% adapter content)

**Step 1: Quality Control**

Run fastqc to inspect raw read quality::

    pip install sequana_fastqc --upgrade
    sequana_fastqc --input-directory raw_data --execute

Open ``summary.html`` to review read length distribution, adapter content, and
base quality. Sequana will show coverage patterns and adapter percentages.

**Step 2: Variant Calling**

If quality is acceptable, run variant calling on the same data::

    pip install sequana_variant_calling --upgrade

Download reference genome and annotation::

    from sequana.snpeff import download_fasta_and_genbank
    download_fasta_and_genbank("K01711.1", "measles")

Initialize pipeline::

    sequana_variant_calling --input-directory raw_data \
                            --reference measles.fa \
                            --annotation measles.gbk \
                            --execute

Output: ``summary.html`` contains MultiQC report with variant calls, coverage,
snpEff annotations, and per-sample variant summaries.

**Step 3: Combine and Interpret**

Using Python, load and compare QC metrics with variant results::

    from sequana import FastQC, VCF

    # Load quality metrics from fastqc output
    qc = FastQC("raw_data/fastqc_results.json")
    print(f"Mean quality: {qc.get_quality()}")

    # Load variant calls
    vcf = VCF("variant_calling/variants.vcf")
    variants = list(vcf)
    print(f"Total variants: {len(variants)}")

    # Filter high-confidence variants (snpEff score > 50)
    high_conf = [v for v in variants if float(v.INFO.get("ANN", [0])[0]) > 50]
    print(f"High-confidence: {len(high_conf)}")

See :ref:`tutorial` for more pipeline recipes (RNA-seq, de-novo assembly, taxonomy).
