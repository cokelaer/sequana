def yield_bigwig_by_chromosome(filename):
    """Yield ``(chrom, entries)`` from a binary BigWig file.

    ::

        for chrom, entries in yield_bigwig_by_chromosome("coverage.bw"):
            print(f"Chromosome: {chrom}, {len(entries)} entries")
            print(entries[:3])

    Each ``entries`` is a list of ``(start, value)`` tuples, one per interval
    stored in the BigWig (e.g. one per ``bamCoverage`` bin). ``pyBigWig`` is
    imported lazily.
    """
    try:
        import pyBigWig
    except ImportError:  # pragma: no cover
        raise ImportError(
            "Reading BigWig files requires the 'pyBigWig' package. " "Install it with 'pip install pyBigWig'."
        )

    bw = pyBigWig.open(str(filename))
    try:
        for chrom in bw.chroms():
            intervals = bw.intervals(chrom)
            if not intervals:  # pragma: no cover
                continue
            current_data = [(start, value) for start, end, value in intervals]
            yield chrom, current_data
    finally:
        bw.close()
