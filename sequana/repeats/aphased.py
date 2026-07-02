#!/usr/bin/env python3
"""Detect A-phased repeats (intrinsically bent DNA).

A-phased repeats are clusters of short A-tracts (runs of A/T) repeated in phase
with the ~10.5 bp helical period, which bends the double helix. This reproduces
the *A_Phased_Repeat* output of the non-B gfa / nBMST tool (``gfa_APR.gff``),
ported from its ``findAPR.c`` / ``getAtracts`` routines.

Algorithm (gfa defaults shown):

1. **A-tracts** — every maximal run of A/T whose length is in
   ``[min_tract_len, max_tract_len]`` = ``[3, 9]``. A run is a valid A-tract
   when, on the forward **or** reverse-complement strand, the longest A / AnTn
   pattern minus the longest pure-T run is ``>= min_tract_len``. Each tract is
   summarised by the centre of its longest A-run.
2. **Phasing** — consecutive tract centres separated by ``[min_sep, max_sep]``
   (gfa uses 9.9..11.1, i.e. ~10 bp) are chained; a chain of at least
   ``min_tracts`` (3) tracts is reported as one A-phased repeat.

Validated against gfa on a 30+ Mb *Leishmania* assembly: identical count (889)
and 100% reciprocal positional overlap with ``gfa_APR.gff``.
"""
import re

import pandas as pd
from tqdm import tqdm

from sequana import FastA
from sequana.tools import reverse_complement

_ATRUN = re.compile(r"[AT]+")


def _metric(ctx):
    """Longest A/AnTn pattern and longest pure-T run within ``ctx[1:]``.

    ``ctx[0]`` is the preceding-base context. Mirrors gfa's ``getAtracts`` inner
    loop. Returns (max_AT_len, max_T_len, end_index_of_longest_AT_within_run).
    """
    a_len = t_len = at_len = ta_len = 0
    max_at = max_t = 0
    max_at_end = -1
    for j in range(1, len(ctx)):
        ch = ctx[j]
        if ch == "A":
            t_len = ta_len = 0
            if ctx[j - 1] == "T":
                a_len = at_len = 0
            else:
                a_len += 1
                at_len += 1
        elif ch == "T":
            if ta_len < a_len:  # T following an A-run extends the AnTn pattern
                ta_len += 1
                at_len += 1
            else:  # pure T
                t_len += 1
                ta_len = at_len = a_len = 0
        if max_at < at_len:
            max_at = at_len
            max_at_end = j - 1
        if max_t < t_len:
            max_t = t_len
    return max_at, max_t, max_at_end


def _tract_center(seq, start, length, min_tract_len):
    """Return the A-run centre (float, forward coords) of a valid A-tract, else None."""
    run = seq[start : start + length]
    prev = seq[start - 1] if start > 0 else "N"
    f_at, f_t, f_end = _metric(prev + run)
    nxt = seq[start + length] if start + length < len(seq) else "N"
    r_at, r_t, r_end = _metric(reverse_complement(run + nxt))
    f_diff, r_diff = f_at - f_t, r_at - r_t
    if max(f_diff, r_diff) < min_tract_len:
        return None
    if f_diff >= r_diff:
        return start + (f_end - (f_at - 1) / 2.0)
    # reverse-complement branch: map the centre back to forward coordinates
    return start + (length - 1 - r_end) + (r_at - 1) / 2.0


class APhasedRepeats:
    """Detect A-phased repeats (bent DNA) like gfa ``gfa_APR.gff``.

    :param fasta_file: input FASTA.
    :param min_tract_len: minimum A-tract (A/T run) length (gfa default 3).
    :param max_tract_len: maximum A-tract length (gfa default 9).
    :param min_tracts: minimum number of phased tracts (gfa default 3).
    :param min_sep: minimum centre-to-centre tract separation (gfa 9.9).
    :param max_sep: maximum centre-to-centre tract separation (gfa 11.1).
    """

    def __init__(self, fasta_file, min_tract_len=3, max_tract_len=9, min_tracts=3, min_sep=9.9, max_sep=11.1):
        self.fasta_file = fasta_file
        self.min_tract_len = min_tract_len
        self.max_tract_len = max_tract_len
        self.min_tracts = min_tracts
        self.min_sep = min_sep
        self.max_sep = max_sep
        self.df = pd.DataFrame(columns=["seqid", "start", "end", "length", "tracts", "sequence"])

    def _find(self, seq):
        # collect valid A-tracts: (start, end, centre)
        tracts = []
        for m in _ATRUN.finditer(seq):
            length = m.end() - m.start()
            if not (self.min_tract_len <= length <= self.max_tract_len):
                continue
            center = _tract_center(seq, m.start(), length, self.min_tract_len)
            if center is not None:
                tracts.append((m.start(), m.end(), center))
        tracts.sort(key=lambda x: x[2])

        # chain tracts whose centres are one helical turn apart
        hits = []
        if not tracts:
            return hits
        chain = [tracts[0]]
        for t in tracts[1:]:
            if self.min_sep <= t[2] - chain[-1][2] <= self.max_sep:
                chain.append(t)
            else:
                if len(chain) >= self.min_tracts:
                    hits.append((chain[0][0], chain[-1][1], len(chain)))
                chain = [t]
        if len(chain) >= self.min_tracts:
            hits.append((chain[0][0], chain[-1][1], len(chain)))
        return hits

    def run(self, progress=True):
        fa = FastA(self.fasta_file)
        frames = []
        for seqid in tqdm(fa.names, disable=not progress, desc="A-phased repeats", unit="seq"):
            seq = fa.sequences[fa.names.index(seqid)].upper()
            hits = self._find(seq)
            if not hits:
                continue
            frames.append(
                pd.DataFrame(
                    {
                        "seqid": seqid,
                        "start": [h[0] for h in hits],
                        "end": [h[1] for h in hits],
                        "length": [h[1] - h[0] for h in hits],
                        "tracts": [h[2] for h in hits],
                        "sequence": [seq[h[0] : h[1]] for h in hits],
                    }
                )
            )
        cols = ["seqid", "start", "end", "length", "tracts", "sequence"]
        self.df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=cols)

    def to_bed(self, output_file, mode="w"):
        if self.df.empty:
            raise ValueError("No results. Run `.run()` first.")
        bed = self.df[["seqid", "start", "end"]].copy()
        bed["name"] = self.df["sequence"]
        bed["score"] = 0
        bed["strand"] = "+"
        bed.to_csv(output_file, sep="\t", header=False, index=False, mode=mode)

    def to_gff(self, output_file, source="sequana"):
        """Write a gff3 comparable to gfa ``gfa_APR.gff`` (1-based, inclusive)."""
        if self.df.empty:
            raise ValueError("No results. Run `.run()` first.")
        with open(output_file, "w") as fout:
            fout.write("##gff-version 3\n")
            for row in self.df.itertuples(index=False):
                start = row.start + 1  # gff is 1-based
                attrs = f"ID={row.seqid}_{start}_{row.end}_APR;tracts={row.tracts};sequence={row.sequence.lower()}"
                fout.write(f"{row.seqid}\t{source}\tA_Phased_Repeat\t{start}\t{row.end}\t.\t+\t.\t{attrs}\n")
