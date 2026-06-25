"""Reverse-complement palindrome detection.

A DNA *palindrome* is a sequence equal to its own reverse complement (e.g.
``GAATTC``). Such palindromes always have an even length: the central base of
an odd-length window would have to be its own complement, which never happens
for A/C/G/T. Only even window sizes are therefore scanned.
"""
import pandas as pd
from tqdm import tqdm

from sequana import FastA
from sequana.tools import reverse_complement

_COLS = ["seqid", "start", "end", "length", "sequence"]


class Palindromes:
    def __init__(self, fasta_file, min_len=4, max_len=12):
        self.fasta_file = fasta_file
        self.min_len = min_len
        self.max_len = max_len
        self.df = pd.DataFrame(columns=_COLS)

    def is_palindrome(self, seq):
        # Upper-case first: the complement table preserves case, so a
        # soft-masked (lower-case) palindrome would otherwise be missed.
        seq = seq.upper()
        return seq == reverse_complement(seq)

    def run(self):
        fa = FastA(self.fasta_file)

        results = []
        for name, sequence in zip(tqdm(fa.names), fa.sequences):
            sequence = sequence.upper()
            seq_len = len(sequence)

            # Even sizes only: odd-length windows can never be palindromes.
            start = self.min_len + (self.min_len % 2)
            for size in range(start, self.max_len + 1, 2):
                for i in range(seq_len - size + 1):
                    subseq = sequence[i : i + size]
                    if subseq == reverse_complement(subseq):
                        results.append({"seqid": name, "start": i, "end": i + size, "length": size, "sequence": subseq})

        self.df = pd.DataFrame(results, columns=_COLS)

    def to_bed(self, output_file):
        if self.df.empty:
            raise ValueError("No data. Run `.run()` first.")
        bed = self.df[["seqid", "start", "end"]].copy()
        bed["name"] = self.df["sequence"]
        bed["score"] = 0
        bed["strand"] = "+"
        bed.to_csv(output_file, sep="\t", header=False, index=False)
