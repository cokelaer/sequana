#
#  This file is part of Sequana software
#
#  Copyright (c) 2025 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Multiple sequence alignment parsing and manipulation."""

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import colorlog

from sequana.lazy import numpy as np

logger = colorlog.getLogger(__name__)

__all__ = ["Alignment", "parse_phylip", "parse_stockholm", "parse_nexus"]


@dataclass
class Alignment:
    """Multiple sequence alignment.

    Attributes:
        sequences: dict mapping sequence name -> sequence string
        names: ordered list of sequence names
        format: source format (phylip, stockholm, nexus, fasta)
        annotations: dict of per-sequence metadata
    """

    sequences: Dict[str, str] = field(default_factory=dict)
    format: str = "unknown"
    annotations: Dict[str, Dict] = field(default_factory=dict)

    @property
    def names(self) -> List[str]:
        """Return sequence names in order."""
        return list(self.sequences.keys())

    def __len__(self) -> int:
        """Return number of sequences."""
        return len(self.sequences)

    def length(self) -> int:
        """Return alignment length (sequence length)."""
        if not self.sequences:
            return 0
        return len(next(iter(self.sequences.values())))

    def add_sequence(self, name: str, sequence: str, annotations: Optional[Dict] = None) -> None:
        """Add sequence to alignment.

        Args:
            name: sequence identifier
            sequence: aligned sequence string
            annotations: optional metadata dict
        """
        if name in self.sequences and self.sequences[name] != sequence:
            logger.warning(f"Overwriting sequence {name}")

        self.sequences[name] = sequence
        if annotations:
            self.annotations[name] = annotations

    def get_column(self, pos: int) -> Dict[str, str]:
        """Get all characters at alignment position.

        Args:
            pos: column position (0-indexed)

        Returns:
            dict mapping sequence name -> character
        """
        if pos >= self.length():
            return {}

        return {name: seq[pos] for name, seq in self.sequences.items()}

    def consensus(self, ignore_gaps: bool = True) -> str:
        """Return consensus sequence (most common char per position).

        Args:
            ignore_gaps: exclude gaps/Ns from consensus calculation

        Returns:
            consensus sequence string
        """
        if not self.sequences:
            return ""

        consensus = []
        for pos in range(self.length()):
            col = self.get_column(pos)
            chars = list(col.values())

            if ignore_gaps:
                chars = [c for c in chars if c not in ["-", "N", "n", "."]]

            if not chars:
                consensus.append("-")
            else:
                # Most common character
                most_common = max(set(chars), key=chars.count)
                consensus.append(most_common)

        return "".join(consensus)

    def write_phylip(self, filename: str, interleaved: bool = False) -> None:
        """Write alignment to Phylip format.

        Args:
            filename: output file path
            interleaved: if True, use interleaved format; else sequential
        """
        with open(filename, "w") as f:
            n_seqs = len(self.sequences)
            seq_len = self.length()
            f.write(f"{n_seqs} {seq_len}\n")

            if interleaved:
                self._write_phylip_interleaved(f)
            else:
                self._write_phylip_sequential(f)

    def _write_phylip_sequential(self, f) -> None:
        """Write sequential Phylip format."""
        for name, seq in self.sequences.items():
            # Phylip names: max 10 chars, padded
            name_padded = name[:10].ljust(10)
            f.write(f"{name_padded}{seq}\n")

    def _write_phylip_interleaved(self, f) -> None:
        """Write interleaved Phylip format."""
        chunk_size = 60
        seq_len = self.length()

        for start in range(0, seq_len, chunk_size):
            end = min(start + chunk_size, seq_len)
            if start > 0:
                f.write("\n")

            for name, seq in self.sequences.items():
                name_padded = name[:10].ljust(10)
                chunk = seq[start:end]
                f.write(f"{name_padded}{chunk}\n")

    def write_stockholm(self, filename: str) -> None:
        """Write alignment to Stockholm format.

        Includes sequence features/annotations if available.
        """
        with open(filename, "w") as f:
            f.write("# STOCKHOLM 1.0\n")

            for name, seq in self.sequences.items():
                f.write(f"{name}\t{seq}\n")

            # Write annotations if present
            for name in self.sequences:
                if name in self.annotations:
                    annot = self.annotations[name]
                    for key, val in annot.items():
                        f.write(f"#=GS {name} {key} {val}\n")

            f.write("//\n")

    def to_fasta(self) -> str:
        """Return FASTA format string."""
        lines = []
        for name, seq in self.sequences.items():
            lines.append(f">{name}")
            lines.append(seq)
        return "\n".join(lines)

    def stats(self) -> dict:
        """Return alignment statistics."""
        return {
            "num_sequences": len(self.sequences),
            "length": self.length(),
            "format": self.format,
            "names": self.names,
        }


class PhylipParser:
    """Parse Phylip format alignments."""

    @staticmethod
    def parse(filename: str) -> Alignment:
        """Parse Phylip file (sequential or interleaved).

        Args:
            filename: path to Phylip file

        Returns:
            Alignment object
        """
        with open(filename) as f:
            lines = f.readlines()

        return PhylipParser.parse_string("".join(lines))

    @staticmethod
    def parse_string(content: str) -> Alignment:
        """Parse Phylip format from string.

        Args:
            content: Phylip format content

        Returns:
            Alignment object
        """
        lines = content.strip().split("\n")
        if not lines:
            return Alignment(format="phylip")

        # Parse header: "num_seqs seq_length"
        header = lines[0].split()
        try:
            n_seqs = int(header[0])
            seq_len = int(header[1])
        except (ValueError, IndexError):
            raise ValueError("Invalid Phylip header format")

        alignment = Alignment(format="phylip")

        # Detect format (sequential vs interleaved)
        # Sequential: n_seqs lines of sequences (after header)
        # Interleaved: multiple blocks of n_seqs lines
        data_lines = [l for l in lines[1:] if l.strip()]  # Skip header, skip empty

        if len(data_lines) > n_seqs:
            # More data lines than sequences -> interleaved
            return PhylipParser._parse_interleaved(lines[1:], n_seqs, seq_len)
        else:
            # Exactly n_seqs lines -> sequential
            return PhylipParser._parse_sequential(lines[1:], n_seqs)

    @staticmethod
    def _parse_sequential(lines: List[str], n_seqs: int) -> Alignment:
        """Parse sequential Phylip format."""
        alignment = Alignment(format="phylip")

        for i in range(min(n_seqs, len(lines))):
            line = lines[i]
            parts = line.split(None, 1)  # Split on first whitespace
            if len(parts) == 2:
                name, seq = parts
                seq = seq.replace(" ", "").replace("\n", "")
                alignment.add_sequence(name, seq)

        return alignment

    @staticmethod
    def _parse_interleaved(lines: List[str], n_seqs: int, seq_len: int) -> Alignment:
        """Parse interleaved Phylip format."""
        alignment = Alignment(format="phylip")
        seqs = {}

        i = 0
        while i < len(lines):
            line = lines[i].strip()

            # Skip empty lines
            if not line:
                i += 1
                continue

            # Parse sequence line
            parts = line.split(None, 1)
            if len(parts) == 2:
                name, chunk = parts
                chunk = chunk.replace(" ", "").replace("\n", "")
                if name not in seqs:
                    seqs[name] = ""
                seqs[name] += chunk

            i += 1

        for name, seq in seqs.items():
            alignment.add_sequence(name, seq)

        return alignment


class StockholmParser:
    """Parse Stockholm format alignments."""

    @staticmethod
    def parse(filename: str) -> Alignment:
        """Parse Stockholm file.

        Args:
            filename: path to Stockholm file

        Returns:
            Alignment object
        """
        with open(filename) as f:
            content = f.read()
        return StockholmParser.parse_string(content)

    @staticmethod
    def parse_string(content: str) -> Alignment:
        """Parse Stockholm format from string.

        Args:
            content: Stockholm format content

        Returns:
            Alignment object
        """
        alignment = Alignment(format="stockholm")
        lines = content.split("\n")

        in_alignment = False
        for line in lines:
            line = line.rstrip()

            if line.startswith("# STOCKHOLM"):
                in_alignment = True
                continue

            if line == "//":
                in_alignment = False
                continue

            if not in_alignment or not line.strip() or line.startswith("#"):
                continue

            # Parse sequence line: "name sequence"
            parts = line.split(None, 1)
            if len(parts) == 2 and not parts[0].startswith("#"):
                name, seq = parts
                if name not in alignment.sequences:
                    alignment.add_sequence(name, seq)
                else:
                    # Append to existing sequence (multi-line)
                    alignment.sequences[name] += seq

        return alignment


class NexusParser:
    """Parse Nexus format alignments."""

    @staticmethod
    def parse(filename: str) -> Alignment:
        """Parse Nexus file.

        Args:
            filename: path to Nexus file

        Returns:
            Alignment object
        """
        with open(filename) as f:
            content = f.read()
        return NexusParser.parse_string(content)

    @staticmethod
    def parse_string(content: str) -> Alignment:
        """Parse Nexus format from string.

        Args:
            content: Nexus format content

        Returns:
            Alignment object
        """
        alignment = Alignment(format="nexus")
        lines = content.split("\n")

        in_data = False
        for i, line in enumerate(lines):
            line_lower = line.lower().strip()

            if line_lower.startswith("begin data"):
                in_data = True
                continue

            if line_lower.startswith("end;") or line_lower.startswith("endblock;"):
                in_data = False
                continue

            if not in_data:
                continue

            # Parse matrix block
            if line_lower.startswith("matrix"):
                j = i + 1
                while j < len(lines):
                    seq_line = lines[j].strip()
                    if seq_line.startswith(";") or not seq_line:
                        break
                    if seq_line.lower() in ["end;", "endblock;"]:
                        break

                    parts = seq_line.split(None, 1)
                    if len(parts) == 2:
                        name, seq = parts
                        seq = seq.replace(" ", "").replace("\n", "")
                        if name not in alignment.sequences:
                            alignment.add_sequence(name, seq)
                        else:
                            alignment.sequences[name] += seq
                    j += 1

        return alignment


def parse_phylip(filename: str) -> Alignment:
    """Parse Phylip format file."""
    return PhylipParser.parse(filename)


def parse_stockholm(filename: str) -> Alignment:
    """Parse Stockholm format file."""
    return StockholmParser.parse(filename)


def parse_nexus(filename: str) -> Alignment:
    """Parse Nexus format file."""
    return NexusParser.parse(filename)


def auto_detect_format(filename: str) -> str:
    """Detect alignment format from file content.

    Returns: 'phylip', 'stockholm', 'nexus', 'fasta', or 'unknown'
    """
    with open(filename) as f:
        content = f.read(1000)  # Read first 1000 chars

    content_lower = content.lower()

    if "stockholm" in content_lower:
        return "stockholm"
    if "begin data" in content_lower or "nexus" in content_lower:
        return "nexus"
    if ">" in content:
        return "fasta"

    # Try to detect Phylip: first line is "num_seqs seq_length"
    lines = content.split("\n")
    if lines:
        parts = lines[0].split()
        if len(parts) == 2:
            try:
                int(parts[0])
                int(parts[1])
                return "phylip"
            except ValueError:
                pass

    return "unknown"


def parse_alignment(filename: str, format: Optional[str] = None) -> Alignment:
    """Parse alignment file with automatic or specified format.

    Args:
        filename: path to alignment file
        format: format name ('phylip', 'stockholm', 'nexus', 'fasta')
                if None, auto-detect

    Returns:
        Alignment object
    """
    if format is None:
        format = auto_detect_format(filename)

    if format == "phylip":
        return parse_phylip(filename)
    elif format == "stockholm":
        return parse_stockholm(filename)
    elif format == "nexus":
        return parse_nexus(filename)
    elif format == "fasta":
        from sequana.fasta import FastA

        fasta = FastA(filename)
        alignment = Alignment(format="fasta")
        for record in fasta:
            alignment.add_sequence(record.name, record.sequence)
        return alignment
    else:
        raise ValueError(f"Unknown format: {format}")
