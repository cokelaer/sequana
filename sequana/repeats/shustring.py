#  This file is part of Sequana software
#
#  Copyright (c) 2018-2022 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Find exact repeats in a genome using the *shustring* tool."""
import subprocess

import colorlog

from sequana.fasta import FastA
from sequana.lazy import pandas as pd
from sequana.lazy import pylab

logger = colorlog.getLogger(__name__)


__all__ = ["Repeats"]


class Repeats:
    """Class for finding repeats in DNA or RNA linear sequences.

    Computation is performed each time the :attr:`threshold` is set
    to a new value.

    .. plot::
        :include-source:

        from sequana import sequana_data, Repeats
        rr = Repeats(sequana_data("measles.fa"))
        rr.threshold = 4
        rr.hist_length_repeats()

    .. note:: Works with shustring package from Bioconda (April 2017)
    .. todo:: use a specific sequence (first one by default). Others can be
        selected by name

    """

    def __init__(self, filename_fasta, merge=False, name=None):
        """.. rubric:: Constructor

        Input must be a fasta file with valid DNA or RNA characters

        :param str filename_fasta: a Fasta file, only the first
            sequence is used !
        :param int threshold: Minimal length of repeat to output
        :param str name: if name is provided, scan the Fasta file
            and select the corresponding sequence. if you want to
            analyse all sequences, you need to use a loop by setting
            _header for each sequence with the sequence name found in
            sequence header.


        .. note:: known problems. Header with a > character (e.g. in the
            comment) are left strip and only the comments is kept. Another issue
            is for multi-fasta where one sequence is ignored (last or first ?)

        """
        # used to check everything is fine with the header/name
        self._fasta = FastA(filename_fasta)

        # Define the attributes, and set the header if already provided
        self._threshold = None
        self._df_shustring = None
        self._header = None
        self._length = None
        self._longest_shustring = None
        self._begin_end_repeat_position = None
        self._begin_end_repeat_position_merge = None
        self._filename_fasta = filename_fasta
        self._previous_thr = None
        self._list_len_repeats = None
        self._contig_names = None
        if not isinstance(merge, bool):
            raise TypeError("do_merge must be boolean")
        self._do_merge = merge
        if name is not None:
            self.header = name
        else:
            self.header = self._fasta.names[0]

    def _get_header(self):
        return self._header

    def _set_header(self, name):
        if name not in self._fasta.names:
            raise ValueError("invalid name. Use one of %s" % self._fasta.names)
        self._header = name
        self._df_shustring = None

    header = property(_get_header, _set_header)

    def _get_names(self):
        if self._contig_names is None:
            self._contig_names = self._fasta.names
        return self._contig_names

    names = property(_get_names)

    def _get_shustrings_length(self):
        """Return dataframe with shortest unique substring length at each position
        shortest unique substrings are unique in the sequence and its complement
        Uses shustring tool"""
        if self._df_shustring is None:
            # read fasta
            task_read = subprocess.Popen(["cat", self._filename_fasta], stdout=subprocess.PIPE)

            # shustring command
            # the -l option uses a regular expression
            task_shus = subprocess.Popen(
                # ["shustring", "-r", "-q", "-l", ">{}[\s,\n]*?".format(self.header)],
                ["shustring", "-r", "-q", "-l", r">{}($|\s+)".format(self.header)],
                stdin=task_read.stdout,
                stdout=subprocess.PIPE,
            )

            # read stdout line by line and append to list
            list_df = []
            for line in task_shus.stdout:
                list_df.append(line.decode("utf8").replace("\n", "").split("\t"))
                # df=pd.concat([df, line])
            task_shus.wait()

            # convert to dataframe
            df = pd.DataFrame(list_df[2 : len(list_df)])
            self._df_shustring = df.astype(int)
            self._df_shustring.columns = ["position", "repeat_length"]

            # get input sequence length and longest shustring in the first line
            self._length = int(list_df[0][1])
            self._longest_shustring = int(list_df[0][3].split("<=")[2])
        return self._df_shustring

    df_shustring = property(_get_shustrings_length)

    def _get_genome_length(self):
        if self._df_shustring is None:
            self._get_shustrings_length()
        return self._length

    length = property(_get_genome_length)

    def _get_longest_shustring(self):
        if self._df_shustring is None:
            self._get_shustrings_length()
        return self._longest_shustring

    longest_shustring = property(_get_longest_shustring)

    def _find_begin_end_repeats(self, force=False):
        """Returns position of repeats longer than threshold
        as an ordered list
        """
        if self.df_shustring is None:
            self._get_shustrings_length()

        if self._threshold is None:
            # print("No threshold : please set minimul length of repeats to output")
            raise ValueError("threshold : please set threshold (minimum length of repeats to output)")

        # if there is no result yet, or the threshold has changed
        if (self._begin_end_repeat_position is None) | (self.threshold != self._previous_thr) | force:
            nb_row = self.df_shustring.shape[0]
            i = 0
            step_repeat_seq = []
            be_repeats = []
            e = 0

            # use list because faster
            list_len_shus = list(self.df_shustring.loc[:, "repeat_length"])

            while i < nb_row:
                # begining of repeat
                if list_len_shus[i] > self.threshold:
                    b = i
                    # compute new end of repeat
                    len_repeat = list_len_shus[i]
                    e = b + len_repeat
                    # save (b,e)
                    be_repeats.append((b, e))
                    # update i
                    i = e - 1
                i += 1

            self._begin_end_repeat_position = be_repeats

        self._get_merge_repeats()

    def _get_be_repeats(self):
        self._find_begin_end_repeats()
        return self._begin_end_repeat_position

    begin_end_repeat_position = property(_get_be_repeats)

    def _set_threshold(self, value):
        if value < 0:
            raise ValueError("Threshold must be positive integer")
        if value != self._threshold:
            self._previous_thr = self._threshold
        self._threshold = value
        self._find_begin_end_repeats()
        self._list_len_repeats = [tup[1] - tup[0] for tup in self._begin_end_repeat_position]

    def _get_threshold(self):
        return self._threshold

    threshold = property(_get_threshold, _set_threshold)

    def _get_list_len_repeats(self):
        if self._list_len_repeats is None:
            raise UserWarning("Please set threshold (minimum length of repeats to output)")
        return self._list_len_repeats

    list_len_repeats = property(_get_list_len_repeats)

    def _get_merge_repeats(self):
        if self._do_merge:
            # if there are repeats, merge repeats that are fragmented
            if len(self._begin_end_repeat_position) > 0:
                prev_tup = self._begin_end_repeat_position[0]
                b = prev_tup[0]
                begin_end_repeat_position_merge = []
                for i in range(1, len(self._begin_end_repeat_position)):
                    tup = self._begin_end_repeat_position[i]

                    if tup[0] == prev_tup[1]:
                        # concat
                        e = tup[1]
                        prev_tup = tup
                        if i == (len(self._begin_end_repeat_position) - 1):
                            # last tup : append to result
                            begin_end_repeat_position_merge.append((b, e))

                    else:
                        # real end of repeat : append result and update b, e
                        e = prev_tup[1]
                        begin_end_repeat_position_merge.append((b, e))
                        prev_tup = tup
                        b = prev_tup[0]
                        if i == (len(self._begin_end_repeat_position) - 1):
                            # last tup : append to result
                            begin_end_repeat_position_merge.append(tup)

                self._begin_end_repeat_position = begin_end_repeat_position_merge

    def _get_do_merge(self):
        return self._do_merge

    def _set_do_merge(self, do_merge):
        if not isinstance(do_merge, bool):
            raise TypeError("do_merge must be boolean")
        # if different
        if do_merge != self._do_merge:
            self._do_merge = do_merge
            if self._do_merge:
                # did not merge before, merge now
                if self._begin_end_repeat_position is None:
                    self._find_begin_end_repeats()
            else:
                # data is already merged : need to compute again to un-merge
                self._find_begin_end_repeats(force=True)

    do_merge = property(_get_do_merge, _set_do_merge)

    def hist_length_repeats(
        self,
        bins=20,
        alpha=0.5,
        hold=False,
        fontsize=12,
        grid=True,
        title="Repeat length",
        xlabel="Repeat length",
        ylabel="#",
        logy=True,
    ):
        """Plots histogram of the repeat lengths"""
        # check that user has set a threshold
        if hold is False:
            pylab.clf()
        pylab.hist(self.list_len_repeats, alpha=alpha, bins=bins)
        pylab.title(title)
        pylab.xlabel(xlabel, fontsize=fontsize)
        pylab.ylabel(ylabel, fontsize=fontsize)
        if grid is True:
            pylab.grid(True)
        if logy:
            pylab.semilogy()

    def plot(self, clf=True, fontsize=12):
        if clf:
            pylab.clf()
        pylab.grid(True)
        pylab.plot(self.df_shustring.repeat_length)
        pylab.xlabel("Position (bp)", fontsize=fontsize)
        pylab.ylabel("Repeat lengths (bp)", fontsize=fontsize)
        pylab.ylim(bottom=0)

    def to_wig(self, filename, step=1000):
        """export repeats into WIG format to import in IGV"""
        assert self.threshold and self.df_shustring is not None

        N = len(self.df_shustring)
        with open(filename, "w") as fout:
            fout.write(
                f'track type=wiggle_0 name="Repeat Density" visibility=full fixedStep chrom=1 start=0 step={step} span={step}\n'
            )
            # min_period = 1 to prevent NAN at first position
            rolling_max = self.df_shustring.rolling(step, step=step).max()
            for _, row in rolling_max.iterrows():
                if row.name == 0:
                    continue
                M = row.repeat_length
                start = row.name - step
                stop = row.name
                fout.write(f"{self.header}\t{start}\t{stop}\t{row.repeat_length}\n")

    def get_peak_position_and_length(self, THRESHOLD=3000):

        df = self.df_shustring.copy()
        # Filter for repeats above threshold
        df = df[df["repeat_length"] > THRESHOLD].copy()

        # Sort by genomic position
        df = df.sort_values(by="position").reset_index(drop=True)

        # Identify local maxima: where the value is greater than the left and right neighbors
        df["prev"] = df["repeat_length"].shift(1, fill_value=0)
        df["next"] = df["repeat_length"].shift(-1, fill_value=0)

        # Keep only local maxima
        peaks = df[(df["repeat_length"] > df["prev"]) & (df["repeat_length"] > df["next"])]

        # Select only relevant columns
        final_df = peaks[["position", "repeat_length"]]

        # Print results
        print(f"Final number of peaks after filtering: {len(final_df)}")
        return final_df
