#
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
"IEM class"
import collections
import io
import os
import sys

import colorlog

from sequana.lazy import pandas as pd

logger = colorlog.getLogger(__name__)


__all__ = ["SampleSheet", "IEM", "BCLConvert", "get_sample_sheet_version", "SampleSheetFactory"]


class SampleSheet:
    """Reader and validator of Illumina samplesheets

    The Illumina samplesheet reader and validator verifies the correctness of the sections
    in the samplesheet, which are not case-sensitive and are enclosed within square brackets.

    Following the closing bracket, no additional characters are permitted, except
    for commas and the end-of-line marker. For instance [Data]a prevents the [Data] section
    from being correctly processed.

    The sections then consist of key-value pairs represented as records, with each line consisting
    of precisely two fields.

    An optional [Settings] section can contain key-value pairs, and
    the [Reads] section specifies the number of cycles per read, which
    is exclusively required for MiSeq.

    The [Data] section, which is a table similar to CSV format, is optional.
    However, without [Data] section all reads are sent to a single 'undetermined'
    output file. Sample_ID is highly recommended.

    Example of typical Data section to be used with bcl2fastq::

        [Header]

        [Data]
        Sample_ID,Sample_Name,I7_Index_ID,index,I5_INdex_ID,index2
        A10001,Sample_A,D701,AATACTCG,D501,TATAGCCT
        A10002,Sample_B,D702,TCCGGAGA,D501,TATAGCCT
        A10003,Sample_C,D703,CGCTCATT,D501,TATAGCCT
        A10004,Sample_D,D704,GAGATTCC,D501,TATAGCCT

    Important: altough we have upper case names as specified in the Illumina specs, the
    bcl2fastq does not care about the upper case. This is not intuitive since IEM produces
    keys with upper and lower case names similarly to the specs.

    **Sequana Standalone**

    The standalone application **sequana** contains a subcommand based on this class::

        sequana samplesheet

    that can be used to check the correctness of a samplesheet::

        sequana samplesheet --check SampleSheet.csv

    :references: illumina specifications 970-2017-004.pdf
    """

    expected_headers_fields = [
        "IEMFileVersion",
        "Investigator Name",
        "Instrument Type",
        "Experiment Name",
        "Date",
        "Workflow",
        "Application",
        "Assay",
        "Description",
        "Chemistry",
        "Index Adapters",
    ]

    expected_data_headers = {"SE": [], "PE": []}

    #: name of the section holding the tabular sample data (title-cased). Overridden in
    #: :class:`BCLConvert` where the section is ``[BCLConvert_Data]``.
    data_section = "Data"

    #: name of the settings section (title-cased). Overridden in :class:`BCLConvert`
    #: where the required section is ``[BCLConvert_Settings]``.
    settings_section = "Settings"

    def __init__(self, filename):

        self.filename = filename
        if not os.path.exists(self.filename):
            raise IOError(f"{filename} does not exist")
        # figures out the sections in the sample sheet.
        # we use a try/except so that even in case of failure, we can still use
        # quickfix or attributes.
        try:
            self._scan_sections()
        except Exception as err:  # pragma: no cover
            print(err)

    def _line_cleaner(self, line, line_count):
        # We can get rid of EOL and spaces
        line = line.strip()

        # is it an empty line ?
        if len(line) == 0:
            return line

        # if we are dealing with a section title, we can cleanup the
        # line. A section must start with '[' and ends with ']' but
        # there could be spaces and commands. If it ends with a ; then
        # the section will not be found as expected since this is
        # sympatomatic of further issues
        if line.startswith("["):
            # [Header], ,, ,\n becomes [Header]
            line = line.strip(", ;")  # note the space AND comma

        return line

    def _scan_sections(self):
        # looks for special section Header/Data/Reads/Settings or
        # any other matching section thst looks like [XXX]

        # sections can be of any types of cases (lower, upper, mixed)

        current_section = None
        data = {}
        with open(self.filename, "r") as fin:
            for line_count, line in enumerate(fin.readlines()):
                line = self._line_cleaner(line, line_count + 1)
                if len(line) == 0:
                    continue
                if line.startswith("[") and line.endswith("]"):
                    name = line.lstrip("[").rstrip("]")
                    current_section = name
                    data[current_section] = []  # create empty list just to create the section
                else:
                    if current_section in data:
                        data[current_section] += [line]
                    else:
                        data[current_section] = [line]
        self.sections = data

        # if the sample sheet starts with an incorrect name
        # e.g section not closed by square brackets or empty lines
        # the a None key is present and should be ignored or cast
        # into a string:
        if None in self.sections:
            self.sections["None"] = self.sections[None]
            del self.sections[None]

        # Some cleanup. Since the sections are case insensitive within
        # bcl2fastq, we need to convert everything back to titles
        self.sections = {k.title(): v for k, v in self.sections.items()}

    def _get_df(self):
        if self.data_section in self.sections:
            if self.sections[self.data_section]:
                # cope with the case of comma or semicolon separators.
                df1 = pd.read_csv(io.StringIO("\n".join(self.sections[self.data_section])), index_col=False, sep=",")
                df2 = pd.read_csv(io.StringIO("\n".join(self.sections[self.data_section])), index_col=False, sep=";")

                if len(df1.columns) > len(df2.columns):

                    # rename all columns to lower case
                    df1.columns = [x.lower() for x in df1.columns]
                    return df1
                else:
                    df2.columns = [x.lower() for x in df2.columns]
                    return df2
            else:
                return pd.DataFrame()
        else:  # pragma: no cover
            return pd.DataFrame()

    df = property(_get_df, doc="Returns the [Data] section")

    def _get_samples(self):
        try:
            return self.df["sample_id"].values
        except AttributeError:  # pragma: no cover
            return "No Sample_ID found in the Data header section"

    samples = property(_get_samples, doc="returns the sample identifiers as a list")

    def _get_version(self):
        try:
            return self.header["IEMFileVersion"]
        except KeyError:
            return None

    version = property(_get_version, doc="return the version of the IEM file")

    def checker(self):

        from sequana.utils.checker import Checker

        checks = Checker()

        if "Reads" in self.sections:
            pass
        else:
            checks.results.append({"msg": f"The optional [Reads] section is missing", "status": "Warning"})

        if "Header" in self.sections:
            pass
        else:
            checks.results.append({"msg": f"The optional [Header] section is missing", "status": "Warning"})

        if "Data" in self.sections:
            checks.tryme(self._check_data_section_csv_format)
            checks.tryme(self._check_optional_data_column_names)
            checks.tryme(self._check_mandatory_data_columns)
            checks.tryme(self._check_sample_ID)
            checks.tryme(self._check_sample_project)
            checks.tryme(self._check_unique_sample_ID)
            checks.tryme(self._check_unique_indices)
            checks.tryme(self._check_nucleotide_indices)
            checks.tryme(self._check_sample_lane_number)
            checks.tryme(self._check_alpha_numerical)

            if "sample_name" in self.df.columns:
                checks.tryme(self._check_sample_names)
                checks.tryme(self._check_unique_sample_name)
            else:
                checks.results.append(
                    {
                        "name": "check_sample_names",
                        "msg": f"Column Sample_Name not found in the header of the [Data] section. Recommended.",
                        "status": "Warning",
                    }
                )

            # if data section is incorrect, self.df is not accessible so the following validation
            # will fail.
            try:
                if "index" in self.df.columns:
                    checks.tryme(self._check_homogene_I7_length)

                if "index2" in self.df.columns:
                    checks.tryme(self._check_homogene_I5_length)
                    checks.tryme(self._check_homogene_I5_and_I7_length)
            except pd.errors.ParserError:  # pragma: no cover
                pass
        else:
            checks.results.append({"msg": f"The [Data] section is missing. ", "status": "Error"})

        if "Settings" in self.sections:
            checks.tryme(self._check_settings)
        else:
            checks.results.append({"msg": f"The optional [Settings] section is missing", "status": "Warning"})

        checks.tryme(self._check_semi_column_presence)

        return checks.results

    def _check_settings(self):

        # Detect malformed lines that are not valid 'key,value' pairs. _get_settings
        # tolerates them so downstream callers do not crash; here we surface a clear
        # error instead of letting a parsing exception leak (e.g. sheets mangled by
        # spreadsheet software into ';'-separated values or with stray ';;;').
        for line in self.sections.get(self.settings_section, []):
            if line.strip() and "," not in line:
                return {
                    "name": "check_settings",
                    "msg": (
                        f"The [Settings] section has a line that is not a valid "
                        f"'key,value' pair: '{line}'. This is often caused by stray "
                        f"';' separators (e.g. added by spreadsheet software)."
                    ),
                    "status": "Error",
                }

        for k, v in self.settings.items():

            # checks ACGT content (no acgt allowed)
            if k.lower() in [
                x.lower()
                for x in [
                    "Adapter",
                    "AdapterRead1",
                    "TrimAdapter",
                    "AdapterRead2",
                    "TrimAdapterRead2",
                    "MaskAdapter",
                    "MaskAdapterRead2",
                ]
            ]:
                allowed_chars = set("ACGT")

                def is_valid_string(s):
                    return set(s).issubset(allowed_chars)

                if is_valid_string(v) is False:
                    return {
                        "name": "check_settings",
                        "msg": f"Invalid nucleotide sequence found for {k} (v)",
                        "status": "Error",
                    }
            elif k.lower() in [
                x.lower()
                for x in ["ReverseComplement", "FindAdaptersWithIndels", "TrimUMI", "CreateFastqForIndexReads"]
            ]:
                if v not in ("true", "false", "t", "f", "yes", "no", "y", "n", "1", "0"):
                    return {
                        "name": "check_settings",
                        "msg": f"Invalid valid for {k} ({v}). Must be set to one of : true, false, t, f, yes, no, y, n, 1, 0",
                        "status": "Error",
                    }

            elif k.lower() in [
                x.lower()
                for x in [
                    "Read1StartFromCycle",
                    "Read2StartFromCycle",
                    "Read1EndWithCycle",
                    "Read2EndWithCycle",
                    "Read1UMILength",
                    "Read2UMILength",
                    "Read1UMIStartFromCycle",
                    "Read2UMIStartFromCycle",
                ]
            ]:
                allowed_chars = set("0123456789")

                def is_valid_int(s):
                    return set(s).issubset(allowed_chars)

                if is_valid_int(v) is False or int(v) < 0:
                    return {
                        "name": "check_settings",
                        "msg": f"Invalid value for {k}. Must be positive. You provided: {v}",
                        "status": "Error",
                    }
            elif k.lower() in [x.lower() for x in ["ExcludeTiles", "ExcludeTilesLaneX"]]:

                allowed_chars = set("0123456789")

                def is_valid_int(s):
                    return set(s).issubset(allowed_chars)

                # valid is 1101+1102+1103-1110  (only +-,numbers)
                for item in v.split("+"):
                    for x in item.split("-"):
                        if is_valid_int(x) is False:
                            return {
                                "name": "check_settings",
                                "msg": f"Invalid value for {k} (v). Must be made of integers, + and - signs. e.g. 1101+1105-1110 to exclude 1101 and values in [1105-1110].",
                                "status": "Error",
                            }
        else:
            return {"name": "check_settings", "msg": f"The [Settings] section looks good", "status": "Success"}

    def _check_sample_ID(self):

        if "sample_id" in self.df.columns:  # optional
            # check that names are not in 'all' or 'undetermined'
            if (self.df["sample_id"] == "unknown").sum() or (self.df["sample_id"] == "all").sum():
                return {
                    "name": "check_sample_ID",
                    "msg": "Sample_ID column contains forbidden name ('all' or 'unknown')",
                    "status": "Error",
                }
            else:
                return {
                    "name": "check_sample_ID",
                    "msg": "Sample_ID column (no unknown/undetermined label). Looks correct",
                    "status": "Success",
                }

        else:
            return {
                "name": "check_sample_ID",
                "msg": f"Column Sample_ID not found in the header of the [Data] section. All data will be stored in Undetermined.fastq.gz",
                "status": "Warning",
            }

    def _check_sample_names(self):

        # check that names are not in 'all' or 'undetermined'
        if (self.df["sample_name"] == "unknown").sum() or (self.df["sample_name"] == "undetermined").sum():
            return {
                "name": "check_sample_names",
                "msg": "Sample_Name column contains forbidden name ('all' or 'undetermined')",
                "status": "Error",
            }
        else:
            return {
                "name": "check_sample_names",
                "msg": "Sample_Name column (no unknown/undetermined label). Looks correct",
                "status": "Success",
            }

    def _check_sample_project(self):

        if "sample_project" in self.df.columns:  # optional
            # check that names are not in 'all' or 'undetermined'
            if (self.df["sample_project"] == "all").sum() or (self.df["sample_project"] == "default").sum():
                return {
                    "name": "check_sample_project",
                    "msg": "Sample_Project column contains forbidden name ('all' or 'default')",
                    "status": "Error",
                }
            else:
                return {
                    "name": "check_sample_project",
                    "msg": "Sample_Project column (no all/default label). Looks correct",
                    "status": "Success",
                }
        else:
            return {
                "name": "check_sample_project",
                "msg": f"Column Sample_Project not found in the header of the [Data] section. Recommended.",
                "status": "Warning",
            }

    def _check_mandatory_data_columns(self):
        # In fact, all columns are optional except index and Sample_Name is optional
        # Sample_Project is optional. If provided, fastq are saved in that sub directory.

        for column in ["sample_id", "index"]:
            if column not in self.df.columns:
                return {
                    "name": "check_mandatory_data_columns",
                    "msg": f"Mandatory '{column}' column not found in the header of the [Data] section",
                    "status": "Error",
                }
        return {
            "name": "check_mandatory_data_columns",
            "msg": f"Mandatory columns (index, Sample_ID) found in the header of the [Data] section",
            "status": "Success",
        }

    def _check_unique_sample_name(self):
        # check that sample names are unique and that sample Names are unique too

        if self.df["sample_name"].isnull().sum() > 0:
            return {"name": "check_unique_sample_name", "msg": "Some sample names are empty", "status": "Warning"}

        elif len(self.df.sample_name) != len(self.df.sample_name.unique()):
            duplicated = self.df.sample_name[self.df.sample_name.duplicated()].index
            duplicated = ",".join([str(x + 1) for x in duplicated])
            return {
                "name": "check_unique_sample_name",
                "msg": f"Sample_Name not unique. Duplicated entries on lines: {duplicated}",
                "status": "Warning",
            }
        else:
            return {"name": "check_unique_sample_name", "msg": "Sample name uniqueness", "status": "Success"}

    def _check_unique_sample_ID(self):
        # check that sample names are unique and that sample Names are unique too

        if "sample_id" not in self.df.columns:
            return {
                "name": "check_unique_sample_ID",
                "msg": "Sample ID not found in the header of the [Data] section",
                "status": "Warning",
            }

        if len(self.df["sample_id"]) != len(self.df["sample_id"].unique()):
            duplicated = self.df.sample_id[self.df.sample_id.duplicated()].index
            duplicated = ",".join([str(x + 1) for x in duplicated])
            return {
                "name": "check_unique_sample_ID",
                "msg": f"Sample ID not unique. Duplicated entries on lines: {duplicated}",
                "status": "Error",
            }
        else:
            return {"name": "check_unique_sample_ID", "msg": "Sample ID uniqueness", "status": "Success"}

    def _check_sample_lane_number(self):
        if "sample_lane" in self.df.columns:

            # Define the allowed lanes
            allowed_chars = set("12345678")

            def is_valid_lane(s):
                return set(str(s)).issubset(allowed_chars)

            # Apply the function to the DataFrame
            invalid_lanes = list(self.df[~self.df["sample_lane"].apply(is_valid_lane)].index)

            if len(invalid_lanes):
                invalid_lanes = [x + 1 for x in invalid_lanes]
                return {
                    "name": "check_sample_lane_number",
                    "msg": f"Incorrect lane number in these rows: {invalid_lanes}. Must be in the range 1-8",
                    "status": "Error",
                }

        return {"name": "check_sample_lane_number", "msg": "Correct lane number range", "status": "Success"}

    def _check_nucleotide_indices(self):

        # Define the allowed characters
        allowed_chars = set("ACGTN")

        def is_valid_string(s):
            try:
                return set(s).issubset(allowed_chars)
            except:
                return False

        # Apply the function to the DataFrame
        invalid_rows = list(self.df[~self.df["index"].apply(is_valid_string)].index)

        if "index2" in self.df.columns:
            invalid_rows += list(self.df[~self.df["index2"].apply(is_valid_string)].index)
        if len(invalid_rows):
            invalid_rows = [x + 1 for x in invalid_rows]
            return {
                "name": "check_nucleotide_indices",
                "msg": f"these rows have invalid index with wrong nucleotides {invalid_rows}",
                "status": "Error",
            }

        return {"msg": "Indices are made of A, C, G, T, N.", "status": "Success"}

    def _check_unique_indices(self):
        if "index2" in self.df.columns:
            indices = self.df["index"] + "," + self.df["index2"]
            msg = "You have duplicated index I7/I5."
        elif "index" in self.df.columns:
            indices = self.df["index"]
            msg = "You have duplicated index I7."
        else:
            return {
                "name": "check_unique_indices",
                "msg": f"column 'index' not found in the header of the [Data] section.",
                "status": "Error",
            }

        if indices.duplicated().sum() > 0:
            duplicated = indices[indices.duplicated()].values
            try:
                IDs = self.df[indices.duplicated()].sample_id.values
                IDs = ", ".join([str(x) for x in IDs])
                IDs = f"related to sample IDs: {IDs}"
            except Exception as err:  # pragma: no cover
                IDs = ""

            return {"name": "check_unique_indices", "msg": f"{msg} {duplicated} {IDs}", "status": "Error"}
        else:
            return {"name": "check_unique_indices", "msg": "Indices are unique.", "status": "Success"}

    def _check_optional_data_column_names(self):
        msg = ""
        warnings = []
        # check whether minimal columns are included
        for x in [
            "i7_index_id",
            "sample_project",
            "description",
            "sample_plate",
            "sample_well",
            "lane",
            "index_plate",
            "index_plate_well",
        ]:
            if x not in self.df.columns:
                warnings.append(x)

        if len(warnings):
            warnings = ",".join(warnings)
            msg = f"Some columns are missing in the [Data] section: {warnings}"
            return {"msg": msg, "status": "Warning"}
        else:  # pragma: no cover
            return {"msg": "Columns of the data section looks good", "status": "Success"}

    def _get_data_length(self):
        N = len(set([x.count(",") for x in self.sections[self.data_section]]))
        if len(self.sections[self.data_section]) >= 2 and N == 1:
            return True
        else:
            return False

    def _check_homogene_I7_length(self):
        if self._get_data_length():
            if len(self.df) == 1:
                L = len(self.df["index"].values[0].strip())
                if L:
                    return {"msg": f"Only one sample. Index length is {L}", "status": "Success"}
                else:
                    return {"msg": f"Only one sample. Index length is {L}", "status": "Error"}

            diff = self.df["index"].apply(lambda x: len(x)).std()
            if diff == 0:
                return {"msg": "Indices length in I7 have same lengths", "status": "Success"}
            else:
                return {"msg": "Indices length in I7 have different lengths", "status": "Error"}
        else:
            return {"msg": "Indices length could not be read.", "status": "Warning"}

    def _check_homogene_I5_length(self):
        if self._get_data_length():

            if len(self.df) == 1:
                L = len(self.df["index2"].values[0].strip())
                if L:
                    return {"msg": f"Only one sample. Index length for I5 is {L}", "status": "Success"}
                else:
                    return {"msg": f"Only one sample. Index length for I5 is {L}", "status": "Error"}
            diff = self.df["index2"].apply(lambda x: len(x)).std()
            if diff == 0:
                return {"msg": "Indices length in I5 have same lengths", "status": "Success"}
            else:
                return {"msg": "Indices length in I5 have different lengths", "status": "Error"}
        else:
            return {"msg": "Indices length could not be read. ", "status": "Warning"}

    def _check_homogene_I5_and_I7_length(self):
        # this is a warning only since you may have custom index
        lengths = [len(x) for x in self.df["index"]] + [len(x) for x in self.df["index2"]]
        lengths = list(set(lengths))
        L = lengths[0]

        if len(lengths) == 1:
            return {"msg": f"I5 and I7 have coherent length of {L}", "status": "Success"}
        else:
            return {"msg": "I5 and I7 have different lengths", "status": "Warning"}

        return {"msg": "Indices length could not be read. ", "status": "Success"}

    def _check_csv_format(self, section, name, level="Error"):
        # Generic CSV-consistency check for a tabular section. ``level`` is the
        # status used when the section is malformed (Error for the mandatory data
        # section, Warning for auxiliary sections such as [Cloud_Data]).
        N = len(set([x.count(",") for x in self.sections[section]]))

        if N == 1:  # looks correct
            if len(self.sections[section]) == 1:
                return {
                    "name": name,
                    "msg": f"The [{section}] section CSV format looks empty. Remove if of fill it.",
                    "status": level,
                }
            else:
                return {
                    "name": name,
                    "msg": f"The [{section}] section CSV format looks correct.",
                    "status": "Success",
                }
        elif N == 0:
            return {
                "name": name,
                "msg": f"The [{section}] section CSV format looks empty. Remove it of fill it",
                "status": level,
            }
        else:
            lengths = set([x.count(",") for x in self.sections[section]])
            return {
                "name": name,
                "msg": f"The [{section}] section has lines with different number of entries {lengths}. Probably missing or extra commas in the [{section}] section.",
                "status": level,
            }

    def _check_data_section_csv_format(self):
        return self._check_csv_format(self.data_section, "check_data_section_csv_format", level="Error")

    def _check_semi_column_presence(self):
        with open(self.filename, "r") as fp:
            line_count = 1
            for line in fp.readlines():
                if line.rstrip().endswith(";"):
                    return {
                        "name": "check_semi_column_presence",
                        "msg": f"suspicous ; at the end of line ({line_count})",
                        "status": "Error",
                    }
                line_count += 1
        return {"name": "check_semi_column_presence", "msg": "No extra semi column found.", "status": "Success"}

    def _check_alpha_numerical(self):
        for column in ["sample_id", "sample_name", "sample_project"]:
            if column not in self.df.columns:
                continue
            for i, x in enumerate(self.df[column].values):
                status = str(x).replace("-", "").replace("_", "").isalnum()
                if status is False:
                    msg = f"type error: wrong {column} name in [Data] section (line {i+1}). Must be made of  alpha numerical characters, _, and - characters only. Found {x}"
                    return {"msg": msg, "name": "check_alpha_numerical", "status": "Error"}
        return {
            "name": "check_alpha_numerical",
            "msg": "sample names and ID looks correct in the Sample_ID, Sample_Name, and Project column (alpha numerical and - or _ characters)",
            "status": "Success",
        }

    def validate(self):
        """This method checks whether the sample sheet is correctly formatted

        Checks for:
            * presence of ; at the end of lines indicated an edition with excel that
              wrongly transformed the data into a pure CSV file
            * inconsistent numbers of columns in the [DATA] section, which must be
              CSV-like section
            * Extra lines at the end are ignored
            * special characters are forbidden except - and _
            * checks for Sample_ID column uniqueness
            * checks for index uniqueness (if single index)
            * checks for combo of dual indices uniqueness
            * checks that sample names are unique

        and raise a SystemExit error on the first found error.

        """
        # aggregates all checks
        checks = self.checker()

        # Stop after first error
        for check in checks:
            if check["status"] == "Error":
                sys.exit("\u274C " + str(check["msg"]))

    def _get_settings(self):
        data = {}
        for line in self.sections[self.settings_section]:
            # tolerate malformed / comma-less lines (e.g. Excel-mangled sheets
            # that use ';' separators or trailing ';;;'). We do not crash here so
            # that the meaningful error is reported by _check_semi_column_presence.
            if "," not in line:
                continue
            key, value = line.split(",", 1)
            data[key] = value
        return data

    settings = property(_get_settings)

    def _get_header(self):
        data = {}
        for line in self.sections["Header"]:
            key, value = line.split(",", 1)
            data[key] = value
        return data

    header = property(_get_header)

    def _get_instrument(self):
        try:
            return self.header["Instrument Type"]
        except KeyError:
            return None

    instrument = property(_get_instrument, doc="returns instrument name")

    def _get_adapter_kit(self):
        try:
            return self.header["Index Adapters"]
        except KeyError:
            return None

    index_adapters = property(_get_adapter_kit, doc="returns index adapters")

    def to_fasta(self, adapter_name=""):
        """Extract adapters from [Adapter] section and print them as a fasta file"""
        ar1 = self.settings["Adapter"]
        try:
            ar2 = self.settings["AdapterRead2"]
        except KeyError:
            ar2 = ""

        for name, index in zip(self.df["i7_index_id"], self.df["index"]):
            read = f"{ar1}{index}{ar2}"
            frmt = {"adapter": adapter_name, "name": name, "index": index}
            print(">{adapter}_index_{name}|name:{name}|seq:{index}".format(**frmt))
            print(read)

        if "index2" in self.df.columns:
            for name, index in zip(self.df["i5_index_id"], self.df["index2"]):
                read = f"{ar1}{index}{ar2}"
                frmt = {"adapter": adapter_name, "name": name, "index": index}
                print(">{adapter}_index_{name}|name:{name}|seq:{index}".format(**frmt))
                print(read)

    def quick_fix(self, output_filename):
        """Fix sample sheet

        Tyical error is when users save the samplesheet as CSV file in excel.
        This may add trailing ; characters at the end of section, which raises error
        in bcl2fastq.

        """
        found_data = False
        with open(self.filename) as fin:
            with open(output_filename, "w") as fout:
                for line in fin.readlines():
                    if line.startswith("[Data]"):
                        found_data = True

                    if found_data:
                        line = line.replace(";", ",")
                    else:
                        line = line.strip().rstrip(";")
                        line = line.replace(";", ",")
                        line = line.strip().rstrip(",")
                    fout.write(line.strip("\n") + "\n")


class BCLConvert(SampleSheet):
    """Reader and validator of Illumina v2 (BCL Convert) sample sheets.

    BCL Convert is the replacement for the bcl2fastq software. The sample sheet
    format (a.k.a. v2) is close to the bcl2fastq one (a.k.a. v1) handled by
    :class:`SampleSheet` but differs in a few places:

    - The tabular data lives in a ``[BCLConvert_Data]`` section (required) instead
      of ``[Data]``.
    - Settings live in a ``[BCLConvert_Settings]`` section (required). ``[Settings]``
      may still appear but is optional.
    - ``[Reads]`` uses cycle counts (``Read1Cycles``, ``Index1Cycles``, ...) rather
      than one integer per line.
    - The header carries ``FileFormatVersion,2``.
    - ``Sample_ID`` is compulsory (error if absent) and at least one sample is
      required. ``Sample_Name`` is ignored by BCL Convert (warning if present).
    - Adapter trimming/masking options moved into the settings section, e.g.
      ``Adapter`` became ``AdapterRead1`` and command-line options such as
      ``--minimum-trimmed-read-length`` became ``MinimumTrimmedReadLength``.
    - Barcode mismatches are set in the settings with ``BarcodeMismatchesIndex1``
      and ``BarcodeMismatchesIndex2``.

    Because sections are parsed generically by :meth:`SampleSheet._scan_sections`,
    most of the [Data] checks (mandatory columns, unique Sample_ID, unique/valid
    indices, homogeneous index lengths, ...) are reused as-is by pointing
    :attr:`data_section` and :attr:`settings_section` to the v2 sections.

    :references:
        - https://support.illumina.com/sequencing/sequencing_software/bcl-convert/compatibility.html
        - https://knowledge.illumina.com/software/general/software-general-reference_material-list/000003710
    """

    #: the tabular section is [BCLConvert_Data]; title-cased by the parser
    data_section = "Bclconvert_Data"

    #: the settings section is [BCLConvert_Settings]; title-cased by the parser
    settings_section = "Bclconvert_Settings"

    def checker(self):

        from sequana.utils.checker import Checker

        checks = Checker()

        checks.tryme(self._check_file_format_version)

        # [BCLConvert_Settings] is required in v2
        if self.settings_section in self.sections:
            checks.tryme(self._check_settings)
        else:
            checks.results.append({"msg": "The required [BCLConvert_Settings] section is missing", "status": "Error"})

        # [BCLConvert_Data] is required in v2
        if self.data_section in self.sections:
            checks.tryme(self._check_data_section_csv_format)
            checks.tryme(self._check_mandatory_data_columns)
            checks.tryme(self._check_sample_ID)
            checks.tryme(self._check_unique_sample_ID)
            checks.tryme(self._check_unique_indices)
            checks.tryme(self._check_nucleotide_indices)
            checks.tryme(self._check_sample_lane_number)
            checks.tryme(self._check_alpha_numerical)
            checks.tryme(self._check_sample_name_ignored)

            try:
                if "index" in self.df.columns:
                    checks.tryme(self._check_homogene_I7_length)
                if "index2" in self.df.columns:
                    checks.tryme(self._check_homogene_I5_length)
                    checks.tryme(self._check_homogene_I5_and_I7_length)
            except pd.errors.ParserError:  # pragma: no cover
                pass
        else:
            checks.results.append({"msg": "The required [BCLConvert_Data] section is missing", "status": "Error"})

        # [Cloud_Data] is optional BaseSpace metadata not used by BCL Convert demux,
        # but a malformed row (e.g. missing comma) is still worth flagging (Warning).
        if "Cloud_Data" in self.sections:
            checks.tryme(self._check_cloud_data_section_csv_format)

        checks.tryme(self._check_semi_column_presence)

        return checks.results

    def _check_file_format_version(self):
        version = self.header.get("FileFormatVersion") if "Header" in self.sections else None
        if version is None:
            return {
                "name": "check_file_format_version",
                "msg": "FileFormatVersion not found in the [Header] section. Expected FileFormatVersion,2 for BCL Convert.",
                "status": "Warning",
            }
        if str(version).strip() != "2":
            return {
                "name": "check_file_format_version",
                "msg": f"Unexpected FileFormatVersion ({version}). BCL Convert expects FileFormatVersion,2.",
                "status": "Error",
            }
        return {
            "name": "check_file_format_version",
            "msg": "FileFormatVersion is 2 (BCL Convert).",
            "status": "Success",
        }

    def _check_cloud_data_section_csv_format(self):
        return self._check_csv_format("Cloud_Data", "check_cloud_data_section_csv_format", level="Warning")

    def _check_sample_name_ignored(self):
        # Sample_Name is ignored by BCL Convert; warn if present to avoid confusion
        if "sample_name" in self.df.columns:
            return {
                "name": "check_sample_name_ignored",
                "msg": "Sample_Name column is present but is ignored by BCL Convert. Remove it to avoid confusion.",
                "status": "Warning",
            }
        return {"name": "check_sample_name_ignored", "msg": "No ignored Sample_Name column.", "status": "Success"}

    def quick_fix(self, output_filename):
        """Fix a v2 sample sheet by removing trailing semicolons only.

        Unlike :meth:`SampleSheet.quick_fix`, internal semicolons are preserved
        because they are legitimate v2 separators (e.g. ``OverrideCycles`` uses
        ``R1:Y151;I1:I10;I2:I10;R2:Y151``). Only trailing ``;`` (typical Excel
        artefact) are stripped.
        """
        with open(self.filename) as fin, open(output_filename, "w") as fout:
            for line in fin.readlines():
                fout.write(line.rstrip().rstrip(";").rstrip() + "\n")

    # obsolete v1 settings names (kept for reference):
    # Adapter, TrimAdapter, MaskAdapter, MaskAdapterRead2, Read1StartFromCycle, Read1EndWithCycle,
    # Read1UMIStartFromCycle, Read1UMILength, Read2UMIStartFromCycle, Read2UMILength, Read2StartFromCycle


def get_sample_sheet_version(filename):
    """Return ``"v2"`` for BCL Convert sample sheets, ``"v1"`` otherwise.

    Detection relies on the v2 markers: presence of a ``[BCLConvert_Data]`` or
    ``[BCLConvert_Settings]`` section, or ``FileFormatVersion,2`` in the header.
    """
    ss = SampleSheet(filename)
    if BCLConvert.data_section in ss.sections or BCLConvert.settings_section in ss.sections:
        return "v2"
    try:
        if str(ss.header.get("FileFormatVersion", "")).strip() == "2":
            return "v2"
    except Exception:  # pragma: no cover
        pass
    return "v1"


def SampleSheetFactory(filename):
    """Return a :class:`BCLConvert` or :class:`SampleSheet` instance based on the
    detected sample sheet version (see :func:`get_sample_sheet_version`)."""
    if get_sample_sheet_version(filename) == "v2":
        return BCLConvert(filename)
    return SampleSheet(filename)


class IEM(SampleSheet):
    def __init__(self, filename):
        super().__init__(filename)
        logger.warning("IEM class is deprecated. Use SampleShee instead.")
