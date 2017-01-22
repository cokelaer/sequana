# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2017 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import re

from sequana import logger

import ruamel.yaml

__all__ = ["YamlDocParser"]


class YamlDocParser(object):
    """A simple parser to extract block content in YAML files


    ::

        from sequana import snaketools
        from sequana.iotools import YamlDocParser
        module = snaketools.Module('quality_control')
        r = YamlDocParser(module.config)
        r.sections['fastqc']


    """
    def __init__(self, filename):
        """
        ::

            # main documentation

            # block comment
            section1:
                - item

            # blcok comment
            section2:

            # a comment

            section3:

        Here, section1 and section2 have block comments but not section3

        """
        self.filename = filename
        self.regex_section = re.compile("^[a-z,A-Z,_,0-9]+:")
        self.sections = {}
        self._read_data()
        self._parse_data()

    def _get_expected_sections(self):
        """Get the top level keys in the YAML file

        :return: list of top level sections' names"""
        with open(self.filename, "r") as fh:
            data = ruamel.yaml.load(fh.read(), ruamel.yaml.RoundTripLoader)
        keys = list(data.keys())
        return keys

    def _read_data(self):
        with open(self.filename, "r") as fh:
            self.data = fh.readlines()

    def _parse_data(self):
        """Parse the YAML file to get the block content (comments)
        before each top-level sections. See doc in the constructor

        Removes all # so that the block of comments can be interpreted as
        a docstring in Sequanix
        """

        current_block = []
        current_section = "docstring"

        # if we get a line that starts with #, this is a new comment or
        # part of a block comment. Otherwise, it means the current block
        # comment has ended.

        for this in self.data:
            # Beginning of a new section at top level
            if self.regex_section.findall(this):
                name = self.regex_section.findall(this)[0]
                current_section = name.strip(":")
                self.sections[current_section] = "".join(current_block)
                current_block = []
                current_section = None
            elif this.startswith('#'):    # a comment at top level
                current_block.append(this)
            elif this.strip() == "":      # an empty line
                #this was the main comment, or an isolated comment
                current_block = []
            else:  # a non-empty line to skip
                current_block = []

        for key in self._get_expected_sections():
            if key not in self.sections.keys():
                logger.warning("section %s not dealt by the parsing function" % key)

    def _block2docstring(self, section):
        if section not in self.sections.keys():
            logger.warning("%s not found in the yaml " % section)
            return
        comments = self.sections[section]
        docstring = []
        for line in comments.split("\n"):
            if "#############" in line:
                pass
            else:
                if len(line)<2: # an empty line (to keep)
                    docstring.append("")
                else:
                    docstring.append(line[2:]) # strip the "# "characters
        docstring = "\n".join(docstring)
        return docstring
