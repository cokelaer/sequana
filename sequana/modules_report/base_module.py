# coding: utf-8
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Dimitri Desvillechabrol <dimitri.desvillechabrol@pasteur.fr>,
#          <d.desvillechabrol@gmail.com>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""
DEPRECATED: This module is maintained for backward compatibility.

The base class has been moved to sequana_report package.
Please use: from sequana_report import SequanaBaseModule
"""

import warnings

# Re-export from sequana_report for backward compatibility
from sequana_report.base_module import SequanaBaseModule

warnings.warn(
    "Importing SequanaBaseModule from sequana.modules_report is deprecated. "
    "Please use: from sequana_report import SequanaBaseModule",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = ["SequanaBaseModule"]
