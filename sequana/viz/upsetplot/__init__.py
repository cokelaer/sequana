"""UpSet plots, vendored from UpSetPlot 0.9.0.

UpSetPlot (https://github.com/jnothman/UpSetPlot) is distributed under the New
BSD License, Copyright (c) 2018-2024 Joel Nothman. The LICENSE file sits next to
this module.

The code is vendored (rather than used as a dependency) because the last release
(0.9.0) is not compatible with pandas >= 3: the in-place ``fillna`` calls used to
set the default colours are silently ignored by the copy-on-write semantics,
leaving NaN colours that matplotlib rejects. Local changes are marked with a
``sequana:`` comment.
"""

from .data import from_contents, from_indicators, from_memberships, generate_counts, generate_data, generate_samples
from .plotting import UpSet, plot
from .reformat import query

__all__ = [
    "UpSet",
    "from_contents",
    "from_indicators",
    "from_memberships",
    "generate_counts",
    "generate_data",
    "generate_samples",
    "plot",
    "query",
]
