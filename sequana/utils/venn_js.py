# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Interactive Venn diagram (HTML/SVG) to be embedded in a Sequana report.

The Python side only exports the geometry of the Venn diagrams (see
:func:`sequana.viz.venn.get_layout`) together with a compact table of the gene
memberships. The intersections themselves are recomputed in the browser so that
the thresholds can be changed without regenerating the report.
"""

import json

import colorlog

from sequana.viz.venn import default_colors, get_layout

logger = colorlog.getLogger(__name__)


__all__ = ["DynamicVenn"]


def _to_rgba(color):
    r, g, b, a = color
    return f"rgba({int(255 * r)},{int(255 * g)},{int(255 * b)},{a})"


class DynamicVenn:
    """Build a self-contained interactive Venn diagram.

    :param names: the name of each set (e.g. the comparison names).
    :param rows: a list of rows. Each row has one integer code per set followed
        by the number of genes sharing that combination of codes. A code of 0
        means that the gene is not significant in that set. Otherwise, its sign
        gives the direction of the regulation and its magnitude minus one gives
        the index of the largest fold-change bin the gene belongs to (see
        :attr:`step`).
    :param step: the width of a fold-change bin.
    :param nbins: the number of fold-change bins, that is the number of
        positions of the fold-change slider.
    :param html_id: a unique identifier, in case several diagrams are included
        in a single report.

    The codes are meant to be built by :meth:`encode_log2_fold_change`.
    """

    def __init__(self, names, rows, step=0.25, nbins=17, html_id="dynamic_venn", max_sets=6):
        self.names = list(names)
        self.rows = rows
        self.step = step
        self.nbins = nbins
        self.html_id = html_id
        self.max_sets = max_sets

    @staticmethod
    def encode_log2_fold_change(l2fc, significant, step=0.25, nbins=17):
        """Encode a log2 fold change into a bin index (see :class:`DynamicVenn`).

        :param l2fc: a pandas Series of log2 fold changes.
        :param significant: a boolean Series telling whether the gene is significant.
        """
        import numpy as np

        magnitude = np.minimum((l2fc.abs() / step).fillna(0).astype(int) + 1, nbins)
        code = np.sign(l2fc).fillna(0).astype(int) * magnitude
        return (code * significant.astype(int)).astype(int)

    def _data(self):
        selected = list(range(min(len(self.names), self.max_sets)))
        return {
            "names": self.names,
            "rows": self.rows,
            "step": self.step,
            "nbins": self.nbins,
            "maxSets": self.max_sets,
            "selected": selected,
            "colors": [_to_rgba(c) for c in default_colors],
            "layouts": {n: get_layout(n) for n in range(2, self.max_sets + 1) if n <= len(self.names)},
        }

    def _controls(self):
        uid = self.html_id
        boxes = "\n".join(
            f'<label class="venn-set"><input type="checkbox" data-index="{i}"'
            f'{" checked" if i < self.max_sets else ""}> <span class="venn-swatch" id="{uid}_swatch_{i}"></span>'
            f"{name}</label>"
            for i, name in enumerate(self.names)
        )
        return f"""
<div class="venn-controls" id="{uid}_controls">
  <div>
    <label for="{uid}_slider">log2 fold change threshold:
      <b id="{uid}_threshold">1</b></label>
    <input type="range" id="{uid}_slider" min="0" max="{self.nbins - 1}" step="1" value="{int(round(1 / self.step))}"
           style="vertical-align:middle; width:280px">
  </div>
  <div>
    Genes:
    <label><input type="radio" name="{uid}_direction" value="all" checked> all</label>
    <label><input type="radio" name="{uid}_direction" value="up"> up-regulated</label>
    <label><input type="radio" name="{uid}_direction" value="down"> down-regulated</label>
    <label><input type="radio" name="{uid}_direction" value="split"> up and down separately</label>
  </div>
  <div id="{uid}_sets">Comparisons: {boxes}</div>
</div>
"""

    def to_html(self):
        """Return the HTML (markup, style and javascript) of the diagram."""
        uid = self.html_id
        data = json.dumps(self._data())

        style = """
<style>
.venn-controls { margin-bottom: 10px; }
.venn-controls > div { margin: 4px 0; }
.venn-controls label { font-weight: normal; margin-right: 12px; }
.venn-set { display: inline-block; white-space: nowrap; }
.venn-swatch { display:inline-block; width:12px; height:12px; margin-right:4px;
               border:1px solid #888; vertical-align:middle; }
.venn-message { color: #a00; }
</style>
"""

        script = """
<script>
(function () {
  var D = %(data)s;
  var uid = "%(uid)s";
  var svgNS = "http://www.w3.org/2000/svg";

  function el(name, attrs) {
    var node = document.createElementNS(svgNS, name);
    for (var key in attrs) { node.setAttribute(key, attrs[key]); }
    return node;
  }

  // axes coordinates (y upwards, 0-1) to svg coordinates (y downwards, 0-100)
  function X(x) { return 100 * x; }
  function Y(y) { return 100 * (1 - y); }

  var anchors = {left: "start", center: "middle", right: "end"};
  var baselines = {bottom: "text-after-edge", center: "central", top: "hanging"};

  function selectedSets() {
    var out = [];
    document.querySelectorAll("#" + uid + "_sets input").forEach(function (box) {
      if (box.checked) { out.push(parseInt(box.dataset.index, 10)); }
    });
    return out;
  }

  function direction() {
    return document.querySelector("input[name='" + uid + "_direction']:checked").value;
  }

  // the sets to be drawn. In the 'split' mode, every comparison gives two sets
  // (its up- and its down-regulated genes), which are disjoint by construction.
  function activeSets() {
    var way = direction();
    var out = [];
    selectedSets().forEach(function (i) {
      if (way === "split") {
        out.push({index: i, way: "up", label: D.names[i] + " (up)"});
        out.push({index: i, way: "down", label: D.names[i] + " (down)"});
      } else {
        out.push({index: i, way: way, label: D.names[i]});
      }
    });
    return out;
  }

  // count the genes of every region of the diagram. Regions are identified by
  // a binary key, '0110' being the genes shared by the 2nd and 3rd sets only.
  function counts(sets, threshold) {
    var out = {}, sizes = [];
    var i, j, row, code, keep, key, any;
    for (j = 0; j < sets.length; j++) { sizes.push(0); }
    for (i = 0; i < D.rows.length; i++) {
      row = D.rows[i];
      key = "";
      any = false;
      for (j = 0; j < sets.length; j++) {
        code = row[sets[j].index];
        keep = code !== 0 && Math.abs(code) - 1 >= threshold &&
               (sets[j].way === "all" || (sets[j].way === "up" && code > 0) ||
                (sets[j].way === "down" && code < 0));
        key += keep ? "1" : "0";
        if (keep) { any = true; sizes[j] += row[D.names.length]; }
      }
      if (any) { out[key] = (out[key] || 0) + row[D.names.length]; }
    }
    return {regions: out, sizes: sizes};
  }

  function draw(sets, result) {
    var layout = D.layouts[sets.length];
    // the shapes span 0-100 but the set names stick out on the sides
    var svg = el("svg", {viewBox: "-22 -6 144 112", width: "100%%", style: "max-width:800px"});

    layout.shapes.forEach(function (shape, i) {
      var kind = shape[0], p = shape[1], node;
      if (kind === "ellipse") {
        var cx = X(p.xy[0]), cy = Y(p.xy[1]);
        node = el("ellipse", {cx: cx, cy: cy, rx: 50 * p.width, ry: 50 * p.height,
                              fill: D.colors[i],
                              transform: "rotate(" + (-p.angle) + " " + cx + " " + cy + ")"});
      } else {
        node = el("polygon", {fill: D.colors[i],
                              points: p.xy.map(function (xy) { return X(xy[0]) + "," + Y(xy[1]); }).join(" ")});
      }
      svg.appendChild(node);
    });

    // the 5 and 6-set diagrams have small regions, hence a smaller font
    var fontsize = sets.length > 4 ? 2.2 : 3;

    layout.texts.forEach(function (t) {
      var node = el("text", {x: X(t.x), y: Y(t.y), "text-anchor": "middle",
                             "dominant-baseline": "central", "font-size": fontsize});
      node.textContent = result.regions[t.key] || 0;
      svg.appendChild(node);
    });

    layout.names.forEach(function (t, i) {
      var node = el("text", {x: X(t.x), y: Y(t.y), "text-anchor": anchors[t.ha],
                             "dominant-baseline": baselines[t.va], "font-size": fontsize + 0.2});
      node.textContent = sets[i].label + " (" + result.sizes[i] + ")";
      svg.appendChild(node);
    });

    return svg;
  }

  function update() {
    var target = document.getElementById(uid + "_plot");
    var sets = activeSets();
    var threshold = parseInt(document.getElementById(uid + "_slider").value, 10);
    document.getElementById(uid + "_threshold").textContent = (threshold * D.step).toFixed(2);
    target.innerHTML = "";
    if (sets.length < 2 || sets.length > D.maxSets) {
      var msg = document.createElement("p");
      msg.className = "venn-message";
      msg.textContent = "Please select the comparisons so that 2 to " + D.maxSets +
                        " sets are drawn (" + sets.length + " at the moment" +
                        (direction() === "split" ? ", every comparison counting as two sets" : "") + ").";
      target.appendChild(msg);
      return;
    }
    target.appendChild(draw(sets, counts(sets, threshold)));
  }

  D.names.forEach(function (name, i) {
    var swatch = document.getElementById(uid + "_swatch_" + i);
    if (swatch) { swatch.style.background = D.colors[i %% D.colors.length]; }
  });
  document.getElementById(uid + "_slider").addEventListener("input", update);
  document.querySelectorAll("#" + uid + "_controls input[type=checkbox], " +
                            "#" + uid + "_controls input[type=radio]").forEach(function (box) {
    box.addEventListener("change", update);
  });
  update();
})();
</script>
""" % {
            "data": data,
            "uid": uid,
        }

        return style + self._controls() + f'<div id="{uid}_plot"></div>' + script
