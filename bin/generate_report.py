#!/usr/bin/env python3

# Written by Fernando Duarte with AI assistance, released under the MIT license.
# Generates a self-contained HTML quality report for GenomeQC pipeline results
"""Generate a self-contained HTML quality report for GenomeQC pipeline results."""

import argparse
import re
import sys
from datetime import datetime
from pathlib import Path


def _safe_id(s):
    return re.sub(r'[^a-zA-Z0-9]', '_', s)

# ── Colour palette ────────────────────────────────────────────────────────────

BUSCO_COLORS = {
    "Single":      "#2196F3",
    "Duplicated":  "#4CAF50",
    "Fragmented":  "#FF9800",
    "Missing":     "#F44336",
}

TIDK_PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]

# ── Data parsers ──────────────────────────────────────────────────────────────

def parse_busco_batch_summaries(paths):
    """Return list of row-dicts from BUSCO batch_summary_modified.txt files."""
    rows = []
    for p in paths:
        with open(p) as fh:
            lines = fh.readlines()
        if len(lines) < 2:
            continue
        header = lines[0].strip().split("\t")
        for line in lines[1:]:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            row = dict(zip(header, parts + [""] * max(0, len(header) - len(parts))))
            rows.append(row)
    rows.sort(key=lambda r: r.get("Input_file", ""))
    return rows


def parse_decontam_tsv(paths):
    """Return {species: [row_dict, ...]} from decontamination report TSV/TXT files.

    Handles FCS-GX (##-prefixed metadata + #-prefixed header), FCS-Adaptor
    (#-prefixed header), and Tiara (plain header) formats.
    """
    result = {}
    for p in paths:
        species = Path(p).stem.split(".")[0]
        rows = []
        header = None
        with open(p) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line or line.startswith("##"):
                    continue
                if header is None:
                    header = line.lstrip("#").lstrip().split("\t")
                    continue
                parts = line.split("\t")
                rows.append(dict(zip(header, parts + [""] * max(0, len(header) - len(parts)))))
        if header is not None:
            result[species] = rows
    return result


def parse_tidk_tsvs(paths):
    """Return {species: [row_dict, ...]} from tidk aposteriori search TSV files."""
    result = {}
    for p in paths:
        species = Path(p).stem
        rows = []
        with open(p) as fh:
            lines = fh.readlines()
        if not lines:
            continue
        header = lines[0].strip().split("\t")
        for line in lines[1:]:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            row = dict(zip(header, parts))
            for k in ("window", "forward_repeat_number", "reverse_repeat_number"):
                if k in row:
                    try:
                        row[k] = int(row[k])
                    except ValueError:
                        pass
            rows.append(row)
        result[species] = rows
    return result


def parse_busco_seqs_table(path):
    """Return ({Id: count}, col_label) from ortho_seqs.py TSV output."""
    result = {}
    col_label = "Seqs above threshold"
    with open(path) as fh:
        lines = fh.readlines()
    if not lines:
        return result, col_label
    header = lines[0].strip().split('\t')
    count_col = next((h for h in header if h.startswith('Num_Seqs_Above_')), None)
    if count_col:
        threshold = count_col.replace('Num_Seqs_Above_', '')
        col_label = f"Seqs >{threshold} BUSCOs"
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        parts = line.split('\t')
        row = dict(zip(header, parts))
        sp_id = row.get('Id', '')
        if sp_id and count_col:
            try:
                result[sp_id] = int(row[count_col])
            except (ValueError, KeyError):
                result[sp_id] = '—'
    return result, col_label


# ── AGAT sp_statistics parsing ───────────────────────────────────────────────

# AGAT's own section names aren't consistently spelled (e.g. "lncrna" vs
# "lnc_rna" depending on the source annotation), so sections are grouped under
# a normalised key (lowercased, underscores stripped) and given a display label
# here. Unrecognised sections fall back to a title-cased version of the key.
_AGAT_SECTION_LABELS = {
    "region": "Region", "sequencefeature": "Sequence feature",
    "genefeature": "Gene feature", "guiderna": "Guide RNA",
    "lncrna": "lncRNA", "mirna": "miRNA", "mrna": "mRNA",
    "primarytranscript": "Primary transcript", "rna": "Other RNA (pseudogene)",
    "rrna": "rRNA", "snorna": "snoRNA", "snrna": "snRNA",
    "transcript": "Transcript", "trna": "tRNA",
}


def _agat_section_key(name):
    return re.sub(r"[^a-z0-9]", "", name.lower())


def _agat_section_label(key):
    return _AGAT_SECTION_LABELS.get(key, key.replace("_", " ").capitalize())


_AGAT_SECTION_HEADER_RE = re.compile(r"^-{2,}\s+([A-Za-z][\w]*)\s+-{2,}$")
_AGAT_STAT_LINE_RE = re.compile(r"^(.*\S)\s{2,}(\S+)\s*$")


def parse_agat_stats(path):
    """Parse an AGAT sp_statistics.pl text report.

    Returns {section_key: {"label": str, "stats": {field: value},
    "stats_no_isoforms": {field: value} | None}}. Sections and fields present
    depend entirely on which feature types exist in the annotation - callers
    must not assume any particular section (e.g. "mrna") is present.
    """
    sections = {}
    current = None
    target = None  # points at either "stats" or "stats_no_isoforms" dict
    with open(path) as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")
            stripped = line.strip()
            if not stripped:
                continue

            header = _AGAT_SECTION_HEADER_RE.match(stripped)
            if header:
                key = _agat_section_key(header.group(1))
                current = sections.setdefault(
                    key, {"label": _agat_section_label(key), "stats": {}, "stats_no_isoforms": None}
                )
                target = current["stats"]
                continue

            if current is not None and "have isoforms!" in stripped:
                current["stats_no_isoforms"] = {}
                target = current["stats_no_isoforms"]
                continue

            m = _AGAT_STAT_LINE_RE.match(stripped)
            if m and target is not None:
                target[m.group(1)] = m.group(2)

    return sections


def parse_agat_stats_files(paths):
    """Return {species: {section_key: {...}}} from AGAT sp_statistics.pl files."""
    result = {}
    for p in paths:
        species = Path(p).stem.replace(".stats", "")
        parsed = parse_agat_stats(p)
        if parsed:
            result[species] = parsed
    return result


def agat_gene_count(species_stats):
    """Best-effort protein-coding-like gene count for the Overview summary table.

    Prefers the "mrna" section; some annotations (e.g. AUGUSTUS predictions,
    or AGAT-merged multi-evidence files) label protein-coding transcripts
    "transcript" instead - see ORTHOLOGOUS_CHROMOSOMES's own handling of this
    same ambiguity. Falls back through both before giving up.
    """
    for key in ("mrna", "transcript", "primarytranscript"):
        section = species_stats.get(key)
        if section:
            count = section["stats"].get("Number of gene")
            if count is not None:
                return count
    return None


# ── RepeatMasker .tbl parsing ──────────────────────────────────────────────────

_RM_HEADER_RES = {
    "sequences":    re.compile(r'^sequences:\s*([\d,]+)'),
    "total_length": re.compile(r'^total length:\s*([\d,]+)\s*bp\s*\(([\d,]+)\s*bp excl N/X-runs\)'),
    "gc_level":     re.compile(r'^GC level:\s*([\d.]+)\s*%'),
    "bases_masked": re.compile(r'^bases (?:masked|covered):\s*([\d,]+)\s*bp\s*\(\s*([\d.]+)\s*%\)'),
}

_RM_ROW_RE = re.compile(
    r'^(?P<indent>\s*)(?P<label>[^\d].*?):?\s+'
    r'(?:(?P<elements>\d[\d,]*)\s+)?'
    r'(?P<length>\d[\d,]*)\s*bp\s*\(?\s*(?P<pct>[\d.]+)\s*%\)?\s*$'
)


def parse_repeatmasker_tbl(path):
    """Parse a RepeatMasker (or te_tbl.py) .tbl file.

    Returns (info_dict, rows) where rows is a list of
    {label, indent, elements, length_bp, pct}.
    """
    info = {}
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            stripped = line.strip()
            if not stripped or stripped.startswith(("=", "-", "*", "file name")):
                continue
            if stripped.startswith("number of") or stripped.startswith("elements"):
                continue
            matched = False
            for key, rex in _RM_HEADER_RES.items():
                m = rex.match(stripped)
                if not m:
                    continue
                matched = True
                if key == "sequences":
                    info["sequences"] = m.group(1).replace(",", "")
                elif key == "total_length":
                    info["total_length_bp"] = m.group(1).replace(",", "")
                    info["total_length_excl_bp"] = m.group(2).replace(",", "")
                elif key == "gc_level":
                    info["gc_level_pct"] = m.group(1)
                elif key == "bases_masked":
                    info["bases_masked_bp"] = m.group(1).replace(",", "")
                    info["bases_masked_pct"] = m.group(2)
                break
            if matched:
                continue
            m = _RM_ROW_RE.match(line)
            if m:
                elements = m.group("elements")
                rows.append({
                    "label":     m.group("label").strip(),
                    "indent":    len(m.group("indent")) // 2,
                    "elements":  elements.replace(",", "") if elements else "",
                    "length_bp": m.group("length").replace(",", ""),
                    "pct":       m.group("pct"),
                })
    return info, rows


def parse_repeatmasker_tbls(paths):
    """Return {species: (info_dict, rows)} from RepeatMasker .tbl files."""
    result = {}
    for p in paths:
        species = Path(p).stem
        result[species] = parse_repeatmasker_tbl(p)
    return result


def parse_quast_tsvs(paths):
    """Return ({species: {metric: value}}, metric_order) from QUAST report TSVs.

    Each file is transposed (one metric per line: name, then value), so we
    read them directly rather than treating the first line as a header.
    """
    species_data = {}
    metric_order = []
    for p in sorted(paths, key=lambda x: Path(x).name):
        species = Path(p).stem.replace(".quast", "")
        data = {}
        with open(p) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                metric, value = parts[0], parts[1]
                data[metric] = value
                if metric not in metric_order:
                    metric_order.append(metric)
        species_data[species] = data
    return species_data, metric_order


# ── SVG chart generators ──────────────────────────────────────────────────────

def busco_stacked_bar_svg(rows, mode_label=""):
    if not rows:
        return "<p><em>No BUSCO data available.</em></p>"

    bar_h     = 30
    label_w   = 220
    chart_w   = 550
    pad       = 12
    legend_h  = 36
    row_gap   = 6
    header_h  = 55
    total_h   = header_h + len(rows) * (bar_h + row_gap) + legend_h + pad
    svg_w     = label_w + chart_w + pad * 2

    bars = []
    y = header_h

    for row in rows:
        species = row.get("Input_file", "Unknown")
        display = (species[:35] + "…") if len(species) > 36 else species
        try:
            single = float(row.get("Single", 0))
            dup    = float(row.get("Duplicated", 0))
            frag   = float(row.get("Fragmented", 0))
            miss   = float(row.get("Missing", 0))
        except ValueError:
            continue

        bars.append(
            f'<text x="{label_w - 6}" y="{y + bar_h // 2 + 4}" '
            f'text-anchor="end" font-size="11" font-family="sans-serif" fill="#333">'
            f'<title>{species}</title>{display}</text>'
        )

        x = label_w + pad
        for val, color, lbl in [
            (single, BUSCO_COLORS["Single"],     f"Complete Single: {single:.1f}%"),
            (dup,    BUSCO_COLORS["Duplicated"],  f"Complete Dup: {dup:.1f}%"),
            (frag,   BUSCO_COLORS["Fragmented"],  f"Fragmented: {frag:.1f}%"),
            (miss,   BUSCO_COLORS["Missing"],     f"Missing: {miss:.1f}%"),
        ]:
            w = val / 100 * chart_w
            if w > 0.5:
                bars.append(
                    f'<rect x="{x:.1f}" y="{y}" width="{w:.1f}" height="{bar_h}" '
                    f'fill="{color}" rx="2"><title>{lbl}</title></rect>'
                )
            x += w

        complete = single + dup
        bars.append(
            f'<text x="{label_w + pad + complete / 100 * chart_w + 5}" '
            f'y="{y + bar_h // 2 + 4}" font-size="10" fill="#555">{complete:.1f}%</text>'
        )
        y += bar_h + row_gap

    # Axis grid lines and tick labels
    axis = []
    for pct in range(0, 101, 20):
        xpos = label_w + pad + pct / 100 * chart_w
        axis.append(
            f'<text x="{xpos}" y="46" text-anchor="middle" font-size="10" fill="#888">{pct}%</text>'
            f'<line x1="{xpos}" y1="50" x2="{xpos}" y2="{total_h - legend_h - pad}" '
            f'stroke="#ddd" stroke-width="1" stroke-dasharray="3,3"/>'
        )

    # Legend
    ley = total_h - legend_h + 4
    legend = []
    lx = label_w + pad
    for key, color in BUSCO_COLORS.items():
        legend.append(
            f'<rect x="{lx}" y="{ley}" width="14" height="14" fill="{color}" rx="2"/>'
            f'<text x="{lx + 18}" y="{ley + 11}" font-size="11" fill="#555">{key}</text>'
        )
        lx += 110

    title_lbl = "BUSCO Completeness" + (f" — {mode_label}" if mode_label else "")
    title = (
        f'<text x="{svg_w // 2}" y="22" text-anchor="middle" '
        f'font-size="14" font-weight="600" font-family="sans-serif" fill="#222">'
        f'{title_lbl}</text>'
    )

    return (
        f'<div style="overflow-x:auto">'
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_w}" height="{total_h}" '
        f'style="max-width:100%;height:auto;display:block">'
        f'{title}{"".join(axis)}{"".join(bars)}{"".join(legend)}'
        f'</svg></div>'
    )


def tidk_line_svg(species, rows, width=700, height=220, plot_id=None):
    """Return an HTML block (optional chromosome dropdown + SVG waveform).

    Mirrors the original tidk SVG plot: one path per chromosome, y encodes total
    repeat count (fwd+rev) with high density plotted toward the top.
    A chromosome selector dropdown is shown when the species has >1 sequence.
    """
    if not rows:
        return ""

    chroms = {}
    for r in rows:
        cid = r.get("id", "chr")
        chroms.setdefault(cid, []).append(r)

    # plot_id lets the caller differentiate IDs when both modes share a page
    sp_id = _safe_id(plot_id if plot_id else species)

    # Use fwd+rev sum as the signal value
    max_sum = max(
        (r["forward_repeat_number"] + r["reverse_repeat_number"]
         for rs in chroms.values() for r in rs),
        default=1,
    ) or 1
    max_window = max(
        (r["window"] for rs in chroms.values() for r in rs),
        default=1,
    ) or 1

    pl, pr, pt, pb = 50, 15, 30, 50
    pw = width - pl - pr
    ph = height - pt - pb
    baseline_y = pt + ph  # SVG y of zero-count baseline (bottom of plot)

    def tx(w):
        return pl + w / max_window * pw

    def ty(s):
        return baseline_y - s / max_sum * ph

    paths_svg = []
    for i, (cid, chrom_rows) in enumerate(list(chroms.items())[:10]):
        color = TIDK_PALETTE[i % len(TIDK_PALETTE)]
        cid_safe = _safe_id(cid)
        s0 = chrom_rows[0]["forward_repeat_number"] + chrom_rows[0]["reverse_repeat_number"]
        pts = [f"M{pl:.1f},{ty(s0):.1f}"]
        for r in chrom_rows:
            s = r["forward_repeat_number"] + r["reverse_repeat_number"]
            pts.append(f"L{tx(r['window']):.1f},{ty(s):.1f}")
        paths_svg.append(
            f'<path id="tp-{sp_id}-{cid_safe}" d="{" ".join(pts)}" fill="none" '
            f'stroke="{color}" stroke-width="1.5" opacity="0.85">'
            f'<title>{cid}</title></path>'
        )

    # Baseline
    baseline = (
        f'<line x1="{pl}" y1="{baseline_y}" x2="{pl + pw}" y2="{baseline_y}" '
        f'stroke="#ccc" stroke-width="1"/>'
    )

    # Horizontal grid lines
    y_grid = []
    for pct in (25, 50, 75, 100):
        ypos = ty(max_sum * pct / 100)
        y_grid.append(
            f'<line x1="{pl}" y1="{ypos:.1f}" x2="{pl + pw}" y2="{ypos:.1f}" '
            f'stroke="#eee" stroke-width="1" stroke-dasharray="3,3"/>'
            f'<text x="{pl - 4}" y="{ypos + 4:.1f}" text-anchor="end" '
            f'font-size="9" fill="#aaa">{int(max_sum * pct / 100)}</text>'
        )

    # X-axis ticks
    x_grid = []
    for pct in range(0, 101, 25):
        xpos = tx(max_window * pct / 100)
        mbp = max_window * pct / 100 / 1_000_000
        label = f"{mbp:.2f}Mb" if mbp >= 0.1 else f"{int(max_window * pct / 100 / 1000)}k"
        x_grid.append(
            f'<line x1="{xpos:.1f}" y1="{baseline_y}" x2="{xpos:.1f}" y2="{baseline_y + 4}" '
            f'stroke="#bbb" stroke-width="1"/>'
            f'<text x="{xpos:.1f}" y="{baseline_y + 14}" text-anchor="middle" '
            f'font-size="9" fill="#888">{label}</text>'
        )

    # Chromosome legend (below x-axis)
    leg = []
    n_chroms = min(len(chroms), 10)
    leg_item_w = min(80, pw // max(n_chroms, 1))
    for i, cid in enumerate(list(chroms.keys())[:10]):
        color = TIDK_PALETTE[i % len(TIDK_PALETTE)]
        lx = pl + i * leg_item_w
        if lx + leg_item_w > pl + pw:
            break
        leg.append(
            f'<line x1="{lx}" y1="{baseline_y + 28}" x2="{lx + 12}" y2="{baseline_y + 28}" '
            f'stroke="{color}" stroke-width="2"/>'
            f'<text x="{lx + 15}" y="{baseline_y + 32}" font-size="9" fill="#555">{cid[:10]}</text>'
        )

    border = (
        f'<rect x="{pl}" y="{pt}" width="{pw}" height="{ph}" '
        f'fill="none" stroke="#ccc" stroke-width="1"/>'
    )
    title_svg = (
        f'<text x="{width // 2}" y="18" text-anchor="middle" '
        f'font-size="12" font-weight="600" font-family="sans-serif" fill="#333">'
        f'{species}</text>'
    )
    y_axis_label = (
        f'<text x="{pl - 38}" y="{pt + ph // 2}" text-anchor="middle" '
        f'font-size="9" fill="#888" '
        f'transform="rotate(-90 {pl - 38} {pt + ph // 2})">Repeat density</text>'
    )

    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'style="max-width:100%;height:auto;display:block">'
        f'{title_svg}{border}{baseline}'
        f'{"".join(y_grid)}{"".join(x_grid)}'
        f'{y_axis_label}'
        f'{"".join(paths_svg)}'
        f'{"".join(leg)}'
        f'</svg>'
    )

    # Chromosome selector dropdown — only rendered when there are multiple sequences
    if len(chroms) > 1:
        chrom_items = []
        for i, cid in enumerate(list(chroms.keys())[:10]):
            cid_safe = _safe_id(cid)
            color = TIDK_PALETTE[i % len(TIDK_PALETTE)]
            chrom_items.append(
                f'<label>'
                f'<input type="checkbox" value="{cid_safe}" checked '
                f'onchange="tidkUpdate(\'{sp_id}\')">'
                f'<span class="chrom-swatch" style="background:{color}"></span>'
                f'{cid}'
                f'</label>'
            )
        dropdown = (
            f'<div class="chrom-select-wrap">'
            f'<button class="chrom-btn" onclick="tidkToggleMenu(\'{sp_id}\',event)">'
            f'Sequences ▾</button>'
            f'<div class="chrom-menu" id="cmenu-{sp_id}">'
            f'<div class="chrom-menu-actions">'
            f'<button onclick="tidkSelectAll(\'{sp_id}\')">All</button>'
            f'<button onclick="tidkSelectNone(\'{sp_id}\')">None</button>'
            f'</div>'
            f'{"".join(chrom_items)}'
            f'</div>'
            f'</div>'
        )
    else:
        dropdown = ""

    return f'{dropdown}{svg}'


# ── HTML helpers ──────────────────────────────────────────────────────────────

def _th(cells):
    return "<tr>" + "".join(f"<th>{c}</th>" for c in cells) + "</tr>"


def _td(cells):
    return "<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>"


def busco_table_html(rows):
    cols = ["Input_file", "Dataset", "Complete", "Single", "Duplicated",
            "Fragmented", "Missing", "n_markers", "Scaffold N50", "Number of scaffolds"]
    available = [c for c in cols if any(c in r for r in rows)]
    header = _th(available)
    body = "\n".join(_td([row.get(c, "") for c in available]) for row in rows)
    return f'<table class="table">{header}{body}</table>'


def busco_panel_html(tab_id, rows, mode_label):
    """Build a BUSCO tab panel (chart + table) for the given mode."""
    chart = busco_stacked_bar_svg(rows, mode_label=mode_label)
    table = busco_table_html(rows)
    return (
        f'<div id="{tab_id}" class="tab-panel">'
        f'<div class="card"><h2>Completeness chart</h2>{chart}</div>'
        f'<div class="card"><h2>Completeness table</h2>'
        f'<p style="margin-bottom:10px;font-size:12px;color:#888">Hover chart bars for tooltips. '
        f'C(S) = Complete single-copy &nbsp;·&nbsp; C(D) = Complete duplicated &nbsp;·&nbsp; '
        f'F = Fragmented &nbsp;·&nbsp; M = Missing</p>'
        f'{table}</div></div>'
    )


def quast_table_html(species_data, metric_order):
    """Wide QUAST table: species as rows, metrics as columns."""
    if not species_data:
        return "<p><em>No QUAST statistics available.</em></p>"
    species_list = sorted(species_data.keys())
    header = _th(["Assembly"] + metric_order)
    body = "\n".join(
        _td([sp] + [species_data[sp].get(m, "NA") for m in metric_order])
        for sp in species_list
    )
    return f'<div style="overflow-x:auto"><table class="table">{header}{body}</table></div>'


def decontam_table_html(data, no_hit_msg="No contamination detected."):
    """Render per-species tables for FCS-GX or FCS-Adaptor data."""
    if not data:
        return "<p><em>No data available.</em></p>"
    parts = []
    for species in sorted(data.keys()):
        rows = data[species]
        if not rows:
            parts.append(
                f'<h3 style="margin:16px 0 4px">{species} '
                f'<span class="badge badge-green">Clean</span></h3>'
                f'<p style="color:#888;font-size:12px;margin-bottom:16px">{no_hit_msg}</p>'
            )
            continue
        cols = list(rows[0].keys())
        body = "\n".join(_td([r.get(c, "") for c in cols]) for r in rows)
        n = len(rows)
        badge = f'<span class="badge badge-orange">{n} hit{"s" if n != 1 else ""}</span>'
        parts.append(
            f'<h3 style="margin:16px 0 8px">{species} {badge}</h3>'
            f'<div style="overflow-x:auto;margin-bottom:8px">'
            f'<table class="table">{_th(cols)}{body}</table></div>'
        )
    return "".join(parts)


def tiara_summary_html(data):
    """Render a per-species classification summary table for Tiara data."""
    if not data:
        return "<p><em>No Tiara data available.</em></p>"
    classes = ["eukarya", "bacteria", "archaea", "prokarya", "mitochondria", "plastid", "unknown"]
    header = _th(["Assembly", "Status"] + [c.capitalize() for c in classes] + ["Total"])
    rows_html = []
    for species in sorted(data.keys()):
        rows = data[species]
        counts = {c: 0 for c in classes}
        for r in rows:
            cls = r.get("class_fst_stage", "unknown").lower()
            counts[cls] = counts.get(cls, 0) + 1
        total = len(rows)
        non_euk = total - counts.get("eukarya", 0)
        status = (
            f'<span class="badge badge-orange">{non_euk} non-eukaryote</span>'
            if non_euk > 0 else
            f'<span class="badge badge-green">All eukaryote</span>'
        )
        cells = [species, status] + [str(counts.get(c, 0)) for c in classes] + [str(total)]
        rows_html.append(_td(cells))
    body = "\n".join(rows_html)
    return f'<div style="overflow-x:auto"><table class="table">{header}{body}</table></div>'


def agat_feature_keys(agat_data):
    """Union of feature-type keys across all species, as {key: label}."""
    features = {}
    for sections in agat_data.values():
        for key, section in sections.items():
            features.setdefault(key, section["label"])
    return features


def agat_feature_has_isoforms(agat_data, key):
    return any(
        sections.get(key, {}).get("stats_no_isoforms") is not None
        for sections in agat_data.values()
    )


# AGAT reports dozens of metrics per feature type, most of them rarely useful
# at a glance. Default tables show a compact subset - AGAT's own first few
# counts (which are always the basic "how many X" fields) plus a couple of
# specific metrics worth always keeping - with a toggle to see every column.
_AGAT_RELEVANT_FIRST_N = 4
_AGAT_RELEVANT_EXTRA_METRICS = [
    "mean gene length (bp)", "Number of single exon gene",
    "Number gene overlapping", "Total gene length (bp)",
]


def agat_feature_metrics(agat_data, key, view):
    """Ordered union of metric names for one feature type/view: metrics that
    more species report come first, ties broken by order of first appearance.

    AGAT's field set for a section isn't fully fixed across species - e.g. it
    inserts an extra "Number of pseudogene" count (shifting every later field
    over by one) when a GFF tags pseudogenes, and some sections' labels echo
    the source GFF's own feature-type spelling (e.g. "lnc_rna" vs "lncrna").
    Ordering by first-appearance alone would make the default columns depend
    on whichever species happens to sort first, rather than on what most
    species actually share.
    """
    stat_key = "stats" if view == "all" else "stats_no_isoforms"
    order = []
    seen = set()
    counts = {}
    for sp in sorted(agat_data.keys()):
        stats = (agat_data[sp].get(key) or {}).get(stat_key) or {}
        for m in stats:
            counts[m] = counts.get(m, 0) + 1
            if m not in seen:
                seen.add(m)
                order.append(m)
    rank = {m: i for i, m in enumerate(order)}
    return sorted(order, key=lambda m: (-counts[m], rank[m]))


def agat_relevant_metrics(metrics):
    relevant = list(metrics[:_AGAT_RELEVANT_FIRST_N])
    for m in _AGAT_RELEVANT_EXTRA_METRICS:
        if m in metrics and m not in relevant:
            relevant.append(m)
    return relevant


def agat_feature_table_html(agat_data, key, view, metrics):
    """Wide table for one feature type: species as rows, given metrics as
    columns. A species missing this feature, or missing a particular metric
    within it, shows "NA" rather than being left out of the table.
    """
    stat_key = "stats" if view == "all" else "stats_no_isoforms"
    species_list = sorted(agat_data.keys())

    header = _th(["Assembly"] + metrics)
    body = "\n".join(
        _td([sp] + [((agat_data[sp].get(key) or {}).get(stat_key) or {}).get(m, "NA") for m in metrics])
        for sp in species_list
    )
    return f'<div style="overflow-x:auto"><table class="table">{header}{body}</table></div>'


def _agat_view_block(agat_data, key, view, view_id, hidden=False):
    """Build one isoform-view's markup: a compact table plus a hidden full-column
    table, both scoped by a per-key CSS class so a single "show all columns"
    checkbox can toggle them together regardless of which isoform view is active.
    """
    metrics = agat_feature_metrics(agat_data, key, view)
    relevant = agat_relevant_metrics(metrics)
    has_extra = len(relevant) < len(metrics)

    relevant_table = agat_feature_table_html(agat_data, key, view, relevant)
    full_table = (
        f'<div class="agat-cols-full-{key}" style="display:none">'
        f'{agat_feature_table_html(agat_data, key, view, metrics)}</div>'
        if has_extra else ""
    )
    relevant_div = (
        f'<div class="agat-cols-relevant-{key}">{relevant_table}</div>'
        if has_extra else relevant_table
    )
    style = ' style="display:none"' if hidden else ""
    return f'<div id="{view_id}"{style}>{relevant_div}{full_table}</div>', has_extra


def agat_panel_html(agat_data):
    """Build the AGAT stats tab body: a feature-type dropdown showing one wide
    table (species x metric) at a time, so species are easy to compare for a
    given feature type instead of clicking through species one at a time."""
    if not agat_data:
        return "<p><em>No AGAT statistics available.</em></p>"

    features = agat_feature_keys(agat_data)
    feature_keys = sorted(features.keys(), key=lambda k: features[k])
    feature_options = "\n".join(
        f'<option value="{key}">{features[key]}</option>' for key in feature_keys
    )

    feature_panels = []
    for i, key in enumerate(feature_keys):
        display = "" if i == 0 else ' style="display:none"'
        has_iso = agat_feature_has_isoforms(agat_data, key)

        all_block, has_extra_all = _agat_view_block(agat_data, key, "all", f"agat-table-all-{key}")
        if has_iso:
            collapsed_block, has_extra_collapsed = _agat_view_block(
                agat_data, key, "collapsed", f"agat-table-collapsed-{key}", hidden=True
            )
        else:
            collapsed_block, has_extra_collapsed = "", False

        iso_toggle = (
            f'<label style="font-size:12px;color:#555;margin-bottom:8px;display:inline-block;margin-right:16px">'
            f'<input type="checkbox" onchange="agatToggleIsoforms(\'{key}\',this.checked)"> '
            f'Collapse isoforms (one transcript per gene)</label>'
            if has_iso else ""
        )
        cols_toggle = (
            f'<label style="font-size:12px;color:#555;margin-bottom:8px;display:inline-block">'
            f'<input type="checkbox" onchange="agatToggleColumns(\'{key}\',this.checked)"> '
            f'Show all columns</label>'
            if (has_extra_all or has_extra_collapsed) else ""
        )

        feature_panels.append(
            f'<div id="agat-feature-{key}" class="agat-feature-panel"{display}>'
            f'{iso_toggle}{cols_toggle}'
            f'{all_block}'
            f'{collapsed_block}'
            f'</div>'
        )

    return (
        f'<div style="margin-bottom:14px">'
        f'<label for="agat-feature-select" style="font-size:13px;color:#555;margin-right:6px">'
        f'Feature type:</label>'
        f'<select id="agat-feature-select" onchange="agatShowFeature(this.value)">{feature_options}</select>'
        f'</div>'
        f'{"".join(feature_panels)}'
    )


def _fmt_int(v):
    try:
        return f"{int(v):,}"
    except (ValueError, TypeError):
        return v or "—"


def repeatmasker_table_html(info, rows):
    """Render a RepeatMasker .tbl as a summary table + repeat-class breakdown table."""
    header_rows = [
        ("Sequences",                          _fmt_int(info.get("sequences"))),
        ("Total length (bp)",                  _fmt_int(info.get("total_length_bp"))),
        ("Total length excl. N/X-runs (bp)",   _fmt_int(info.get("total_length_excl_bp"))),
        ("GC level",                           f'{info.get("gc_level_pct", "—")} %'),
        ("Bases masked (bp)",                  _fmt_int(info.get("bases_masked_bp"))),
        ("Bases masked",                       f'{info.get("bases_masked_pct", "—")} %'),
    ]
    info_html = "\n".join(
        f'<tr><td><strong>{k}</strong></td><td>{v}</td></tr>' for k, v in header_rows
    )

    body_rows = []
    for r in rows:
        label = (
            f'<span style="padding-left:{r["indent"] * 16}px;display:inline-block">'
            f'{r["label"]}</span>'
        )
        elements = _fmt_int(r["elements"]) if r["elements"] else "—"
        body_rows.append(_td([label, elements, f'{_fmt_int(r["length_bp"])} bp', f'{r["pct"]} %']))

    table = (
        f'<table class="table">'
        f'{_th(["Category", "Number of elements", "Length occupied", "% of sequence"])}'
        f'{"".join(body_rows)}'
        f'</table>'
    )

    return (
        f'<table class="table" style="max-width:420px;margin-bottom:16px">{info_html}</table>'
        f'<div style="overflow-x:auto">{table}</div>'
    )


def _busco_complete_pct(row):
    """Format a BUSCO row's Complete value as a percentage, or em-dash if absent."""
    if not row:
        return "—"
    try:
        return f"{float(row.get('Complete', 0)):.1f}%"
    except ValueError:
        return row.get("Complete", "—")


def summary_table_html(busco_rows, tidk_data, busco_seqs_data=None, busco_seqs_col=None,
                       busco_prot_rows=None, agat_data=None):
    """Cross-tool summary table shown on the Overview tab."""
    has_prot = bool(busco_prot_rows)
    # Collect unique species from all data sources
    species_set = {r.get("Input_file", "") for r in busco_rows} | set(tidk_data.keys())
    if has_prot:
        species_set |= {r.get("Input_file", "") for r in busco_prot_rows}
    if agat_data is not None:
        species_set |= set(agat_data.keys())
    species_list = sorted(s for s in species_set if s)

    busco_by_species      = {r.get("Input_file", ""): r for r in busco_rows}
    busco_prot_by_species = {r.get("Input_file", ""): r for r in (busco_prot_rows or [])}

    # When protein BUSCO is present, disambiguate the genome column and add a
    # protein completeness column right after it.
    genome_col = "BUSCO genome complete (%)" if has_prot else "BUSCO complete (%)"
    cols = ["Assembly", genome_col]
    if has_prot:
        cols.append("BUSCO proteins complete (%)")
    cols += ["BUSCO lineage", "Scaffold N50", "# scaffolds", "Telomeric repeat"]
    if busco_seqs_data is not None:
        cols.append(busco_seqs_col or "Seqs above threshold")
    if agat_data is not None:
        cols.append("Genes (AGAT)")
    header = _th(cols)
    rows_html = []
    for sp in species_list:
        br = busco_by_species.get(sp, {})
        n50 = br.get("Scaffold N50", "—") or "—"
        n_scaffolds = br.get("Number of scaffolds", "—") or "—"
        lineage = br.get("Dataset", "—") or "—"
        tidk_rows = tidk_data.get(sp, [])
        repeat = tidk_rows[0].get("telomeric_repeat", "—") if tidk_rows else "—"
        cells = [sp, _busco_complete_pct(br)]
        if has_prot:
            cells.append(_busco_complete_pct(busco_prot_by_species.get(sp, {})))
        cells += [lineage, n50, n_scaffolds, repeat]
        if busco_seqs_data is not None:
            cells.append(str(busco_seqs_data.get(sp, "—")))
        if agat_data is not None:
            cells.append(agat_gene_count(agat_data.get(sp, {})) or "—")
        rows_html.append(_td(cells))
    body = "\n".join(rows_html)
    return f'<table class="table">{header}{body}</table>'


# ── Inline CSS + JS ───────────────────────────────────────────────────────────

CSS = """
*{box-sizing:border-box;margin:0;padding:0}
body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;font-size:14px;color:#222;background:#f5f6fa}
header{background:#1565C0;color:#fff;padding:14px 24px;display:flex;align-items:center;gap:12px}
header h1{font-size:20px;font-weight:600}
header span{font-size:12px;opacity:.7}
.container{max-width:1200px;margin:24px auto;padding:0 16px}
nav.tabs{display:flex;gap:0;border-bottom:2px solid #e0e0e0;margin-bottom:20px;flex-wrap:wrap}
nav.tabs button{background:none;border:none;padding:10px 20px;cursor:pointer;font-size:14px;color:#555;border-bottom:3px solid transparent;margin-bottom:-2px;transition:color .15s,border-color .15s}
nav.tabs button:hover{color:#1565C0}
nav.tabs button.active{color:#1565C0;border-bottom-color:#1565C0;font-weight:600}
.tab-panel{display:none}
.tab-panel.active{display:block}
.card{background:#fff;border-radius:8px;box-shadow:0 1px 4px rgba(0,0,0,.1);padding:20px 24px;margin-bottom:20px}
.card h2{font-size:16px;font-weight:600;margin-bottom:14px;color:#1565C0}
.card h3{font-size:14px;font-weight:600;margin:16px 0 8px;color:#333}
table.table{width:100%;border-collapse:collapse;font-size:13px}
table.table th{background:#e3f2fd;color:#1565C0;padding:8px 10px;text-align:left;font-weight:600;white-space:normal;vertical-align:bottom;overflow-wrap:break-word}
table.table td{padding:7px 10px;border-bottom:1px solid #f0f0f0}
table.table tr:hover td{background:#fafafa}
.badge{display:inline-block;padding:2px 8px;border-radius:99px;font-size:11px;font-weight:600}
.badge-blue{background:#e3f2fd;color:#1565C0}
.badge-green{background:#e8f5e9;color:#2e7d32}
.badge-orange{background:#fff3e0;color:#e65100}
.tag{font-size:11px;color:#888}
.tidk-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(660px,1fr));gap:16px}
.tidk-item{background:#fff;border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,.08);padding:12px}
.tidk-item h3{font-size:13px;font-weight:600;margin-bottom:8px;color:#333}
.chrom-select-wrap{position:relative;display:inline-block;margin-bottom:8px}
.chrom-btn{background:#f0f4ff;border:1px solid #c5cae9;border-radius:4px;padding:4px 10px;cursor:pointer;font-size:12px;color:#3949ab;line-height:1.4}
.chrom-btn:hover{background:#e8eaf6}
.chrom-menu{position:absolute;top:calc(100% + 4px);left:0;background:#fff;border:1px solid #ddd;border-radius:6px;box-shadow:0 4px 12px rgba(0,0,0,.12);padding:8px;z-index:100;min-width:180px;max-height:260px;overflow-y:auto;display:none}
.chrom-menu.open{display:block}
.chrom-menu label{display:flex;align-items:center;gap:6px;padding:3px 4px;font-size:12px;cursor:pointer;white-space:nowrap;border-radius:3px}
.chrom-menu label:hover{background:#f5f5f5}
.chrom-menu input[type=checkbox]{cursor:pointer}
.chrom-swatch{display:inline-block;width:10px;height:10px;border-radius:2px;flex-shrink:0}
.chrom-menu-actions{display:flex;gap:6px;margin-bottom:6px;padding-bottom:6px;border-bottom:1px solid #eee}
.chrom-menu-actions button{flex:1;background:#f0f4ff;border:1px solid #c5cae9;border-radius:3px;padding:2px 6px;font-size:11px;cursor:pointer;color:#3949ab}
.chrom-menu-actions button:hover{background:#e8eaf6}
.tidk-mode-toggle{display:flex;gap:4px;margin-bottom:8px;flex-wrap:wrap}
.tidk-mode-toggle button{background:#f5f5f5;border:1px solid #ddd;border-radius:4px;padding:4px 12px;cursor:pointer;font-size:12px;color:#555;transition:background .1s,border-color .1s}
.tidk-mode-toggle button.active{background:#e3f2fd;border-color:#90caf9;color:#1565C0;font-weight:600}
.tidk-mode-toggle button:hover:not(.active){background:#ebebeb}
.tidk-mode-toggle code{font-size:11px;background:rgba(0,0,0,.06);padding:1px 4px;border-radius:3px}
footer{text-align:center;padding:24px;font-size:12px;color:#aaa}
"""

TAB_JS = """
document.querySelectorAll('nav.tabs button').forEach(btn => {
  btn.addEventListener('click', () => {
    document.querySelectorAll('nav.tabs button').forEach(b => b.classList.remove('active'));
    document.querySelectorAll('.tab-panel').forEach(p => p.classList.remove('active'));
    btn.classList.add('active');
    document.getElementById(btn.dataset.tab).classList.add('active');
  });
});

function tidkToggleMenu(spId, event) {
  event.stopPropagation();
  var menu = document.getElementById('cmenu-' + spId);
  menu.classList.toggle('open');
}

function tidkUpdate(spId) {
  var menu = document.getElementById('cmenu-' + spId);
  menu.querySelectorAll('input[type=checkbox]').forEach(function(cb) {
    var path = document.getElementById('tp-' + spId + '-' + cb.value);
    if (path) path.style.display = cb.checked ? '' : 'none';
  });
}

function tidkSelectAll(spId) {
  var menu = document.getElementById('cmenu-' + spId);
  menu.querySelectorAll('input[type=checkbox]').forEach(function(cb) { cb.checked = true; });
  tidkUpdate(spId);
}

function tidkSelectNone(spId) {
  var menu = document.getElementById('cmenu-' + spId);
  menu.querySelectorAll('input[type=checkbox]').forEach(function(cb) { cb.checked = false; });
  tidkUpdate(spId);
}

document.addEventListener('click', function(e) {
  if (!e.target.closest('.chrom-select-wrap')) {
    document.querySelectorAll('.chrom-menu.open').forEach(function(m) { m.classList.remove('open'); });
  }
});

function rmShowSpecies(spId) {
  document.querySelectorAll('.rm-panel').forEach(function(p) {
    p.style.display = (p.id === 'rm-panel-' + spId) ? '' : 'none';
  });
}

function agatShowFeature(key) {
  document.querySelectorAll('.agat-feature-panel').forEach(function(p) {
    p.style.display = (p.id === 'agat-feature-' + key) ? '' : 'none';
  });
}

function agatToggleIsoforms(key, checked) {
  var allView = document.getElementById('agat-table-all-' + key);
  var collapsedView = document.getElementById('agat-table-collapsed-' + key);
  if (allView) allView.style.display = checked ? 'none' : '';
  if (collapsedView) collapsedView.style.display = checked ? '' : 'none';
}

function agatToggleColumns(key, checked) {
  document.querySelectorAll('.agat-cols-relevant-' + key).forEach(function(el) {
    el.style.display = checked ? 'none' : '';
  });
  document.querySelectorAll('.agat-cols-full-' + key).forEach(function(el) {
    el.style.display = checked ? '' : 'none';
  });
}

function tidkSetMode(spId, mode) {
  ['aposteriori', 'apriori'].forEach(function(m) {
    var el = document.getElementById('plot-' + m + '-' + spId);
    if (el) el.style.display = (m === mode) ? '' : 'none';
  });
  var toggle = document.getElementById('mode-' + spId);
  if (toggle) toggle.querySelectorAll('button').forEach(function(btn) {
    btn.classList.toggle('active', btn.dataset.mode === mode);
  });
}
"""

# ── Main HTML assembly ────────────────────────────────────────────────────────

def build_html(busco_rows, tidk_data, tidk_apriori_data=None,
               fcsgx_data=None, fcsadp_data=None, tiara_data=None,
               busco_seqs_data=None, busco_seqs_col=None,
               busco_prot_rows=None, repeatmasker_data=None, agat_data=None,
               quast_data=None):
    tabs = []
    panels = []
    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    # ── Overview tab ──────────────────────────────────────────────────────────
    tabs.append(("overview", "Overview"))
    n_species = len({r.get("Input_file") for r in busco_rows}
                    | {r.get("Input_file") for r in (busco_prot_rows or [])}
                    | set(tidk_data.keys()))
    lineages  = sorted({r.get("Dataset", "") for r in busco_rows if r.get("Dataset")})
    badges = " ".join(
        f'<span class="badge badge-blue">{lg}</span>' for lg in lineages
    )
    summary_tbl = (
        summary_table_html(busco_rows, tidk_data, busco_seqs_data, busco_seqs_col,
                           busco_prot_rows, agat_data)
        if (busco_rows or busco_prot_rows or tidk_data or agat_data) else "<p>No data available.</p>"
    )

    panels.append(
        f'<div id="overview" class="tab-panel active">'
        f'<div class="card">'
        f'<h2>Run summary</h2>'
        f'<p><strong>{n_species}</strong> assemblies analysed &nbsp;·&nbsp; '
        f'BUSCO lineage(s): {badges or "—"} &nbsp;·&nbsp; '
        f'<span class="tag">Generated {now}</span></p>'
        f'</div>'
        f'<div class="card"><h2>Per-assembly overview</h2>'
        f'<div style="overflow-x:auto">{summary_tbl}</div></div>'
        f'</div>'
    )

    # ── BUSCO tab(s) ──────────────────────────────────────────────────────────
    # When protein BUSCO is present (genome + annotation mode), show genome and
    # protein completeness in separate tabs; otherwise a single "BUSCO" tab.
    if busco_rows:
        genome_label = "BUSCO (Genome)" if busco_prot_rows else "BUSCO"
        tabs.append(("busco", genome_label))
        panels.append(busco_panel_html("busco", busco_rows, "Genome" if busco_prot_rows else ""))

    if busco_prot_rows:
        tabs.append(("busco_prot", "BUSCO (Proteins)"))
        panels.append(busco_panel_html("busco_prot", busco_prot_rows, "Proteins"))

    # ── QUAST tab ─────────────────────────────────────────────────────────────
    if quast_data and quast_data[0]:
        quast_species, quast_metrics = quast_data
        tabs.append(("quast", "Assembly stats"))
        panels.append(
            f'<div id="quast" class="tab-panel">'
            f'<div class="card">'
            f'<h2>QUAST — Assembly statistics</h2>'
            f'{quast_table_html(quast_species, quast_metrics)}'
            f'</div></div>'
        )

    # ── Telomeres tab ─────────────────────────────────────────────────────────
    all_tidk_species = sorted(set(tidk_data.keys()) | set((tidk_apriori_data or {}).keys()))
    if all_tidk_species:
        tabs.append(("tidk", "Telomeres"))
        items = []
        for sp in all_tidk_species:
            post_rows = tidk_data.get(sp, [])
            pre_rows  = (tidk_apriori_data or {}).get(sp, [])
            sp_id     = _safe_id(sp)
            has_both  = bool(post_rows and pre_rows)

            repeat = (post_rows or pre_rows)[0].get("telomeric_repeat", "—")
            repeat_badge = (
                f'<span class="badge badge-green">Repeat: {repeat}</span>'
                if repeat != "—" else ""
            )

            if has_both:
                # Repeat badges for each mode
                post_repeat = post_rows[0].get("telomeric_repeat", "—")
                pre_repeat  = pre_rows[0].get("telomeric_repeat", "—")
                toggle = (
                    f'<div class="tidk-mode-toggle" id="mode-{sp_id}">'
                    f'<button class="active" data-mode="aposteriori" '
                    f'onclick="tidkSetMode(\'{sp_id}\',\'aposteriori\')">'
                    f'A posteriori'
                    f'{f" &middot; <code>{post_repeat}</code>" if post_repeat != "—" else ""}'
                    f'</button>'
                    f'<button data-mode="apriori" '
                    f'onclick="tidkSetMode(\'{sp_id}\',\'apriori\')">'
                    f'A priori'
                    f'{f" &middot; <code>{pre_repeat}</code>" if pre_repeat != "—" else ""}'
                    f'</button>'
                    f'</div>'
                )
                plot_content = (
                    f'{toggle}'
                    f'<div id="plot-aposteriori-{sp_id}">'
                    f'{tidk_line_svg(sp, post_rows, plot_id=sp + "__post")}</div>'
                    f'<div id="plot-apriori-{sp_id}" style="display:none">'
                    f'{tidk_line_svg(sp, pre_rows, plot_id=sp + "__pre")}</div>'
                )
            elif post_rows:
                plot_content = tidk_line_svg(sp, post_rows)
            else:
                plot_content = tidk_line_svg(sp, pre_rows)

            items.append(
                f'<div class="tidk-item">'
                f'<h3>{sp} {repeat_badge}</h3>'
                f'{plot_content}</div>'
            )

        panels.append(
            f'<div id="tidk" class="tab-panel">'
            f'<div class="card">'
            f'<h2>Telomeric repeat analysis</h2>'
            f'<p style="margin-bottom:14px;font-size:12px;color:#888">'
            f'Total telomere repeat density (forward + reverse) per 10 kb window.</p>'
            f'<div class="tidk-grid">{"".join(items)}</div>'
            f'</div></div>'
        )

    # ── Repeats tab ───────────────────────────────────────────────────────────
    if repeatmasker_data:
        tabs.append(("repeats", "Repeats"))
        species_list = sorted(repeatmasker_data.keys())
        options = "\n".join(
            f'<option value="{_safe_id(sp)}">{sp}</option>' for sp in species_list
        )
        rm_panels = []
        for i, sp in enumerate(species_list):
            info, rows = repeatmasker_data[sp]
            sp_id = _safe_id(sp)
            display = "" if i == 0 else ' style="display:none"'
            rm_panels.append(
                f'<div id="rm-panel-{sp_id}" class="rm-panel"{display}>'
                f'{repeatmasker_table_html(info, rows)}</div>'
            )
        panels.append(
            f'<div id="repeats" class="tab-panel">'
            f'<div class="card">'
            f'<h2>Repeat content (RepeatMasker)</h2>'
            f'<div style="margin-bottom:14px">'
            f'<label for="rm-species-select" style="font-size:13px;color:#555;margin-right:6px">'
            f'Assembly:</label>'
            f'<select id="rm-species-select" onchange="rmShowSpecies(this.value)">{options}</select>'
            f'</div>'
            f'{"".join(rm_panels)}'
            f'</div></div>'
        )

    # ── AGAT annotation stats tab ─────────────────────────────────────────────
    if agat_data:
        tabs.append(("agat", "Annotation stats"))
        panels.append(
            f'<div id="agat" class="tab-panel">'
            f'<div class="card">'
            f'<h2>AGAT annotation statistics</h2>'
            f'<p style="margin-bottom:14px;font-size:12px;color:#888">'
            f'Per-assembly gene/transcript statistics, broken down by feature type. '
            f'Sections shown depend on which feature types are present in each annotation.</p>'
            f'{agat_panel_html(agat_data)}'
            f'</div></div>'
        )

    # ── FCS-GX tab ────────────────────────────────────────────────────────────
    if fcsgx_data is not None:
        tabs.append(("fcsgx", "FCS-GX"))
        content = decontam_table_html(fcsgx_data, "No foreign sequences detected.")
        panels.append(
            f'<div id="fcsgx" class="tab-panel">'
            f'<div class="card"><h2>FCS-GX — Foreign sequence detection</h2>'
            f'<p style="margin-bottom:14px;font-size:12px;color:#888">'
            f'Sequences flagged for removal due to cross-species contamination.</p>'
            f'{content}</div></div>'
        )

    # ── FCS-Adaptor tab ───────────────────────────────────────────────────────
    if fcsadp_data is not None:
        tabs.append(("fcsadp", "FCS-Adaptor"))
        content = decontam_table_html(fcsadp_data, "No adaptor contamination detected.")
        panels.append(
            f'<div id="fcsadp" class="tab-panel">'
            f'<div class="card"><h2>FCS-Adaptor — Adaptor contamination</h2>'
            f'<p style="margin-bottom:14px;font-size:12px;color:#888">'
            f'Sequences flagged for adaptor trimming or removal.</p>'
            f'{content}</div></div>'
        )

    # ── Tiara tab ─────────────────────────────────────────────────────────────
    if tiara_data is not None:
        tabs.append(("tiara", "Tiara"))
        content = tiara_summary_html(tiara_data)
        panels.append(
            f'<div id="tiara" class="tab-panel">'
            f'<div class="card"><h2>Tiara — Sequence classification</h2>'
            f'<p style="margin-bottom:14px;font-size:12px;color:#888">'
            f'Deep-learning classification of sequences by taxonomic origin.</p>'
            f'{content}</div></div>'
        )

    # ── Assemble page ─────────────────────────────────────────────────────────
    tab_nav = "\n".join(
        f'<button data-tab="{tid}" class="{"active" if i == 0 else ""}">{label}</button>'
        for i, (tid, label) in enumerate(tabs)
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>GenomeQC Report</title>
<style>{CSS}</style>
</head>
<body>
<header>
  <svg width="28" height="28" viewBox="0 0 28 28" fill="none" xmlns="http://www.w3.org/2000/svg">
    <circle cx="14" cy="14" r="13" stroke="white" stroke-width="2"/>
    <path d="M7 14 Q10 7 14 14 Q18 21 21 14" stroke="white" stroke-width="2" fill="none"/>
    <circle cx="14" cy="14" r="2.5" fill="white"/>
  </svg>
  <div>
    <h1>GenomeQC Report</h1>
    <span>nf-core/genomeqc &nbsp;·&nbsp; {now}</span>
  </div>
</header>
<div class="container">
  <nav class="tabs">{tab_nav}</nav>
  {"".join(panels)}
</div>
<footer>Generated by <strong>nf-core/genomeqc</strong> — <a href="https://github.com/nf-core/genomeqc" style="color:#888">github.com/nf-core/genomeqc</a></footer>
<script>{TAB_JS}</script>
</body>
</html>
"""


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Generate a self-contained HTML report for GenomeQC results."
    )
    parser.add_argument(
        "--busco_tables", nargs="*", default=[],
        metavar="TSV",
        help="BUSCO (genome) batch_summary_modified.txt files (one per species or combined)",
    )
    parser.add_argument(
        "--busco_prot_tables", nargs="*", default=[],
        metavar="TSV",
        help="BUSCO (proteins) batch_summary_modified.txt files (one per species or combined)",
    )
    parser.add_argument(
        "--tidk_tsvs", nargs="*", default=[],
        metavar="TSV",
        help="tidk aposteriori search TSV files (one per species)",
    )
    parser.add_argument(
        "--tidk_apriori_tsvs", nargs="*", default=[],
        metavar="TSV",
        help="tidk apriori search TSV files (one per species)",
    )
    parser.add_argument(
        "--fcsgx_reports", nargs="*", default=None,
        metavar="TXT",
        help="FCS-GX *.fcs_gx_report.txt files (one per species, omit to hide tab)",
    )
    parser.add_argument(
        "--fcsadp_reports", nargs="*", default=None,
        metavar="TXT",
        help="FCS-Adaptor *.fcs_adaptor_report.txt files (one per species, omit to hide tab)",
    )
    parser.add_argument(
        "--tiara_reports", nargs="*", default=None,
        metavar="TXT",
        help="Tiara *.txt classification files (one per species, omit to hide tab)",
    )
    parser.add_argument(
        "--busco_seqs_table", default=None,
        metavar="TSV",
        help="ortho_seqs.py output TSV (sequences above BUSCO threshold)",
    )
    parser.add_argument(
        "--repeatmasker_tbls", nargs="*", default=[],
        metavar="TBL",
        help="RepeatMasker *.tbl files (one per species, omit to hide tab)",
    )
    parser.add_argument(
        "--agat_stats", nargs="*", default=[],
        metavar="TXT",
        help="AGAT sp_statistics.pl *.stats.txt files (one per species, omit to hide tab)",
    )
    parser.add_argument(
        "--quast_tsvs", nargs="*", default=[],
        metavar="TSV",
        help="QUAST report.tsv files (one per species, omit to hide tab)",
    )
    parser.add_argument(
        "--output", default="genomeqc_report.html",
        metavar="HTML",
        help="Output HTML file path (default: genomeqc_report.html)",
    )
    args = parser.parse_args()

    # Parse BUSCO (genome and protein)
    busco_rows      = parse_busco_batch_summaries(args.busco_tables)      if args.busco_tables      else []
    busco_prot_rows = parse_busco_batch_summaries(args.busco_prot_tables) if args.busco_prot_tables else []

    # Parse tidk TSVs
    tidk_data         = parse_tidk_tsvs(args.tidk_tsvs)         if args.tidk_tsvs         else {}
    tidk_apriori_data = parse_tidk_tsvs(args.tidk_apriori_tsvs) if args.tidk_apriori_tsvs else None

    # Parse decontamination reports (None = tool not run → tab hidden)
    fcsgx_data  = parse_decontam_tsv(args.fcsgx_reports)  if args.fcsgx_reports  is not None else None
    fcsadp_data = parse_decontam_tsv(args.fcsadp_reports) if args.fcsadp_reports is not None else None
    tiara_data  = parse_decontam_tsv(args.tiara_reports)  if args.tiara_reports  is not None else None

    busco_seqs_data, busco_seqs_col = parse_busco_seqs_table(args.busco_seqs_table) if args.busco_seqs_table else (None, None)

    repeatmasker_data = parse_repeatmasker_tbls(args.repeatmasker_tbls) if args.repeatmasker_tbls else None

    agat_data = parse_agat_stats_files(args.agat_stats) if args.agat_stats else None

    quast_data = parse_quast_tsvs(args.quast_tsvs) if args.quast_tsvs else None

    if not busco_rows and not busco_prot_rows and not tidk_data and not tidk_apriori_data and not agat_data and not quast_data:
        print("WARNING: no input data found; generating empty report.", file=sys.stderr)

    html = build_html(busco_rows, tidk_data, tidk_apriori_data,
                      fcsgx_data=fcsgx_data, fcsadp_data=fcsadp_data, tiara_data=tiara_data,
                      busco_seqs_data=busco_seqs_data, busco_seqs_col=busco_seqs_col,
                      busco_prot_rows=busco_prot_rows, repeatmasker_data=repeatmasker_data,
                      agat_data=agat_data, quast_data=quast_data)

    Path(args.output).write_text(html)
    print(f"Report written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
