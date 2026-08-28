#!/usr/bin/env python3

# Written by Fernando Duarte with AI assistance, released under the MIT license.
# Generates an Excel workbook with raw data tables from GenomeQC pipeline results
"""Generate an Excel workbook with raw data tables from GenomeQC pipeline results.

Each tool's output is placed in a separate sheet. Multiple per-species files of
the same type are concatenated into one sheet (with the species name prepended
as an extra column where the format allows).

No external dependencies — XLSX is written using only the Python standard library.
"""

import argparse
import io
import re
import sys
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape as _xe


# ── Minimal XLSX writer ───────────────────────────────────────────────────────

class _XLSX:
    """Write a simple xlsx workbook with no external dependencies."""

    def __init__(self):
        self._sheets = []   # [(safe_name, rows)]

    def add_sheet(self, name, rows):
        safe = name[:31].translate(str.maketrans(r'/\[]*?:', '-------'))
        self._sheets.append((safe, list(rows)))

    def save(self, path):
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("[Content_Types].xml", self._content_types())
            zf.writestr("_rels/.rels",          self._root_rels())
            zf.writestr("xl/workbook.xml",       self._workbook())
            zf.writestr("xl/_rels/workbook.xml.rels", self._wb_rels())
            zf.writestr("xl/styles.xml",         self._styles())
            for i, (_, rows) in enumerate(self._sheets, 1):
                zf.writestr(f"xl/worksheets/sheet{i}.xml", self._sheet(rows))
        Path(path).write_bytes(buf.getvalue())

    # ── XML fragments ─────────────────────────────────────────────────────────

    def _content_types(self):
        overrides = "\n".join(
            f'<Override PartName="/xl/worksheets/sheet{i}.xml" '
            f'ContentType="application/vnd.openxmlformats-officedocument'
            f'.spreadsheetml.worksheet+xml"/>'
            for i in range(1, len(self._sheets) + 1)
        )
        return (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
            '<Default Extension="xml"  ContentType="application/xml"/>'
            '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
            '<Override PartName="/xl/styles.xml"   ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
            f'{overrides}</Types>'
        )

    def _root_rels(self):
        return (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" '
            'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
            'Target="xl/workbook.xml"/>'
            '</Relationships>'
        )

    def _workbook(self):
        sheets = "".join(
            f'<sheet name="{_xe(name)}" sheetId="{i}" r:id="rId{i}"/>'
            for i, (name, _) in enumerate(self._sheets, 1)
        )
        return (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
            'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
            f'<sheets>{sheets}</sheets></workbook>'
        )

    def _wb_rels(self):
        rels = "".join(
            f'<Relationship Id="rId{i}" '
            f'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
            f'Target="worksheets/sheet{i}.xml"/>'
            for i in range(1, len(self._sheets) + 1)
        )
        return (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            f'{rels}</Relationships>'
        )

    def _styles(self):
        return (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
            '<fonts><font><sz val="11"/><name val="Calibri"/></font></fonts>'
            '<fills>'
            '<fill><patternFill patternType="none"/></fill>'
            '<fill><patternFill patternType="gray125"/></fill>'
            '</fills>'
            '<borders><border><left/><right/><top/><bottom/><diagonal/></border></borders>'
            '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
            '<cellXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/></cellXfs>'
            '</styleSheet>'
        )

    @staticmethod
    def _col_letter(idx):
        """0-based column index → Excel letter (A, B, …, Z, AA, …)."""
        s = ""
        idx += 1
        while idx:
            idx, r = divmod(idx - 1, 26)
            s = chr(65 + r) + s
        return s

    def _sheet(self, rows):
        rows_xml = []
        for r_i, row in enumerate(rows, 1):
            cells = []
            for c_i, val in enumerate(row):
                ref = f"{self._col_letter(c_i)}{r_i}"
                txt = "" if val is None else str(val)
                cells.append(
                    f'<c r="{ref}" t="inlineStr"><is><t>{_xe(txt)}</t></is></c>'
                )
            rows_xml.append(f'<row r="{r_i}">{"".join(cells)}</row>')
        return (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
            f'<sheetData>{"".join(rows_xml)}</sheetData>'
            '</worksheet>'
        )


# ── Summary-table parsers (mirror generate_report.py) ────────────────────────

def _parse_busco_batch_summaries(paths):
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
            rows.append(dict(zip(header, parts + [""] * max(0, len(header) - len(parts)))))
    rows.sort(key=lambda r: r.get("Input_file", ""))
    return rows


def _parse_tidk_tsvs(paths):
    """Return {species: [row_dict, ...]} from tidk aposteriori search TSV files."""
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


def _parse_busco_seqs_table(path):
    """Return ({species: count}, col_label) from ortho_seqs.py TSV output."""
    result = {}
    col_label = "Seqs above threshold"
    with open(path) as fh:
        lines = fh.readlines()
    if not lines:
        return result, col_label
    header = lines[0].strip().split("\t")
    count_col = next((h for h in header if h.startswith("Num_Seqs_Above_")), None)
    if count_col:
        threshold = count_col.replace("Num_Seqs_Above_", "")
        col_label = f"Seqs >{threshold} BUSCOs"
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        row = dict(zip(header, parts))
        sp_id = row.get("Id", "")
        if sp_id and count_col:
            try:
                result[sp_id] = int(row[count_col])
            except (ValueError, KeyError):
                result[sp_id] = "—"
    return result, col_label


# ── RepeatMasker .tbl parsing (mirrors generate_report.py) ───────────────────

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


def _parse_repeatmasker_tbl(path):
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


def _repeatmasker_rows(paths):
    """Combine RepeatMasker .tbl files into a single table (one sheet, all species)."""
    rows = [["assembly", "category", "number_of_elements", "length_occupied_bp", "percent_of_sequence"]]
    for p in sorted(paths, key=lambda x: Path(x).name):
        species = Path(p).stem
        info, cat_rows = _parse_repeatmasker_tbl(p)
        header_metrics = [
            ("Sequences",                        info.get("sequences", ""),         ""),
            ("Total length (bp)",                info.get("total_length_bp", ""),   ""),
            ("Total length excl. N/X-runs (bp)", info.get("total_length_excl_bp", ""), ""),
            ("GC level (%)",                     "",                                 info.get("gc_level_pct", "")),
            ("Bases masked (bp)",                info.get("bases_masked_bp", ""),   info.get("bases_masked_pct", "")),
        ]
        for label, length, pct in header_metrics:
            rows.append([species, label, "", length, pct])
        for r in cat_rows:
            label = ("  " * r["indent"]) + r["label"]
            rows.append([species, label, r["elements"], r["length_bp"], r["pct"]])
    return rows


def _busco_complete_pct(row):
    """Format a BUSCO row's Complete value as a percentage, or em-dash if absent."""
    if not row:
        return "—"
    try:
        return f"{float(row.get('Complete', 0)):.1f}%"
    except (ValueError, TypeError):
        return row.get("Complete", "—") or "—"


def _summary_rows(busco_rows, tidk_data, busco_seqs_data=None, busco_seqs_col=None,
                  busco_prot_rows=None):
    """Build the summary table rows (header + data) mirroring the HTML overview."""
    has_prot = bool(busco_prot_rows)
    species_set = {r.get("Input_file", "") for r in busco_rows} | set(tidk_data.keys())
    if has_prot:
        species_set |= {r.get("Input_file", "") for r in busco_prot_rows}
    species_list = sorted(s for s in species_set if s)
    busco_by_species      = {r.get("Input_file", ""): r for r in busco_rows}
    busco_prot_by_species = {r.get("Input_file", ""): r for r in (busco_prot_rows or [])}

    genome_col = "BUSCO genome complete (%)" if has_prot else "BUSCO complete (%)"
    cols = ["Assembly", genome_col]
    if has_prot:
        cols.append("BUSCO proteins complete (%)")
    cols += ["BUSCO lineage", "Scaffold N50", "# scaffolds", "Telomeric repeat"]
    if busco_seqs_data is not None:
        cols.append(busco_seqs_col or "Seqs above threshold")
    rows = [cols]
    for sp in species_list:
        br = busco_by_species.get(sp, {})
        n50        = br.get("Scaffold N50", "—") or "—"
        n_scaffolds = br.get("Number of scaffolds", "—") or "—"
        lineage    = br.get("Dataset", "—") or "—"
        tidk_sp    = tidk_data.get(sp, [])
        repeat     = tidk_sp[0].get("telomeric_repeat", "—") if tidk_sp else "—"
        cells = [sp, _busco_complete_pct(br)]
        if has_prot:
            cells.append(_busco_complete_pct(busco_prot_by_species.get(sp, {})))
        cells += [lineage, n50, n_scaffolds, repeat]
        if busco_seqs_data is not None:
            cells.append(str(busco_seqs_data.get(sp, "—")))
        rows.append(cells)
    return rows


# ── TSV/TXT readers ───────────────────────────────────────────────────────────

def _read_tsv(path, skip_comment=True):
    """Yield rows (list of str) from a TSV, stripping leading # from header."""
    with open(path) as fh:
        first = True
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if skip_comment and line.startswith("##"):
                continue  # metadata lines (e.g. FCS-GX JSON header) — always skip
            if first and skip_comment:
                line = line.lstrip("#").lstrip()
                first = False
            elif skip_comment and line.startswith("#"):
                continue
            yield line.split("\t")


def _combine_tsv_files(paths, add_species_col=True):
    """Concatenate multiple TSV files into a single list of rows.

    When add_species_col is True, a 'species' column is prepended and the
    header is written once (from the first file).
    """
    rows = []
    header_written = False
    for p in sorted(paths, key=lambda x: Path(x).name):
        species = Path(p).stem.split(".")[0]  # strip extra suffixes
        file_rows = list(_read_tsv(p))
        if not file_rows:
            continue
        if not header_written:
            header = (["assembly"] + file_rows[0]) if add_species_col else file_rows[0]
            rows.append(header)
            header_written = True
            data_rows = file_rows[1:]
        else:
            data_rows = file_rows[1:]   # skip repeated header
        if data_rows:
            for r in data_rows:
                rows.append(([species] + r) if add_species_col else r)
        elif add_species_col:
            rows.append([species, "no contamination detected"] + [""] * (len(file_rows[0]) - 1))
    return rows


# AGAT's report has no "key: value" syntax - fields are a label, then 2+
# spaces, then the value (e.g. "Number of gene<spaces>9934"), grouped under
# "--- sectionname ---" headers. The same label (e.g. "Number of gene")
# repeats under every section with a different value each time, so the
# section name has to be carried along as its own column, not folded into
# the metric label.
_AGAT_SECTION_HEADER_RE = re.compile(r"^-{2,}\s+([A-Za-z][\w]*)\s+-{2,}$")
_AGAT_STAT_LINE_RE = re.compile(r"^(.*\S)\s{2,}(\S+)\s*$")


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


def _parse_agat_stats_nested(paths):
    """Parse AGAT stats files into {section_key: {species: {"all": {...}, "collapsed": {...}|None}}}."""
    sections = {}
    for p in sorted(paths, key=lambda x: Path(x).name):
        species = Path(p).stem.replace(".stats", "")
        entry = None
        target = None
        with open(p) as fh:
            for raw_line in fh:
                stripped = raw_line.strip()
                if not stripped:
                    continue
                header = _AGAT_SECTION_HEADER_RE.match(stripped)
                if header:
                    section_key = _agat_section_key(header.group(1))
                    entry = sections.setdefault(section_key, {}).setdefault(
                        species, {"all": {}, "collapsed": None}
                    )
                    target = entry["all"]
                    continue
                if entry is not None and "have isoforms!" in stripped:
                    entry["collapsed"] = {}
                    target = entry["collapsed"]
                    continue
                m = _AGAT_STAT_LINE_RE.match(stripped)
                if m and target is not None:
                    target[m.group(1)] = m.group(2)
    return sections


def _agat_ordered_metrics(species_map):
    """Union of metric names across species/views, metrics more species report
    coming first (ties broken by first appearance).

    AGAT's field set for a section isn't fully fixed across species - e.g. it
    inserts an extra "Number of pseudogene" count (shifting every later field
    over by one) when a GFF tags pseudogenes, and some sections' labels echo
    the source GFF's own feature-type spelling (e.g. "lnc_rna" vs "lncrna").
    Ordering by first-appearance alone would make the column order depend on
    whichever species happens to sort first, rather than on what most species
    actually share.
    """
    order = []
    seen = set()
    counts = {}
    for entry in species_map.values():
        for view in ("all", "collapsed"):
            for m in (entry.get(view) or {}):
                counts[m] = counts.get(m, 0) + 1
                if m not in seen:
                    seen.add(m)
                    order.append(m)
    rank = {m: i for i, m in enumerate(order)}
    return sorted(order, key=lambda m: (-counts[m], rank[m]))


def _agat_feature_sheets(paths):
    """Build a single AGAT sheet with one wide table per feature type, stacked
    vertically and separated by a title row naming the feature: species as
    rows, metrics as columns (NA where a species lacks that feature or
    metric). Feature types that report isoform-collapsed numbers get an
    extra "isoforms" column (all/collapsed) rather than a separate table,
    since Excel has no toggle.

    Every species passed in gets a row in every table, even ones that don't
    have that feature type at all - filled entirely with NA - rather than
    silently disappearing from the table.
    """
    all_species = sorted({Path(p).stem.replace(".stats", "") for p in paths})
    sections = _parse_agat_stats_nested(paths)
    combined_rows = []
    for key, species_map in sorted(sections.items()):
        has_iso = any(v["collapsed"] is not None for v in species_map.values())
        empty_entry = {"all": {}, "collapsed": None}

        metrics = _agat_ordered_metrics(species_map)

        header = (["assembly", "isoforms"] if has_iso else ["assembly"]) + metrics
        table_rows = [header]
        for species in all_species:
            entry = species_map.get(species, empty_entry)
            views = [("all", entry["all"])] + ([("collapsed", entry["collapsed"] or {})] if has_iso else [])
            for view_name, stats in views:
                row = [species] + ([view_name] if has_iso else [])
                row += [stats.get(m, "NA") for m in metrics]
                table_rows.append(row)

        label = _AGAT_SECTION_LABELS.get(key, key.replace("_", " ").capitalize())
        if combined_rows:
            combined_rows.append([])
        combined_rows.append([label])
        combined_rows.extend(table_rows)

    return {"Annotation_AGAT": combined_rows} if combined_rows else {}


def _quast_rows(paths):
    """Combine per-species QUAST report TSVs.

    QUAST TSV is transposed (metric rows, species column). Each file has
    two columns: metric name and value. We pivot to species-as-columns.
    """
    # Read all files first
    species_data = {}  # species -> {metric: value}
    metric_order = []
    for p in sorted(paths, key=lambda x: Path(x).name):
        species = Path(p).stem.replace(".quast", "")
        data = {}
        file_rows = list(_read_tsv(p, skip_comment=False))
        for row in file_rows:
            if len(row) >= 2:
                metric = row[0]
                val    = row[1]
                data[metric] = val
                if metric not in metric_order:
                    metric_order.append(metric)
        species_data[species] = data

    if not species_data:
        return []

    species_list = sorted(species_data.keys())
    header = ["metric"] + species_list
    rows = [header]
    for m in metric_order:
        rows.append([m] + [species_data[sp].get(m, "") for sp in species_list])
    return rows


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Generate an Excel workbook with raw GenomeQC result tables."
    )
    ap.add_argument("--busco_tables",         nargs="*", default=[], metavar="TSV",
                    help="BUSCO (genome) batch_summary_modified.txt files")
    ap.add_argument("--busco_prot_tables",    nargs="*", default=[], metavar="TSV",
                    help="BUSCO (proteins) batch_summary_modified.txt files")
    ap.add_argument("--busco_seqs_table",     default=None, metavar="TSV",
                    help="ortho_seqs.py output TSV (sequences above BUSCO threshold)")
    ap.add_argument("--quast_tsvs",           nargs="*", default=[], metavar="TSV",
                    help="Per-species QUAST report TSV files")
    ap.add_argument("--agat_stats",           nargs="*", default=[], metavar="TXT",
                    help="AGAT spstatistics *.txt files")
    ap.add_argument("--tidk_tsvs",            nargs="*", default=[], metavar="TSV",
                    help="tidk aposteriori search TSV files")
    ap.add_argument("--fcs_gx_reports",       nargs="*", default=[], metavar="TXT",
                    help="FCS-GX *.fcs_gx_report.txt files (optional)")
    ap.add_argument("--fcs_adaptor_reports",  nargs="*", default=[], metavar="TXT",
                    help="FCS-Adaptor *.fcs_adaptor_report.txt files (optional)")
    ap.add_argument("--tiara_reports",        nargs="*", default=[], metavar="TXT",
                    help="Tiara classification *.txt files (optional)")
    ap.add_argument("--repeatmasker_tbls",    nargs="*", default=[], metavar="TBL",
                    help="RepeatMasker *.tbl files (optional)")
    ap.add_argument("--output", default="genomeqc_tables.xlsx", metavar="XLSX",
                    help="Output Excel file (default: genomeqc_tables.xlsx)")
    args = ap.parse_args()

    wb = _XLSX()
    added = 0

    # ── Summary sheet (first) ─────────────────────────────────────────────────
    busco_rows_dicts      = _parse_busco_batch_summaries(args.busco_tables)      if args.busco_tables      else []
    busco_prot_rows_dicts = _parse_busco_batch_summaries(args.busco_prot_tables) if args.busco_prot_tables else []
    tidk_data_dicts  = _parse_tidk_tsvs(args.tidk_tsvs)               if args.tidk_tsvs  else {}
    busco_seqs_data, busco_seqs_col = (
        _parse_busco_seqs_table(args.busco_seqs_table)
        if args.busco_seqs_table else (None, None)
    )
    if busco_rows_dicts or busco_prot_rows_dicts or tidk_data_dicts:
        summary = _summary_rows(busco_rows_dicts, tidk_data_dicts, busco_seqs_data, busco_seqs_col,
                                busco_prot_rows_dicts)
        wb.add_sheet("Summary", summary)
        added += 1

    if args.busco_tables:
        rows = _combine_tsv_files(args.busco_tables, add_species_col=False)
        if rows:
            wb.add_sheet("BUSCO_Genome" if args.busco_prot_tables else "BUSCO", rows)
            added += 1

    if args.busco_prot_tables:
        rows = _combine_tsv_files(args.busco_prot_tables, add_species_col=False)
        if rows:
            wb.add_sheet("BUSCO_Proteins", rows)
            added += 1

    if args.quast_tsvs:
        rows = _quast_rows(args.quast_tsvs)
        if rows:
            wb.add_sheet("Assembly_QUAST", rows)
            added += 1

    if args.agat_stats:
        for sheet_name, rows in _agat_feature_sheets(args.agat_stats).items():
            wb.add_sheet(sheet_name, rows)
            added += 1

    if args.tidk_tsvs:
        rows = _combine_tsv_files(args.tidk_tsvs, add_species_col=True)
        if rows:
            wb.add_sheet("Telomeres_tidk", rows)
            added += 1

    if args.fcs_gx_reports:
        rows = _combine_tsv_files(args.fcs_gx_reports, add_species_col=True)
        if rows:
            wb.add_sheet("FCS_GX", rows)
            added += 1

    if args.fcs_adaptor_reports:
        rows = _combine_tsv_files(args.fcs_adaptor_reports, add_species_col=True)
        if rows:
            wb.add_sheet("FCS_Adaptor", rows)
            added += 1

    if args.tiara_reports:
        rows = _combine_tsv_files(args.tiara_reports, add_species_col=True)
        if rows:
            wb.add_sheet("Tiara", rows)
            added += 1

    if args.repeatmasker_tbls:
        rows = _repeatmasker_rows(args.repeatmasker_tbls)
        if rows:
            wb.add_sheet("Repeats_RepeatMasker", rows)
            added += 1

    if added == 0:
        print("WARNING: no input data found; writing empty workbook.", file=sys.stderr)
        wb.add_sheet("empty", [["no data"]])

    wb.save(args.output)
    print(f"Workbook written to {args.output} ({added} sheet(s))", file=sys.stderr)


if __name__ == "__main__":
    main()
