#!/usr/bin/env python3

# Written by Chris Wyatt with AI assistance, released under the MIT license.
# Parses FCS-GX genome reports into a combined TSV for tree-summary pie/ring plotting
"""
Parse one or more NCBI FCS-GX `*.fcs_gx_report.txt` genome reports and produce
a combined TSV suitable for pie-chart/ring plotting in the tree summary.

Each report's first line is a `##`-prefixed JSON header containing
`run-info.agg-cvg`, the fraction of the genome's total bases covered by
putative-contaminant hits - this is used directly as the contamination
percentage rather than summing the per-sequence hit rows below it, since
those rows only cover flagged sequences (not the whole genome) and an empty
report (nothing flagged) has no rows at all. An empty/missing report means
0% contamination.

Output columns:
  species  contaminant_pct  non_contaminant_pct

Values are percentages (0-100) summing to 100 per row.

Usage:
  fcsgx_report_2_table.py species1.fcs_gx_report.txt species2.fcs_gx_report.txt -o fcs_summary.tsv
"""

import argparse
import json
import os
import re
import sys

_SUFFIX_RE = re.compile(r"\.fcs_gx_report\.txt$")


def species_name(path: str) -> str:
    """Strip the `.fcs_gx_report.txt` suffix from a filename to get the species name."""
    return _SUFFIX_RE.sub("", os.path.basename(path))


def parse_report(path: str) -> float:
    """Return the contaminant percentage (0-100) for one FCS-GX report."""
    if os.path.getsize(path) == 0:
        return 0.0

    with open(path) as fh:
        header = fh.readline().strip()

    if not header.startswith("##"):
        raise ValueError(f"unexpected header line: {header[:80]!r}")

    data = json.loads(header[2:])
    agg_cvg = data[1]["run-info"]["agg-cvg"]
    return round(agg_cvg * 100, 4)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("report_files", nargs="+", help="FCS-GX *.fcs_gx_report.txt files (one per species)")
    ap.add_argument("-o", "--output", default="-", help="Output TSV (default: stdout)")
    args = ap.parse_args()

    rows = []
    for path in args.report_files:
        try:
            contaminant_pct = parse_report(path)
        except Exception as exc:
            print(f"WARNING: skipping {path}: {exc}", file=sys.stderr)
            continue
        non_contaminant_pct = round(100.0 - contaminant_pct, 4)
        rows.append([species_name(path), str(contaminant_pct), str(non_contaminant_pct)])

    out = open(args.output, "w") if args.output != "-" else sys.stdout
    try:
        print("\t".join(["species", "contaminant_pct", "non_contaminant_pct"]), file=out)
        for row in rows:
            print("\t".join(row), file=out)
    finally:
        if args.output != "-":
            out.close()


if __name__ == "__main__":
    main()
