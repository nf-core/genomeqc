#!/usr/bin/env python3

# Written by Fernando Duarte with AI assistance, released under the MIT license.
# Parses RepeatMasker .tbl files into a combined TSV for tree-summary pie charts
"""
Parse one or more RepeatMasker .tbl files and produce a combined TSV
suitable for pie-chart plotting in the tree summary (scatterpie).

Output columns:
  species  SINE  LINE  LTR  Penelope  DNA  Rolling_Circle  Unclassified  Other  Non_Repeat

All values are percentages (0-100) summing to 100 per row.
'Other' aggregates Small RNA, Satellites, Simple repeats, Low complexity,
and any unassigned retroelements.

Usage:
  te_tbl_2_table.py species1.tbl species2.tbl -o te_summary.tsv
"""

import argparse
import os
import re
import sys

# ---------------------------------------------------------------------------
# Label → output-column mapping (matched against stripped line prefix).
# Order matters: sub-class patterns must come before their parent category.
# ---------------------------------------------------------------------------

_PATTERNS = [
    # Retroelement sub-classes
    (re.compile(r"SINEs?:",                     re.IGNORECASE), "SINE"),
    (re.compile(r"LINEs?:",                     re.IGNORECASE), "LINE"),
    (re.compile(r"LTR elements?:",              re.IGNORECASE), "LTR"),
    (re.compile(r"Penelope:",                   re.IGNORECASE), "Penelope"),
    (re.compile(r"Other retroelements?:?",      re.IGNORECASE), "Other"),
    # Top-level categories
    (re.compile(r"DNA transposons",             re.IGNORECASE), "DNA"),
    (re.compile(r"Rolling.circles",             re.IGNORECASE), "Rolling_Circle"),
    (re.compile(r"Unclassified:?",              re.IGNORECASE), "Unclassified"),
    # Aggregated into "Other"
    (re.compile(r"Small RNA:?",                 re.IGNORECASE), "Other"),
    (re.compile(r"Satellites?:?",               re.IGNORECASE), "Other"),
    (re.compile(r"Simple repeats?:?",           re.IGNORECASE), "Other"),
    (re.compile(r"Low complexity:?",            re.IGNORECASE), "Other"),
]

_TE_COLS = ["SINE", "LINE", "LTR", "Penelope",
            "DNA", "Rolling_Circle", "Unclassified", "Other"]

_ALL_COLS = _TE_COLS + ["Non_Repeat"]

_PCT_RE = re.compile(r"(\d+\.\d+)\s*%\s*$")


def parse_tbl(path: str) -> dict:
    """Return {column: percentage} for one .tbl file."""
    totals = {c: 0.0 for c in _TE_COLS}

    with open(path) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            m_pct = _PCT_RE.search(stripped)
            if m_pct is None:
                continue
            pct = float(m_pct.group(1))
            for pattern, col in _PATTERNS:
                if pattern.match(stripped):
                    totals[col] += pct
                    break  # first match only – avoids double-counting sub-classes

    total_te = sum(totals.values())
    totals["Non_Repeat"] = max(0.0, round(100.0 - total_te, 4))
    return totals


def species_name(path: str) -> str:
    """Strip the final extension from a filename to get the species name."""
    return os.path.splitext(os.path.basename(path))[0]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("tbl_files", nargs="+", help="RepeatMasker .tbl files (one per species)")
    ap.add_argument("-o", "--output", default="-", help="Output TSV (default: stdout)")
    args = ap.parse_args()

    rows = []
    for path in args.tbl_files:
        try:
            data = parse_tbl(path)
        except Exception as exc:
            print(f"WARNING: skipping {path}: {exc}", file=sys.stderr)
            continue
        sp = species_name(path)
        rows.append([sp] + [str(round(data[c], 4)) for c in _ALL_COLS])

    out = open(args.output, "w") if args.output != "-" else sys.stdout
    try:
        print("\t".join(["species"] + _ALL_COLS), file=out)
        for row in rows:
            print("\t".join(row), file=out)
    finally:
        if args.output != "-":
            out.close()


if __name__ == "__main__":
    main()
