#!/usr/bin/env python3

# Written by Fernando Duarte with AI assistance, released under the MIT license.
# Converts famdb EMBL output to RepeatMasker-compatible FASTA
"""Convert famdb EMBL output to RepeatMasker-compatible FASTA.

Input:  famdb families -f embl  (on stdin)
Output: FASTA with headers >Accession#Type/SubType  (on stdout)

The #Type/SubType suffix is read from the RepeatMasker Annotations
block in each EMBL record's CC lines.
"""
import sys
import re


def convert(fh, out):
    accession = rm_type = rm_subtype = ''
    seq_lines = []
    in_seq = in_rm = False

    for line in fh:
        if line.startswith('ID '):
            accession = line.split()[1].rstrip(';')
            rm_type = rm_subtype = ''
            seq_lines = []
            in_seq = in_rm = False

        elif line.startswith('CC'):
            if 'RepeatMasker Annotations:' in line:
                in_rm = True
            elif in_rm:
                field = line[2:].strip()
                if field.startswith('Type:') and 'SubType' not in field:
                    rm_type = field.split('Type:', 1)[1].strip()
                elif field.startswith('SubType:'):
                    rm_subtype = field.split('SubType:', 1)[1].strip()

        elif line.startswith('XX'):
            in_rm = False

        elif line.startswith('SQ'):
            in_seq = True

        elif line.startswith('//'):
            if accession and seq_lines:
                if rm_type and rm_subtype:
                    header = f'>{accession}#{rm_type}/{rm_subtype}'
                elif rm_type:
                    header = f'>{accession}#{rm_type}'
                else:
                    header = f'>{accession}'
                seq = ''.join(seq_lines).upper()
                out.write(header + '\n')
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i+60] + '\n')
            in_seq = False

        elif in_seq:
            seq_lines.append(re.sub(r'[\s0-9]', '', line))


def main():
    try:
        convert(sys.stdin, sys.stdout)
    except BrokenPipeError:
        sys.stderr.close()


if __name__ == '__main__':
    main()
