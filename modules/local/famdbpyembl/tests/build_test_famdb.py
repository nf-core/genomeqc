#!/usr/bin/env python3

# Written by Chris Wyatt with AI assistance, released under the MIT license.
#
# Builds test.0.h5, the tiny synthetic famdb fixture used by this module's
# nf-test. Real Dfam famdb databases are large (tens of MB minimum for just
# the root/taxonomy partition), so this constructs a minimal, schema-valid
# database directly via famdb.py's own internal classes rather than
# downloading a real one.
#
# famdb.py's on-disk schema is versioned and has changed across releases
# (see https://github.com/Dfam-consortium/FamDB); this script targets famdb
# format version 2.0.5, the version bundled with the pipeline's pinned
# RepeatMasker container (quay.io/biocontainers/repeatmasker:4.2.2--pl5321hdfd78af_0).
# If that container is ever upgraded to a RepeatMasker/famdb.py release with
# a different schema, this script (and the checked-in test.0.h5) will need
# regenerating against the new famdb_classes.py/famdb_helper_classes.py.
#
# Usage (run inside the RepeatMasker container, which already has the
# bundled famdb modules and a matching h5py on its Python install):
#   docker run --rm -v "$(pwd):/data" \
#     quay.io/biocontainers/repeatmasker:4.2.2--pl5321hdfd78af_0 \
#     python /data/build_test_famdb.py
#
# Verified against the real famdb.py CLI:
#   famdb.py -i <dir> info
#   famdb.py -i <dir> lineage -ad 3
#   famdb.py -i <dir> families --descendants 3 -f embl

import sys
import os

sys.path.insert(0, "/usr/local/share/RepeatMasker")
from famdb_classes import FamDBRoot, FamDB  # noqa: E402
from famdb_helper_classes import TaxNode, Family  # noqa: E402

OUT_DIR = os.path.dirname(os.path.abspath(__file__))

# Minimal 3-node taxonomy: root(1) -> order(2) -> species(3).
# The test family below is attached to taxon 3.
TAX_DB = {
    1: TaxNode(1, None),
    2: TaxNode(2, 1),
    3: TaxNode(3, 2),
}
for tax_id, node in TAX_DB.items():
    if tax_id != 1:
        node.parent_node = TAX_DB[node.parent_id]
        node.parent_node.children.append(node)
    node.names = [["scientific name", f"Testus taxon{tax_id}"]]

family = Family()
family.name = "TESTFAM1"
family.accession = "TEST0000001"
family.version = 1
family.consensus = "ACGTACGTACGTACGTACGTACGTACGTACGT"
family.length = 32
family.classification = "root;Interspersed_Repeat;Transposable_Element;Class_I;LTR;Ty1_Copia"
family.clades = [3]

FILE_INFO = {
    "meta": {"uuid": "test-uuid-0001", "db_version": "TESTV1", "db_date": "2026-01-01"},
    "file_map": {
        "0": {
            "T_root": 1,
            "filename": "test.0.h5",
            "F_roots": [],
            "T_root_name": "root",
            "F_roots_names": [],
        }
    },
}

# Single-partition layout: all taxa live in partition 0.
PARTITION_NODES = {0: [1, 2, 3]}

with FamDBRoot(os.path.join(OUT_DIR, "test.0.h5"), "w") as db:
    db.set_metadata("0", FILE_INFO, "Test Dfam", "TESTV1", "2026-01-01", "Test copyright")
    # write_full_taxonomy builds the parent/child tree structure;
    # write_taxonomy separately creates the Lookup/ByTaxon group entries
    # that add_family() requires to already exist before it will link a
    # family to a taxon (it silently no-ops the link otherwise).
    db.write_full_taxonomy(TAX_DB, PARTITION_NODES)
    db.write_taxonomy(TAX_DB.keys())
    db.add_family(family)
    db.finalize()

print(f"Built {OUT_DIR}/test.0.h5")

# families --descendants queries walk a separate "pruned" tree (the
# Val_Children/Val_Parent datasets), which is only populated by reopening
# the database and running build_pruned_tree() as a distinct step.
db_dir = FamDB(OUT_DIR, "r+")
db_dir.build_pruned_tree()
print("Pruned tree built")
