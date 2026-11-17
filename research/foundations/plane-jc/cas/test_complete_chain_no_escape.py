#!/usr/bin/env python3
"""Regression test for the triple-root complete-chain obstruction."""

from complete_chain_no_escape import (
    Corner,
    Edge,
    complete_chains,
    maximum_chain_length,
    triple_root_no_escape_audit,
)


audit = triple_root_no_escape_audit()
assert audit.translated_initial_edge == (("8", 40), ("8", 12))
assert audit.maximum_length == 3
assert audit.open_chain_counts == (1, 6, 3, 0)
assert audit.complete_chain_count == 0
assert audit.published_fixture_reproduced
assert audit.excluded

fixture = Edge(Corner(8, 1, 40), Corner(8, 1, 28))
assert maximum_chain_length(fixture) == 3
fixture_chains, _ = complete_chains(fixture)
assert any(chain[-1].end == Corner(11, 4, 7) for chain in fixture_chains)

print("PASS: the published repeated tail reaches final corner (11/4,7)")
print("PASS: the triple-root translation has maximum chain length three")
print("PASS: permissive open-chain counts are exactly 1,6,3,0")
print("PASS: the triple-root repeated-tail row has no complete-chain escape")
