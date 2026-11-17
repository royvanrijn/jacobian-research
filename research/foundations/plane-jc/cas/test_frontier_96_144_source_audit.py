#!/usr/bin/env python3
"""Regression tests for the repeated-tail source reconciliation."""

from frontier_96_144_source_audit import repeated_tail_96_144_audit


audit = repeated_tail_96_144_audit()

assert audit.degree_pair == (96, 144)
assert audit.final_edge_weight == (4, -1)
assert audit.forcing_fraction == (3, 4)
assert audit.q1 == audit.d0_upper_bound == audit.d0 == 4
assert audit.primitive_vertical_factor_endpoints == ((2, 7), (2, 10))
assert audit.residual_factor_degree == 3

by_partition = {item.partition: item for item in audit.root_partitions}
assert set(by_partition) == {(3,), (2, 1), (1, 1, 1)}
assert by_partition[(2, 1)].excluded_by_forbidden_corner
assert by_partition[(1, 1, 1)].excluded_by_forbidden_corner
assert not by_partition[(3,)].excluded_by_forbidden_corner
assert audit.surviving_partitions == ((3,),)
assert audit.surviving_factor == "R=kappa*x^2*y^7*(y-lambda)^3"
assert audit.triple_root_open_chain_counts == (1, 6, 3, 0)
assert audit.triple_root_complete_chain_count == 0
assert audit.repeated_tail_excluded
assert audit.status.startswith("excluded")

print("PASS: the common tail forces q1=d0=4")
print("PASS: the (8,40) vertical factor has residual degree three")
print("PASS: every partition with a simple root yields the forbidden corner (8,4)")
print("PASS: the sole triple-root branch has no complete-chain escape")
print("PASS: the repeated-tail (96,144) row is excluded")
