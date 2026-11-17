#!/usr/bin/env python3
"""Turn every independent signed-lift class over an S4 quotient into a fibre.

Bit coordinates 0,1 are the affine plane F2^2. Higher bits are the radical.
The common alternating form is B(x,y)=x0*y1+x1*y0, of rank two.
"""

import argparse
import json
from pathlib import Path

from signed_lifts import contains, nullspace, relator_rows, row_basis, solve


def polar(x, y):
    return (((x & 1) * ((y >> 1) & 1)) ^ (((x >> 1) & 1) * (y & 1))) & 1


def apply_generator(x, direction, offset):
    return x ^ direction if polar(direction, x) ^ offset else x


def build(data):
    assert data["degree"] == 4, "the affine-plane tower uses a specified S4 quotient"
    result = solve(data)
    edges = [tuple(e) for e in data["transpositions"]]
    s = len(edges)
    rows = relator_rows(4, edges, data["relators"])
    lifts = nullspace(rows, s)
    cuts = [sum(((p == v) ^ (q == v)) << i for i, (p, q) in enumerate(edges)) for v in range(4)]
    spanning = list(row_basis(cuts).values())
    quotient_basis = []
    for a in lifts:
        if not contains(row_basis(spanning), a):
            quotient_basis.append(a)
            spanning.append(a)
    d = len(quotient_basis)
    assert d == result["lift_space_dimension"] - 3
    labels = [sum(((a >> i) & 1) << j for j, a in enumerate(quotient_basis)) for i in range(s)]
    directions = [(p ^ q) | (labels[i] << 2) for i, (p, q) in enumerate(edges)]
    offsets = [1 ^ polar(p ^ q, p) for p, q in edges]
    output = {"label": data.get("label"), "meridional_rank_lower_bound": d + 3,
              "dimension": d + 2, "sheets": 1 << (d + 2), "form_rank": 2,
              "radical_dimension": d, "image_order": 24 * 8 ** d,
              "image": f"Aff(F2^2,F2^{d}) semidirect S4",
              "quotient_basis_bits": quotient_basis, "radical_edge_labels_bits": labels,
              "directions_bits": directions, "offsets": offsets,
              "conditional_on_presentation_and_meridian_provenance": True}
    if d <= 12:
        seen, stack = {0}, [0]
        while stack:
            x = stack.pop()
            for u, a in zip(directions, offsets):
                y = apply_generator(x, u, a)
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        assert len(seen) == output["sheets"]
        output["independent_orbit_size"] = len(seen)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("presentation", type=Path)
    args = parser.parse_args()
    print(json.dumps(build(json.loads(args.presentation.read_text())), indent=2))


if __name__ == "__main__":
    main()
