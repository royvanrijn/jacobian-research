#!/usr/bin/env python3
"""Export the proved minimal curve and its 26 transported independent points."""
import argparse
from pathlib import Path
import certify_new_compact_rank26_v2 as proof


def content():
    data = proof.cert.read(proof.OUTPUT)
    proof.verify(data)
    lines = [
        '# New compact R17 curve, stable ID new-20260905-37.',
        '# Rank >=26, not an exact-rank assertion. Global minimal model.',
        '# Exact independence proof: new_compact_rank26_proof_v1.json.',
        '# Replay with elliptic-curves/cas/certify_new_compact_rank26_v2.py.',
        'from sage.all import EllipticCurve, QQ',
        'E = EllipticCurve(QQ, ['+', '.join(data['minimal_curve'])+'])',
        'point_coordinates = [',
        *['    ('+repr(x)+', '+repr(y)+'),' for x, y in data['points']],
        ']',
        'points = [E(QQ(x), QQ(y)) for x, y in point_coordinates]',
        'assert len(points) == 26',
        'assert E.c4().gcd(E.c6()) == 27',
        'print("Loaded global minimal curve and 26 certified independent points:", E)',
        '',
    ]
    return '\n'.join(lines)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=proof.ART/'new_compact_rank26_curve.sage')
    p.add_argument('--check', action='store_true')
    a = p.parse_args(); text = content()
    if a.check:
        if a.output.read_text() != text:
            raise ArithmeticError('Sage export differs')
    else:
        with a.output.open('x') as stream:
            stream.write(text)
