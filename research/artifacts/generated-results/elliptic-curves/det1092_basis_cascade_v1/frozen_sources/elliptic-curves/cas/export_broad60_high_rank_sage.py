#!/usr/bin/env python3
"""Executable equations and exactly transported points for the unmatched broader higher MW16 fibres."""
import argparse
import certify_compact_r17_candidates as cert
import certify_broad60_high_rank_minimal as proof

OUT = proof.ART / 'new_broad60_high_rank_curves.sage'

def main(check):
    data = proof.expected()
    if cert.read(proof.OUT) != data:
        raise ArithmeticError('exact minimal-model certificate differs')
    lines = [f'# {len(data["curves"])} distinct higher MW16 curves with certified rank lower bounds.',
             '# Proof: broad60_high_rank_minimal_proof_v1.json.',
             '# These are not exact-rank or universal-novelty claims.',
             'from sage.all import QQ, EllipticCurve', 'curves = {}']
    for row in data['curves']:
        lines.extend([f"E = EllipticCurve(QQ, {row['minimal_curve']!r})",
                      f"points = [E([QQ(x), QQ(y)]) for x, y in {row['points']!r}]",
                      f"assert len(points) == {row['rank_lower_bound']}",
                      f"curves[{row['id']!r}] = (E, points)",
                      f"print({row['id']!r}, E, 'rank >= {row['rank_lower_bound']}')"])
    text = '\n'.join(lines) + '\n'
    if check:
        if OUT.read_text() != text:
            raise ArithmeticError('executable equation/point export differs')
    else:
        with OUT.open('x') as stream:
            stream.write(text)
    print('BROADER HIGHER MW16 EQUATION/POINT EXPORTS PASS', flush=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    main(parser.parse_args().check)
