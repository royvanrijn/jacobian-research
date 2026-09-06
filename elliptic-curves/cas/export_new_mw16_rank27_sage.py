#!/usr/bin/env python3
"""Exact standalone export of the newly certified MW16 rank27 point subgroup."""
import argparse
import certify_compact_r17_candidates as cert
from certify_new_mw16_rank27_minimal import ART,OUT
PATH=ART/'new_mw16_rank27_curve_a1_01.sage'
def main(check):
    d=cert.read(OUT)
    if d['status']!='PASS' or len(d['curves'])!=1:raise ArithmeticError('fixed single new27 proof required')
    r=d['curves'][0]
    if r['id']!='a1-fibration-01-052' or r['rank_lower_bound']!=27 or r['icarm_matches'] or r['previous_matches']:raise ArithmeticError('certified unmatched27 roster differs')
    text='# Certified rank lower bound27, not an exact-rank or universal-novelty claim.\n# Minimal model and independent-point proof: new_mw16_rank27_minimal_proof_v1.json.\nfrom sage.all import QQ, EllipticCurve\nE = EllipticCurve(QQ, '+repr(r['minimal_curve'])+')\npoints = [E([QQ(x), QQ(y)]) for x, y in '+repr(r['points'])+']\nassert len(points) == 27\nprint(E)\nprint("27 exactly certified independent rational points")\n'
    if check:
        if PATH.read_text()!=text:raise ArithmeticError('point export differs')
    else:
        with PATH.open('x') as f:f.write(text)
    print('EXACT NEW MW16 RANK27 SAGE EXPORT',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();main(a.check)
