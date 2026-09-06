#!/usr/bin/env python3
"""Exact standalone export of the newly certified higher-height R17 rank26 point subgroup."""
import argparse
import certify_compact_r17_candidates as cert
from certify_higher26_minimal import ART,OUT
PATH=ART/'new_higher_rank26_curve_11952.sage'
def main(check):
    d=cert.read(OUT)
    if d['status']!='PASS' or len(d['curves'])!=1:raise ArithmeticError('fixed single new26 proof required')
    r=d['curves'][0]
    if r['id']!='11952-069' or r['rank_lower_bound']!=26 or r['icarm_matches'] or r['previous_matches']:raise ArithmeticError('certified unmatched26 roster differs')
    text='# Certified rank lower bound26, not an exact-rank or universal-novelty claim.\n# Minimal model and independent-point proof: higher26_minimal_proof_v1.json.\nfrom sage.all import QQ, EllipticCurve\nE = EllipticCurve(QQ, '+repr(r['minimal_curve'])+')\npoints = [E([QQ(x), QQ(y)]) for x, y in '+repr(r['points'])+']\nassert len(points) == 26\nprint(E)\nprint("26 exactly certified independent rational points")\n'
    if check:
        if PATH.read_text()!=text:raise ArithmeticError('point export differs')
    else:
        with PATH.open('x') as f:f.write(text)
    print('EXACT NEW higher-height R17 RANK26 SAGE EXPORT',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();main(a.check)
