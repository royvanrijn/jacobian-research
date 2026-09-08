#!/usr/bin/env python3
"""Deterministic Sage point files for the three catalogue-unmatched MW16 rank26s."""
import argparse
import certify_compact_r17_candidates as cert
from certify_extended20_mw16_minimal import ART,OUT

def main(check):
    d=cert.read(OUT)
    if d['status']!='PASS' or len(d['curves'])!=3:raise ArithmeticError('fixed three certified curves required')
    for r in d['curves']:
        if r['rank_lower_bound']!=26 or r['icarm_matches'] or r['previous_matches']:raise ArithmeticError('unmatched26 roster differs')
        path=ART/('new_mw16_rank26_'+r['id'].replace('-','_')+'.sage')
        text='# Certified rank lower bound; no exact rank or universal novelty claim.\n# No Q-isomorphism match in the pinned593-equation catalogue or382 prior equations.\nfrom sage.all import QQ, EllipticCurve\nE = EllipticCurve(QQ, '+repr(r['minimal_curve'])+')\npoints = [E([QQ(x), QQ(y)]) for x, y in '+repr(r['points'])+']\nassert len(points) == 26\nprint(E)\nprint("26 exactly certified independent points; see extended20_mw16_high_rank_minimal_v1.json")\n'
        if check:
            if path.read_text()!=text:raise ArithmeticError('Sage export differs')
        else:
            with path.open('x') as f:f.write(text)
    print('EXACT THREE MW16 SAGE EXPORTS',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();main(a.check)
