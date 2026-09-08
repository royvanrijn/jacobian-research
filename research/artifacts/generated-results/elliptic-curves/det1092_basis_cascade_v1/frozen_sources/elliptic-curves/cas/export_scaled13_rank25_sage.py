#!/usr/bin/env python3
"""Exact standalone export of the newly certified 13-scaled R17 rank25 point subgroup."""
import argparse
import certify_compact_r17_candidates as cert
from certify_scaled13_24_rank25_minimal import ART,OUT
PATH=ART/'new_scaled13_rank25_curve_11952.sage'
def main(check):
    d=cert.read(OUT)
    if d['status']!='PASS' or len(d['curves'])!=1:raise ArithmeticError('fixed single new25 proof required')
    r=d['curves'][0]
    if r['id']!='11952-300' or r['rank_lower_bound']!=25 or r['icarm_matches'] or r['previous_matches']:raise ArithmeticError('certified unmatched25 roster differs')
    text='# Certified rank lower bound25, not an exact-rank or universal-novelty claim.\n# Minimal model and independent-point proof: scaled13_24_rank25_minimal_proof_v1.json.\nfrom sage.all import QQ, EllipticCurve\nE = EllipticCurve(QQ, '+repr(r['minimal_curve'])+')\npoints = [E([QQ(x), QQ(y)]) for x, y in '+repr(r['points'])+']\nassert len(points) == 25\nprint(E)\nprint("25 exactly certified independent rational points")\n'
    if check:
        if PATH.read_text()!=text:raise ArithmeticError('point export differs')
    else:
        with PATH.open('x') as f:f.write(text)
    print('EXACT NEW 13-scaled R17 RANK25 SAGE EXPORT',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();main(a.check)
