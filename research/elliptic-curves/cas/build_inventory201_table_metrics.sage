#!/usr/bin/env sage-python
"""Minimal models, transported points and ICARM-compatible table metrics."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
from sage.all import EllipticCurve, QQ, ZZ, RealField, pari
from sage.version import version

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert
ART = ROOT/'artifacts/generated-results/elliptic-curves'
INPUT = ART/'new_high_rank_curve_index_v22.json'
CAT = ROOT/'artifacts/local/elliptic-curves/conductor-record-screen-v1/database.json'
AUDIT = ART/'inventory201_conductor_bounds_v1.json'
SOURCE = ROOT/'artifacts/local/elliptic-curves/icarm-table-source'
OUT = ART/'inventory201_table_metrics_v1.json'


def run(check=False):
    audit = {r['id']:r for r in cert.read(AUDIT)['rows']}
    public = {r['id']:r for r in cert.read(CAT)['curves']}
    R = RealField(160)
    rows = []
    for r in cert.read(INPUT)['curves']:
        prior = audit[r['id']]
        original = EllipticCurve(QQ,[QQ(x) for x in r['curve']])
        raw,change = pari.ellinit([QQ(x) for x in prior['integral_curve']]).ellminimalmodel()
        model = [QQ(x) for x in raw[:5]]
        E = EllipticCurve(QQ,model)
        iso = original.isomorphism_to(E)
        c4,c6 = E.c_invariants()
        if (c4,c6,E.discriminant()) != tuple(QQ(raw[i]) for i in (9,10,11)):
            raise ArithmeticError('independent exact invariants disagree')
        g = ZZ(c4).gcd(ZZ(c6)).abs()
        if not 0 < g <= 10000000:
            raise ArithmeticError('minimality GCD exceeds finite trial-factor gate')
        primes = sorted(set(g.prime_divisors())|{ZZ(2),ZZ(3)})
        local = []
        for p in primes:
            data = E.local_data(p,algorithm='generic',proof=True)
            if data.discriminant_valuation() != E.discriminant().valuation(p):
                raise ArithmeticError('model is not minimal')
            local.append({'prime':str(p),'minimal_discriminant_valuation':int(data.discriminant_valuation())})
        points = []
        for pair in r['points']:
            P = original([QQ(x) for x in pair])
            Q = iso(P)
            if Q.is_zero() or Q not in E:
                raise ArithmeticError('point transport failed')
            points.append([str(Q[0]),str(Q[1])])
        naive = R(max(abs(c4)**3,c6**2)).log()
        h = E.faltings_height(stable=False,prec=160)
        low = E.faltings_height(stable=False,prec=96)
        if abs(R(low)-h) > R('1e-22'):
            raise ArithmeticError('period-area height precision check failed')
        matches = prior['catalogue_matches']
        for i in matches:
            q = public[i]
            if not cert.isomorphic(list(map(str,model)),q['ainvs']) or str(E.discriminant()) != q['discriminant']:
                raise ArithmeticError('public minimal-model comparison failed')
            if abs(float(naive)-q['naive_height']) > 1e-9 or abs(float(h)-q['faltings_height']) > 1e-9:
                raise ArithmeticError('ICARM metric convention differs')
        rows.append({'id':r['id'],'ainvs':list(map(str,model)),
            'c4':str(c4),'c6':str(c6),'discriminant':str(E.discriminant()),
            'naive_height':float(naive),'faltings_height':float(h),
            'log_abs_discriminant':float(R(abs(E.discriminant())).log()),
            'original_to_minimal_isomorphism':list(map(str,iso.tuple())),
            'points':points,'invariant_gcd':str(g),'minimality_local_checks':local,
            'catalogue_matches':matches})
    result = {'schema':'elliptic-curves.inventory201-table-metrics.v1','status':'PASS',
        'rows':rows,'point_memberships_and_transports':sum(len(r['points']) for r in rows),
        'public_metric_matches':sum(len(r['catalogue_matches']) for r in rows),
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__),INPUT,AUDIT,CAT,SOURCE/'src/verify.ts']},
        'icarm_source_commit':subprocess.check_output(['git','-C',str(SOURCE),'rev-parse','HEAD'],text=True).strip(),
        'sage_version':version,'precision_bits':[96,160],
        'definitions':{'log':'natural logarithm','naive_height':'log max(abs(c4)^3,c6^2) on global minimal model',
            'faltings_height':'-1/2 log(period lattice area) on global minimal model, matching ICARM; Sage stable=False',
            'minimality':'Generic local minimization at2,3 and every prime dividing gcd(c4,c6), whose factorization is bounded by10^7. Outside these primes at least one invariant is a unit, precluding a smaller integral model.'},
        'claim_boundary':'Exact minimal equations, discriminants and point transports. Height and logarithm columns are numerical approximations with precision checks; ICARM period-area convention is used, not Sage stable=True. Independence is inherited through the exact isomorphism from the frozen rank certificates; no new rank claim or parameter search.'}
    if check:
        if cert.read(OUT) != result:
            raise ArithmeticError('table metric replay differs')
    elif OUT.exists():
        if cert.read(OUT) != result:
            raise FileExistsError('preserve table metrics')
    else:
        cert.write(OUT,result)
    print('TABLE METRICS PASS',len(rows),'curves;',result['point_memberships_and_transports'],'point transports;',result['public_metric_matches'],'ICARM matches')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    run(parser.parse_args().check)
