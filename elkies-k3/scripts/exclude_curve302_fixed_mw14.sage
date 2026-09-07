#!/usr/bin/env sage-python
"""Exact projective inverse test of the certified fixed3A1/MW14 fibration.

One equation, one prime, no parameter height bound; 60-second time cap.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import signal
import sys
from sage.all import EllipticCurve,GF,PolynomialRing,QQ

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'artifacts/local/elkies-k3/fixed-reverse-3a1-rr-qq.json'
PUBLIC=ROOT/'elliptic-curves/cas/icarm_curve302.py'
OUT=ROOT/'artifacts/generated-results/elkies-k3-curve302-fixed-3a1-mw14-exclusion-v1.json'
sys.path.insert(0,str(PUBLIC.parent))
from icarm_curve302 import GENERAL_WEIERSTRASS_COEFFICIENTS


def digest(p):return sha256(p.read_bytes()).hexdigest()


def build():
    source=json.loads(SOURCE.read_text());assert source['status']=='PASS_EXACT_QQ_FIXED_REVERSE_3A1_RR_JACOBIAN'
    child=source['child'];assert child['ADE']=='3A1' and child['MW_rank_if_rho19']==14
    assert child['root_rank']==3 and child['degrees_A_B_Delta']==[8,12,24]
    p=1009;F=GF(p);R=PolynomialRing(F,'t');t=R.gen()
    a=R([F(QQ(c)) for c in child['minimal_A_coefficients_low_to_high']])
    b=R([F(QQ(c)) for c in child['minimal_B_coefficients_low_to_high']])
    assert a.degree()==8 and b.degree()==12
    target=EllipticCurve(QQ,list(map(QQ,GENERAL_WEIERSTRASS_COEFFICIENTS)))
    c4,c6=target.c_invariants();ta=F(-27*c4);tb=F(-54*c6)
    comparison=a**3*tb**2-ta**3*b**2
    assert comparison.degree()==24 and comparison.gcd(t**p-t)==1
    values=[int(comparison(v)) for v in F];infinity=int(comparison[24]);assert infinity and all(values)
    return {'schema':'curve302.fixed-mw14-inverse.v1','status':'PASS_NO_RATIONAL_PARAMETER_ON_FIXED_3A1_MW14',
        'input_sha256':{str(f.relative_to(ROOT)):digest(f) for f in [Path(__file__),SOURCE,PUBLIC]},
        'parent_rank_authority':'EC-K3-H3-FIXED-REVERSE-3A1-QQ',
        'family':'fixed3A1/MW14 on the source-identified determinant948 K3',
        'limits':{'equations':1,'prime':p,'seconds':60},
        'A_mod_p_coefficients_low_to_high':list(map(int,a.list())),
        'B_mod_p_coefficients_low_to_high':list(map(int,b.list())),
        'target_short_A_B_mod_p':[int(ta),int(tb)],
        'comparison_mod_p_coefficients_low_to_high':list(map(int,comparison.list())),
        'projective_degree':24,'finite_values_in_residue_order':values,'infinity_value':infinity,
        'boundary':'This one certified MW14 fibration has no rational parameter with j equal to302. Fixed homogeneous degree24 retains infinity and parameters whose denominator is divisible by1009. No exclusion of other MW14,MW15,MW16 fibrations or all parents is asserted.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');args=parser.parse_args();signal.alarm(60)
    result=build()
    if args.check:assert result==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'])
