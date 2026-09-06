#!/usr/bin/env python3
"""Check inclusive numerator/denominator endpoints against independent tiny boxes."""
import argparse,subprocess
from fractions import Fraction as F
from math import gcd,isqrt
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];OUTPUT=ROOT/'artifacts/generated-results/elliptic-curves/pari_height_boundary_regression_v1.json'

def exact(coefficients,height):
    points=set()
    for d in range(1,height+1):
        for n in range(-height,height+1):
            if gcd(n,d)!=1:continue
            x=F(n,d);f=sum(v*x**i for i,v in enumerate(coefficients))
            if f>=0 and isqrt(f.numerator)**2==f.numerator and isqrt(f.denominator)**2==f.denominator:points.add(x)
    return sorted(map(str,points))

def main(check=False):
    if check:
        data=cert.read(OUTPUT)
        for row in data['rows']:
            xs=sorted(set(line[2:] for line in row['stdout'].splitlines() if line.startswith('X|')))
            if xs!=row['x_coordinates'] or xs!=exact(row['coefficients'],row['height']):raise ArithmeticError('boundary regression differs')
        print('REPLAYED FOUR INCLUSIVE BOUNDARY BOXES');return
    if OUTPUT.exists():raise FileExistsError('preserve boundary test')
    rows=[]
    for coefficients in ((-56,0,0,0,1),(-1,3,0,0,81)):
        for h in (2,3):
            poly='+'.join(f'({a})*x^{i}' for i,a in enumerate(coefficients));program=f'R=hyperellratpoints({poly},{h});for(i=1,#R,print("X|",R[i][1]));quit\n';p=subprocess.run(['/usr/bin/gp','-q'],input=program,text=True,capture_output=True,timeout=3)
            if p.returncode or '***' in p.stderr:raise ArithmeticError('GP boundary test failed')
            xs=sorted(set(line[2:] for line in p.stdout.splitlines() if line.startswith('X|')))
            if xs!=exact(coefficients,h):raise ArithmeticError('inclusive boundary mismatch')
            rows.append({'coefficients':list(coefficients),'height':h,'program':program,'stdout':p.stdout,'stderr':p.stderr,'returncode':p.returncode,'x_coordinates':xs})
    checkpoint(OUTPUT,{'schema':'elliptic-curves.pari-inclusive-height-boundaries.v1','status':'PASS','gp_binary_sha256':cert.hashed(Path('/usr/bin/gp')),'sources':{str(Path(__file__).resolve().relative_to(ROOT)):cert.hashed(Path(__file__).resolve())},'rows':rows,'claim_boundary':'Pinned PARI includes |numerator|=H and denominator=H in these independently enumerated nonsingular quartic boxes. This resolves the wording of its short help on these boundary regressions; production finite completeness still trusts the pinned implementation.'})
    print('FOUR INCLUSIVE BOUNDARY BOXES PASS')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');main(p.parse_args().check)
