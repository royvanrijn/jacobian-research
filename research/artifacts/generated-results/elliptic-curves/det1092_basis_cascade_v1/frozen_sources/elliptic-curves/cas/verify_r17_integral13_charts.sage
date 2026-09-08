#!/usr/bin/env sage-python
"""Independent exact symbolic substitution checks for the determinant13 charts."""
import sys
from pathlib import Path
from sage.all import PolynomialRing,QQ,matrix,ZZ
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert
import audit_r17_integral13_charts as audit

def main():
    d=cert.read(audit.OUT);families={f['family']:f for f in cert.read(audit.INPUT)['families']};R=PolynomialRing(QQ,names=('u','v'));u,v=R.gens()
    if d['status']!='PASS' or len(d['rows'])!=84:raise ArithmeticError('complete84-cell coefficient audit required')
    for r in d['rows']:
        f=families[r['family']];m=matrix(ZZ,2,2,r['matrix']);n=m[0,0]*u+m[0,1]*v;den=m[1,0]*u+m[1,1]*v
        if m.det()!=13:raise ArithmeticError('determinant13 required')
        for label,degree,scale in [('A',8,4),('B',12,6)]:
            original=sum(QQ(c)*n**i*den**(degree-i) for i,c in enumerate(f[label+'_coefficients_low_to_high']));new=sum(QQ(c)*u**i*v**(degree-i) for i,c in enumerate(r[label+'_coefficients_low_to_high']))
            if original!=13**scale*new:raise ArithmeticError('symbolic binary-form identity differs')
        integral=all(QQ(c).denominator()==1 for label in ('A','B') for c in r[label+'_coefficients_low_to_high'])
        if integral!=r['integral_after_curve_scale13']:raise ArithmeticError('coefficient integrality differs')
    print('INDEPENDENT SYMBOLIC84 BINARY-FORM TRANSPORTS PASS',flush=True)
if __name__=='__main__':main()
