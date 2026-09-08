#!/usr/bin/env sage-python
"""Exact coefficient support for complete residual binary-form collapse."""
import json,hashlib,sys
from pathlib import Path
from sage.all import ZZ,PolynomialRing,gcd
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
from research_runtime.store import checkpoint
D=ROOT/'artifacts/local/elliptic-curves/det1092-parameter-hessian-support-v1';SOURCE=ROOT/'artifacts/local/elliptic-curves/det1092-integral-parameter-model-v1/result.json'
def hashed(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=json.loads((D/'protocol.json').read_bytes());assert p['source_sha256']==hashed(SOURCE) and p['script_sha256']==hashed(Path(__file__))
 s=json.loads(SOURCE.read_bytes());assert s['status']=='PASS';R=PolynomialRing(ZZ,'t');F=R(s['A_coefficients']);G=R(s['B_coefficients'])
 H4=8*7*F*F.derivative(2)-7**2*F.derivative()**2;H6=12*11*G*G.derivative(2)-11**2*G.derivative()**2;J=12*F.derivative()*G-8*F*G.derivative()
 contents=[abs(gcd(list(h))) for h in [H4,H6,J]];N=gcd(contents);assert N>0
 data=dict(status='FACTORING',source_sha256=hashed(SOURCE),hessian_contents=list(map(str,contents)),joint_content=str(N));checkpoint(D/'result.json',data);fac=list(N.factor(proof=True));assert all(q.is_prime(proof=True) for q,e in fac)
 data.update(status='PASS',prime_factors=[[str(q),int(e)] for q,e in fac],candidate_primes=sorted({2,3,5,7,11,13}|{int(q) for q,e in fac}),scope='If the two binary forms reduce to powers of the same linear form, both Hessians and their cross-Jacobian vanish. Their coefficient gcd gives necessary prime support; no converse, local/global minimality or complete partial-collapse support is claimed.');checkpoint(D/'result.json',data);print('JOINT CONTENT bits',N.nbits(),'FACTORS',data['prime_factors'],flush=True)
if __name__=='__main__':main()
