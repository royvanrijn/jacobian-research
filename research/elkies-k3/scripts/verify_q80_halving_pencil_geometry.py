#!/usr/bin/env python3
"""Finite premises for the written ruled-quotient halving-pencil argument."""
import hashlib,json
from pathlib import Path
import sympy as S
root=Path(__file__).resolve().parents[2]
source=root/'artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisections-cheapest-1024-v1.json'
control=root/'artifacts/generated-results/elkies-k3-q80-one-cover-halving-control-v1/result.json'
b=json.loads(source.read_text())['bisections'][0]
t=S.symbols('t')
def poly(cs):return S.Poly(sum(S.Rational(c)*t**i for i,c in enumerate(cs)),t)
lift=b['lifted_section'];x=poly(lift['x1_coefficients']);y=poly(lift['y1_coefficients']);q=poly(b['branch']['numerator_coefficients'])
assert q.degree()==2 and S.gcd(q,q.diff()).degree()==0
assert S.gcd(x,y).degree()==0
# In the infinity K3 chart X=t^-4*x,Y=t^-6*y, w/t has nonzero square.
assert x.degree()==3 and y.degree()<=5 and q.LC()!=0
assert poly(lift['x0_coefficients']).degree()<=4
assert poly(lift['y0_coefficients']).degree()<=6
assert sum((a%2)<<i for i,a in enumerate(b['section_basis_w']))==65538
# F1: C^2=-1, C.F=1, F^2=0; branch R=4C+6F, K=-2C-3F.
def intersect(a,b):return -a[0]*b[0]+a[0]*b[1]+a[1]*b[0]
r=(4,6);k=(-2,-3);f=(0,1);c=(1,0)
assert intersect(r,f)==4 and intersect(r,c)==2
assert 1+(intersect(r,r)+intersect(r,k))//2==9
# K_R-2f is the restriction of 2C+F; all its three sections vanish on C.
assert (r[0]+k[0],r[1]+k[1]-2)==(2,1)
assert sum(max(1-j+1,0) for j in range(3))==3
out={'schema':'q80-halving-pencil-geometry-v1','status':'PASS','label':b['label'],
 'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
 'prior_control_sha256':hashlib.sha256(control.read_bytes()).hexdigest(),
 'finite_premises':{'x1_y1_gcd_degree':0,'x1_degree':3,'y1_degree':y.degree(),'branch_degree':2,'trace_parity':65538,'branch_curve_genus':9,'fixed_divisor_degree':2,'adjoint_space_dimension':3},
 'boundary':'Checks smoothness premises for the retained bisection image and ruled-surface intersection arithmetic. The quotient geometry, completeness of adjoint sections and uniqueness of the degree-four pencil are written arguments, not formally checked. No new branch point or new section is constructed.'}
p=root/'artifacts/generated-results/elkies-k3-q80-halving-pencil-geometry-v1/result.json'
if '--write' in __import__('sys').argv:p.write_text(json.dumps(out,indent=2)+'\n')
else:assert out==json.loads(p.read_text())
print('PASS: embedded rational bisection; F1 branch genus9; adjoint fixed divisor degree2')
