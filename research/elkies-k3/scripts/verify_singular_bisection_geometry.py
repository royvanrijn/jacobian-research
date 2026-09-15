#!/usr/bin/env python3
"""Check symbolic intersection premises; geometric/finiteness proofs are written."""
import json
from pathlib import Path
import sympy as S
root=Path(__file__).resolve().parents[2]
a,z=S.symbols('a z',integer=True)
# NS=U+(-MW), B=2O+(z+4)F+v_T, O^2=-2.
hT=4*z+10-2*a
assert S.expand(-8+4*(z+4)-hT-(2*a-2))==0
hQ=8+2*z
assert S.expand(4*hQ-2*hT-4*(a+3))==0
# Primitive-trace quotient has irreducible branch -2K on F_n;
# its intersection with C0 is 4-2n, so n<=2.
rows=[]
for n in range(3):
 # Section C0+kF must have square1 for arithmetic genus2 pullback.
 k=S.Rational(n+1,2)
 rows.append({'n':n,'required_k':str(k),'integral':bool(k.is_Integer)})
assert [r['n'] for r in rows if r['integral']]==[1]
def product(x,y,n):return -n*x[0]*y[0]+x[0]*y[1]+x[1]*y[0]
r=(4,6);K=(-2,-3);line=(1,1)
assert product(line,line,1)==1
assert product(r,line,1)==6
assert 1+(product(r,r,1)+product(r,K,1))//2==9
# Local delta of y^2=u^m at a branch tangency.
contacts=[{'contact_order':m,'delta':m//2} for m in range(1,7)]
assert [v['contact_order'] for v in contacts if v['delta']==1]==[2,3]
out={'schema':'singular-bisection-geometry-v1','status':'PASS',
 'anti_trace_height':'4*(arithmetic_genus+3)',
 'quotient_candidates':rows,'plane_branch_degree':6,'branch_normalization_genus':9,
 'local_contact_delta':contacts,
 'scope':'Symbolic height/intersection identities and finite ruled-surface cases only. Quotient extension, monodromy, tangent correspondence and Faltings finiteness are written arguments; no arithmetic point list or positive common cover is supplied.'}
p=root/'artifacts/generated-results/elkies-k3-singular-bisection-geometry-v1/result.json'
if '--write' in __import__('sys').argv:p.parent.mkdir(exist_ok=True);p.write_text(json.dumps(out,indent=2)+'\n')
else:assert out==json.loads(p.read_text())
print('PASS: anti-trace height identity; unique F1 case; contact orders2/3')
