#!/usr/bin/env python3
"""Finite exact frame controls for QEP1; normal-form exhaustion is written."""
from pathlib import Path
import hashlib,json,tempfile
import sympy as sp
HERE=Path(__file__).resolve().parent
x,q,h=sp.symbols('x q h')
r=x*(2-3*x*q)
b=(3*x*x,2*(1-3*x*q))
d=((1+3*x*q)/2,-3*q*q)
def act(v,f):return sp.expand(v[0]*sp.diff(f,x)+v[1]*sp.diff(f,q))
assert act(b,r)==0
assert act(d,r)==1
assert act((b[0]+1,b[1]),r)!=0  # A corrupted tangential coefficient is rejected.
assert sp.det(sp.Matrix([b,d]))==-1
assert all(sp.expand(act(d,b[i])-act(b,d[i])-sp.Rational(9,2)*q*b[i])==0 for i in range(2))
assert tuple(sp.expand(2*(1-3*x*q)*d[i]+3*q*q*b[i]) for i in range(2))==(1,0)
assert tuple(sp.expand(-3*x*x*d[i]+(1+3*x*q)*b[i]/2) for i in range(2))==(0,1)
# Two disjoint coefficient-fiber components, stable under the relative derivation.
factor=2-3*x*q
assert act(b,x)==3*x*x and sp.expand(act(b,factor)+3*r)==0
id1=factor/2;id2=3*x*q/2
assert sp.expand(id1+id2)==1
assert sp.expand(id1*id2-sp.Rational(3,4)*q*r)==0
assert sp.expand(id1*id1-id1+sp.Rational(3,4)*q*r)==0
assert b[1].subs(x,0)==2
assert act(b,1/x)==-3
# The generic coefficient chart is exact, not merely a fraction-field inclusion.
v,z=sp.symbols('v z')
Q=(2*v-z*v*v)/3
assert sp.expand(r.subs({x:1/v,q:Q})-z)==0
assert sp.cancel(1/x-(2-3*x*q)/r)==0
names=['EXACT_POSITION_CENTRALIZER_AND_COMPLETION.md','check_exact_position_centralizer.py','CENTRALIZER_NORMAL_FORM.md','QUINTIC_PHASE_AND_BRANCH_OBSTRUCTION.md']
receipt={'claim':'QEP1-global-Ore','status':'PASS','checked':['unimodular differential frame','both exact inverse vector-field formulas','R and frame commutators','Lie bracket [D0,B]=(9h/2)QB','Chinese-remainder fiber idempotents and preserved ideals','canonical first-component pair and Laurent second-component pair','exact generic coefficient chart'],'hashes':{n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},'boundary':'Unique finite PBW normal form, all-operator centralizer exhaustion and component algebra identifications use the written symbol and Ore arguments; finite identities do not establish the missing Long mates.','independent_replay':False,'formal_verification':False,'external_review':False}
out=Path(tempfile.mkdtemp(prefix='qep1-global-'))/'exact_position_controls.json';out.write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'claim':receipt['claim'],'status':'PASS','receipt':str(out)}))
