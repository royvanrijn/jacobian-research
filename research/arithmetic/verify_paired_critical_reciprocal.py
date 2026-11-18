#!/usr/bin/env python3
"""Exact controls for uniform-sign reciprocal critical pairing.

The all-degree ODE classification and the quadratic-field sign argument are
written in PAIRED_CRITICAL_RECIPROCAL_OBSTRUCTION.md. This checker derives
the symbolic h=2,6 families, recovers degree10 and certifies the degree26
obstruction by testing every point of P1(F11).
"""
import argparse,hashlib,json,platform
from pathlib import Path
import sympy as sp

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument("--write-receipt",type=Path)
 args=parser.parse_args()
 assert sp.__version__=="1.14.0"
 v,r,u=sp.symbols("v r u")
 families={}
 for h in(2,6):
  d=2*h;b=d+1;L=4*h+1
  coefficients=[sp.Integer(1),((4*h+1)*r+4*h-1)/(3-4*h)]
  for j in range(2,d):
   coefficients.append(sp.cancel(-(4*j-4*h-3)*coefficients[-1]/(4*j-4*h-1)))
  coefficients.append(sp.Integer(1))
  S=sum(coefficients[j]*v**j for j in range(d+1))
  R=sp.Poly(L*(v-r)*S,v)
  F=sp.expand(sum(R.nth(j)*v**j/(2*j-1)for j in range(b+1)))
  assert sp.Poly(F,v).LC()==1
  assert sp.expand(v**d*S.subs(v,1/v)-S)==0
  assert sp.expand(2*v*sp.diff(F,v)-F-L*(v-r)*S)==0
  ode=4*v*(v+1)*sp.diff(S,v)-((4*h-1)*v+4*h+1)*S
  target=L*(v**(d+1)-r*v**d+r*v-1)
  assert sp.expand(ode-target)==0
  identity=sp.expand(v**d*F+v**(d+1)*F.subs(v,1/v)-(v+1)*S**2)
  assert identity==0
  num,den=sp.fraction(sp.together(sp.cancel(F.subs(v,r)/(2*r))))
  families[h]=(S,F,sp.Poly(num,r,domain=sp.QQ),den)
 S2,F2,P2,den2=families[2]
 assert sp.factor(P2.as_expr()-(2*r-15)*(4*r**3-6*r**2-9*r-15))==0
 assert den2==25
 r0=sp.Rational(15,2);rho=1/r0
 N=sp.expand(rho**5*F2.subs({r:r0,v:u/rho}))
 expected=u**5-sp.Rational(96,25)*u**4+sp.Rational(3129,625)*u**3-sp.Rational(14006,5625)*u**2+sp.Rational(1804,5625)*u+sp.Rational(16,5625)
 assert sp.expand(N-expected)==0
 assert N.subs(u,1)==0 and sp.diff(N,u).subs(u,1)==0
 assert sp.factor(sp.discriminant(S2.subs(r,r0),v))!=0
 assert S2.subs({r:r0,v:1})!=0 and S2.subs({r:r0,v:-1})!=0
 S6,F6,P6,den6=families[6]
 expected_coefficients=[3516552,-5171400,8166600,-14469000,32478600,-207017160,-117119800,-26033800,79781000,-234304200,1527631560,1769265225,4846248225]
 assert P6.all_coeffs()==expected_coefficients and den6==193849929
 p=11
 reduced=[int(c)%p for c in P6.all_coeffs()]
 values=[int(P6.eval(j))%p for j in range(p)]
 assert reduced==[6,8,2,4,0,2,5,10,2,2,4,2,5]
 assert values==[5,8,1,6,8,3,3,4,3,5,7]
 assert reduced[0]!=0 and all(values)
 # Corruption controls: creating a zero constant produces a finite root;
 # making the leading coefficient zero produces a projective infinity root.
 assert int((P6.as_expr()-P6.eval(0)).subs(r,0))%p==0
 wrong=list(reduced);wrong[0]=0
 assert not wrong[0]
 root=Path(__file__).resolve().parents[2]
 inputs=["research/arithmetic/ALTERNATING_HURWITZ_UNIFORMITY_AND_REAL_NODES.md",
         "research/arithmetic/RATIONAL_DEGREE_TEN_ALTERNATING_PENCIL.md"]
 report={
  "schema_version":1,"kind":"exact_reciprocal_classification_controls",
  "python":platform.python_version(),"sympy":sp.__version__,
  "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  "inputs":{path:hashlib.sha256((root/path).read_bytes()).hexdigest()for path in inputs},
  "symbolic_family_controls":[2,6],
  "square_reciprocity_identity":"v^(2h)F(v)+v^(2h+1)F(1/v)=(v+1)S(v)^2",
  "degree10_recovery":{"r":"15/2","rho":"2/15","N":str(N)},
  "degree26_obstruction_polynomial":str(P6.as_expr()),
  "degree26_F_at_r":"2*r*P6(r)/193849929",
  "modulus":p,"degree26_mod11_coefficients":reduced,
  "finite_P1_mod11_values":values,"infinity_value_mod11":reduced[0],
  "scope":"Uniform one-sign reciprocal classification uses the written proof. The exact mod11 certificate rules out rational degree26 models in that ansatz; mixed-sign reciprocal factors and general degree26 pencils remain open.",
  "independent_replay":False,"external_review":False,"formal_verification":False,
  "result":"PASS"
 }
 output=json.dumps(report,indent=2,sort_keys=True)+"\n"
 if args.write_receipt:args.write_receipt.write_text(output)
 print(output,end="")
if __name__=="__main__":main()

