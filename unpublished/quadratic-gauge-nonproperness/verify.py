"""Finite exact controls for the nonproperness formulas; not a general proof."""
import json
from pathlib import Path
import sympy as s

x,y,z,t,q,pi,B,C,S=s.symbols("x y z t q pi B C S")
g1,g2,g3=s.symbols("g1 g2 g3", nonzero=True)
t0=1+x*y;q0=t0**2*z+g1/g3*y**2*(1+3*t0)
Pi=t0*q0
B0=y+3*g3/g1*x*q0+2*g2/g1*t0*q0
C0=x*(5-3*t0)-g3/g1*x**3*z
assert s.factor(s.Matrix([Pi,B0,C0]).jacobian((x,y,z)).det())==-2
assert s.simplify(B0.subs(z,-g1/g3*y**2*(1+3*t0)/t0**2)-y)==0
cq=s.factor(C0.subs(z,-g1/g3*y**2*(1+3*t0)/t0**2))
assert s.factor(cq-(2*x/t0-y*(x/t0)**2))==0
assert s.factor(1-y*x/t0-1/t0)==0
assert s.simplify(q0.subs(y,-1/x)-g1/g3/x**2)==0
assert s.simplify(B0.subs(y,-1/x)-2/x)==0
assert s.diff(C0.subs(y,-1/x),z)==-g3/g1*x**3
records=[]
for n in range(3,8):
    Gpi=S+pi*(S**2+S**3)+sum(pi**k*S**k for k in range(4,n+1))
    E=Gpi-(B*S**2+C)/2
    if n==3:
        D=s.Poly(s.discriminant(E,S),pi)
        assert s.factor(D.as_expr().subs(pi,0)-B**2*(1-B*C)/4)==0
        records.append({"degree":n,"exact_global_order":0,"symbolic_boundary_coefficient":True})
    else:
        targets=[(2,3),(3,2),(1,-1)];checks=[]
        if n==4:
            D=s.Poly(s.discriminant(E,S),pi)
            order=min((k[0] for k,v in D.terms() if v!=0))
            assert order==n*n-3*n-2
            lead=s.factor(D.coeff_monomial(pi**order))
            assert s.factor(lead**2-(B**2*(1-B*C)/4)**2)==0
        for b,c in targets:
            D=s.Poly(s.discriminant(E.subs({B:b,C:c}),S),pi)
            order=min(k[0] for k,v in D.terms() if v!=0)
            assert order==n*n-3*n-2
            r=n-3;expected=s.Rational(r**r,4)*b*b*(1-b*c)
            lead=D.coeff_monomial(pi**order)
            assert lead**2==expected**2
            checks.append({"B":b,"C":c,"order":order})
        records.append({"degree":n,"targets":checks,"symbolic_boundary_coefficient":n==4})
report={"scope":"Generic cubic boundary identities and selected degree-three-through-seven discriminant controls; not the all-degree Newton proof, Jelonek criterion, or Lean verification.","checks":records}
Path(__file__).with_name("symbolic_checks.json").write_text(json.dumps(report,indent=2)+"\n")
print(json.dumps(report,indent=2))
