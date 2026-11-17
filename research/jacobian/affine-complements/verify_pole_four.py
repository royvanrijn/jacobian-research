"""Exact controls for PF2RULED3's necessary normal form and one excluded lift family."""
from pathlib import Path
from functools import reduce
from itertools import combinations
from time import monotonic
import hashlib,json,platform
import sympy as S

HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def zero(name,expression):
    check(name,S.cancel(expression)==0)
def wedge(F,G,u,v):
    return S.cancel(S.diff(F,u)*S.diff(G,v)-S.diff(F,v)*S.diff(G,u))

start=monotonic()
T,A,H,C,k,U,D=S.symbols("T A H C k U D")
P=T**3+T*A+H
Q=T**4+S.Rational(4,3)*T*T*A+S.Rational(4,3)*T*H+S.Rational(2,9)*A*A+k*T*T+C
expected=[
 ((T,A),S.Rational(4,9)*A*A-S.Rational(4,3)*T*H-2*k*T*T),
 ((T,H),-S.Rational(4,3)*T*A-S.Rational(4,3)*H-2*k*T),
 ((T,C),3*T*T+A),((A,H),-S.Rational(4,9)*A),((A,C),T),((H,C),1)]
for pair,value in expected:
    zero("normal-form exterior coefficient "+str(pair),wedge(P,Q,*pair)-value)
coreP=T**3+U
coreQ=T**4+S.Rational(4,3)*T*U+D
zero("first residual normalization core",wedge(coreP,coreQ,T,U)+S.Rational(4,3)*U)
zero("first residual normalization dT dD",wedge(coreP,coreQ,T,D)-3*T*T)

w=S.symbols("w")
tau,a2,d2=[S.Function(name)(w) for name in ("tau","a2","d2")]
pole6=(-S.Rational(8,3)*S.diff(tau,w)*a2*a2
       +S.Rational(4,3)*tau*a2*S.diff(a2,w)
       +6*tau*tau*S.diff(tau,w)*d2-3*tau**3*S.diff(d2,w))
zero("pole-six integrating identity",
     pole6+3*tau**5*S.diff(d2/tau**2-S.Rational(2,9)*a2*a2/tau**4,w))

alpha,eta,gamma=[S.Function(name)(w) for name in ("alpha","eta","gamma")]
def W(f,g):
    return S.diff(f,w)*g-f*S.diff(g,w)
pole5=((S.Rational(4,9)*alpha**2-S.Rational(4,3)*tau*eta-2*k*tau**2)*W(tau,alpha)
       -S.Rational(4,3)*tau*alpha*W(tau,eta)+3*tau**2*W(tau,gamma))
primitive=3*gamma/tau-S.Rational(4,3)*alpha*eta/tau**2+S.Rational(4,27)*alpha**3/tau**3-2*k*alpha/tau
zero("pole-five cubic integrating identity",pole5+tau**4*S.diff(primitive,w))

s,t=S.symbols("s t")
symbols=[s**3,3*s*s*t,3*s*t*t,S.Rational(8,3)*t**3]
tr,ar,hr,cr=symbols
zero("base-point-free cubic leading relation",4*ar**3-36*tr*ar*hr+81*tr*tr*cr)
check("cubic leading data have no common zero",
      reduce(S.gcd,[S.Poly(f,s,t,domain=S.QQ) for f in symbols]).total_degree()==0)

z,w=S.symbols("z w")
t0,a0,b0,c0,ell=S.symbols("t0 a0 b0 c0 ell")
g1=w/z
g2=w*w/z
g3=(z+w*w)*w/z
g4=(z+w*w)**2/z
functions={T:g4+t0,A:3*g3+a0,H:3*g2+b0}
functions[C]=S.Rational(8,3)*g1+S.Rational(2,3)*k*functions[A]+ell*functions[T]+c0
# Reconstruct through exterior coefficients, independently of the producer's
# direct differentiation of the large P,Q expressions.
volume=0
for pair,value in expected:
    u,v=pair
    volume+=S.sympify(value).subs(functions)*wedge(functions[u],functions[v],z,w)
volume=S.cancel(volume)
coefficient=S.cancel((z**4*volume)).subs(z,0)
formula=w**7*((-4*b0-3*ell+8)*w**3-6*k*w*w+12*t0*w+4)
zero("independent exterior reconstruction of all lift parameters",coefficient-formula)
check("unavoidable pole-four coefficient is four",
      S.Poly(coefficient,w).coeff_monomial(w**7)==4)
payload=json.loads((HERE/"pole4_cubic_control.json").read_text())
check("control producer hash",
      payload["source_sha256"]==hashlib.sha256((HERE/"explore_pole4_cubic_control.py").read_bytes()).hexdigest())
zero("saved control coefficient matches independent reconstruction",
     S.sympify(payload["pole_four_coefficient"])-formula)
equations=list(map(S.sympify,payload["equations"]))
constant_index=next(i for i,v in enumerate(equations) if v!=0 and not v.free_symbols)
multiplier=1/equations[constant_index]
zero("exact unit certificate for the full six-parameter control family",
     multiplier*equations[constant_index]-1)
report=dict(as_of="2026-11-14",outcome="PASS",check_count=len(checks),checks=checks,
            written_result="Every reduced (3,4) candidate has a four-V1-section normal form; the specified six-parameter cubic-symbol lift family is excluded",
            remaining_scalar_parameter_count=21,
            scope="PF2RULED3 normal form and one control; complete pole-four exclusion is PF2RULED4",
            excluded_control_unit_certificate=dict(equation_index=constant_index,multiplier=str(multiplier)),
            assurance=dict(full_pole_four_layer_excluded=False,finite_map_constructed=False,
                           actual_ruled_compatibility_proved=False,formal_verification=False,
                           independent_full_reproof=False),
            source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                           for p in [Path(__file__),HERE/"POLE_FOUR_FRONTIER.md",
                                     HERE/"explore_pole4_cubic_control.py"]},
            python=platform.python_version(),sympy=S.__version__,
            elapsed_seconds=round(monotonic()-start,6))
(HERE/"pole_four_receipt.json").write_text(json.dumps(report,indent=2)+"\n")
print(f"PASS {len(checks)} exact controls; 21-parameter normal form and control verified; complete exclusion is PF2RULED4")

