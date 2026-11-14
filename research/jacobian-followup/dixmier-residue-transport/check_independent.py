#!/usr/bin/env python3
"""Separate SymPy replay of the scalar period and native-support failure.

Does not import the transport implementation or the legacy sparse engine.
This checks core scalar identities, not all 945 producer columns.
"""
from pathlib import Path
from math import comb
import hashlib
import json
import sympy as s

HERE=Path(__file__).resolve().parent
X,Q,Z,L=s.symbols("X Q Z lambda")
vs=(X,Q,Z,L)
packet=json.loads((HERE/"inputs.json").read_text())
receipt=json.loads((HERE/"receipt.json").read_text())
assert hashlib.sha256((HERE/"inputs.json").read_bytes()).hexdigest()==receipt["inputs_sha256"]
def read(rows):
    return s.Add(*(s.Rational(c)*s.prod(v**e for v,e in zip(vs,m)) for m,c in rows))
S,T,u,v=(read(packet[k]) for k in ("S","T","s2","t2"))
lam={tuple(m[:3]):s.Rational(c) for m,c in packet["lambda"]}
def delta(p): return 3*X**2*s.diff(p,X)+(2-6*X*Q)*s.diff(p,Q)
def iterate(p,n):
    for _ in range(n): p=s.expand(delta(p))
    return p
def pi(p,q,n=1):
    return s.expand(sum((-1)**j*comb(n,j)*iterate(s.diff(p,Z,n-j),j)*s.diff(iterate(q,n-j),Z,j) for j in range(n+1)))
def evaluate(p):
    p=s.Poly(s.expand(p),X,Q,Z)
    return sum(c*p.coeff_monomial(X**i*Q**j*Z**k) for (i,j,k),c in lam.items())
assert pi(S,T)==1
assert s.expand(pi(u,T)+pi(S,v)+pi(S,T,3)/24)==0
O=pi(u,v)+pi(u,T,3)/24+pi(S,v,3)/24+pi(S,T,5)/1920
value=evaluate(O)
assert value==-s.Rational(47547660815739,190658)
assert evaluate(pi(X**12*Q**6*Z,T))==0
assert evaluate(pi(X**16*Q**6,T))==432
assert s.expand(evaluate(pi(X**12*Q**6*(Z-3*L*X**4),T)))==-1296*L
H=L*X**3/3
assert s.expand(iterate(H,3)/24)==s.Rational(45,2)*L*X**6
assert s.expand(iterate(H,5)/1920)==s.Rational(1701,16)*L*X**8
# Independently check every listed dual coefficient against inverse
# substitution on its input monomial, with no sparse transform import.
dual={tuple(m):s.Rational(c) for m,c in packet["transported_lambda"]}
triples={m[:3] for m in dual}
for i,j,k in triples:
    expected=evaluate(X**i*Q**j*(Z-3*L*X**4)**k)
    actual=sum(c*L**m[3] for m,c in dual.items() if m[:3]==(i,j,k))
    assert s.expand(expected-actual)==0
# All six ambient generator relations for the Weyl shear. Since these
# images are affine in momenta, their Moyal commutators equal h*Poisson.
x,q,p,z=s.symbols("x q p z")
images=(x,q,p+L*x**2,z)
src=(x,q,p,z)
def pb(a,b): return s.expand(s.diff(a,p)*s.diff(b,x)-s.diff(a,x)*s.diff(b,p)+s.diff(a,z)*s.diff(b,q)-s.diff(a,q)*s.diff(b,z))
for i in range(4):
    for j in range(i+1,4): assert pb(images[i],images[j])==pb(src[i],src[j])
out={"status":"PASS","scalar_period":str(value),"native_failure":"-1296*lambda", "original_extra_direction":432,
     "symbolic_dual_coefficients_checked":len(dual),"ambient_generator_relations":6,
     "inputs_sha256":receipt["inputs_sha256"],"checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     "sympy":s.__version__,"scope":"Independent scalar identities, dual coefficient values and six ambient shear relations. No independent all-degree proof or full 945-column replay."}
(HERE/"independent.json").write_text(json.dumps(out,indent=2)+"\n")
print("PASS: independent seed, lower equation, scalar residue, native failure, quantum coefficients, dual table and six ambient shear relations")
