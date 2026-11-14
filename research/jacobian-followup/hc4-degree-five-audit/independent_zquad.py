"""One complete unary-top / quadratic-passive-block gate; no source verifier import.

Generate every permitted coefficient independently with FLINT, and saturate
by det(A1) and by the nonzero constant full determinant. Bounded externally.
"""
from pathlib import Path
import itertools
import json
import sys
import time

sys.path.insert(0,'/tmp/jacobian-hc4-sage-deps')
from sage.all import PolynomialRing, QQ, prod
import flint

HERE=Path(__file__).resolve().parent
start=time.monotonic()
monomials=[]
# General base terms of degrees 2..4; the x^4 term has already been removed
# by translation of x in a unary-top potential.
for degree in range(2,5):
    for i in range(degree+1):
        if (i,degree-i)==(4,0):continue
        monomials.append((i,degree-i,0,0))
# Terms linear in either passive variable, with base degree 1..3.
for passive in (2,3):
    for degree in range(1,4):
        for i in range(degree+1):
            e=[i,degree-i,0,0];e[passive]=1;monomials.append(tuple(e))
# Complete constant and x-linear symmetric passive Hessian blocks.
for xpower in (0,1):
    for j in range(3):monomials.append((xpower,0,j,2-j))
assert len(set(monomials))==len(monomials)==35
names=['x','y','z','w']+[f'c{i}' for i in range(len(monomials))]
ctx=flint.fmpq_mpoly_ctx.get(names,'lex');g=ctx.gens();X=g[:4]
F=X[0]**5
for coefficient,exponent in zip(g[4:],monomials):
    term=coefficient
    for variable,power in zip(X,exponent):term*=variable**power
    F+=term
H=[[F.derivative(i).derivative(j) for j in range(4)] for i in range(4)]
D=ctx.from_dict({})
for permutation in itertools.permutations(range(4)):
    sign=(-1)**sum(permutation[i]>permutation[j] for i in range(4) for j in range(i+1,4))
    term=ctx.constant(sign)
    for i,j in enumerate(permutation):term*=H[i][j]
    D+=term
R=PolynomialRing(QQ,names[4:]+['T','U'],order='degrevlex');rg=R.gens()
bucket={}
for exponent,coefficient in D.to_dict().items():
    term=QQ(str(coefficient))*prod(rg[i]**int(power) for i,power in enumerate(exponent[4:]) if power)
    bucket[tuple(exponent[:4])]=bucket.get(tuple(exponent[:4]),R(0))+term
constant=bucket.pop((0,0,0,0))
coefficient_map={exponent:rg[i] for i,exponent in enumerate(monomials)}
detA1=4*coefficient_map[1,0,2,0]*coefficient_map[1,0,0,2]-coefficient_map[1,0,1,1]**2
equations=list(bucket.values())+[rg[-2]*detA1-1,rg[-1]*constant-1]
payload=dict(monomials=monomials,variables=list(R.variable_names()),equations=list(map(str,equations)),constant_determinant=str(constant),detA1=str(detA1),coverage='All degree-five unary-top potentials quadratic in (z,w), with f_zz=A0+x*A1 and no pure x^4 term; no normalized rank-two representative used')
(HERE/'zquad-input.json').write_text(json.dumps(payload,indent=2)+'\n')
print('INDEPENDENT_INPUT',len(monomials),'coefficients',len(equations),'equations',flush=True)
G=R.ideal(equations).groebner_basis(algorithm='libsingular:slimgb')
result=dict(unit_ideal=list(G)==[1],basis_size=len(G),elapsed_seconds=round(time.monotonic()-start,3),basis=list(map(str,G)))
(HERE/'zquad-result.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result),flush=True)
