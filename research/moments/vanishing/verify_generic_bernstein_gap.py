"""Exact controls for the proposed all-degree proof; not a theorem prover."""
from pathlib import Path
import hashlib, json
import sympy as S

root=Path(__file__).resolve().parents[3]
draft=Path(__file__).with_name('GENERIC_BERNSTEIN_EQUATIONS_AND_FORMAL_GAP_CLOSURE.md')
x,y,e,s,z=S.symbols('x y e s z')

def action(poly, phase, symbol):
    out=0
    for (a,b),c in S.Poly(symbol,x,y).terms():
        out += c*S.diff(poly,x,a,y,b)
    return S.expand(out)

# Two genuinely different spatial directions, including a mixed cubic phase.
P=x*x*y+y*y+x
operators=[lambda f:S.diff(f,x)+e*S.diff(P,x)*f,
           lambda f:S.diff(f,y)+e*S.diff(P,y)*f,
           lambda f:-x*f,lambda f:-y*f,
           lambda f:S.diff(f,e)+P*f]
f=x*x*y+e*x*y*y+e*e*y
weyl=[]
for i,j,want in [(2,0,f),(3,1,f),(0,1,0),(4,0,0),(4,1,0),
                 (4,2,0),(4,3,0),(2,1,0),(3,0,0)]:
    assert S.expand(operators[i](operators[j](f))-
                    operators[j](operators[i](f))-want)==0
    weyl.append([i,j])

# Coefficientwise conjugation: symbol has both low and high participating jets.
lam=x*x+x*y*y+y**5
time_order=5
expP=sum(e**j*S.expand(P**j)/S.factorial(j) for j in range(time_order+1))
v=expP
conjugations=[]
for k in range(1,4):
    v=action(v,P,lam)
    out=0
    for j in range(time_order+1):
        for q in range(j+1):
            out+=e**j*(-P)**q/S.factorial(q)*S.expand(v).coeff(e,j-q)
    def xi_apply(g,i):
        return S.expand(S.diff(g,[x,y][i])+e*S.diff(P,[x,y][i])*g)
    w=S.Integer(1)
    for _ in range(k):
        accum=0
        for (a,b),c in S.Poly(lam,x,y).terms():
            term=w
            for __ in range(a):term=xi_apply(term,0)
            for __ in range(b):term=xi_apply(term,1)
            accum+=c*term
        w=S.expand(accum)
        w=sum(w.coeff(e,j)*e**j for j in range(time_order+1))
    assert all(S.expand(w-out).coeff(e,j)==0 for j in range(time_order+1))
    conjugations.append(k)

# Arbitrary-unary-unit all-index formula: highest polynomial degree cannot cancel.
# Here use a lacunary-jet unit, with only the degree-relevant finite jet applied.
unary=[]
for D,r in [(1,2),(2,3),(3,4),(4,7)]:
    for k in range(1,5):
        tau=(r*k+D-1)//D
        jet=S.Poly(S.expand((x**r*(1+x+x**4+x**11))**k),x)
        for m in [max(0,tau-1),tau]:
            val=sum(c*S.diff(x**(D*m),x,a[0]) for a,c in jet.terms()
                    if a[0]<=D*m)
            assert (S.expand(val)!=0)==(m>=tau)
        assert S.prod(D*tau-r*k-j for j in range(D))==0
        unary.append(dict(D=D,r=r,k=k,tau=tau))

# Finite negative-time specialization: the correct exponent is -tau-1.
# P=x^D, lambda=xi^r gives leading nonzero coefficient known exactly.
pole=[]
for D,r,k in [(2,3,2),(3,4,2),(3,4,3)]:
    tau=(r*k+D-1)//D
    vk=sum(S.diff(x**(D*j),x,r*k)*e**j/S.factorial(j)
           for j in range(tau,tau+3))
    unit=S.expand(vk/e**tau)
    assert S.expand(unit).coeff(e,0)!=0
    assert S.expand(S.diff(unit,e)).as_leading_term(e).as_powers_dict().get(e,0)>=0
    assert S.expand(vk/e**(tau+1)).as_leading_term(e).as_powers_dict().get(e,0)==-1
    pole.append([D,r,k,tau])

# Exact abstract controls: monotonicity, eventuality, all-index relation are separate.
counter=[]
for n in range(3,9):
    before=2*(n*n-1); after=n*n+n
    assert before>after
    for k,tau in [(n*n,after),(n*n-1,before)]:
        assert ((tau-k)**2-k)*(tau-2*k)==0
    counter.append(n)
assert ((2+2-2)**2-2)!=0  # ceil(sqrt(2)) is not an algebraic equality root.
q=S.Poly((s+z-s-1),s)  # bounded excess polynomial, zero only at fixed z=1.
assert q.as_expr()==z-1

receipt=dict(completed=True,kind='exact illustrative controls, not formal proof',
             canonical_proof_sha256=hashlib.sha256(draft.read_bytes()).hexdigest(),
             checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             weyl_relations=weyl,mixed_cubic_conjugation_k=conjugations,
             unary_first_survival_checks=unary,negative_time_pole_checks=pole,
             nonmonotone_polynomial_constraint_square_controls=counter,
             generic_Bernstein_import_computationally_verified=False,
             universal_integer_sequence_lemma_verified_by='written proof only',
             formal_verification=False,independent_replay=False,external_review=False)
out=Path(__file__).with_name('generic_bernstein_gap_verification.json')
out.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(dict(completed=True,weyl=len(weyl),conjugations=len(conjugations),
                      unary=len(unary),poles=len(pole),squares=len(counter))))
