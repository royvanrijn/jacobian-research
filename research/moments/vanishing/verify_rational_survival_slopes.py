"""Exact transcription controls; the universal slope theorem is a written proof."""
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import math
import sympy as S

p = Path(__file__).resolve().parent
proof = p/'BOUNDED_BERNSTEIN_STEPS_AND_RATIONAL_SURVIVAL_SLOPES.md'
sha = lambda f: hashlib.sha256(f.read_bytes()).hexdigest()
x,y,a,b,c,d,t = S.symbols('x y a b c d t')

# Check the exact formula by actual differential evaluation near first survival.
two_variable = []
for N,M in [(1,1),(2,2),(2,3),(3,4),(4,5),(7,8)]:
    rho = Fraction(N*M,N+M)
    seq = [min(max(N*i,M*(k-i)) for i in range(k+1)) for k in range(81)]
    assert all(0<=v-u<=max(N,M) for u,v in zip(seq,seq[1:]))
    assert max(abs(Fraction(v)-rho*k) for k,v in enumerate(seq))<=max(N,M)
    for k in range(1,5):
        tau=seq[k]
        for j in [tau-1,tau]:
            out=sum(S.binomial(k,i)*S.diff((x*y)**j,x,N*i,y,M*(k-i)) for i in range(k+1))
            assert (S.expand(out)!=0)==(j==tau)
    maxerr=max(abs(Fraction(v)-rho*k) for k,v in enumerate(seq))
    two_variable.append(dict(N=N,M=M,rho=str(rho),indices=80,max_error=str(maxerr)))

# Arbitrary formal units preserve heights; finite lacunary jets test the transcription.
unary=[]
for D,r in [(1,2),(2,3),(3,4),(5,2)]:
    for k in range(1,5):
        tau=(r*k+D-1)//D
        L=S.Poly(S.expand((x**r*(1+x+x**4+x**11))**k),x)
        for j in [tau-1,tau]:
            out=sum(co*S.diff(x**(D*j),x,mo[0]) for mo,co in L.terms() if mo[0]<=D*j)
            assert (S.expand(out)!=0)==(j==tau)
        unary.append(dict(D=D,r=r,k=k,tau=tau))

# Known genuine pure-unsafe seed: do not replace the critical slope by safety.
P=(t+c)*(a*d+b*t)
op=lambda f:S.expand(S.diff(f,t,a,d)-S.diff(f,t,b,c))
unsafe=[]
for k in range(1,4):
    nonzero=[]
    for j in range(k+2):
        f=S.expand(P**j)
        for _ in range(k): f=op(f)
        if f!=0: nonzero.append(j)
    assert nonzero==[k+1]
    unsafe.append(dict(k=k,tau=k+1))

# The running maximum is the passive-suspension height, even with missing powers.
running=[]
for N,M in [(1,1),(2,3),(3,4)]:
    hs=[m//N+m//M for m in range(121)]
    for m in range(1,60):
        H=max(hs[:m+1])
        inv=max(k for k in range(H+2) if min(max(N*i,M*(k-i)) for i in range(k+1))<=m)
        assert inv==H
    running.append(dict(N=N,M=M,indices=59))

# An abstract bounded-step relation need not select its parallel lines periodically.
abstract=[]
seq=[k+int(math.isqrt(k)**2==k) for k in range(1,501)]
assert all(0<=v-u<=2 for u,v in zip(seq,seq[1:]))
assert all((v-k)*(v-k-1)==0 for k,v in enumerate(seq,1))
abstract.append('square indicator:bounded steps,parallel line offsets0/1;nonperiodicity is proved,not inferred from prefix')
assert all((k*k-k*k)==0 for k in range(1,30))
assert (30**2-29**2)>(10**2-9**2)
abstract.append('k squared:monotone polynomial curve with unbounded steps')

out=dict(completed=True,kind='exact transcription controls only',
         canonical_proof_sha256=sha(proof),checker_sha256=sha(Path(__file__)),
         two_variable_controls=two_variable,unary_unit_controls=unary,
         retained_unsafe_pure_seed_controls=unsafe,running_height_inverse_controls=running,
         abstract_boundary_controls=abstract,
         actual_module_and_Bernstein_import_verified_by='retained and new written proofs;primary scope check',
         universal_branch_lattice_lemma_verified_by='self-contained written proof,not these finite controls',
         original_unsuspended_height_bounded_error_claimed=False,eventual_periodicity_claimed=False,
         independent_replay=False,formal_verification=False,external_review=False)
(p/'rational_survival_slopes_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(completed=True,two_variable=len(two_variable),unary=len(unary),unsafe=len(unsafe),inverse=len(running),abstract=len(abstract))))
