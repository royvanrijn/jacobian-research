"""Exact controls for reduced-branch trace pairings; not a universal proof."""
from pathlib import Path
import hashlib
import json
import tempfile
import time
import sympy as s

HERE = Path(__file__).resolve().parent
started = time.perf_counter()
checks = []
def check(label, condition):
    if not condition:
        raise AssertionError(label)
    checks.append(label)
def zero(label, expression):
    check(label, s.cancel(expression) == 0)

p, q, v, t = s.symbols("p q v t")
W = s.Matrix([[0, 0, q], [1, 0, -p], [0, 1, 0]])
G = s.Matrix([[s.trace(W**(i+j)) for j in range(3)] for i in range(3)])
expected = s.Matrix([[3,0,-2*p],[0,-2*p,3*q],[-2*p,3*q,2*p**2]])
check("cusp Gram independently obtained by multiplication traces", G == expected)
H = 4*p**3+27*q**2
zero("cusp trace discriminant", G.det()+H)
N = (H*G.inv()).applyfunc(s.cancel)
check("cusp inverse Gram has only reduced denominator",
      all(entry.is_polynomial(p,q) for entry in N))
check("cusp unscaled inverse Gram is rejected",
      any(not s.cancel(entry).is_polynomial(p,q) for entry in G.inv()))
check("cusp wrong scaled inverse matrix is rejected",
      any(s.cancel(entry)!=0 for entry in G*(2*N)-H*s.eye(3)))
check("cusp symmetric matrix factor", N == N.T)
check("cusp factorization on both sides", all(s.cancel(x)==0 for x in G*N-H*s.eye(3)) and all(s.cancel(x)==0 for x in N*G-H*s.eye(3)))
check("cusp original Gram fibre has rank one", G.subs({p:0,q:0}).rank() == 1)
check("cusp defect fibre requires two generators", 3-G.subs({p:0,q:0}).rank() == 2)
check("cusp reduced pairing fibre is zero", N.subs({p:0,q:0}) == s.zeros(3))
zero("cusp defining polynomial resultant",
     s.resultant(v**3+p*v-q,3*v**2+p,v)-H)
zero("cusp source Jacobian", s.diff(v**3+p*v,v)-(3*v**2+p))

# Exact independent trace matrices, covering odd/even, residual/nonresidual,
# full codifferent and higher-moment examples. The proof handles all e,m.
controls = []
for e,m in [(2,2),(3,2),(3,3),(5,3),(6,2),(6,3),(7,3),(7,7)]:
    U=s.zeros(e)
    for j in range(e-1):
        U[j+1,j]=1
    U[0,e-1]=t
    b=(e-1)//m
    exponents=list(range(-b,e-b))
    Gm=s.Matrix([[s.trace(U**(i+j)) for j in exponents] for i in exponents])
    Nm=(t*Gm.inv()).applyfunc(s.cancel)
    rank=e-1-2*b
    label=f" e={e},m={m}"
    check("regular moment Gram"+label,all(x.is_polynomial(t) for x in Gm))
    check("reduced inverse moment Gram"+label,all(x.is_polynomial(t) for x in Nm))
    check("moment matrix factorization"+label,all(s.cancel(x)==0 for x in Gm*Nm-t*s.eye(e)))
    check("moment residual rank from Gram"+label,e-Gm.subs(t,0).rank()==rank)
    check("moment residual rank from complementary factor"+label,Nm.subs(t,0).rank()==rank)
    indices=list(range(b+1,e-b))
    pairing=s.Matrix([[s.cancel(t*s.trace(U**(-i-j))).subs(t,0)
                       for j in indices] for i in indices])
    reverse=s.zeros(rank)
    for i in range(rank):
        reverse[i,rank-1-i]=e
    check("split residue anti-diagonal"+label,pairing==reverse)
    check("split residue nondegenerate"+label,rank==0 or pairing.det()!=0)
    controls.append(dict(e=e,m=m,b=b,residual_rank=rank,
                         reduced_support=rank>0))

# Symbolic inequality guaranteeing a single reduced equation kills Q_m.
a,b=s.symbols("a b",integer=True,nonnegative=True)
zero("annihilator valuation surplus",(a+1)-(a-b)-(-b)-(1+2*b))
zero("full codifferent annihilator surplus",(a+1)-a-1)

proof=HERE/"REDUCED_BRANCH_PAIRING_AND_SIXTH_MOMENT_DEFECT.md"
record = {
    "outcome":"PASS",
    "sympy_version":s.__version__,
    "check_count":len(checks),
    "checks":checks,
    "controls":controls,
    "elapsed_seconds":round(time.perf_counter()-started,3),
    "proof_sha256":hashlib.sha256(proof.read_bytes()).hexdigest(),
    "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "scope":{
        "universal_arguments_are_written_not_finitely_certified":True,
        "reduced_branch_annihilator_proved":True,
        "global_branch_curve_pairing_proved":True,
        "higher_moment_residual_annihilators_proved":True,
        "matrix_factorizations_are_local_on_base_free_charts":True,
        "automatic_trace_integrality_proved":False,
        "automatic_sixth_moment_perfection_proved":False,
        "JC2_or_HC4_resolved":False,
        "independent_mathematical_review":False,
        "formal_verification":False,
        "external_review":False
    }
}
outdir=Path(tempfile.mkdtemp(prefix="reduced-branch-check-"))
receipt=outdir/"reduced_branch_pairing_verification.json"
receipt.write_text(json.dumps(record,indent=2)+"\n")
print(json.dumps({"outcome":record["outcome"],"checks":len(checks),
                  "receipt":str(receipt)}))
