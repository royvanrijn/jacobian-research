#!/usr/bin/env sage
# Recover exact mu for Q80 orbit 1222 from factors of
# R = 3*C'*S - 2*C*S', using certified exact C,S.
#
# No resultant, no maximal orders, no new modular fibers.

import json
import re
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, QuadraticField

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
ORBIT_DATA = (
    REPO_ROOT
    / "elkies-k3"
    / "data"
    / "fibrations"
    / "q80-orbit1222-char0"
)
MODULAR_DATA = (
    REPO_ROOT
    / "artifacts"
    / "generated-results"
    / "q80-orbit1222-char0"
    / "modular"
)
CSFILE = ORBIT_DATA / "q80_char0_orbit1222_cs.json"

if not CSFILE.exists():
    raise SystemExit(f"missing certified C,S cache: {CSFILE}")

data = json.loads(CSFILE.read_text())

K = QuadraticField(-3, "j")
j = K.gen()
KV = PolynomialRing(K, "V")
V = KV.gen()


def row_to_K(row):
    a = QQ(row[0])
    s = QQ(row[1])
    b = QQ(row[2])
    sj = QQ(row[3])
    if s != 0 or sj != 0:
        raise ArithmeticError(
            f"C/S cache unexpectedly uses sqrt(-6) or sqrt(18): {row}"
        )
    return K(a + b*j)


C = KV([row_to_K(row) for row in data["C"]])
S = KV([row_to_K(row) for row in data["S"]])

assert C.degree() == 8 and C[8] == 1
assert S.degree() == 12 and S[12] == 1

C3 = C**3
S2 = S**2
Rcrit = 3*C.derivative()*S - 2*C*S.derivative()

print(
    "Q80CRITFACT|"
    f"C_degree={C.degree()}|S_degree={S.degree()}|"
    f"R_degree={Rcrit.degree()}|stage=READY",
    flush=True,
)

t0 = time.monotonic()
fac = tuple(Rcrit.factor())
print(
    "Q80CRITFACT|"
    f"stage=FACTOR_DONE|seconds={time.monotonic()-t0:.3f}|"
    f"factor_count={len(fac)}|"
    f"degrees={','.join(str(g.degree())+('x'+str(m) if m != 1 else '') for g,m in fac)}",
    flush=True,
)

groups = {}
nonconstant = []
infinite = []

for index, (g, multiplicity) in enumerate(fac, start=1):
    A = C3.mod(g)
    B = S2.mod(g)

    if B == 0:
        infinite.append((index, int(g.degree()), int(multiplicity)))
        print(
            f"Q80CRITFACT|factor={index}|degree={g.degree()}|mult={multiplicity}|"
            "critical_value=INFINITY|status=POLE_CRITICAL_FACTOR",
            flush=True,
        )
        continue

    pivot = next(i for i in range(B.degree()+1) if B[i] != 0)
    lam = K(A[pivot] / B[pivot])

    if A != lam*B:
        nonconstant.append((index, int(g.degree()), int(multiplicity)))
        print(
            f"Q80CRITFACT|factor={index}|degree={g.degree()}|mult={multiplicity}|"
            "critical_value=NONCONSTANT_IN_QUOTIENT|status=OTHER_CRITICAL_FACTOR",
            flush=True,
        )
        continue

    contribution = int(g.degree()) * int(multiplicity)
    groups.setdefault(lam, []).append(
        {
            "factor_index": index,
            "degree": int(g.degree()),
            "multiplicity_in_R": int(multiplicity),
            "contribution": contribution,
            "factor": str(g),
        }
    )

    print(
        f"Q80CRITFACT|factor={index}|degree={g.degree()}|mult={multiplicity}|"
        f"critical_value={lam}|status=PASS_CONSTANT_CRITICAL_VALUE",
        flush=True,
    )

summary = []
for lam, entries in groups.items():
    degree_sum = sum(e["contribution"] for e in entries)
    summary.append((degree_sum, lam, entries))
    print(
        "Q80CRITFACT|"
        f"group_value={lam}|degree_sum={degree_sum}|"
        f"factors={','.join(str(e['factor_index']) for e in entries)}",
        flush=True,
    )

hits = [(deg, lam, entries) for deg, lam, entries in summary if deg == 15]
if len(hits) != 1:
    print(
        "Q80CRITFACT|"
        f"status=AMBIGUOUS_DEGREE15_GROUPS|count={len(hits)}|"
        f"all_group_degrees={','.join(str(x[0]) for x in summary)}",
        flush=True,
    )
    raise SystemExit(2)

_, mu, mu_entries = hits[0]
mu = K(mu)

print(
    f"Q80CRITFACT|mu={mu}|"
    f"mu_factors={','.join(str(e['factor_index']) for e in mu_entries)}|"
    "stage=EXACT_MU",
    flush=True,
)

D0 = C3 - mu*S2
assert D0.degree() == 24

g = D0
gcd_chain = []
while True:
    g = g.gcd(g.derivative()).monic()
    gcd_chain.append(int(g.degree()))
    if g.degree() == 0:
        break
    if len(gcd_chain) > 20:
        raise ArithmeticError("unexpectedly long discriminant gcd chain")

expected_chain = [15, 10, 8, 6, 4, 2, 0]
if gcd_chain != expected_chain:
    raise ArithmeticError(
        f"wrong discriminant gcd chain {gcd_chain}; expected {expected_chain}"
    )

print(
    "Q80CRITFACT|"
    f"D_degree={D0.degree()}|"
    f"gcd_chain={','.join(map(str,gcd_chain))}|"
    "repeated_fibers=2I7+3I2|stage=PASS_EXACT_MULTIPLICITIES",
    flush=True,
)


def parse_expected_kernel(path):
    text = path.read_text()
    match = re.search(
        r"expected_kernel\s*=\s*vector\(\s*finite\s*,\s*\[(.*?)\]\s*,?\s*\)",
        text,
        re.S,
    )
    if not match:
        raise RuntimeError(f"could not parse expected_kernel from {path}")
    values = [int(x) for x in re.findall(r"-?\d+", match.group(1))]
    if len(values) != 50:
        raise RuntimeError(f"{path}: expected 50 values, got {len(values)}")
    return values


def pad(values, n):
    values = list(values)
    return values + [0]*(n-len(values))


def pmul(a,b,p,n=25):
    out=[0]*n
    for i,x in enumerate(a):
        for k,y in enumerate(b):
            if i+k<n:
                out[i+k]=(out[i+k]+x*y)%p
    return out


def ppow(a,e,p,n=25):
    out=[1]+[0]*(n-1)
    base=pad(a,n)
    while e:
        if e&1:
            out=pmul(out,base,p,n)
        e//=2
        if e:
            base=pmul(base,base,p,n)
    return out


def proot(target,e,p,degree):
    target=[x%p for x in target]
    root=[0]*25
    root[0]=1
    ie=pow(e,-1,p)
    for n in range(1,25):
        known=ppow(root,e,p,25)[n]
        root[n]=(target[n]-known)*ie%p
    if any(root[n] for n in range(degree+1,25)):
        raise ArithmeticError("root degree overflow")
    return root[:degree+1]


def kernel_mu(kernel,p):
    N=[x%p for x in kernel[:25]]
    D=[x%p for x in kernel[25:]]
    C0=proot(N,3,p,8)
    Q=[(N[i]-1728*D[i])%p for i in range(25)]
    m0=Q[0]
    S0=proot([x*pow(m0,-1,p)%p for x in Q],2,p,12)
    return m0*S0[12]**2*pow(C0[8]**3,-1,p)%p


orig = parse_expected_kernel(
    REPO_ROOT / "elkies-k3/scripts/reconstruct_q80_third_q12_jacobian_gf73.sage"
)
conj = parse_expected_kernel(
    REPO_ROOT / "elkies-k3/scripts/analyze_q80_third_q12_galois_descent_gf73.sage"
)
up = kernel_mu(orig,73)
um = kernel_mu(conj,73)

F73 = GF(73)
mu_plus = F73(QQ(mu[0])) + F73(17)*F73(QQ(mu[1]))
mu_minus = F73(QQ(mu[0])) - F73(17)*F73(QQ(mu[1]))

if int(mu_plus) != up or int(mu_minus) != um:
    raise ArithmeticError(
        f"p73 mismatch: exact -> {int(mu_plus)},{int(mu_minus)}; "
        f"pinned -> {up},{um}"
    )

print(
    f"Q80CRITFACT|p73_plus={up}|p73_minus={um}|"
    "status=PASS_P73_INDEPENDENT",
    flush=True,
)

fast_checked = 0
for path in sorted(MODULAR_DATA.glob("q80_mu_fast_p*.json")):
    record = json.loads(path.read_text())
    p = int(record["prime"])
    Fp = GF(p)
    jr = Fp(int(record["j_positive_root"]))
    exact_plus = Fp(QQ(mu[0])) + jr*Fp(QQ(mu[1]))
    exact_minus = Fp(QQ(mu[0])) - jr*Fp(QQ(mu[1]))

    if int(exact_plus) != int(record["mu_plus_j"]):
        raise ArithmeticError(f"{path.name}: +j mu mismatch")
    if int(exact_minus) != int(record["mu_minus_j"]):
        raise ArithmeticError(f"{path.name}: -j mu mismatch")
    fast_checked += 1

print(
    f"Q80CRITFACT|fast_mu_files_checked={fast_checked}|"
    "status=PASS_ALL_FAST_MU",
    flush=True,
)

OUT_SAGE = ORBIT_DATA / "q80_char0_orbit1222_mu_critical_factors.sage"
OUT_JSON = ORBIT_DATA / "q80_char0_orbit1222_mu_critical_factors.json"
OUT_NOTE = ORBIT_DATA / "Q80_CHAR0_ORBIT1222_CRITICAL_FACTOR_MU.md"

sage_lines = [
    "#!/usr/bin/env sage",
    "from sage.all import PolynomialRing, QuadraticField",
    'K = QuadraticField(-3, "j")',
    "j = K.gen()",
    'R = PolynomialRing(K, "V")',
    "V = R.gen()",
    f"C = {C}",
    f"S = {S}",
    f"mu = {mu}",
    "Delta = (C^3-mu*S^2)/1728",
    "jmap = C^3/Delta",
    "assert C.degree() == 8 and C[8] == 1",
    "assert S.degree() == 12 and S[12] == 1",
    "assert Delta.degree() == 24",
    'print(f"Q80ORBIT1222EXACTINVARIANTS|mu={mu}|Delta_degree={Delta.degree()}|status=PASS_EXACT_CRITICAL_FACTOR_MU")',
]
OUT_SAGE.write_text("\n".join(sage_lines) + "\n")

payload = {
    "version": 1,
    "method": "factor critical polynomial and group constant critical values",
    "field": "QQ(sqrt(-3))",
    "mu": str(mu),
    "critical_factor_degrees": [
        {"degree": int(g.degree()), "multiplicity": int(m)}
        for g,m in fac
    ],
    "mu_factor_indices": [e["factor_index"] for e in mu_entries],
    "mu_degree_sum": 15,
    "discriminant_gcd_chain": gcd_chain,
    "p73_plus": up,
    "p73_minus": um,
    "fast_mu_files_checked": fast_checked,
}
OUT_JSON.write_text(
    json.dumps(payload, indent=2, sort_keys=True, default=int) + "\n"
)

OUT_NOTE.write_text(
    "# Q80 orbit 1222 — exact mu from critical factors\n\n"
    "Status: **PASS_EXACT_CRITICAL_FACTOR_MU**\n\n"
    "The degree-18 critical polynomial `R=3*C'*S-2*C*S'` was factored "
    "over `QQ(sqrt(-3))`. For each irreducible factor `g`, the reductions "
    "of `C^3` and `S^2` modulo `g` were tested for constant proportionality. "
    "The unique critical value whose factor-degrees total 15 is the exact mu "
    "for the repeated fibers `2 I7 + 3 I2`.\n\n"
    f"- exact mu: `{mu}`\n"
    f"- discriminant gcd chain: `{gcd_chain}`\n"
    f"- fast modular mu files checked: `{fast_checked}`\n"
    "- both independent p=73 Galois embeddings validate\n\n"
    "No further modular primes are needed for mu.\n"
)

print(
    f"Q80CRITFACT|out_sage={OUT_SAGE}|out_json={OUT_JSON}|"
    f"out_note={OUT_NOTE}|status=PASS_EXACT_CRITICAL_FACTOR_MU",
    flush=True,
)
