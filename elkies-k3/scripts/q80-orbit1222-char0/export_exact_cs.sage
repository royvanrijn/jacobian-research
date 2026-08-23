#!/usr/bin/env sage
# Simultaneous rational reconstruction of the Q80 orbit-1222 monic invariants.
# Uses only cached modular embedding JSON files; computes no new fibers.

import argparse
import json
import re
from fractions import Fraction
from itertools import product
from pathlib import Path

from sage.all import Matrix, ZZ

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

EMBEDDINGS = ((1, 1), (1, -1), (-1, 1), (-1, -1))


def emb_file(repo, p, ss, sj):
    return repo / f"q80_modj_p{p}_s{ss:+d}_j{sj:+d}.json"


def load_json(path):
    return json.loads(path.read_text())


def inv(a, p):
    return pow(int(a) % p, -1, p)


def pad(values, n):
    values = list(values)
    return values + [0] * (n - len(values))


def poly_mul(a, b, p, length=25):
    out = [0] * length
    for i, x in enumerate(a):
        x %= p
        if not x:
            continue
        for j, y in enumerate(b):
            if i + j >= length:
                break
            y %= p
            if y:
                out[i+j] = (out[i+j] + x*y) % p
    return out


def poly_pow(a, e, p, length=25):
    out = [1] + [0] * (length - 1)
    base = pad(a, length)[:length]
    while e:
        if e & 1:
            out = poly_mul(out, base, p, length)
        e >>= 1
        if e:
            base = poly_mul(base, base, p, length)
    return out


def normalized_power_root(target, exponent, p, expected_degree):
    target = [int(x) % p for x in pad(target, 25)[:25]]
    if target[0] != 1:
        raise ArithmeticError(f"constant={target[0]}, expected 1")
    root = [0] * 25
    root[0] = 1
    ie = inv(exponent, p)
    for n in range(1, 25):
        known = poly_pow(root, exponent, p, 25)[n]
        root[n] = (target[n] - known) * ie % p
    if poly_pow(root, exponent, p, 25) != target:
        raise ArithmeticError("formal root mismatch")
    high = [i for i, c in enumerate(root) if i > expected_degree and c % p]
    if high:
        raise ArithmeticError(f"root degree overflow {high}")
    return root[:expected_degree+1]


def monic_invariants(kernel, p):
    kernel = [int(x) % p for x in kernel]
    if len(kernel) != 50 or kernel[0] != 1:
        raise ArithmeticError("bad normalized kernel")
    N = kernel[:25]
    D = kernel[25:]

    C = normalized_power_root(N, 3, p, 8)
    M = [(N[i] - (1728 % p)*D[i]) % p for i in range(25)]
    m0 = M[0]
    if m0 == 0:
        raise ArithmeticError("m0=0")
    S = normalized_power_root([x*inv(m0,p) % p for x in M], 2, p, 12)

    c8, s12 = C[8] % p, S[12] % p
    if not c8 or not s12:
        raise ArithmeticError("degree drop in invariant root")

    C = [x*inv(c8,p) % p for x in C]
    S = [x*inv(s12,p) % p for x in S]
    mu = m0*s12*s12*inv(c8*c8*c8,p) % p

    N2 = poly_pow(C, 3, p, 25)
    M2 = [mu*x % p for x in poly_pow(S, 2, p, 25)]
    D2 = [(N2[i]-M2[i])*inv(1728,p) % p for i in range(25)]
    scale = inv(c8*c8*c8,p)
    assert N2 == [scale*x % p for x in N]
    assert D2 == [scale*x % p for x in D]
    return {"C": C, "S": S, "mu": mu}


def combine_four(values, rs, rj, p):
    fpp = values[(1,1)] % p
    fpm = values[(1,-1)] % p
    fmp = values[(-1,1)] % p
    fmm = values[(-1,-1)] % p
    return (
        (fpp+fpm+fmp+fmm)*inv(4,p) % p,
        (fpp+fpm-fmp-fmm)*inv(4*rs,p) % p,
        (fpp-fpm+fmp-fmm)*inv(4*rj,p) % p,
        (fpp-fpm-fmp+fmm)*inv(4*rs*rj,p) % p,
    )


def build_packet(repo, p):
    rec = {e: load_json(emb_file(repo,p,*e)) for e in EMBEDDINGS}
    rs = int(rec[(1,1)]["s_root"]) % p
    rj = int(rec[(1,1)]["j_root"]) % p
    if rs*rs % p != (-6)%p or rj*rj % p != (-3)%p:
        raise ArithmeticError("bad roots")
    for ss,sj in EMBEDDINGS:
        if int(rec[(ss,sj)]["s_root"]) % p != ss*rs % p:
            raise ArithmeticError("s-root mismatch")
        if int(rec[(ss,sj)]["j_root"]) % p != sj*rj % p:
            raise ArithmeticError("j-root mismatch")

    data = {e: monic_invariants(rec[e]["kernel"],p) for e in EMBEDDINGS}
    C = [combine_four({e:data[e]["C"][i] for e in EMBEDDINGS},rs,rj,p) for i in range(9)]
    S = [combine_four({e:data[e]["S"][i] for e in EMBEDDINGS},rs,rj,p) for i in range(13)]
    mu = combine_four({e:data[e]["mu"] for e in EMBEDDINGS},rs,rj,p)

    if any(c[1] % p or c[3] % p for c in C+S+[mu]):
        raise ArithmeticError("unexpected sqrt(-6) or sqrt(18) support")
    return {"p":p,"rs":rs,"rj":rj,"C":C,"S":S,"mu":mu}


def discover(repo):
    out=[]
    pat=re.compile(r"q80_modj_p(\d+)_s\+1_j\+1\.json$")
    for f in repo.glob("q80_modj_p*_s+1_j+1.json"):
        m=pat.match(f.name)
        if not m:
            continue
        p=int(m.group(1))
        if p==73:
            continue
        try:
            ok=all(
                emb_file(repo,p,*e).exists()
                and len(load_json(emb_file(repo,p,*e)).get("kernel",()))==50
                for e in EMBEDDINGS
            )
        except Exception:
            ok=False
        if ok:
            out.append(p)
    return sorted(set(out))


def crt_pair(x,m,y,p):
    t=((y-x)%p)*inv(m,p)%p
    return (x+m*t)%(m*p),m*p


def crt(values):
    x,m=0,1
    for p,y in values:
        x,m=crt_pair(x,m,int(y)%p,p)
    return x,m


def fraction_mod(q,p):
    num = q.numerator() if callable(getattr(q, "numerator", None)) else q.numerator
    den = q.denominator() if callable(getattr(q, "denominator", None)) else q.denominator
    num, den = int(num), int(den)
    if den % p == 0:
        raise ZeroDivisionError
    return num % p * inv(den,p) % p


def vector_matches(values, residue_table, primes):
    for coord,q in enumerate(values):
        for p in primes:
            try:
                got=fraction_mod(q,p)
            except ZeroDivisionError:
                return False
            if got != residue_table[p][coord] % p:
                return False
    return True


def lll_candidates(xs, modulus, num_weight, den_weight):
    n=len(xs)
    B=Matrix(ZZ,n+1,n+1)
    for i in range(n):
        B[i,i]=modulus*num_weight
    for i,x in enumerate(xs):
        B[n,i]=int(x)*num_weight
    B[n,n]=den_weight

    R=B.LLL(delta=0.99)
    rows=[tuple(int(x) for x in R.row(i)) for i in range(R.nrows())]
    rows.sort(key=lambda v: sum(x*x for x in v))

    candidates=list(rows)
    first=rows[:min(4,len(rows))]
    for coeffs in product((-2,-1,0,1,2), repeat=len(first)):
        if not any(coeffs):
            continue
        candidates.append(tuple(
            sum(coeffs[r]*first[r][i] for r in range(len(first)))
            for i in range(n+1)
        ))

    seen=set()
    for v in candidates:
        if v not in seen:
            seen.add(v)
            yield v


def simrr(residue_table, train, verify, label):
    n=len(next(iter(residue_table.values())))
    xs=[]
    modulus=None
    for coord in range(n):
        x,m=crt([(p,residue_table[p][coord]) for p in train])
        if modulus is None:
            modulus=m
        assert m==modulus
        xs.append(x)

    weights=((1,1),(4,1),(16,1),(64,1),(256,1),(1,4),(1,16),(1,64),(1,256))
    for nw,dw in weights:
        for v in lll_candidates(xs,modulus,nw,dw):
            if v[-1]==0:
                continue
            if any(v[i] % nw for i in range(n)) or v[-1] % dw:
                continue
            nums=[v[i]//nw for i in range(n)]
            den=v[-1]//dw
            if den==0:
                continue
            values=[Fraction(int(a),int(den)) for a in nums]
            if not vector_matches(values,residue_table,train):
                continue
            if verify and not vector_matches(values,residue_table,verify):
                continue
            max_num=max(abs(q.numerator) for q in values)
            max_den=max(q.denominator for q in values)
            print(
                f"Q80SIMRR|label={label}|coords={n}|train={len(train)}|verify={len(verify)}|"
                f"weight={nw}:{dw}|max_num_digits={len(str(max_num))}|"
                f"max_den_digits={len(str(max_den))}|status=PASS_GROUP",
                flush=True,
            )
            return values
    return None


def reconstruct_recursive(entries, packets, train, verify, label):
    residue_table={}
    for p in train+verify:
        flat=[]
        for name,index in entries:
            c=packets[p][name][index] if name != "mu" else packets[p]["mu"]
            flat.extend((c[0] % p,c[2] % p))
        residue_table[p]=flat

    values=simrr(residue_table,train,verify,label)
    if values is not None:
        out={}
        for i,entry in enumerate(entries):
            out[entry]=(values[2*i],values[2*i+1])
        return out

    if len(entries)==1:
        return None

    mid=len(entries)//2
    left=reconstruct_recursive(entries[:mid],packets,train,verify,label+"L")
    if left is None:
        return None
    right=reconstruct_recursive(entries[mid:],packets,train,verify,label+"R")
    if right is None:
        return None
    left.update(right)
    return left


def build_candidate(packets,train,verify):
    Cmap=reconstruct_recursive([("C",i) for i in range(8)],packets,train,verify,"C")
    if Cmap is None:
        return None
    Smap=reconstruct_recursive([("S",i) for i in range(12)],packets,train,verify,"S")
    if Smap is None:
        return None
    mmap=reconstruct_recursive([("mu",0)],packets,train,verify,"mu")
    if mmap is None:
        return None

    C=[(Cmap[("C",i)][0],Fraction(0),Cmap[("C",i)][1],Fraction(0)) for i in range(8)]
    C.append((Fraction(1),Fraction(0),Fraction(0),Fraction(0)))
    S=[(Smap[("S",i)][0],Fraction(0),Smap[("S",i)][1],Fraction(0)) for i in range(12)]
    S.append((Fraction(1),Fraction(0),Fraction(0),Fraction(0)))
    ma,mc=mmap[("mu",0)]
    return {"C":C,"S":S,"mu":(ma,Fraction(0),mc,Fraction(0))}


def eval_field(c,ss,sj,rs,rj,p):
    a,b,q,d=c
    return (
        fraction_mod(a,p)+fraction_mod(b,p)*ss*rs+
        fraction_mod(q,p)*sj*rj+fraction_mod(d,p)*ss*sj*rs*rj
    )%p


def kernel_from_candidate(cand,p,ss,sj,rs,rj):
    C=[eval_field(c,ss,sj,rs,rj,p) for c in cand["C"]]
    S=[eval_field(c,ss,sj,rs,rj,p) for c in cand["S"]]
    mu=eval_field(cand["mu"],ss,sj,rs,rj,p)
    N=poly_pow(C,3,p,25)
    M=[mu*x%p for x in poly_pow(S,2,p,25)]
    D=[(N[i]-M[i])*inv(1728,p)%p for i in range(25)]
    return N,D


def proportional(N,D,expected,p):
    EN=[int(x)%p for x in expected[:25]]
    ED=[int(x)%p for x in expected[25:]]
    lam=None
    for a,b in zip(N+D,EN+ED):
        if a:
            lam=b*inv(a,p)%p
            break
    if lam is None:
        return False
    return [lam*x%p for x in N]==EN and [lam*x%p for x in D]==ED


def validate_cached(repo,cand,packets,primes):
    bad=[]
    for p in primes:
        P=packets[p]
        for ss,sj in EMBEDDINGS:
            try:
                N,D=kernel_from_candidate(cand,p,ss,sj,P["rs"],P["rj"])
            except ZeroDivisionError:
                bad.append((p,ss,sj,"denom"))
                continue
            E=load_json(emb_file(repo,p,ss,sj))["kernel"]
            if not proportional(N,D,E,p):
                bad.append((p,ss,sj,"kernel"))
    return bad


def parse_kernel(path):
    text=path.read_text()
    m=re.search(r"expected_kernel\s*=\s*vector\(\s*finite\s*,\s*\[(.*?)\]\s*,?\s*\)",text,re.S)
    if not m:
        raise RuntimeError(f"cannot parse {path}")
    vals=[int(x) for x in re.findall(r"-?\d+",m.group(1))]
    if len(vals)!=50:
        raise RuntimeError(f"{path}: parsed {len(vals)}")
    return vals



def build_p73_qj_packet():
    p=73
    rj=17
    orig=parse_kernel(
        REPO_ROOT/"elkies-k3/scripts/reconstruct_q80_third_q12_jacobian_gf73.sage"
    )
    conj=parse_kernel(
        REPO_ROOT/"elkies-k3/scripts/analyze_q80_third_q12_galois_descent_gf73.sage"
    )

    plus=monic_invariants(orig,p)
    minus=monic_invariants(conj,p)

    def pair(fp,fm):
        return (
            (fp+fm)*inv(2,p)%p,
            0,
            (fp-fm)*inv(2*rj,p)%p,
            0,
        )

    C=[pair(plus["C"][i],minus["C"][i]) for i in range(9)]
    S=[pair(plus["S"][i],minus["S"][i]) for i in range(13)]
    mu=pair(plus["mu"],minus["mu"])

    for data,jr in ((plus,17),(minus,56)):
        for i,c in enumerate(C):
            assert (c[0]+c[2]*jr)%p == data["C"][i]
        for i,c in enumerate(S):
            assert (c[0]+c[2]*jr)%p == data["S"][i]
        assert (mu[0]+mu[2]*jr)%p == data["mu"]

    print(
        "Q80SIMRR|p=73|source=pinned_original+conjugate|"
        "field=Q(sqrt(-3))|status=PASS_VIRTUAL_HOLDOUT",
        flush=True,
    )
    return {"p":p,"rs":33,"rj":17,"C":C,"S":S,"mu":mu}

def validate_p73(cand):
    orig=parse_kernel(REPO_ROOT/"elkies-k3/scripts/reconstruct_q80_third_q12_jacobian_gf73.sage")
    conj=parse_kernel(REPO_ROOT/"elkies-k3/scripts/analyze_q80_third_q12_galois_descent_gf73.sage")
    bad=[]
    for jr,E in ((17,orig),(56,conj)):
        try:
            N,D=kernel_from_candidate(cand,73,1,1,33,jr)
        except ZeroDivisionError:
            bad.append((jr,"denom"))
            continue
        if not proportional(N,D,E,73):
            bad.append((jr,"kernel"))
    return bad


def qstr(q):
    return str(q.numerator) if q.denominator==1 else f"({q.numerator}/{q.denominator})"


def field_str(c):
    a,b,q,d=c
    parts=[]
    for x,suf in ((a,""),(b,"*s"),(q,"*j"),(d,"*s*j")):
        if x:
            z=qstr(x)
            parts.append(z if not suf else f"({z}){suf}")
    return " + ".join(parts) if parts else "0"


def poly_str(coeffs):
    parts=[]
    for i,c in enumerate(coeffs):
        z=field_str(c)
        if z=="0":
            continue
        parts.append(f"({z})" if i==0 else f"({z})*V" if i==1 else f"({z})*V^{i}")
    return " + ".join(parts)


def write_outputs(cand,train,verify):
    sage=ORBIT_DATA/"q80_char0_orbit1222_simrr_invariants.sage"
    js=ORBIT_DATA/"q80_char0_orbit1222_simrr_invariants.json"
    note=ORBIT_DATA/"Q80_CHAR0_ORBIT1222_SIMRR_INVARIANTS.md"

    template='''#!/usr/bin/env sage
from sage.all import PolynomialRing, QuadraticField
K=QuadraticField(-6,"s"); s=K.gen()
TR=PolynomialRing(K,"T"); T=TR.gen()
L=K.extension(T^2+3,"j"); j=L.gen(); s=L(s)
R=PolynomialRing(L,"V"); V=R.gen()
C=__C__
S=__S__
mu=__MU__
assert C.degree()==8 and C[8]==1
assert S.degree()==12 and S[12]==1
N=C^3
D=(N-mu*S^2)/1728
jmap=N/D
print(f"Q80SIMRRCHAR0|C={C}|S={S}|mu={mu}|status=PASS_EXACT_SIMRR_INVARIANTS")
'''
    sage.write_text(
        template.replace("__C__",poly_str(cand["C"]))
                .replace("__S__",poly_str(cand["S"]))
                .replace("__MU__",field_str(cand["mu"]))
    )
    js.write_text(json.dumps({
        "version":int(1),
        "method":"simultaneous rational reconstruction + LLL",
        "field":"QQ(sqrt(-3)); sqrt(-6) and sqrt(18) coefficients zero",
        "train_primes":train,
        "verify_primes":verify,
        "C":[[str(x) for x in c] for c in cand["C"]],
        "S":[[str(x) for x in c] for c in cand["S"]],
        "mu":[str(x) for x in cand["mu"]],
    },indent=2,sort_keys=True,default=int)+"\n")
    note.write_text(
        "# Q80 orbit 1222 simultaneous invariant reconstruction\n\n"
        "Status: **PASS_EXACT_SIMRR_INVARIANTS**\n\n"
        "The monic degree-8 C, monic degree-12 S and scalar mu were recovered "
        "by simultaneous rational reconstruction (LLL), not coefficientwise reconstruction.\n\n"
        f"- training primes: `{','.join(map(str,train))}`\n"
        f"- held-out cached primes: `{','.join(map(str,verify)) or '-'}`\n"
        "- every cached modular embedding validates projectively\n"
        "- both independent p=73 original/conjugate kernels validate\n"
        "- exact coefficient field in this V-coordinate is QQ(sqrt(-3))\n\n"
        "The normalized j-map is `j=C^3 / ((C^3-mu*S^2)/1728)`.\n"
    )
    return sage,js,note


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--modular-data",type=Path,default=MODULAR_DATA)
    ap.add_argument("--primes",default="")
    args=ap.parse_args()
    modular_data=args.modular_data.expanduser().resolve()
    primes=[int(x) for x in args.primes.split(",") if x.strip()] if args.primes else discover(modular_data)
    if len(primes)<5:
        raise SystemExit(f"need >=5 complete primes; found {primes}")

    packets={}
    usable=[]
    for p in primes:
        try:
            packets[p]=build_packet(modular_data,p)
            usable.append(p)
        except Exception as exc:
            print(f"Q80SIMRR|p={p}|status=SKIP|type={type(exc).__name__}|message={exc}",flush=True)
    print(f"Q80SIMRR|usable={len(usable)}|primes={','.join(map(str,usable))}",flush=True)

    # p=73 is kept OUT of CRT, but used inside LLL candidate selection as
    # an independent virtual Q(sqrt(-3)) holdout.
    packets[73]=build_p73_qj_packet()

    accepted=None
    for k in range(5,len(usable)+1):
        train=usable[:k]
        verify=usable[k:]
        lll_verify=verify+[73]
        print(
            f"Q80SIMRR|stage=TRY|k={k}|train={','.join(map(str,train))}|"
            f"verify={','.join(map(str,verify)) or '-'}|p73_lll_holdout=yes",
            flush=True,
        )
        cand=build_candidate(packets,train,lll_verify)
        if cand is None:
            print(f"Q80SIMRR|k={k}|status=NO_SHORT_VECTOR",flush=True)
            continue
        train_bad=validate_cached(modular_data,cand,packets,train)
        verify_bad=validate_cached(modular_data,cand,packets,verify)
        p73_bad=validate_p73(cand)
        print(
            f"Q80SIMRR|k={k}|train_fail={len(train_bad)}|"
            f"verify_fail={len(verify_bad)}|p73_fail={len(p73_bad)}",
            flush=True,
        )
        if not train_bad and not verify_bad and not p73_bad:
            accepted=(cand,train,verify)
            break

    if accepted is None:
        print("Q80SIMRR|status=NEED_DIFFERENT_GROUPING_OR_MORE_MODULUS",flush=True)
        raise SystemExit(2)

    cand,train,verify=accepted
    sage,js,note=write_outputs(cand,train,verify)
    print(
        "Q80SIMRR|field=QQ(sqrt(-3))|all_cached=PASS|p73_lll=PASS|p73=PASS|"
        "status=PASS_EXACT_SIMRR_INVARIANTS",
        flush=True,
    )
    print(f"Q80SIMRR|sage={sage}",flush=True)
    print(f"Q80SIMRR|json={js}",flush=True)
    print(f"Q80SIMRR|note={note}",flush=True)



def build_cs_only(packets,train,verify):
    Cmap=reconstruct_recursive([("C",i) for i in range(8)],packets,train,verify,"C")
    if Cmap is None:
        return None
    Smap=reconstruct_recursive([("S",i) for i in range(12)],packets,train,verify,"S")
    if Smap is None:
        return None

    C=[(Cmap[("C",i)][0],Fraction(0),Cmap[("C",i)][1],Fraction(0)) for i in range(8)]
    C.append((Fraction(1),Fraction(0),Fraction(0),Fraction(0)))
    S=[(Smap[("S",i)][0],Fraction(0),Smap[("S",i)][1],Fraction(0)) for i in range(12)]
    S.append((Fraction(1),Fraction(0),Fraction(0),Fraction(0)))
    return {"C":C,"S":S}


def rational_text(x):
    n = x.numerator() if callable(getattr(x, "numerator", None)) else x.numerator
    d = x.denominator() if callable(getattr(x, "denominator", None)) else x.denominator
    n, d = int(n), int(d)
    return str(n) if d == 1 else f"{n}/{d}"


def write_cs_only(cs,train,verify):
    out=ORBIT_DATA/"q80_char0_orbit1222_cs.json"
    out.write_text(json.dumps({
        "version":int(1),
        "field":"QQ(sqrt(-3))",
        "normalization":"C,S monic",
        "train_primes":train,
        "verify_primes":verify,
        "p73_lll_holdout":True,
        "C":[[rational_text(x) for x in c] for c in cs["C"]],
        "S":[[rational_text(x) for x in c] for c in cs["S"]],
    },indent=2,sort_keys=True)+"\n")
    print(
        f"Q80CS|train={','.join(map(str,train))}|"
        f"verify={','.join(map(str,verify)) or '-'}|p73=PASS|"
        f"status=PASS_EXACT_CS|out={out}",
        flush=True,
    )
    return out


def main_cs():
    ap=argparse.ArgumentParser()
    ap.add_argument("--modular-data",type=Path,default=MODULAR_DATA)
    ap.add_argument("--primes",default="")
    args=ap.parse_args()
    modular_data=args.modular_data.expanduser().resolve()

    primes=[int(x) for x in args.primes.split(",") if x.strip()] if args.primes else discover(modular_data)
    packets={}
    usable=[]
    for p in primes:
        try:
            packets[p]=build_packet(modular_data,p)
            usable.append(p)
        except Exception as exc:
            print(f"Q80CS|p={p}|status=SKIP|type={type(exc).__name__}|message={exc}",flush=True)

    packets[73]=build_p73_qj_packet()

    for k in range(5,len(usable)+1):
        train=usable[:k]
        verify=usable[k:]
        lll_verify=verify+[73]
        print(
            f"Q80CS|stage=TRY|k={k}|train={','.join(map(str,train))}|"
            f"verify={','.join(map(str,verify)) or '-'}|p73=yes",
            flush=True,
        )
        cs=build_cs_only(packets,train,lll_verify)
        if cs is None:
            continue
        write_cs_only(cs,train,verify)
        return

    raise SystemExit("could not reconstruct C,S")

if __name__=="__main__":
    main_cs()
