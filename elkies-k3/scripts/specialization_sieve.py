from __future__ import annotations

from pathlib import Path
import argparse
import json
import math
import heapq


def polyval(coeffs,t,p):
    """coefficients low-to-high."""
    z=0
    for a in reversed(coeffs):
        z=(z*t+a)%p
    return z


def legendre(a,p):
    a%=p
    if a==0:
        return 0
    z=pow(a,(p-1)//2,p)
    return -1 if z==p-1 else 1


def count_curve(model,t,p):
    """
    General Weierstrass:
      y^2 + a1*x*y + a3*y
        = x^3 + a2*x^2 + a4*x + a6

    For fixed x, quadratic in y has discriminant

      D=(a1*x+a3)^2
        +4*(x^3+a2*x^2+a4*x+a6)

    # solutions in y = 1 + chi(D).
    """

    a=[
        polyval(model[k],t,p)
        for k in ("a1","a2","a3","a4","a6")
    ]

    a1,a2,a3,a4,a6=a

    chi_sum=0

    for x in range(p):
        lhs=(a1*x+a3)%p

        rhs=(
            x*x%p*x
            +a2*x*x
            +a4*x
            +a6
        )%p

        D=(lhs*lhs+4*rhs)%p
        chi_sum += legendre(D,p)

    # p affine base x values + point at infinity.
    N=p+1+chi_sum
    return N


def primes_up_to(n):
    sieve=[True]*(n+1)
    sieve[:2]=[False,False]

    for q in range(2,int(n**0.5)+1):
        if sieve[q]:
            sieve[q*q:n+1:q]=[False]*(
                ((n-q*q)//q)+1
            )

    return [
        q for q in range(3,n+1)
        if sieve[q]
    ]


ap=argparse.ArgumentParser()
ap.add_argument("--model",required=True)
ap.add_argument("--prime-bound",type=int,default=200)
ap.add_argument("--num-bound",type=int,default=5000)
ap.add_argument("--den-bound",type=int,default=100)
ap.add_argument("--top",type=int,default=500)
ap.add_argument("--out",required=True)
args=ap.parse_args()

model=json.loads(Path(args.model).read_text())

required=("a1","a2","a3","a4","a6")
for k in required:
    if k not in model:
        raise SystemExit(f"missing model coefficient {k}")

primes=primes_up_to(args.prime_bound)

print("primes =",len(primes),flush=True)

# ------------------------------------------------------------
# score[p][r] = log(N_p(E_r)/p)
# ------------------------------------------------------------

tables={}

for p in primes:
    T=[]

    for r in range(p):
        N=count_curve(model,r,p)

        # Singular fibers should eventually be filtered with
        # discriminant(t), but avoid log(0) regardless.
        if N<=0:
            s=1000.0
        else:
            s=math.log(N/p)

        T.append(s)

    tables[p]=T

    print("precomputed p =",p,flush=True)

# ------------------------------------------------------------
# Rational t = a/b.
#
# Smaller score is better according to the logarithmic
# Mestre/BSD heuristic as written in Elkies' notes.
# ------------------------------------------------------------

best=[]

for b in range(1,args.den_bound+1):

    for a in range(-args.num_bound,args.num_bound+1):

        if math.gcd(a,b)!=1:
            continue

        score=0.0
        used=0

        for p in primes:
            if b%p==0:
                continue

            r=(a*pow(b,-1,p))%p

            score += tables[p][r]
            used += 1

        if not used:
            continue

        # Use raw accumulated score; candidate comparison with
        # varying skipped primes can also use score/used.
        key=score

        item=(-key,a,b,score,used)

        if len(best)<args.top:
            heapq.heappush(best,item)
        elif item>best[0]:
            heapq.heapreplace(best,item)

results=sorted(
    [
        (score,a,b,used)
        for _,a,b,score,used in best
    ],
    key=lambda x:x[0]
)

out=Path(args.out)
out.parent.mkdir(parents=True,exist_ok=True)

with out.open("w") as f:
    f.write("score\tnum\tden\tused_primes\n")

    for score,a,b,used in results:
        f.write(
            f"{score:.12f}\t{a}\t{b}\t{used}\n"
        )

print()
print("TOP 25")

for x in results[:25]:
    print(x)

print("saved",out)
