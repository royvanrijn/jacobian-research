#!/usr/bin/env python3
from math import gcd, isqrt
import argparse

ap = argparse.ArgumentParser(description="Enumerate imaginary quadratic orders optimally embeddable in the Eichler order (D,N)=(6,79).")
ap.add_argument("--disc-max", type=int, default=20000, help="search negative order discriminants Delta with |Delta| <= this")
ap.add_argument("--embedding-count", type=int, default=4, help="highlight this many optimal embedding classes; Elkies t=2 orbit has four points")
ap.add_argument("--show-all", action="store_true")
ap.add_argument("--top", type=int, default=100)
args = ap.parse_args()

DQUAT = (2,3)
LEVEL = (79,)

def factorint(n):
    n=abs(int(n)); out=[]
    p=2
    while p*p<=n:
        if n%p==0:
            e=0
            while n%p==0:
                n//=p; e+=1
            out.append((p,e))
        p=3 if p==2 else p+2
    if n>1: out.append((n,1))
    return out

def fundamental_disc_and_conductor(Delta):
    # Delta = D_K f^2, D_K fundamental negative discriminant.
    assert Delta < 0 and Delta % 4 in (0,1)
    A=abs(Delta)
    # Try f descending; remaining discriminant must be fundamental.
    for f in range(isqrt(A),0,-1):
        if A%(f*f): continue
        dk=Delta//(f*f)
        if is_fundamental_discriminant(dk):
            return dk,f
    raise ValueError(Delta)

def is_squarefree(n):
    return all(e==1 for _,e in factorint(n))

def is_fundamental_discriminant(d):
    if d>=0: return False
    if d%4==1:
        return is_squarefree(abs(d))
    if d%4==0:
        m=d//4
        # m squarefree and m == 2 or 3 mod 4
        return is_squarefree(abs(m)) and (m%4 in (2,3))
    return False

def kronecker_prime(a,p):
    # Kronecker (a/p) for p prime, sufficient for p=2,3,79.
    a=int(a)
    if p==2:
        if a%2==0: return 0
        r=a%8
        return 1 if r in (1,7) else -1
    r=a%p
    if r==0: return 0
    return 1 if pow(r,(p-1)//2,p)==1 else -1

def eichler_symbol(Delta,p):
    dk,f=fundamental_disc_and_conductor(Delta)
    # Standard modified Eichler symbol for quadratic order Delta=dk*f^2:
    # if p divides conductor, symbol is 1; otherwise ordinary Kronecker.
    if f%p==0:
        return 1
    return kronecker_prime(dk,p)

def class_number_order(Delta):
    # Count primitive reduced positive-definite binary quadratic forms [a,b,c]
    # of discriminant Delta. For Delta<-4 this is the proper class number.
    assert Delta<0 and Delta%4 in (0,1)
    # Reduced form bound a <= sqrt(|D|/3).
    Amax=isqrt(abs(Delta)//3)+2
    h=0
    reps=[]
    for aa in range(1,Amax+1):
        # -a < b <= a is a convenient unique convention after boundary rules.
        for bb in range(-aa,aa+1):
            num=bb*bb-Delta
            den=4*aa
            if num%den: continue
            cc=num//den
            if aa>cc: continue
            if gcd(gcd(aa,abs(bb)),cc)!=1: continue
            # Reduced boundary convention: if |b|=a or a=c, require b>=0.
            if (abs(bb)==aa or aa==cc) and bb<0: continue
            h+=1; reps.append((aa,bb,cc))
    return h,reps

def embedding_count(Delta):
    h,_=class_number_order(Delta)
    facD=[]
    facN=[]
    e=h
    for p in DQUAT:
        s=eichler_symbol(Delta,p)
        facD.append((p,s,1-s))
        e*=1-s
    for q in LEVEL:
        s=eichler_symbol(Delta,q)
        facN.append((q,s,1+s))
        e*=1+s
    return h,e,facD,facN

def admissible_basic(Delta):
    # Necessary field embedding condition: ramified quaternion primes must not split,
    # and optimality forces conductor coprime to D.
    dk,f=fundamental_disc_and_conductor(Delta)
    if gcd(f,6)!=1: return False
    return all(kronecker_prime(dk,p)!=1 for p in DQUAT)

rows=[]
for n in range(3,args.disc_max+1):
    Delta=-n
    if Delta%4 not in (0,1): continue
    try:
        dk,f=fundamental_disc_and_conductor(Delta)
    except ValueError:
        continue
    if not admissible_basic(Delta):
        continue
    h,e,fd,fn=embedding_count(Delta)
    if e==0: continue
    rows.append((Delta,dk,f,h,e,fd,fn))

# Primary rank: exact target count first, then small |Delta| / class number.
rows.sort(key=lambda r:(r[4]!=args.embedding_count, abs(r[4]-args.embedding_count), abs(r[0]), r[3]))

print(f"EMBED|stage=start|D=6|N=79|disc_max={args.disc_max}|target={args.embedding_count}")
matches=[r for r in rows if r[4]==args.embedding_count]
print(f"EMBED|stage=summary|admissible_nonzero={len(rows)}|target_matches={len(matches)}")

shown=rows if args.show_all else rows[:args.top]
for i,(Delta,dk,f,h,e,fd,fn) in enumerate(shown,1):
    ds=",".join(f"{p}:{s}:{factor}" for p,s,factor in fd)
    ns=",".join(f"{p}:{s}:{factor}" for p,s,factor in fn)
    flag="TARGET" if e==args.embedding_count else ""
    print(f"EMBED|rank={i}|Delta={Delta}|field_disc={dk}|conductor={f}|h={h}|e={e}|Dlocal={ds}|Nlocal={ns}|{flag}")

print("")
print("TARGET MATCHES")
for Delta,dk,f,h,e,fd,fn in matches:
    ds=",".join(f"{p}:{s}:{factor}" for p,s,factor in fd)
    ns=",".join(f"{p}:{s}:{factor}" for p,s,factor in fn)
    print(f"TARGET|Delta={Delta}|field_disc={dk}|conductor={f}|h={h}|e={e}|Dlocal={ds}|Nlocal={ns}")

print("")
print("NOTE|e counts optimal embedding classes on the Shimura curve before quotienting by extra Atkin-Lehner identifications.")
print("NOTE|Therefore e=4 is a strong fingerprint for a four-point CM orbit, not by itself a proof that the orbit is Elkies t=2.")
