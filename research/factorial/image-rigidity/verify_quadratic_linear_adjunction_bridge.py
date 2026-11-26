#!/usr/bin/env python3
"""Finite native identities for the private reversible linear-adjunction proof."""
if not __debug__:
    raise SystemExit("Optimized Python is refused.")
from fractions import Fraction as F
from math import factorial
import json

def require(ok, msg):
    if not ok:
        raise ValueError(msg)

def compositions(n,d):
    if d==1:
        yield (n,)
    else:
        for i in range(n+1):
            for rest in compositions(n-i,d-1):
                yield (i,)+rest

def native(n,s):
    out={}
    for e in compositions(n,3+s):
        a,b,c=e[:3]
        # The extra variables have exponents e[3:]; their factorials
        # cancel the corresponding multinomial denominators exactly.
        out[e]=F(factorial(n)*factorial(2*a+b),
                 factorial(a)*factorial(b)*factorial(c))
    return out

def add(p,q):
    out=p.copy()
    for e,c in q.items():
        out[e]=out.get(e,F(0))+c
        if not out[e]:
            del out[e]
    return out

def shift_last(p,c):
    return {e[:-1]+(e[-1]+1,):c*v for e,v in p.items()}

def embed(p):
    return {e+(0,):c for e,c in p.items()}

def polyadd(p,q):
    out=[F(0)]*max(len(p),len(q))
    for i,c in enumerate(p): out[i]+=c
    for i,c in enumerate(q): out[i]+=c
    while len(out)>1 and not out[-1]: out.pop()
    return out

def polymul(p,q):
    out=[F(0)]*(len(p)+len(q)-1)
    for i,c in enumerate(p):
        for j,d in enumerate(q): out[i+j]+=c*d
    while len(out)>1 and not out[-1]: out.pop()
    return out

def fixed_moment(n,C):
    out=[F(0)]*(n+1)
    for (a,b,c),v in native(n,0).items():
        out[a]+=v*C**c
    return out

controls=[]
negative=[]
for s in (1,2):
    for n in range(1,9):
        rhs=add(shift_last(native(n-1,s),n),embed(native(n,s-1)))
        require(rhs==native(n,s),"Native linear-adjunction identity")
        require(native(n,s)[(0,)*(2+s)+(n,)]==factorial(n),"Leading coefficient")
        controls.append((s,n))
        wrong=add(shift_last(native(n-1,s),n+1),embed(native(n,s-1)))
        require(wrong!=native(n,s),"Incorrect recurrence multiplier accepted")
        negative.append((s,n))
P=fixed_moment(2,-4)
Q=fixed_moment(3,-4)
require(P==[10,-4,24] and Q==[-34,24,72,720],"Native initial polynomials")
left=polyadd(polymul([5577,4032,-10980],P),polymul([-232,366],Q))
require(left==[63658],"Initial-four-window literal Bezout identity")
require(polyadd(left,[-1])!=[63658],"Corrupted Bezout constant accepted")
require(sum(fixed_moment(2,-3))==29,"Initial-three-window native value")
result={"status":"PASS","scope":"Finite native polynomial identities only; root-existence and universal equivalences are proved in the written source.","adjunction_identities":len(controls),"literal_initial_Bezout":True,"native_initial_three_value":29,"rejected_corruptions":len(negative)+1}
import argparse,hashlib
from pathlib import Path
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument("--output",type=Path)
args=ap.parse_args()
result["checker_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
encoded=json.dumps(result,indent=2,sort_keys=True)+"\n"
if args.output:
    with args.output.open("x",encoding="utf-8") as f:
        f.write(encoded)
print(encoded,end="")
