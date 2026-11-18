#!/usr/bin/env python3
"""Exact Rouche continuation and an explicit degree-eight Hurwitz braid word.

Default verification uses only rational arithmetic and a retained certificate.
--generate uses mpmath to propose centers; every proposal is checked exactly.
"""
from __future__ import annotations
import argparse
from collections import deque
from fractions import Fraction as F
import gzip,hashlib,json,math
from pathlib import Path
import platform

HERE=Path(__file__).resolve().parent
CERT=HERE/"degree_eight_braid_word_certificate.json.gz"
ZERO=(F(0),F(0))
ONE=(F(1),F(0))
BITS=70
D=1<<BITS
ANCHOR=F(1,1<<40)
SQRT_BITS=80
SQRT_FLOOR=math.isqrt(329*(1<<(2*SQRT_BITS)))
ERE=F(-1009,10584)
EIM=F(94*(2*SQRT_FLOOR+1),10584*(1<<(SQRT_BITS+1)))
E=(ERE,EIM)
EERR=F(94,10584*(1<<(SQRT_BITS+1)))
A=F(-625,21952)
CENTERS=[(F(-1323,10000),F(752,10000)),
         (F(-811,10000),F(405,10000)),
         (F(811,10000),F(-405,10000)),
         (F(1323,10000),F(-752,10000))]
HALF=F(1,500)

def req(ok,message):
    if not ok:raise AssertionError(message)
def add(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[0],-a[1]
def sub(a,b):return a[0]-b[0],a[1]-b[1]
def mul(a,b):return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def scale(a,b):return a[0]*b,a[1]*b
def norm1(a):return abs(a[0])+abs(a[1])
def normsq(a):return a[0]*a[0]+a[1]*a[1]
def lower(a):return max(abs(a[0]),abs(a[1]))
def cpow(a,n):
    out=ONE
    for _ in range(n):out=mul(out,a)
    return out
def dyadic(exp):return F(1<<exp) if exp>=0 else F(1,1<<(-exp))
def root_center(record):return F(record[0],D),F(record[1],D)
def taylor(z,s):
    powers=[ONE]
    for _ in range(8):powers.append(mul(powers[-1],z))
    coeff=[]
    for k in range(9):
        v=scale(powers[8-k],math.comb(8,k))
        if k<=6:v=add(v,scale(powers[6-k],math.comb(6,k)))
        if k<=2:v=add(v,scale(mul(E,powers[2-k]),math.comb(2,k)))
        if k==0:v=add(sub(v,mul(s,z)),(A,F(0)))
        if k==1:v=sub(v,s)
        coeff.append(v)
    return coeff
def residual_bounds(z,s,r,ds):
    coeff=taylor(z,s)
    error0=EERR*norm1(cpow(z,2))
    error1=2*EERR*norm1(z)
    constant=norm1(coeff[0])+error0+norm1(ds)*norm1(z)
    derivative=lower(coeff[1])-error1-norm1(ds)
    higher=sum((norm1(coeff[k])+(EERR if k==2 else 0))*r**k for k in range(2,9))
    return constant+higher,derivative*r
def rouche(z,s,r,ds=ZERO):
    l,rhs=residual_bounds(z,s,r,ds)
    return l<rhs
def allowance(z,s,r):
    coeff=taylor(z,s)
    residual=norm1(coeff[0])+EERR*norm1(cpow(z,2))
    deriv=lower(coeff[1])-2*EERR*norm1(z)
    higher=sum((norm1(coeff[k])+(EERR if k==2 else 0))*r**k for k in range(2,9))
    return (deriv*r-residual-higher)/(norm1(z)+r)
def separated(centers,radii):
    return all(normsq(sub(centers[i],centers[j]))>(radii[i]+radii[j])**2
               for i in range(len(centers)) for j in range(i))
def vertices(center):
    x,y=center;h=HALF
    east=(x+h,y)
    return [ONE,east,(x+h,y+h),(x-h,y+h),(x-h,y-h),(x+h,y-h),east,ONE]
def q_taylor(z):
    # Q=s^4+B*s^2+C with coefficients affine in the exact e.
    b=scale(sub(scale(E,3246152),(F(22375),F(0))),F(1,19765032))
    c=scale(add(scale(E,6143242504),(F(991745875),F(0))),F(-625,3515908409649216))
    bzerr=F(3246152,19765032)*EERR
    cerr=F(625*6143242504,3515908409649216)*EERR
    powers=[ONE]
    for _ in range(4):powers.append(mul(powers[-1],z))
    values=[]
    errors=[]
    for k in range(5):
        v=scale(powers[4-k],math.comb(4,k))
        err=F(0)
        if k<=2:
            v=add(v,scale(mul(b,powers[2-k]),math.comb(2,k)))
            err=bzerr*norm1(powers[2-k])*math.comb(2,k)
        if k==0:v=add(v,c);err+=cerr
        values.append(v);errors.append(err)
    return values,errors
def verify_geometry():
    req(SQRT_FLOOR**2<329*(1<<(2*SQRT_BITS))<(SQRT_FLOOR+1)**2,"sqrt enclosure")
    radius=F(1,1000)
    for center in CENTERS:
        vals,errs=q_taylor(center)
        left=norm1(vals[0])+errs[0]+sum((norm1(vals[k])+errs[k])*radius**k for k in range(2,5))
        right=(lower(vals[1])-errs[1])*radius
        req(left<right,"one branch value per rational disk")
        req(radius<HALF,"branch disk lies inside square")
    req(separated(CENTERS,[2*HALF]*4),"disjoint branch squares")
    # A segment from the base to the east side of its square must avoid every
    # other square. Exact interval clipping detects all intersections.
    def meets_box(a,b,center):
        lo,hi=F(0),F(1)
        for j in range(2):
            delta=b[j]-a[j];bottom=center[j]-HALF;top=center[j]+HALF
            if delta==0:
                if not bottom<=a[j]<=top:return False
            else:
                q0,q1=(bottom-a[j])/delta,(top-a[j])/delta
                if q0>q1:q0,q1=q1,q0
                lo=max(lo,q0);hi=min(hi,q1)
                if lo>hi:return False
        return True
    spokes=[vertices(c)[1] for c in CENTERS]
    for i,spoke in enumerate(spokes):
        req(spoke[0]<1,"base lies right of branch squares")
        for j,center in enumerate(CENTERS):
            if i!=j:req(not meets_box(ONE,spoke,center),"spoke avoids other square")
    slopes=[v[1]/(1-v[0]) for v in spokes]
    req(all(slopes[i]>slopes[i+1] for i in range(3)),
        "counterclockwise spoke order")
    for i,a in enumerate(spokes):
        for b in spokes[:i]:
            aa,bb=sub(a,ONE),sub(b,ONE)
            req(aa[0]*bb[1]-aa[1]*bb[0]!=0,"spokes only meet at base")

def comp(a,b):return tuple(a[b[j]] for j in range(8))
def inv(a):
    r=[0]*8
    for i,v in enumerate(a):r[v]=i
    return tuple(r)
def perm(*cycles):
    p=list(range(8))
    for cycle in cycles:
        for a,b in zip(cycle,cycle[1:]+cycle[:1]):p[a]=b
    return tuple(p)
def prod(values):
    p=tuple(range(8))
    for value in values:p=comp(p,value)
    return p
def cycles(value):
    seen=set();out=[]
    for i in range(8):
        if i not in seen:
            c=[];j=i
            while j not in seen:seen.add(j);c.append(j);j=value[j]
            out.append(c)
    return out
def hurwitz(values,letter):
    j=abs(letter)-1;a,b=values[j:j+2]
    pair=(b,comp(comp(b,a),b)) if letter>0 else (comp(comp(a,b),a),a)
    return values[:j]+pair+values[j+2:]
def braid_to_seed(values):
    knot=(perm((0,1),(2,3)),perm((1,4),(2,3)),perm((0,5),(1,6)),perm((1,2),(6,7)))
    target=prod(knot);actual=prod(values)
    ac=next(c for c in cycles(actual) if len(c)==7)
    tc=next(c for c in cycles(target) if len(c)==7)
    af=next(c[0] for c in cycles(actual) if len(c)==1)
    tf=3
    relabel=[0]*8
    for a,b in zip(ac,tc):relabel[a]=b
    relabel[af]=tf;relabel=tuple(relabel)
    normalized=tuple(comp(comp(relabel,p),inv(relabel)) for p in values)
    req(prod(normalized)==target,"normalized product")
    parents={normalized:(None,0)};queue=deque([normalized])
    while queue and knot not in parents:
        here=queue.popleft()
        for letter in (1,-1,2,-2,3,-3):
            there=hurwitz(here,letter)
            if there not in parents:parents[there]=(here,letter);queue.append(there)
    req(knot in parents,"model tuple reaches knot tuple")
    word=[];here=knot
    while parents[here][0] is not None:
        here0,letter=parents[here];word.append(letter);here=here0
    word.reverse()
    check=normalized
    for letter in word:check=hurwitz(check,letter)
    req(check==knot,"braid word replay")
    return {"sheet_relabeling":list(relabel),"normalized_model_tuple":[list(p) for p in normalized],
            "word":word,"knot_tuple":[list(p) for p in knot]}

def generate():
    import mpmath as mp
    mp.mp.dps=70
    def mpc(v):return mp.mpc(mp.mpf(v[0].numerator)/v[0].denominator,mp.mpf(v[1].numerator)/v[1].denominator)
    ee=mpc(E);aa=mp.mpf(A.numerator)/A.denominator
    def quantize(v):return [int(mp.nint(mp.re(v)*D)),int(mp.nint(mp.im(v)*D))]
    def newton(z,s):
        for _ in range(14):
            value=z**8+z**6+ee*z**2-s*z+aa
            deriv=8*z**7+6*z**5+2*ee*z-s
            delta=value/deriv;z-=delta
            if abs(delta)<mp.mpf("1e-58"):break
        return z
    base=sorted(mp.polyroots([1,0,1,0,0,0,ee,-1,aa],maxsteps=500),key=lambda z:(mp.re(z),mp.im(z)))
    initial=[quantize(z) for z in base]
    root_initial=[root_center(v) for v in initial]
    req(separated(root_initial,[2*ANCHOR]*8),"initial separation")
    req(all(rouche(z,ONE,2*ANCHOR) for z in root_initial),"base root anchors")
    loops=[]
    for loopid,center in enumerate(CENTERS):
        roots=initial;steps=[];verts=vertices(center)
        for edge in range(len(verts)-1):
            start,end=verts[edge:edge+2];direction=sub(end,start);t=F(0)
            while t<1:
                current=add(start,scale(direction,t));zz=[root_center(v) for v in roots]
                exps=[];allow=[]
                for i,z in enumerate(zz):
                    minimum=min(abs(complex(float(z[0]-w[0]),float(z[1]-w[1]))) for j,w in enumerate(zz) if i!=j)
                    exp=math.floor(math.log2(minimum/16))
                    while allowance(z,current,dyadic(exp))<=0:exp-=1
                    exps.append(exp);allow.append(allowance(z,current,dyadic(exp)))
                limit=min(allow)/4;expstep=0
                while norm1(direction)*dyadic(expstep)>limit:expstep-=1
                dt=min(1-t,dyadic(expstep));nextt=t+dt
                endpoint=add(start,scale(direction,nextt));ds=sub(endpoint,current)
                nextroots=[quantize(newton(mpc(z),mpc(endpoint))) for z in zz]
                radii=[dyadic(v) for v in exps]
                req(separated(zz,radii),"tracking disks separated")
                for z,new,r in zip(zz,nextroots,radii):
                    req(rouche(z,current,r,ds),"segment root disk")
                    req(rouche(root_center(new),endpoint,ANCHOR),"endpoint root anchor")
                    req(norm1(sub(root_center(new),z))+ANCHOR<r,"endpoint contained in tracking disk")
                steps.append({"edge":edge,"t":str(nextt),"roots":nextroots,"radii":exps})
                roots=nextroots;t=nextt
        permutation=[]
        for z in map(root_center,roots):
            matches=[j for j,w in enumerate(root_initial) if norm1(sub(z,w))+ANCHOR<2*ANCHOR]
            req(len(matches)==1,"closed-loop label");permutation.append(matches[0])
        req(sorted(map(len,cycles(tuple(permutation))))==[1,1,1,1,2,2],"double transposition")
        loops.append({"steps":steps,"permutation":permutation})
        print("generated loop",loopid,"steps",len(steps),"permutation",permutation,flush=True)
    values=tuple(tuple(l["permutation"]) for l in loops)
    req(sorted(map(len,cycles(prod(values))))==[1,7],"finite product seven-cycle")
    braid=braid_to_seed(values)
    output={"schema_version":1,"model":"quadratic degree8; sqrt(-329)=+i sqrt329",
            "root_center_bits":BITS,"sqrt_bits":SQRT_BITS,"base_roots":initial,
            "branch_square_centers":[[str(x),str(y)] for x,y in CENTERS],
            "branch_square_halfside":str(HALF),"loops":loops,"braid":braid}
    payload=json.dumps(output,separators=(",",":")).encode()
    CERT.write_bytes(gzip.compress(payload,mtime=0))
    print("wrote",CERT,"word",braid["word"],flush=True)

def verify():
    verify_geometry()
    with gzip.open(CERT,"rt",encoding="utf-8") as f:data=json.load(f)
    req(data["schema_version"]==1,"certificate schema")
    req(data["root_center_bits"]==BITS and data["sqrt_bits"]==SQRT_BITS,"precision binding")
    req(data["branch_square_centers"]==[[str(x),str(y)] for x,y in CENTERS]
        and data["branch_square_halfside"]==str(HALF),"path geometry binding")
    initial=data["base_roots"];zz0=list(map(root_center,initial))
    req(len(zz0)==8 and separated(zz0,[2*ANCHOR]*8),"initial root count and separation")
    req(all(rouche(z,ONE,2*ANCHOR) for z in zz0),"initial roots")
    req(all(zz0[i][0]+2*ANCHOR<zz0[i+1][0]-2*ANCHOR for i in range(7)),
        "base labels ordered by real part")
    req(not rouche(zz0[0],ONE,dyadic(-200)),"too-small root-disk rejection")
    req(len(data["loops"])==4,"four loops")
    values=[]
    total=0
    for loopid,loop in enumerate(data["loops"]):
        verts=vertices(CENTERS[loopid]);current=ONE;roots=zz0;edge=0;t=F(0)
        for step in loop["steps"]:
            if t==1:edge+=1;t=F(0)
            req(step["edge"]==edge and edge<7,"ordered polygon edges")
            nextt=F(step["t"]);req(t<nextt<=1,"positive segment progression")
            endpoint=add(verts[edge],scale(sub(verts[edge+1],verts[edge]),nextt))
            ds=sub(endpoint,current);radii=[dyadic(v) for v in step["radii"]]
            req(len(radii)==8 and separated(roots,radii),"disjoint tracking disks")
            newroots=list(map(root_center,step["roots"]));req(len(newroots)==8,"eight endpoint anchors")
            for z,w,r in zip(roots,newroots,radii):
                req(ANCHOR<r and rouche(z,current,r,ds),"uniform segment Rouche inequality")
                req(rouche(w,endpoint,ANCHOR),"endpoint root anchor")
                req(norm1(sub(w,z))+ANCHOR<r,"endpoint-root inclusion")
            roots=newroots;current=endpoint;t=nextt;total+=1
        req(edge==6 and t==1 and current==ONE,"closed full polygon loop")
        permutation=[]
        for z in roots:
            matches=[j for j,w in enumerate(zz0) if norm1(sub(z,w))+ANCHOR<2*ANCHOR]
            req(len(matches)==1,"exact endpoint matching");permutation.append(matches[0])
        req(permutation==loop["permutation"],"recorded branch permutation")
        req(sorted(map(len,cycles(tuple(permutation))))==[1,1,1,1,2,2],"DT branch permutation")
        values.append(tuple(permutation))
    req(sorted(map(len,cycles(prod(values))))==[1,7],"seven-cycle product")
    braid=data["braid"];relabel=tuple(braid["sheet_relabeling"])
    req(sorted(relabel)==list(range(8)),"sheet relabeling")
    current=tuple(comp(comp(relabel,p),inv(relabel)) for p in values)
    req([list(p) for p in current]==braid["normalized_model_tuple"],"normalized tuple")
    braid_states=[[list(v) for v in current]]
    for letter in braid["word"]:
        req(letter in (-3,-2,-1,1,2,3),"braid generator")
        current=hurwitz(current,letter)
        braid_states.append([list(v) for v in current])
    knot=(perm((0,1),(2,3)),perm((1,4),(2,3)),perm((0,5),(1,6)),perm((1,2),(6,7)))
    req(current==knot and braid["knot_tuple"]==[list(v) for v in knot],"literal knot equality")
    # Invalid braid replay must not satisfy the claimed final equality.
    altered=tuple(comp(comp(relabel,p),inv(relabel)) for p in values)
    for letter in braid["word"][:-1]:altered=hurwitz(altered,letter)
    req(altered!=knot,"truncated braid rejection")
    return {"schema_version":1,"kind":"exact_rational_Rouche_and_braid_certificate",
            "python":platform.python_version(),"result":"PASS","model":data["model"],
            "parameter_segments":total,"branch_permutations":[list(v) for v in values],
            "sheet_relabeling":list(relabel),"normalized_model_tuple":braid["normalized_model_tuple"],
            "braid_word":braid["word"],"braid_word_length":len(braid["word"]),
            "braid_states":braid_states,
            "negative_controls":["too-small root disk rejected","truncated braid rejected"],
            "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "certificate_sha256":hashlib.sha256(CERT.read_bytes()).hexdigest(),
            "arithmetic":"exact rational complex arithmetic; Rouche bounds for every point of every segment",
            "scope":"Literal quadratic-field model to knot tuple, with positive meridians on rational polygon paths. Numerical centers are proposals only; all root counts, label transitions, branch disks and the braid word are checked exactly.",
            "independent_replay":False,"formal_verification":False,"external_review":False}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generate",action="store_true")
    parser.add_argument("--write-receipt",type=Path)
    args=parser.parse_args()
    if args.generate:generate()
    report=verify();text=json.dumps(report,indent=2)+"\n"
    if args.write_receipt:args.write_receipt.write_text(text)
    print(text,end="")
if __name__=="__main__":main()
