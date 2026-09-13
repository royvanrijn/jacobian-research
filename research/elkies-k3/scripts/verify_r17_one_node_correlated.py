#!/usr/bin/env python3
"""Portable independent replay of the frozen one-node halving experiment.

The verifier proves the negative rational-root decisions by retained finite
projective witnesses, not by trusting Sage's factorization flags. --export
creates those witnesses from the retained exact polynomials, without searching
for rational points or repeating the constructor.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import itertools
import json
from math import gcd, isqrt, lcm
from pathlib import Path
import time
import zipfile

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/"artifacts/generated-results/elkies-k3-r17-one-node-correlated-v1"
SPEC=importlib.util.spec_from_file_location("pencil_fraction_replay",ROOT/"elkies-k3/scripts/verify_r17_correlated_genus_one_pencils.py")
V=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V)
require=V.require


def evaluate(a,x):
    result=Q(0)
    for c in reversed(a): result=result*x+c
    return result


def divide_exact(a,b):
    a=list(a); result=[Q(0)]*max(0,len(a)-len(b)+1)
    while len(a)>=len(b):
        n,c=len(a)-len(b),a[-1]/b[-1]
        result[n]=c
        for j,x in enumerate(b): a[n+j]-=c*x
        a=V.trim(a)
    require(not a,"nonexact polynomial division")
    return V.trim(result)


def primitive(a):
    scale=lcm(*(x.denominator for x in a))
    integers=[int(x*scale) for x in a]
    content=gcd(*integers)
    require(content,"zero root-test polynomial")
    return [x//content for x in integers]


def has_projective_root(coefficients,p):
    require(p>=2 and all(p%d for d in range(2,isqrt(p)+1)),"nonprime root witness")
    if coefficients[-1]%p==0: return True
    for x in range(p):
        y=0
        for c in reversed(coefficients): y=(y*x+c)%p
        if not y: return True
    return False


def root_witness(a):
    integer=primitive(a)
    for p in range(2,1000):
        if all(p%d for d in range(2,isqrt(p)+1)) and not has_projective_root(integer,p):
            return p
    raise ValueError("No root witness below1000; retain UNKNOWN instead of a no-root conclusion")


def special(coeffs,r,weight,p):
    if r is None: return int(coeffs[weight])%p if len(coeffs)>weight else 0
    out=0
    for c in reversed(coeffs): out=(out*r+int(c))%p
    return out


def sum_ff(P,Q0,a,p):
    x,y=P;u,v=Q0
    if x==u:
        if (y+v)%p==0:return None
        m=(3*x*x+a)*pow(2*y,-1,p)%p
    else:m=(v-y)*pow((u-x)%p,-1,p)%p
    xx=(m*m-x-u)%p
    return xx,(m*(x-xx)-y)%p


def finite_tables(packet):
    source=packet["generic_source"]
    A,B=V.poly(source["A"]),V.poly(source["B"])
    basis=[(V.poly(r["x"]),V.poly(r["y"])) for r in source["basis"]]
    require(all(x.denominator==1 for c in [A,B]+[c for point in basis for c in point] for x in c),"nonintegral finite-model input")
    output=[]
    for p in packet["limits"]["primes"]:
        values=[]
        # Cache the four-power tuples for the tangency quartic criterion.
        powers=[(m,pow(m,2,p),pow(m,4,p)) for m in range(p)]
        for r in list(range(p))+[None]:
            a,b=special(A,r,8,p),special(B,r,12,p)
            if (4*a**3+27*b*b)%p==0:
                values.append(None);continue
            points=[(special(x,r,4,p),special(y,r,6,p)) for x,y in basis]
            require(all((y*y-x**3-a*x-b)%p==0 for x,y in points),"finite generic section identity")
            mask=0
            for k,row in enumerate(packet["traces"]):
                i,j=[i for i,n in enumerate(row["word"]) if n]
                U=points[j]
                if row["word"][j]<0:U=(U[0],-U[1]%p)
                T=sum_ff(points[i],U,a,p)
                if T is None:continue  # O always has the rational half O.
                x,y=T
                if not any((m4-6*x*m2-8*y*m-3*x*x-4*a)%p==0 for m,m2,m4 in powers):mask|=1<<k
            values.append(str(mask))
        output.append({"prime":p,"rejected_trace_masks_by_projective_residue":values})
    return output


def read_traces(path):
    manifest=json.loads((path/"trace-records-manifest.json").read_text())
    archive=path/"trace-records.zip"
    require(sha256(archive.read_bytes()).hexdigest()==manifest["archive_sha256"],"trace archive hash")
    result=[]
    with zipfile.ZipFile(archive) as z:
        require(z.namelist()==[r["path"] for r in manifest["members"]],"trace archive coverage")
        for row in manifest["members"]:
            raw=z.read(row["path"])
            require(len(raw)==row["bytes"] and sha256(raw).hexdigest()==row["sha256"],"trace member hash")
            result.append(json.loads(raw))
    require(len(result)==67,"trace record count")
    return result


def q_channels(trace,A):
    h,Nx,Ny,M=[V.poly(trace[k]) for k in ("h","Nx","Ny","M0")]
    q0=divide_exact(V.sub(V.sub(V.sub(V.sub(V.power(M,4),V.scale(V.mul(V.power(M,2),Nx),6)),V.scale(V.mul(M,Ny),8)),V.scale(V.power(Nx,2),3)),V.scale(V.mul(A,V.power(h,4)),4)),V.power(h,6))
    q1=divide_exact(V.sub(V.sub(V.scale(V.power(M,3),4),V.scale(V.mul(M,Nx),12)),V.scale(Ny,8)),V.power(h,4))
    q2=V.scale(divide_exact(V.sub(V.power(M,2),Nx),V.power(h,2)),6)
    return [q0,q1,q2,V.scale(M,4),V.power(h,2)]


def node_polynomial(channels,a,b):
    q=[]
    for j,f in enumerate(channels):q=V.add(q,V.mul(f,V.power([a,b],j)))
    return q


def identify_split_branches(trace,carrier,basis):
    """Identify both rational branches by exact generic group-law identities."""
    h,Nx,Ny,M0=[V.poly(trace[k]) for k in ("h","Nx","Ny","M0")]
    M=V.add(M0,V.mul([Q(carrier["a"]),Q(carrier["b"])],V.power(h,2)))
    q=V.poly(carrier["q"])
    r=Q(-21,4)
    sqrt_q=V.scale(V.mul([-r,Q(1)],[Q(348480,94393),Q(430623,94393),Q(1)]),
                   isqrt(442567226524897748498605056))
    require(V.power(sqrt_q,2)==q,"split sextic square root")
    x0=divide_exact(V.sub(V.power(M,2),Nx),V.scale(V.power(h,2),2))
    y0=divide_exact(V.sub(V.mul(M,V.sub(V.mul(x0,V.power(h,2)),Nx)),Ny),V.power(h,3))
    lifted=[(V.add(x0,V.scale(V.mul(h,sqrt_q),sign*Q(1,2))),
             V.add(y0,V.scale(V.mul(M,sqrt_q),sign*Q(1,2)))) for sign in (1,-1)]
    expected=[]
    for i,si,j,sj in ((0,1,2,1),(0,-1,15,-1)):
        x,y=basis[i];u,v=basis[j];y=V.scale(y,si);v=V.scale(v,sj)
        D,N=V.sub(u,x),V.sub(v,y)
        xx=divide_exact(V.sub(V.power(N,2),V.mul(V.add(x,u),V.power(D,2))),V.power(D,2))
        yy=divide_exact(V.sub(V.mul(N,V.sub(x,xx)),V.mul(y,D)),D)
        expected.append((xx,yy))
    require(sorted(lifted)==sorted(expected),"split inherited section identities")
    return ["P1+P3","-P1-P16"]


def verify(path,export=False):
    packet=json.loads((path/"input.json").read_text())
    result=json.loads((path/"result.json").read_text())
    boundary=json.loads((path/"boundary.json").read_text())
    raw_tables=json.loads((path/"finite-halving-tables.json").read_text())
    traces=read_traces(path)
    require(packet["schema"]=="r17-one-node-correlated-input-v1","input schema")
    require(result["input_sha256"]==boundary["input_sha256"]==sha256((path/"input.json").read_bytes()).hexdigest(),"input binding")
    source=packet["generic_source"]
    V.source_projection(source)
    A,B=V.poly(source["A"]),V.poly(source["B"])
    basis=[(V.poly(r["x"]),V.poly(r["y"])) for r in source["basis"]]
    G=[[Q(x) for x in r] for r in source["gram"]]
    words=[]
    for i,j in itertools.combinations(range(17),2):
        for sign in (1,-1):
            if G[i][i]+G[j][j]+2*sign*G[i][j]==6:
                w=[0]*17;w[i],w[j]=1,sign;words.append(w);break
    words.sort()
    require(len(words)==67 and [t["word"] for t in packet["traces"]]==words,"trace-word selection")
    channels=[]
    for index,(word,trace) in enumerate(zip(words,packet["traces"])):
        require(trace["index"]==index,"trace order")
        h,Nx,Ny,M=[V.poly(trace[k]) for k in ("h","Nx","Ny","M0")]
        require(len(h)==2 and h[-1]==1 and len(M)<=2,"regular norm-six frame")
        require(len(V.gcd_poly(Nx,h))==1,"trace pole numerator")
        i,j=[i for i,n in enumerate(word) if n]
        x,y=basis[i];u,v=basis[j];v=V.scale(v,word[j])
        D,N=V.sub(u,x),V.sub(v,y)
        xx=V.sub(V.power(N,2),V.mul(V.add(x,u),V.power(D,2)))
        yy=V.sub(V.mul(N,V.sub(V.mul(x,V.power(D,2)),xx)),V.mul(y,V.power(D,3)))
        require(V.mul(Nx,V.power(D,2))==V.mul(xx,V.power(h,2)),"trace x identity")
        require(V.mul(Ny,V.power(D,3))==V.mul(yy,V.power(h,3)),"trace y identity")
        require(not V.remainder(V.add(V.mul(M,Nx),Ny),V.power(h,2)),"trace RR congruence")
        channels.append(q_channels(trace,A))
    fibres=sorted({Q(a,b) for a in range(-32,33) for b in range(1,33) if gcd(a,b)==1})
    poles=sorted({-Q(r["h"][0])/Q(r["h"][1]) for r in packet["traces"]})
    panel=[str(x) for x in sorted(set(fibres)|set(poles))]+[None]
    require(packet["height_bound"]==32 and packet["height_box_fibres"]==list(map(str,fibres)),"fibre box")
    require(packet["automatic_pole_controls"]==list(map(str,poles)) and packet["finite_fibres"]==panel[:-1] and packet["include_infinity"],"projective fibre panel")
    require(packet["limits"]["primes"]==[101,103,107,109,113,127,131,137],"finite prime roster")
    tables=finite_tables(packet)
    require(tables==raw_tables,"independent tangency/finite-doubling table comparison")
    require(sha256((path/"finite-halving-tables.json").read_bytes()).hexdigest()==result["finite_halving_tables_sha256"],"finite table binding")
    witnesses={} if export else json.loads((path/"root-witnesses.json").read_text())["witnesses"]
    used=set(); counts=Counter();deferred=set();positive=[]
    def exclude(key,poly):
        if export:witnesses[key]=root_witness(poly)
        require(key in witnesses and not has_projective_root(primitive(poly),witnesses[key]),"rational root exclusion witness")
        used.add(key)
    for index,record in enumerate(traces):
        trace=packet["traces"][index]
        require(record["trace_index"]==index and record["word"]==trace["word"],"trace-record attachment")
        require([r["r"] for r in record["records"]]==panel,"per-trace fibre coverage")
        h,Nx,Ny=[V.poly(trace[k]) for k in ("h","Nx","Ny")]
        local=Counter()
        for row in record["records"]:
            status=row["status"];local[status]+=1;counts[status]+=1
            rr=row["r"];r=None if rr is None else Q(rr)
            key=f"finite:{index}:{rr}"
            if status=="NO_RATIONAL_HALF_BY_GOOD_REDUCTION":
                tab=next(t for t in tables if t["prime"]==row["prime"]);p=tab["prime"]
                residue=p if r is None or r.denominator%p==0 else r.numerator*pow(r.denominator,-1,p)%p
                mask=tab["rejected_trace_masks_by_projective_residue"][residue]
                require(mask is not None and (int(mask)>>index)&1,"finite half exclusion")
            elif status=="DEFERRED_TRACE_POLE_OR_INFINITY":
                require(r is None or not evaluate(h,r),"unjustified deferred chart")
                deferred.add((index,rr))
            else:
                require(r is not None and evaluate(h,r),"regular fibre chart")
                ar,br=evaluate(A,r),evaluate(B,r)
                require(4*ar**3+27*br*br,"smooth rational fibre")
                x,y=evaluate(Nx,r)/evaluate(h,r)**2,evaluate(Ny,r)/evaluate(h,r)**3
                F=[-3*x*x-4*ar,-8*y,-6*x,Q(0),Q(1)]
                require(V.poly(row["halving_slope_quartic"])==F,"tangency quartic identity")
                if status=="EXACT_NO_RATIONAL_HALF":exclude(key,F)
                else:
                    require(status=="RATIONAL_HALF_FOUND" and index==54 and r==Q(-21,4),"unexpected positive case needs review")
                    require(len(row["rational_slopes"])==len(row["carriers"])==1,"positive branch count")
                    slope=Q(row["rational_slopes"][0]);require(not evaluate(F,slope),"rational tangency slope")
                    quotient=divide_exact(F,[-slope,Q(1)])
                    exclude(key+":residual-cubic",quotient)
                    c=row["carriers"][0];xq,yq=map(Q,c["half"])
                    require(xq==(slope*slope-x)/2 and yq==slope*(xq-x)-y,"rational half coordinates")
                    require(yq*yq==xq**3+ar*xq+br,"half lies on curve")
                    tangent=(3*xq*xq+ar)/(2*yq)
                    dx=tangent*tangent-2*xq;dy=tangent*(xq-dx)-yq
                    require((dx,dy)==(x,y),"exact doubling of half")
                    q=node_polynomial(channels[index],Q(c["a"]),Q(c["b"]))
                    require(q==V.poly(c["q"]),"one-node carrier identity")
                    d=divide_exact(q,V.power([-r,Q(1)],2))
                    require(d==V.poly(c["d"]),"one-node normalization identity")
                    # The observed quartic is a rational square; certify this by
                    # its literal monic quadratic factor and rational scalar.
                    f=[Q(348480,94393),Q(430623,94393),Q(1)]
                    constant=Q(442567226524897748498605056)
                    require(d==V.scale(V.power(f,2),constant),"split quartic factor identity")
                    require(constant.denominator==1 and isqrt(constant.numerator)**2==constant.numerator,"split constant squareclass")
                    identities=identify_split_branches(trace,c,basis)
                    positive.append({"trace":index,"r":rr,"normalization":"split over Q(t)","inherited_branches":identities,"new_directions":0})
        require(dict(local)==record["counts"]==result["trace_summaries"][index]["counts"],"per-trace count reconciliation")
    require(dict(counts)==result["totals"] and sum(counts.values())==91187,"total case reconciliation")
    require({(r["trace_index"],r["r"]) for r in boundary["records"]}==deferred and len(deferred)==len(boundary["records"])==70,"boundary coverage")
    for row in boundary["records"]:
        k,rr=row["trace_index"],row["r"]
        if rr is None:
            equation=[]
            for j,f in enumerate(channels[k]):
                equation.append(f[6-j] if 0<=6-j<len(f) else Q(0))
        else:
            r=Q(rr)
            equation=[evaluate(f,r) for f in channels[k]]
        equation=V.trim(equation)
        require(equation==V.poly(row["equation"]),"projective boundary coefficient identity")
        require(row["status"]=="NO_NODE_IN_THIS_CHART" and row["carriers"]==[],"boundary conclusion")
        exclude(f"boundary:{k}:{rr}",equation)
    require(set(witnesses)==used,"root witness coverage")
    if export:
        with (path/"root-witnesses.json").open("x") as f:
            json.dump({"schema":"r17-one-node-projective-root-witnesses-v1","witnesses":witnesses},f,indent=2,sort_keys=True);f.write("\n")
    return {"status":"PASS_FIXED_ONE_NODE_PANEL_NO_NONSPLIT_CARRIER","traces":67,"projective_fibres":1361,"cases":91187,
            "finite_halving_rejections":counts["NO_RATIONAL_HALF_BY_GOOD_REDUCTION"],
            "exact_regular_root_exclusions":counts["EXACT_NO_RATIONAL_HALF"],"boundary_exclusions":70,
            "rational_root_witnesses":len(used),"split_positive_control":positive,
            "scope":"The frozen fibre panel supplies no nonsplit one-node genus-one cover for these 67 traces. This is not an exhaustion of rational node positions or of quadratic covers."}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir",type=Path,default=DEFAULT)
    parser.add_argument("--export",action="store_true")
    parser.add_argument("--output",type=Path)
    args=parser.parse_args();started=time.monotonic()
    result=verify(args.input_dir,args.export)
    result["elapsed_seconds"]=round(time.monotonic()-started,3)
    if args.output:
        with args.output.open("x") as f:json.dump(result,f,indent=2,sort_keys=True);f.write("\n")
    print(json.dumps(result,sort_keys=True))


if __name__=="__main__":main()
