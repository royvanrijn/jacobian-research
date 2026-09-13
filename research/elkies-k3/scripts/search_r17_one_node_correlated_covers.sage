#!/usr/bin/env sage-python
"""Bounded generic-only halving construction of one-node quadratic carriers.

The finite fibre panel is not an exhaustive rational-point calculation on the
trace-halving curves. Boundary cases and singular fibres remain explicit.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import importlib.util
from importlib.machinery import SourceFileLoader
import itertools
import json
from math import gcd
from pathlib import Path
import resource
import time

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, matrix, vector

ROOT = Path(__file__).resolve().parents[2]
PRIOR = ROOT / "artifacts/generated-results/elkies-k3-r17-correlated-genus-one-v1"
DEFAULT = ROOT / "artifacts/generated-results/elkies-k3-r17-one-node-correlated-v1"
PRIMES = [101,103,107,109,113,127,131,137]


def write_new(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        json.dump(data,stream,indent=2,sort_keys=True)
        stream.write("\n")


def coeffs(p):
    return [str(x) for x in p.list()]


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def chord_module():
    path = ROOT / "elkies-k3/scripts/construct_elkies_2026_bisections.sage"
    loader = SourceFileLoader("one_node_chord", str(path))
    spec = importlib.util.spec_from_loader(loader.name,loader)
    out = importlib.util.module_from_spec(spec)
    loader.exec_module(out)
    return out


def prepare(output):
    source = json.loads((PRIOR/"input.json").read_text())
    G = matrix(QQ,source["gram"])
    words = []
    for i,j in itertools.combinations(range(17),2):
        for sign in (1,-1):
            if G[i,i]+G[j,j]+2*sign*G[i,j] == 6:
                w=[0]*17
                w[i],w[j]=1,sign
                words.append(w)
                break
    words.sort()
    assert len(words) == 67
    R=PolynomialRing(QQ,"t")
    K=R.fraction_field()
    A,B=R(source["A"]),R(source["B"])
    E=EllipticCurve(K,[A,B])
    basis=[E(K(R(p["x"])),K(R(p["y"]))) for p in source["basis"]]
    chord=chord_module()
    traces=[]
    for index,word in enumerate(words):
        T=sum((int(n)*p for n,p in zip(word,basis) if n),E(0))
        frame=chord.trace_chord_frame(T[0],T[1],R)
        traces.append({"index":index,"word":word,
                       **{k:coeffs(frame[k]) for k in ("h","Nx","Ny","M0")}})
    fibres=sorted({Fraction(a,b) for a in range(-32,33) for b in range(1,33)
                   if gcd(a,b)==1})
    # Pole fibres are automatic controls, selected from generic traces only.
    poles=sorted({Fraction(str(-QQ(r["h"][0])/QQ(r["h"][1]))) for r in traces if len(r["h"])==2})
    panel=sorted(set(fibres)|set(poles))
    packet={"schema":"r17-one-node-correlated-input-v1",
            "generic_source":source,"prior_input_sha256":digest(PRIOR/"input.json"),
            "target_pencils_sha256":digest(PRIOR/"pencils.json"),
            "traces":traces,"height_bound":32,
            "height_box_fibres":[str(r) for r in fibres],
            "automatic_pole_controls":[str(r) for r in poles],
            "finite_fibres":[str(r) for r in panel],"include_infinity":True,
            "limits":{"cpu_seconds":120,"address_space_gib":4,"primes":PRIMES},
            "selection":"All67 parity-distinct height-six two-term generic-basis words, rational projective height<=32, and automatic trace-pole controls; no exceptional-point data.",
            "positive_gates":["rational half of trace at smooth fibre in regular chord chart",
                              "squarefree quartic normalization and exact lifted section",
                              "same rational squareclass as a smooth target-pencil member",
                              "heights16 and20 or other exact independence proof",
                              "covering base with certified infinitely many rational points"]}
    write_new(output/"input.json",packet)
    print(json.dumps({"traces":67,"box_fibres":len(fibres),"automatic_pole_controls":len(poles),"finite_panel":len(panel)}),flush=True)


def ff_add(P,Q,a,p):
    if P is None: return Q
    if Q is None: return P
    x,y=P; z,w=Q
    if x==z:
        if (y+w)%p==0: return None
        m=(3*x*x+a)*pow(2*y,-1,p)%p
    else:
        m=(w-y)*pow((z-x)%p,-1,p)%p
    u=(m*m-x-z)%p
    return u,(m*(x-u)-y)%p


def evaluate(coefficients,t,p,weight):
    if t is None:
        return coefficients[weight]%p if weight<len(coefficients) else 0
    out=0
    for c in reversed(coefficients): out=(out*t+c)%p
    return out


def screen(packet):
    source=packet["generic_source"]
    def integers(values):
        out=[Fraction(x) for x in values]
        assert all(x.denominator==1 for x in out)
        return [int(x) for x in out]
    A,B=integers(source["A"]),integers(source["B"])
    basis=[(integers(r["x"]),integers(r["y"])) for r in source["basis"]]
    words=[r["word"] for r in packet["traces"]]
    tables=[]
    for p in packet["limits"]["primes"]:
        roots={}
        for y in range(p): roots.setdefault(y*y%p,[]).append(y)
        row=[]
        for value in list(range(p))+[None]:
            a,b=evaluate(A,value,p,8),evaluate(B,value,p,12)
            if (4*a*a*a+27*b*b)%p==0:
                row.append(None)
                continue
            doubles={None}
            for x in range(p):
                for y in roots.get((x*x*x+a*x+b)%p,[]):
                    doubles.add(ff_add((x,y),(x,y),a,p))
            points=[(evaluate(x,value,p,4),evaluate(y,value,p,6)) for x,y in basis]
            for x,y in points: assert (y*y-x*x*x-a*x-b)%p==0
            bits=0
            for k,word in enumerate(words):
                indices=[i for i,n in enumerate(word) if n]
                i,j=indices
                P,Q=points[i],points[j]
                if word[j]<0: Q=(Q[0],-Q[1]%p)
                if ff_add(P,Q,a,p) not in doubles: bits|=1<<k
            row.append(str(bits))
        tables.append({"prime":p,"rejected_trace_masks_by_projective_residue":row})
    return tables


def run(output):
    resource.setrlimit(resource.RLIMIT_CPU,(120,125))
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    started=time.process_time()
    packet=json.loads((output/"input.json").read_text())
    assert digest(PRIOR/"input.json")==packet["prior_input_sha256"]
    assert digest(PRIOR/"pencils.json")==packet["target_pencils_sha256"]
    tables=screen(packet)
    write_new(output/"finite-halving-tables.json",tables)
    print(json.dumps({"stage":"finite-halving-tables","cpu_seconds":time.process_time()-started}),flush=True)
    R=PolynomialRing(QQ,"t"); t=R.gen(); K=R.fraction_field()
    A,B=R(packet["generic_source"]["A"]),R(packet["generic_source"]["B"])
    Delta=4*A**3+27*B**2
    E=EllipticCurve(K,[A,B])
    slopes=PolynomialRing(QQ,"m"); m=slopes.gen()
    target_matrices=[matrix(QQ,r["branch_matrix"]) for r in json.loads((PRIOR/"pencils.json").read_text())]
    inverses=[M.inverse() for M in target_matrices]
    summaries=[]
    for k,trace in enumerate(packet["traces"]):
        h,Nx,Ny,M0=[R(trace[name]) for name in ("h","Nx","Ny","M0")]
        X,Y=K(Nx)/h**2,K(Ny)/h**3
        base_slope=K(M0)/h
        records=[]
        for text_r in packet["finite_fibres"]+[None]:
            r=None if text_r is None else QQ(text_r)
            witness=None
            for tab in tables:
                p=tab["prime"]
                residue=p if r is None or r.denominator()%p==0 else int(r.numerator()%p)*pow(int(r.denominator()%p),-1,p)%p
                bits=tab["rejected_trace_masks_by_projective_residue"][residue]
                if bits is not None and (int(bits)>>k)&1:
                    witness=p
                    break
            if witness is not None:
                records.append({"r":text_r,"status":"NO_RATIONAL_HALF_BY_GOOD_REDUCTION","prime":witness})
                continue
            if r is None or h(r)==0:
                records.append({"r":text_r,"status":"DEFERRED_TRACE_POLE_OR_INFINITY"})
                continue
            if Delta(r)==0:
                records.append({"r":text_r,"status":"OUTSIDE_SMOOTH_FIBRE_DOMAIN"})
                continue
            xt,yt,ar,br=X(r),Y(r),A(r),B(r)
            F=m**4-6*xt*m**2-8*yt*m-3*xt**2-4*ar
            factorization=F.factor()
            rational_slopes=[-f[0]/f[1] for f,e in factorization if f.degree()==1]
            record={"r":text_r,"halving_slope_quartic":coeffs(F),
                    "factorization":[{"coefficients":coeffs(f),"multiplicity":int(e)} for f,e in factorization],
                    "rational_slopes":[str(v) for v in rational_slopes],"carriers":[]}
            record["status"]="EXACT_NO_RATIONAL_HALF" if not rational_slopes else "RATIONAL_HALF_FOUND"
            for slope in rational_slopes:
                curve=EllipticCurve(QQ,[ar,br])
                xq=(slope**2-xt)/2
                yq=slope*(xq-xt)-yt
                Q=curve(xq,yq)
                assert 2*Q==curve(xt,yt)
                denom=4*slope**3-12*xt*slope-8*yt
                if not denom:
                    record["carriers"].append({"status":"DEFERRED_TANGENCY_DERIVATIVE","slope":str(slope)})
                    continue
                slope_derivative=(6*X.derivative()(r)*slope**2+8*Y.derivative()(r)*slope+6*xt*X.derivative()(r)+4*A.derivative()(r))/denom
                L=(slope-base_slope(r))/h(r)
                b=(slope_derivative-base_slope.derivative()(r)-L*h.derivative()(r))/h(r)
                a=L-b*r
                M=M0+(a+b*t)*h**2
                numerator=M**4-6*M*M*Nx-8*M*Ny-3*Nx*Nx-4*A*h**4
                q,rem=numerator.quo_rem(h**6)
                assert not rem and q(r)==0 and q.derivative()(r)==0
                d,rem=q.quo_rem((t-r)**2)
                assert not rem
                carrier={"slope":str(slope),"half":[str(xq),str(yq)],"a":str(a),"b":str(b),
                         "q":coeffs(q),"d":coeffs(d),"matches":[]}
                if d.degree() not in (3,4) or d.gcd(d.derivative()).degree()!=0:
                    carrier["status"]="OUTSIDE_SMOOTH_GENUS_ONE_NORMALIZATION"
                else:
                    carrier["status"]="GENUS_ONE_NORMALIZATION"
                    x0=R((M*M-Nx)//(2*h**2)); x1=h*(t-r)/2
                    y0=K(M)/h*(x0-K(Nx)/h**2)-K(Ny)/h**3
                    assert y0.denominator()==1
                    y0=R(y0); y1=M*(t-r)/2
                    assert y0*y0+y1*y1*d==x0**3+3*x0*x1*x1*d+A*x0+B
                    assert 2*y0*y1==3*x0*x0*x1+x1**3*d+A*x1
                    carrier.update({"x0":coeffs(x0),"x1":coeffs(x1),"y0":coeffs(y0),"y1":coeffs(y1),
                                    "branch_disjoint_from_Delta":bool(d.gcd(Delta).degree()==0)})
                    for j,inverse in enumerate(inverses):
                        v=inverse*vector(QQ,[d[i] for i in range(5)])
                        if not all(v[u]*v[vv+1]==v[u+1]*v[vv] for u,vv in itertools.combinations(range(4),2)):
                            continue
                        if not v[0]:
                            carrier["matches"].append({"target":j,"status":"SINGULAR_TARGET_AT_INFINITY"})
                            continue
                        lam=v[1]/v[0]
                        assert v==v[0]*vector(QQ,[1,lam,lam**2,lam**3,lam**4])
                        carrier["matches"].append({"target":j,"lambda":str(lam),"constant_ratio":str(v[0]),
                                                    "same_quadratic_extension":bool(v[0].is_square())})
                record["carriers"].append(carrier)
            records.append(record)
        counts={}
        for row in records: counts[row["status"]]=counts.get(row["status"],0)+1
        payload={"trace_index":k,"word":trace["word"],"records":records,"counts":counts}
        write_new(output/"traces"/f"{k:03d}.json",payload)
        summaries.append({"trace_index":k,"counts":counts})
        if k%10==0: print(json.dumps({"completed_traces":k+1,"counts":counts,"cpu_seconds":time.process_time()-started}),flush=True)
    totals={}
    for row in summaries:
        for key,value in row["counts"].items(): totals[key]=totals.get(key,0)+value
    result={"schema":"r17-one-node-correlated-result-v1","input_sha256":digest(output/"input.json"),
            "finite_halving_tables_sha256":digest(output/"finite-halving-tables.json"),
            "trace_summaries":summaries,"totals":totals,"cpu_seconds":time.process_time()-started,
            "boundary":"Only the frozen finite fibre panel and67 traces are tested. Pole/infinity charts are deferred. No complete rational-point theorem for halving curves is claimed; no positive two-section result follows without the remaining common-cover, independence and rational-base gates."}
    write_new(output/"result.json",result)
    print(json.dumps({"completed_traces":len(summaries),"totals":totals,"cpu_seconds":result["cpu_seconds"]}),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode",choices=("prepare","run"))
    parser.add_argument("--output",type=Path,default=DEFAULT)
    args=parser.parse_args()
    (prepare if args.mode=="prepare" else run)(args.output)
