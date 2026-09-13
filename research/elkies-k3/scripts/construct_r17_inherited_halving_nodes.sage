#!/usr/bin/env sage-python
"""Construct node positions from generic intersections T=2Q.

This is a fixed finite set of generic words, not an enlarged rational fibre
box or a rational-point exhaustion of the trace-halving curves. It reads no
exceptional points. Every trace is checkpointed before the next is attempted.
"""
import argparse
from hashlib import sha256
import importlib.util
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import resource
import time

from sage.all import EllipticCurve, PolynomialRing, QQ, matrix

ROOT=Path(__file__).resolve().parents[2]
PRIOR=ROOT/"artifacts/generated-results/elkies-k3-r17-one-node-correlated-v1"
PENCILS=ROOT/"artifacts/generated-results/elkies-k3-r17-correlated-genus-one-v1/pencils.json"
DEFAULT=ROOT/"artifacts/generated-results/elkies-k3-r17-inherited-halving-nodes-v1"


def write_new(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("x") as f:json.dump(data,f,indent=2,sort_keys=True);f.write("\n")


def coefficients(f):return [str(x) for x in f.list()]


def load_boundary():
    loader=SourceFileLoader("inherited_node_inspector",str(ROOT/"elkies-k3/scripts/close_r17_one_node_boundary.sage"))
    spec=importlib.util.spec_from_loader(loader.name,loader)
    module=importlib.util.module_from_spec(spec);loader.exec_module(module)
    return module


def prepare(output):
    prior=json.loads((PRIOR/"input.json").read_text())
    source=prior["generic_source"]
    assert all(QQ(source["gram"][i][i])==4 for i in range(17))
    packet={"schema":"r17-inherited-halving-node-input-v1",
            "generic_source":source,"traces":prior["traces"],
            "probe_words":[{"index":i,"sign":s} for i in range(17) for s in (1,-1)],
            "source_input_sha256":sha256((PRIOR/"input.json").read_bytes()).hexdigest(),
            "target_pencils_sha256":sha256(PENCILS.read_bytes()).hexdigest(),
            "selection":"All67 frozen norm-six traces and all34 signed generic published basis sections Q. Take every rational finite zero of the pole polynomial of T-2Q. Node positions are consequences of generic group-law identities, not desired specialized points.",
            "limits":{"trace_probe_pairs":2278,"cpu_seconds":120,"address_space_gib":4},
            "scope":"The finite word bank does not exhaust rational halves or node positions; infinity is already closed for these traces by the retained fixed-panel proof."}
    write_new(output/"input.json",packet)
    print(json.dumps({"prepared_pairs":2278}),flush=True)


def run(output):
    resource.setrlimit(resource.RLIMIT_CPU,(120,125))
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    started=time.process_time()
    packet=json.loads((output/"input.json").read_text())
    assert sha256(PENCILS.read_bytes()).hexdigest()==packet["target_pencils_sha256"]
    R=PolynomialRing(QQ,"t");t=R.gen();K=R.fraction_field()
    source=packet["generic_source"]
    A,B=R(source["A"]),R(source["B"]);E=EllipticCurve(K,[A,B])
    basis=[E(K(R(p["x"])),K(R(p["y"]))) for p in source["basis"]]
    doubles=[2*P for P in basis]
    inspect=load_boundary().inspect_carrier
    inverses=[matrix(QQ,p["branch_matrix"]).inverse() for p in json.loads(PENCILS.read_text())]
    summaries=[];all_carriers=[];seen=set()
    for trace in packet["traces"]:
        index=trace["index"]
        h,Nx,Ny,M0=[R(trace[k]) for k in ("h","Nx","Ny","M0")]
        X,Y=K(Nx)/h**2,K(Ny)/h**3;T=E(X,Y);base=K(M0)/h
        records=[];carriers=[]
        for probe in packet["probe_words"]:
            j,sign=probe["index"],probe["sign"]
            S=T-sign*doubles[j]
            den=S[0].denominator()
            pole=den.sqrt().monic()
            factors=pole.factor()
            roots=sorted(-f[0]/f[1] for f,e in factors if f.degree()==1)
            row={"probe":probe,"pole_polynomial":coefficients(pole),
                 "factorization":[{"coefficients":coefficients(f),"multiplicity":int(e)} for f,e in factors],
                 "rational_roots":[str(r) for r in roots]}
            records.append(row)
            for r in roots:
                if not h(r) or not (4*A**3+27*B**2)(r):
                    row.setdefault("outside_regular_smooth_chart",[]).append(str(r));continue
                xq,yq=R(basis[j][0])(r),sign*R(basis[j][1])(r)
                Er=EllipticCurve(QQ,[A(r),B(r)])
                assert 2*Er(xq,yq)==Er(X(r),Y(r)) and yq
                slope=(3*xq*xq+A(r))/(2*yq)
                identity=(index,str(r),str(slope))
                if identity in seen:continue
                seen.add(identity)
                slope_prime=(6*X.derivative()(r)*slope**2+8*Y.derivative()(r)*slope+6*X(r)*X.derivative()(r)+4*A.derivative()(r))/(4*slope**3-12*X(r)*slope-8*Y(r))
                L=(slope-base(r))/h(r)
                b=(slope_prime-base.derivative()(r)-L*h.derivative()(r))/h(r)
                a=L-b*r;M=M0+(a+b*t)*h*h
                q,rem=(M**4-6*M*M*Nx-8*M*Ny-3*Nx*Nx-4*A*h**4).quo_rem(h**6)
                assert not rem and not q(r) and not q.derivative()(r)
                candidate=inspect(q,r,h,Nx,Ny,M,A,B,inverses,R)
                candidate.update({"trace_index":index,"r":str(r),"probe":probe,
                                  "half":[str(xq),str(yq)],"slope":str(slope),"a":str(a),"b":str(b)})
                carriers.append(candidate)
        write_new(output/"traces"/f"{index:03d}.json",{"trace_index":index,"records":records,"carriers":carriers})
        summary={"trace_index":index,"pairs":len(records),"rational_root_incidence":sum(len(r["rational_roots"]) for r in records),
                 "distinct_carriers":len(carriers),"genus_one_carriers":sum(c["status"]=="GENUS_ONE_NORMALIZATION" for c in carriers)}
        summaries.append(summary);all_carriers.extend(carriers)
        print(json.dumps(summary),flush=True)
    positive=[c for c in all_carriers if c["status"]=="GENUS_ONE_NORMALIZATION"]
    result={"schema":"r17-inherited-halving-nodes-result-v1","input_sha256":sha256((output/"input.json").read_bytes()).hexdigest(),
            "trace_summaries":summaries,"pairs":sum(s["pairs"] for s in summaries),
            "rational_root_incidence":sum(s["rational_root_incidence"] for s in summaries),
            "distinct_carriers":len(all_carriers),"genus_one_carriers":positive,
            "cpu_seconds":time.process_time()-started,
            "scope":"Only rational finite intersections T=2Q for the frozen67-by34 word bank. This is not an exhaustion of rational node positions or quadratic covers."}
    write_new(output/"result.json",result)
    print(json.dumps({k:result[k] for k in ("pairs","rational_root_incidence","distinct_carriers","cpu_seconds")}|{"genus_one_carriers":len(positive)}),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("mode",choices=("prepare","run"));p.add_argument("--output",type=Path,default=DEFAULT)
    a=p.parse_args();prepare(a.output) if a.mode=="prepare" else run(a.output)


if __name__=="__main__":main()
