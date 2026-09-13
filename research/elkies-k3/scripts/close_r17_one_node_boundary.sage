#!/usr/bin/env sage-python
"""Resolve all retained pole/infinity cases of the fixed one-node experiment."""
import argparse
from hashlib import sha256
import itertools
import json
from pathlib import Path
import resource
import time

from sage.all import PolynomialRing, QQ, matrix, vector

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/"artifacts/generated-results/elkies-k3-r17-one-node-correlated-v1"
TARGET=ROOT/"artifacts/generated-results/elkies-k3-r17-correlated-genus-one-v1/pencils.json"


def coeffs(q):
    return [str(x) for x in q.list()]


def write_new(path,data):
    with path.open("x") as stream:
        json.dump(data,stream,indent=2,sort_keys=True)
        stream.write("\n")


def factor_record(f):
    return [{"coefficients":coeffs(g),"multiplicity":int(e)} for g,e in f.factor()]


def inspect_carrier(q,r,h,Nx,Ny,M,A,B,inverses,R):
    if r is None:
        assert q.degree()<=4
        d=q
        linear=R(1)
    else:
        linear=R.gen()-r
        d,rem=q.quo_rem(linear**2)
        assert not rem
    out={"q":coeffs(q),"d":coeffs(d),"d_factorization":factor_record(d),"matches":[]}
    if d.degree() not in (3,4) or d.gcd(d.derivative()).degree()!=0:
        out["status"]="OUTSIDE_SMOOTH_GENUS_ONE_NORMALIZATION"
        out["split_over_original_field"]=bool(d.is_square())
        return out
    out["status"]="GENUS_ONE_NORMALIZATION"
    K=R.fraction_field()
    x0=R((M*M-Nx)//(2*h*h));x1=h*linear/2
    y0=K(M)/h*(x0-K(Nx)/h**2)-K(Ny)/h**3
    assert y0.denominator()==1
    y0=R(y0);y1=M*linear/2
    assert y0*y0+y1*y1*d==x0**3+3*x0*x1*x1*d+A*x0+B
    assert 2*y0*y1==3*x0*x0*x1+x1**3*d+A*x1
    out.update({"x0":coeffs(x0),"x1":coeffs(x1),"y0":coeffs(y0),"y1":coeffs(y1),
                "branch_disjoint_from_Delta":bool(d.gcd(4*A**3+27*B**2).degree()==0)})
    for j,inverse in enumerate(inverses):
        v=inverse*vector(QQ,[d[i] for i in range(5)])
        if not all(v[a]*v[b+1]==v[a+1]*v[b] for a,b in itertools.combinations(range(4),2)):
            continue
        if not v[0]:
            out["matches"].append({"target":j,"status":"SINGULAR_TARGET_AT_INFINITY"})
        else:
            lam=v[1]/v[0]
            assert v==v[0]*vector(QQ,[1,lam,lam**2,lam**3,lam**4])
            out["matches"].append({"target":j,"lambda":str(lam),"constant_ratio":str(v[0]),
                                   "same_quadratic_extension":bool(v[0].is_square())})
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir",type=Path,default=DEFAULT)
    args=parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(60,65))
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    started=time.process_time()
    packet=json.loads((args.input_dir/"input.json").read_text())
    assert sha256(TARGET.read_bytes()).hexdigest()==packet["target_pencils_sha256"]
    R=PolynomialRing(QQ,"t");t=R.gen()
    S=PolynomialRing(QQ,names=("a","b"));a,b=S.gens()
    U=PolynomialRing(S,"t");u=U.gen()
    L=PolynomialRing(QQ,"z");z=L.gen()
    A,B=R(packet["generic_source"]["A"]),R(packet["generic_source"]["B"])
    inverses=[matrix(QQ,p["branch_matrix"]).inverse() for p in json.loads(TARGET.read_text())]
    records=[]
    for trace in packet["traces"]:
        index=trace["index"]
        prior=json.loads((args.input_dir/"traces"/f"{index:03d}.json").read_text())
        deferred=[r for r in prior["records"] if r["status"]=="DEFERRED_TRACE_POLE_OR_INFINITY"]
        h,Nx,Ny,M0=[R(trace[key]) for key in ("h","Nx","Ny","M0")]
        hu,nxu,nyu,mu=U(h),U(Nx),U(Ny),U(M0)+(a+b*u)*U(h)**2
        qsym,rem=(mu**4-6*mu*mu*nxu-8*mu*nyu-3*nxu*nxu-4*U(A)*hu**4).quo_rem(hu**6)
        assert not rem and qsym.degree()<=6
        for item in deferred:
            r=None if item["r"] is None else QQ(item["r"])
            row={"trace_index":index,"r":item["r"],"carriers":[]}
            if r is None:
                equation=L(qsym[6](0,z))
                assert all(exponents[0]==0 for exponents in qsym[6].dict())
                roots=equation.roots(QQ,multiplicities=False)
                pairs=[]
                for bv in roots:
                    other=L(qsym[5](z,bv))
                    assert other.degree()==1
                    pairs.append((-other[0]/other[1],bv))
                row.update({"equation_variable":"b","equation":coeffs(equation),"factorization":factor_record(equation)})
            else:
                assert h(r)==0 and (4*A**3+27*B**2)(r)!=0
                at_r=qsym(r)
                equation=L(at_r(z,0))
                # q(r) depends only on a+b*r; verify this coefficientwise.
                assert at_r==S(sum(equation[i]*(a+b*r)**i for i in range(equation.degree()+1)))
                pairs=[]
                for lv in equation.roots(QQ,multiplicities=False):
                    other=L(qsym.derivative()(r)(lv-z*r,z))
                    assert other.degree()==1
                    bv=-other[0]/other[1]
                    pairs.append((lv-bv*r,bv))
                row.update({"equation_variable":"a+b*r","equation":coeffs(equation),"factorization":factor_record(equation)})
            for av,bv in pairs:
                M=M0+(av+bv*t)*h*h
                q=R([c(av,bv) for c in qsym])
                result=inspect_carrier(q,r,h,Nx,Ny,M,A,B,inverses,R)
                result.update({"a":str(av),"b":str(bv)})
                row["carriers"].append(result)
            row["status"]="NO_NODE_IN_THIS_CHART" if not pairs else "EXACT_BOUNDARY_CARRIERS"
            records.append(row)
    positive=[(r["trace_index"],r["r"],c["matches"]) for r in records for c in r["carriers"] if c["status"]=="GENUS_ONE_NORMALIZATION"]
    result={"schema":"r17-one-node-boundary-v1","input_sha256":sha256((args.input_dir/"input.json").read_bytes()).hexdigest(),
            "records":records,"resolved_cases":len(records),"genus_one_carriers":positive,
            "cpu_seconds":time.process_time()-started,
            "boundary":"These equations close only the70 previously deferred cases in the same fixed finite fibre panel and finite regular slope plane. They do not enumerate all rational node positions."}
    write_new(args.input_dir/"boundary.json",result)
    print(json.dumps({"resolved_cases":len(records),"genus_one_carriers":positive,"cpu_seconds":result["cpu_seconds"]}),flush=True)


if __name__=="__main__":
    main()
