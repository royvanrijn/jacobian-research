#!/usr/bin/env python3
"""Bare exact replay of complete graded Long cap23 homogeneous tail ranks."""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial,isqrt
import argparse,hashlib,importlib.util,json,platform,tempfile

HERE=Path(__file__).resolve().parent
CAP=23
def load(path,name):
    specification=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module
def basis(weights,m):
    out=[]
    for i,w in enumerate(weights):
        for b in range(CAP+1):
            for c in range(CAP+1-b):
                for d in range(CAP+1-b-c):
                    a=b+4*c+2*d+3*m-w
                    if a>=0 and a+b+c+d<=CAP:
                        out.append((i,(a,b,c,d)))
    return sorted(out)
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path)
    args=parser.parse_args()
    cert=load(HERE/"check_graded_long_cap23.py","graded_column_engine")
    alg=load(HERE/"check_long_graph_curvature.py","independent_polynomial_engine")
    source=HERE.parent/"dixmier-audit/verify_long_pbw.py"
    producer=load(source,"canonical_input")
    images=[{tuple(e):Q(c) for e,c in poly.items()}
            for poly in producer.long_images()]  # Polynomial INPUT only.
    witnesspath=HERE/"graded_long_tail_rigidity_witness.json"
    witness=json.loads(witnesspath.read_text())
    assert [record["m"] for record in witness]==list(range(3,10))
    p=1000003
    assert all(p%d for d in range(2,isqrt(p)+1))
    expected_counts=[[34,54,71,47],[20,34,47,29],[10,20,29,16],
                     [4,10,16,8],[1,4,8,3],[0,1,3,1],[0,0,1,0]]
    expected_ranks=[205,130,75,38,16,5,1]
    summary=[]
    for record,counts,rank in zip(witness,expected_counts,expected_ranks):
        m=record["m"]
        unknowns=basis(cert.WEIGHTS,m)
        assert [sum(i==j for i,e in unknowns) for j in range(4)]==counts
        columns=cert.linear_columns(alg,images,unknowns)
        rows=sorted(set().union(*(set(column) for column in columns)))
        assert record["counts"]==counts
        assert record["unknowns"]==len(unknowns)
        assert record["rows"]==len(rows)
        assert record["nonzeros"]==sum(map(len,columns))
        assert record["prime"]==p
        det=cert.minor(columns,{},rows,record["rank_minor"],p)
        assert len(record["rank_minor"]["rows"])==rank
        if m==3:
            r=images[0]
            generator=cert.power(alg,r,6)
            kernel=[{},{},alg.scale(cert.power(alg,r,5),-6),{}]
            assert kernel==[alg.bidifferential(generator,z,1) for z in images]
            assert kernel[2]
            assert not cert.product(columns,cert.vector(unknowns,kernel))
            assert len(unknowns)==rank+1
        else:
            assert len(unknowns)==rank
            assert record["rank_minor"]["columns"]==list(range(len(unknowns)))
        summary.append({"m":m,"basis_counts":counts,"unknowns":len(unknowns),
                        "rows":len(rows),"rank_Q":rank,
                        "kernel_dimension":len(unknowns)-rank,"minor_det_modp":det})
    assert not basis(cert.WEIGHTS,10)
    # For all m>=10, a>=3m-max(w_i)>=26>23: no family is admitted.
    assert 3*10-max(cert.WEIGHTS)>CAP

    # The target h*R^3 shear only changes F2D through F4 under R1=R2=0.
    # Its other possible contributions are the h^1,h^2 coefficients of R#R.
    r=images[0]
    assert all(e[2]+e[3]==0 for e in r)
    assert not alg.bidifferential(r,r,1)
    assert not alg.bidifferential(r,r,2)
    r3=cert.power(alg,r,3)
    assert [alg.bidifferential(r3,z,1) for z in images]==[
        {},{},alg.scale(cert.power(alg,r,2),-3),{}]
    assert max(map(sum,cert.power(alg,r,2)))==6

    # Independently enumerate the universal Moyal coefficient indices.
    # The two unknown-F_(n-1) Poisson terms are excluded from the RHS.
    def rhs_indices(n):
        result=[]
        for k in range(1,n+1,2):
            for a in range(n-k+1):
                b=n-k-a
                if k==1 and (a==n-1 or b==n-1):
                    continue
                result.append((a,b,k,Q(1,2**(k-1)*factorial(k))))
        return result
    assert rhs_indices(4)==[
        (1,2,1,Q(1)),(2,1,1,Q(1)),
        (0,1,3,Q(1,24)),(1,0,3,Q(1,24))]
    assert rhs_indices(5)==[
        (1,3,1,Q(1)),(2,2,1,Q(1)),(3,1,1,Q(1)),
        (0,2,3,Q(1,24)),(1,1,3,Q(1,24)),(2,0,3,Q(1,24)),
        (0,0,5,Q(1,1920))]
    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    receipt={
        "claim":"Complete graded Long cap23 tail: one F3 kernel, injective F4..F9, no coefficients beyond F9",
        "status":"PASS",
        "environment":{"python":platform.python_version(),
                       "arithmetic":"Python integers/Fraction; no symbolic or nullspace library"},
        "source_order":["X","Q","P","Z"],"source_weights":[-1,1,4,2],
        "image_order":list(cert.LABELS),"image_weights":list(cert.WEIGHTS),
        "h_weight":3,"cap":CAP,"prime":p,"homogeneous_ranks":summary,
        "F3_kernel":"G=R^6; variation=(0,0,-6R^5,0)",
        "all_m_ge10_empty":True,
        "u_target_normalization":{
            "formula":"D -> D-3u h^2 R^2; R,T,S fixed",
            "hypotheses":"R1=R2=0, R0 positional",
            "conclusion":"Only F2D changes through F4; degree6",
            "boundary":"No all-orders cap-preservation or R3=0 assertion"},
        "order_bookkeeping":{"h4_rhs":[[a,b,k,str(q)] for a,b,k,q in rhs_indices(4)],
                             "h5_rhs":[[a,b,k,str(q)] for a,b,k,q in rhs_indices(5)]},
        "parameter_consequence":"Every compatible normalized full tuple is determined by (a,b,c,u,v,w); if the separate h4 family imposes a=544c/27, by (b,c,u,v,w). No existence or higher compatibility is inferred.",
        "hashes":{"GRADED_LONG_TAIL_RIGIDITY.md":sha(HERE/"GRADED_LONG_TAIL_RIGIDITY.md"),
                  "check_graded_long_tail_rigidity.py":sha(Path(__file__)),
                  "graded_long_tail_rigidity_witness.json":sha(witnesspath),
                  "check_graded_long_cap23.py":sha(HERE/"check_graded_long_cap23.py"),
                  "check_long_graph_curvature.py":sha(HERE/"check_long_graph_curvature.py"),
                  "../dixmier-audit/verify_long_pbw.py":sha(source)},
        "coverage_boundary":"Complete homogeneous coefficient maps F3..F9 and finite termination. No inhomogeneous h4/h5 solution, compatibility sufficiency, bounded normalization, full lift or DC1/DC2 claim.",
        "independent_replay":False,"formal_verification":False,"external_review":False}
    output=args.output_dir or Path(tempfile.mkdtemp(prefix="graded-long-tail-rigidity-"))
    output.mkdir(parents=True,exist_ok=True)
    path=output/"graded_long_tail_rigidity_controls.json"
    path.write_text(json.dumps(receipt,indent=2)+"\n")
    print("PASS: F3 rank205/206 with kernel R^6; F4..F9 injective; m>=10 empty")
    print("Minor determinants:",[x["minor_det_modp"] for x in summary],"mod",p)
    print("Receipt:",path)
if __name__=="__main__":
    main()
