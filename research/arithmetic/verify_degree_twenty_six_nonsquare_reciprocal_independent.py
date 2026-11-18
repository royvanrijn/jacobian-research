#!/usr/bin/env python3
"""Independent direct-ODE binding and Singular verification of all 18 cases.

Run with Sage Python. No retained coefficient generator is imported.
Default replay rebuilds every polynomial over F37^2, verifies every native
input and runs ordinary exact Singular slimgb on all 18 specializations.
"""
from __future__ import annotations
import argparse,hashlib,json,shutil,subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from sage.all import GF,PolynomialRing,binomial

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
SYSTEM=HERE/"degree26_balanced_rational_quotient_system.json"
INPUTS=HERE/"degree_twenty_six_reciprocal_global_inputs"

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def run():
    data=json.loads(SYSTEM.read_text())
    assert data["variables"]==["x0","x1","x2","z0","z1","z2","eta"]
    names=data["variables"][:-1]
    F=GF(37)
    R7=PolynomialRing(F,names=data["variables"],order="degrevlex")
    canonical=[R7(text)for text in data["strong_primitive_integer_equations"]]
    assert len(canonical)==13
    R6=PolynomialRing(F,names=names,order="degrevlex")
    xs=R6.gens()
    nonres=[ev for ev in range(1,37)if pow(ev,18,37)==36]
    assert len(nonres)==18
    records=[];pins={SYSTEM.relative_to(ROOT).as_posix():digest(SYSTEM),
                     Path(__file__).relative_to(ROOT).as_posix():digest(Path(__file__))}
    for ev in nonres:
        specialization=R7.hom(list(xs)+[R6(ev)],R6)
        expected=[specialization(poly)for poly in canonical]
        field_polynomials=PolynomialRing(F,"tt")
        tt=field_polynomials.gen()
        K=GF(37**2,name="rr",modulus=tt*tt-ev)
        rr=K.gen()
        RK=PolynomialRing(K,names=names,order="degrevlex")
        xx=RK.gens()
        W=PolynomialRing(RK,"ww");ww=W.gen()
        ps=[(xx[i]+rr*xx[i+3])/K(2)for i in range(3)]
        qs=[(-1)**(3+i)*(xx[i]-rr*xx[i+3])/K(2)for i in range(3)]
        P=ww**3+sum(ps[i]*ww**i for i in range(3))
        Q=ww**3+sum(qs[i]*ww**i for i in range(3))
        s=P*Q
        sigma=[RK(0)for j in range(13)]
        for k in range(7):
            for ell in range(k+1):
                sigma[6+k-2*ell]+=s[k]*binomial(k,ell)
        coefficients=[]
        for j in range(14):
            numerator=(sigma[j-1]if j else RK(0))-rr*(
                       sigma[j]if j<13 else RK(0))
            coefficients.append(K(25)/K(2*j-1)*numerator)
        assert coefficients[-1]==1
        C=[W(1),ww-1]
        for j in range(1,12):C.append(ww*C[-1]-C[-2])
        X=W(coefficients[0])+sum(W(coefficients[j])*C[j-1]
                                      for j in range(1,14))
        remainder=X%(P*P)
        def split(poly):
            even={};odd={}
            for exponent,coefficient in RK(poly).dict().items():
                entries=coefficient.polynomial()
                if entries[0]:even[exponent]=F(entries[0])
                if entries[1]:odd[exponent]=F(entries[1])
            return R6(even),R6(odd)
        rows=[]
        for j in range(6):rows.extend(split(remainder[j]))
        G=sum(coefficients[j]*rr**j for j in range(14))/(K(2)*rr)
        ge,go=split(G);assert go==0;rows.append(ge)
        scales=[]
        for actual,target in zip(rows,expected):
            assert actual and target
            exponent=next(iter(actual.dict()))
            scale=target.dict()[exponent]/actual.dict()[exponent]
            assert scale and scale*actual==target
            scales.append(int(scale))
        path=INPUTS/f"eta_{ev}.sing"
        text=path.read_text()
        ideal_text=text.split("ideal I=",1)[1].split(";",1)[0]
        native_polynomials=[R6(piece.strip())for piece in ideal_text.split(",")]
        assert native_polynomials==expected
        output=INPUTS/f"eta_{ev}.stdout"
        assert output.read_text().strip()==f"ETA {ev} UNIT"
        for source in [path,output]:
            pins[source.relative_to(ROOT).as_posix()]=digest(source)
        records.append(dict(eta=ev,
             independently_reconstructed_all13_rows=True,
             primitive_row_scalars=scales,
             canonical_native_input_exact_match=True,
             input=path.relative_to(ROOT).as_posix(),
             output=output.relative_to(ROOT).as_posix(),
             canonical_specialization_term_counts=[len(poly.dict())for poly in expected]))
    singular=shutil.which("Singular")
    if singular is None:raise RuntimeError("Singular is required")
    def verify(record):
        process=subprocess.run([singular,"-q",str(ROOT/record["input"])],
                         check=True,capture_output=True,text=True,timeout=120)
        assert not process.stderr.strip()
        assert process.stdout.strip()==f'ETA {record["eta"]} UNIT'
        assert process.stdout.encode()==(ROOT/record["output"]).read_bytes()
        return record["eta"]
    with ThreadPoolExecutor(max_workers=2)as executor:
        verified=list(executor.map(verify,records))
    assert verified==nonres
    return dict(schema_version=1,status="PASS",prime=37,
        fresh_independent_singular_replay=True,
        method="Independent direct Euler/Chebyshev reconstruction over F37[r]/(r^2-eta), followed by ordinary exact Singular slimgb; no retained generator imported.",
        nonresidues=nonres,records=records,sha256=pins,
        formal_verification=False,independent_external_review=False)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    result=run();payload=json.dumps(result,indent=2,sort_keys=True)+"\n"
    if args.output:args.output.write_text(payload)
    print(payload,end="")
if __name__=="__main__":main()
