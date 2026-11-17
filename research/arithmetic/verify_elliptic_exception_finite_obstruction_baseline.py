#!/usr/bin/env python3
"""Exact baseline certificate; optional arithmetic on a supplied weight table."""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import gcd, isqrt, prod
from pathlib import Path
import argparse
import json
import platform
import runpy

HERE=Path(__file__).resolve().parent
DEPENDENCY="verify_elliptic_exception_one_prime_positivity.py"


def classify_weight_table(table):
    """EFO2 arithmetic only: FQS image/weight certification is an external input."""
    primes=table["primes"]
    cells=table["cells"]
    a=table["class"]
    order=table["quotient_order"]
    if len(set(primes))!=len(primes) or len(cells)!=order or order<=0:
        raise ValueError("inconsistent finite weight table")
    if len({c["id"] for c in cells})!=order:
        raise ValueError("cell identifiers must be distinct")
    counts=Counter(c["class"] for c in cells)
    classes={1,2,13545865,27091730}
    if set(counts)!=classes or any(4*counts[x]!=order for x in classes) or a not in classes:
        raise ValueError("the full quotient must retain the uniform four-class map")
    rows=[]
    for cell in cells:
        weights={p:F(cell["weights"][str(p)]) for p in primes}
        if any(w<0 or w>1 for w in weights.values()):
            raise ValueError("local weights must be in [0,1]")
        if cell["class"]==a:
            rows.append((cell["id"],weights,{p for p,w in weights.items() if not w}))
    coefficient=sum((prod(row[1].values()) for row in rows),F(0))/order
    cover=all(row[2] for row in rows)
    critical={p:[name for name,weights,zeros in rows if zeros=={p}] for p in primes}
    deleted={p:sum((prod(w for q,w in weights.items() if q!=p)
                    for name,weights,zeros in rows),F(0))/order for p in primes}
    return {"coefficient":str(coefficient),
            "zero_weight_cover":{name:sorted(zeros) for name,weights,zeros in rows},
            "obstruction_if_supplied_input_is_certified":bool(cover),
            "minimal_if_supplied_input_is_certified":bool(cover and all(critical.values())),
            "unique_zero_cells":{str(p):ids for p,ids in critical.items()},
            "one_prime_deletion_coefficients":{str(p):str(w) for p,w in deleted.items()},
            "full_rational_image_and_local_weights_certified_by_this_script":False,
            "required_input_scope":"Exact full quotient and local weights certified as in FQS2-FQS3"}


def replay():
    # Import exact arithmetic helpers without invoking the dependency's replay.
    dep=runpy.run_path(str(HERE/DEPENDENCY),run_name="baseline_arithmetic_dependency")
    cert=json.loads((HERE/"elliptic_exception_finite_obstruction_baseline.json").read_text())
    checks=[]
    def check(label,condition):
        if not condition:
            raise AssertionError(label)
        checks.append(label)
    a,b=cert["curve"]["a2"],cert["curve"]["a4"]
    points={name:tuple(F(q) for q in point) for name,point in cert["points"].items()}
    for name,point in points.items():
        check("exact curve point "+name,dep["on_curve"](point,a,b))
    constructed=dep["add"](dep["add"](dep["mul"](points["P"],2,a),
                                     dep["mul"](points["Q"],3,a),a),points["T"],a)
    check("literal U group-law certificate",constructed==points["U"])
    data={}
    for name in ["Q","U"]:
        point=points[name]
        c=isqrt(point[0].denominator)
        aa=point[0].numerator
        bb=point[1]*c**3
        check("square denominator "+name,c*c==point[0].denominator)
        check("cube denominator "+name,bb.denominator==1)
        bb=int(bb)
        check("primitive coordinates "+name,gcd(aa,c)==gcd(bb,c)==1)
        alpha=cert["classes"][name]
        check("exact global squareclass "+name,dep["square"](point[0]/alpha))
        q,h,D,Delta=dep["factors_at"](point)
        check("irreducible quadratic "+name,D<0)
        check("separable cubic "+name,Delta!=0)
        check("complete fibre factors coprime "+name,dep["gcd_degree"](q,h)==0)
        check("matching resolvents "+name,dep["square"](Delta/D))
        p=cert["cubic_irreducible_primes"][name]
        check("irreducibility witness prime "+name,dep["prime"](p))
        check("integral monic irreducibility reduction "+name,all(z.denominator%p for z in h))
        roots=[t for t in range(p) if sum(dep["modq"](coef,p)*pow(t,i,p)
                                         for i,coef in enumerate(h))%p==0]
        check("complete cubic residue-root exclusion "+name,not roots)
        val=dep["valuation"]
        if c%2:
            check("two-adic translated component "+name,val(point[0],2)>=9)
        else:
            check("two-adic permitted formal depth "+name,val(c,2)%3==0)
        check("five-adic integral coordinate "+name,val(c,5)==0)
        check("three-adic integral coordinate "+name,val(c,3)==0)
        check("three-adic soluble residue and depth "+name,
              aa%3==2 and val(bb,3)>=2 or aa%3==1 and val(bb,3)>=3)
        data[name]={"class":alpha,"irreducible_cubic_prime":p,
                    "c_v2":val(c,2),"c_v3":val(c,3),"c_v5":val(c,5),
                    "a_mod3":aa%3,"b_v3":val(bb,3)}
    expected=cert["U_reduced_coordinate_claims"]
    check("literal U local coordinate claims",all(data["U"][k]==v for k,v in expected.items()))
    wrong=(points["U"][0]+1,points["U"][1])
    check("corrupted U coordinate rejected",not dep["on_curve"](wrong,a,b))
    return {"passed":True,"exact_controls":len(checks),"checks":checks,
            "baseline_witnesses":data,
            "actual_finite_obstruction_found":False,
            "actual_obstruction_weight_table_supplied":False,
            "weight_table_routine_scope":"Exact arithmetic assuming independently certified FQS inputs",
            "finite_checks_establish_unrestricted_obstruction_criterion":False,
            "all_finite_sets_in_other_two_classes_claimed":False,
            "formal_verification":False,"external_review":False,"independent_replay":False,
            "python":platform.python_version(),"arithmetic":"standard-library exact integers and fractions"}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weights",type=Path,help="classify an externally certified full finite weight table")
    args=parser.parse_args()
    if args.weights:
        print(json.dumps(classify_weight_table(json.loads(args.weights.read_text())),indent=2))
    else:
        report=replay()
        names=["ELLIPTIC_EXCEPTION_FINITE_OBSTRUCTION_CERTIFICATES.md",Path(__file__).name,
               "elliptic_exception_finite_obstruction_baseline.json",DEPENDENCY]
        report["input_sha256"]={name:sha256((HERE/name).read_bytes()).hexdigest() for name in names}
        (HERE/"elliptic_exception_finite_obstruction_baseline_checks.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
        print(json.dumps({"passed":True,"exact_controls":report["exact_controls"],
                          "baseline_witnesses":report["baseline_witnesses"],
                          "input_sha256":report["input_sha256"]},indent=2))
