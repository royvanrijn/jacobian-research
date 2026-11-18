#!/usr/bin/env python3
"""Complete balanced degree26 Q37 obstruction, including nonsquare eta.

Fresh exact msolve checks use both sparse and sparse/dense linear algebra.
The square-eta case replays the existing complete labelled F37 census.
A separately retained Sage/Singular replay binds all 18 inputs directly
to the Euler/Chebyshev construction and verifies their unit ideals.
"""
from __future__ import annotations
import argparse,hashlib,importlib.util,json,re,shutil,subprocess,tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import sympy as sp

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
SYSTEM=HERE/"degree26_balanced_rational_quotient_system.json"
INDEPENDENT=HERE/"degree_twenty_six_nonsquare_reciprocal_independent_verification.json"
VARIABLE_NAMES=["x0","x1","x2","z0","z1","z2","eta"]
VARIABLES=sp.symbols(" ".join(VARIABLE_NAMES))
PRIME=37

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def polynomial_text(terms):
    out=[]
    for exponent,coefficient in sorted(terms.items(),reverse=True):
        coefficient%=37
        if not coefficient:continue
        factors=[str(coefficient)]
        for name,power in zip(VARIABLE_NAMES[:-1],exponent):
            if power:factors.append(name+("^"+str(power)if power>1 else""))
        out.append("*".join(factors))
    return"+".join(out)or"0"

def parse_safe_polynomial(text):
    """Independent roundtrip for our positive, explicit-multiplication format."""
    if text=="0":return{}
    index={name:i for i,name in enumerate(VARIABLE_NAMES[:-1])}
    result={}
    for monomial in text.split("+"):
        factors=monomial.split("*");coefficient=int(factors[0])
        assert 0<coefficient<37
        exponent=[0]*6
        for factor in factors[1:]:
            pieces=factor.split("^")
            assert len(pieces)<=2 and pieces[0]in index
            power=int(pieces[1])if len(pieces)==2 else 1
            assert 0<power<37
            exponent[index[pieces[0]]]+=power
        key=tuple(exponent)
        assert key not in result
        result[key]=coefficient
    return result

def load_inputs():
    data=json.loads(SYSTEM.read_text())
    assert data["k"]==3 and data["h"]==6 and data["degree"]==26
    assert data["variables"]==VARIABLE_NAMES
    assert data["strong_equation_count"]==13
    polys=[sp.Poly(sp.sympify(text),*VARIABLES,domain=sp.ZZ)
           for text in data["strong_primitive_integer_equations"]]
    assert len(polys)==13
    eta_zero=sp.expand(polys[-1].as_expr().subs(VARIABLES[-1],0))
    assert eta_zero==6692786100 and int(eta_zero)%37==30
    nonres=[ev for ev in range(1,37)if pow(ev,18,37)==36]
    assert nonres==[2,5,6,8,13,14,15,17,18,19,20,22,23,24,29,31,32,35]
    cases=[]
    for ev in nonres:
        rows=[];term_counts=[]
        for poly in polys:
            specialized={}
            for exponent,coefficient in poly.terms():
                key=exponent[:-1]
                specialized[key]=(specialized.get(key,0)+
                     int(coefficient)*pow(ev,exponent[-1],37))%37
            specialized={key:co for key,co in specialized.items()if co}
            text=polynomial_text(specialized)
            assert parse_safe_polynomial(text)==specialized
            rows.append(text);term_counts.append(len(specialized))
        source=",".join(VARIABLE_NAMES[:-1])+"\n37\n"+",\n".join(rows)+"\n"
        assert max(map(int,re.findall(r"\d+",source)))<=37
        cases.append(dict(eta=ev,source=source,term_counts=term_counts,
             input_sha256=hashlib.sha256(source.encode()).hexdigest()))
    return nonres,cases

def verify_independent_record(nonres):
    report=json.loads(INDEPENDENT.read_text())
    assert report["status"]=="PASS"and report["prime"]==37
    assert report["fresh_independent_singular_replay"]is True
    assert report["nonresidues"]==nonres
    assert [item["eta"]for item in report["records"]]==nonres
    for item in report["records"]:
        assert item["independently_reconstructed_all13_rows"]
        assert item["canonical_native_input_exact_match"]
        assert (ROOT/item["output"]).read_text().strip()==f'ETA {item["eta"]} UNIT'
    for path,expected in report["sha256"].items():
        source=(ROOT/path).resolve()
        assert source.is_relative_to(ROOT.resolve())
        assert digest(source)==expected,path
    return dict(status="PASS",all18_direct_ODE_bindings=True,
         all18_independent_Singular_unit_ideals=True,
         retained_record_sha256=digest(INDEPENDENT))

def square_case():
    path=HERE/"verify_degree_twenty_six_balanced_local_square.py"
    spec=importlib.util.spec_from_file_location("retained_balanced_square",path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    report=module.run()
    assert report["status"]=="PASS"and report["prime"]==37
    assert report["whole_labelled_strong_closure_has_Q37_point"]is False
    assert all(case["points"]==[]for case in report["censuses"])
    return dict(status="PASS",fresh_complete_labelled_census_replay=True,
         censuses=report["censuses"],
         exact_coefficient_binding_controls=len(report["polynomial_control_bindings"]))

def run():
    nonres,cases=load_inputs()
    independent=verify_independent_record(nonres)
    compiler=shutil.which("msolve")
    if compiler is None:raise RuntimeError("msolve is required")
    with tempfile.TemporaryDirectory(prefix="reciprocal26-nonsquare-")as folder:
        paths={}
        for case in cases:
            path=Path(folder)/f'eta_{case["eta"]}.ms'
            path.write_text(case["source"]);paths[case["eta"]]=path
        def solve(job):
            case,variant=job;ev=case["eta"]
            output=Path(folder)/f"eta_{ev}_l{variant}.gb"
            process=subprocess.run([compiler,"-f",str(paths[ev]),"-o",str(output),
                  "-t","2","-v","0","-l",str(variant),"-g","2"],
                   check=True,capture_output=True,text=True,timeout=120)
            text=output.read_text()
            body="".join(line.strip()for line in text.splitlines()
                          if line.strip()and not line.lstrip().startswith("#"))
            assert body=="[1]:",(ev,variant,body[:200])
            assert "error"not in process.stderr.lower()
            return dict(eta=ev,exact_linear_algebra_variant=variant,
                  unit_basis=[1],output_sha256=digest(output))
        jobs=[(case,variant)for case in cases for variant in [2,1]]
        with ThreadPoolExecutor(max_workers=2)as executor:
            solved=list(executor.map(solve,jobs))
    square=square_case()
    paths=[
       SYSTEM.relative_to(ROOT).as_posix(),
       INDEPENDENT.relative_to(ROOT).as_posix(),
       "research/arithmetic/RECIPROCAL_PARAMETER_FINITENESS_AND_INTEGRALITY.md",
       "research/arithmetic/BALANCED_RECIPROCAL_RATIONAL_QUOTIENT.md",
       "research/arithmetic/DEGREE_TWENTY_SIX_BALANCED_LOCAL_SQUARE_OBSTRUCTION.md",
       "research/arithmetic/verify_degree_twenty_six_balanced_local_square.py",
       "research/arithmetic/verify_degree_twenty_six_balanced_mod37.cpp",
       "research/arithmetic/DEGREE_TWENTY_SIX_RECIPROCAL_GLOBAL_OBSTRUCTION.md",
       Path(__file__).relative_to(ROOT).as_posix()]
    return dict(schema_version=1,status="PASS",prime=37,
       eta_zero_exclusion=dict(integer_constant=6692786100,residue=30),
       nonsquare_residues=nonres,
       nonsquare_cases=[dict(eta=case["eta"],input_sha256=case["input_sha256"],
                            equation_term_counts=case["term_counts"])
                         for case in cases],
       fresh_exact_msolve_eliminations=solved,
       independent_Singular_record=independent,
       square_case=square,
       msolve_binary_sha256=digest(Path(compiler)),
       all_emitted_integer_literals_at_most_37=True,
       discarded_unreduced_large_integer_experiment_used=False,
       whole_strong_quotient_has_Q37_point=False,
       whole_strong_quotient_has_point_over_residue_degree_one_Q37_extension=False,
       normalized_good_reciprocal_degree26_rational_model_exists=False,
       full_reciprocal_conclusion_uses_prior_uniform_and_unbalanced_obstructions=True,
       general_nonreciprocal_degree26_existence_claim=False,
       formal_verification=False,external_human_review=False,
       sha256={path:digest(ROOT/path)for path in paths})

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    result=run();payload=json.dumps(result,indent=2,sort_keys=True)+"\n"
    if args.output:args.output.write_text(payload)
    print(payload,end="")
if __name__=="__main__":main()
