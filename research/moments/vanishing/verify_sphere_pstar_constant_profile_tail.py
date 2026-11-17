"""Exact global index-seven and finite-tail certificate on a cubic sphere slice.

No Singular installation, Groebner computation, numerical moment inference,
or parameter localization is used. Literal polynomial witnesses are replayed.
The companion proof supplies the unrestricted mixed-cutoff tail theorem.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from math import factorial
from pathlib import Path
import platform
import sympy as sp
from sympy.polys.rings import ring
from sympy.polys.domains import QQ

if not __debug__:
    raise RuntimeError("Run without -O: exact assertions must remain enabled")

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
u,v,e,Z,Y,W=sp.symbols("u v e Z Y W")


def polynomial(text):
    value=sp.sympify(text,locals={"u":u,"v":v,"e":e})
    return sp.Poly(value,u,v,e,domain=QQ).as_expr()


def zero(value):
    assert sp.Poly(sp.expand(value),u,v,e,domain=QQ).is_zero


def phase_and_perturbation():
    r0=-u/3
    r1=-(2*u*u+3*v)/15
    r2=-4*u*(31*u*u+54*v)/945
    r3=-sp.Rational(1,2)-sp.Rational(2312,14175)*u**4-sp.Rational(188,525)*u*u*v-sp.Rational(8,75)*v*v
    phase=W+(u+v*Z)*Y**2+(Z+e*Z**2)*Y+r0+r1*Z+r2*Z**2+r3*Z**3
    return phase,phase-(W+Z*Y-Z**3/2)


def original_moments_through_eight():
    """Regenerate the actual sphere moments by sparse multinomial contraction."""
    R,U,V,E=ring("u,v,e",QQ)
    coefficients=[U,V,R.one,E,-U/3,-(2*U**2+3*V)/15,
                  -4*U*(31*U**2+54*V)/945,
                  -QQ(1,2)-QQ(2312,14175)*U**4-QQ(188,525)*U**2*V-QQ(8,75)*V**2]
    weights=[1,2,2,3,1,2,3,4]
    z_degree=[0,1,1,2,0,1,2,3]
    y_degree=[2,2,1,1,0,0,0,0]
    powers=[[coefficient**k for k in range(9)] for coefficient in coefficients]
    def vectors(n,index=0,prefix=()):
        if index==8:
            if n==0:
                yield prefix
            return
        for count in range(n//weights[index]+1):
            yield from vectors(n-count*weights[index],index+1,prefix+(count,))
    answer=[]
    for n in range(1,9):
        value=R.zero
        for counts in vectors(n):
            z=sum(a*b for a,b in zip(counts,z_degree))
            y=sum(a*b for a,b in zip(counts,y_degree))
            if y%2:
                continue
            scalar=(sp.factorial(n)*sp.factorial2(y-1)
                    /sp.factorial2(2*z+y+1))
            term=R.one
            for index,count in enumerate(counts):
                scalar/=sp.factorial(count)
                term*=powers[index][count]
            value+=QQ.convert(scalar)*term
        answer.append(value.as_expr())
    return answer


def check_lifts(primitive,small,lift):
    for j in range(7):
        zero(sum(primitive[i]*lift[i,j] for i in range(4))-small[j])


def checks(certificate):
    assert certificate["format"]=="sphere-pstar-global-nu7-literal-witnesses-v1"
    assert certificate["parameter_ring"]=="Q[u,v,e]"
    assert certificate["N"]==8 and certificate["nu"]==7
    for source in certificate["sources"]:
        path=ROOT/source["path"]
        assert hashlib.sha256(path.read_bytes()).hexdigest()==source["sha256"]
    moments=original_moments_through_eight()
    assert all(value==0 for value in moments[:4])
    for actual,stored in zip(moments[4:],certificate["raw_moments_5_8"]):
        zero(actual-polynomial(stored))
    primitive=list(map(polynomial,certificate["primitive_moments_5_8"]))
    for raw,scaled,scale in zip(moments[4:],primitive,certificate["primitive_over_raw_scalars"]):
        zero(scaled-sp.Rational(scale)*raw)
    small=list(map(polynomial,certificate["small_intermediate_generators"]))
    lift=sp.Matrix([[polynomial(value) for value in values]
                    for values in certificate["small_generators_from_primitive_moments"]])
    assert lift.shape==(4,7)
    check_lifts(primitive,small,lift)
    corrupted=lift.copy()
    corrupted[0,0]+=1
    try:
        check_lifts(primitive,small,corrupted)
    except AssertionError:
        pass
    else:
        raise AssertionError("The deliberately corrupted polynomial lift was accepted")

    uf=-e+sp.Rational(1184,155925)*e**5
    vf=-e**2/2+sp.Rational(6022904,638512875)*e**6
    J=[u-uf,v-vf,e**7]
    assert list(map(polynomial,certificate["curvilinear_generators"]))==J
    second=sp.Matrix([[polynomial(value) for value in values]
                      for values in certificate["curvilinear_from_small_generators"]])
    assert second.shape==(7,3)
    for j in range(3):
        zero(sum(small[i]*second[i,j] for i in range(7))-J[j])
    reverse=sp.Matrix([[polynomial(value) for value in values]
                       for values in certificate["primitive_moments_from_curvilinear_generators"]])
    assert reverse.shape==(4,3)
    for i in range(4):
        zero(primitive[i]-sum(reverse[i,j]*J[j] for j in range(3)))

    phase,H=phase_and_perturbation()
    assert sp.expand(phase-sp.sympify(certificate["phase"],locals={"u":u,"v":v,"e":e,"Z":Z,"Y":Y,"W":W}))==0
    # Independent quotient moments: inverse shear, then finite nilpotent binomial.
    G=sp.expand(H.subs(Y,Y+Z**2).subs({u:uf,v:vf}))
    R,ee,zz,yy=ring("e,Z,Y",QQ)
    def trunc(value):
        return R.from_dict({powers:coefficient for powers,coefficient in value.items() if powers[0]<7})
    g=trunc(R.from_expr(G))
    assert all(powers[0]>=1 and powers[1]<=5 for powers in g)
    power=R.one
    remainders=[R.zero for _ in range(61)]
    sizes=[]
    max_possible_moment=0
    for j in range(7):
        if j:
            power=trunc(power*g)
        sizes.append(len(power))
        for (ep,a,b),coefficient in power.items():
            n=a+j
            max_possible_moment=max(max_possible_moment,n)
            if n==0 or n>60 or b%2:
                continue
            scalar=(sp.binomial(n,j)*sp.factorial(a)*sp.factorial2(b-1)
                    /sp.factorial2(2*a+b+1))
            remainders[n]+=coefficient*QQ.convert(scalar)*ee**ep
    assert trunc(power*g)==0
    assert all(value==0 for value in remainders[1:])
    assert sizes==certificate["independent_finite_memberships"]["nilpotent_power_term_counts"]
    assert max_possible_moment<=36
    return {
        "status":"passed","functional":"normalized S2 sphere",
        "source":"literal original sparse multinomial moments1..8",
        "polynomial_lift_identities":14,"corrupted_lift_rejected":True,
        "global_ideal":"I8=(u+e-1184e5/155925,v+e2/2-6022904e6/638512875,e7)",
        "quotient":"Q[e]/(e7)","exact_nilpotence_index":7,
        "intervening_moment_memberships":[9,60],"all_remainders_zero":True,
        "nilpotent_power_term_counts":sizes,"largest_possible_quotient_moment":max_possible_moment,
        "scope":"Exact polynomial witnesses prove the global ideal equality and length/index seven on this specified graph slice. Original finite moments and nilpotent-binomial controls are checked; unrestricted tail conclusions use the written proof.",
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    parser.add_argument("--checks-only",action="store_true")
    args=parser.parse_args()
    witness=HERE/"sphere_pstar_constant_profile_tail_certificate.json"
    result=checks(json.loads(witness.read_text(encoding="utf-8")))
    if args.checks_only:
        print(json.dumps(result,indent=2))
        return
    proof=HERE/"MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md"
    receipt={"status":"passed","checked_on":"2026-11-17",
             "environment":{"python":platform.python_version(),"sympy":sp.__version__},
             "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             "certificate_sha256":hashlib.sha256(witness.read_bytes()).hexdigest(),
             "proof_sha256":hashlib.sha256(proof.read_bytes()).hexdigest(),
             "checks":result,"independent_replay":False,"formal_verification":False,"external_review":False}
    if args.output:
        args.output.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(receipt,indent=2))


if __name__=="__main__":
    main()
