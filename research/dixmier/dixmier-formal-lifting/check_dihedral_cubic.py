#!/usr/bin/env python3
"""QCS1 binary-dihedral quantum normal-closure controls; no Long identification."""
from __future__ import annotations
import argparse
from pathlib import Path
from math import comb
import hashlib
import json
import platform
import tempfile
import sympy as s

HERE=Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    h,x,t,z,omega=s.symbols("h x t z omega")
    def clean(op):
        return {i:s.factor(c) for i,c in op.items() if s.factor(c)!=0}
    def add(*ops):
        out={}
        for op in ops:
            for i,c in op.items():
                out[i]=out.get(i,0)+c
        return clean(out)
    def scale(op,c):
        return clean({i:c*a for i,a in op.items()})
    def mul(left,right,var):
        out={}
        for i,a in left.items():
            for j,b in right.items():
                for k in range(i+1):
                    n=i+j-k
                    out[n]=out.get(n,0)+a*comb(i,k)*h**k*s.diff(b,var,k)
        return clean(out)
    def comm(left,right,var):
        return add(mul(left,right,var),scale(mul(right,left,var),-1))
    def same(left,right):
        return all(s.factor(left.get(i,0)-right.get(i,0))==0
                   for i in set(left)|set(right))
    def omega_zero(expr):
        num,den=s.fraction(s.cancel(expr))
        n=s.rem(s.expand(num),omega**2+omega+1,omega)
        d=s.rem(s.expand(den),omega**2+omega+1,omega)
        assert s.factor(d)!=0
        return s.factor(n)==0

    # Exact original-field identities in finite x-normal order.
    P={1:s.Integer(1)}
    X={0:x}
    torig=mul(P,{0:1/x},x)
    Zorig={0:x*x}
    assert comm(torig,Zorig,x)=={0:2*h}
    assert mul(torig,X,x)==P
    assert same(add(mul(mul(torig,torig,x),Zorig,x),scale(torig,-h)),mul(P,P,x))

    # On the central fixed field, represent Z=-2A, [A,t]=h.
    Z={1:s.Integer(-2)}
    zt=-s.I*t
    b=add(scale(Z,-s.I/2),{0:h/(4*zt)})
    assert comm(b,{0:zt},t)=={0:h}
    jZ=add(mul({0:t*t},Z,t),{0:-h*t})
    jzt=1/zt
    jb=add(scale(jZ,-s.I/2),{0:h/(4*jzt)})
    assert same(jb,scale(b,-zt*zt))
    braw=scale(Z,-s.I/2)
    jbraw=scale(jZ,-s.I/2)
    raw_defect=add(jbraw,scale(braw,zt*zt))
    assert raw_defect=={0:s.I*h*t/2} and raw_defect!={}
    assert comm(jb,{0:jzt},t)=={0:h}

    # Finite dicyclic group, including its exact index-three subgroup.
    def times(a,b):
        return ((a[0]+(-1 if a[1] else 1)*b[0]+(3 if a[1] and b[1] else 0))%6,
                (a[1]+b[1])%2)
    identity=(0,0); g=(1,0); j=(0,1)
    elements={(a,b) for a in range(6) for b in range(2)}
    assert len(elements)==12 and all(times(a,b) in elements for a in elements for b in elements)
    assert all(times(times(a,b),c)==times(a,times(b,c))
               for a in elements for b in elements for c in elements)
    def power(a,n):
        v=identity
        for _ in range(n):v=times(v,a)
        return v
    H={power(j,n) for n in range(4)}
    center={identity,(3,0)}
    assert len(H)==4 and power(j,2)==power(g,3)
    cosets=set(frozenset(times(a,b) for b in H) for a in elements)
    assert len(cosets)==3
    cosets=sorted(cosets,key=lambda q:sorted(q))
    perms={}
    for a in elements:
        perms[a]=tuple(cosets.index(frozenset(times(a,b) for b in q)) for q in cosets)
    assert len(set(perms.values()))==6
    assert {a for a,perm in perms.items() if perm==(0,1,2)}==center

    # Coefficient-field quotient maps and exact ordered cotangent lifts.
    yy=-3*z/(z*z+z+1)
    uu=27*z**3/(z**3-1)**2
    fy=-s.Symbol("Y")**3/(s.Symbol("Y")+1)
    assert s.factor(uu-fy.subs(s.Symbol("Y"),yy))==0
    yp,up=s.diff(yy,z),s.diff(uu,z)
    eta={1:-1/yp}
    zeta={1:-1/up}
    assert comm({0:yy},eta,z)=={0:h}
    assert comm({0:uu},zeta,z)=={0:h}
    fprime=s.diff(fy,s.Symbol("Y")).subs(s.Symbol("Y"),yy)
    assert same(eta,mul({0:fprime},zeta,z))
    assert s.factor(yy.subs(z,1/z)-yy)==0
    assert s.factor(uu.subs(z,1/z)-uu)==0
    assert omega_zero(uu.subs(z,z/omega)-uu)
    assert not omega_zero(yy.subs(z,z/omega)-yy)
    def transform_first(op,phi):
        result={}
        for order,coefficient in op.items():
            assert order in (0,1)
            factor=s.diff(phi,z)**(-order)
            result[order]=s.factor(coefficient.subs(z,phi)*factor)
        return result
    assert same(transform_first(eta,1/z),eta)
    assert same(transform_first(zeta,1/z),zeta)
    assert all(omega_zero(transform_first(zeta,z/omega).get(i,0)-zeta.get(i,0))
               for i in set(zeta))
    deg_y=max(s.degree(q,z) for q in s.fraction(s.cancel(yy)))
    deg_u=max(s.degree(q,z) for q in s.fraction(s.cancel(uu)))
    assert deg_y==2 and deg_u==6
    Y,U=s.symbols("Y U")
    cubic_discriminant=s.factor(s.discriminant(Y**3+U*Y+U,Y))
    quadratic=Y*z*z+(Y+3)*z+Y
    assert s.factor(quadratic.subs(Y,yy))==0
    quadratic_discriminant=s.factor(s.discriminant(quadratic,z))
    assert s.factor(quadratic_discriminant+3*(Y-3)*(Y+1))==0

    output=args.output or Path(tempfile.mkdtemp(prefix="qcs1-dihedral-"))/"dihedral_cubic_controls.json"
    output.parent.mkdir(parents=True,exist_ok=True)
    receipt={
        "claim":"QCS1","scope":"Explicit binary-dihedral quantum S3 normal closure",
        "status":"PASS","environment":{"python":platform.python_version(),"sympy":s.__version__},
        "original_relations":{"t":"p*x^-1","Z":"x^2","commutator_t_Z":"2h","p_squared":"t^2*Z-h*t"},
        "shifted_generators":{"z":"-i*t","b":"-i*Z/2+h/(4*z)","commutator_b_z":"h"},
        "projective_action":{"g":["omega^-1*z","omega*b"],"j":["1/z","-z^2*b"]},
        "unshifted_generator_cotangent_defect":"-h*z/2",
        "binary_dihedral_order":12,"subgroup_order":4,"subgroup_index":3,
        "coset_action_order":6,"coset_action_kernel":"{1,-I}",
        "coefficient_invariants":{"y":str(yy),"u":str(uu),"u_equals_f_y":True,
                                  "degree_z_over_y":2,"degree_z_over_u":6,
                                  "quadratic_normal_closure_discriminant":str(quadratic_discriminant),
                                  "cubic_discriminant":str(cubic_discriminant)},
        "quantum_invariants":{"eta":"-y_prime(z)^-1*b","zeta":"-u_prime(z)^-1*b",
                              "commutator_y_eta":"h","commutator_u_zeta":"h",
                              "eta_equals_fprime_y_zeta":True,
                              "reflection_eta_invariance":True,
                              "full_S3_zeta_invariance":True},
        "fixed_field_and_index_scope":"Written finite differential-field basis and division spanning proof; averaging acts on coefficient basis only. D^H/D^G has left and right degree three; D^{+-I}/D^G is the explicit S3 normal closure.",
        "classical_bridge":"This is exactly the QCS cotangent left column. No Poisson equivalence to the Long right column or compatible quintic root is supplied.",
        "literature":{"primary_general_statement":"https://doi.org/10.1080/00927879708825943",
                      "general_original_proof_inspected":False,
                      "author_exposition_inspected":"https://lmbp.uca.fr/~fdumas/fichiers/B-Mat60.pdf",
                      "exact_exposition_location":"Section 5.3.2, theorem/proof and binary-dihedral Example 2, printed page 64 (PDF page 69); heading printed page 63",
                      "cached_file":"sources/dumas-2011-noncommutative-invariants.pdf",
                      "cached_sha256":sha(HERE/"sources/dumas-2011-noncommutative-invariants.pdf"),
                      "theorem_scope":"Finite subgroups of Aut A1(C) extending to D1(C), not arbitrary finite subgroups of Aut D1(C)",
                      "proof_premise":False},
        "hashes":{name:sha(HERE/name) for name in ["CLASSICAL_QUINTIC_COVER_SQUARE.md","check_dihedral_cubic.py"]},
        "coverage_boundary":"Exact rational Weyl identities, h-correction, finite group/coset enumeration, coefficient invariants and ordered mate identities. Fixed division field equalities, universal degree statements and the Poisson comparison are written arguments. No actual Long quantum identification, compatible fifth root, global polynomial admission, bounded lift or DC1/DC2 resolution.",
        "independent_replay":False,"formal_verification":False,"external_review":False
    }
    output.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    print("QCS1 binary-dihedral controls PASS")
    print("Verification output:",output)

if __name__=="__main__":
    main()
