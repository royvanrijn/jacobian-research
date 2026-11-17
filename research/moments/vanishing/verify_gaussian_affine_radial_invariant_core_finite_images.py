#!/usr/bin/env python3
"""Bounded original-source controls for affine radial invariant cores.

Recorded finite identities do not prove the all-model volume bounds or
complete analytic comparison in the accompanying written proof.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sympy as sp

if not __debug__:
    raise RuntimeError("Run without -O: exact assertions must remain enabled")

def affine_radial_controls():
    Z,W,Y,S,t,eps,alpha=sp.symbols('Z W Y S t eps alpha')
    order=4
    def tr(value,n=order):
        return sp.expand(sum(co*t**int(e[0]) for e,co in sp.Poly(sp.expand(value),t).terms() if int(e[0])<=n))
    def ey(value):
        return sp.expand(sum(co*sp.factorial2(k-1) for (k,),co in sp.Poly(sp.expand(value),Y).terms() if k%2==0))
    def wick(value):
        return sp.expand(sum(co*sp.factorial(i)*sp.factorial2(k-1)
                             for (i,j,k),co in sp.Poly(sp.expand(value),Z,W,Y).terms()
                             if i==j and k%2==0))
    def q_inverse(H,k,n=order):
        if k==0:
            return sp.Integer(1)
        return sp.expand(sum(sp.binomial(k+j-1,j)*t**j*H**j for j in range(n+1)))
    def exp_G(G,n=order):
        return sp.expand(sum(t**j*G**j/sp.factorial(j) for j in range(n+1)))
    fixtures=[]
    for H,G in [(Y+1,Y**3+2*Y+3),((Y-1)**3,Y**4+Y),
                (Y**2+Y+1,Y**2-Y+3),(sp.Integer(2),2*Y+3),
                (sp.Integer(0),2*Y+3)]:
        phase=Z*W*H+G
        powers=[sp.Integer(1)]
        for _ in range(order):
            powers.append(sp.expand(powers[-1]*phase))
        E=exp_G(G)
        def actual(marker):
            return sp.expand(sum(wick(marker*powers[n])*t**n/sp.factorial(n) for n in range(order+1)))
        for m in range(3):
            for j in range(3):
                original=actual(Z**m*W**m*Y**j)
                predicted=sp.factorial(m)*ey(tr(Y**j*E*q_inverse(H,m+1)))
                assert sp.expand(original-predicted)==0
        assert actual(Z)==0 and actual(W)==0
        # q times the actual unit recovers the endpoint Gaussian source.
        assert tr(actual(1)-t*actual(H)-ey(E))==0
        for k in range(3):
            F=Y*q_inverse(H,k)
            spatial=tr(sp.diff(F,Y)-Y*F+t*sp.diff(G,Y)*F)
            assert ey(tr(spatial*E))==0
            value=ey(tr(F*E))
            derivative=tr(sp.diff(value,t),order-1)
            inside=tr(k*Y*H*q_inverse(H,k+1)+G*F,order-1)
            expected=ey(tr(inside*E,order-1))
            assert sp.expand(derivative-expected)==0
        fixtures.append({'H':str(H),'G':str(G),'actual_unit_jet':str(actual(1)),
                         'balanced_markers_m_j':[0,2],'unbalanced_markers':['Z','W'],
                         'tracked_endpoint_q_identity':True,
                         'localized_Y_Stein_orders':[0,order],
                         'localized_parameter_orders':[0,order-1]})
    residues=[]
    for m in range(1,7):
        H=(Y-alpha)**m
        q=1-eps**(-m)*H
        residue=sp.cancel(1/sp.diff(q,Y).subs(Y,alpha+eps))
        assert sp.cancel(residue+eps/m)==0
        assert sp.gcd(q,sp.diff(q,Y),Y)==1
        residues.append({'multiplicity':m,'root':'alpha+eps','r':'eps^m',
                         'exact_log_residue':str(residue),'positive_eps_order':1})
    return {'actual_source_fixtures':fixtures,'original_moment_orders':[0,order],
            'moving_log_residue_controls':residues,
            'scope':'Bounded exact source and residue identities; all-model volume orders and comparison are written arguments.'}

def finite_controls():
    Z,W,Y,t=sp.symbols("Z W Y t")
    order=6
    def truncate(value):
        return sp.expand(sum(co*t**int(power[0]) for power,co in sp.Poly(sp.expand(value),t).terms()
                             if int(power[0])<=order))
    def wick(value):
        return sp.expand(sum(co*sp.factorial(i)*sp.factorial2(k-1)
            for (i,j,k),co in sp.Poly(sp.expand(value),Z,W,Y).terms() if i==j and k%2==0))
    def exp_series(value):
        answer=1
        power=1
        for n in range(1,order+1):
            power=truncate(power*value)
            answer+=power/sp.factorial(n)
        return truncate(answer)
    R=Z**5*W**3*Y**4+Z**2
    rows=[]
    for c,b in [(0,0),(3,0),(3,2)]:
        P=c+b*Y+R
        powers=[sp.Integer(1)]
        for _ in range(order):
            powers.append(sp.expand(powers[-1]*P))
        factor=exp_series(c*t+b*b*t*t/2)
        markers=[]
        for Q in [sp.Integer(1),W**2,Y*W**2,W**4]:
            D=sp.Poly(Q,Z,W,Y).total_degree()
            normalized=sp.expand(sum(t**n*wick(Q.subs(Y,Y+b*t)*R.subs(Y,Y+b*t)**n)/sp.factorial(n)
                                     for n in range(D//2+1)))
            original=sum(wick(Q*powers[n])*t**n/sp.factorial(n) for n in range(order+1))
            assert sp.expand(original-truncate(factor*normalized))==0
            if Q==1:
                assert normalized==1
            if Q==W**2:
                assert sp.expand(normalized-(362*t+720*b*b*t**3+120*b**4*t**5))==0
            bound=D//2 if b==0 else D+5*(D//2)
            assert normalized==0 or sp.degree(normalized,t)<=bound
            markers.append({"marker":str(Q),"normalized_finite_formula":str(normalized),
                            "sufficient_degree_bound":bound})
        rows.append({"H":"0","G":str(c+b*Y),"R":"Z^5 W^3 Y^4+Z^2",
                     "gamma":2,"degree_Y_R":4,"original_moment_orders":[0,order],"markers":markers})
    return rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,help="Write the exact JSON receipt to this path")
    args=parser.parse_args()
    checks={"actual_source":affine_radial_controls(),"finite_branches":finite_controls()}
    here=Path(__file__).resolve()
    repository=here.parents[3]
    canonical=[here]+[here.with_name(name) for name in [
        "GAUSSIAN_AFFINE_RADIAL_INVARIANT_CORE_FINITE_IMAGES.md",
        "GAUSSIAN_ANGULAR_CORE_FULL_IMAGE_INVARIANCE.md",
        "GAUSSIAN_AFFINE_ISOTROPIC_CHART.md",
        "GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md",
        "GAUSSIAN_SAFE_RECOMPOSITION_CLASSIFICATION.md",
        "MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md",
    ]]
    base=repository/"unpublished/gaussian-boundary"
    imports=[base/name for name in ["seed_period_comparison.tex","analytic_core.tex",
                                   "polynomial_stationary_balls.tex","sources/laurent_primary_manifest.json"]]
    manifest=json.loads((base/"sources/laurent_primary_manifest.json").read_text(encoding="utf-8"))
    primary=[]
    for entry in manifest["sources"]:
        path=base/entry["file"]
        assert hashlib.sha256(path.read_bytes()).hexdigest()==entry["sha256"],entry["file"]
        primary.append(path)
    checks["primary_manifest_consistency"]={"status":"passed","files":[e["file"] for e in manifest["sources"]]}
    paths=canonical+imports+primary
    receipt={
        "status":"passed","checked_on":"2026-11-17",
        "environment":{"python":platform.python_version(),"sympy":sp.__version__},
        "sha256":{p.relative_to(repository).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        "checks":checks,
        "independent_replay":False,"formal_verification":False,"external_review":False,
        "unrestricted_result_from_finite_prefix":False,
        "scope":"Bounded exact original Gamma elimination, localized Gaussian source and moving-log residue controls. All-model volume bounds and complete moving-puncture comparison are the pinned written proof with retained review limits.",
    }
    encoded=json.dumps(receipt,indent=2)+"\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(encoded,encoding="utf-8")
    print(encoded,end="")

if __name__=="__main__":
    main()
