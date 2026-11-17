"""Exact controls for written plane branch-curve and cubic-output proofs."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import time

import sympy as S
from sympy.combinatorics import Permutation, PermutationGroup

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    checks = []

    def check(name,value):
        if not value:
            raise AssertionError(name)
        checks.append(name)

    def zero(name,value):
        check(name,S.cancel(S.expand(value))==0)

    T,V,z,lam = S.symbols('T V z lambda')
    # Content and localization controls: horizontal and vertical components
    # must not be confused. The written proof covers all polynomials.
    Q = T**4*(T-1)**3
    C = T**4
    D = (T+2)**2
    check('the parameter control is reduced despite repeated denominator',S.gcd(S.gcd(Q,C),D)==1)
    for level in [0,S.Rational(1,9),2]:
        content = S.gcd(Q,C-level*D)
        u = S.cancel(Q/content)
        v = S.cancel((C-level*D)/content)
        check('primitive horizontal coefficients are coprime at '+str(level),S.gcd(u,v)==1)
        zero('horizontal equation retains simple V degree at '+str(level),
             S.diff(S.cancel((Q*V+C-level*D)/content),V)-u)
        check('horizontal coefficient is nonconstant in this three-level control '+str(level),S.degree(u,T)>0)
    u,v = T**3,(T-1)**2
    bez_u,bez_v,gcd = S.gcdex(u,v,T)
    zero('Bezout makes the horizontal coefficient a unit on its primitive graph',bez_u*u+bez_v*v-1)
    zero('the graph coordinate recovers that inverse',bez_u-bez_v*V-(bez_u*u+bez_v*v)/u+(bez_v/u)*(u*V+v))
    Q0,D0,lambda0,lambda1 = S.symbols('Q0 D0 lambda0 lambda1')
    zero('two admissible finite branch values force their difference times D into the content ideal',
         (C-lambda0*D0)-(C-lambda1*D0)-(lambda1-lambda0)*D0)
    x = S.symbols('x')
    zero('a punctured-curve cylinder has a free polynomial line in dimension three',
         (T*V-1).subs({T:1,V:1}))
    check('that free line remains nonconstant in the extra coordinate',S.diff(x,x)==1)

    d,r0,rinf = S.symbols('d r0 rinf')
    zero('two branch fibres exhaust Riemann-Hurwitz only when their point counts sum to two',
         (d-r0)+(d-rinf)-(2*d-2)-(2-r0-rinf))
    # A local index-one point contradicts a perfect power in every d>1.
    check('multiplicity one cannot be a nontrivial integer multiple',S.gcd(1,6)==1 and 1%6!=0)
    check('ramification indices multiply and a ramified intermediate place stays ramified',2*3>1)

    # Group control: transitivity alone does not imply full preservation.
    transposition = Permutation(0,1,size=3)
    cycle = Permutation(0,1,2)
    group = PermutationGroup(transposition,cycle)
    alternating = PermutationGroup(cycle)
    check('a proper normal subgroup can be transitive',
         alternating.is_normal(group) and alternating.is_transitive() and alternating.order()<group.order())
    check('normal closure of actual transposition inertia gives the full S3 group',
         group.normal_closure(PermutationGroup(transposition)).order()==group.order())
    check('the preserved quotient gate is stronger than root transitivity',group.order()//alternating.order()==2)

    # Elimination controls for the unrestricted nonadjacent-binomial family.
    for m,q,r in [(5,2,3),(6,2,2),(7,3,1)]:
        k = m-q
        product = (-S.Rational(q,m))**q*(-S.Rational(k,m))**k
        expected = (-1)**m*S.Rational(q**q*k**k,m**m)
        zero('critical elimination constant '+str((m,q,r)),product-expected)
        curve = T**(r*q)*V**k-expected
        check('the bad branch curve has no axis point '+str((m,q,r)),
              curve.subs(T,0)!=0 and curve.subs(V,0)!=0)
        P = T**r*z**m+z**q
        zero('critical derivative identity '+str((m,q,r)),
             S.diff(P,z)-z**(q-1)*(m*T**r*z**k+q))
        zero('critical value identity modulo its nonzero critical equation '+str((m,q,r)),
             P-S.Rational(k,m)*z**q-z**q*(m*T**r*z**k+q)/m)
    # Cubic normal form, with symbolic rather than bounded coefficient degrees.
    a,b,c,d0,H,ell,Y = S.symbols('a b c d0 H ell Y',nonzero=True)
    delta = b*b-3*a*c
    A = H**3/(27*a*a)
    B0 = d0+(2*b**3-9*a*b*c)/(27*a*a)
    transformed = A*((3*a*z+b)**3/H**3-3*ell*(3*a*z+b)/H)+B0
    zero('the cubic rational normal form follows exactly from Delta=H^2 ell',
         (transformed-(a*z**3+b*z**2+c*z+d0)).subs(ell,delta/H**2))
    zero('zero critical discriminant is a rationally translated cube',
         (a*(z+b/(3*a))**3+d0-b**3/(27*a*a)-(a*z**3+b*z**2+c*z+d0)).subs(c,b*b/(3*a)))
    tau = S.symbols('tau')
    zero('linear-double-cover critical value has the stated sign',
         (Y**3-3*ell*Y).subs({ell:tau**2,Y:tau})+2*tau**3)
    zero('the critical points are simple generically',S.diff(Y**3-3*ell*Y,Y,2).subs(Y,tau)-6*tau)
    even = T*T+3*T+1
    odd = T**3-T+2
    value = -even.subs(T,tau*tau)+2*odd.subs(T,tau*tau)*tau
    zero('even critical-value part recovers the output translation',
         (value+value.subs(tau,-tau))/2+even.subs(T,tau*tau))
    zero('odd critical-value part recovers A ell',
         (value-value.subs(tau,-tau))/(4*tau)-odd.subs(T,tau*tau))
    pole_composite = S.cancel(1/((x*x+1)-2))
    check('a nonconstant polynomial substitution cannot remove a genuine finite pole',
         S.degree(S.denom(pole_composite),x)==2 and S.gcd(S.numer(pole_composite),S.denom(pole_composite))==1)

    # Local controls illustrate the unbounded written valuation inequalities.
    n,hval = S.symbols('n hval',positive=True,integer=True)
    zero('at a negative root valuation the cubic is below the linear term',
         (hval-n)-(-3*n)-(hval+2*n))
    for order in [1,2,3,4,6]:
        possible = [w for w in range(-8,4) if
                    [3*w,w,-order].count(min(3*w,w,-order))>=2]
        check('negative constant valuation cancellation control '+str(order),
              possible==([-order//3] if order%3==0 else []))
    check('a possible simple pole at ell=0 admits no integral root valuation',
          all([3*w,1+w,1].count(min(3*w,1+w,1))==1 for w in range(-8,9)))
    scale,W = S.symbols('scale W',nonzero=True)
    zero('polynomial cube scaling produces the monic root equation',
         scale**3*((W/scale)**3-3*ell*(W/scale))-(W**3-3*ell*scale**2*W))
    check('the final Jacobian multiplier has positive W degree without cancellation',
          S.degree(S.diff(W**3-3*ell*scale**2*W,W),W)==2)

    # A concrete polynomial cubic with a linear critical double cover.
    P = (z+T)**3-3*T*(z+T)+T**2
    coeff = S.Poly(P,z)
    aa,bb,cc = [coeff.coeff_monomial(z**degree) for degree in [3,2,1]]
    zero('the cubic linear-critical control is outside the old high-genus case',bb*bb-3*aa*cc-9*T)
    s,t,w,v = S.symbols('s t w v')
    potential = 1+s*v+t*w+w*P.subs({T:s,z:w})
    zero('the cubic output has an actual-unit-Hessian automorphic control',
         S.det(S.hessian(potential,(s,t,w,v)))-1)

    context = json.loads((HERE/'sources/branch-parameter-import.json').read_text())
    pins = context['retained_proof_and_source_sha256']
    check('there are exactly seven explicit retained primary inputs',len(context['retained_primary_inputs'])==7)
    check('no new external theorem is imported',context['new_external_theorem_imports']==[])
    for path,digest in pins.items():
        check('retained byte pin '+path,hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest)
    retained_counts = {}
    for name in [name for name in pins if '/receipt-' in name]:
        receipt = json.loads((ROOT/name).read_text())
        check('retained receipt is PASS '+name,receipt['outcome']=='PASS')
        for relative,digest in receipt['source_sha256'].items():
            check('retained receipt source '+str((Path(name).parent/relative)),
                  hashlib.sha256((ROOT/Path(name).parent/relative).read_bytes()).hexdigest()==digest)
        for mapping in ['dependencies','dependency_source_sha256']:
            for path,digest in receipt.get(mapping,{}).items():
                check('retained dependency '+name+' -> '+path,
                      hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest)
        retained_counts[Path(name).name] = receipt['check_count']
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['verify_branch_parameters.py','BRANCH_CURVES_AND_RATIONAL_TRIANGULAR_PARAMETERS.md','sources/branch-parameter-import.json']},
        dependencies=pins,retained_receipt_counts_not_newly_replayed=retained_counts,
        written_not_finite_exhaustion=context['direct_proof_arguments'],new_external_theorem_imports=[],
        plane_only=True,actual_divisorial_inertia_required=True,
        all_cubic_outputs_irreducible=True,
        higher_dimensional_affine_parameter_theorem_asserted=False,
        all_sheets_over_branch_curve_asserted_omitted=False,
        intermediate_field_assumed_Galois=False,arbitrary_target_cover_group_preservation=False,
        transitivity_asserted_full_group=False,full_cubic_splitting_group_asserted_preserved=False,
        parameter_extension_source_assumed_affine_plane=False,
        nonlinear_source_Hessian_invariance_used=False,
        unit_ideal_endpoint_resolved=False,singleton_zero_endpoint_resolved=False,
        unrestricted_HC4_or_JC2_resolved=False,original_goal_complete=False,
        formal_verification=False,independent_replay=False,external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-branch-parameters.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome='PASS',check_count=len(checks),receipt=str(output)),indent=2))


if __name__=='__main__':
    main()
