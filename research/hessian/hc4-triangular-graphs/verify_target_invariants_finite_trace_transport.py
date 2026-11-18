#!/usr/bin/env python3
"""Exact controls for target invariants and finite trace transport.

The universal proof is written separately. Branched cyclic controls are
finite algebras, not actual noninverse plane Keller normalizations.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
import re
import time
from pathlib import Path
import sympy as S

BASE=Path(__file__).resolve().parent
PROOF=BASE/'TARGET_INVARIANTS_AND_FINITE_TRACE_TRANSPORT_OBSTRUCTION.md'
RECEIPT=BASE/'target_invariants_finite_trace_transport_verification.json'
CHECKS=[]

def check(name,condition):
    ok=bool(condition)
    CHECKS.append({'name':name,'passed':ok})
    if not ok:
        raise AssertionError(name)

def eq(name,lhs,rhs=0):
    check(name,S.cancel(S.sympify(lhs)-S.sympify(rhs))==0)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def cyclic_trace(power,e,h):
    return e*h**(power//e) if power%e==0 else S.Integer(0)

def cyclic_gram(e,h):
    return S.Matrix(e,e,lambda i,j:cyclic_trace(i+j,e,h))

def run():
    x,y,u,p,q,z=S.symbols('x y u p q z')
    # Minimality and reduced pullback have distinct roles.
    for n in range(2,12):
        h=y**n
        target=q
        check(f'critical fibre repeated in nonetale example degree{n}',
              S.degree(S.gcd(h,S.diff(h,y)),y)==n-1)
        eq(f'nonetale example Jacobian degree{n}',
           S.det(S.Matrix([[S.diff(x,x),S.diff(x,y)],
                           [S.diff(h,x),S.diff(h,y)]])),n*y**(n-1))
        # Identity source: choosing h=y^n is not least degree.
        check(f'target power has a lower invariant factor degree{n}',
              S.degree(q**n,q)>S.degree(q,q)>0)
        check(f'repeated single target factor is not reduced degree{n}',
              S.gcd(q**n,S.diff(q**n,q))==q**(n-1))
    for h in [y*y+y,y**3-3*y,y**4+y*y,y**5-y]:
        derivative=S.diff(h,y)
        roots=S.roots(derivative,y)
        check('every nonlinear polynomial has a critical root '+str(h),
              sum(roots.values())==S.degree(derivative,y)>0)
        for root in roots:
            shifted=S.expand(h-h.subs(y,root))
            eq('critical shifted value vanishes '+str((h,root)),shifted.subs(y,root))
            eq('critical shifted derivative vanishes '+str((h,root)),
               S.diff(shifted,y).subs(y,root))
    # Automorphisms are actual positive Keller controls, with reduced pulls.
    P=x+y*y
    Q=y+P*P
    J=S.diff(P,x)*S.diff(Q,y)-S.diff(P,y)*S.diff(Q,x)
    eq('actual nonlinear positive Keller control Jacobian',J,1)
    y_inverse=q-p*p
    eq('source invariant is actually target polynomial',
       y_inverse.subs({p:P,q:Q}),y)
    eq('invariant implies first velocity in target P',S.diff(P,x),
       S.diff(y_inverse,q).subs({p:P,q:Q}))
    eq('invariant implies first velocity in target Q',S.diff(Q,x),
       -S.diff(y_inverse,p).subs({p:P,q:Q}))
    for H in [p,q,p*q,p*p-q,q*q-p,p*p+q*q-1]:
        pull=S.expand(H.subs({p:P,q:Q}))
        check('reduced actual Keller pullback '+str(H),
              S.gcd(S.gcd(pull,S.diff(pull,x)),S.diff(pull,y))==1)
    # Degree comparison for all represented source replicas and generators.
    for degree in range(1,9):
        f=(y+1)*x**degree+x**(degree-1)+y
        for replica in [S.Integer(1),y,y*y+1,(y-1)**3]:
            df=replica*S.diff(f,x)
            eq(f'replica lowers x degree d{degree}q{replica}',
               S.degree(df,x),degree-1)
            check(f'noninvariant principal ideal not stable d{degree}q{replica}',
                  S.rem(df,f,x)!=0)
    for h in [y,y*y+1,(y-1)**3]:
        eq('invariant ideal genuinely stable '+str(h),S.diff(h,x),0)
    # Cyclic finite algebra Gram and connection, exact in the coefficient field.
    for e in range(2,9):
        for h in [x+y*y,y-1,x*y]:
            G=cyclic_gram(e,h)
            derivative=S.diff(h,x)
            M=S.diag(*[S.Rational(j,e)*derivative/h for j in range(e)])
            check(f'cyclic trace connection matrix e{e}h{h}',
                  G.diff(x)==(M*G+G*M.T).applyfunc(S.cancel))
            disc=S.factor(G.det())
            expected=(-1)**((e-1)*(e-2)//2)*e**e*h**(e-1)
            eq(f'cyclic discriminant e{e}h{h}',disc,expected)
            eq(f'determinant connection identity e{e}h{h}',
               S.diff(disc,x),2*S.trace(M)*disc)
            eq(f'cyclic connection trace e{e}h{h}',S.trace(M),
               S.Rational(e-1,2)*derivative/h)
            if h==y-1:
                check(f'invariant branch admits regular replica lift e{e}',M==S.zeros(e))
            else:
                check(f'moving branch has polar regular-lift obstruction e{e}h{h}',
                      not S.cancel(M[1,1]).is_polynomial(x,y))
            # Nonconstant local basis changes test the derivative correction.
            U=S.eye(e)
            U[1,0]=x+y
            if e>2:
                U[2,0]=x*x
                U[2,1]=y*x
            G2=U*G*U.T
            M2=(U.diff(x)+U*M)*U.inv()
            check(f'changed basis trace connection e{e}h{h}',
                  (G2.diff(x)-M2*G2-G2*M2.T).applyfunc(S.cancel)==S.zeros(e))
            eq(f'changed basis determinant retains discriminant e{e}h{h}',G2.det(),disc)
        # Without LND the discriminant can be stable and noninvariant.
        h=x
        disc=S.factor(cyclic_gram(e,h).det())
        eq(f'nonLND Euler lift preserves nonconstant discriminant e{e}',
           x*S.diff(disc,x),(e-1)*disc)
        eq(f'nonLND Euler regular lift relation e{e}',
           e*u**(e-1)*(u/e),u**e)
    # Cyclic trace-defect modules have the REDUCED annihilator h,
    # although their discriminant has exponent e-1. A connection must
    # preserve every defining relation h*generator=0.
    for e in range(2,10):
        for h in [x+y*y,y-1,x*y]:
            G=cyclic_gram(e,h)
            check(f'cyclic trace defect constant unit summand e{e}h{h}',
                  G[0,0]==e and all(G[0,j]==G[j,0]==0 for j in range(1,e)))
            reduced=(G[1:,1:]/h).applyfunc(S.cancel)
            check(f'cyclic reduced support block invertible e{e}h{h}',
                  reduced.det()!=0 and not reduced.det().has(x,y))
            for replica in [S.Integer(1),y,y*y+1]:
                dh=replica*S.diff(h,x)
                relation=S.rem(dh,h,x) if S.degree(h,x)>0 else S.rem(dh,h,y)
                if h==y-1:
                    eq(f'invariant cyclic module connection relation e{e}q{replica}',relation,0)
                else:
                    check(f'moving cyclic module connection obstructed e{e}h{h}q{replica}',
                          relation!=0)
    # Nontrivial reduced support with two independently supported packets.
    h1=x+y*y;h2=x-y
    support=h1*h2
    eq('packet annihilator is reduced lcm',S.lcm(h1,h2),support)
    check('packet connection relation fails at first component',
          S.rem(S.diff(support,x),h1,x)!=0)
    check('packet connection relation fails at second component',
          S.rem(S.diff(support,x),h2,x)!=0)
    # General matrix identity, independently from cyclic formulas.
    for size in range(2,6):
        G=S.diag(*[x+j+1 for j in range(size)])
        U=S.eye(size)
        U[1,0]=x*y
        G2=U*G*U.T
        M=(G2.diff(x)*G2.inv())/2
        check(f'generic symmetric connection size{size}',
              (G2.diff(x)-M*G2-G2*M.T).applyfunc(S.cancel)==S.zeros(size))
        eq(f'generic symmetric determinant size{size}',
           S.diff(G2.det(),x),2*S.trace(M)*G2.det())
    # Finite nontrivial trace transport: Laurent Keller controls.
    for e in range(2,8):
        for k in range(-4,8):
            for ell in range(4):
                n=k+(e-1)*ell
                lhs=k*q**ell*cyclic_trace(n-1,e,p)
                af=e*q**ell*cyclic_trace(n+e-1,e,p)
                bf=-(e-1)*q**(ell+1)*cyclic_trace(n-1,e,p)
                rhs=S.diff(af,p)+S.diff(bf,q)
                eq(f'exact finite trace divergence e{e}k{k}ell{ell}',lhs,rhs)
        check(f'trace does not commute with source translation e{e}',
              cyclic_trace(0,e,p)!=0 and cyclic_trace(1,e,p)==0)
        # Missing u=0 is essential: u^e is a unit in the Laurent source.
        eq(f'Laurent critical fibre has unit inverse e{e}',u**e*u**(-e),1)
    # Trace divergence for actual positive polynomial Keller pairs.
    a=S.diff(P,x);b=S.diff(Q,x)
    DP=lambda f:S.diff(Q,y)*S.diff(f,x)-S.diff(Q,x)*S.diff(f,y)
    DQ=lambda f:-S.diff(P,y)*S.diff(f,x)+S.diff(P,x)*S.diff(f,y)
    eq('actual source translation has zero target divergence',DP(a)+DQ(b),0)
    for f in [1,x,y,x*y,x*x+y**3,P*Q]:
        eq('exact divergence identity actual positive control '+str(f),
           S.diff(f,x),DP(a*f)+DQ(b*f))
    return True

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write-receipt',action='store_true')
    args=parser.parse_args()
    start=time.monotonic()
    run()
    text=PROOF.read_text()
    for phrase in ['EVERY polynomial source frame','B intersect C[y]=C',
                   'R tensor_B (Lambda/S)','T=R times T_other',
                   'D_R(G)=M*G+G*M^t','The lift itself was NOT assumed locally nilpotent',
                   'The two requests are different','No automatic ordinary-codifferent membership',
                   'Tr(delta f)=delta Tr(f)','regular-base projective-dimension formula',
                   'Ann_R(M)=H(P,Q)*R','M admits NO D_R-connection',
                   'It need not be locally nilpotent']:
        check('written scope retains '+phrase,phrase in text)
    for link in re.findall(r'\]\(([^)]+)\)',text):
        if '://' not in link:
            check('local proof link exists '+link,(BASE/link).is_file())
    parents=['DIRECT_REPLICA_INDEX_FROM_A_DESCENDED_DIFFERENTIAL.md',
             'CANONICAL_BOUNDARY_AND_STABLE_COMPLETE_INTERSECTION_OBSTRUCTION.md',
             'TRACE_LINKING_HALF_RAMIFICATION_AND_CUBIC_MULTIPLICATION.md',
             'CODIFFERENT_PREIMAGES_UNDER_TARGET_DERIVATIONS.md',
             'PLINTH_IDEALS_AND_INVARIANT_REPLICA_EXTENSION.md',
             'REDUCED_BRANCH_PAIRING_AND_SIXTH_MOMENT_DEFECT.md',
             'LOGARITHMIC_TRACE_DEFECT_AND_NON_DENSE_DERIVATIVE_DOMAINS.md']
    receipt={'claim':'Target invariants and finite trace-transport obstruction',
             'outcome':'PASS','check_count':len(CHECKS),'checks':CHECKS,
             'proof_sha256':sha(PROOF),'checker_sha256':sha(Path(__file__)),
             'retained_inputs':[{'path':p,'sha256':sha(BASE/p)} for p in parents],
             'environment':{'python':platform.python_version(),'sympy':S.__version__},
             'scope':{'actual_target_polynomial_invariant_forces_inverse':True,
                      'finite_source_basechange_trace_and_diagonal_factor':True,
                      'regular_trace_connection_preserves_discriminant':True,
                      'noninverse_no_regular_source_replica_lift_to_basechange':True,
                      'source_basechanged_trace_defect_admits_no_replica_connection':True,
                      'lift_itself_assumed_LND':False,
                      'original_S_replica_extension_identified_with_basechange':False,
                      'automatic_codifferent_membership_proved':False,
                      'surviving_fold_family_excluded':False},
             'assurance':'Internal written proof and exact controls. Branched covers and Laurent maps are scoped controls, not plane Keller counterexamples. No formal verification or external human review.'}
    if args.write_receipt:
        RECEIPT.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'outcome':'PASS','check_count':len(CHECKS),
                      'elapsed_seconds':time.monotonic()-start,
                      'receipt_written':args.write_receipt},indent=2))

if __name__=='__main__':
    main()
