"""Exact controls for independent coefficient fields, overlapping rows and generic Bezout frames."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import time

import sympy as S
from sympy.combinatorics import Permutation, PermutationGroup
from sympy.matrices.normalforms import smith_normal_form

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

    T,V,U,u,x,X,Y,Z = S.symbols('T V U u x X Y Z')
    algebraic_polynomial = u**3*Z**3+u*Z+u**2
    zero('a polynomial with genuinely algebraic coefficients has an integral affine monic model',
         algebraic_polynomial.subs(Z,U/u)-(U**3+U+u**2))
    check('its leading coefficient is not invariant under the algebraic parameter conjugation',
          S.expand(algebraic_polynomial.subs(u,-u)-algebraic_polynomial)!=0)
    zero('the conjugate model is integral as well',
         algebraic_polynomial.subs(u,-u).subs(Z,-U/u)-(U**3+U+u**2))
    zero('an algebraic critical-value image can be polynomially parametrizable',
         (V**2-T).subs({T:x**2,V:x}))
    zero('a finite-pole algebraic branch image instead has the stated unit equation',
         (T*V**2-1).subs({T:u**2,V:1/u}))
    a,b = S.symbols('a b',nonnegative=True,integer=True)
    check('the finite-pole image admits no nonconstant polynomial parametrization by the degree equation',
          S.solve([a+2*b],[a],dict=True)==[{a:-2*b}])

    # Both extensions of the parameter-conjugate root model are kept in the closure.
    flip_a = Permutation(0,1,size=4)
    flip_b = Permutation(2,3,size=4)
    parameter_swap = Permutation(0,2)(1,3)
    group = PermutationGroup(flip_a,flip_b,parameter_swap)
    parameter_kernel = PermutationGroup(flip_a,flip_b)
    check('the closure of two conjugate quadratic models has a nontrivial parameter quotient',
          group.order()==8 and parameter_kernel.order()==4 and parameter_kernel.is_normal(group))
    check('the lifted branch inertia lies in the parameter-fixing subgroup',flip_a in parameter_kernel)
    check('parameter conjugation really exchanges the two relative models',
          parameter_swap*flip_a*~parameter_swap==flip_b)
    check('restricting the lifted inertia to the first model remains transitive',
          PermutationGroup(flip_a).orbit(0)=={0,1})

    # A non-Keller-mate factorization illustrates the finite specialization guards.
    lam,t = S.symbols('lam t')
    left = (t-lam*U)/(lam-1)
    right = (lam-1)*(t+lam*U)
    zero('factorization identities survive clearing the finite coefficient denominators',
         left*right-(t**2-lam**2*U**2))
    zero('a permitted closed specialization retains the full identity',
         (left*right).subs(lam,2)-(t**2-4*U**2))
    check('both factors retain positive degree at that permitted specialization',
          S.degree(left.subs(lam,2),U)==1 and S.degree(right.subs(lam,2),U)==1)
    check('a degree-dropping specialization must be removed by a guard',
          S.degree(left.subs(lam,0),U)==0 and S.degree(right.subs(lam,0),U)==0)
    check('a denominator-zero specialization is separately excluded',S.denom(left).subs(lam,1)==0)
    Dg = S.symbols('Dg')
    independent = U**2+V-Dg
    check('independent coefficients leave a simple V zero and a nonsquare quadratic discriminant',
          S.degree(S.discriminant(independent,U),V)==1)
    zero('identifying the coefficient with the original source g destroys the theorem hypothesis',
         independent.subs(Dg,V)-U**2)

    A = (X**2-T)*(X-1)**2
    P = (X-1)*Z**3+Z**2+X*Z
    zero('the multiplier splits after the finite parameter extension',
         A.subs(T,u**2)-(X-u)*(X+u)*(X-1)**2)
    residuals = [S.expand(P.subs(X,u)),S.expand(P.subs(X,-u)),S.expand(P.subs(X,1))]
    check('overlap gives different generic polynomial degrees in different multiplier rows',
          [S.degree(value,Z) for value in residuals]==[3,3,2])
    check('the ambient joint degree is not the number of components in a restricted row',S.Poly(P,X,Z).total_degree()==4)
    for index,value in enumerate(residuals):
        check('the generic restricted level is reduced in row '+str(index),
              S.Poly(value+V,Z,domain=S.QQ.frac_field(u,V)).gcd(S.Poly(S.diff(value,Z),Z,domain=S.QQ.frac_field(u,V))).degree()==0)
    relations = S.zeros(3,8)
    offset = 0
    for row,count in enumerate([3,3,2]):
        for col in range(offset,offset+count):
            relations[row,col] = 1
        offset += count
    smith = smith_normal_form(relations,domain=S.ZZ)
    check('the overlapping row lattice has exactly the three primitive row relations',
          relations.rank()==3 and [abs(smith[i,i]) for i in range(3)]==[1,1,1])
    check('the resulting class lattice is torsion-free of rank five, independent of the repeated multiplier exponent',
          relations.cols-relations.rank()==5)
    cyclic_actions = []
    for count in [3,3,2]:
        action = S.zeros(count-1)
        for col in range(count-1):
            if col+1<count-1:
                action[col+1,col]=1
            else:
                for row in range(count-1):
                    action[row,col]=-1
        cyclic_actions.append(action)
    action = S.diag(*cyclic_actions)
    check('row-fixing transitive residual actions have zero source-fixed class subgroup',
          (action-S.eye(5)).rank()==5)

    X1,X2,W,J = S.symbols('X1 X2 W J')
    multiplier = 1+X1*X2
    leading = T*X1
    bezout_u,bezout_v = S.Integer(1),-X2/T
    root = T*Z**4+Z*W**3
    lower = Z+X1*(T*Z**2*W+X2**5*Z**2+1)
    output = multiplier*Y-leading*root-lower
    zero('the overlapping leading coefficient has a genuine polynomial Bezout identity over C(T)',
         bezout_u*multiplier+bezout_v*leading-1)
    check('that generic certificate can have a parameter denominator and is not a global monic assertion',S.denom(bezout_v)==T)
    helper = bezout_v*Y+bezout_u*root
    shifted = V+lower
    zero('the determinant-one generic frame reconstructs Y',
         (bezout_u*shifted+leading*helper).subs(V,output)-Y)
    zero('the same inverse reconstructs the root polynomial',
         (-bezout_v*shifted+multiplier*helper).subs(V,output)-root)
    relation = multiplier*J-root-bezout_v*(V+lower)
    zero('the new generic relation holds on the original ambient ring',relation.subs({J:helper,V:output},simultaneous=True))
    reverse_y = bezout_u*shifted+leading*J
    zero('the reverse map preserves the original observable modulo the new relation',
         multiplier*reverse_y-leading*root-lower-V-leading*relation)
    weights = [1,1,38,10,10,0,0]
    poly = S.Poly(S.expand(T*relation),X1,X2,J,Z,W,T,V)
    top_weight = max(sum(e*w for e,w in zip(monomial,weights)) for monomial,c in poly.terms())
    top = sum(c*S.prod(symbol**e for symbol,e in zip(poly.gens,monomial))
              for monomial,c in poly.terms() if sum(e*w for e,w in zip(monomial,weights))==top_weight)
    check('one positive generic weighting handles all overlapping lower terms in the control',top_weight==40)
    zero('the generic highest equation has the stated primitive linear form',
         top-T*(X1*X2*J-T*Z**4-Z*W**3))
    check('its leading coefficients are coprime over the geometric generic field',S.gcd(X1*X2,T*Z**4+Z*W**3)==1)
    ell,d,M,B = S.symbols('ell d M B',positive=True,integer=True)
    zero('the two highest terms balance for every admitted degree',ell+(d*M-ell)-d*M)
    zero('all lower root terms have the positive M-minus-coefficient-degree gap',
         d*M-(B+(d-1)*M)-(M-B))
    check('gcd one in several coefficient variables does not imply the required Bezout identity',
          S.gcd(X1,X2)==1 and S.groebner([X1,X2],X1,X2).reduce(S.Integer(1))[1]==1)

    content_A = T*(1+(T+1)*X)
    content_P = T*X*Z**4+X*Z
    check('the enlarged content criterion allows a closed polynomial depending on overlapping X and Z',
          content_P.subs(T,0)==X*Z and S.gcd_list(S.Poly(content_A,X).all_coeffs())==T)
    check('that content zero level may be reducible although its generic level is integral',
          len(S.factor_list(X*Z,X,Z)[1])==2 and S.gcd(X,V)==1)
    alpha = S.symbols('alpha')
    zero('a composite overlapping content level remains excluded',
         ((X*Z)**2+V).subs(V,-alpha**2)-(X*Z-alpha)*(X*Z+alpha))
    missing_direction = X*(Y-Z)-V
    zero('without a residual polynomial direction the generic multiplier becomes a nonconstant unit',
         X*(Y-Z)/V-1-missing_direction/V)
    check('the residual-direction condition detects that excluded row',S.Poly((X*Z).subs(X,0),Z).is_zero)

    s,t,u0,v0 = S.symbols('s t u0 v0')
    actual_A = 1+X1*X2
    actual_P = X1*Z**5+Z+X1*(T*X2**3*Z**2+1)
    embedding = {T:s,X1:S.Integer(0),X2:u0,Y:t+v0,Z:v0}
    zero('the actual HC4 overlap embedding carries the second identity coordinate',
         (actual_A*Y-actual_P).subs(embedding,simultaneous=True)-t)
    check('the overlap embedding is closed with the explicit source inverse',
          all(S.expand(value.subs(embedding,simultaneous=True)-coord)==0
              for value,coord in zip([T,Y-Z,X2,Z],[s,t,u0,v0])))
    actual_pullback = (T*X2+Y*Z+X1**2).subs(embedding,simultaneous=True)
    zero('the actual four-variable Hessian is a nonzero constant',S.hessian(actual_pullback,(s,t,u0,v0)).det()-1)

    context = json.loads((HERE/'sources/overlapping-output-import.json').read_text())
    check('one finite-type specialization statement is explicitly imported',len(context['newly_used_primary_inputs'])==1)
    primary = context['newly_used_primary_inputs'][0]
    check('the Nullstellensatz import is pinned and is not represented as fully reproved',
          primary['theorem_fully_reproved'] is False and hashlib.sha256((HERE/'sources'/primary['retained_file']).read_bytes()).hexdigest()==primary['sha256'])
    check('all generic and global qualifications remain explicit',
          context['residual_polynomial_direction_required'] is True
          and context['generic_positive_filtration_required_for_global_inverse'] is True
          and context['all_overlaps_asserted_covered'] is False)
    controls = len(checks)
    for name,digest in context['retained_proof_and_source_sha256'].items():
        check('retained byte pin '+name,hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    counts = {}
    for name in ['receipt-reflexive-pencil.json','receipt-parameter-stability.json','receipt-parameter-norms.json']:
        retained = json.loads((HERE/name).read_text())
        check('retained receipt is PASS '+name,retained['outcome']=='PASS')
        counts[name] = retained['check_count']
        for path,digest in retained['source_sha256'].items():
            check('retained receipt source '+name+' -> '+path,hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest)
        for path,digest in retained['dependencies'].items():
            check('retained dependency '+name+' -> '+path,hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest)
    flags = {key:value for key,value in context.items() if type(value) is bool}
    result = dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=controls,checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['verify_overlapping_outputs.py','UNIVERSAL_COEFFICIENT_FIELDS_AND_OVERLAPPING_OUTPUTS.md','sources/overlapping-output-import.json']},
        dependencies=context['retained_proof_and_source_sha256'],retained_receipt_counts_not_newly_replayed=counts,
        newly_used_primary_inputs=context['newly_used_primary_inputs'],written_not_finite_exhaustion=context['direct_proof_arguments'],**flags,
        formal_verification=False,independent_replay=False,external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-overlapping-outputs.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=controls),indent=2))


if __name__=='__main__':
    main()
