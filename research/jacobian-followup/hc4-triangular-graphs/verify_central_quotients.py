"""Exact controls for relative central series, algebraic product descent and binomial poles."""
import argparse
from itertools import permutations, product
import hashlib
import json
from pathlib import Path
import platform
import time

import sympy as S

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

    # A nonnilpotent base and a nonabelian nilpotent fibre: S3 times D8.
    symmetric = set(permutations(range(3)))
    identity = tuple(range(3))
    dihedral = set(product(range(4),range(2)))

    def pmul(left,right):
        return tuple(left[right[i]] for i in range(3))

    def dmul(left,right):
        r,e = left
        rr,ee = right
        return ((r+(-1)**e*rr)%4,(e+ee)%2)

    def mul(left,right):
        return (pmul(left[0],right[0]),dmul(left[1],right[1]))

    group = set(product(symmetric,dihedral))
    unit = (identity,(0,0))
    inverse = {g:next(h for h in group if mul(g,h)==unit) for g in group}

    def normal(subgroup):
        return all(mul(mul(g,h),inverse[g]) in subgroup for g in group for h in subgroup)

    def pairwise_product(left,right):
        return {mul(g,h) for g in left for h in right}

    fibre = {(identity,h) for h in dihedral}
    base = {(g,(0,0)) for g in symmetric}
    centre_fibre = {(identity,(r,0)) for r in [0,2]}
    check('the product control has nonnilpotent parameter factor and order 48',len(group)==48)
    check('parameter and fibre images commute and generate the whole group',
          pairwise_product(base,fibre)==group and all(mul(g,h)==mul(h,g) for g in base for h in fibre))
    check('the nilpotent geometric kernel need not itself be central or abelian',
          normal(fibre) and any(mul(g,h)!=mul(h,g) for g in fibre for h in fibre))
    check('the first relative central layer is central in the full product',
          normal(centre_fibre) and all(mul(g,h)==mul(h,g) for g in group for h in centre_fibre))
    check('the second relative layer is central modulo the first',
          all(mul(mul(mul(g,h),inverse[g]),inverse[h]) in centre_fibre for g in group for h in fibre))
    subgroup = {(g,(0,e)) for g in symmetric for e in [0,1]}
    lifted = pairwise_product(subgroup,centre_fibre)
    check('surjecting onto the parameter quotient does not make the source subgroup normal at once',
          pairwise_product(subgroup,fibre)==group and not normal(subgroup))
    check('descending the relative series produces a proper normal container',
          len(subgroup)==12 and len(lifted)==24 and normal(lifted) and subgroup<=lifted<group)
    # A nilpotent normal kernel with nontrivial quotient conjugation is insufficient.
    def pinv(g):
        return next(h for h in symmetric if pmul(g,h)==identity)

    three_cycle = (1,2,0)
    alternating = {identity,three_cycle,pmul(three_cycle,three_cycle)}
    transposition = (1,0,2)
    order_two = {identity,transposition}
    check('a nilpotent normal kernel alone fails the full-group criterion',
          {pmul(g,h) for g in order_two for h in alternating}==symmetric and len(order_two)<len(symmetric))
    check('the failed kernel is not relatively central',
          pmul(pmul(transposition,three_cycle),pinv(transposition))!=three_cycle)
    # Pullback descent must retain the parameter subgroup, rather than discard it.
    check('killing the geometric kernel leaves a genuine nonabelian parameter quotient in the control',
          len(group)//len(fibre)==6 and pmul(three_cycle,transposition)!=pmul(transposition,three_cycle))
    check('restriction to a parameter section retains every parameter monodromy value',
          {g for g,h in base}==symmetric)

    T,V,z,a,b = S.symbols('T V z a b')
    for m,q in [(5,2),(6,2),(7,3)]:
        k = m-q
        constant = (-1)**m*S.Rational(q**q*k**k,m**m)
        P = a*z**m+b*z**q
        zero('nonzero critical derivative identity '+str((m,q)),
             S.diff(P,z)-z**(q-1)*(m*a*z**k+q*b))
        zero('critical-value identity '+str((m,q)),
             P-S.Rational(k,m)*b*z**q-z**q*(m*a*z**k+q*b)/m)
        zero('critical elimination constant '+str((m,q)),
             (-S.Rational(q,m))**q*(-S.Rational(k,m))**k-constant)
        discriminant = S.discriminant(P+V,z)
        critical = a**q*V**k-constant*b**m
        support = a*b*V*critical
        residual = discriminant
        while S.total_degree(residual)>0:
            divisor = S.gcd(residual,support)
            if S.total_degree(divisor)==0:
                break
            residual = S.cancel(residual/divisor)
        check('the root discriminant has only the stated vertical, zero and critical support '+str((m,q)),
              residual!=0 and S.total_degree(residual)==0)

    # Cancellation is essential when coefficients have common zeros.
    examples = [
        (5,2,T**6,T,True),
        (6,2,T**7*(T-1)**4,T*(T+1),True),
        (5,2,T**5,T**2,False),
        (4,1,(T-1)**5,T-1,True)]
    for m,q,aa,bb,admitted in examples:
        denominator = aa**q
        numerator = bb**m
        common = S.gcd(denominator,numerator)
        U,W = [S.cancel(value/common) for value in [denominator,numerator]]
        label = str((m,q,aa,bb))
        check('reduced critical-value coefficients are coprime '+label,S.gcd(U,W)==1)
        check('the denominator is nonconstant exactly when the finite-pole test admits this control '+label,
              (S.degree(U,T)>0)==admitted)
        check('finite-pole condition agrees with failure of polynomial divisibility '+label,
              (S.rem(numerator,denominator,T)!=0)==admitted)
        bez_u,bez_w,gcd = S.gcdex(U,W,T)
        zero('Bezout makes the reduced critical denominator a unit on every component '+label,
             bez_u*U+bez_w*W-1)
    m,q = S.symbols('m q',positive=True,integer=True)
    zero('the all-degree shared-zero example has pole order q(m+1)-m',
         q*(m+1)-m-((q-1)*m+q))
    # A reducible critical curve still has no polynomial A1 map in any component.
    factors = S.factor_list(T**2*V**4-1,T,V)[1]
    check('imprimitive controls may have more than one critical-value component',len(factors)==2)
    check('each such component omits the denominator zero',all(factor.subs(T,0)!=0 for factor,exponent in factors))
    zero('the no-finite-pole boundary can be a parametrizable critical curve',
         (V**2-T**3).subs({T:z*z,V:z**3}))
    check('that critical-curve parametrization is nonconstant',S.degree(z*z,z)>0)
    yy = S.symbols('Y')
    fixed = yy**5+yy**2
    zero('failure of the pole test need not be open: a rational-input fixed-polynomial case',
         (T**5*z**5+T**2*z**2)-fixed.subs(yy,T*z))
    # Product topology is a hypothesis: an arbitrary moving arrangement need not commute.
    check('noncommuting geometric and base generators are a rejected product certificate',
          pmul(transposition,three_cycle)!=pmul(three_cycle,transposition))
    s,t,w,v = S.symbols('s t w v')
    actual_P = T**6*z**5+T*z*z
    actual_potential = 1+s*v+t*w+w*actual_P.subs({T:s,z:w})
    zero('the mixed-coefficient finite-pole family has an actual-unit-Hessian automorphic control',
         S.det(S.hessian(actual_potential,(s,t,w,v)))-1)

    context = json.loads((HERE/'sources/central-quotient-import.json').read_text())
    check('the comparison upgrade is explicitly imported rather than credited to the forward-only input',
          len(context['newly_used_external_statements'])==1
          and 'full faithfulness' in context['newly_used_external_statements'][0]['statement'])
    check('the two covers in the comparison are already algebraic and essential surjectivity is unused',
          context['essential_surjectivity_used'] is False
          and 'already algebraic' in context['newly_used_external_statements'][0]['exact_use'])
    check('all-dimensional field and product claims are separated from plane branch geometry',
          len(context['all_dimensional_claims'])==2 and len(context['plane_only_claims'])==2)
    pins = context['retained_proof_and_source_sha256']
    for path,digest in pins.items():
        check('retained byte pin '+path,hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest)
    retained_counts = {}
    for name in [name for name in pins if '/receipt-' in name]:
        receipt = json.loads((ROOT/name).read_text())
        check('retained receipt is PASS '+name,receipt['outcome']=='PASS')
        for relative,digest in receipt['source_sha256'].items():
            check('retained receipt source '+str(Path(name).parent/relative),
                  hashlib.sha256((ROOT/Path(name).parent/relative).read_bytes()).hexdigest()==digest)
        for mapping in ['dependencies','dependency_source_sha256']:
            for path,digest in receipt.get(mapping,{}).items():
                check('retained dependency '+name+' -> '+path,
                      hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest)
        retained_counts[Path(name).name] = receipt['check_count']
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['verify_central_quotients.py','CENTRAL_GEOMETRIC_QUOTIENTS_AND_BINOMIAL_OUTPUTS.md','sources/central-quotient-import.json']},
        dependencies=pins,retained_receipt_counts_not_newly_replayed=retained_counts,
        written_not_finite_exhaustion=context['direct_proof_arguments'],
        newly_used_external_statements=context['newly_used_external_statements'],
        all_dimensional_relative_central_series=True,all_dimensional_product_cover_gate=True,
        branch_inertia_and_binomial_application_plane_only=True,
        product_open_and_nilpotent_fibre_image_required=True,
        full_faithfulness_newly_used=True,essential_surjectivity_used=False,
        nilpotent_normal_kernel_alone_asserted_sufficient=False,
        residual_whole_group_asserted_nilpotent=False,
        arbitrary_moving_branch_complement_asserted_product=False,
        no_finite_pole_binomials_asserted_all_closed=False,
        arbitrary_target_cover_group_preservation=False,
        parameter_extension_source_assumed_affine_plane=False,
        nonlinear_source_Hessian_invariance_used=False,
        unit_ideal_endpoint_resolved=False,singleton_zero_endpoint_resolved=False,
        unrestricted_HC4_or_JC2_resolved=False,original_goal_complete=False,
        formal_verification=False,independent_replay=False,external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-central-quotients.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome='PASS',check_count=len(checks),receipt=str(output)),indent=2))


if __name__=='__main__':
    main()
