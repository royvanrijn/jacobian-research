#!/usr/bin/env python3
"""Elementary exact rational point-orbit and reciprocal-exponent decision.

The algorithm returns the entire nonnegative exponent set. It does not use
an Orbit oracle, floating point, or an unproved search cutoff. The written
proof establishes termination; these controls exercise every branch.
"""
from __future__ import annotations
if not __debug__:
    raise RuntimeError('Active assertions are required; do not use -O or -OO.')

import argparse
import hashlib
import json
from math import comb, gcd, lcm
from pathlib import Path
import platform
import sympy as sp

T = sp.Symbol('T')


def poly(expr):
    return sp.Poly(expr, T, domain=sp.QQ)


def encoded(p):
    return [str(c) for c in p.all_coeffs()]


def primitive(p):
    scale = lcm(*(int(c.q) for c in p.all_coeffs()))
    values = [int(c*scale) for c in p.all_coeffs()]
    content = 0
    for value in values:
        content = gcd(content,abs(value))
    values = [value//content for value in values]
    if values[0] < 0:
        values = [-value for value in values]
    return sp.Poly.from_list(values,T,domain=sp.ZZ)


def ppow(p,n,f):
    assert n >= 0 and f.degree() >= 1
    value = poly(1).rem(f)
    p = p.rem(f)
    while n:
        if n&1:
            value = (value*p).rem(f)
        p = (p*p).rem(f)
        n >>= 1
    return value


def multiplication_matrix(p,f):
    columns = []
    for i in range(f.degree()):
        v = (p*poly(T**i)).rem(f)
        columns.append(sp.Matrix([v.nth(j) for j in range(f.degree())]))
    return sp.Matrix.hstack(*columns)


def vector(p,f):
    p = p.rem(f)
    return sp.Matrix([p.nth(i) for i in range(f.degree())])


def empty():
    return {'kind':'empty'}


def singleton(j):
    return {'kind':'singleton','exponent':int(j)}


def progression(first,step):
    assert first >= 0 and step > 0
    return {'kind':'progression','first':int(first),'step':int(step)}


def contains(solution,j):
    if solution['kind']=='empty':
        return False
    if solution['kind']=='singleton':
        return j == solution['exponent']
    return j >= solution['first'] and (j-solution['first'])%solution['step']==0


def intersect_tail(solution,onset):
    if solution['kind']=='empty':
        return solution
    if solution['kind']=='singleton':
        return solution if solution['exponent'] >= onset else empty()
    first,step = solution['first'],solution['step']
    if first < onset:
        first += ((onset-first+step-1)//step)*step
    return progression(first,step)


def torsion_order(h):
    """h is irreducible monic over Q; a root of unity has degree phi(order)."""
    d = h.degree()
    if primitive(h).LC() != 1:
        return None
    for n in range(1,2*d*d+1):
        if sp.totient(n)==d and poly(sp.cyclotomic_poly(n,T))==h:
            return n
    return None


def integral_growth_bound(h,q):
    """Exact Kronecker coefficient-growth witness and an exclusive j bound."""
    assert h.LC()==1 and all(c.q==1 for c in h.all_coeffs())
    assert torsion_order(h) is None and h.eval(0) != 0
    d = h.degree()
    matrix = multiplication_matrix(poly(T),h)
    power = sp.eye(d)
    m = 0
    witness = None
    while witness is None:
        m += 1
        power = power*matrix
        coefficients = sp.Poly(power.charpoly(T).as_expr(),T,domain=sp.ZZ)
        for k in range(1,d+1):
            value = abs(int(coefficients.nth(d-k)))
            if value > comb(d,k):
                witness = (k,value,comb(d,k))
                break
    k,value,binomial = witness
    denominator = 2
    while binomial*sp.Rational(denominator+1,denominator)**(m*k) >= value:
        denominator *= 2
    lower = sp.Rational(denominator+1,denominator)
    root_upper = 1+max(abs(h.nth(i)) for i in range(d))
    target_upper = max(sp.Integer(1),sum(abs(q.nth(i))*root_upper**i
                                       for i in range(q.degree()+1)))
    lo,hi = 0,1
    while lower**hi <= target_upper:
        hi *= 2
    while lo+1 < hi:
        mid = (lo+hi)//2
        if lower**mid > target_upper:
            hi = mid
        else:
            lo = mid
    assert lower**(hi-1) <= target_upper < lower**hi
    return hi,{
        'branch':'integral_nontorsion_growth','factor':encoded(h),
        'power_index':m,'coefficient_index':k,'coefficient_absolute_value':value,
        'unit_disc_binomial_bound':binomial,'conjugate_modulus_lower':str(lower),
        'all_conjugates_modulus_upper':str(root_upper),'target_modulus_upper':str(target_upper),
        'exclusive_exponent_bound':hi,
    }


def unit_power_solutions(f,q):
    """Solve T^j=q in Q[T]/f for j>=0, with f(0)!=0."""
    assert f.LC()==1 and f.eval(0) != 0
    q = q.rem(f)
    if sp.gcd(q,f).degree() > 0:
        return empty(),{'branch':'unit_power_cannot_equal_nonunit'}
    factors = [(h.monic(),int(e)) for h,e in sp.factor_list(f)[1]]
    repeated = next(((h,e) for h,e in factors if e >= 2),None)
    if repeated is not None:
        h,e = repeated
        inv = sp.invert(q,h)
        # Differentiate the full target BEFORE reducing modulo h.
        residue = (poly(T)*q.diff()*inv).rem(h)
        evidence = {'branch':'repeated_factor_derivative','factor':encoded(h),
                    'multiplicity':e,'forced_exponent_residue':encoded(residue)}
        if residue.degree() > 0:
            return empty(),evidence
        candidate = residue.nth(0)
        if candidate.q != 1 or candidate < 0:
            return empty(),evidence
        j = int(candidate)
        evidence['candidate'] = j
        good = ppow(poly(T),j,f)==q
        evidence['full_modulus_check'] = good
        return (singleton(j) if good else empty()),evidence
    orders = [(h,torsion_order(h)) for h,e in factors]
    exceptional = next(((h,n) for h,n in orders if n is None),None)
    if exceptional is None:
        d = lcm(*(n for h,n in orders))
        value = poly(1).rem(f)
        found = None
        for j in range(d):
            if value==q:
                assert found is None
                found = j
            value = (value*poly(T)).rem(f)
        assert value==poly(1).rem(f)
        evidence = {'branch':'squarefree_torsion','factor_orders':[n for h,n in orders],
                    'period':d,'tested_exponents':d}
        return (progression(found,d) if found is not None else empty()),evidence
    h,_ = exceptional
    qh = q.rem(h)
    assert not qh.is_zero
    if primitive(h).LC() > 1:
        denominator = lcm(*(int(c.q) for c in qh.all_coeffs()))
        degree = int(qh.degree())
        bound = degree+h.degree()*(denominator.bit_length()-1)
        exclusive = bound+1
        evidence = {'branch':'nonintegral_valuation_bound','factor':encoded(h),
                    'primitive_factor':encoded(primitive(h)),
                    'target_residue':encoded(qh),'target_degree':degree,
                    'coefficient_denominator':denominator,
                    'inclusive_exponent_bound':bound}
    else:
        exclusive,evidence = integral_growth_bound(h,qh)
    value = poly(1).rem(f)
    found = []
    for j in range(exclusive):
        if value==q:
            found.append(j)
        value = (value*poly(T)).rem(f)
    assert len(found) <= 1
    evidence['tested_exponents'] = exclusive
    return (singleton(found[0]) if found else empty()),evidence


def quotient_power_solutions(f,q):
    """Solve T^j=q modulo a nonconstant monic rational f, including T|f."""
    assert f.degree() >= 1 and f.LC()==1
    q = q.rem(f)
    zero_order = 0
    unit_part = f
    while unit_part.eval(0)==0:
        unit_part = unit_part.exquo(poly(T))
        zero_order += 1
    evidence = {'minimal_polynomial':encoded(f),'target_polynomial':encoded(q),
                'zero_primary_exponent':zero_order}
    value = poly(1).rem(f)
    for j in range(zero_order):
        if value==q:
            evidence['branch']='early_zero_primary_witness'
            return singleton(j),evidence
        value = (value*poly(T)).rem(f)
    if zero_order:
        if not q.rem(poly(T**zero_order)).is_zero:
            evidence['branch']='zero_primary_obstruction'
            return empty(),evidence
        if unit_part.degree()==0:
            evidence['branch']='nilpotent_zero_tail'
            return progression(zero_order,1),evidence
    solution,inner = unit_power_solutions(unit_part,q.rem(unit_part))
    evidence['unit_part'] = inner
    return intersect_tail(solution,zero_order),evidence


def orbit_solutions(matrix,start,target):
    """Exact complete exponent set for a rational matrix and two vectors."""
    matrix,start,target = sp.Matrix(matrix),sp.Matrix(start),sp.Matrix(target)
    n = matrix.rows
    assert n >= 1 and matrix.cols==n and start.shape==target.shape==(n,1)
    assert all(entry.is_Rational for obj in (matrix,start,target) for entry in obj)
    if start==sp.zeros(n,1):
        return (progression(0,1) if target==start else empty()),{'branch':'zero_start'}
    columns = [start]
    next_vector = matrix*start
    while True:
        basis = sp.Matrix.hstack(*columns)
        if basis.row_join(next_vector).rank()==len(columns):
            coefficients = basis.gauss_jordan_solve(next_vector)[0]
            d = len(columns)
            f = poly(T**d-sum(coefficients[i]*T**i for i in range(d)))
            break
        assert len(columns) < n
        columns.append(next_vector)
        next_vector = matrix*next_vector
    evidence = {'ambient_dimension':n,'cyclic_dimension':len(columns),
                'cyclic_minimal_polynomial':encoded(f)}
    if basis.row_join(target).rank()!=len(columns):
        evidence['branch']='target_outside_cyclic_span'
        return empty(),evidence
    coefficients = basis.gauss_jordan_solve(target)[0]
    q = poly(sum(coefficients[i]*T**i for i in range(len(columns))))
    solution,inner = quotient_power_solutions(f,q)
    evidence['quotient'] = inner
    return solution,evidence


def expect_rejected(name,action):
    try:
        action()
    except (AssertionError,ValueError):
        return name
    raise AssertionError('Negative control accepted: '+name)


def equal(a,b):
    assert a==b


def verify():
    records = []
    prefix_checks = 0
    cases = [
        ('jet_singleton',(T-1)**2,1+7*(T-1),singleton(7)),
        ('torsion_value_wrong_jet',(T+1)**2,-1,empty()),
        ('higher_jet_singleton',(T-1)**3,T**7,singleton(7)),
        ('higher_jet_corruption',(T-1)**3,T**7+(T-1)**2,empty()),
        ('quadratic_torsion_jet',(T*T+1)**2,T*T,singleton(2)),
        ('number_field_jet',(T*T-2)**2,T**3,singleton(3)),
        ('number_field_value_only_false',(T*T-2)**2,2,empty()),
        ('dyadic_denominator_hit',T-sp.Rational(1,2),sp.Rational(1,128),singleton(7)),
        ('dyadic_denominator_miss',T-sp.Rational(1,2),sp.Rational(1,129),empty()),
        ('unit_circle_nonintegral_hit',T*T-sp.Rational(6,5)*T+1,T**7,singleton(7)),
        ('unit_circle_nonintegral_miss',T*T-sp.Rational(6,5)*T+1,2,empty()),
        ('integral_unit_growth_hit',T*T-T-1,T**7,singleton(7)),
        ('integral_unit_growth_norm_obstruction',T*T-T-1,T+3,empty()),
        ('mixed_torsion_orders',(T*T+1)*(T*T+T+1),T**5,progression(5,12)),
        ('torsion_target_miss',T*T+1,T+1,empty()),
        ('nilpotent_early',T**3,T*T,singleton(2)),
        ('nilpotent_tail',T**3,0,progression(3,1)),
        ('nilpotent_identity',T**3,1,singleton(0)),
        ('nilpotent_nonzero_miss',T**3,2,empty()),
        ('mixed_zero_torsion_odd',T*T*(T+1),T**3,progression(3,2)),
        ('mixed_zero_torsion_even',T*T*(T+1),T*T,progression(2,2)),
        ('mixed_zero_early',T*T*(T+1),T,singleton(1)),
        ('mixed_zero_no_hit',T*T*(T+1),0,empty()),
        ('mixed_zero_growth',T*T*(T-2),T**5,singleton(5)),
        ('irrational_forced_jet',(T*T+1)**2,T+2,empty()),
        ('negative_forced_jet',(T-1)**2,2-T,empty()),
        ('fractional_forced_jet',(T-1)**2,1+(T-1)/2,empty()),
        ('nonunit_target_in_unit_algebra',(T-1)*(T-2),T-1,empty()),
        ('unit_witness_below_zero_threshold',T**3*(T-2),T**3/4,empty()),
    ]
    for name,ff,qq,expected in cases:
        f = poly(ff).monic();q=poly(qq).rem(f)
        matrix = multiplication_matrix(poly(T),f)
        start,target = vector(poly(1),f),vector(q,f)
        solution,evidence = orbit_solutions(matrix,start,target)
        assert solution==expected,(name,solution,expected)
        value = start
        for j in range(80):
            assert (value==target)==contains(solution,j),(name,j)
            value = matrix*value
            prefix_checks += 1
        records.append({'name':name,'solution':solution,'evidence':evidence})
    # An original nonreduced jet algebra need not give a repeated cyclic
    # minimal polynomial: a constant multiplier kills its nilpotent direction.
    f = poly((T-1)**2)
    matrix = multiplication_matrix(poly(1+(T-1)**2),f)
    solution,evidence = orbit_solutions(matrix,vector(poly(1),f),vector(poly(1),f))
    assert solution==progression(0,1) and evidence['cyclic_dimension']==1
    records.append({'name':'nonreduced_ambient_constant_multiplier','solution':solution,'evidence':evidence})
    for name,start,target,expected in [
        ('zero_start_equal',sp.zeros(2,1),sp.zeros(2,1),progression(0,1)),
        ('zero_start_unequal',sp.zeros(2,1),sp.Matrix([1,0]),empty()),
        ('outside_cyclic_span',sp.Matrix([1,0]),sp.Matrix([0,1]),empty()),
    ]:
        solution,evidence=orbit_solutions(sp.diag(2,3),start,target)
        assert solution==expected
        records.append({'name':name,'solution':solution,'evidence':evidence})
    # Rational restriction of scalars for Q(sqrt2)[epsilon]/epsilon^2.
    matrix=sp.Matrix([[0,2,0,0],[1,0,0,0],[1,0,0,2],[0,1,1,0]])
    start=sp.Matrix([1,0,0,0]);target=sp.Matrix([0,2,6,0])
    solution,evidence=orbit_solutions(matrix,start,target)
    assert solution==singleton(3)
    assert evidence['quotient']['unit_part']['branch']=='repeated_factor_derivative'
    records.append({'name':'number_field_full_jet_restriction','solution':solution,'evidence':evidence})
    change=sp.Matrix([[1,2,0,0],[0,1,1,0],[0,0,1,3],[0,0,0,1]])
    changed,_=orbit_solutions(change*matrix*change.inv(),change*start,change*target)
    assert changed==solution
    # Verify the finite cyclotomic index bound on a broad exact control range.
    for n in range(1,513):
        assert n <= 2*int(sp.totient(n))**2
    lookup = {record['name']:record for record in records}
    negatives=[
        expect_rejected('jet test falsely accepts torsion point value',lambda:equal(lookup['torsion_value_wrong_jet']['solution'],progression(1,2))),
        expect_rejected('higher jets discarded after forced exponent',lambda:equal(lookup['higher_jet_corruption']['solution'],singleton(7))),
        expect_rejected('unit-circle nonintegral root falsely periodic',lambda:equal(lookup['unit_circle_nonintegral_hit']['solution']['kind'],'progression')),
        expect_rejected('nilpotent tail starts too early',lambda:equal(lookup['nilpotent_tail']['solution'],progression(0,1))),
        expect_rejected('target outside cyclic span accepted',lambda:equal(lookup['outside_cyclic_span']['solution'],singleton(0))),
        expect_rejected('inexact input accepted',lambda:orbit_solutions([[sp.Float(1)]],[1],[1])),
    ]
    return {
        'algorithm_cases':len(records),'independent_prefix_comparisons':prefix_checks,
        'cyclotomic_index_controls':512,'rational_basis_change_controls':1,
        'negative_controls':negatives,'cases':records,
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root',type=Path,default=Path(__file__).resolve().parents[3])
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    repo=args.repo_root.resolve()
    directory=Path('research/foundations/extended-geometry')
    ledger_path=directory/'elementary_reciprocal_exponents_source_ledger.json'
    ledger=json.loads((repo/ledger_path).read_text())
    for source in ledger['sources']+ledger.get('context_sources',[]):
        actual=hashlib.sha256((repo/source['path']).read_bytes()).hexdigest()
        assert actual==source['sha256'],source['path']
    data=verify()
    result={
        'status':'PASS','python':platform.python_version(),'sympy':sp.__version__,
        'arithmetic':'Exact rational polynomial and matrix arithmetic; no floating point or Orbit oracle.',
        'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'source_pins_checked':len(ledger['sources'])+len(ledger.get('context_sources',[])),
        'packet_pins':{str(path):hashlib.sha256((repo/path).read_bytes()).hexdigest()
                       for path in (directory/'ELEMENTARY_RECIPROCAL_EXPONENT_DECISION.md',
                                    directory/Path(__file__).name,ledger_path)},
        'algorithm':'Implemented rational point-to-point exponent decision with complete solution descriptions.',
        'verification':data,
        'scope':'Finite exact controls for the implemented general algorithm. Written proof establishes universal termination; no polynomial-time claim or end-to-end fibre implementation.',
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','cases':data['algorithm_cases'],'prefix_comparisons':data['independent_prefix_comparisons'],'negative_controls':len(data['negative_controls']),'output':str(args.output)},sort_keys=True))


if __name__=='__main__':
    main()
