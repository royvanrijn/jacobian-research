"""Replay universal disjoint-image exclusion without a Groebner solver."""
from pathlib import Path
import gzip
import hashlib
import json
import platform
import sys
import time
import sympy as s
import flint
from flint import fmpq, fmpq_poly

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE.parent/'jc2-r7-image-family'))
import family as F
import critical_trace as C

CHECKS=[]


def check(label,condition):
    if not condition:
        raise AssertionError(label)
    CHECKS.append(label)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def modular_identity(polys,record,mutate=False):
    prime=record['prime']
    assert prime>1 and all(prime%i for i in range(2,int(prime**0.5)+1))
    first=list(record['first'])
    if mutate:
        first[0]=(first[0]+1)%prime
    result=[0]*(max(len(polys[0])+len(first),len(polys[1])+len(record['second']))-1)
    for poly,multiplier in zip(polys,[first,record['second']]):
        for i,a in enumerate(poly):
            for j,b in enumerate(multiplier):
                result[i+j]=(result[i+j]+a*b)%prime
    return result[0]==1 and all(co==0 for co in result[1:])


def source_h(chart):
    U,V,z=s.symbols('U V z')
    if chart=='L':
        p=s.expand(F.p.subs(F.line));q=s.expand(F.q.subs(F.line));parameter=F.d
    else:
        def scaled(poly,degree):
            raw=s.Poly(s.expand(poly.subs(F.laurent).subs(F.t,F.t/F.x)*F.x**degree),F.x,F.t)
            assert all(m[0]%3==0 for m,co in raw.terms())
            return s.expand(sum(co*z**(m[0]//3)*F.t**m[1] for m,co in raw.terms()))
        p,q=scaled(F.p,6),scaled(F.q,10);parameter=z
    h=s.resultant(p-U,V-q,F.t)
    check(chart+' reconstructed norm monic V6',s.Poly(h,V).LC()==1 and s.degree(h,V)==6)
    return C.decode(C.encode_expression(h,U,V,parameter))


def replay(chart):
    start=time.monotonic()
    path=HERE/f'critical_trace_{chart}.json.gz'
    data=json.loads(gzip.decompress(path.read_bytes()))
    check(chart+' certificate chart',data['chart']==chart)
    H=source_h(chart)
    gradients=[(H[0].derivative(i),H[1]) for i in range(2)]
    generators=[C.decode(row) for row in data['generators']]
    lifts={(row['i'],row['j']):C.decode(row['value']) for row in data['lifts']}
    check(chart+' six certificate generators and distinct lift entries',len(generators)==6 and len(lifts)==len(data['lifts']))
    for j,G in enumerate(generators):
        terms=[(1,[G])]+[(-1,[gradients[i],lifts[i,j]]) for i in range(2) if (i,j) in lifts]
        check(f'{chart} literal gradient membership{j}',C.sum_products_zero(terms))
    leading=[C.leading_UV(G) for G in generators]
    bU=min(i for i,j in leading if j==0)
    bV=min(j for i,j in leading if i==0)
    standard={(i,j) for i in range(bU) for j in range(bV) if not any(i>=a and j>=b for a,b in leading)}
    basis=[tuple(m) for m in data['basis']]
    check(chart+' full24-element leading staircase',len(basis)==24 and len(set(basis))==24 and set(basis)==standard)
    entries={(row['i'],row['j']):C.decode(row['value']) for row in data['entries']}
    quotients={(row['i'],row['j']):C.decode(row['value']) for row in data['column_quotients']}
    check(chart+' distinct matrix entries in range',len(entries)==len(data['entries']) and all(0<=i<24 and 0<=j<24 for i,j in entries))
    check(chart+' distinct quotient entries in range',len(quotients)==len(data['column_quotients']) and all(0<=i<6 and 0<=j<24 for i,j in quotients))

    def column_identity(j,matrix):
        terms=[(1,[H,C.monomial(*basis[j])])]
        terms+=[(-1,[matrix[i,j],C.monomial(*basis[i])]) for i in range(24) if (i,j) in matrix]
        terms+=[(-1,[quotients[k,j],G]) for k,G in enumerate(generators) if (k,j) in quotients]
        return C.sum_products_zero(terms)

    for j in range(24):
        check(f'{chart} literal multiplication column{j}',column_identity(j,entries))
    if chart=='L':
        changed=dict(entries)
        i,j=next(iter(changed))
        num,den=changed[i,j]
        changed[i,j]=(num+C.CTX.constant(1),den)
        check('matrix-entry mutation rejected',not column_identity(j,changed))
    traces={key:C.decode(row) for key,row in data['traces'].items()}
    actual=C.trace_powers(entries,24)
    for j,(num,den) in enumerate(actual,1):
        check(f'{chart} direct trace power{j}',C.same_univariate((num,den),C.univariate(traces[f'TR{j}'])))
        tnum,tden=C.univariate(traces[f'TR{j}'])
        check(f'{chart} trace{j} polynomial degree{30*j}',tnum.degree()==30*j and tden.degree()==0)
    check(chart+' moment2 identity',C.sum_products_zero([
        (1,[traces['M2']]),(-12,[traces['TR2']]),(1,[traces['TR1'],traces['TR1']])]))
    check(chart+' moment3 identity',C.sum_products_zero([
        (1,[traces['M3']]),(-144,[traces['TR3']]),(1,[traces['TR1'],traces['TR1'],traces['TR1']])]))
    polys=[[int(co) for co in row] for row in data['moment_primitive_polynomials']]
    degrees=[]
    for key,coefficients in zip(['M2','M3'],polys):
        numerator,denominator=C.univariate(traces[key])
        primitive=fmpq_poly(coefficients)
        quotient,remainder=divmod(numerator,primitive)
        check(chart+key+' primitive moment scalar and constant denominator',
              not remainder and quotient.degree()==0 and bool(quotient) and denominator.degree()==0)
        degrees.append(primitive.degree())
        check(chart+key+' degree preserved modulo prime',coefficients[-1]%data['modular_bezout']['prime']!=0)
    check(chart+' moment degrees60,90',degrees==[60,90])
    check(chart+' literal modular Bezout identity',modular_identity(polys,data['modular_bezout']))
    check(chart+' Bezout-coefficient mutation rejected',not modular_identity(polys,data['modular_bezout'],mutate=True))
    print(json.dumps(dict(chart=chart,outcome='PASS',seconds=round(time.monotonic()-start,6))),flush=True)
    return dict(chart=chart,parameter=data['parameter'],certificate_sha256=digest(path),
        certificate_bytes=path.stat().st_size,gradient_memberships=len(generators),multiplication_columns=24,
        dimension=24,moment_degrees=degrees,primitive_moment_polynomials=polys,
        modular_bezout=data['modular_bezout'],all_allowed_parameters_exclude_disjoint_images=True,
        basis_denominator_zeros_removed=False)


def main():
    started=time.monotonic()
    rows=[replay(chart) for chart in ['L','K']]
    ledger=dict(as_of='2026-11-14',claim='PF2OMIT17',charts=rows,
        theorem='Every additional irreducible polynomially normalized plane image meets every complete gap21 first image C',
        global_trace_extension='Constant total critical Milnor length24 and local intersection conservation give finite flat critical algebras near every allowed parameter; generic polynomial traces specialize everywhere',
        all_nonnodal_and_nodal_parameters_covered=True,scaled_K_z_zero_included=False,
        actual_additional_image_degrees_restricted=False,all_residual_nonproperness_supports_connected=True,
        first_image_is_unique_rational_fibre_of_its_equation=True,nonzero_level_normalization_genus_at_least=1,
        singleton_forced_groupings=['A_H=C','A_H=B','C'],shared_A_omitted_lengths=[2,2,1],
        cyclic_C_omitted_lengths=[2,3],retained_length_at_every_C_B_intersection=1,
        smooth_connected_Euler_theorem_applied_to_singular_union=False,
        actual_B_equation=None,actual_H_equation=None,actual_exceptional_values=None,
        actual_global_normalization=None,actual_source_realization=None,
        unrestricted_singleton_zero_properness_proved=False,full_goal_complete=False)
    ledger_path=HERE/'critical_trace_ledger.json'
    ledger_path.write_text(json.dumps(ledger,indent=2)+'\n')
    inputs=['jacobian-followup/jc2-omitted-divisors/'+name for name in [
        'UNIVERSAL_IMAGE_INTERSECTION.md','critical_trace.py','critical_trace_L.json.gz','critical_trace_K.json.gz',
        'DISJOINT_IMAGE_GATE.md','RESIDUE_LEDGER.md','NONIMMERSED_LEDGER.md','SINGLETON_INVERSE_VALUATIONS.md']]
    inputs+=['jacobian-followup/jc2-r7-image-family/family.py',
        'archive/non-elliptic/scripts/verify_f2_degree_6_10_gap.py']
    receipt=dict(as_of='2026-11-14',outcome='PASS',check_count=len(CHECKS),checks=CHECKS,
        checker_sha256=digest(Path(__file__)),input_sha256={p:digest(ROOT/p) for p in inputs},ledger_sha256=digest(ledger_path),
        runtime_seconds=round(time.monotonic()-started,6),python=platform.python_version(),sympy=s.__version__,flint=flint.__version__,
        proof_boundary='Solver-free literal gradient memberships,24 leading monomials, all48 multiplication-column identities, exact trace powers1..3, moment identities and good-prime univariate Bezout certificates replayed. Rank24 Milnor geometry and finite-flat trace extension across all parameter denominators remain written arguments. No global cover, source identification, formalization or external review.',
        formal_verification=False,external_review=False,full_goal_complete=False)
    (HERE/'critical_trace_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['outcome','check_count','runtime_seconds']}))


if __name__=='__main__':
    main()
