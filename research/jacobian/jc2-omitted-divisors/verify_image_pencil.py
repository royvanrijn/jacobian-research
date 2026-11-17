"""Exact group and universal genus/omission controls for PF2OMIT18."""
from pathlib import Path
import hashlib
import json
import platform
import time
from itertools import combinations,permutations
import sympy as s
from verify_global_monodromy import closure,cycle,conjugates,sign,orbits

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CHECKS=[]


def check(label,condition):
    if not condition:
        raise AssertionError(label)
    CHECKS.append(label)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def partitions(total,least=1):
    if total==0:
        yield []
    for first in range(least,total+1):
        for rest in partitions(total-first,first):
            yield [first]+rest


def permutation_of_parts(parts):
    out=list(range(sum(parts)))
    offset=0
    for length in parts:
        letters=list(range(offset,offset+length))
        for a,b in zip(letters,letters[1:]+letters[:1]):
            out[a]=b
        offset+=length
    return tuple(out)


def main():
    started=time.monotonic()
    symmetric=set(permutations(range(6)))
    alternating={p for p in symmetric if sign(p)==1}
    stabilizer={p for p in symmetric if p[5]==5}
    even_stabilizer=stabilizer&alternating
    transposition=cycle(6,[0,1])
    three=cycle(6,[0,1,2])
    check('B transposition normal closure S6',closure(6,conjugates(transposition,symmetric))==symmetric)
    check('B three-cycle normal closure A6',closure(6,conjugates(three,symmetric))==alternating)
    check('geometric groups orders720,360',len(symmetric)==720 and len(alternating)==360)
    check('source geometric stabilizer orders120,60',len(stabilizer)==120 and len(even_stabilizer)==60)
    check('both generic curve covers degree6',len(symmetric)//len(stabilizer)==len(alternating)//len(even_stabilizer)==6)
    check('S5 acts nontrivially on sign constants',{sign(p) for p in stabilizer}=={-1,1})
    for name,group in [('S6',symmetric),('A6',alternating)]:
        check(name+' transitive',orbits(6,group)==[list(range(6))])
        pairs={(p[0],p[1]) for p in group}
        check(name+' two-transitive',len(pairs)==30)
        blocks=[]
        for length in range(1,7):
            for others in combinations(range(1,6),length-1):
                block={0,*others}
                if all(({p[i] for i in block}==block or not({p[i] for i in block}&block)) for p in group):
                    blocks.append(sorted(block))
        check(name+' no proper curve-cover block',blocks==[[0],list(range(6))])
    b,u,n,e=s.symbols('b u n e',integer=True)
    genus=(140+(e-1)*b-n)/2
    punctures=n+b+u
    chi=-138-e*b-u
    check('universal nonlinear-pencil RH identity',s.expand(2*genus-2-6*(2*12-2)-(e-1)*b-(6-n))==0)
    check('universal source genus/puncture/Euler identity',s.expand(2-2*genus-punctures-chi)==0)
    check('universal finite-fibre Euler integration',s.expand(6*(-23)-(e*b+u)-chi)==0)
    check('universal different vs removed-point identity',s.expand((e*b+u)-(b+u)-(e-1)*b)==0)
    check('universal affine-normal Euler identity',s.expand(2-2*genus-n+138+(e-1)*b)==0)
    m=s.symbols('m',integer=True,positive=True)
    for parity,bvalue,nvalues in [(0,2*m,[2,4,6]),(1,2*m+1,[1,3,5])]:
        options={s.expand(genus.subs({e:2,b:bvalue,n:nvalue})-s.ceiling(bvalue/2)) for nvalue in nvalues}
        check(f'e2 universal parity{parity} genus set',options=={s.Integer(67),s.Integer(68),s.Integer(69)})
    check('e3 universal genus set',{s.expand(genus.subs({e:3,n:nvalue})-b) for nvalue in [2,4,6]}=={s.Integer(67),s.Integer(68),s.Integer(69)})
    infinity=[]
    for parts in partitions(6):
        perm=permutation_of_parts(parts)
        check('infinity partition'+str(parts)+' parity',sign(perm)==(-1)**(6-len(parts)))
        infinity.append(dict(cycles=parts,geometric_point_count=len(parts),different=6-len(parts),
            sign=sign(perm),allowed_in_A6=perm in alternating,geometric_realization_claimed=False))
    check('all11 geometric infinity partitions',len(infinity)==11)
    check('exactly6 even infinity partitions',sum(r['allowed_in_A6'] for r in infinity)==6)
    rows=[]
    for row,other_index,u_formula in [('A_distinct',2,'deg h|H>=6'),('A_H=C',2,'0'),('A_H=B',2,'b'),
                                      ('B',2,'0'),('N',2,'0'),('C',3,'0')]:
        rows.append(dict(row=row,e_S=other_index,B_intersection_degree='b>=6',
            unramified_omitted_count_on_generic_level=u_formula,
            geometric_generic_monodromy='S6' if other_index==2 else 'A6',
            Galois_closure_constant_field='C(h)' if other_index==2 else 'C(h)(sqrt(h))',
            Galois_closure_geometric_component_count=1 if other_index==2 else 2,
            source_constant_field='C(h)',source_geometric_component_count=1,
            geometric_curve_stabilizer='S5' if other_index==2 else 'A5',
            proper_intermediate_generic_curve_fields=[],
            compact_source_genus_formula='2g=140+(e_S-1)*b-n_h',
            source_punctures_formula='n_h+b+u',source_Euler_formula='-138-e_S*b-u',
            actual_b=None,actual_u=None,actual_infinity_partition=None,actual_generic_curve_cover=None))
    for other_index,extra,expected_genus,expected_punctures,expected_Euler in [(2,0,75,14,-162),(2,12,75,26,-174),(3,0,81,14,-174)]:
        values={e:other_index,b:12,n:2,u:extra}
        check(f'numeric control e{other_index} u{extra}',genus.subs(values)==expected_genus and punctures.subs(values)==expected_punctures and chi.subs(values)==expected_Euler)
    check('unramified omission changes punctures but not compact RH',s.diff(genus,u)==0 and s.diff(punctures,u)==1 and s.diff(chi,u)==-1)
    check('minimum genus e2/e3',s.ceiling(s.Integer(6)/2)+67==70 and 6+67==73)
    check('quintic minimum genus e2/e3',s.ceiling(s.Integer(10)/2)+67==72 and 10+67==77)
    prior=json.loads((HERE/'critical_trace_ledger.json').read_text())
    check('universal intersection and unique rational fibre input',prior['first_image_is_unique_rational_fibre_of_its_equation'] and prior['all_residual_nonproperness_supports_connected'])
    ledger=dict(as_of='2026-11-14',claim='PF2OMIT18',
        hypotheses='Four oriented PF2OMIT13 degree6/index2/gap21 rows, actual full ramification and source field; no degree/contact restriction on B/H',
        target_pencil='h=H_C',generic_compact_target_genus=12,generic_affine_target_Euler=-23,
        all_source_generic_fibres_geometrically_integral=True,rows=rows,
        infinity_geometric_cycle_options=infinity,
        infinity_count_is_general_target_line_count=False,
        normal_closure_constant_quadratic_field_is_inside_source=False,
        minimum_compact_source_genus_e2=70,minimum_compact_source_genus_e3=73,
        numeric_controls_are_global_Keller_realizations=False,
        central_h_zero_fibre_accounted_by_generic_formulas=False,
        actual_B_equation=None,actual_H_equation=None,actual_global_cover=None,actual_source_realization=None,
        unrestricted_singleton_zero_properness_proved=False,full_goal_complete=False)
    ledger_path=HERE/'image_pencil_ledger.json'
    ledger_path.write_text(json.dumps(ledger,indent=2)+'\n')
    inputs=['jacobian/jc2-omitted-divisors/'+name for name in ['IMAGE_PENCIL_MONODROMY.md',
        'GLOBAL_MONODROMY.md','verify_global_monodromy.py','NONIMMERSED_LEDGER.md','UNIVERSAL_IMAGE_INTERSECTION.md',
        'critical_trace_ledger.json','DISJOINT_IMAGE_GATE.md','GENERIC_LINE_ACCOUNTING.md','SIMPLE_RAMIFICATION_AND_RADIAL.md']]
    receipt=dict(as_of='2026-11-14',outcome='PASS',check_count=len(CHECKS),checks=CHECKS,
        checker_sha256=digest(Path(__file__)),input_sha256={p:digest(ROOT/p) for p in inputs},ledger_sha256=digest(ledger_path),
        runtime_seconds=round(time.monotonic()-started,6),python=platform.python_version(),sympy=s.__version__,
        proof_boundary='Inertia fixes the algebraic pencil constants, regular-field base change and sign-field intersection are written field arguments. Exact S6/A6 closures, stabilizers, blocks, all infinity partitions, symbolic RH/Euler identities and unramified-deletion controls are replayed. No actual curve cover, infinity realization, source identification or central-fibre calculation is claimed.',
        formal_verification=False,external_review=False,full_goal_complete=False)
    (HERE/'image_pencil_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['outcome','check_count','runtime_seconds']}))


if __name__=='__main__':
    main()
