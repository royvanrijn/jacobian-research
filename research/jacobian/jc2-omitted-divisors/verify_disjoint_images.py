"""Exact pole-module, Milnor-catalogue and critical-trace controls."""
from pathlib import Path
import hashlib
import json
import platform
import sys
import time
import sympy as s

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE.parent/'jc2-r7-image-family'))
import family as F

CHECKS=[]


def check(label,condition):
    if not condition:
        raise AssertionError(label)
    CHECKS.append(label)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def monomial(exponents,U,V):
    return U**exponents[0]*V**exponents[1]


def critical_witness(chart):
    U,V=s.symbols('U V')
    seed=F.line if chart=='L' else F.laurent
    parameter=F.d if chart=='L' else F.x
    p=s.expand(F.p.subs(seed).subs(parameter,1))
    q=s.expand(F.q.subs(seed).subs(parameter,1))
    h=s.resultant(p-U,V-q,F.t)
    h=s.expand(h/s.Poly(h,V).LC())
    gradients=[s.diff(h,U),s.diff(h,V)]
    gb=s.groebner(gradients,U,V,order='grevlex',domain=s.QQ)
    check(chart+' gradient generators reduce to zero',all(gb.reduce(f)[1]==0 for f in gradients))
    leading=[f.LM(order=gb.order).exponents for f in gb.polys]
    for i,f in enumerate(gb.polys):
        for j in range(i):
            g=gb.polys[j]
            common=tuple(max(a,b) for a,b in zip(leading[i],leading[j]))
            spoly=(monomial(tuple(c-a for c,a in zip(common,leading[i])),U,V)*f.as_expr()/f.LC(order=gb.order)
                  -monomial(tuple(c-a for c,a in zip(common,leading[j])),U,V)*g.as_expr()/g.LC(order=gb.order))
            check(f'{chart} Buchberger pair{i},{j}',gb.reduce(s.expand(spoly))[1]==0)
    bound_U=min(a for a,b in leading if b==0)
    bound_V=min(b for a,b in leading if a==0)
    basis=[(i,j) for i in range(bound_U) for j in range(bound_V)
           if not any(i>=a and j>=b for a,b in leading)]
    check(chart+' exact standard basis dimension24',len(basis)==24)
    positions={m:j for j,m in enumerate(basis)}
    matrix=s.zeros(24)
    for j,m in enumerate(basis):
        rem=s.Poly(gb.reduce(s.expand(h*monomial(m,U,V)))[1],U,V)
        check(f'{chart} multiplication column{j} in standard basis',all(mon in positions for mon,c in rem.terms()))
        for mon,coefficient in rem.terms():
            matrix[positions[mon],j]=coefficient
    tr1=s.trace(matrix)
    tr2=sum(matrix[i,j]*matrix[j,i] for i in range(24) for j in range(24))
    moment=12*tr2-tr1**2
    check(chart+' nonzero disjoint-fibre obstruction',moment!=0)
    check(chart+' norm retains image degrees',s.total_degree(h)==10 and s.degree(h,V)==6)
    return dict(chart=chart,parameter=str(parameter),parameter_value=1,p=str(p),q=str(q),
        h=str(h),groebner_order='grevlex',groebner_basis=[str(f.as_expr()) for f in gb.polys],
        ideal_equality_boundary='Gradient ideal contained in checked Groebner ideal; written global Milnor dimension24 and checked quotient dimension24 imply equality.',
        standard_monomials=[list(m) for m in basis],
        multiplication_matrix_sparse=[[i,j,str(matrix[i,j])] for i in range(24) for j in range(24) if matrix[i,j]],
        trace_h=str(tr1),trace_h_squared_operator=str(tr2),moment_12tr2_minus_tr1_squared=str(moment),
        disjoint_polynomial_image_excluded=True,global_Keller_cover_claimed=False)


def main():
    started=time.monotonic()
    manifest_path=HERE/'sources/disjoint_image_manifest.json'
    manifest=json.loads(manifest_path.read_text())
    for source in manifest['sources']:
        path=manifest_path.parent/source['path']
        check('source pin '+source['path'],path.stat().st_size==source['bytes'] and digest(path)==source['sha256'])
    degrees=[0,10,20,9,19,29]
    apery=[min(d for d in degrees if d%6==r) for r in range(6)]
    check('module basis distinct residues',len({d%6 for d in degrees})==6)
    check('Apéry representatives',[0,19,20,9,10,29]==apery)
    check('module codimension12',(sum(apery)-sum(range(6)))//6==12)
    gaps=[n for n in range(max(apery)+1) if n<apery[n%6]]
    check('complete semigroup gap list',gaps==[1,2,3,4,5,7,8,11,13,14,17,23])
    check('conductor24',max(gaps)+1==24 and all(n>=apery[n%6] for n in range(24,30)))
    generated={6*i+9*j+10*k for i in range(40) for j in range(40) for k in range(40)}
    check('module degrees agree with generated semigroup through conductor and one residue block',
          all((n>=apery[n%6])==(n in generated) for n in range(36)))
    geometry=json.loads((HERE.parent/'jc2-r7-complement/geometry.json').read_text())
    G=s.sympify(geometry['closest_quintics'][0]['equation'],locals={'U':s.Symbol('U'),'V':s.Symbol('V'),'d':F.d}).subs(F.d,0)
    p=s.expand(F.p.subs(F.line).subs(F.d,0))
    q=s.expand(F.q.subs(F.line).subs(F.d,0))
    g=s.Poly(s.expand(G.subs({s.Symbol('U'):p,s.Symbol('V'):q})),F.t)
    check('exact degree9 ring generator control',g.degree()==9 and g.LC()==-s.Rational(1,6561))
    check('root derivative pole charge29',4*10-11==29)
    check('shift quadratic correction is lower order',40-2*29==-18)
    check('generic genus and Milnor budget',(10-1)*(10-2)//2-24==12 and 1-(2-2*12-1)==24)
    catalogue=json.loads((HERE/'generic_line_ledger.json').read_text())
    milnor_weights=dict(nodes=1,triples=4,tacnodes=3,cusps=2,A5=5)
    catalogue_checks=[]
    for row in catalogue['first_image_catalogue']:
        mu=sum(milnor_weights[k]*v for k,v in row['singularity_counts'].items())
        check(row['chart']+row['stratum']+' nodal/disjoint budget',mu==12 if row['open_stratum'] else mu>12)
        catalogue_checks.append(dict(chart=row['chart'],stratum=row['stratum'],total_affine_Milnor=mu,
            disjoint_polynomial_image_excluded_by_budget=not row['open_stratum']))
    check('nonnodal union268',catalogue['first_image_union']['nonnodal_coefficient_points']==268)
    witnesses=[critical_witness(chart) for chart in ['L','K']]
    # Abstract trace control: necessary, never a converse over C.
    lam=s.symbols('lambda',nonzero=True)
    check('two nodal fibre trace identity',s.expand(12*(12*lam**2)-(12*lam)**2)==0)
    check('three disjoint polynomial fibres exceed Milnor budget',3*12>24)
    check('a cusped polynomial fibre cannot be in a disjoint pair',13+12>24)
    check('complex zero moment need not equal values',sum([1,-1,s.I,-s.I])==0 and sum(z*z for z in [1,-1,s.I,-s.I])==0)
    ledger=dict(as_of='2026-11-14',claim='PF2OMIT16',
        disjoint_polynomial_curve_rule='Equal degrees; H_B=alpha*H_C+beta with alpha,beta nonzero',
        only_possible_disjoint_degree_B=10,only_possible_disjoint_k=2,
        affine_ring_module_basis=['1','q','q^2','g','q*g','q^2*g'],module_basis_degrees=degrees,
        affine_pole_semigroup=[6,9,10],affine_pole_gaps=gaps,affine_pole_conductor=24,
        local_infinity_branch_semigroup=[4,10,41],affine_delta=12,
        every_shifted_fibre_irreducible_one_infinity_place=True,shifted_fibre_infinity_delta=24,
        generic_shifted_fibre_genus=12,total_affine_critical_Milnor_budget=24,
        disjoint_pair_forces='Two twelve-node fibres, no other critical values',
        support_connected_unconditionally_in_groupings=['A_distinct','B','N'],
        only_potentially_disconnected_groupings=['A_H=C','A_H=B','C'],
        disconnected_support_forces='Exactly two distinct degree10 twelve-node fibres of H_C, trace moment0',
        nonnodal_first_image_points_excluding_disjointness=268,catalogue_Milnor_checks=catalogue_checks,
        critical_trace_witnesses=witnesses,
        generic_nodal_open_in_each_family_excluding_disjoint_images=True,
        remaining_disjoint_candidate_parameters='Finite uncomputed nodal parameter set satisfying trace moment=0 (and any unresolved specialization poles); k2 and fibre-translation form required',
        trace_moment_zero_is_sufficient=False,
        current_quintic_intersection_degree_options=[10,12,16,18,20,22,24,26,28,30],
        singleton_forced_groupings=['A_H=C','A_H=B','C'],
        singleton_exists_if='k!=2, or C nonnodal, or C in trace-certified nodal open',
        shared_A_omitted_lengths=[2,2,1],cyclic_C_omitted_lengths=[2,3],retained_length_at_each_C_B_intersection=1,
        connected_singular_nonproperness_set_excluded=False,
        smooth_connected_Euler_theorem_applied_to_singular_union=False,
        actual_B_equation=None,actual_H_equation=None,actual_intersection_values=None,
        explicit_finite_trace_zero_set=None,actual_source_realization=None,
        unrestricted_singleton_zero_properness_proved=False,full_goal_complete=False)
    ledger_path=HERE/'disjoint_image_ledger.json'
    ledger_path.write_text(json.dumps(ledger,indent=2)+'\n')
    inputs=['jacobian/jc2-omitted-divisors/'+n for n in [
        'DISJOINT_IMAGE_GATE.md','GENERIC_LINE_ACCOUNTING.md','generic_line_ledger.json','ODD_IMAGE_INTERSECTIONS.md',
        'RESIDUE_LEDGER.md','NONIMMERSED_LEDGER.md','sources/disjoint_image_manifest.json']]
    inputs+=['jacobian/'+n for n in [
        'jc2-r7-complement/PROOF.md','jc2-r7-complement/geometry.json','jc2-r7-image-family/family.py',
        'jc2-r7-quintic-normalizations/PROOF.md','jc2-r7-quintic-incidence/PROOF.md']]
    inputs+=['foundations/scripts/verify_f2_degree_6_10_gap.py']
    inputs+=['jacobian/jc2-omitted-divisors/sources/'+row['path'] for row in manifest['sources']]
    receipt=dict(as_of='2026-11-14',outcome='PASS',check_count=len(CHECKS),checks=CHECKS,
        checker_sha256=digest(Path(__file__)),input_sha256={p:digest(ROOT/p) for p in inputs},ledger_sha256=digest(ledger_path),
        runtime_seconds=round(time.monotonic()-started,6),python=platform.python_version(),sympy=s.__version__,
        proof_boundary='Disjoint-curve lemma, module codimension, shifted infinity branch, global Milnor dimension24, ideal dimension comparison and openness are written arguments. Exact pins, semigroup arithmetic, catalogue Milnor charges, Buchberger pairs, gradient containment and two full rational trace matrices are replayed. No generic solver printout, incomplete symbolic trace or converse to the complex moment is assumed.',
        formal_verification=False,external_review=False,full_goal_complete=False)
    (HERE/'disjoint_image_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['outcome','check_count','runtime_seconds']}))


if __name__=='__main__':
    main()
