"""PF2RULED70: integral cone quotients, marked budget and unrestricted sector proofs."""
from pathlib import Path
import hashlib,importlib.util,itertools,json,math,os,platform
import sympy as S
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('base69',HERE/'verify_seifert_base_genus.py')
V=importlib.util.module_from_spec(spec);spec.loader.exec_module(V);E=V.E

def group(orders=(),free=0):
    orders=[int(x) for x in orders if x>1]
    if not orders:return dict(free=free,torsion=[])
    out=E.h1(S.diag(*orders).tolist());out['free']+=free;return out

def direct_sum(groups):
    return group([x for g in groups for x in g['torsion']],sum(g['free'] for g in groups))

def cone_group(orders,genus=0):
    if not orders:return group(free=2*genus)
    assert all(x>=1 for x in orders)
    mat=S.diag(*orders).col_join(S.ones(1,len(orders)))
    out=E.h1(mat.tolist());out['free']+=2*genus;return out

def dim(group,prime):return group['free']+sum(x%prime==0 for x in group['torsion'])

def lengths(pi):return sorted(len(x) for x in V.orbits([pi],range(len(pi))))

def piece_quotients(p,q,r,row):
    companion=[];cable=[]
    for data in V.companion_bases(row,p,q,r):
        a,b,M=data['induced'];la,lb=lengths(a),lengths(b)
        assert all(p%x==0 for x in la) and all(q%x==0 for x in lb)
        cones=[p//x for x in la]+[q//x for x in lb]
        companion.append(dict(cones=cones,genus=data['genus'],quotient=cone_group(cones,data['genus'])))
    for data in V.cable_bases(row):
        z,B,tau=data['induced'];lz=lengths(z);assert all(r%x==0 for x in lz)
        cones=[r//x for x in lz]
        cable.append(dict(cones=cones,genus=data['genus'],quotient=cone_group(cones,data['genus'])))
    return dict(companion=companion,cable=cable,
                companion_group=direct_sum([x['quotient'] for x in companion]),
                cable_group=direct_sum([x['quotient'] for x in cable]),
                protected=direct_sum([x['quotient'] for x in companion+cable]))

def table_groups(e,o,r,tag):
    assert e%2==0 and o%2==1 and o%3==0 and math.gcd(e,o)==1
    if tag==(24,24):return group([e//2,e,o//3]),group([r//3])
    if tag==(12,6):return group([e//2,o//3]),group([r//2])
    if tag==(6,6):return group([e//2]),group()
    assert tag==(24,2);return group([o//3]),group([r])

def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    preds=['marked_end_lattice_receipt.json','binomial_periodic_receipt.json','cyclic_companion_receipt.json',
           'cyclic_parity_receipt.json','seifert_base_genus_receipt.json']
    for name in preds:
        old=json.loads((HERE/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for p,h in old['source_sha256'].items():check('sealed predecessor '+name+' '+p,sha(HERE/p)==h)
    manifest=HERE/'sources/peripheral-cone-budget/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        p=HERE/row['path'];check('primary bytes '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    # The diagonal presentation is exact over Z, not just an order calculation.
    for orders in [[],[1],[9],[3,3],[4,6],[6,6,6],[5,5,5,5],[2,2,1,3,3],[7,14,14,1,3]]:
        out=cone_group(orders)
        for ell in [2,3,5,7]:
            check('prime dimension of entire cone quotient '+str((orders,ell)),dim(out,ell)==max(sum(a%ell==0 for a in orders)-1,0))
    a,beta,h=S.symbols('alpha beta h')
    check('killing regular fibre removes twisting coefficients',(a*S.Symbol('x')+beta*h).subs(h,0)==a*S.Symbol('x'))
    # Primitive canonical marking: kill eR in the rank2 source presentation.
    for mat in [[[-3,1],[1,-2]],[[-1,3],[3,-12]],[[-7,2],[2,-4]]]:
        out=E.h1(mat,[[1,0]])
        check('quartic marked quotient cyclic '+str(mat),out['free']<=1 and len(out['torsion'])<=1)
    check('nonprimitive vector does not improve budget at its bad prime',
          E.h1([[2,0],[0,2]],[[2,0]])==group([2,2]) and dim(group([2,2]),2)==2)
    check('individual divisor mark improves budget even at nonprimitive canonical prime',
          E.h1([[2,0],[0,2]],[[1,0]])==group([2]))
    check('two individual divisor marks remove both source generators',
          E.h1([[2,0],[0,2]],[[1,0],[0,1]])==group())
    # Complete finite permutation-type audit. Actual alpha_i remain unbounded.
    census=json.loads((HERE/'binomial_periodic_receipt.json').read_text())['census']['complete_cells']
    expected_types={(24,24):([1,1,2],[1,3],[1,3]),(12,6):([2,2],[1,3],[1,1,2]),
                    (6,6):([1,1,2],[1,3],[1,3]),(24,2):([4],[1,3],[1,1,1,1])}
    type_counts={tag:0 for tag in expected_types}
    for key,cell in census.items():
        p0,q0,r0=map(int,key.split(','));p=p0 if p0>=2 else p0+12;q=q0 if q0>=2 else q0+12;r=r0 if r0>=2 else r0+12
        for i,row in enumerate(cell['orbits']):
            tag=row['companion_order'],row['cable_order']
            if row['cycle_rank'] or tag not in expected_types:continue
            if sum(x['genus'] for x in V.companion_bases(row,p,q,r)):continue
            a,b,B,tau=map(tuple,row['rho']);even,odd=(a,b) if p%2==0 else (b,a);e,o=(p,q) if p%2==0 else (q,p)
            root=E.mul(E.inv(tau),B)
            check('entire planar noncyclic cycle type '+str((key,i)),
                  (lengths(even),lengths(odd),lengths(root))==expected_types[tag])
            check('noncyclic companion fibre identity '+str((key,i)),E.power(a,p)==tuple(range(4)))
            actual=piece_quotients(p,q,r,row);cg,pg=table_groups(e,o,r,tag)
            check('cone table from actual fibre-orbit cover '+str((key,i)),actual['companion_group']==cg and actual['cable_group']==pg)
            type_counts[tag]+=1
    check('all272 planar noncyclic periodic classes audited',type_counts=={(24,24):64,(12,6):96,(6,6):64,(24,2):48})
    # Exact algebra controls for arbitrary-parameter group identities proved in the note.
    for e,o in [(2,3),(4,3),(10,3),(14,15),(26,9),(8,15)]:
        check('unrestricted S4 companion cone form control '+str((e,o)),cone_group([e//2,e,e,o//3,o])==group([e//2,e,o//3]))
        check('unrestricted A4 companion cone form control '+str((e,o)),cone_group([e//2,e//2,o//3,o])==group([e//2,o//3]))
        if e>2:
            ell=min(S.factorint(e//2));check('every e greater2 forces noncyclic protected companion '+str((e,o)),dim(group([e//2,e,o//3]),int(ell))>=2)
    controls=[]
    cases=[(14,3,9,(24,24),[7,168],[7,42],True),
           (26,3,9,(24,24),[13,312],[13,78],True),
           (2,15,45,(24,24),[5,120],[5,30],True),
           (2,9,27,(24,24),[3,72],[3,18],True),
           (4,3,16,(12,6),[4,8],[4,8],True),
           (10,3,10,(12,6),[5,10],[5,10],True),
           (4,9,3,(24,2),[3,3],[3,3],True),
           (2,3,9,(24,24),[24],[6],False),
           (2,3,26,(12,6),[26],[26],False),
           (4,3,3,(24,2),[3],[3],False)]
    for p,q,r,tag,h1,hq,excluded in cases:
        key=','.join(str(x%12) for x in [p,q,r]);row=next(x for x in census[key]['orbits'] if
              (x['companion_order'],x['cable_order'])==tag and not x['cycle_rank'] and
              not sum(y['genus'] for y in V.companion_bases(x,p,q,r)))
        actual=piece_quotients(p,q,r,row);full=V.full_satellite(p,q,r,row);T=actual['protected']
        check('full ordinary satellite homology '+str((p,q,r,tag)),full==dict(H1=group(h1),quotient=group(hq)))
        check('marked protected quotient exclusion matches control '+str((p,q,r,tag)),(len(T['torsion'])>1)==excluded)
        for ell in [2,3,5,7,13]:
            check('full marked quotient retains protected prime dimension '+str((p,q,r,tag,ell)),dim(full['quotient'],ell)>=dim(T,ell))
        determinant=(r*(p*q-p-q)+1) if r%2==0 else r*q
        check('new controls pass earlier target determinant nine gate '+str((p,q,r,tag)),determinant%9==0)
        controls.append(dict(p=p,q=q,r=r,images=list(tag),protected=T,full=full,excluded=excluded))
    # Recover the known cyclic gates without importing an Alexander calculation.
    cyclic=[]
    for p,q,r,tag,want in [(3,5,32,(3,24),group([5,5])),(4,3,3,(4,24),group([3,3,3])),
                          (2,3,3,(4,24),group([3])),(14,3,3,(4,24),group([3]))]:
        key=','.join(str(x%12) for x in [p,q,r]);row=next(x for x in census[key]['orbits'] if
                    (x['companion_order'],x['cable_order'])==tag and not x['cycle_rank'])
        out=piece_quotients(p,q,r,row)
        check('cyclic predecessor gate recovered from cones '+str((p,q,r,tag)),out['protected']==want)
        cyclic.append(dict(p=p,q=q,r=r,images=list(tag),protected=out['protected']))
    passing=next(x for x in controls if (x['p'],x['q'],x['r'])==(2,3,9))
    check('protected quotient does NOT assert whole homology or zero R mark',passing['protected']==group([2,3]) and passing['full']['H1']!=passing['full']['quotient'])
    inputs=[Path(__file__),HERE/'PERIPHERAL_CONE_TORSION_AND_MARKED_BUDGET.md',manifest]+[HERE/n for n in preds]+[HERE/r['path'] for r in sources['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
        general=dict(all_finite_field_degrees=True,all_cabling_lengths=True,
          orientable_selected_Seifert_bases_required=True,arbitrary_boundary_gluings_and_fillings=True,
          protected_group_is_quotient_not_direct_summand=True,all_filling_core_marks_killed=True,
          canonical_mark_must_be_filling_core_combination=True,
          smooth_disjoint_entire_A1_source_hypotheses_required=True,
          individual_omitted_prime_marks_give_b_minus_s_bound=True,
          individual_mark_bound_needs_no_canonical_primitivity=True,
          all_actual_ramified_omitted_primes_are_filling_cores=True,
          unramified_omitted_primes_not_automatically_filling_cores=True,
          primitive_canonical_vector_required_for_canonical_only_all_prime_improvement=True,
          canonical_only_nonprimitive_bad_primes_keep_b_bound=True,
          cone_group='(direct sum Z/alpha_i)/<(1,...,1)>',
          prime_bound='b1(G)+sum[2g_v+max(N_v(ell)-1,0)]<=b-s; s tracked omitted filling-core knots',
          canonical_only_bound='b-1 if canonical vector primitive or nonzero modulo ell',
          no_cyclic_companion_assumption=True,actual_cone_orders_not_period12=True),
        quartic=dict(marked_end_quotient_cyclic=True,all_planar_noncyclic_cycle_classes_audited=272,
          S4_S4_all_even_first_exponents_greater2_excluded=True,
          S4_S4_remaining_even_exponent=2,S4_S4_remaining_gcd='gcd(o/3,r/3)=1',
          A4_S3_necessary_gcd='gcd(e*o/6,r/2)=1',S4_C2_necessary_gcd='gcd(o/3,r)=1',
          S3_S3_new_restriction=False,all_S4_S4_representations_excluded=False,
          no_whole_target_key_exclusion_from_one_representation=True,
          cyclic_predecessor_generator_gates_recovered_without_Alexander=True),
        periodic_cycle_type_counts={str(k):v for k,v in type_counts.items()},full_controls=controls,cyclic_controls=cyclic,
        assurance=dict(written_unbounded_quotient_prime_budget_and_cone_list_proofs=True,
          primary_presentation_visually_read=True,finite_controls_not_universal_proof=True,
          independent_spanning_tree_full_cover_predecessor_reused=True,
          protected_quotient_identified_with_H1=False,ramification_mark_declared_zero=False,
          entire_noncylclic_sector_closed=False,full_linking_matched=False,
          full_quartic_exclusion=False,higher_degree_exclusion=False,finite_Keller_map=False,
          PF63_polynomial_curves_declared_nonexistent=False,formal_verification=False,
          independent_full_reproof=False,external_review=False,novelty_priority=False,research_goal_completed=False),
        source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'peripheral_cone_budget_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; all-degree marked cone budget; S4/S4 forces even first exponent2')
if __name__=='__main__':main()
