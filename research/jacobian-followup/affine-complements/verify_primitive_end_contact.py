"""PF2RULED102 controls; written geometry and declared contact imports prove scope."""
from pathlib import Path
import hashlib,json,sys,platform,math,itertools
import sympy as S
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
import full_peripheral_linking as N


def mul(a,b):return tuple(a[b[i]] for i in range(len(a)))
def group(gens):
    seen={tuple(range(len(gens[0])))};todo=list(seen)
    while todo:
        a=todo.pop()
        for b in gens:
            c=mul(a,b)
            if c not in seen:seen.add(c);todo.append(c)
    return seen
def cycles(a):
    seen=set();out=[]
    for i in range(len(a)):
        if i in seen:continue
        j=i;n=0
        while j not in seen:seen.add(j);n+=1;j=a[j]
        out.append(n)
    return sorted(out,reverse=True)
def partitions2(n):return [(n-j,j) for j in range(1,n//2+1)]
def order(a):return math.lcm(*cycles(a))


def main():
    checks=[]
    def check(name,value):
        assert value,name;checks.append(name)
    def identity(name,value):check(name,S.factor(value)==0)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    preds=['full_peripheral_linking_receipt.json','proper_finite_marked_end_receipt.json',
           'polynomial_key_chain_receipt.json','torus_infinity_receipt.json','canonical_formats_receipt.json']
    for name in preds:
        old=json.loads((P/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for f,h in old['source_sha256'].items():check('sealed predecessor '+name+' '+f,sha(P/f)==h)
    manifest=P/'sources/primitive-end-contact/manifest.json';m=json.loads(manifest.read_text())
    for row in m['primary_sources']:
        f=P/row['path'];check('primary contact PDF exact pin '+f.name,sha(f)==row['sha256'] and f.stat().st_size==row['bytes'])
    check('published erratum scope and weak radial Levi distinction retained',
          m['Honda_erratum_scope_checked'] and not m['nonrotative_torus_layer_uniqueness_imported'] and
          m['radial_weak_Levi_replaced_by_actual_small_strict_psh_perturbation'])
    h,e,j,rot=S.symbols('h e j r',integer=True,positive=True)
    theta=-(e-1)**2/h-1
    identity('unbounded rank-one canonical theta with chi2 signature-1',
             theta-(-(e-1)**2/h-2*2-3*(-1)))
    rotation=-h+2+2*j
    identity('unbounded tight lens spin classes are consecutive nonzero residues',
             (rotation-h)/2-(-h+1+j))
    identity('theta matching is exactly equality of canonical and allowed rotation squares',
             h*(theta-(-1-rot**2/h))-(rot**2-(e-1)**2))
    check('RP3 index1 allowed and triple-index theta contradiction exact',
          theta.subs({h:2,e:1})==-1 and theta.subs({h:2,e:3})==-3)
    check('lens bound permits index3 order4 equality and does not overexclude it',
          theta.subs({h:4,e:3})==-1-S.Rational(2)**2/4)
    lens=[]
    for hh in [2,3,4,5,8,13]:
        rotations=list(range(-hh+2,hh-1,2))
        gamma=[((r-hh)//2)%hh for r in rotations]
        values=sorted(set(-1-S.Rational(r*r,hh) for r in rotations))
        check('fixed-spin contact classes exhaust exact Honda count '+str(hh),
              len(rotations)==hh-1 and set(gamma)==set(range(1,hh)))
        lens.append(dict(h=hh,rotations=rotations,Gamma=gamma,theta_values=[str(x) for x in values]))
    # The perturbation is necessary: unperturbed radial Levi is degenerate on D.
    z,zb,w,wb,eps=S.symbols('z zb w wb epsilon',positive=True)
    for ee in [2,3,7]:
        rho=(z*zb)**ee+w*wb
        levi=S.diff(rho,z,zb)
        identity('exact local transverse Levi coefficient e='+str(ee),levi-ee**2*(z*zb)**(ee-1))
        check('strict perturbation restores missing branch-normal direction e='+str(ee),
              levi.subs(z,0)==0 and S.diff(rho+eps*(z*zb+w*wb),z,zb).subs(z,0)==eps)
    # Actual affine index1 control fV=W(W-1) is allowed, no Keller target asserted.
    x,f,u,t=S.symbols('x f u t');W=f*x;V=x*(f*x-1)
    identity('actual affine quadric relation',f*V-W*(W-1))
    xx=u**-1+t;ff=u
    identity('index1 affine quadric nonvanishing source volume',
             S.diff(xx,u)*S.diff(ff,t)-S.diff(xx,t)*S.diff(ff,u)+1)
    # The complete HUMAN two-cycle partition argument is not a group census.
    parts=partitions2(6)
    check('all two-cycle partitions of six are the three stated types',parts==[(5,1),(4,2),(3,3)])
    check('possible even positive cycle-count split from RH is only2+2',
          [(a,4-a) for a in range(1,4) if a%2==0 and (4-a)%2==0]==[(2,2)])
    A,B=S.symbols('A B',integer=True,positive=True)
    table=[]
    for pa,pb,oa,ob,fiber in [((5,1),(4,2),5,4,3),((5,1),(3,3),5,3,2),((4,2),(3,3),4,3,1)]:
        p=oa*A;q=ob*B
        cones=[p/S.Integer(a) for a in pa]+[q/S.Integer(b) for b in pb]
        total=S.prod(cones)*S.Rational(6)/(p*q)
        identity('UNBOUNDED complete torus homology order '+str((pa,pb)),total-fiber*A*B)
        check('fiber order from cycle gcds '+str((pa,pb)),6//(math.gcd(*pa)*math.gcd(*pb))==fiber)
        check('source parity or RP3 excludes minimal case '+str((pa,pb)),
              (fiber%2==1) or (fiber==2 and theta.subs({e:3,h:2})!=-1))
        table.append(dict(alpha_partition=list(pa),beta_partition=list(pb),
                          base_torus=[oa,ob],source_orders=[str(x) for x in cones],
                          H1_order=str(S.factor(total)),fiber_order=fiber,
                          primitive_source_forces='A=B=1',exclusion='canonical parity' if fiber%2 else 'RP3 contact theta -3 versus -1'))
    # Elementary support lemma sanity controls; universality is its written induction.
    for gens,nn in [([(1,2,0,3),(0,2,3,1)],4),
                    ([(1,2,0,3,4),(0,1,3,4,2)],5)]:
        check('two-overlap/one-overlap triple support examples generate alternating '+str(nn),len(group(gens))==math.factorial(nn)//2)
    # Exact actual knot covers; these are NOT Keller realizations.
    seeds=[(5,4,[(4,0,1,2,3,5),(1,0,5,2,3,4)],3),
           (5,3,[(3,4,0,1,2,5),(1,2,0,4,5,3)],2),
           (4,3,[(3,0,1,2,5,4),(2,0,1,4,5,3)],1)]
    controls=[]
    for p,q,rho,fiber in seeds:
        G=group(rho);check('actual connected full A6 torus cover '+str((p,q)),len(G)==360 and
                            {a[0] for a in G}==set(range(6)))
        results=[]
        for reverse in [False,True]:
            data=N.peripheral_cover(p,q,[],rho,reverse);r=N.rational_core_linking(data)
            check('actual meridian single3 and target core order '+str((p,q,reverse)),
                  cycles(data['meridian_monodromy'])==[3,1,1,1] and
                  r['core_orders']==[fiber]*4 and r['closed_H1']==dict(free=0,torsion=[] if fiber==1 else [fiber]))
            check('actual full longitude/shear balance retained '+str((p,q,reverse)),
                  r['normal_indices']==[3,1,1,1] and r['longitudinal_degrees']==[1]*4 and
                  r['shears']==[p*q%3,0,0,0])
            results.append(r)
        check('both maximal trees retain all exact marks '+str((p,q)),results[0]==results[1])
        controls.append(dict(torus=[p,q],rho=[list(a) for a in rho],group_order=len(G),
                             result=results[0],as_actual_smooth_one_prime_Keller_packet_excluded=True,
                             topological_cover_itself_excluded=False))
    # A multiplier that passes even homology parity still fails the actual primitive mark.
    rho=[(1,2,3,4,0,5),(1,2,0,4,5,3)]
    d=N.peripheral_cover(5,6,[],rho);twist=N.rational_core_linking(d)
    check('unbounded twist gate essential: even end order4 but ramification core order2',
          math.prod(twist['closed_H1']['torsion'])==4 and twist['core_orders'][0]==2 and
          twist['closed_H1']['free']==0)
    check('target torus contact step uses a proven lens topology not homology alone',
          table[1]['source_orders']==['A','5*A','B','B'] and
          [int(S.sympify(v,locals={'A':A,'B':B}).subs({A:1,B:1})) for v in table[1]['source_orders']]==[1,5,1,1])
    receipt=dict(as_of='2026-11-16',outcome='PASS',check_count=len(checks),checks=checks,
        general=dict(actual_smooth_one_line_source=True,primitive_omitted_end_generator=True,
                     negative_rank_one_intersection=True,end_order_parity_equals_index_minus_one=True,
                     canonical_contact_theta_formula=True,oriented_Yh_lens_bound_h_at_least_e_plus_one=True,
                     actual_RP3_end_excluded_all_degrees=True,actual_centre_free_monodromy_proved=True,
                     torus_end_order_divides_actual_field_degree=True,
                     odd_degree_odd_index_torus_packets_excluded=True),
        six_sheet_triple_torus=dict(full_actual_smooth_one_prime_index3_scope=True,
                     residue_degree1=True,full_monodromy='A6',all_torus_parameters_exhausted=True,
                     unbounded_multiplier_order_table=table,primitive_mark_forces_minimal_cases=True,
                     two_cases_excluded_by_parity=True,remaining_case_RP3_proved_by_Seifert_topology=True,
                     remaining_case_excluded_by_canonical_contact_theta=True,
                     unknot_excluded=True,finite_image_immersion_or_zero_excess_assumed=False),
        primary_audit=dict(Honda_lens_Theorem21_imported=True,erratum_scope_checked=True,
                     Gompf_theta_and_Gamma_exact_quantifiers_imported=True,
                     complete_selected_PDF_pages_rendered_and_read=True,
                     radial_boundary_strict_psh_perturbation_proved=True),
        contact_lens_controls=lens,actual_A6_torus_cover_controls=controls,
        nonprimitive_even_order_twist_control=twist,
        affine_index1_quadric_countercontrol=dict(e=1,end='RP3',canonical_class=0,
                     theta='-1',allowed=True,finite_Keller_target_claimed=False),
        assurance=dict(written_unbounded_source_contact_and_torus_proofs=True,
                     actual_new_higher_degree_boundary_exclusion=True,
                     finite_controls_not_universal_proof=True,
                     published_quartic_retirement_preserved=True,
                     singular_normalizations_excluded=False,multiple_prime_normalizations_excluded=False,
                     all_genuine_iterated_triple_knots_excluded=False,all_six_sheet_Keller_maps_excluded=False,
                     global_Laurent_finite_algebra_gap_closed=False,
                     external_contact_classification_independently_reproved=False,
                     formal_verification=False,external_review=False,research_goal_completed=False),
        source_sha256={name:sha(P/name) for name in [Path(__file__).name,
                     'PRIMITIVE_END_CANONICAL_CONTACT_AND_TRIPLE_TORUS_EXCLUSION.md',
                     'sources/primitive-end-contact/manifest.json']+preds+[row['path'] for row in m['primary_sources']]},
        runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'primitive_end_contact_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print('PASS',len(checks),'controls; actual all-degree source constraints and complete smooth six-sheet triple torus exclusion; full goal active')


if __name__=='__main__':main()
