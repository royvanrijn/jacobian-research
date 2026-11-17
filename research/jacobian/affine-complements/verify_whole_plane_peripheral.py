"""PF111 exact controls: universal scope rests on written Thom/linking/cable proofs."""
from pathlib import Path
import hashlib,json,math,platform
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
import whole_plane_peripheral as W
import full_peripheral_linking as F

P=Path(__file__).resolve().parent


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def ident(name,value):check(name,S.factor(value)==0)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    for predecessor in ['planar_circle_bundle_receipt.json','full_peripheral_linking_receipt.json','fibre_transitive_spin_receipt.json']:
        old=json.loads((P/predecessor).read_text());check('predecessor PASS '+predecessor,old['outcome']=='PASS')
        for name,h in old['source_sha256'].items():check('predecessor seal '+predecessor+' '+name,sha(P/name)==h)
    manifest=json.loads((P/'sources/whole-plane-peripheral/manifest.json').read_text())
    for item in manifest['sources']:
        f=P/item['path'];check('exact context '+f.name,sha(f)==item['sha256'] and f.stat().st_size==item['bytes'])
    check('no new deep classification or false fresh PDF audit',not manifest['new_deep_classification_imported'] and not manifest['new_primary_PDF_audit_claimed'])
    check('six-sheet application honestly credited as prior maximal-index exclusion',manifest['six_sheet_application_already_excluded_by_retained_maximal_index_gate'])
    # Algebraic identities in the arbitrary-rank source/contact proof.
    n,z=S.symbols('n z',positive=True,integer=True)
    ident('UNBOUNDED normal source canonical theta',-z-2*(1+n)-3*(-n)-(n-2-z))
    b,M,t,s,u,v=S.symbols('b M t s u v',integer=True)
    r=1+M*t;k=b+M*s
    ident('UNBOUNDED Bezout meridian residue',(u-b*v-1)-(r*u-k*v-1)+M*(t*u-s*v))
    ident('UNBOUNDED longitude residue',r*(b-k)+r*M*s)
    R,T,a,e=S.symbols('R T a e',nonzero=True)
    ident('UNBOUNDED shear recursion composes',T*T*a*(1-R*R)/e+a*(1-T*T)/e-a*(1-R*R*T*T)/e)
    ident('UNBOUNDED mutual linking rescales',T*T*(R*R*z)-(R*T)**2*z)
    check('old root sector fails on the six-sheet ramified orbit',math.gcd(5,1*1-1)==5)
    check('new b0 root is fibre transitive on both orbits',math.gcd(1,0)==1 and math.gcd(5,-1)==1)
    # Compact exceptional Picard classes vanish at the end; strict D still generates.
    Q=S.Matrix([[-10,-5],[-5,-3]]);exceptional=Q.inv()*S.Matrix([0,1]);j=S.Matrix([2,1])
    check('formal negative rank-two lattice retains compact exceptional class',Q[0,0]<0 and Q.det()==5 and exceptional==S.Matrix([1,-2]) and Q*exceptional==S.Matrix([0,1]))
    check('formal canonical vector satisfies exceptional adjunction and positive energy',(j.T*exceptional)[0]==0 and (j.T*(-Q).inv()*j)[0]==S.Rational(2,5))
    normal_module=W.linking_module([[S.Rational(3,5)]],generating_cores_certified=True)
    check('one strict omitted knot detects full discriminant after exceptionals vanish',normal_module['normal_filling_H1']==dict(free=0,torsion=[5]) and normal_module['row_denominators']==[5])
    # Verify the seed representation and its group by finite closure, not a label guess.
    rho=[W.projective_action(0,-1,1,0),W.projective_action(0,-1,1,1)]
    check('exact six projective images',rho==[(5,4,2,3,1,0),(4,2,3,1,5,0)])
    identity=tuple(range(6));group={identity};queue=[identity]
    while queue:
        h=queue.pop()
        for g in rho:
            q=tuple(g[h[i]] for i in range(6))
            if q not in group:group.add(q);queue.append(q)
    check('projective generator closure has exact order60',len(group)==60)
    seeds=[]
    for reverse in [False,True]:
        data=W.six_sheet_cover(reverse=reverse);out=F.rational_core_linking(data);gen=W.core_generation(data)
        bs=data['boundaries'];rows=data['relation_rows'];ms=[x['meridian'] for x in bs];ls=[x['longitude'] for x in bs]
        check('actual full peripheral types '+str(reverse),[(x['normal_index'],x['longitudinal_degree'],x['shear']) for x in bs]==[(1,1,0),(5,1,1)])
        check('actual seed rational matrix '+str(reverse),out['matrix']==[['-5','1'],['1','0']])
        check('actual seed row defects '+str(reverse),out['core_orders']==[2,2] and out['rational_row_denominators']==[1,1] and out['cleared_relation_torsion_orders']==[2,2])
        check('actual seed group and non-generation '+str(reverse),gen['closed_H1']==dict(free=0,torsion=[4]) and gen['exterior_H1']==dict(free=2,torsion=[2]) and gen['unseen_quotient']==dict(free=0,torsion=[2]) and not gen['all_cores_generate'])
        for label,extra,want in [('exterior',[],dict(free=2,torsion=[2])),('closed',ms,dict(free=0,torsion=[4])),('cores killed',ms+ls,dict(free=0,torsion=[2]))]:
            check('independent unreduced Smith '+str((reverse,label)),W.group_of_columns(S.Matrix(rows+extra).T)==want)
        seed=dict(rho=rho,meridian_monodromy=data['meridian_monodromy'],longitude_monodromy=data['longitude_monodromy'],
                  peripheral=out,generation=gen,integer_relation_rows=rows,integer_meridian_rows=ms,integer_longitude_rows=ls)
        seeds.append(seed)
    check('both maximal trees agree on complete marked seed',seeds[0]['peripheral']==seeds[1]['peripheral'] and seeds[0]['generation']==seeds[1]['generation'])
    seed=seeds[0];mu=seed['meridian_monodromy'];lam=seed['longitude_monodromy']
    check('longitude equals inverse meridian',lam==F.inverse(mu))
    # The full integral module is reconstructible only after generation is proved.
    try:W.linking_module(S.Matrix(seed['peripheral']['matrix']),generating_cores_certified=False)
    except ValueError:check('genuine non-generating seed cannot invoke the complete-F reconstruction',True)
    else:raise AssertionError('uncertified core generation accepted')
    abstract=dict(relation_rows=[[0]],boundaries=[dict(meridian=[2],longitude=[1],normal_index=2,longitudinal_degree=1,shear=1)])
    abstract_gen=W.core_generation(abstract);abstract_F=F.rational_core_linking(abstract)
    check('same half-linking generating control certified',abstract_gen['all_cores_generate'] and abstract_gen['closed_H1']==dict(free=0,torsion=[2]) and abstract_gen['exterior_H1']==dict(free=1,torsion=[]))
    module=W.linking_module([[S.Rational(1,2)]],generating_cores_certified=True)
    check('fractional meridian-lattice index retained',module['basis']==S.Matrix([[S.Rational(1,2)]]) and module['normal_filling_H1']==dict(free=0,torsion=[2]))
    for b0,c0 in [(1,0),(0,1),(1,1),(1,-2),(3,1),(-2,3)]:
        row=[b0*2+c0]
        check('arbitrary integral filling from F matches exact peripheral control '+str((b0,c0)),W.further_filling(module,[(b0,c0)])==F.reduced_h1(abstract['relation_rows'],[row]))
    # Ordinary integer linking, zero closed homology and fillings with free summands.
    integral=W.linking_module([[-5,1],[1,0]],generating_cores_certified=True)
    check('integral F in generating regime forces ZERO closed homology',integral['normal_filling_H1']==dict(free=0,torsion=[]))
    check('same integral F cannot encode actual non-generating Z4 seed',integral['normal_filling_H1']!=seed['generation']['closed_H1'])
    unlink=W.linking_module(S.zeros(2),generating_cores_certified=True)
    check('longitude filling can retain free homology',W.further_filling(unlink,[(0,1),(0,1)])==dict(free=2,torsion=[]))
    cable_controls=[]
    for cables in [[(6,5)],[(11,5)],[(6,5),(11,5)],[(16,85)]]:
        data=W.six_sheet_cover(cables);out=F.rational_core_linking(data);gen=W.core_generation(data);theory=W.six_sheet_transfer(cables)
        check('new root sector preserves full monodromy '+str(cables),data['meridian_monodromy']==mu and data['longitude_monodromy']==lam)
        check('all-parameter predicted F agrees with direct cable lift '+str(cables),S.Matrix(out['matrix'])==S.Matrix(theory['matrix']))
        check('exact filled group/core transfer agrees '+str(cables),out['closed_H1']==theory['closed_H1'] and out['core_orders']==theory['core_orders'])
        check('outside-prime obstruction persists through full cable lift '+str(cables),gen['exterior_H1']==theory['exterior_H1'] and gen['unseen_quotient']==theory['unseen_quotient'] and gen['required_outside_prime_count']==1)
        check('six-sheet known ramification leaves no repair budget '+str(cables),6-1-5*1==theory['remaining_Euler_budget']==0 and theory['required_outside_prime_count']>theory['remaining_Euler_budget'])
        cable_controls.append(dict(cables=cables,direct_peripheral=out,direct_generation=gen,unbounded_formula=theory))
    even=next(x for x in cable_controls if x['cables']==[(6,5)])
    check('zero individual defects DO NOT certify generation',even['direct_peripheral']['cleared_relation_torsion_orders']==[1,1] and even['direct_generation']['unseen_quotient']==dict(free=0,torsion=[4]))
    for bad in [[(7,5)],[(6,6)],[(6,10)],[(1,5)]]:
        try:W.six_sheet_transfer(bad)
        except ValueError:check('out-of-sector cable rejected '+str(bad),True)
        else:raise AssertionError('invalid cable accepted')
    # Existing PF101 genuine S4 control remains a counterexample to unrestricted reconstruction.
    old_data=F.peripheral_cover(2,3,[],[(1,0,2,3),(0,2,3,1)])
    old_out=F.rational_core_linking(old_data);old_gen=W.core_generation(old_data)
    check('same F1/2 has a distinct genuine integral-cover group',old_out['matrix']==[['1/2']] and old_gen['closed_H1']==dict(free=0,torsion=[8]) and not old_gen['all_cores_generate'])
    check('torsion-exterior/unseen quotient orders agree in both genuine seeds',math.prod(old_gen['exterior_H1']['torsion'])==math.prod(old_gen['unseen_quotient']['torsion'])==2)
    text=' '.join((P/'WHOLE_PLANE_THOM_LATTICE_AND_COMPLETE_PERIPHERAL_OBSTRUCTION.md').read_text().lower().split())
    for phrase in ['integral isomorphism','exact rank','generate all','compact resolution exceptionals','all-omitted-','torsion-free','partial coverage','at least g','no new deep','any finite list','r_j=1 modulo5','k_j=0 modulo5','zero excess','extra unramified','not every six-sheet','remain open']:
        check('written scope '+phrase,phrase in text)
    files=['WHOLE_PLANE_THOM_LATTICE_AND_COMPLETE_PERIPHERAL_OBSTRUCTION.md','whole_plane_peripheral.py','verify_whole_plane_peripheral.py','sources/whole-plane-peripheral/manifest.json','planar_circle_bundle_receipt.json','full_peripheral_linking_receipt.json','fibre_transitive_spin_receipt.json']+[x['path'] for x in manifest['sources']]
    result=dict(claim='PF2RULED111',outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={f:sha(P/f) for f in dict.fromkeys(files)},
        theorem=dict(all_field_degrees=True,all_normal_actual_sources=True,whole_marked_plane_required=True,
                     integral_boundary_Picard_Chern_isomorphism=True,exact_resolved_rank_r_plus_m=True,
                     all_normal_source_intersection_lattices_negative=True,all_actual_end_H1_finite=True,
                     strict_omitted_end_knots_generate_all_H1=True,
                     all_omitted_knots_covered_exterior_torsion_free=True,
                     complete_F_recovers_integral_peripheral_module_if_cores_generate=True,
                     partial_coverage_outside_prime_budget=True,
                     arbitrary_b_fibre_transitive_cable_sector_written=True,
                     complete_six_sheet_colored_iterated_family_excluded=True,
                     arbitrary_cable_length_and_stated_parameters=True),
        exact_seed=seed,cable_controls=cable_controls,
        assurance=dict(written_unbounded_proof=True,finite_sweep_proves_universal_theorem=False,
                       new_deep_classification_imported=False,all_source_boundary_components_retained=True,
                       arbitrary_cover_F_determines_full_torsion=False,
                       zero_row_defects_certify_generation=False,coverage_hypothesis_omitted=False,
                       additional_unramified_images_assumed_absent=False,
                       smooth_source_assumed_for_six_sheet_exclusion=False,
                       six_sheet_application_is_new_exclusion_class=False,
                       all_six_sheet_Keller_maps_excluded=False,uncolored_knots_alone_excluded=False,
                       all_Lp1_source_controls_excluded=False,
                       global_finite_algebra_closed=False,formal_verification=False,external_review=False,
                       external_review_deferred=True,research_goal_completed=False),
        environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'whole_plane_peripheral_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',len(checks),'PF111 controls; complete normal peripheral obstruction and six-sheet iterated-family exclusion, full goal active')


if __name__=='__main__':main()
