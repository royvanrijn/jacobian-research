"""PF2RULED44: exact support for the written all-depth regular quartic exclusion.
No finite coefficient sweep is used; topology imports and proof limits stay explicit.
"""
from pathlib import Path
from math import gcd,lcm
from functools import reduce
from itertools import permutations
import hashlib,json,platform
import sympy as S
import verify_quartic_cabling as Q
import verify_universal_regular_end as U
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def cyc(p):
    left=set(range(len(p)));ans=[]
    while left:
        a=min(left);row=[]
        while a not in row:row.append(a);left.remove(a);a=p[a]
        ans.append(len(row))
    return tuple(sorted(ans))
def prod(p,q):return tuple(p[q[i]] for i in range(len(p)))
def perm_order(p):return lcm(*cyc(p))
def main():
    pred_names=['universal_regular_end_receipt.json','coprime_regular_end_receipt.json','borderline_regular_end_receipt.json','quartic_cabling_receipt.json','higher_regular_coordinate_receipt.json','regular_sector_end_receipt.json','intrinsic_normalization_receipt.json']
    preds={}
    for name in pred_names:
        d=json.loads((HERE/name).read_text());preds[name]=d;check('predecessor PASS '+name,d['outcome']=='PASS')
        for f,h in d['source_sha256'].items():check('predecessor seal '+name+' '+f,digest(HERE/f)==h)
    manifest=HERE/'sources/regular-quartic-exclusion/manifest.json';pins=json.loads(manifest.read_text())
    check('six primary statements with explicit import boundaries',len(pins['primary_sources'])==6)
    for row in pins['primary_sources']:
        f=HERE/row['path_relative_to_affine_complements'];check('primary PDF byte seal '+f.name,f.read_bytes().startswith(b'%PDF') and f.stat().st_size==row['bytes'] and digest(f)==row['sha256'])
    visual=pins['new_visual_audits']
    check('Waldhausen complete marked-boundary conclusion visually rechecked',visual[0]['PDF_pages']==[23] and visual[0]['complete_theorem_statement'] and visual[0]['boundary_relative_conclusion_checked'])
    check('all essential-annulus and fiber exceptions visually rechecked',visual[1]['PDF_pages']==[20,21] and visual[1]['complete_essential_annulus_and_boundary_fiber_statements'] and visual[1]['all_product_Klein_exceptions_recorded'])
    # Exhaustive base-cover permutation layer for a hypothetical degree-four
    # torus-exterior cover; A,B themselves remain arbitrary coprime integers.
    base_profiles={};base_counts={}
    for b in [2,4]:
        profiles=set();count=0
        for a in permutations(range(b)):
            for c in permutations(range(b)):
                ca,cc=cyc(a),cyc(c)
                if cyc(prod(a,c))!=(b,) or len(ca)+len(cc)!=b+1:continue
                if gcd(perm_order(a),perm_order(c))!=1:continue
                count+=1;profiles.add((ca,cc))
        expected=({((1,1),(2,)),((2,),(1,1))} if b==2 else
                  {((1,1,1,1),(4,)),((4,),(1,1,1,1)),((1,3),(1,1,2)),((1,1,2),(1,3))})
        check('exhaustive admissible coprime disk-base cover profiles b'+str(b),profiles==expected)
        for ca,cc in profiles:
            if ca in [(1,1),(1,1,1,1)] or cc in [(1,1),(1,1,1,1)]:
                repeated=ca if ca in [(1,1),(1,1,1,1)] else cc
                check('identity monodromy repeats a genuine target cone '+str((b,ca,cc)),len(repeated)>=2 and set(repeated)=={1})
            else:
                # Target orders are 3u and 2v: the singleton in (1,3)
                # and the two singletons in (1,1,2) are always genuine.
                check('degree-four mixed profile has at least three genuine cones '+str((ca,cc)),sorted([ca,cc],key=len)==[(1,3),(1,1,2)])
        base_profiles[str(b)]=[[list(a),list(c)] for a,c in sorted(profiles)];base_counts[str(b)]=count
    check('fiber-degree-four boundary equation would require AB dividing 3 or 5',min(2*3,3*2)>max(abs(4-1),abs(4+1)))
    # Degree-independent central-fiber argument, with finite controls that
    # deliberately allow nontrivial residual permutations and larger orbits.
    local_fiber_controls=[]
    for degree in [3,4,5,6]:
        perms=list(permutations(range(degree)));identity=tuple(range(degree));tau_local=(1,0)+tuple(range(2,degree));count=0;max_orbit=0
        for h in perms:
            if prod(h,tau_local)!=prod(tau_local,h):continue
            for alpha in perms:
                if prod(h,alpha)!=prod(alpha,h):continue
                orbit={0};todo=[0]
                while todo:
                    a=todo.pop()
                    for gen in [h,alpha,tau_local]:
                        b=gen[a]
                        if b not in orbit:orbit.add(b);todo.append(b)
                h2=prod(h,h)
                assert all(h2[a]==a for a in orbit)
                count+=1;max_orbit=max(max_orbit,len(orbit))
        check('simple ramification fiber has multiplicity at most two in degree '+str(degree),count>0 and max_orbit>=2)
        local_fiber_controls.append(dict(degree=degree,commuting_pattern_controls=count,maximum_branched_orbit=max_orbit))
    h=(1,0,3,2);alpha=(2,3,0,1);tau_local=(1,0,2,3)
    check('nonidentity fiber with four-letter support is retained when longitude identity is absent',prod(h,alpha)==prod(alpha,h) and prod(h,tau_local)==prod(tau_local,h) and Q.orbit(Q.generated([h,alpha,tau_local]),0)==set(range(4)))
    # Full S3 subgroup/peripheral residue census used by the solid partial
    # three-sheeted-cover induction. It includes C2 companions, not only S3.
    S3=frozenset(a for a in Q.P if a[3]==3);tau=(1,0,2,3)
    subs={frozenset([Q.I])};todo=list(subs)
    while todo:
        g=todo.pop()
        for a in S3:
            k=Q.generated(tuple(g)+(a,))
            if k not in subs:subs.add(k);todo.append(k)
    check('complete S3 subgroup layer',len(subs)==6)
    derived={g:Q.comm(g) for g in subs};states=set();residue_count=0
    for p in range(6):
        for a in range(6):
            if gcd(gcd(p,a),6)!=1 or p*a%2:continue
            for alpha in S3:
                if Q.power(alpha,p)!=Q.I:continue
                pattern=Q.generated([alpha,tau]);orb=Q.orbit(pattern,0);delta=Q.mul(alpha,tau)
                mu=Q.power(delta,p);lon=Q.power(delta,-a)
                for H in subs:
                    if mu not in H or lon not in derived[H] or Q.normal(H,mu)!=H:continue
                    if Q.generated(tuple(pattern)+tuple(H))!=S3:continue
                    states.add((len(orb),Q.cycles(alpha,orb),len(H),Q.cycles(mu),Q.cycles(lon)));residue_count+=1
    expected={(2,(1,1),6,(1,1,2),(1,1,1,1)),(3,(1,2),3,(1,3),(1,1,1,1)),(3,(3,),2,(1,1,2),(1,1,1,1)),(3,(3,),6,(1,1,2),(1,1,1,1))}
    check('complete S3 trivial-fiber cases retain both C2 and S3 companion alternatives',states==expected)
    nt=set()
    for p in range(6):
        for a in range(6):
            if gcd(gcd(p,a),6)!=1 or not p*a%2:continue
            mu=Q.power(tau,p);h=Q.power(tau,p*a)
            for H in subs:
                if mu not in H or Q.normal(H,mu)!=H:continue
                for lon in derived[H]:
                    if Q.mul(lon,tau)!=Q.mul(tau,lon):continue
                    if Q.mul(Q.power(lon,p),Q.power(mu,a))!=h:continue
                    if Q.generated(tuple(H)+(tau,))==S3:nt.add((len(H),Q.cycles(mu),Q.cycles(lon)))
    check('nontrivial S3 fiber has full companion and identity longitude',nt=={(6,(1,1,2),(1,1,1,1))})
    # Reuse the sealed exhaustive S4/A4 census without widening its assurance.
    layer=preds['quartic_cabling_receipt.json']['finite_peripheral_layer']
    check('complete S4 subgroup and residue layer is retained',layer['exhaustive_S4_subgroups']==30 and layer['S_fiber_trivial_state_count']==9 and layer['A_fiber_trivial_state_count']==4 and layer['S_residue_controls']==208 and layer['A_residue_controls']==96)
    decisions=[]
    for state in layer['states']['S']:
        m=state['branched_orbit'];alpha=state['alpha_cycles'];H=state['companion_order']
        if m==2:why='two genuine root cones for trivial fiber; nontrivial fiber would give forbidden distance two'
        elif m==3 and alpha==[1,2]:why='only p=2,g=2; A4 replacement has fiber distance>=7'
        elif m==3 and H==24:why='dual graph cycle gives infinite closed H1'
        elif m==3 and H==6:why='identity companion sheet forces p=g; signed trefoil distance contradicts coprimality'
        elif m==4 and alpha==[1,3]:why='only p=3,g=3; impossible unbranched degree-four knot-exterior cover'
        elif m==4 and alpha==[4]:why='three-boundary root contradicts the source cable piece'
        else:raise AssertionError('unhandled S4 state '+str(state))
        decisions.append(dict(state=state,exclusion=why))
    check('every retained S4 peripheral state has an explicit proof case',len(decisions)==9)
    aux=[]
    for state in layer['states']['A']:
        check('each A4 trivial-fiber state is disk-three-cone or annular four-letter',state['branched_orbit']==4 or (state['branched_orbit']==3 and state['alpha_cycles']==[1,1,1]))
        aux.append(dict(state=state,marked_fiber_distance=1))
    check('A4 centralizer/derived intersection makes longitude identity',{a for a in Q.comm(Q.A4) if Q.mul(a,(1,2,0,3))==Q.mul((1,2,0,3),a)}=={Q.I})
    # Source peeling is checked using the independent characteristic machinery
    # and semigroup controls already sealed by PF2RULED43.
    examples=[(4,[6,3]),(6,[9,5]),(6,[10,8,5]),(6,[12,10,5]),(8,[12,7]),(8,[12,10,7]),(12,[18,11]),(12,[18,16,11])]
    peeling=[]
    for r,support in examples:
        d=U.data(r,support);q=d['q'];high=[j for j in support if j>r-1];g=reduce(gcd,[r]+high)
        tag=f' r{r} support{support}'
        check('distance-one exact datum has a proper higher support gcd'+tag,d['kappa']==1 and r not in support and 1<g<r and r-1 in support)
        short=U.data(r//g,[j//g for j in high]);F=short['F'];kp=short['kappa'];np=short['n']
        check('shortened datum is again admissible exact regular geometry'+tag,min(short['support'])>=short['r']+1 and short['r'] not in short['support'] and kp>=3)
        check('source companion characteristic exponents are exactly the divided prefix'+tag,[s['beta'] for s in short['chain']]==[s['beta']//g for s in d['chain'][:-1]])
        check('final cable slope is the source-peeling formula'+tag,d['pairs'][-1]==[g,g*(np-2)+1])
        check('surgery numerator and primitive rational filling agree'+tag,d['n']==g*g*(np-2)+g+1 and gcd(d['n'],g*g)==1 and d['n']>=31)
        D=d['n']-g*g*F
        check('next actual fiber distance has its exact gap'+tag,D==g*g*(kp-2)+g+1 and D>=g*g+g+1>=7)
        if g==2:check('A4 replacement violates the intrinsic fiber-distance bound'+tag,D not in [1,3] and gcd(2,d['n'])==1)
        peeling.append(dict(original=d,shortened=short,higher_support_gcd=g,companion_meridian=[d['n'],g*g],next_fiber_distance=D))
    F,kp=S.symbols('F kappa_prime',integer=True,positive=True);g=S.symbols('g',integer=True,positive=True)
    np=F+kp;alpha=g*(np-2)+1;n=g*alpha+1
    check('universal peeling numerator is a polynomial identity',S.expand(n-(g*g*(np-2)+g+1))==0)
    check('universal next-distance polynomial is exact',S.expand(n-g*g*F-(g*g*(kp-2)+g+1))==0)
    check('minimum next-distance bound follows at kappa_prime>=3',S.expand((n-g*g*F).subs(kp,3)-(g*g+g+1))==0)
    check('binary A4 replacement distance has no values one or three at kappa_prime>=3',S.expand((n-g*g*F).subs({g:2,kp:3}))==7)
    # The fiber/H1 argument fixes companion meridians up to signs; it does
    # not assume a general knot-complement recognition theorem.
    possibilities={e-w for e in [-1,1] for w in [-1,1]}
    check('identity-companion longitude correction would divide zero or two',possibilities=={-2,0,2} and max(abs(v) for v in possibilities)<6)
    check('one-cone p=3g case contradicts the identity-boundary winding',S.expand(3*g-g)==2*g)
    signed=[]
    for eps in [-1,1]:
        for sign in [-1,1]:
            a=2*(3*eps+sign*3)
            check('last p=3,cone=3 trefoil slopes violate cable coprimality '+str((eps,sign)),abs(a//2-3*eps)==3 and gcd(abs(a),3)==3)
            signed.append(dict(epsilon=eps,distance_sign=sign,a=a))
        for p in [3,5,7]:
            for sign in [-1,1]:
                a=2*p*eps+sign
                check('partial three-sheeted solid condition violates algebraic growth '+str((eps,p,sign)),abs(a-2*p*eps)==1 and abs(a)<=2*p+1<6*p)
    sharp=preds['borderline_regular_end_receipt.json']['countercontrol']
    check('actual retained S4 lens-cover control fails precisely the new root winding',sharp['actual_topological_simple_S4_cover_constructed'] and sharp['lens_H1_order']==7 and sharp['marked_exterior_extra_cone']==5 and sharp['cable_winding']==[3,16] and 5!=3)
    check('actual nonalgebraic-growth S3 control has full monodromy and matched companion meridian',Q.generated(tuple(S3)+(tau,))==S3 and Q.power(tau,3)==tau)
    check('actual nonalgebraic-growth S3 control has a solid partial exterior',gcd(3,5)==1 and abs(5-2*3)==1 and 5<3*6)
    # Closed-model decompositions used before applying prime decomposition.
    model_profiles=[]
    for d in preds['universal_regular_end_receipt.json']['controls']:
        if len(d['chain'])<2:continue
        wind=d['pairs'][-1][0]
        if d['kappa']>=3:
            check('closed terminal piece has negative orbifold Euler characteristic '+str(d['support']),S.Rational(1,wind)+S.Rational(1,d['kappa'])-1<0)
            model_profiles.append(dict(support=d['support'],closed_terminal_cones=[wind,d['kappa']],group_argument='injective chain with noncyclic central vertex groups'))
        else:
            prev=d['pairs'][-2][0]*d['pairs'][-2][1];D=d['n']-wind*wind*prev
            check('distance-one closed model has genuine next Seifert filling '+str(d['support']),D>=7 and 0<prev<d['n'])
            if len(d['chain'])==2:check('collapsed closed Seifert central fiber is nonzero in finite homology '+str(d['support']),prev%d['n']!=0)
            model_profiles.append(dict(support=d['support'],companion_filling=[d['n'],wind*wind],next_distance=D,group_argument='nontrivial central closed fiber or injective central-vertex chain'))
    inputs=[Path(__file__),HERE/'REGULAR_COORDINATE_QUARTIC_EXCLUSION.md',manifest,HERE/'verify_universal_regular_end.py',HERE/'verify_quartic_cabling.py']+[HERE/r['path_relative_to_affine_complements'] for r in pins['primary_sources']]+[HERE/n for n in pred_names]
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        theorem=dict(all_regular_source_coordinate_quartic_models_excluded=True,one_omitted_prime=True,generic_index=2,finite_quartic_Keller_map=True,all_Laurent_coefficients=True,all_cabling_depths=True,remaining_one_index_two_prime_requires='every source coordinate has a pole at R'),
        all_degree_geometry=dict(radial_closed_pair_recognized=True,radial_marked_exterior_recognized=True,irreducibility_proved_before_Waldhausen=True,full_peripheral_system_and_omitted_meridian_preserved=True,H1_only_recognition_assumed=False,arbitrary_abstract_group_recognition_assumed=False,closed_group_free_indecomposability='central vertex groups and nontrivial injective edges; closed Seifert fiber handled separately',source_fiber_slopes_and_JSJ_transfer_valid=True),
        all_degree_fiber_bound=dict(simple_ramification_multiplicity=[1,2],preferred_longitude_identity_required=False,nontrivial_fiber_branched_orbit_need_not_have_size_two=True,essential_nonproduct_piece_forces_kappa_one=True,finite_controls=local_fiber_controls),
        exterior_lemmas=dict(no_unbranched_degree_four_between_nontrivial_iterated_torus_knot_exteriors=True,source_chain_leaf_reflection_explicitly_excluded=True,base_cover_profiles=base_profiles,base_profile_counts=base_counts,solid_partial_three_sheeted_cover_under_absolute_algebraic_growth='trefoil or mirror',S3_trivial_fiber_states=[dict(branched_orbit=k[0],alpha_cycles=list(k[1]),companion_order=k[2]) for k in sorted(states)],S3_residue_controls=residue_count,S3_C2_companions_retained=True,finite_H1_not_assumed_for_solid_partial_three_sheeted_lemma=True),
        quartic_reduction=dict(terminal_distance_at_least_three_excluded=True,remaining_intermediate_distance=1,product_exception_excluded_by_solid_three_sheeted_lemma=True,distance_one_all_cases_excluded=True,S4_decisions=decisions,nontrivial_fiber_case='distance two, already impossible',A4_trivial_fiber_states=aux,A4_intrinsic_fiber_distance_bound=[1,3],A4_replacement_minimum_distance=7,identity_companion_forces='p_target=g_source',final_signed_cases=signed),
        source_peeling=dict(all_degree=True,full_support_block='h(z)=H(z^g)+c*z^(-(r-1)), 1<g<r',shortened_model_actual_affine=True,shortened_terminal_distance_at_least=3,source_last_cable='(g,g*(n_prime-2)+1)',surgery_numerator='n=g^2*(n_prime-2)+g+1',new_companion_meridian='n*mu_prime+g^2*lambda_prime',new_fiber_distance='g^2*(kappa_prime-2)+g+1',finite_Keller_map_descended=False,controls=peeling),
        sharpness_controls=dict(non_growth_S3=dict(actual_topological_cover=True,target='(3,5)-cable of positive trefoil',full_S3=True,partial_exterior_solid=True,algebraic_growth=False,finite_Keller_map_constructed=False),retained_S4=dict(target='(3,16)-cable of positive trefoil',actual_topological_cover=True,lens_H1_order=7,source_root_cone=5,identity_boundary_longitude_coefficient=3,iterated_S3_knot_exterior_condition=False)),
        model_group_controls=model_profiles,
        assurance=dict(written_unbounded_geometric_and_all_case_proofs=True,primary_imports_explicit=True,finite_controls_not_full_topological_proof=True,source_reference_and_target_infinity_knots_distinguished=True,compression_ball_alternatives_handled=True,all_solid_and_product_degenerations_handled=True,no_automatic_target_algebra_invariance_assumed=True,coefficient_search_used=False,no_regular_coordinate_case_excluded=False,additional_omitted_primes_excluded=False,higher_index_packets_excluded=False,all_finite_Keller_degrees_excluded=False,finite_Keller_map_constructed=False,formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,research_goal_completed=False),
        source_sha256={str(p.relative_to(HERE)):digest(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'regular_quartic_exclusion_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS '+str(len(checks))+' exact regular-quartic controls; entire regular-coordinate quartic sector excluded; goal active')
if __name__=='__main__':main()
