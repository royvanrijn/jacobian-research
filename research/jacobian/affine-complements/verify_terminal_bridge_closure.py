"""PF2RULED95: terminal bridge laws, S3 closure and still-solid A4 closure.

Finite exact controls support the written unbounded proof. They do not
recognize affine realizations or independently prove imported topology.
"""
from pathlib import Path
import sys,json,hashlib,math,os,platform
import sympy as S
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
import exact_poly_cluster as Q
import composed_peripheral_cover as C
import verify_cyclic_solid_cap_genus as H
import verify_nonunimodular_a4 as N
A=H.A;K=H.K;E=H.E;G=H.G;B=H.B;I=tuple(range(4))
def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def identity(name,value):check(name,S.expand(value)==0)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    preds=['cyclic_solid_cap_genus_receipt.json','integer_grid_common_receipt.json',
        'recursive_linking_endpoints_receipt.json','nonunimodular_a4_receipt.json',
        'genus_final_s3_receipt.json','cyclic_pattern_a4_receipt.json',
        'persistent_jsj_terminal_modulus_receipt.json','trefoil_complementary_fusion_receipt.json',
        'terminal_cone_twists_receipt.json','iterated_peripheral_tail_receipt.json']
    for name in preds:
        old=json.loads((P/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for path,digest in old['source_sha256'].items():check('sealed predecessor '+name+' '+path,sha(P/path)==digest)
    manifest=P/'sources/terminal-bridge-closure/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        f=P/row['path'];check('primary pin '+f.name,sha(f)==row['sha256'] and f.stat().st_size==row['bytes'])
    fresh=sources['primary_sources'][0]
    check('new primary characteristic-difference audit keeps signed splice convention separate',len(fresh['complete_HTML_sections_read'])==2 and fresh['new_deep_theorem_import'] is False)
    t,kap,kp,h,dd,qq,w,f,Di,Ga=S.symbols('t kappa kappa_prime h delta_prime q w f D_in gamma',integer=True)
    din=t*kp-(t+kap-1);dout=1-kp
    identity('terminal adjacent residue exact',t*dout+(din+kap-1))
    cp=(f*t+1)*dd-h*(f*(t+kap-1)+1);cm=(f*t-1)*dd-h*(f*(t+kap-1)-1)
    identity('positive signed divisor elimination',cp.subs(dd,h*kp)-h*(dout-f*din)*(-1))
    identity('negative signed divisor elimination',cm.subs(dd,h*kp)-h*(dout+f*din))
    strip=lambda delta,hh,tt,kk:S.Rational(delta+hh*(tt-1)*(tt+kk),tt*tt)
    check('terminal2 free4 strips determinant2 to EXACT5',strip(2,3,2,4)==5)
    check('both earlier reset pairs strip to EXACT8',strip(2,3,2,8)==strip(2,3,4,10)==8)
    def constants(hh,tt,kk,dp,ff):return [int(z.subs({h:hh,t:tt,kap:kk,dd:dp,f:ff})) for z in [cp,cm]]
    pairs={'still_solid':constants(3,2,4,5,2),'reset_2_8':constants(3,2,8,8,1),'reset_4_10':constants(3,4,10,8,1)}
    check('complete exact terminal divisor constants',pairs=={'still_solid':[-8,-12],'reset_2_8':[-6,-16],'reset_4_10':[-2,-12]})
    def divisors(vals):return sorted({int(x) for val in vals for x in S.divisors(abs(val)) if x>=2 and x%3})
    check('still-solid tail has NO odd genuine winding unit3',not any(x%2 for x in divisors(pairs['still_solid'])))
    lists=[divisors(pairs['reset_2_8']),divisors(pairs['reset_4_10'])]
    check('complete post-reset first winding divisor lists',lists==[[2,4,8,16],[2,4]])
    check('zero divisor constant branch explicitly unbounded',constants(3,2,4,6,1)[0]==0)
    for ww in [11,17,101]:
        kpv=ww+2
        check('zero-constant necessary congruences admit arbitrarily chosen large w controls '+str(ww),
            (6-3*kpv)%ww==0 and ((1-kpv)-(2*kpv-5))%ww==0)
    # Same complete blowup rows as the old rational-expression helper.
    z,U=S.symbols('z U');shadow=[]
    for sig,support in [(2,{10:1,5:1}),(4,{14:1,13:1}),(4,{8:2,6:1,5:3})]:
        new=Q.local_cluster(sig,support)
        old=Q.C.local_cluster(sum(c*z**-j for j,c in support.items())+z**(1-sig)*U,z**sig,0,z,U)
        check('independent symbolic and polynomial rational full-row cluster agreement '+str((sig,support)),new==old)
        shadow.append(dict(sigma=sig,support=support,blowups=len(new['rows']),all_rows_identical=True))
    # Full integer graphs, including non-self-inverse internal cone residues.
    graphs=[]
    for sig,support in [(2,{10:1,5:1}),(10,{50:1,18:1,13:1}),
        (22,{110:1,38:1,25:1}),(20,{60:1,50:1,46:1,23:1})]:
        row=N.contact(sig,support);graph=Q.graph(sig,support)
        drops=row['drops'];js=[row['p']]+[int(x[2]) for x in drops];gs=[int(x[1]) for x in drops]
        want=[(js[i]-js[i+1])//gs[i] for i in range(len(gs))]
        check('all exact full graph bridge differences '+str((sig,support)),[x['D'] for x in graph['bridges']]==want)
        check('full graph source contact diagonal exact '+str((sig,support)),graph['source']==[[-3,row['p']],[row['p'],-row['nD']]])
        for i,node in enumerate(graph['ordered_nodes'][1:-1],1):
            order=row['windings'][i];leg=graph['legs'][node][0];inv=pow(want[i-1],-1,order)
            check('complete internal cone cofactor is inverse OUTER bridge '+str((sig,support,node)),leg['order']==order and leg['cofactor']%order in {inv,(-inv)%order})
        tt=row['windings'][-1];kk=row['kappa'];mid=js[-2];kpv=mid//tt+1-sig//tt
        check('terminal inward bridge EXACT t*kappa_prime minus contact gap '+str((sig,support)),graph['bridges'][-1]['D']==tt*kpv-(tt+kk-1))
        graph['contact']=row;graphs.append(graph)
    check('completed delta2 source charts pass contact but new still-solid bridge condition excludes',
        all(x['contact']['delta']==2 and x['contact']['kappa']==4 for x in graphs[1:3]) and
        all(x['contact']['windings'][-2] not in divisors(pairs['still_solid']) for x in graphs[1:3]))
    # S3 ALL nonempty ordinary solid histories reduce to a finite divisor.
    solid_cases=[]
    for eps in [-1,1]:
        vals=[int(3*eps-sgn*(1-eps)) for sgn in [-1,1]]
        choices=sorted({int(x) for val in vals for x in S.divisors(abs(val)) if x>=3 and x%2})
        check('complete first ordinary S3 loop divisor choices '+str(eps),choices==([5] if eps==-1 else [3]))
        for ww in choices:
            kval=2*ww+eps;inner=6*ww-kval
            obstruction='W divisible3 forbidden PF94' if ww%3==0 else 'terminal bridge not2mod3'
            check('every first-loop divisor candidate exactly excluded '+str((eps,ww)),ww%3==0 or inner%3!=2)
            solid_cases.append(dict(epsilon=eps,w=ww,k=kval,D_in=inner,exclusion=obstruction))
    # Degree-zero source terminal strip is not an arbitrary parameter cutoff.
    kbase=[kk for kk in range(2,18,2) if math.gcd(3,kk)==1]
    kept=[kk for kk in kbase if kk//2-1>0 and (kk//2-1)%9==6]
    check('complete depth-zero r3 infinity bound and terminal strip forceONLYk14',kept==[14])
    check('depth-zero k14 inward fibre fails terminal3 residue',(18-14)%3!=2)
    check('depth-zero product roots already fail exact source terminal strip',all((kk//2-1)%9!=6 for kk in [4,8]))
    unit=[]
    for ss in [5,7,11,17,25]:
        for sign in [-1,1]:
            rr=3*ss;kk=6*ss+2*sign;v=pow(-kk,-1,rr);co=(-v)%ss;inner=6*rr-kk
            check('actual signed S3 unit-cap normalized cone and inner bridge '+str((ss,sign)),(2*sign*co)%ss==1 and inner==12*ss-2*sign)
            check('all unit-cap S3 root survivors excluded by inward residue or divisor '+str((ss,sign)),inner%3!=2 or (ss%3!=0 and 3%ss!=0 and 9%ss!=0))
            unit.append(dict(s=ss,zeta=sign,r=rr,k=kk,cone_coefficient=co,D_in=inner))
    # Actual full target groups that pass old homology but fail the NEW local law.
    controls=[]
    oldS=next(x for x in G.torus(2,3) if len(G.groups[x[0],x[1]])==6 and A.order(G.perms[x[2]])==2)
    seedS=A.normalized(2,3,[],[G.perms[oldS[0]],G.perms[oldS[1]]]);MS,_=A.periphery(2,3,[],seedS)
    def marked(cables,rho):
        data=C.peripheral_cover(2,3,cables,rho);bd=next(x for x in data['boundaries'] if x['degree']==2)
        rows=data['relation_rows']+[x['meridian'] for x in data['boundaries']]
        return dict(H1=C.reduced_h1(rows),quotient=C.reduced_h1(rows,[bd['longitude']]))
    for older,aa,rr,kk in [([],1,3,14),([],1,51,100),([(5,11)],25,3,62)]:
        pref=seedS+[E.power(MS,tk%2) for tr,tk in older]
        zz,tau,_,_=next(x for x in A.roots(MS,I,rr,kk,2,6) if E.generated(pref+[x[0]])==24)
        tail=[(5,3)];ans=marked(older+[(rr,kk)]+tail,pref+[zz,E.power(tau,1)])
        d0=kk//2-aa*(rr//3);order=abs(d0);coef=(rr//3)*kk*5
        check('actual S3 old-prime-passing countercontrol still has exact marked cyclic H1 '+str((older,rr,kk)),
            order%3==0 and order%9!=0 and ans==dict(H1=B.group([order]),quotient=B.group([math.gcd(order,coef)])))
        controls.append(dict(kind='S3_newly_excluded',older=older,arrival=[rr,kk],tail=tail,actual=ans,excluded_by_written_terminal_law=True,not_affine_realization=True))
    oldA=next(x for x in G.torus(2,3) if len(G.groups[x[0],x[1]])==12 and A.order(G.perms[x[2]])==3)
    seedA=A.normalized(2,3,[],[G.perms[oldA[0]],G.perms[oldA[1]]]);MA,_=A.periphery(2,3,[],seedA)
    for steps in [[],[(5,1)],[(5,-1)]]:
        hist=H.history(steps);rho=list(seedA);mm=MA
        for tr,tk in hist['cables']:
            rho.append(E.power(mm,(tk*pow(tr,-1,3))%3));mm=E.power(mm,pow(tr,-1,3))
        kval=hist['A']+(2 if not steps else -2)
        zz,tau,_,_=next(x for x in A.roots(mm,I,2,3*kval,2,6) if E.generated(rho+[x[0]])==24)
        tail=[(11,3),(5,3)];ans=marked(hist['cables']+[(2,3*kval)]+tail,rho+[zz]+[E.power(tau,1)]*2)
        check('actual A4 still-solid two-tail countercontrol retains old marked H1 '+str(steps),ans==dict(H1=B.group([2]),quotient=B.group([2])))
        controls.append(dict(kind='A4_still_solid_newly_excluded',steps=steps,arrival=[2,3*kval],tail=tail,actual=ans,excluded_by_written_terminal_law=True,not_affine_realization=True))
    inputs=[Path(__file__),P/'TERMINAL_BRIDGES_AND_S3_SOLID_A4_CLOSURE.md',P/'exact_poly_cluster.py',P/'composed_peripheral_cover.py',manifest]+[P/f for f in preds]+[P/x['path'] for x in sources['primary_sources']]
    result=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
        general=dict(all_depth_characteristic_bridge_and_inverse_cone_law=True,
            terminal_adjacent_signed_residue_and_inward_modulus=True,
            terminal_strip_two_constant_divisor_law_no_degree_bound=True,
            zero_constant_branch_not_claimed_bounded=True,
            primitive_local_stripped_chart_not_assumed_affine=True),
        quartic=dict(ALL_incoming_S3_BOTH_pencils_ALL_finite_lengths_closed=True,
            S3_nonempty_solid_histories_closed_by_FIRST_ordinary_loop=True,
            S3_depthzero_r3_and_unitcap_rge9_closed=True,
            ALL_still_solid_A4_arrivals_BOTH_pencils_ALL_lengths_closed=True,
            A4_remainingONLY_earlier_nonsolid_reset_terminal2_8_or4_10=True,
            first_postreset_C3_winding_complete_finite_lists=True,
            no_postreset_loop_and_later_reset_tails_NOT_claimed_closed=True,
            earlier_common_class_and_all_other_incoming_closures_retained=True),
        divisor_constants=pairs,postreset_first_winding_lists=lists,
        symbolic_old_helper_full_row_controls=shadow,complete_integer_source_graphs=graphs,
        complete_first_S3_solid_loop_divisors=solid_cases,
        S3_depthzero_r3_complete_k=kbase,S3_unitcap_controls=unit,full_actual_target_controls=controls,
        assurance=dict(written_unbounded_terminal_slope_and_divisor_proofs=True,
            actual_oriented_slopes_with_independent_signs_allowed=True,
            exact_full_graph_and_composed_target_controls=True,
            primary_characteristic_convention_and_retained_deep_topology_explicit=True,
            finite_controls_not_universal_proof=True,
            all_complementary_A4_resets_closed=False,full_quartic_exclusion=False,
            higher_degree_packet_classification=False,degenerate_stage_exclusion=False,
            finite_Keller_map=False,formal_verification=False,external_review=False,
            independent_full_reproof=False,novelty_priority=False,research_goal_completed=False),
        source_sha256={os.path.relpath(f,P):sha(f) for f in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'terminal_bridge_closure_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; ALL S3 and still-solid A4 tails closed BOTH pencils; earlier reset A4 frontier OPEN')
if __name__=='__main__':main()
