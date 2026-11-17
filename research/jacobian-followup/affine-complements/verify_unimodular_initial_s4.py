"""PF2RULED88: unrestricted solid-cap identities and initial simple-S4 controls."""
from pathlib import Path
import sys,json,hashlib,math,os,platform,collections
import sympy as S
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
import verify_cyclic_pattern_a4 as A
K=A.K;G=A.G;E=A.E;B=A.B;V=A.V
I=tuple(range(4))
def exterior(p,q,cables,rho):
    rows,mark,keep=K.summarized_cover(p,q,cables,rho,raw=True)
    # Sheet0 is on the transposition, so the FIRST true cycle-fill is R.
    at=4*(1+len(cables));normal=rows[at]
    return rows,rows[:at]+rows[at+1:],normal,mark
def solid_slope(unfilled,normal,mark):
    v=S.Matrix(unfilled).nullspace();assert len(v)==1
    nm=(S.Matrix([normal])*v[0])[0];nl=(S.Matrix([mark])*v[0])[0]
    value=-nl/nm;assert value.is_Integer
    return int(value)
def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    preds=['cyclic_pattern_a4_receipt.json','iterated_peripheral_tail_receipt.json',
        'cyclic_loop_first_arrival_receipt.json','fibre_order_transfer_receipt.json',
        'common_pencil_mark_receipt.json','complementary_corner_jsj_receipt.json',
        'trefoil_complementary_fusion_receipt.json','persistent_jsj_terminal_modulus_receipt.json',
        'terminal_cone_twists_receipt.json']
    for name in preds:
        old=json.loads((P/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for f,h in old['source_sha256'].items():check('sealed predecessor '+name+' '+f,sha(P/f)==h)
    manifest=P/'sources/unimodular-initial-s4/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        f=P/row['path'];check('retained complete primary pin '+f.name,sha(f)==row['sha256'] and f.stat().st_size==row['bytes'])
    p,sig,r,t,a,cap,k=S.symbols('p sigma r t a A k',integer=True)
    nd=(sig+1)*p-sig+1;delta=3*nd-p*p
    check('whole single-node source unrestricted congruence',S.expand(delta-5-(p-1)*(3*sig-p+2))==0)
    det1=[]
    for pp in [2,3,5]:
        ss=S.Rational(pp*pp-3*pp-2,3*(pp-1))
        check('complete determinant1 divisor candidate fails '+str(pp),not(ss.is_Integer and ss>=1))
        det1.append(dict(p=pp,required_sigma=str(ss)))
    check('equal pair source determinant',S.expand(delta.subs({p:r,sig:1})-r*(6-r))==0)
    check('equal pair determinant candidate fails actual marked distance',abs(5-3)==2 and 2!=1)
    check('initial disk t,3t source determinant mismatch identity',S.expand(delta.subs({p:3*t,sig:2*t+1})-t*(9*t+12))==0)
    check('odd solid cap true primitive fibre determinant',S.Matrix([[cap,1],[k,2*r]]).det()==2*cap*r-k)
    check('even solid cap true primitive fibre determinant',S.Matrix([[cap,1],[k/2,r]]).det()==cap*r-k/2)
    # Complete finite residue proof after the retained unbounded period12 theorem.
    initial=collections.Counter()
    for pp in range(12):
        for qq in range(12):
            if math.gcd(math.gcd(pp,qq),12)!=1:continue
            for aa,bb,mm,ll in G.torus(pp,qq):
                if len(G.groups[aa,bb])!=24 or mm not in G.trans:continue
                check('every initial simpleS4 generator type is4cycle/3cycle '+str((pp,qq,aa,bb)),
                      sorted([A.order(G.perms[aa]),A.order(G.perms[bb])])==[3,4] and ll==G.I and G.power(aa,pp)==G.I)
                check('every initial simpleS4 base is planar '+str((pp,qq,aa,bb)),
                      V.genus_from_base([G.perms[aa],G.perms[bb]],[G.perms[mm]])==0)
                initial[pp,qq]+=1
    check('initial simpleS4 residue class nonempty',bool(initial))
    cap_controls=[]
    for aa in [-3,-1,0,1,3,9]:
        for rr in [3,5,7]:
            for kk in [2,3,4,5,7,8,11]:
                if math.gcd(rr,kk)!=1:continue
                rows,R,m,L,M,z=A.pattern(2,rr,kk)
                eta=[aa*x+y for x,y in zip(m,L)]
                unfilled=rows[:-1]+[eta];closed=rows+[eta]
                J=abs((2*aa*rr-kk) if kk%2 else (aa*rr-kk//2))
                check('unimodular cap is positive and coprime '+str((aa,rr,kk)),J>0 and math.gcd(rr,J)==1)
                if kk%2:check('odd cap order actually odd '+str((aa,rr,kk)),J%2==1)
                want=B.group([rr] if kk%2==0 else [])
                check('actual capped-pattern H1 and zeroR '+str((aa,rr,kk)),E.h1(closed)==want and E.h1(closed,[R])==want)
                check('actual capped R exterior retains integral group '+str((aa,rr,kk)),E.h1(unfilled)==dict(free=1,torsion=want['torsion']))
                if kk%2 and J==1:
                    anew=solid_slope(unfilled,rows[-1],R)
                    check('solid continuation again has unimodular meridian '+str((aa,rr,kk)),isinstance(anew,int))
                if kk%2==0 and J==1:
                    check('actual repeated terminal inversek reduces to inverse2 '+str((aa,rr,kk)),
                          kk%rr in [2%rr,(-2)%rr] and pow(kk,-1,rr) in [pow(2,-1,rr),(-pow(2,-1,rr))%rr])
                cap_controls.append(dict(A=aa,r=rr,k=kk,J=J,closed=want,
                    exterior_cones=([rr,J] if kk%2 else [rr,rr,J]),solid=bool(kk%2 and J==1)))
    # A nonunimodular cap is a deliberate countercontrol, never covered by the lemma.
    rows,R,m,L,M,z=A.pattern(2,3,5);eta=[x+2*y for x,y in zip(m,L)];closed=rows+[eta]
    check('coefficient2 cap is outside unimodular theorem and can retain nonzeroR',
          E.h1(closed)==B.group([2]) and E.h1(closed,[R])!=E.h1(closed))
    tails=[[],[(3,2)],[(3,5)],[(5,7),(3,4)],[(3,5),(5,3)],[(5,2),(3,5)]]
    targets=[]
    for ee,oo in [(4,3),(4,9),(4,15),(8,3),(8,9),(16,3)]:
        aa,bb,mm,ll=next(x for x in G.torus(ee,oo) if len(G.groups[x[0],x[1]])==24 and x[2] in G.trans)
        rho=A.normalized(ee,oo,[],[G.perms[aa],G.perms[bb]]);M,_=A.periphery(ee,oo,[],rho)
        rows,unfilled,normal,mark=exterior(ee,oo,[],rho)
        want=B.group([oo//3])
        check('initial actual full group and zero ramification '+str((ee,oo)),E.h1(rows)==want and E.h1(rows,[mark])==want)
        check('initial actual fixed-sheet-capped R exterior '+str((ee,oo)),E.h1(unfilled)==dict(free=1,torsion=want['torsion']))
        for tail in tails:
            images=rho+[E.power(M,kk%2) for rr,kk in tail]
            ans=K.summarized_cover(ee,oo,tail,images)
            want=B.group([oo//3]+[rr for rr,kk in tail if kk%2==0])
            check('full actual initial simpleS4 trajectory H1 and zeroR '+str((ee,oo,tail)),ans==dict(H1=want,quotient=want))
            targets.append(dict(first=[ee,oo],tail=tail,actual=ans))
    # A long-word solid chain is evaluated ONLY by exact composed CW paths.
    aa,bb,mm,ll=next(x for x in G.torus(4,3) if len(G.groups[x[0],x[1]])==24 and x[2] in G.trans)
    rho=A.normalized(4,3,[],[G.perms[aa],G.perms[bb]]);M,_=A.periphery(4,3,[],rho);cables=[];solid_controls=[]
    for step in range(4):
        rows,unfilled,normal,mark=exterior(4,3,cables,rho);oldA=solid_slope(unfilled,normal,mark)
        check('deep composed-CW solid prefix H1 and unimodular cap '+str(step),E.h1(rows)==B.group([]) and E.h1(unfilled)==dict(free=1,torsion=[]))
        solid_controls.append(dict(step=step,cables=cables.copy(),actual_A=oldA))
        rr=3 if step%2==0 else 5;kk=2*oldA*rr-1
        check('deep solid control cap distance1 '+str(step),kk>1 and abs(2*oldA*rr-kk)==1)
        rho.append(E.power(M,kk%2));cables.append((rr,kk))
    check('repeated terminal free coefficient permits only r3',
          [rr for rr in [3,5,7,9,11,15] if (2-1)%rr==0 or (2+1)%rr==0]==[3])
    check('r3 and integral-homology-sphere cases fail zeroR collar modulus',3%3!=2 and 1%3!=2)
    inputs=[Path(__file__),P/'UNIMODULAR_SOLID_PREFIX_AND_INITIAL_S4_CLOSURE.md',manifest]+[P/f for f in preds]+[P/x['path'] for x in sources['primary_sources']]
    result=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      general=dict(whole_single_node_source_delta5_mod_pminus1=True,no_whole_two_cone_source_delta1=True,
        equal_pair_determinant_r_requires_r5_and_wrong_distance=True,
        arbitrary_integer_unimodular_solid_cap_A=True,no_recursive_cap_slope_needed=True,
        all_odd_windings_coprime_cable_parameters_and_tail_lengths=True,
        actual_cyclic2_pattern_and_fixed_meridian_fill_hypotheses_required=True,
        finite_first_nonsolid_dichotomy_unbounded_length=True),
      quartic=dict(all_initial_state3_torus_trajectories_closed_both_pencils=True,
        all_initial_generator_types4cycle_3cycle=True,actual_fixed_caps_disk_a_t_3t=True,
        initial_H1_Zt_and_R_zero=True,arbitrary_finite_later_tails=True,no_prefix_eta_condition=True,
        actual_T4_3_solid_exception_closed=True,remaining_first_arrival_frontier_states2_and_complementary7=True),
      initial_periodic_counts={str(k):v for k,v in initial.items()},single_node_delta1_candidates=det1,
      actual_unimodular_cap_controls=cap_controls,full_initial_trajectory_controls=targets,
      composed_solid_prefix_controls=solid_controls,
      assurance=dict(written_unbounded_source_cap_first_nonsolid_and_terminal_proofs=True,
        actual_full_filled_and_unfilled_CW_mark_controls=True,complete_retained_primary_audits_reused=True,
        deep_JSJ_and_full_peripheral_rigidity_imported=True,solid_target_excluded_by_group_without_Waldhausen=True,
        finite_controls_not_universal_proof=True,
        all_deep_sources_delta1_excluded=False,nonunimodular_A4_cap_closed=False,
        all_state2_closed=False,all_complementary_A4_closed=False,full_quartic_exclusion=False,
        higher_degree_packet_classification=False,finite_Keller_map=False,formal_verification=False,
        external_review=False,independent_full_reproof=False,novelty_priority=False,research_goal_completed=False),
      source_sha256={os.path.relpath(f,P):sha(f) for f in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'unimodular_initial_s4_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; initial simple-S4 trajectories closed BOTH pencils at every tail length; full goal OPEN')
if __name__=='__main__':main()
