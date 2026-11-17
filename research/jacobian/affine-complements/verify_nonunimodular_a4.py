"""PF2RULED89: coefficient2 cap descent and an exact all-depth terminal list."""
from pathlib import Path
import sys,json,hashlib,math,os,platform
import sympy as S
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
import verify_cyclic_pattern_a4 as A
K=A.K;G=A.G;E=A.E;B=A.B
I=tuple(range(4))
def contact(sigma,support):
    support=sorted(support,reverse=True);p=support[0]
    assert all(j>=sigma+1 for j in support)
    g=math.gcd(sigma,p);first=g;total=(sigma-g)*p;wind=[p//g];last=p;drops=[]
    for j in support[1:]:
        h=math.gcd(g,j)
        if h<g:
            total+=(g-h)*(2*p-j);drops.append((g,h,j));wind.append(g//h);g=h;last=j
    assert g==1
    n=total+2*p-sigma+1
    return dict(sigma=sigma,p=p,support=support,S=total,nD=n,delta=3*n-p*p,
        windings=[x for x in wind if x>1],kappa=last+1-sigma,root_distance=abs(p//first-3*sigma//first),drops=drops)
def companion(p,q,cables,rho):
    rows,mark,keep=K.summarized_cover(p,q,cables,rho,raw=True);at=4*(1+len(cables))
    return rows,rows[:at]+rows[at+1:],rows[at],mark
def cap_A(unfilled,m,L):
    basis=S.Matrix(unfilled).nullspace();assert len(basis)==1
    value=-2*(S.Matrix([L])*basis[0])[0]/(S.Matrix([m])*basis[0])[0]
    assert value.is_Integer and int(value)%2
    return int(value)
def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    preds=['unimodular_initial_s4_receipt.json','cyclic_pattern_a4_receipt.json',
        'iterated_peripheral_tail_receipt.json','persistent_jsj_terminal_modulus_receipt.json',
        'terminal_cone_twists_receipt.json','fibre_order_transfer_receipt.json',
        'complementary_corner_jsj_receipt.json','trefoil_complementary_fusion_receipt.json']
    for name in preds:
        old=json.loads((P/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for f,h in old['source_sha256'].items():check('sealed predecessor '+name+' '+f,sha(P/f)==h)
    manifest=P/'sources/nonunimodular-a4/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        f=P/row['path'];check('retained complete primary pin '+f.name,sha(f)==row['sha256'] and f.stat().st_size==row['bytes'])
    a,kap,np,pp,nr,q=S.symbols('a kappa nprime pprime nR q',integer=True)
    n=a*a*(np-1)+a-(a-1)*kap
    check('exact arbitrary nR q last-gcd determinant strip identity',
          S.expand(nr*n-q*q*(a*pp)**2-a*a*(nr*np-q*q*pp*pp)+nr*(a-1)*(a+kap))==0)
    check('two residual strip deltas exactly8',all(S.Rational(2+3*(aa-1)*(aa+kk),aa*aa)==8 for aa,kk in [(2,8),(4,10)]))
    check('internal cone2 makes2,4 stripped determinant parity impossible',S.Rational(2+3*(2-1)*(2+4),4)==5 and 5%2==1)
    candidates=[]
    for aa in S.divisors(8):
        if aa<2:continue
        for kk in sorted(set(S.divisors(3*aa-2)+S.divisors(3*aa+2))):
            if kk<2 or math.gcd(aa,kk)!=2:continue
            if (2-3*(aa-(aa-1)*kk))%(aa*aa):continue
            candidates.append([int(aa),int(kk)])
    check('complete finite DIVISOR list is exactly three ordered pairs',candidates==[[2,4],[2,8],[4,10]])
    b,g,QQ,p,sig=S.symbols('b g Q p sigma',integer=True,positive=True)
    lower=3*((sig+1)*p-sig+1)-p*p
    check('all-depth root minus-branch contact lower bound identity',
          S.expand(lower.subs({p:g*QQ,sig:g*b})-(g*g*QQ*(3*b-QQ)+3*g*(QQ-b)+3))==0)
    check('rho3 zero-mark source determinant2 has no square solution',3 not in {x*x%7 for x in range(7)})
    check('rho2 rootcone2 source has strict determinant lower bound',S.expand(lower.subs({p:2*g,sig:g})-(2*g*g+3*g+3))==0)
    check('one simple tail normalized cone and edge conditions force bminus2 or4',
          S.expand(2*(4*b+1)-b-2*(3*b+2)-(b-2))==0 and
          S.expand(2*(4*b+1)+b-3*(3*b+2)+4)==0)
    # Exact contact arithmetic controls replace expensive full symbolic clusters.
    source_controls=[]
    for sig0,support in [(2,[10,5]),(6,[22,9]),(8,[28,26,15]),(16,[80,28,25]),(28,[140,48,46,35]),(14,[70,24,23])]:
        row=contact(sig0,support)
        check('all-depth exact contact diagonal even '+str((sig0,support)),row['nD']%2==0)
        if row['drops']:
            aa,hh,jj=row['drops'][-1];assert hh==1
            stripped=contact(sig0//aa,[j//aa for j in support if j>jj])
            check('actual primitive local chart exact diagonal strip '+str((sig0,support)),
                  row['nD']==aa*aa*(stripped['nD']-1)+aa-(aa-1)*row['kappa'])
            check('actual primitive local chart determinant strip '+str((sig0,support)),
                  row['delta']==aa*aa*stripped['delta']-3*(aa-1)*(aa+row['kappa']))
            row['stripped']=stripped
        source_controls.append(row)
    hard=next(x for x in source_controls if x['p']==140)
    check('delta2 internal cone2 and rootdistance2 are NOT a universal obstruction',
          hard['delta']==2 and hard['root_distance']==2 and hard['windings']==[5,7,2,2] and hard['kappa']==8)
    # Reuse one genuinely completed small whole cluster to cross-check the new local edge formula.
    import verify_trefoil_complementary_fusion as F
    cluster=F.cluster_control(2,2,{10:1,5:1})
    check('actual full small cluster has predicted root5 and terminal2,4',sorted(cluster['legs'].values())==[[2,4],[5]])
    check('actual full small cluster fibre edge determinant is pminusj5',len(cluster['bridges'])==1 and abs(cluster['bridges'][0]['transported_fibre'][1])==5)
    # A non-self-inverse root coefficient distinguishes inverse(b) from b.
    zc,wc,Tc,Uc=S.symbols('z w T U');C=F.C;J=F.J
    rr=J.retag(C.local_cluster(zc**-1+Tc,zc**2,0,zc,Tc),'R','Hx')
    dd=J.retag(C.local_cluster(wc**-14+wc**-13+wc**-3*Uc,wc**4,0,wc,Uc),'D','Hy')
    raw=C.assemble([rr,dd]);mat=raw['full'];vertices=set(range(mat.rows))-set(raw['index'])
    arrow=raw['labels'][rr['rows'][-1]['parents'][0]];spine=J.rooted_nodes(mat,vertices,arrow)
    root_leg=spine['legs'][spine['root']][0];legpath=root_leg['vertices']
    co=int(abs(mat.extract(legpath[1:],legpath[1:]).det()))
    check('full actual first-cone coefficient is INVERSE b up to sign',root_leg['order']==7 and co in [pow(2,-1,7),(-pow(2,-1,7))%7])
    check('full root coefficient control excludes wrong b direction',co not in [2,(-2)%7])
    root_coefficient_control=dict(sigma=4,p=14,support=[14,13],alpha=7,b=2,cofactor=co,weights=[int(mat[i,i]) for i in legpath])
    # Actual initial A4 cap and several solid cyclic3 repetitions use COMPOSED paths only.
    aa,bb,mm,ll=next(x for x in G.torus(2,3) if len(G.groups[x[0],x[1]])==12 and A.order(G.perms[x[2]])==3)
    seed=A.normalized(2,3,[],[G.perms[aa],G.perms[bb]]);M,_=A.periphery(2,3,[],seed)
    rows,unfilled,m,L=companion(2,3,[],seed);initialA=cap_A(unfilled,m,L)
    check('actual initial coefficient2 cap is odd and solid',initialA in [1,7] and E.h1(unfilled)==dict(free=1,torsion=[]))
    solid_controls=[];rho=seed.copy();cables=[];currentM=M
    for step in range(3):
        rows,unfilled,m,L=companion(2,3,cables,rho);oldA=cap_A(unfilled,m,L)
        check('actual solid A4 chain retains H1Z2 and cap coefficient2 '+str(step),E.h1(rows)==B.group([2]) and E.h1(unfilled)==dict(free=1,torsion=[]))
        rr=5 if step%2==0 else 7;kk=(3*oldA*rr-1)//2
        check('actual unrestricted solid loop cap distance1 '+str(step),kk>=2 and math.gcd(rr,kk)==1 and abs(3*oldA*rr-2*kk)==1)
        solid_controls.append(dict(step=step,cables=cables.copy(),actual_A=oldA,next_loop=[rr,kk]))
        rho.append(E.power(currentM,(kk*pow(rr,-1,3))%3));currentM=E.power(currentM,pow(rr,-1,3));cables.append((rr,kk))
    targets=[]
    # Positive/zero cap controls at the final arrival while the old A4 companion is solid.
    for dK in [0,2,4]:
        Kval=initialA+dK;fr,fk=2,3*Kval;z,tau,lam,beta=next(A.roots(M,I,fr,fk,2,6))
        for tail in [[],[(5,3)],[(5,3),(11,3)]]:
            images=seed+[z]+[E.power(tau,kk%2) for rr,kk in tail]
            ans=K.summarized_cover(2,3,[(fr,fk)]+tail,images)
            check('actual solid-at-arrival target full marked H1 '+str((dK,tail)),ans==dict(H1=B.group([2]),quotient=B.group([2])))
            targets.append(dict(kind='solid_at_arrival',arrival=[fr,fk],J_final=2*abs(dK),tail=tail,actual=ans))
    for rr,JJ in [(2,4),(2,8),(8,2),(4,10),(10,4)]:
        kk=(3*initialA*rr+JJ)//2
        check('actual nonsolid reset primitive cap '+str((rr,JJ)),math.gcd(rr,kk)==1 and kk%3 and abs(3*initialA*rr-2*kk)==JJ)
        z=E.power(M,(kk*pow(rr,-1,3))%3);newM=E.power(M,pow(rr,-1,3));rho0=seed+[z]
        for fr,fk in [(2,3),(10,3)]:
            fz,tau,lam,beta=next(A.roots(newM,I,fr,fk,2,6))
            for tail in [[],[(5,3)],[(5,3),(11,3)]]:
                full=K.summarized_cover(2,3,[(rr,kk),(fr,fk)]+tail,rho0+[fz]+[E.power(tau,tk%2) for tr,tk in tail])
                want=B.group([2,fr//2]);check('actual reset arrival/tail Rzero and full H1 '+str((rr,JJ,fr,tail)),full==dict(H1=want,quotient=want))
                targets.append(dict(kind='earlier_reset',reset=[rr,kk],terminal=[rr,JJ],arrival=[fr,fk],tail=tail,actual=full))
    # Integral cap section congruence controls include even r,J, without illegal division by2.
    cap_controls=[]
    for oldA in [1,3,7]:
        for rr,JJ in [(2,4),(2,8),(8,2),(4,10),(10,4)]:
            kk=(3*oldA*rr+JJ)//2;cc,dd,gg=S.gcdex(3*rr,kk);assert gg==1;cc=int(cc);dd=-int(dd)
            beta=2*cc-oldA*dd
            check('actual signed cap coefficient integral congruence '+str((oldA,rr,JJ)),(3*rr*beta-2)%JJ==0)
            cap_controls.append(dict(A=oldA,r=rr,k=kk,J=JJ,beta=int(beta%JJ)))
    inputs=[Path(__file__),P/'NONUNIMODULAR_A4_CAPS_AND_THREE_TERMINAL_PAIRS.md',manifest]+[P/f for f in preds]+[P/x['path'] for x in sources['primary_sources']]
    result=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      general=dict(exact_last_gcd_local_chart_strip_any_nR_q=True,all_source_depths_and_no_support_bound=True,
        stripped_chart_not_assumed_separate_affine_packet=True,unramified_chart_diagonal_even=True,
        all_depth_contact_lower_bound=True,odd_coefficient2_solid_caps_left_symbolic=True,
        rootdistance_and_normalized_edge_coefficient_comparison_retained=True),
      quartic=dict(all_final_stage_A4_arrivals_closed_both_pencils=True,no_prefix_eta_condition=True,
        arbitrary_A4_companion_depth=True,remaining_A4_requires_later_simple_tails=True,
        complete_terminal_pairs_2_4_2_8_4_10=True,solid_at_arrival_requires2_4_and_at_least2_simple_tails=True,
        earlier_reset_requires2_8_or4_10_and_internal_arrival2=True,
        one_simple_tail_no_intervening_C3_loop_reset_closed=True,
        final_simple_tail_winding5mod6=True),
      ordered_terminal_pairs=candidates,exact_source_contact_controls=source_controls,
      actual_solid_C3_controls=solid_controls,full_target_controls=targets,cap_coefficient_controls=cap_controls,
      completed_small_cluster_control=dict(matrix=[[int(x) for x in row] for row in cluster['matrix'].tolist()],bridges=cluster['bridges'],legs=cluster['legs']),
      actual_root_coefficient_control=root_coefficient_control,
      assurance=dict(written_unbounded_cap_strip_root_and_three_pair_proofs=True,
        actual_full_CW_and_exact_contact_controls=True,complete_retained_primary_audits_reused=True,
        deep_JSJ_and_full_peripheral_rigidity_imported=True,finite_divisors_follow_proved_unbounded_congruences=True,
        local_countercontrols_not_affine_or_full_topological_matches=True,finite_controls_not_universal_proof=True,
        all_A4_complementary_closed=False,all_determinant2_or8_sources_excluded=False,
        all_state2_closed=False,full_quartic_exclusion=False,higher_degree_packet_classification=False,
        finite_Keller_map=False,formal_verification=False,external_review=False,independent_full_reproof=False,
        novelty_priority=False,research_goal_completed=False),
      source_sha256={os.path.relpath(f,P):sha(f) for f in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'nonunimodular_a4_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; all final-stage A4 arrivals closed; later tails reduced to three ordered terminal types; full goal OPEN')
if __name__=='__main__':main()
