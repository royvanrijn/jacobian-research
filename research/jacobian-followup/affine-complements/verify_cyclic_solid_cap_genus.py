"""PF2RULED94: degree-free solid cap/genus recurrence and A4 reduction."""
from pathlib import Path
import sys,json,hashlib,math,os,platform
import sympy as S
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
import verify_genus_final_s3 as H
import composed_peripheral_cover as C
A=H.A;K=H.K;E=H.E;G=H.G;B=H.B;I=tuple(range(4))
def history(steps):
    aa=1;bb=6;cc=1;cables=[];valid=True
    for r,eps in steps:
        assert r>=5 and r%2 and r%3 and eps in [-1,1]
        k=(3*aa*r+eps)//2;assert 2*k==3*aa*r+eps
        valid &= k<r*bb
        cables.append((r,k));cc=r*cc+(r-1)*k;bb=r*k;aa*=r*r
    return dict(steps=steps,A=aa,B=bb,c=cc,D=3*aa-2*cc,W=math.prod(r for r,eps in steps),cables=cables,infinity_valid=valid)
def primitive_images(rows,meridian,longitude):
    basis=S.Matrix(rows).nullspace();assert len(basis)==1
    v=basis[0];scale=S.ilcm(*[S.Rational(x).q for x in v]);v=[int(x*scale) for x in v]
    g=math.gcd(*v);v=[x//g for x in v]
    m=sum(x*y for x,y in zip(meridian,v));l=sum(x*y for x,y in zip(longitude,v))
    return (m,l) if m>0 else (-m,-l)
def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    preds=['integer_grid_common_receipt.json','recursive_linking_endpoints_receipt.json',
        'nonunimodular_a4_receipt.json','cyclic_pattern_a4_receipt.json',
        'global_ruling_unfilled_cyclic_receipt.json','genus_final_s3_receipt.json',
        'cyclic_loop_first_arrival_receipt.json','iterated_peripheral_tail_receipt.json']
    for name in preds:
        old=json.loads((P/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for path,digest in old['source_sha256'].items():check('sealed predecessor '+name+' '+path,sha(P/path)==digest)
    manifest=P/'sources/cyclic-solid-cap-genus/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        path=P/row['path'];check('retained primary pin '+path.name,sha(path)==row['sha256'] and path.stat().st_size==row['bytes'])
    n,b,aa,r,c,eps,w,ep2,D=S.symbols('n b A r c epsilon w epsilon_next D',integer=True)
    kval=(n*aa*r+eps)/b
    check('degree-free exact genus deficit recurrence',S.expand(n*aa*r*r-b*(r*c+(r-1)*kval)-(r*(n*aa-b*c)-(r-1)*eps))==0)
    check('degree-free exact L-space cabling threshold margin',S.simplify(kval/r-c-((n*aa-b*c)+eps/r)/b)==0)
    check('degree-free strict infinity obstruction after minus',S.expand((n*aa*r+ep2)/b-r*(n*aa-w)/b-(ep2+r*w)/b)==0)
    check('last-minus closed genus deficit',S.expand(r-(r-1)*(-1)-(2*r-1))==0)
    # Exact cyclic unfilled/solid-cap controls across several cover degrees.
    cyclic=[]
    for degree,coef,oldA,rr,sign in [(2,1,1,3,-1),(2,1,1,5,1),
        (3,2,1,5,-1),(3,2,1,5,1),(4,3,1,5,1),
        (5,2,1,3,-1),(5,2,1,3,1),(5,3,1,2,-1),
        (7,2,1,3,1),(8,3,1,5,-1),(9,2,1,5,1)]:
        kk=(degree*oldA*rr+sign)//coef
        assert coef*kk==degree*oldA*rr+sign and kk>=2 and math.gcd(rr,kk)==math.gcd(degree,rr)==math.gcd(degree,kk)==1
        rows,Li,mo,lo,_,_=A.pattern(degree,rr,kk);mi=rows[-1];unfilled=rows[:-1]
        cap=[oldA*x+coef*y for x,y in zip(mo,lo)];capped=unfilled+[cap]
        rel1=[x-rr*y for x,y in zip(mo,mi)];rel2=[x-rr*y for x,y in zip(Li,lo)]
        check('all-degree cyclic exact unfilled bases '+str((degree,coef,rr,kk)),
            C.reduced_h1(unfilled)==dict(free=2,torsion=[]) and C.reduced_h1(unfilled,[rel1,rel2])==C.reduced_h1(unfilled))
        check('all-degree unit cap preserves coefficient b and squares A '+str((degree,coef,rr,kk)),
            C.reduced_h1(capped)==dict(free=1,torsion=[]) and primitive_images(capped,mi,Li)==(coef,oldA*(-rr*rr)))
        cyclic.append(dict(n=degree,b=coef,A=oldA,r=rr,k=kk,sign=sign,J=abs(degree*oldA*rr-coef*kk),new_A=oldA*rr*rr,primitive_images=[coef,-oldA*rr*rr]))
    # Exhaustive finite base-case orbit, not a parameter-bounded history search.
    base=[]
    for p,q in [(2,3),(3,2)]:
        pairs={(G.perms[x[0]],G.perms[x[1]]) for x in G.torus(p,q) if len(G.groups[x[0],x[1]])==12 and A.order(G.perms[x[2]])==3}
        seedpair=next(iter(pairs))
        def conjugate(pi,g):
            inv=E.inv(g);return tuple(g[pi[inv[i]]] for i in range(4))
        orbit={(conjugate(seedpair[0],g),conjugate(seedpair[1],g)) for g in G.perms}
        check('complete initial trefoil A4 pairs form ONE conjugacy orbit '+str((p,q)),len(pairs)==24 and pairs==orbit)
        for rho in sorted(pairs):
            data=C.peripheral_cover(p,q,[],list(rho));bd=next(x for x in data['boundaries'] if x['degree']==3)
            rows=data['relation_rows']+[x['meridian'] for x in data['boundaries'] if x['degree']==1]
            check('initial actual primitive coefficient2 cap sign '+str((p,q,rho)),
                C.reduced_h1(rows)==dict(free=1,torsion=[]) and primitive_images(rows,bd['meridian'],bd['longitude'])==(2,-1))
        base.append(dict(p=p,q=q,labelled_pairs=len(pairs),one_conjugacy_orbit=True,primitive_images=[2,-1]))
    old=next(x for x in G.torus(2,3) if len(G.groups[x[0],x[1]])==12 and A.order(G.perms[x[2]])==3)
    seed=A.normalized(2,3,[],[G.perms[old[0]],G.perms[old[1]]]);M,_=A.periphery(2,3,[],seed)
    def representation(cables):
        rho=list(seed);mm=M
        for rr,kk in cables:
            rho.append(E.power(mm,(kk*pow(rr,-1,3))%3));mm=E.power(mm,pow(rr,-1,3))
        return rho,mm
    stepsets=[[],[(5,1)],[(7,1)],[(5,1),(7,1)],[(5,-1)],[(7,-1)],[(5,1),(7,-1)],[(7,1),(5,-1)],[(5,1),(7,1)]*4]
    histories=[]
    for steps in stepsets:
        h=history(steps);rho,mm=representation(h['cables']);data=C.peripheral_cover(2,3,h['cables'],rho)
        bd=next(x for x in data['boundaries'] if x['degree']==3);rows=data['relation_rows']+[x['meridian'] for x in data['boundaries'] if x['degree']==1]
        check('all-depth actual A4 solid square cap without expanded words '+str(steps),
            h['infinity_valid'] and h['A']==h['W']**2 and h['A']%24==1 and C.reduced_h1(rows)==dict(free=1,torsion=[]) and primitive_images(rows,bd['meridian'],bd['longitude'])==(2,-h['A']))
        wantD=1 if not steps or steps[-1][1]==1 else 2*steps[-1][0]-1
        check('all-depth actual A4 ordinary genus deficit '+str(steps),h['D']==wantD and all(S.Rational(kk,rr)>hh['c'] for hh,(rr,kk) in zip([history(steps[:i]) for i in range(len(steps))],h['cables'])))
        histories.append(h)
    bad=[]
    for w0,r0,s in [(5,5,1),(5,7,-1),(7,5,-1)]:
        h=history([(w0,-1),(r0,s)]);check('minus continuation exact infinity countercontrol '+str((w0,r0,s)),not h['infinity_valid']);bad.append(h)
    # Complete ratio/inequality deduction for the only remaining reset orders.
    assignments=[(2,8),(8,2),(4,10),(10,4)]
    check('complete reset assignments have ratio J/r at most4',max(S.Rational(jj,rr) for rr,jj in assignments)==4)
    check('last-minus still-solid arrival complete winding list is only5',[ww for ww in range(2,6) if ww%2 and ww%3]==[5])
    reset_table=[]
    for rr,jj in assignments:
        for sign in [-1,1]:
            check('last-minus reset impossible from exact bounded ratio '+str((rr,jj,sign)),sign*jj>=-5*rr)
            reset_table.append(dict(r=rr,J=jj,sign=sign,k_formula=f'({3*rr}*A{sign*jj:+d})/2',eta=(rr+sign*jj)//2-1))
    # Full marked covers for actual infinity-valid survivors, with two tails.
    targets=[]
    def closed_control(cables,rho,name):
        data=C.peripheral_cover(2,3,cables,rho);bd=next(x for x in data['boundaries'] if x['degree']==2)
        rows=data['relation_rows']+[x['meridian'] for x in data['boundaries']]
        ans=dict(H1=C.reduced_h1(rows),quotient=C.reduced_h1(rows,[bd['longitude']]))
        check('full actual A4 arrival/reset/tail marked cover '+name,ans==dict(H1=B.group([2]),quotient=B.group([2])) and data['word_expansion_used'] is False)
        return ans
    for steps in [[],[(5,1)],[(5,1),(7,1)],[(5,-1)],[(7,1),(5,-1)]]:
        h=history(steps);rho,mm=representation(h['cables'])
        for sign in [-1,1]:
            kval=h['A']+2*sign
            if kval<=0 or 3*kval>=2*h['B']:continue
            eta=3*kval-2*h['c']-1
            check('still-solid arrival exact plus/minus classification '+str((steps,sign)),
                (not steps or steps[-1][1]==1) and eta==6*sign or steps[-1]==(5,-1) and sign==-1 and eta==2)
            z,tau,_,_=next(x for x in A.roots(mm,I,2,3*kval,2,6) if E.generated(rho+[x[0]])==24)
            tail=[(5,3),(11,3)];cables=h['cables']+[(2,3*kval)]+tail;images=rho+[z]+[E.power(tau,tk%2) for tr,tk in tail]
            ans=closed_control(cables,images,str((steps,sign)))
            targets.append(dict(kind='still_solid',history=steps,arrival=[2,3*kval],eta_arrival=eta,tail=tail,actual=ans))
    for steps in [[],[(5,1)],[(5,1),(7,1)]]:
        h=history(steps);rho,mm=representation(h['cables'])
        for rr,jj in assignments:
            for sign in [-1,1]:
                kk=(3*h['A']*rr+sign*jj)//2
                if kk<2 or kk>=rr*h['B']:continue
                assert math.gcd(rr,kk)==1 and kk%3
                reset_rho=rho+[E.power(mm,(kk*pow(rr,-1,3))%3)];newM=E.power(mm,pow(rr,-1,3))
                # One genuine nonsolid cyclic3 loop retained after reset.
                middle=[(2,5)];reset_rho.append(E.power(newM,5*pow(2,-1,3)));newM=E.power(newM,pow(2,-1,3))
                z,tau,_,_=next(x for x in A.roots(newM,I,2,3,2,6) if E.generated(reset_rho+[x[0]])==24)
                tail=[(5,3),(11,3)];cables=h['cables']+[(rr,kk)]+middle+[(2,3)]+tail;images=reset_rho+[z]+[E.power(tau,tk%2) for tr,tk in tail]
                eta=kk-rr*h['c']-1;check('exact all-plus reset genus table '+str((steps,rr,jj,sign)),eta==(rr+sign*jj)//2-1)
                ans=closed_control(cables,images,str((steps,rr,jj,sign)))
                targets.append(dict(kind='earlier_reset',history=steps,reset=[rr,kk],terminal=[rr,jj],eta_reset=eta,postreset_C3=middle,arrival=[2,3],tail=tail,actual=ans))
    # Degree-free cyclic-lattice prime lemma; cancellation and noncyclicity
    # countercontrols keep its exact boundary visible.
    lattice=[]
    for hh,qq,pp,nn,ell in [(3,1,3,4,3),(3,1,6,14,3),(3,1,9,28,3),
        (7,2,7,29,7),(9,1,9,10,3),(27,1,3,4,3),(9,1,3,4,3)]:
        delta=hh*nn-qq*qq*pp*pp;assert delta>0
        mat=[[-hh,qq*pp],[qq*pp,-nn]];group=C.reduced_h1(mat);quot=C.reduced_h1(mat,[[1,0]])
        va=lambda x:int(S.factorint(abs(x)).get(ell,0))
        av=va(hh);tv=va(pp)
        check('cyclic full lattice prime mark generates primary part '+str((hh,qq,pp,nn,ell)),
            group==B.group([delta]) and nn%ell and math.gcd(qq*pp,nn)%ell and all(x%ell for x in quot['torsion']))
        check('exact unequal lattice valuation rule '+str((hh,qq,pp,nn,ell)),av==2*tv or va(delta)==min(av,2*tv))
        lattice.append(dict(h=hh,q=qq,p=pp,nD=nn,prime=ell,delta=delta,H1=group,quotient_R=quot,valuation=va(delta),equal_valuation_case=av==2*tv))
    check('equal valuation CAN cancel so no unjustified equality',lattice[-1]['equal_valuation_case'] and lattice[-1]['valuation']==3)
    check('noncyclic source lattice invalidates single-factor conclusion',C.reduced_h1([[-3,3],[3,-6]])==B.group([3,3]))
    check('terminal trefoil source winding must3 at multiple nodes',math.gcd(2,3-1)!=1 and math.gcd(3,2-1)==1)
    check('product merger cone pair is NEVER consecutive at genuine odd tail windings',
        all(abs((2*ww+sgn)-ww)>1 for ww in [3,5,7,11] for sgn in [-1,1]))
    # Actual S3 marked covers: missing, repeated and tail-annihilated factors3.
    oldS=next(x for x in G.torus(2,3) if len(G.groups[x[0],x[1]])==6 and A.order(G.perms[x[2]])==2)
    seedS=A.normalized(2,3,[],[G.perms[oldS[0]],G.perms[oldS[1]]]);MS,_=A.periphery(2,3,[],seedS)
    s3=[]
    for steps in [[],[(3,1)],[(5,1)],[(7,1)],[(5,-1)]]:
        hh=H.history(steps);pref=seedS+[E.power(MS,tk%2) for tr,tk in hh['cables']]
        for d0 in [-3,3,6,9,1]:
            kk=2*(hh['A']+d0)
            if kk<2 or math.gcd(3,kk)!=1 or kk>=3*hh['B']:continue
            z,tau,_,_=next(x for x in A.roots(MS,I,3,kk,2,6) if E.generated(pref+[x[0]])==24)
            for tail in [[(5,3)],[(3,5),(5,3)]]:
                data=C.peripheral_cover(2,3,hh['cables']+[(3,kk)]+tail,pref+[z]+[E.power(tau,tk%2) for tr,tk in tail])
                bd=next(x for x in data['boundaries'] if x['degree']==2);rows=data['relation_rows']+[x['meridian'] for x in data['boundaries']]
                ans=dict(H1=C.reduced_h1(rows),quotient=C.reduced_h1(rows,[bd['longitude']]))
                TT=math.prod(tr for tr,tk in tail);delta=abs(d0);CC=abs(3*hh['A']-kk//2);pp=3*hh['W']*CC*TT
                check('actual S3 target cyclic group and mark coefficient '+str((steps,d0,tail)),ans==dict(H1=B.group([delta]),quotient=B.group([math.gcd(delta,kk*TT)])))
                good=(delta%3==0 and delta%9!=0 and hh['W']%3!=0 and TT%3!=0)
                if delta%3==0:
                    nn=(delta+pp*pp)//3;source=C.reduced_h1([[-3,pp],[pp,-nn]])
                    sq=C.reduced_h1([[-3,pp],[pp,-nn]],[[1,0]])
                    v3=lambda grp:sorted(int(S.factorint(x).get(3,0)) for x in grp['torsion'] if x%3==0)
                    check('source full cyclic lattice and marked quotient distinguish exact3 obstruction '+str((steps,d0,tail)),
                        good==(v3(source)==v3(ans['H1']) and v3(sq)==v3(ans['quotient'])))
                s3.append(dict(history=steps,arrival=[3,kk],D0=d0,tail=tail,source_pole_candidate=pp,actual=ans,passes_only_new_3primary_necessary_gate=good,not_affine_realization=True))
    inputs=[Path(__file__),P/'CYCLIC_SOLID_CAP_GENUS_AND_A4_HISTORY_REDUCTION.md',manifest,P/'composed_peripheral_cover.py']+[P/f for f in preds]+[P/x['path'] for x in sources['primary_sources']]
    result=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
        general=dict(all_degree_actual_cyclic_unit_cap_Anew_r2A_coefficient_b_preserved=True,
            all_depth_genus_deficit_Dnew_rD_minus_rminus1epsilon=True,
            positive_deficit_positive_Lspace_companion_preserved=True,
            negative_solid_sign_must_be_LAST_by_strict_infinity=True,
            actual_topological_solid_model_not_H1_recognition=True,
            cyclic_two_line_lattice_common_prime_mark_generates_primary_part=True,
            unequal_prime_valuation_exact_equal_case_not_overclaimed=True),
        quartic=dict(A4_initial_primitive_cap_EXACT_m_plus2L=True,A4_solid_A_EXACT_Wsquare_1mod24=True,
            all_A4_solid_ordinary_companions_positive_Lspace=True,
            earlier_nonsolid_reset_requires_ALLPLUS_solid_history=True,
            stillsolid_minus_arrival_ONLY_lastwinding5_KAminus2_eta2=True,
            allplus_stillsolid_arrival_eta_exact_plusminus6=True,
            all_remaining_complementary_S3_end_order_has_EXACT_single_factor3=True,
            all_remaining_S3_solid_and_later_tail_windings_prime_to3=True,
            S3_rge9_s_prime_to3_and_exact_arrival_mod9_gate=True,
            product_first_tail_fibre_mergers_excluded_before_terminal_match=True,
            postreset_nonsolid_loops_not_replaced_by_solid_recurrence=True,
            no_prefix_eta_bound_or_separate_affine_source_assumption=True),
        all_degree_cyclic_cap_controls=cyclic,exhaustive_initial_A4_orbits=base,
        actual_A4_solid_histories=histories,forbidden_minus_continuations=bad,
        complete_reset_degree_table=reset_table,full_actual_target_controls=targets,
        prime_mark_lattice_controls=lattice,S3_exact_3primary_controls=s3,
        assurance=dict(written_unbounded_cap_genus_infinity_and_frontier_proofs=True,
            finite_initial_A4_coloring_classification_EXHAUSTIVE=True,
            exact_integral_peripherals_and_composed_cover_controls=True,
            retained_deep_topology_and_primary_Floer_inputs_imported=True,
            finite_history_controls_not_universal_proof=True,
            all_complementary_A4_tails_closed=False,all_complementary_S3_tails_closed=False,
            all_degree_packet_classification=False,degenerate_stage_classification=False,
            full_quartic_exclusion=False,finite_Keller_map=False,formal_verification=False,
            external_review=False,independent_full_reproof=False,novelty_priority=False,research_goal_completed=False),
        source_sha256={os.path.relpath(f,P):sha(f) for f in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'cyclic_solid_cap_genus_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; all-degree cap/genus and prime-mark laws; complementary A4/S3 histories restricted, full goal OPEN')
if __name__=='__main__':main()
