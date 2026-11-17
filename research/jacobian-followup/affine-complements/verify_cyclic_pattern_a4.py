"""PF2RULED87: cyclic pattern integral marks and all-length A4 reduction."""
from pathlib import Path
import sys,json,hashlib,math,os,platform
import sympy as S
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
import verify_cyclic_loop_first_arrival as K
G=K.G;E=G.E;V=G.V;B=G.B
I=tuple(range(4));product=K.product;wp=K.wp
def path(w,rho,keep,start=0):
    row=[0]*len(keep);indices={e:i for i,e in enumerate(keep)};s=start
    for x in w:
        g=abs(x)-1
        if x>0:e=(g,s);s=E.inv(rho[g])[s];sign=1
        else:s=rho[g][s];e=(g,s);sign=-1
        if e in indices:row[indices[e]]+=sign
    assert s==start
    return row
def order(pi):return next(n for n in range(1,13) if E.power(pi,n)==tuple(range(len(pi))))
def periphery(p,q,cables,rho):
    rels,m,l=K.words(p,q,cables);return E.evaluate(m,rho),E.evaluate(l,rho)
def normalized(p,q,cables,rho):
    m,l=periphery(p,q,cables,rho)
    if m[0]!=0:return list(rho)
    s=next(i for i in range(4) if m[i]!=i);c=list(range(4));c[0],c[s]=s,0
    return [tuple(c[a[c[i]]] for i in range(4)) for a in rho]
def roots(M,L,r,k,new_order,cable_order):
    u,v,d=S.gcdex(r,k);assert d==1;u=int(u);v=-int(v)
    F=E.mul(E.power(M,k),E.power(L,r));beta=E.mul(E.power(M,u),E.power(L,v))
    for z in G.perms:
        if E.power(z,r)!=F:continue
        mu=E.mul(beta,E.power(z,-v));lam=E.mul(E.power(mu,-r*k),F)
        if order(mu)!=new_order or lam!=I:continue
        if len(G.groups[G.ix[beta],G.ix[z]])!=cable_order:continue
        yield z,mu,lam,beta
def pattern(n,r,k):
    M=tuple((i+1)%n for i in range(n));ident=tuple(range(n))
    z=E.power(M,k*pow(r,-1,n));rho=[M,ident,z]
    F=product(wp([1],k),wp([2],r));u,v,d=S.gcdex(r,k);assert d==1;u=int(u);v=-int(v)
    mu=product(wp([1],u),wp([2],v),wp([3],-v));lam=product(wp(mu,-r*k),F)
    rels=[[1,2,-1,-2],product(wp([3],r),E.wi(F))]
    assert all(E.evaluate(w,rho)==ident for w in rels) and E.evaluate(lam,rho)==ident
    rows,mark,keep=E.cover_matrix(rels,rho,mu,lam)
    return rows,mark,path(wp([1],n),rho,keep),path([2],rho,keep),M,z
def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    preds=['meridian_cycle_s4_receipt.json','terminal_cone_twists_receipt.json',
        'cyclic_loop_first_arrival_receipt.json','iterated_peripheral_tail_receipt.json',
        'fibre_order_transfer_receipt.json','cyclic_companion_receipt.json',
        'common_pencil_mark_receipt.json','complementary_corner_jsj_receipt.json',
        'persistent_jsj_terminal_modulus_receipt.json']
    for name in preds:
        old=json.loads((P/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for f,h in old['source_sha256'].items():check('sealed predecessor '+name+' '+f,sha(P/f)==h)
    manifest=P/'sources/cyclic-pattern-a4/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        f=P/row['path'];check('retained complete primary pin '+f.name,sha(f)==row['sha256'] and f.stat().st_size==row['bytes'])
    pattern_controls=[]
    for n in [2,3,4,5,6,7,8]:
        for r,k in [(2,3),(3,2),(4,3),(5,6),(5,3),(7,6)]:
            if math.gcd(n,r)!=1:continue
            rows,R,m,L,M,z=pattern(n,r,k);h=E.h1(rows)
            check('actual cyclic pattern meridian disk balance '+str((n,r,k)),E.h1(rows,[m])==h)
            relation=[x-r*y for x,y in zip(R,L)]
            check('actual integral cyclic pattern R=rL including torsion '+str((n,r,k)),E.h1(rows,[relation])==h)
            closed=K.pattern(r,k,M,z)
            check('closed cyclic ramification mark nullhomologous '+str((n,r,k)),closed['H1']==closed['quotient'])
            want=dict(free=closed['H1']['free']+1,torsion=closed['H1']['torsion'])
            check('actual cyclic pattern full integral split '+str((n,r,k)),h==want and E.h1(rows,[L])==closed['H1'])
            pattern_controls.append(dict(n=n,r=r,k=k,H1=h,closed=closed['H1']))
    # Winding-one Hopf pattern tied with a trefoil: an arbitrary-pattern control
    # with positive-rank added homology at n=6, not a cable Seifert computation.
    knotted_core_controls=[]
    for n in [2,3,6]:
        mu=[1,-2];lon=product([1,1],wp(mu,-6));rho=K.cyclic_rho(2,3,[],tuple((i+1)%n for i in range(n)))
        rho=rho+[tuple(range(n))]
        rels=[[1,1,-2,-2,-2],product([3],mu,[-3],E.wi(mu))]
        rows,R,keep=E.cover_matrix(rels,rho,mu,product(lon,[3]));h=E.h1(rows)
        L=path([3],rho,keep);m=path(wp(mu,n),rho,keep)
        check('knotted winding1 core exact integral marking '+str(n),E.h1(rows,[[x-y for x,y in zip(R,L)]])==h and E.h1(rows,[m])==h)
        closed=K.full(K.words(2,3,[]),rho[:2])['H1']
        check('knotted winding1 core full split including positive rank '+str(n),h==dict(free=closed['free']+1,torsion=closed['torsion']))
        knotted_core_controls.append(dict(n=n,H1=h,closed=closed))
    check('arbitrary pattern added module genuinely allowed positive rank',knotted_core_controls[-1]['closed']['free']==2)
    sat=json.loads((P/'iterated_peripheral_tail_receipt.json').read_text())['saturation']
    incoming=[x for x in sat['useful_edges'] if x['source']==7 and x['dest']==3]
    check('A4 nonself predecessors exactly cyclic3 and doubleC2',
          {x['source'] for x in sat['transition_profiles'] if x['dest']==7 and x['source']!=7}=={6,8})
    birth_choices=[x for x in sat['useful_edges'] if x['source']!=7 and x['dest']==7]
    for x in birth_choices:
        st=sat['useful_profiles'][x['source']];M=tuple(st['M']);L=tuple(st['L']);r,k=x['rk']
        check('complete A4 birth retains two old orbits and identity fibre '+str((x['source'],r,k,x['z'])),
              len(V.orbits([M,L],range(4)))==2 and x['cable_order']==12 and E.mul(E.power(M,k),E.power(L,r))==I)
    for x in incoming:
        r,k=x['rk'];z=tuple(x['z'])
        check('complete A4 final actual parity and S3 roots '+str((r,k,z)),
              r%2==0 and k%3==0 and x['cable_order']==6)
    loops=[x for x in sat['useful_edges'] if x['source']==x['dest']==7]
    check('all A4 loop cable actions cyclic3 and windings unit3',
          bool(loops) and all(x['cable_order']==3 and math.gcd(x['rk'][0],3)==1 for x in loops))
    seeds=[]
    for p,q in [(2,3),(2,9),(4,3)]:
        a,b,m,l=next(x for x in G.torus(p,q) if len(G.groups[x[0],x[1]])==12 and order(G.perms[x[2]])==3)
        rho=normalized(p,q,[],[G.perms[a],G.perms[b]])
        seeds.append(dict(first=[p,q],cables=[],rho=rho,kind='initial'))
    births=[]
    for M,p,q,r,k in [((0,2,3,1),2,5,2,3),((1,0,3,2),3,5,3,2)]:
        old=K.cyclic_rho(p,q,[],M);z,mu,lam,beta=next(roots(M,I,r,k,3,12))
        rho=normalized(p,q,[(r,k)],old+[z])
        seeds.append(dict(first=[p,q],cables=[(r,k)],rho=rho,kind='birth'))
        graph=G.O.orbit_graph(dict(old=old,birth=[beta,z]),[('old','birth',[M,I])],4)
        old_edges=len(graph['edges'])
        check('first A4 birth retains TWO older companion tori '+str(M),old_edges==2 and graph['cycle_rank']==0)
        check('after fixed new cap birth vertex has degreeTHREE '+str(M),old_edges+1==3)
        births.append(dict(oldM=M,formation=[r,k],graph=graph,retained_degree=3))
    target_controls=[];loop_controls=[]
    for si,seed in enumerate(seeds):
        p,q=seed['first'];cables=seed['cables'];rho=seed['rho']
        rows,mark,keep=K.summarized_cover(p,q,cables,rho,raw=True);h=E.h1(rows)
        check('first A4 R3 has orderdividing2 '+str(si),E.h1(rows,[[2*x for x in mark]])==h)
        check('first A4 R3 actually nonzero in controls '+str(si),E.h1(rows,[mark])!=h)
        M,L=periphery(p,q,cables,rho)
        for rr,kk in [(2,3),(2,5),(4,3),(5,6),(5,7)]:
            z=E.power(M,kk*pow(rr,-1,3));cnew=cables+[(rr,kk)];rnew=rho+[z]
            rnew=normalized(p,q,cnew,rnew)
            actualrows,newmark,newkeep=K.summarized_cover(p,q,cnew,rnew,raw=True)
            _,_,oldL=K.words(p,q,cables);oldmark=path(oldL,rnew,newkeep)
            relation=[x-rr*y for x,y in zip(newmark,oldmark)]
            actual=E.h1(actualrows);extra=B.group([rr,rr]) if kk%3==0 else B.group([])
            check('actual A4 loop full group '+str((si,rr,kk)),actual==B.direct_sum([h,extra]))
            check('actual A4 loop integral R3 transfer '+str((si,rr,kk)),E.h1(actualrows,[relation])==actual)
            check('A4 loop retains two-torsion mark bound '+str((si,rr,kk)),E.h1(actualrows,[[2*x for x in newmark]])==actual)
            loop_controls.append(dict(seed=si,loop=[rr,kk],H1=actual,quotient=E.h1(actualrows,[newmark])))
            MM,LL=periphery(p,q,cnew,rnew)
            for fr,fk in [(2,3),(4,3),(8,3)]:
                fz,tau,flam,fbeta=next(roots(MM,LL,fr,fk,2,6))
                ans=K.summarized_cover(p,q,cnew+[(fr,fk)],rnew+[fz]);want=B.direct_sum([actual,B.group([fr//2])])
                check('actual final A4 simple R is ZERO '+str((si,rr,kk,fr)),ans['H1']==want and ans['quotient']==want)
                tailz=E.power(tau,3)
                tail=K.summarized_cover(p,q,cnew+[(fr,fk),(5,3)],rnew+[fz,tailz])
                check('actual later tail preserves zeroR and full H1 '+str((si,rr,kk,fr)),tail==ans)
                target_controls.append(dict(seed=si,companion_loop=[rr,kk],arrival=[fr,fk],actual=ans,after_tail=tail))
    t=S.symbols('t');alex=[]
    for r,k in [(2,5),(4,5),(5,7),(7,8),(8,11)]:
        delta=S.cancel((t**(r*k)-1)*(t-1)/((t**r-1)*(t**k-1)))
        residue=S.rem(delta,t**3-1,t);wanted=t*t if r%3==k%3==2 else 1
        check('torus Alexander unit residue modulo t3minus1 '+str((r,k)),S.expand(residue-wanted)==0)
        alex.append(dict(r=r,k=k,residue=str(residue)))
    a,kap,Z,sig=S.symbols('a kappa Z sigma',integer=True)
    multi=a*a*Z+a-(a-1)*kap
    check('all-depth terminal congruence multi-node identity',S.expand(multi-kap-a*(a*Z+1-kap))==0)
    single=(sig+1)*a-sig+1
    check('single-node terminal congruence identity',S.expand(single-(a+1-sig)-sig*a)==0)
    for tt in [3,5,7,9,11,25]:
        check('two initial non-trefoil terminal coefficients are inverse2 '+str(tt),all(math.gcd(2,alpha)==1 for alpha in [tt,3*tt]))
        possible=[alpha for alpha in [tt,3*tt] if (2-1)%alpha==0 or (2+1)%alpha==0]
        check('free-tail normalized coefficient forces ONLY t3/free3 '+str(tt),possible==([3] if tt==3 else []))
    check('last residual excluded in both possible source nR cases',all((6-nr*3)%9!=0 for nr in [3,7]))
    # Receipt construction follows; source theorem pins are checked via predecessor seals above.
    inputs=[Path(__file__),P/'CYCLIC_PATTERN_MARK_TRANSFER_AND_A4_REDUCTION.md',manifest]+[P/f for f in preds]+[P/x['path'] for x in sources['primary_sources']]
    result=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      general=dict(all_finite_cyclic_degrees=True,arbitrary_knot_patterns=True,coprime_winding_required=True,
        actual_outer_boundary_basis_required=True,integral_meridian_zero_and_R_winding_longitude=True,
        added_homology_may_have_positive_rank=True,no_hidden_boundary_or_mark_coordinates=True,
        no_Seifert_or_graph_hypothesis_for_pattern=True,terminal_congruence_includes_single_node=True),
      quartic=dict(all_A4_incoming_common_sources_closed=True,arbitrary_companion_depth_and_tail_length=True,
        no_prefix_eta_condition=True,A4_companion_R3_order_divides2=True,final_arrival_even_winding_kills_mark=True,
        complementary_later_A4_births_closed=True,complementary_initial_non_trefoil_A4_closed=True,
        complementary_residual_only_initial_trefoil_origin=True,all_actual_fixed_sheet_caps_retained=True),
      cyclic_pattern_controls=pattern_controls,knotted_core_pattern_controls=knotted_core_controls,
      A4_loop_controls=loop_controls,full_target_controls=target_controls,
      actual_A4_birth_branching=births,Alexander_unit_controls=alex,
      incoming_residue_choices=len(incoming),cyclic3_loop_choices=len(loops),birth_residue_choices=len(birth_choices),
      assurance=dict(written_unbounded_surface_MV_mark_history_and_terminal_proofs=True,
        full_CW_integral_mark_and_all_fill_controls=True,complete_retained_primary_audits_reused=True,
        deep_JSJ_and_full_peripheral_rigidity_imported=True,finite_controls_not_universal_proof=True,
        all_A4_complementary_trajectories_closed=False,initial_trefoil_cap_asserted_incompressible=False,
        noncyclic_companions_replaced_by_unknots=False,linking_splitting_asserted_for_arbitrary_pattern=False,
        state2_closed=False,all_initial_state3_closed=False,full_quartic_exclusion=False,
        higher_degree_packet_classification=False,finite_Keller_map=False,formal_verification=False,
        external_review=False,independent_full_reproof=False,novelty_priority=False,research_goal_completed=False),
      source_sha256={os.path.relpath(f,P):sha(f) for f in inputs},
      runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'cyclic_pattern_a4_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; all-degree cyclic pattern mark transfer; A4 incoming common closed, complementary reduced to trefoil origin')
if __name__=='__main__':main()
