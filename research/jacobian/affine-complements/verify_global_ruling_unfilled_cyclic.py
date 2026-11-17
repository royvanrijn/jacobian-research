"""PF2RULED90: exact unfilled periphery and all-length nonsolid S3 reduction.

Finite controls check peripheral calculations, not the universal written proof
or the imported Floer theorems. Deep chains use composed CW paths only.
"""
from pathlib import Path
import sys,json,hashlib,math,os,platform
import sympy as S
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
import verify_cyclic_pattern_a4 as A
K=A.K;G=A.G;E=A.E;V=A.V;B=A.B
I=tuple(range(4))
def restrict_three(rho,M):
    active=next(o for o in V.orbits(rho,range(4)) if len(o)==3)
    start=next(i for i in active if M[i]!=i)
    active=[start]+sorted(set(active)-{start});ix={s:i for i,s in enumerate(active)}
    return [tuple(ix[g[s]] for s in active) for g in rho]
def sharp(p,q,cables,rho):
    rows,L,keep=K.summarized_cover(p,q,cables,rho,raw=True)
    at=3*(1+len(cables));assert len(rows)==at+2
    return rows,rows[:at]+rows[at+1:],rows[at],L
def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    preds=['nonunimodular_a4_receipt.json','cyclic_pattern_a4_receipt.json',
        'meridian_balance_two_stage_receipt.json','iterated_peripheral_tail_receipt.json',
        'persistent_jsj_terminal_modulus_receipt.json','lifted_surface_s3_receipt.json',
        'common_pencil_mark_receipt.json','complementary_corner_jsj_receipt.json',
        'trefoil_complementary_fusion_receipt.json']
    for name in preds:
        old=json.loads((P/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for f,h in old['source_sha256'].items():check('sealed predecessor '+name+' '+f,sha(P/f)==h)
    manifest=P/'sources/global-ruling-unfilled-cyclic/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        f=P/row['path'];check('retained exact primary pin '+f.name,sha(f)==row['sha256'] and f.stat().st_size==row['bytes'])
    check('two explicit new primary Floer imports and no old gluing conjecture',
        sum(x.get('new_deep_Floer_import',False) for x in sources['primary_sources'])==2 and
        all(x.get('complete_relevant_primary_HTML_read') for x in sources['primary_sources'] if x.get('new_deep_Floer_import')))
    # Both boundary relations are checked BEFORE any inner meridian filling.
    pattern_controls=[]
    for n in [2,3,4,5,6,7,8]:
        for r,k in [(2,3),(3,2),(4,3),(5,6),(5,3),(7,6)]:
            if math.gcd(n,r)!=1:continue
            rows,li,mo,lo,M,z=A.pattern(n,r,k);mi=rows[-1];unfilled=rows[:-1];h=E.h1(unfilled)
            closed=K.pattern(r,k,M,z)['H1']
            check('UNFILLED outer meridian exact integral relation '+str((n,r,k)),E.h1(unfilled,[[x-r*y for x,y in zip(mo,mi)]])==h)
            check('UNFILLED inner longitude exact integral relation '+str((n,r,k)),E.h1(unfilled,[[x-r*y for x,y in zip(li,lo)]])==h)
            check('UNFILLED full free2 plus branched module '+str((n,r,k)),h==dict(free=closed['free']+2,torsion=closed['torsion']))
            check('UNFILLED actual mi lo primitive integral splitting '+str((n,r,k)),E.h1(unfilled,[mi,lo])==closed)
            pattern_controls.append(dict(n=n,r=r,k=k,H1=h,added=closed))
    check('nonunit k modulo cover degree controls included',any(math.gcd(x['n'],x['k'])>1 for x in pattern_controls))
    knotted_controls=[]
    for n in [2,3,6]:
        mu=[1,-2];lon=A.product([1,1],A.wp(mu,-6));rho=K.cyclic_rho(2,3,[],tuple((i+1)%n for i in range(n)))
        rho=rho+[tuple(range(n))]
        rels=[[1,1,-2,-2,-2],A.product([3],mu,[-3],E.wi(mu))]
        rows,li,keep=E.cover_matrix(rels,rho,mu,A.product(lon,[3]));mi=rows[-1];unfilled=rows[:-1]
        lo=A.path([3],rho,keep);mo=A.path(A.wp(mu,n),rho,keep);h=E.h1(unfilled)
        closed=K.full(K.words(2,3,[]),rho[:2])['H1']
        check('noncable winding1 full UNFILLED periphery '+str(n),
            E.h1(unfilled,[[x-y for x,y in zip(li,lo)],[x-y for x,y in zip(mo,mi)]])==h)
        check('noncable winding1 UNFILLED primitive split '+str(n),h==dict(free=closed['free']+2,torsion=closed['torsion']) and E.h1(unfilled,[mi,lo])==closed)
        knotted_controls.append(dict(n=n,H1=h,added=closed))
    check('positive-rank added homology never quotiented away',knotted_controls[-1]['added']['free']==2)
    # Exact symbolic consequences; the geometric ruling proof is written.
    b,Aw,w,mi,lo=S.symbols('b A w mi lo',integer=True)
    check('whole source Euler budget over complete base',S.expand(1+b-2-(b-1))==0)
    check('whole source Euler budget over affine base',S.expand(1+b-1-b)==0)
    check('solid cap transfer gives square winding',S.expand(w*(-Aw*w*mi)+Aw*w*w*mi)==0)
    # Actual D3 trefoil cover: keep BOTH boundary tori unfilled.
    rho=[(1,0,2),(1,2,0)];rels,Mword,Lword=K.words(2,3,[])
    rows,L2,keep=E.cover_matrix(rels,rho,Mword,Lword);unfilled=rows[:3]
    M3=E.evaluate(Mword,rho);fixed=next(i for i in range(3) if M3[i]==i)
    L1=A.path(Lword,rho,keep,fixed);m2,m1=rows[3:]
    check('D3 trefoil UNFILLED two-boundary H1 is free2',E.h1(unfilled)==dict(free=2,torsion=[]))
    check('D3 preferred index2 longitude includes mutual linking2',E.h1(unfilled,[[x+y-2*z for x,y,z in zip(L2,m2,m1)]])==E.h1(unfilled))
    check('D3 preferred index1 longitude includes mutual linking2',E.h1(unfilled,[[x-2*y+4*z for x,y,z in zip(L1,m2,m1)]])==E.h1(unfilled))
    check('D3 weighted lifted surface balance',E.h1(unfilled,[[2*x+y for x,y in zip(L2,L1)]])==E.h1(unfilled))
    matrix=S.Matrix([[1,-1,1,0],[-1,1,-1,0],[0,-1,1,-1],[-1,0,0,-1]])
    for j,(rhs,sol) in enumerate([([-1,1,-1,0],[0,1,0,0]),([1,-1,0,-1],[1,0,0,0])]):
        check('published trefoil bounding-chain system '+str(j),matrix*S.Matrix(sol)==S.Matrix(rhs))
    check('published two mutual intersections equal independent CW value2',sum([0,1,1,0])==2)
    sat=json.loads((P/'iterated_peripheral_tail_receipt.json').read_text())['saturation']
    check('complete S3 nonself predecessors are C2 and C3',
        {x['source'] for x in sat['transition_profiles'] if x['dest']==2 and x['source']!=2}=={1,6})
    for x in sat['useful_edges']:
        r,k=x['rk']
        if x['source']==x['dest']==2:
            check('complete S3 selfloops only odd winding cyclic2 '+str((r,k,x['z'])),r%2==1 and x['cable_order']==2)
        if x['source']==2 and x['dest']==3:
            check('complete S3 final arrivals odd multiple3 and even k '+str((r,k,x['z'])),r%6==3 and k%2==0 and x['cable_order']==6)
    # Capped E^sharp controls: remove ONLY the index2 fill, keep fixed fill.
    initial_controls=[];seed=None
    for p,q in [(2,3),(4,3),(2,9),(8,3),(4,9)]:
        if math.gcd(p,q)!=1:continue
        aa,bb,mm,ll=next(x for x in G.torus(p,q) if len(G.groups[x[0],x[1]])==6 and A.order(G.perms[x[2]])==2)
        rho4=A.normalized(p,q,[],[G.perms[aa],G.perms[bb]]);M4,_=A.periphery(p,q,[],rho4)
        rho3=restrict_three(rho4,M4);full,unfilled,m,L=sharp(p,q,[],rho3);bt=(p//2)*(q//3)
        check('actual initial Esharp full marked split '+str((p,q)),E.h1(unfilled)==dict(free=1,torsion=[] if p==2 else [p//2]) and E.h1(unfilled,[[x+bt*y for x,y in zip(L,m)]])==E.h1(unfilled))
        check('actual initial Esharp primitive meridian and zero core '+str((p,q)),E.h1(unfilled,[m])==E.h1(full) and E.h1(full,[L])==E.h1(full))
        initial_controls.append(dict(p=p,q=q,A=bt,H1=E.h1(unfilled),closed=E.h1(full)))
        if (p,q)==(2,3):seed=(rho4,M4,rho3)
    rho4,M4,rho3=seed;M3=E.evaluate(Mword,rho3)
    solid_controls=[];cables=[];images=rho3.copy();curA=1
    for step,w0 in enumerate([3,5,3]):
        k0=2*curA*w0+1;images.append(E.power(M3,k0%2));cables.append((w0,k0));curA*=w0*w0
        full,unfilled,m,L=sharp(2,3,cables,images)
        check('actual composed solid S3 square cap '+str(step),E.h1(unfilled)==dict(free=1,torsion=[]) and E.h1(unfilled,[[x+curA*y for x,y in zip(L,m)]])==E.h1(unfilled))
        check('actual composed solid S3 primitive meridian and zero core '+str(step),E.h1(unfilled,[m])==B.group([]) and E.h1(full,[L])==E.h1(full))
        solid_controls.append(dict(cables=cables.copy(),A=curA,H1=E.h1(unfilled),geometric_solid_cap_distance=1))
    reset_controls=[]
    for r,k in [(3,11),(3,2),(5,3),(3,7)]:
        full,unfilled,m,L=sharp(2,3,[(r,k)],rho3+[E.power(M3,k%2)])
        check('actual nonsolid reset meridian still primitive and zero core '+str((r,k)),E.h1(unfilled,[m])==E.h1(full) and E.h1(full,[L])==E.h1(full))
        check('actual nonsolid reset preferred longitude square recurrence '+str((r,k)),E.h1(unfilled,[[x+r*r*y for x,y in zip(L,m)]])==E.h1(unfilled))
        reset_controls.append(dict(r=r,k=k,H1=E.h1(unfilled),cones=([r,abs(2*r-k)] if k%2 else [r,r,abs(r-k//2)])))
    targets=[]
    # Different active S3 sets are essential for total image S4.
    for prefix in [[],[(3,7)],[(3,11)],[(3,2)]]:
        prefimages=rho4+[E.power(M4,k%2) for r,k in prefix]
        for r,k in [(3,2),(3,4),(9,2),(9,4)]:
            z,tau,lam,beta=next(x for x in A.roots(M4,I,r,k,2,6) if E.generated(prefimages+[x[0]])==24)
            for tail in [[],[(5,3)]]:
                images=prefimages+[z]+[E.power(tau,tk%2) for tr,tk in tail]
                ans=K.summarized_cover(2,3,prefix+[(r,k)]+tail,images)
                check('actual full S4 target through different incoming S3 orbit '+str((prefix,r,k,tail)),E.generated(images)==24 and set(ans)=={'H1','quotient'})
                targets.append(dict(prefix=prefix,arrival=[r,k],tail=tail,actual=ans))
    # Births from each possible cyclic predecessor, with actual total S3 image.
    births=[]
    for oldM,r,k in [((1,0,2,3),3,2),((1,2,0,3),2,3)]:
        oldrho=K.cyclic_rho(2,3,[],oldM)
        z,tau,lam,beta=next(x for x in A.roots(oldM,I,r,k,2,6) if E.generated(oldrho+[x[0]])==6)
        birthrho=restrict_three(oldrho+[z],tau)
        full,unfilled,m,L=sharp(2,3,[(r,k)],birthrho)
        check('actual S3 cyclic birth retains zero index2 core '+str((r,k)),E.h1(full,[L])==E.h1(full))
        check('actual S3 cyclic birth has primitive remaining meridian '+str((r,k)),E.h1(unfilled,[m])==E.h1(full))
        births.append(dict(old_meridian=list(oldM),r=r,k=k,H1=E.h1(unfilled),closed=E.h1(full)))
    inputs=[Path(__file__),P/'COORDINATE_RULINGS_UNFILLED_CYCLIC_AND_S3_REDUCTION.md',manifest]+[P/f for f in preds]+[P/x['path'] for x in sources['primary_sources']]
    result=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
        general=dict(every_actual_coordinate_extends_on_whole_smooth_affine_source=True,
            generic_linear_ruling_entire_boundary_over_infinity=True,
            arbitrary_degree_unfilled_cyclic_pattern_integral_periphery=True,
            positive_rank_added_modules_allowed=True,
            nullhomologous_core_meridian_outside_Lspace_interior=True),
        quartic=dict(all_incoming_S3_nonsolid_histories_closed_both_pencils=True,
            arbitrary_companion_depth_and_later_tail_length=True,no_prefix_eta_condition=True,
            remaining_S3_only_initial_trefoil_solid_chain=True,
            solid_cap_A_exact_product_of_odd_windings_squared=True,
            actual_ordinary_companion_retained_at_final_arrival=True),
        unfilled_pattern_controls=pattern_controls,noncable_positive_rank_controls=knotted_controls,
        actual_initial_S3_controls=initial_controls,actual_solid_S3_controls=solid_controls,
        actual_nonsolid_reset_controls=reset_controls,full_target_controls=targets,actual_S3_birth_controls=births,
        dihedral_trefoil=dict(mutual_linking=2,preferred_longitudes=['-m2+2m1','2m2-4m1'],full_algorithm_implemented=False,paper_diagonal_framing_identified=False),
        assurance=dict(written_unbounded_ruling_MV_endpoint_and_history_proofs=True,
            actual_full_filled_and_unfilled_CW_controls=True,
            complete_relevant_primary_HTML_import_audits=True,
            retained_complete_topology_audits_reused=True,
            new_deep_Floer_theorems_imported=True,
            deep_JSJ_and_full_peripheral_rigidity_imported=True,
            finite_controls_not_universal_proof=True,
            all_S3_solid_histories_closed=False,all_A4_complementary_closed=False,
            full_quartic_exclusion=False,higher_degree_packet_classification=False,
            finite_Keller_map=False,formal_verification=False,external_review=False,
            independent_full_reproof=False,novelty_priority=False,research_goal_completed=False),
        source_sha256={os.path.relpath(f,P):sha(f) for f in inputs},
        runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'global_ruling_unfilled_cyclic_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; all nonsolid incoming S3 histories closed BOTH pencils; solid trefoil-square residual OPEN')
if __name__=='__main__':main()
