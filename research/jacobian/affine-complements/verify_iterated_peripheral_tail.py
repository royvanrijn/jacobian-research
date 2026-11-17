"""PF2RULED82: complete finite saturation and actual iterated-tail controls."""
from pathlib import Path
import sys,math,itertools,collections,json,time
P=Path(__file__).resolve().parent; R=P.parents[2]
sys.path.insert(0,str(P))
import verify_general_two_stage as G
V=G.V
I=G.I; pp=G.perms; mul=G.mul; power=G.power
pairs=[(r,k) for r,k in itertools.product(range(12),repeat=2) if math.gcd(math.gcd(r,k),12)==1]
bez={(r,k):next((u,v) for u,v in itertools.product(range(12),repeat=2) if (r*u-k*v)%12==1) for r,k in pairs}
partcache={};groupcache={};topcache={};conjcache={};canoncache={}
def group(gs):
    gs=tuple(sorted(set(gs)))
    if gs not in groupcache:groupcache[gs]=tuple(sorted(G.generated(gs)))
    return groupcache[gs]
def part(gs):
    gs=tuple(sorted(set(gs)))
    if gs not in partcache:partcache[gs]=V.orbits([pp[x] for x in gs],range(4))
    return partcache[gs]
def canonical(h,m,l):
    key=(h,m,l)
    if key not in canoncache:
        variants=[]
        for c in range(24):
            def cj(x):return mul[mul[G.inverse[c]][x]][c]
            variants.append((tuple(sorted(cj(x) for x in h)),cj(m),cj(l)))
        ans=min(variants)
        for val in variants:canoncache[val]=ans
    return canoncache[key]
def cable_top(h,m,l,beta,z,mu,F):
    key=(h,m,l,beta,z,mu,F)
    if key not in topcache:
        ch=group([beta,z]);nh=group(h+(z,))
        cyc=len(part([m,l]))-len(part(h))-len(part(ch))+len(part(nh))
        genera=[]
        for block in part(ch):
            (zz,bb,tt),degree=V.fibre_quotient([pp[z],pp[beta],pp[mu]],pp[F],block)
            genera.append(V.genus_from_base([zz],[bb,tt]))
        topcache[key]=(cyc,sum(genera),ch,nh)
    return topcache[key]
initial={};initial_rows=0
for p,q in pairs:
    for a,b,m,l in G.torus(p,q):
        H=power(a,p);ok=True
        for block in part([a,b]):
            (aa,bb,mm),degree=V.fibre_quotient([pp[a],pp[b],pp[m]],pp[H],block)
            if V.genus_from_base([aa,bb],[mm]):ok=False;break
        if ok:
            key=canonical(group([a,b]),m,l)
            initial_rows+=1
            initial.setdefault(key,dict(first=[p,q,a,b],steps=[]))
states=dict(initial);todo=collections.deque(initial);edges={};reject=collections.Counter();started=time.time()
while todo:
    state=todo.popleft();h,m,l=state;out=[]
    for r,k in pairs:
        u,v=bez[r,k];F=mul[power(m,k)][power(l,r)];beta=mul[power(m,u)][power(l,v)]
        for z in G.roots[r][F]:
            mu=mul[beta][power(z,-v)];lam=mul[power(mu,-r*k)][F]
            cyc,gen,ch,nh=cable_top(h,m,l,beta,z,mu,F)
            if cyc or gen:reject['cycle' if cyc else 'genus']+=1;continue
            dest=canonical(nh,mu,lam)
            out.append((dest,(r,k,z),len(ch)))
            if dest not in states:
                states[dest]=dict(first=states[state]['first'],steps=states[state]['steps']+[(r,k,z)])
                todo.append(dest)
    edges[state]=out
    if len(edges)%20==0:print('saturated',len(edges),'of',len(states),'states',round(time.time()-started,2),flush=True)
final=[s for s in states if len(s[0])==24 and s[1] in G.trans and s[2]==I]
reverse=collections.defaultdict(set)
for s,rows in edges.items():
    for t,_,_ in rows:reverse[t].add(s)
useful=set(final);q=collections.deque(final)
while q:
    t=q.popleft()
    for s in reverse[t]:
        if s not in useful:useful.add(s);q.append(s)
report=dict(initial_states=len(initial),initial_rows=initial_rows,all_states=len(states),edges=sum(map(len,edges.values())),
    useful_states=len(useful),final_states=len(final),reject=dict(reject),seconds=time.time()-started,
    useful_profiles=[dict(order=len(h),M=pp[m],L=pp[l],shortest_depth=1+len(states[h,m,l]['steps'])) for h,m,l in sorted(useful)],
    useful_edges=[dict(source=sorted(useful).index(s),dest=sorted(useful).index(t),rk=[r,k],z=pp[z],cable_order=ch)
        for s in sorted(useful) for t,(r,k,z),ch in edges[s] if t in useful],
    final_profiles=[dict(order=len(h),M=pp[m],L=pp[l]) for h,m,l in sorted(final)])

import hashlib,os,platform
import sympy as S
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
preds=['meridian_balance_two_stage_receipt.json','general_two_stage_receipt.json',
       'lifted_surface_s3_receipt.json','fibre_order_transfer_receipt.json',
       'orbit_tree_receipt.json','seifert_base_genus_receipt.json',
       'rational_end_orderability_receipt.json','fold_schur_s4_receipt.json',
       'common_pencil_mark_receipt.json']
for name in preds:
    old=json.loads((P/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
    for path,digest in old['source_sha256'].items():
        check('sealed predecessor '+name+' '+path,sha(P/path)==digest)
manifest=P/'sources/iterated-peripheral-tail/manifest.json'
sources=json.loads(manifest.read_text())
for row in sources['primary_sources']:
    path=P/row['path']
    check('retained primary byte pin '+path.name,sha(path)==row['sha256'] and path.stat().st_size==row['bytes'])
check('complete mod12 parameter pairs',len(pairs)==96)
check('complete labelled planar torus starting rows',initial_rows==3648)
check('ten initial states saturate without a new state',len(initial)==len(states)==10 and set(initial)==set(states))
check('all transitions and rejections exhausted',report['edges']==1596 and reject==dict(cycle=184,genus=60))
profiles=collections.Counter((sorted(useful).index(s),sorted(useful).index(t),ch)
    for s in useful for t,_,ch in edges[s] if t in useful)
check('all29 transition profiles retained',len(profiles)==29 and len(useful)==10)
check('unique final full simple S4 state',len(final)==1 and len(final[0][0])==24)
check('nontrivial longitude state retained',sum(l!=I for h,m,l in states)==1 and
      all(l==power(m,2) and len(h)==24 for h,m,l in states if l!=I))
for i,state in enumerate(sorted(states)):
    h,m,l=state
    check('state boundary pair commutes '+str(i),mul[m][l]==mul[l][m] and m in h and l in h)
    check('state complete canonical transition closure '+str(i),all(t in states for t,_,_ in edges[state]))
tailstate=final[0];h,m,l=tailstate
check('absorbing final state has64 passing roots',len(edges[tailstate])==64)
for r,k in pairs:
    u,v=bez[r,k];F=mul[power(m,k)][power(l,r)];beta=mul[power(m,u)][power(l,v)]
    good=[]
    for z in G.roots[r][F]:
        mu=mul[beta][power(z,-v)];lam=mul[power(mu,-r*k)][F]
        cyc,gen,ch,nh=cable_top(h,m,l,beta,z,mu,F)
        if not cyc and not gen:good.append((z,mu,lam,ch))
    check('exact absorbing parity/root for residue '+str((r,k)),
          good==[(power(m,k),m,I,group([m]))] if r%2 else not good)
# Exact positive representatives exercise integer gcd beyond modular conditions.
for r,k in [(3,2),(3,5),(5,2),(5,3),(9,2),(9,5),(13,14),(15,14)]:
    check('actual unit transfer pair '+str((r,k)),math.gcd(r,k)==1 and
          (k%2 or math.gcd(k//2,r)==1))
    check('order4 half linking survives odd winding '+str(r),(r*r-1)%8==0)

E=G.E;B=G.B
# Independent agreement with every class of the retained two-stage census.
oldcells=json.loads((P/'general_two_stage_receipt.json').read_text())['census']['complete_cells']
oldcounts=collections.Counter()
for key,cell in oldcells.items():
    p,q,r,k=map(int,key.split(','))
    for row in cell['orbits']:
        a,b,z=[G.ix[tuple(x)] for x in row['rho']]
        mm,ll,beta,F,mu=[G.ix[tuple(row[x])] for x in ['M','L','beta','h','mu']]
        cyc,gen,ch,nh=cable_top(group([a,b]),mm,ll,beta,z,mu,F)
        oldgen=sum(x['genus'] for x in row['companion_bases'])+sum(x['genus'] for x in row['cable_bases'])
        check('whole graph agrees retained two-stage class '+str((key,a,b,z)),cyc==row['cycle_rank'])
        if not cyc:
            oldcounts['genus' if oldgen else 'accepted']+=1
            if not oldgen:
                src=canonical(group([a,b]),mm,ll)
                dest=canonical(nh,mu,I)
                check('retained planar two-stage path accepted '+str((key,a,b,z)),
                      src in initial and dest==tailstate and any(t==dest and rk[:2]==(r,k) for t,rk,_ in edges[src]))
        else:oldcounts['cycle']+=1
check('all2752 old classes agree forest and planar gates',oldcounts==dict(cycle=448,genus=128,accepted=2176))
for a in [2,3,5,7]:
    balance=[2,3,5,7,0];i=balance.index(a)
    caps=[[int(c==j) for c in range(5)] for j in range(4) if j!=i]
    check('nonunit many-boundary cap retains torsion '+str(a),E.h1([balance]+caps)==B.group([a],free=1))
    meridian=[int(c==i) for c in range(5)]
    check('balanced last filling retains infinite H1 '+str(a),E.h1([balance]+caps+[meridian])==B.group(free=1))
def reduce_word(w):
    stack=[]
    for x in w:
        if stack and stack[-1]==-x:stack.pop()
        else:stack.append(x)
    return stack
def product(*ws):return reduce_word([x for w in ws for x in w])
def wp(w,n):return reduce_word(E.wp(w,n))
def ordinary(p,q,cables,rho,reverse=False,mark_multiple=1):
    uu,vv,d=S.gcdex(q,p);assert d==1
    M=product(wp([1],int(uu)),wp([2],int(vv)))
    L=product(wp([1],p),wp(M,-p*q));rels=[product(wp([1],p),wp([2],-q))]
    for j,(r,k) in enumerate(cables,3):
        uu,vv,d=S.gcdex(r,k);assert d==1;u=int(uu);v=-int(vv)
        F=product(wp(M,k),wp(L,r));rels.append(product(wp([j],r),E.wi(F)))
        M=product(wp(M,u),wp(L,v),wp([j],-v));L=product(wp(M,-r*k),F)
    assert all(E.evaluate(w,rho)==pp[I] for w in rels) and E.evaluate(L,rho)==pp[I]
    assert E.evaluate(M,rho) in [pp[x] for x in G.trans]
    rows,mark,_=E.cover_matrix(rels,rho,M,L,reverse_tree=reverse)
    return dict(H1=E.h1(rows),quotient=E.h1(rows,[[mark_multiple*x for x in mark]]))

controls=[]
# Zero marks, order4 marks, nonzero S3 marks and free homology; both parities.
cases=[(2,3,4,3),(4,3,3,2),(2,9,9,8),(2,3,3,2),(2,3,9,2),(2,9,2,3)]
for p,q,r,k in cases:
    cell=G.cell(p,q,r,k)
    rows=[x for x in cell['orbits'] if not x['cycle_rank'] and
          not sum(b['genus'] for b in G.base_records(x,p)[0])]
    for ci,row in enumerate(rows):
        rho=[tuple(x) for x in row['rho']];seed=ordinary(p,q,[(r,k)],rho)
        for tail in [[(3,2)],[(3,5)],[(3,2),(5,3)]]:
            image=list(rho);cables=[(r,k)];added=[];prod=1
            for rr,kk in tail:
                # The boundary meridian permutation is unchanged by every odd tail.
                image.append(E.power(tuple(row['mu']),kk));cables.append((rr,kk));prod*=rr
                if kk%2==0:added.append(rr)
            actual=ordinary(p,q,cables,image)
            wantH=B.direct_sum([seed['H1'],B.group(added)])
            # Compare with the genuinely multiplied seed class, including
            # annihilation when winding and mark order are not coprime.
            scaled=ordinary(p,q,[(r,k)],rho,mark_multiple=prod)
            wantQ=B.direct_sum([scaled['quotient'],B.group(added)])
            check('full ordinary tail mark '+str((p,q,r,k,ci,tail)),actual['quotient']==wantQ)
            check('full ordinary tail H1 '+str((p,q,r,k,ci,tail)),actual['H1']==wantH)
            controls.append(dict(seed=[p,q,r,k],class_index=ci,tail=tail,seed_groups=seed,actual=actual))
            if ci==0 and len(tail)==2:
                check('independent full-cover spanning tree '+str((p,q,r,k)),ordinary(p,q,cables,image,True)==actual)

report['transition_profiles']=[dict(source=a,dest=b,cable_order=c,count=n) for (a,b,c),n in sorted(profiles.items())]
inputs=[Path(__file__),P/'ITERATED_PERIPHERAL_STATES_AND_ABSORBING_TAILS.md',manifest]+[P/n for n in preds]+[P/x['path'] for x in sources['primary_sources']]
receipt=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
    saturation=report,ordinary_full_cover_controls=controls,
    general=dict(all_finite_degrees=True,all_finite_cabling_lengths=True,
        exact_group_exponent_period=True,whole_component_graph_contraction=True,
        nonidentity_fibre_quotients_retained=True,multi_boundary_nonunit_balance_lemma=True,
        every_ordered_partial_filling_nonsolid_incompressible_required=True),
    quartic=dict(exactly_ten_states=True,exactly29_transition_profiles=True,
        final_state_absorbing=True,all_later_windings_odd=True,tail_image_C2=True,
        even_k_adds_Zr_odd_k_no_new_torsion=True,actual_R_mark_multiplied_by_winding=True,
        zero_mark_and_order4_half_linking_persist=True,positive_betti_persists=True,
        non_L_space_tail_requires_actual_seed_R_exterior_irreducible_incompressible_graph=True,
        seed_R_exterior_hypotheses_reproved_after_each_stage=True,
        arbitrary_tails_of_PF81_negative_even_seeds_excluded_for_common_source=True),
    assurance=dict(written_unbounded_state_period_and_tail_proofs=True,
        complete_finite_saturation_after_period_proof=True,retained_primary_page_audits_reused=True,
        actual_all_meridian_cycle_fillings_in_full_cover_controls=True,
        deep_topology_imported=True,finite_controls_not_universal_proof=True,
        graph_states_asserted_to_determine_integral_end=False,
        irreducibility_after_drilling_inferred_from_closed_irreducibility=False,
        complementary_tails_excluded=False,all_longer_target_sequences_excluded=False,
        higher_degree_packets_classified=False,full_quartic_exclusion=False,
        finite_Keller_map=False,formal_verification=False,external_review=False,
        independent_full_reproof=False,novelty_priority=False,research_goal_completed=False),
    source_sha256={os.path.relpath(path,P):sha(path) for path in inputs},
    runtime=dict(python=platform.python_version(),sympy=S.__version__))
(P/'iterated_peripheral_tail_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(f'PASS {len(checks)} controls; ten states at all lengths; actual odd-winding common-source tail transfer')
