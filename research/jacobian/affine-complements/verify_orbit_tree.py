"""PF2RULED62: exact orbit graphs and exhaustive exponent-period certificates."""
from pathlib import Path
import hashlib,importlib.util,itertools,json,os,platform
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('sealed_end61',HERE/'verify_marked_end_lattice.py')
G=importlib.util.module_from_spec(spec);spec.loader.exec_module(G)
def partition(gens,degree):
    blocks=[];seen=set()
    for s in range(degree):
        if s in seen:continue
        block={s};todo=[s]
        while todo:
            t=todo.pop()
            for g in gens:
                if g[t] not in block:block.add(g[t]);todo.append(g[t])
        blocks.append(sorted(block));seen|=block
    return blocks
def orbit_graph(vertex_gens,edges,degree):
    """Each edge is (left,right,edge generators) in the same based sheet action."""
    parts={v:partition(gs,degree) for v,gs in vertex_gens.items()}
    vertices=[(v,i) for v in parts for i in range(len(parts[v]))]
    lifted=[]
    for tag,(left,right,gens) in enumerate(edges):
        for block in partition(gens,degree):
            a=[i for i,p in enumerate(parts[left]) if set(block)<=set(p)]
            b=[i for i,p in enumerate(parts[right]) if set(block)<=set(p)]
            assert len(a)==len(b)==1,'edge orbits require actual consistent based inclusions'
            lifted.append(dict(base_edge=tag,letters=block,left=[left,a[0]],right=[right,b[0]]))
    parent={v:v for v in vertices}
    def root(v):
        while parent[v]!=v:v=parent[v]
        return v
    for e in lifted:
        a=tuple(e['left']);b=tuple(e['right']);parent[root(a)]=root(b)
    components=len({root(v) for v in vertices})
    return dict(vertex_orbits=parts,edges=lifted,vertices=len(vertices),edge_count=len(lifted),components=components,cycle_rank=len(lifted)-len(vertices)+components)
def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    preds=['sign_square_receipt.json','peripheral_zero_receipt.json','marked_end_lattice_receipt.json']
    for name in preds:
        r=json.loads((HERE/name).read_text());check('predecessor PASS '+name,r['outcome']=='PASS')
        for p,h in r['source_sha256'].items():check('sealed predecessor '+name+' '+p,sha(HERE/p)==h)
    manifest=HERE/'sources/orbit-tree/manifest.json';m=json.loads(manifest.read_text())
    for row in m['primary_sources']:
        p=HERE/row['path'];check('primary retained pin '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    perms=list(itertools.permutations(range(4)));I=tuple(range(4));trans=[p for p in perms if sum(p[i]!=i for i in range(4))==2]
    check('exact exponent12 of S4',all(G.power(p,12)==I for p in perms) and any(G.power(p,3)!=I for p in perms) and any(G.power(p,4)!=I for p in perms))
    def orbit_key(gs):return min(tuple(tuple(G.inv(c)[g[c[i]]] for i in range(4)) for g in gs) for c in perms)
    # Orbit graph controls, including a global cycle missed by all one-edge tests.
    swap=(1,0);identity=(0,1)
    path=orbit_graph(dict(A=[swap],B=[identity],C=[swap]),[('A','B',[identity]),('B','C',[identity])],2)
    check('multistage four-cycle control',path['components']==1 and path['cycle_rank']==1 and (path['vertices'],path['edge_count'])==(4,4))
    for side in ['A','C']:
        one=orbit_graph({side:[swap],'B':[identity]},[(side,'B',[identity])],2)
        check('each separate edge can pass '+side,one['components']==1 and one['cycle_rank']==0)
    # Degree2 through7 regular actions verify the generic abelian semiregularity control.
    for degree in range(2,8):
        shift=tuple((i+1)%degree for i in range(degree))
        check('regular cyclic action semiregular degree'+str(degree),all(len({len({G.power(shift,j*k)[s] for j in range(degree)}) for s in range(degree)})==1 for k in range(degree)))
    for a,b in itertools.product(perms,repeat=2):
        if G.mul(a,b)!=G.mul(b,a) or len(partition([a,b],4))!=1:continue
        # A commuting transitive pair generates a regular abelian subgroup.
        subgroup={G.mul(G.power(a,i),G.power(b,j)) for i in range(12) for j in range(12)}
        assert len(subgroup)==4 and all(p==I or all(p[s]!=s for s in range(4)) for p in subgroup)
    check('all transitive commuting pairs in S4 act regularly',True)
    # Every subgroup of A4 is generated incrementally; double transpositions cannot normally generate a transitive image.
    even=[p for p in perms if sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))%2==0]
    def closure(gs):
        out={I};todo=[I]
        while todo:
            a=todo.pop()
            for b in gs:
                c=G.mul(a,b)
                if c not in out:out.add(c);todo.append(c)
        return frozenset(out)
    subgroups={frozenset([I])};todo=list(subgroups)
    while todo:
        h=todo.pop()
        for p in even:
            k=closure(list(h)+[p])
            if k not in subgroups:subgroups.add(k);todo.append(k)
    transitive_even=[h for h in subgroups if len(partition(list(h),4))==1]
    check('transitive subgroups of A4 exactly V4 and A4',sorted(map(len,transitive_even))==[4,12])
    for h in transitive_even:
        for p in h:
            if p==I or G.power(p,2)!=I:continue
            normal=closure([G.mul(G.mul(g,p),G.inv(g)) for g in h])
            check('double transposition fails normal generation in order'+str(len(h)),normal!=h)
    expected={0:(0,0,0),1:(0,0,0),2:(72,3,0),3:(72,3,0),4:(0,0,0),5:(0,0,0),6:(0,0,0),7:(0,0,0),8:(168,7,3),9:(144,6,0),10:(0,0,0),11:(0,0,0)}
    retained61=json.loads((HERE/'marked_end_lattice_receipt.json').read_text());control61={}
    for w in ['3','8']:
        control61[w]={orbit_key([r['rho'][i] for i in [0,1,2,4]]):r for r in retained61['control'][w]['orbits']}
    certificates={}
    for residue in range(12):
        classes={};total=0
        for a,b in itertools.product(perms,repeat=2):
            if G.power(a,2)!=G.power(b,3):continue
            M=G.mul(a,G.inv(b));L=G.mul(G.power(a,2),G.power(M,-6))
            for beta in perms:
                if G.power(beta,residue)!=M or G.power(beta,-residue-1)!=L:continue
                for tau in trans:
                    if G.power(G.mul(tau,beta),residue)!=I or G.generated([a,b,beta,tau])!=24:continue
                    total+=1;key=orbit_key([a,b,beta,tau])
                    if key not in classes:
                        graph=orbit_graph(dict(companion=[a,b],cable=[beta,tau]),[('companion','cable',[M,L])],4)
                        assert graph['components']==1
                        classes[key]=dict(a=a,b=b,beta=beta,tau=tau,orbit_size=0,companion_order=G.generated([a,b]),cable_order=G.generated([beta,tau]),graph=graph)
                    classes[key]['orbit_size']+=1
        check('complete residue census '+str(residue),(total,len(classes),sum(r['graph']['cycle_rank']>0 for r in classes.values()))==expected[residue])
        for i,(key,row) in enumerate(classes.items()):
            check(f'residue{residue} orbit{i} size24',row['orbit_size']==24)
            if str(residue) in control61:
                check(f'residue{residue} orbit{i} agrees with full H1',row['graph']['cycle_rank']==control61[str(residue)][key]['H1']['free'])
            if row['graph']['cycle_rank']>0:
                check(f'residue{residue} orbit{i} global A4 S4 obstruction',residue==8 and row['companion_order']==12 and row['cable_order']==24 and row['graph']['cycle_rank']==1)
            elif row['graph']['vertices']==2:
                M=G.mul(row['a'],G.inv(row['b']))
                check(f'residue{residue} orbit{i} transitive pieces require 4-cycle',len(partition([M],4))==1)
        certificates[str(residue)]=dict(total_representations=total,conjugacy_classes=len(classes),tree_classes=sum(r['graph']['cycle_rank']==0 for r in classes.values()),orbits=list(classes.values()))
    # Exact CRT is modulo36: no finite-winding extrapolation.
    allowed=[]
    for w in range(36):
        peripheral=certificates[str(w%12)]['tree_classes']>0
        determinant_ok=((w+1)%9==0) if w%2==0 else ((3*w)%9==0)
        if peripheral and determinant_ok:allowed.append(w)
        check('CRT equivalent residue'+str(w),(peripheral and determinant_ok)==(w%6==3 or w%18==8))
    check('even positive Betti obstruction repeats at every w8mod12',certificates['8']['tree_classes']==4 and certificates['2']['tree_classes']==3)
    inputs=[Path(__file__),HERE/'MONODROMY_ORBIT_TREE_AND_UNBOUNDED_WINDING.md',manifest]+[HERE/n for n in preds]+[HERE/row['path'] for row in m['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
        graph_theorem=dict(all_finite_degrees=True,all_finite_torus_decomposition_trees=True,all_cabling_lengths=True,all_sheets_and_meridian_cycles=True,cover_need_not_be_regular=True,no_minimal_JSJ_hypothesis=True,free_graph_homology_quotient=True,smooth_affine_source_requires_tree=True,smoothness_explicit_in_higher_degrees=True,tree_sufficient_for_affine_cover=False),
        consequences=dict(transitive_pieces_force_transitive_abelian_peripheral=True,meridian_semiregularity_all_degrees=True,transitive_quartic_pieces_force_companion_4cycle=True,whole_graph_not_separate_edge_counts=True,only_dcycles_asserted_in_higher_degrees=False),
        periodic=dict(primitive_actual_binomial_trefoil_key_required=True,quartic_eta_zero_forces_parameter_w_plus_one=True,period=12,period_from_exact_group_exponent=True,all_windings_at_least_two=True,residue_certificates=certificates,determinant_square_gate=9,combined_period=36,allowed_mod36=allowed,allowed_odd='w=3 mod6',allowed_even='w=8 mod18',w8_mod12_cycle_classes=3,even_w8_mod36_tree_classes=4,even_w26_mod36_tree_classes=3,all_passing_infinity_types_realized=False),
        assurance=dict(written_unbounded_proofs=True,exhaustive_finite_group_certificate_after_period_reduction=True,finite_search_inferred_unbounded_theorem=False,symplectic_filling_imported=False,named_boundary_manifold_identification_used=False,full_quartic_exclusion=False,higher_degree_exclusion=False,finite_Keller_map=False,affine_monodromy_extension_supplied=False,formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,research_goal_completed=False),
        source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=G.S.__version__))
    (HERE/'orbit_tree_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; all-degree orbit tree and exact unbounded winding reduction')
if __name__=='__main__':main()
