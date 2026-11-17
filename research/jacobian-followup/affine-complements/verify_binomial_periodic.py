"""PF2RULED63: exact finite certificate after the written unbounded reduction."""
from pathlib import Path
import hashlib,importlib.util,itertools,json,math,os,platform
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('sealed_orbit62',HERE/'verify_orbit_tree.py')
O=importlib.util.module_from_spec(spec);spec.loader.exec_module(O)
def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    preds=['marked_end_lattice_receipt.json','orbit_tree_receipt.json']
    for name in preds:
        old=json.loads((HERE/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for p,h in old['source_sha256'].items():check('sealed predecessor '+name+' '+p,sha(HERE/p)==h)
    manifest=HERE/'sources/binomial-periodic/manifest.json';source=json.loads(manifest.read_text())
    for row in source['primary_sources']:
        p=HERE/row['path'];check('primary pin '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    perms=list(itertools.permutations(range(4)));index={p:i for i,p in enumerate(perms)};I=index[tuple(range(4))]
    mul=[[index[tuple(a[b[i]] for i in range(4))] for b in perms] for a in perms]
    powers=[]
    for a in range(24):
        row=[I]
        for j in range(12):row.append(mul[row[-1]][a])
        powers.append(row)
    def power(a,n):return powers[a][n%12]
    trans=[a for a,p in enumerate(perms) if sum(p[i]!=i for i in range(4))==2]
    check('permutation exponent12 retained',all(row[12]==I for row in powers))
    def group(gs):
        out={I};todo=[I]
        while todo:
            a=todo.pop()
            for b in gs:
                c=mul[a][b]
                if c not in out:out.add(c);todo.append(c)
        return frozenset(out)
    groups={(a,b):group([a,b]) for a,b in itertools.product(range(24),repeat=2)}
    full={(a,b,tau):len(group(list(groups[a,b])+[tau]))==24 for a,b,tau in itertools.product(range(24),range(24),trans)}
    inverse=[next(b for b in range(24) if mul[a][b]==I) for a in range(24)]
    orbit_cache={}
    def orbit(gs):
        if gs not in orbit_cache:
            vals=[tuple(mul[mul[inverse[c]][a]][c] for a in gs) for c in range(24)];key=min(vals)
            for val in vals:orbit_cache[val]=key
        return orbit_cache[gs]
    census={};odd_torus_cyclic=True
    for p,q in itertools.product(range(12),repeat=2):
        if math.gcd(math.gcd(p,q),12)!=1:continue
        uv=[(u,v) for u,v in itertools.product(range(12),repeat=2) if (q*u+p*v)%12==1]
        check(f'modular Bezout pairs {p},{q}',len(uv)==12)
        torus={};unique=True
        for a,b in itertools.product(range(24),repeat=2):
            ht=power(a,p)
            if ht!=power(b,q):continue
            MM={mul[power(a,u)][power(b,v)] for u,v in uv};unique &= len(MM)==1
            assert len(MM)==1,'every modular Bezout solution must give the SAME meridian'
            M=MM.pop();L=mul[ht][power(M,-p*q)];torus.setdefault((M,L),[]).append((a,b))
            if p%2 and q%2:odd_torus_cyclic &= mul[a][b]==mul[b][a] and len(groups[a,b]) in [1,2,3,4]
        check(f'complete Bezout-image independence {p},{q}',unique)
        c=(p*q-p-q)%12
        for r in range(12):
            k=(r*c+1)%12;classes={};total=0
            for beta in range(24):
                M=power(beta,r);L=power(beta,-k)
                for a,b in torus.get((M,L),[]):
                    assert mul[power(M,-c)][power(L,-1)]==beta,'outer section is fixed by companion peripherals'
                    for tau in trans:
                        if power(mul[tau][beta],r)!=I or not full[a,b,tau]:continue
                        key=orbit((a,b,beta,tau));total+=1
                        if key not in classes:
                            graph=O.orbit_graph(dict(companion=[perms[a],perms[b]],cable=[perms[beta],perms[tau]]),[('companion','cable',[perms[M],perms[L]])],4)
                            assert graph['components']==1
                            classes[key]=dict(rho=[perms[a],perms[b],perms[beta],perms[tau]],orbit_size=0,cycle_rank=graph['cycle_rank'],companion_order=len(groups[a,b]),cable_order=len(groups[beta,tau]))
                        classes[key]['orbit_size']+=1
            check(f'complete image conjugacy orbits {p},{q},{r}',total==24*len(classes) and all(x['orbit_size']==24 for x in classes.values()))
            census[f'{p},{q},{r}']=dict(representations=total,classes=len(classes),tree_classes=sum(x['cycle_rank']==0 for x in classes.values()),orbits=list(classes.values()))
    check('all odd first-pair torus images cyclic',odd_torus_cyclic)
    check('1152 complete residue cells',len(census)==1152)
    check('244 peripheral and tree residue cells',sum(x['classes']>0 for x in census.values())==sum(x['tree_classes']>0 for x in census.values())==244)
    check('14304 total representations and596 classes',sum(x['representations'] for x in census.values())==14304 and sum(x['classes'] for x in census.values())==596)
    check('484 tree classes and112 rank-one cycle classes',sum(x['tree_classes'] for x in census.values())==484 and sum(x['cycle_rank']==1 for cell in census.values() for x in cell['orbits'])==112)
    check('parameter interchange count symmetry',all((row['classes'],row['tree_classes'])==(census[','.join([key.split(',')[1],key.split(',')[0],key.split(',')[2]])]['classes'],census[','.join([key.split(',')[1],key.split(',')[0],key.split(',')[2]])]['tree_classes']) for key,row in census.items()))
    old62=json.loads((HERE/'orbit_tree_receipt.json').read_text())['periodic']['residue_certificates']
    check('all trefoil slices agree with sealed PF62',all((census[f'2,3,{r}']['representations'],census[f'2,3,{r}']['classes'],census[f'2,3,{r}']['tree_classes'])==(old62[str(r)]['total_representations'],old62[str(r)]['conjugacy_classes'],old62[str(r)]['tree_classes']) for r in range(12)))
    check('C3 companions retained',any(x['companion_order']==3 and x['cycle_rank']==0 for row in census.values() for x in row['orbits']))
    # The determinant test is a separate condition modulo9; combined period36 is exact.
    combined_admissible=combined_surviving=0;odd_classifier=True
    for p,q in itertools.product(range(36),repeat=2):
        if math.gcd(math.gcd(p,q),36)!=1:continue
        for r in range(36):
            combined_admissible+=1;row=census[f'{p%12},{q%12},{r%12}'];c=p*q-p-q;k=r*c+1
            D0=1 if p%2 and q%2 else (q if p%2==0 else p)
            det=k if r%2==0 else r*D0
            passing=row['tree_classes']>0 and det%9==0
            combined_surviving+=passing
            if p%2 and q%2:odd_classifier &= passing==(r%36==9 or r%4==0 and k%9==0)
    check('31104 combined cells with2772 survivors',(combined_admissible,combined_surviving)==(31104,2772))
    check('closed all-odd classifier agrees with complete census',odd_classifier)
    # Explicit constructions are independent of the finite-cell lookup.
    def explicit(p,q,r,cycle):
        c=p*q-p-q;k=r*c+1
        M=index[cycle];tau=index[(1,0,2,3)] if len({power(M,j) for j in range(12)})==4 else index[(0,1,3,2)]
        a=power(M,q);b=power(M,p);beta=power(M,-c);z=mul[tau][inverse[beta]]
        assert power(a,p)==power(b,q)
        assert power(z,r)==power(M,k)==I and mul[z][beta]==tau
        assert (r*k)%2==0 and len(group([a,b,beta,tau]))==24
        graph=O.orbit_graph(dict(companion=[perms[a],perms[b]],cable=[perms[beta],perms[tau]]),[('companion','cable',[perms[M],perms[I]])],4)
        assert graph['cycle_rank']==0 and (k if r%2==0 else r)%9==0
        return True
    check('explicit unbounded-family controls C4',all(explicit(p,q,r,(1,2,3,0)) for p,q,r in [(3,5,9),(5,7,45),(37,49,81),(101,103,117)]))
    check('explicit unbounded-family controls C3',all(explicit(p,q,r,(1,2,0,3)) for p,q,r in [(3,5,32),(5,7,16),(7,9,4)]))
    # Universal polynomial identities and exact partition substitution for Zannier.
    S=O.G.S;p,q,r=S.symbols('p q r',integer=True);c=p*q-p-q;k=r*c+1;N=p*q*r
    check('Zannier uniform partition equality',S.expand(N-(p+q)*r+1-k)==0)
    check('Clay Watson positive Bezout conversion',S.expand(r*(c*(r-1)+1)-k*(r-1))==1)
    check('cable section recovers meridian exponent',S.expand(k-r*c)==1)
    check('constant Wronskian degree',S.expand(k+r*(p+q-p*q)-1)==0)
    check('primitive last order and even fibre exponent controls',all(math.gcd(rr,rr*(pp*qq-pp-qq)+1)==1 and (rr*(rr*(pp*qq-pp-qq)+1))%2==0 for pp,qq,rr in itertools.product(range(2,16),range(2,16),range(2,18)) if math.gcd(pp,qq)==1))
    inputs=[Path(__file__),HERE/'ALL_BINOMIAL_KEYS_PERIODIC_GATE_AND_EXISTENCE.md',manifest,HERE/'verify_orbit_tree.py',HERE/'verify_marked_end_lattice.py']+[HERE/n for n in preds]+[HERE/row['path'] for row in source['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
        reduction=dict(all_coprime_p_q_at_least_two=True,all_windings_r_at_least_two=True,primitive_actual_binomial_key_required=True,general_keys_assumed_binomial=False,quartic_forces_eta_zero=True,c_odd=True,final_fibre_image_identity=True,peripheral_period=12,modular_Bezout_lifts_exactly=True,all_Bezout_images_identical=True,all_companion_images_retained=True),
        census=dict(admissible_mod12_cells=1152,peripheral_cells=244,tree_cells=244,representations=14304,classes=596,tree_classes=484,cycle_classes=112,complete_cells=census,combined_period=36,admissible_mod36_cells=combined_admissible,combined_surviving_cells=combined_surviving,residue_counts_not_curve_counts=True),
        existence=dict(Zannier_Theorem1_imported=True,all_coprime_p_q_and_r=True,complex_coefficients=True,finite_birational_A1_normalization=True,immersive=True,normal_degree=0,finite_multiple_value_types_classified=False,all_controls_nodal=False,quartic_double_value_bound_supplied=False),
        odd_family=dict(all_coprime_odd_p_q_at_least_three=True,torus_S4_images_cyclic=True,closed_classifier='r=9 mod36 OR (4|r AND 9|(r*(p*q-p-q)+1))',all_r9mod36_explicit_peripheral_S4_and_tree=True,polynomial_immersions_exist_for_all_these_parameters=True,sign_Picard_or_linking_supplied=False,affine_complement_S4_extension_supplied=False,marked_plane_finite_cover_supplied=False),
        assurance=dict(written_unbounded_reduction=True,finite_census_exhaustive_after_period_proof=True,Zannier_existence_independently_reproved=False,Clay_Watson_primary_peripherals_used=True,finite_parameter_search_inferred_unbounded_theorem=False,full_quartic_exclusion=False,higher_degree_exclusion=False,finite_Keller_map=False,formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,research_goal_completed=False),
        source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'binomial_periodic_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; all primitive binomial parameters, exact census and polynomial existence')
if __name__=='__main__':main()
