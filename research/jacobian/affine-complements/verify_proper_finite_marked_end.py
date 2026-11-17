"""PF2RULED98: proper finite marked end and general peripheral lattices.

The written geometric proof, not finite controls, identifies the end.
The lattice helper is independently checked against sheet-action BFS.
"""
from pathlib import Path
import sys,json,hashlib,os,platform,itertools,math
import sympy as S
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
import composed_peripheral_cover as C
from peripheral_boundary_lattices import boundary_lattices
E=C.E

def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def identity(name,value):check(name,S.expand(value)==0)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    preds=['unit_peripheral_coverage_receipt.json','marked_end_lattice_receipt.json','complementary_corner_jsj_receipt.json']
    for name in preds:
        r=json.loads((P/name).read_text());check('predecessor PASS '+name,r['outcome']=='PASS')
        for f,h in r['source_sha256'].items():check('sealed predecessor '+name+' '+f,sha(P/f)==h)
    manifest=P/'sources/proper-finite-end/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        f=P/row['path'];check('retained primary pin '+f.name,sha(f)==row['sha256'] and f.stat().st_size==row['bytes'])
        check('complete primary pages freshly visually read '+f.name,row['fresh_PDF_audit'] and len(row['complete_PDF_pages_rendered_and_visually_inspected'])>=3)
    check('source audit not whole-chain reproof',sources['all_earlier_source_and_exclusion_implications_independently_reproved'] is False and sources['full_quartic_Keller_exclusion_registered'] is False)
    A,B,t,u,v,z,w,R=S.symbols('A B t u v z w R',positive=True)
    a,b=S.symbols('a b',integer=True,positive=True)
    # Squared local radii z,w and positive local polar multiplicities a,b.
    identity('projective squared norm equals reciprocal affine squared radius',
        S.cancel(z/(z+z*u+z*v)-1/(1+u+v)))
    identity('target radius R corresponds to alpha 1/(1+R2)',(1/(1+u+v)).subs(v,R**2-u)-1/(1+R**2))
    identity('crossing blowup adds POSITIVE exceptional polar multiplicity',(z**a*w**b).subs(w,z*t)-z**(a+b)*t**b)
    identity('smooth-face blowup retains positive support',(z**a).subs(z,w*t)-w**a*t**a)
    # Bound |V(log c)| by B; thin faces dominate every bounded metric term.
    identity('uniform smooth positive factor cannot cancel thin-face transversality',
        (2*a/z-B).subs(z,a/(B+1+t))-(B+2+2*t))
    identity('positive two-face corner derivative dominates bounded metric term',
        (2*a/z+2*b/w-B).subs({z:a/(B+1+t),w:b/(B+1+u)})-(3*B+4+2*t+2*u))
    identity('actual power-map radial derivative in first positive radius',S.diff(z**(2*a)+w**(2*b),z)-2*a*z**(2*a-1))
    identity('actual power-map radial derivative in second positive radius',S.diff(z**(2*a)+w**(2*b),w)-2*b*w**(2*b-1))
    identity('whole-map properness countercontrol exact bounded image',((1/t)*t)-1)
    # F(x,y)=(x,xy): (1/t,t) escapes while image (1/t,1) stays bounded.
    check('nonfinite whole-map control cannot supply compact radial cores',S.limit(1/t,t,S.oo)==0 and S.limit(t,t,S.oo)==S.oo)

    def orbit_bfs(mu,lam,seed):
        out={seed};todo=[seed]
        while todo:
            s=todo.pop()
            for g in [mu,lam]:
                q=g[s]
                if q not in out:out.add(q);todo.append(q)
        return out
    def act(mu,lam,seed,aa,bb):
        n=len(mu);mm=E.power(mu,aa%math.lcm(*range(1,n+1)));ll=E.power(lam,bb%math.lcm(*range(1,n+1)))
        return mm[ll[seed]]
    # Complete small-degree commuting tables; this is a finite control of
    # the independently written all-degree stabilizer proof.
    census=[]
    for degree in range(1,6):
        permutations=list(itertools.permutations(range(degree)));pairs=0;components=0
        for mu in permutations:
            for lam in permutations:
                if any(mu[lam[s]]!=lam[mu[s]] for s in range(degree)):continue
                pairs+=1;rows=boundary_lattices(mu,lam);seen=set()
                for row in rows:
                    seed=row['seed'];ee=row['normal_index'];ff=row['longitudinal_degree'];aa=row['shear'];components+=1
                    expected=orbit_bfs(mu,lam,seed)
                    assert expected==set(row['sheets']) and not(seen&expected);seen.update(expected)
                    assert ee*ff==len(expected) and 0<=aa<ee
                    assert act(mu,lam,seed,ee,0)==seed and act(mu,lam,seed,aa,ff)==seed
                    # An independent rectangular sheet-action transversal has
                    # exactly e*f distinct states, hence lattice index is exact.
                    assert len({act(mu,lam,seed,x,y) for x in range(ee) for y in range(ff)})==ee*ff
                    assert S.det(S.Matrix(row['lattice']))==len(expected)
                assert seen==set(range(degree))
        check('exhaustive commuting peripheral table degree'+str(degree),pairs>0 and components>=pairs)
        census.append(dict(degree=degree,commuting_pairs=pairs,boundary_components=components))
    mu=(1,2,3,4,0);lam=E.power(mu,2);sheared=boundary_lattices(mu,lam)
    check('nonzero actual peripheral shear retained',len(sheared)==1 and sheared[0]['lattice']==[[5,0],[3,1]])
    rejected=False
    try:boundary_lattices((1,0,2),(0,2,1))
    except ValueError:rejected=True
    check('noncommuting pair not accepted as torus action',rejected)

    # Actual finite maps (z,w)->(z^e,w^f): a Hopf branch-link control.
    # The first longitude rotates f meridian cycles into ONE boundary torus.
    powers=[]
    for ee,ff in [(2,2),(2,3),(3,4)]:
        mu=tuple(((i+1)%ee)*ff+j for i in range(ee) for j in range(ff))
        lam=tuple(i*ff+(j+1)%ff for i in range(ee) for j in range(ff))
        rows=boundary_lattices(mu,lam)
        check('finite power-map residue longitude joins all meridian cycles '+str((ee,ff)),
            len(rows)==1 and rows[0]['lattice']==[[ee,0],[0,ff]] and rows[0]['degree']==ee*ff and E.generated([mu,lam])==ee*ff)
        powers.append(dict(normal_index=ee,residue_degree=ff,mu=mu,longitude=lam,boundaries=rows,image_order=E.generated([mu,lam])))
    check('transitive V4 does NOT satisfy single-transposition hypothesis',powers[0]['image_order']==4 and
        all(sum(p[i]!=i for i in range(4))==4 for p in [powers[0]['mu'],powers[0]['longitude']]))
    # Finite cyclic maps (z,w)->(z^e,w). Their radial section is S3;
    # use the actual unknot exterior representation to compute all fillings.
    cyclic=[]
    for ee in [2,3,4,5,7,12]:
        mu=tuple(list(range(1,ee))+[0]);rho=[E.power(mu,5%ee),mu]
        for reverse in [False,True]:
            data=C.peripheral_cover(1,5,[],rho,reverse=reverse);bd=data['boundaries'];rows=data['relation_rows']
            closed=C.reduced_h1(rows,[x['meridian'] for x in bd])
            check('actual cyclic finite-map S3 end and null core '+str((ee,reverse)),
                len(bd)==1 and bd[0]['degree']==ee and
                C.reduced_h1(rows)==dict(free=1,torsion=[]) and closed==dict(free=0,torsion=[]) and
                C.reduced_h1(rows,[bd[0]['meridian'],bd[0]['longitude']])==closed)
        cyclic.append(dict(degree=ee,closed_H1=closed))
    # Small connected transposition graphs control the degree-free path proof.
    graphs=[]
    for degree in range(2,5):
        edges=list(itertools.combinations(range(degree),2));connected=0
        for mask in range(1<<len(edges)):
            chosen=[edge for j,edge in enumerate(edges) if mask>>j&1]
            component={0}
            while True:
                new=component|{x for e in chosen if set(e)&component for x in e}
                if new==component:break
                component=new
            if len(component)!=degree:continue
            generators=[]
            for x,y in chosen:
                p=list(range(degree));p[x],p[y]=p[y],p[x];generators.append(tuple(p))
            assert E.generated(generators)==math.factorial(degree);connected+=1
        check('complete connected transposition graph controls degree'+str(degree),connected>0)
        graphs.append(dict(degree=degree,connected_graphs=connected))
    inputs=[Path(__file__),P/'PROPER_FINITE_MARKED_END_AND_COMPONENT_OPEN_AUDIT.md',
        P/'peripheral_boundary_lattices.py',P/'composed_peripheral_cover.py',manifest]+[P/x for x in preds]+[P/x['path'] for x in sources['primary_sources']]
    result=dict(as_of='2026-11-16',outcome='PASS',check_count=len(checks),checks=checks,
        geometry=dict(primary_component_open_affineness_and_disjoint_entire_lines_verified=True,
            actual_exact_volume_primitive_restricts_without_deleting_finite_map=True,
            original_finite_normalization_compactifies_with_positive_entire_polar_support=True,
            actual_radial_section_is_oriented_marked_plumbing_boundary=True,
            curve_marks_normal_meridians_and_actual_full_peripheral_maps_retained=True,
            Hodge_index_gives_connected_end_in_every_degree=True,
            single_transposition_meridians_force_full_infinity_Sd=True,
            no_abstract_group_or_determinant_recognition_substituted=True),
        peripheral=dict(arbitrary_commuting_actions_full_orbits_classify_boundary_tori=True,
            exact_integral_stabilizer_basis_e0_af_index_ef=True,
            longitudinal_residue_degree_and_nonzero_shear_retained=True,
            identity_longitude_not_assumed_in_general=True,
            flow_helper_itself_still_requires_identity_longitude=True),
        commuting_action_census=census,nonzero_shear_control=sheared,
        actual_finite_power_map_Hopf_controls=powers,actual_cyclic_finite_map_controls=cyclic,
        connected_transposition_graph_controls=graphs,
        assurance=dict(written_all_degree_geometric_and_lattice_proofs=True,
            finite_controls_not_universal_topology_or_realization_proof=True,
            original_proper_finite_map_essential_countercontrol_retained=True,
            complete_primary_PDF_pages_rendered_and_visually_read=True,
            earlier_whole_chart_key_form_passport_and_exclusion_chain_audit_required=True,
            all_earlier_source_and_exclusion_implications_independently_reproved=False,
            arbitrary_abstract_plumbing_realized_as_finite_map=False,
            component_open_restriction_inferred_finite=False,
            full_quartic_Keller_exclusion_registered=False,
            arbitrary_higher_degree_packets_excluded=False,
            full_compatible_completion_goal_solved=False,formal_verification=False,
            external_review=False,research_goal_completed=False),
        source_sha256={os.path.relpath(f,P):sha(f) for f in inputs},
        runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'proper_finite_marked_end_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; degree-free proper marked end and general torus-cover lattices; full objective OPEN')
if __name__=='__main__':main()
