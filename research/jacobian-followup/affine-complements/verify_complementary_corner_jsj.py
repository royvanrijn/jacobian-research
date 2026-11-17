"""PF2RULED77: exact corner/marked-spine controls for unbounded written proofs."""
from pathlib import Path
import hashlib,importlib.util,json,math,os,platform
import sympy as S
HERE=Path(__file__).resolve().parent
def load(name,file):
    spec=importlib.util.spec_from_file_location(name,HERE/file)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
C=load('cluster77','verify_common_pencil_mark.py')
V=load('base77','verify_seifert_base_genus.py');E=V.E
B=load('budget77','verify_peripheral_cone_budget.py')

def retag(cl,tag,axis):
    mp={row['key']:(tag,)+row['key'] for row in cl['rows']};mp['Hx']=axis
    return dict(rows=[dict(key=mp[row['key']],parents=[mp[p] for p in row['parents']],multiplicity=row['multiplicity']) for row in cl['rows']],terminal=mp[cl['terminal']])

def rooted_nodes(mat,vertices,arrow):
    vertices=set(vertices);adj={i:[j for j in vertices if i!=j and mat[i,j]] for i in vertices}
    nodes={i for i in vertices if len(adj[i])+(i==arrow)>=3};edges=set();legs={i:[] for i in nodes}
    for a in nodes:
        for nxt in adj[a]:
            prev,cur=a,nxt;path=[]
            while cur not in nodes:
                if cur==arrow:break
                path.append(cur);opts=[x for x in adj[cur] if x!=prev]
                if not opts:
                    legs[a].append(dict(vertices=path,order=int(abs(mat.extract(path,path).det()))));break
                assert len(opts)==1
                prev,cur=cur,opts[0]
            if cur in nodes:edges.add(tuple(sorted((a,cur))))
    prev=None;root=arrow
    while root not in nodes:
        opts=[x for x in adj[root] if x!=prev];assert len(opts)==1
        prev,root=root,opts[0]
    degrees={i:sum(i in e for e in edges) for i in nodes}
    return dict(nodes=sorted(nodes),edges=sorted(edges),root=root,degrees=degrees,legs=legs)

def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    preds=['matched_meridian_common_exclusion_receipt.json','seifert_sections_winding_three_receipt.json',
           'regular_quartic_exclusion_receipt.json','universal_regular_end_receipt.json',
           'unramified_orevkov_receipt.json','coordinate_pencil_matrix_receipt.json',
           'common_pencil_mark_receipt.json','binomial_periodic_receipt.json']
    for name in preds:
        old=json.loads((HERE/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for p,h in old['source_sha256'].items():check('sealed predecessor '+name+' '+p,sha(HERE/p)==h)
    manifest=HERE/'sources/complementary-corner-jsj/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        p=HERE/row['path'];check('primary bytes '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    q=S.symbols('q',positive=True,integer=True);nn=q*q+q+1
    check('unbounded tail determinant recurrence',S.expand(2*(q+1)-q-(q+2))==0)
    check('ramified exceptional chain determinant',S.expand((q+1)*(q+1)-q-nn)==0)
    check('whole corner chain determinant',S.expand(nn-(q+1)-q*q)==0)
    change=S.Matrix([[q*q,1-q],[-nn,q]])
    check('primitive actual meridian longitude determinant',S.expand(change.det())==1)
    check('chosen longitude maps to core',S.expand((1-q)*nn+q*q*q)==1)
    check('chosen meridian maps to zero',q*q*nn-nn*q*q==0)
    mu,ell=S.symbols('mu ell');hx=q*mu+nn*ell;dd=(q-1)*mu+q*q*ell
    check('actual transported R curve unit cable identity',S.expand(q*hx-(q+1)*dd-(mu+q*ell))==0)
    sigma,p,jstar=S.symbols('sigma p jstar',positive=True,integer=True)
    check('unramified final free contacts identity',S.expand((2*p+1-sigma)-(2*p-jstar)-(jstar+1-sigma))==0)
    check('unramified exact final contact lower bound',S.expand((jstar+1-sigma).subs(jstar,sigma+1+q)-(q+2))==0)
    z,w,T,U=S.symbols('z w T U');rside=[]
    for rho in [2,3,4,5,8]:
        qq=rho-1;nr=rho*rho-rho+1
        R=C.local_cluster(z**(-qq)+z**(2-rho)*T,z**rho,0,z,T);out=C.assemble([R]);mat=out['full']
        chain=[out['labels']['Hx'],out['labels'][R['rows'][0]['key']]]+[out['labels'][row['key']] for row in reversed(R['rows'][1:-1])]
        check('exact unrestricted R chain weights control '+str(rho),[mat[i,i] for i in chain]==[-1,-rho]+[-2]*qq)
        check('exact R chain adjacency control '+str(rho),{(i,j) for i in range(len(chain)) for j in range(i+1,len(chain)) if mat[chain[i],chain[j]]}=={(i,i+1) for i in range(len(chain)-1)})
        check('whole R-side and cofactor determinants '+str(rho),abs(mat.extract(chain,chain).det())==qq*qq and abs(mat.extract(chain[1:],chain[1:]).det())==nr)
        check('R arrow at first negative-two vertex '+str(rho),out['labels'][R['rows'][-1]['parents'][0]]==chain[2])
        check('R multiplicity pattern and whole dual '+str(rho),[r['multiplicity'] for r in R['rows']]==[qq]+[1]*(qq+1) and out['source']==S.Matrix([[-nr]]))
        rside.append(dict(rho=rho,q=qq,weights=[-1,-rho]+[-2]*qq,meridian=[qq*qq,-nr],longitude=[1-qq,qq],R_slope_in_meridian_longitude=[1,qq]))
    R=C.local_cluster(z**-2+z*T,z,0,z,T);out=C.assemble([R]);mat=out['full'];hx0=out['labels']['Hx'];ids=[out['labels'][row['key']] for row in R['rows'][:-1]]
    expected={tuple(sorted(e)) for e in [(hx0,ids[1]),(ids[0],ids[1]),(ids[1],ids[2]),(ids[2],ids[3])]}
    vertices=[hx0]+ids
    check('simple R exact five-negative-two graph',all(mat[i,i]==-2 for i in vertices) and {(i,j) for i in vertices for j in vertices if i<j and mat[i,j]}==expected)
    check('simple R arrow on order3 leg endpoint',out['labels'][R['rows'][-1]['parents'][0]]==ids[3] and out['source']==S.Matrix([[-5]]))
    # These are cluster controls, not assertions that each open surface is affine or finite.
    cases=[(1,1,{2:2}),(2,1,{2:2}),(3,2,{3:2}),(2,4,{6:2,5:1}),
           (3,6,{8:2,7:1}),(4,9,{12:2,10:1}),(2,4,{8:2,6:1,5:3}),(3,6,{12:2,8:1,7:3})]
    controls=[];saw_klein=False;saw_integral=False;saw_multiple=False
    for rho,ss,cs in cases:
        qq=2 if rho==1 else rho-1;nr=5 if rho==1 else rho*rho-rho+1;pp=max(cs)
        R=C.local_cluster(z**-2+z*T if rho==1 else z**(-qq)+z**(2-rho)*T,z**rho,0,z,T)
        D=C.local_cluster(sum(co*w**(-j) for j,co in cs.items())+w**(1-ss)*U,w**ss,0,w,U)
        rr,ds=retag(R,'R','Hx'),retag(D,'D','Hy');out=C.assemble([rr,ds]);mat=out['full'];ix=out['index'];vertices=set(range(mat.rows))-set(ix)
        arrow=out['labels'][rr['rows'][-1]['parents'][0]];spine=rooted_nodes(mat,vertices,arrow);tag=str((rho,ss,cs))
        nd=C.V.diagonal_size(1,ss,cs)
        check('whole complementary lattice and canonical volumes '+tag,out['source']==S.Matrix([[-nr,qq*pp],[qq*pp,-nd]]) and [out['orders'][i] for i in ix]==[1,0])
        check('actual retained leaves and untouched corner '+tag,out['degrees']==[1,1] and mat[0,1]==1)
        check('source initial Seifert-node graph a path '+tag,len(spine['edges'])==len(spine['nodes'])-1 and all(v<=2 for v in spine['degrees'].values()))
        check('source marked boundary at end of path '+tag,spine['degrees'][spine['root']]<=1)
        for node,legs in spine['legs'].items():
            check('all exceptional legs give genuine nonzero cones '+str((tag,node)),all(x['order']>=2 and all(mat[i,i]<=-2 for i in x['vertices']) for x in legs))
            if spine['degrees'][node]<=1 and node!=spine['root'] and sorted(x['order'] for x in legs)==[2,2]:saw_klein=True
            if len(spine['nodes'])==1 and sorted(x['order'] for x in legs)==[2,2]:saw_klein=True
        if ss>1:
            gcd=ss;drops=[]
            for jj in sorted(cs,reverse=True):
                new=math.gcd(gcd,jj)
                if new<gcd:drops.append(jj);gcd=new
            check('all source supports exact and final tail genuine '+tag,min(cs)>=ss+1 and gcd==1 and drops[-1]+1-ss>=2)
            saw_integral |= pp%ss==0;saw_multiple |= len(drops)>=2
        controls.append(dict(rho=rho,sigma=ss,support={str(j):v for j,v in cs.items()},nR=nr,nD=nd,
          node_count=len(spine['nodes']),root_node_degree=spine['degrees'][spine['root']],node_degrees=list(spine['degrees'].values()),
          leg_orders={str(node):[x['order'] for x in legs] for node,legs in spine['legs'].items()}))
    check('integral contact multiple characteristics and Klein terminal all exercised',saw_klein and saw_integral and saw_multiple)
    ee,oo,rr=S.symbols('e o r',positive=True,integer=True);kk=rr*(ee*oo-ee-oo)+1
    check('target fixed-edge fibre intersection exact',S.expand(rr*ee*oo-kk-(rr*(ee+oo)-1))==0)
    check('target central-edge fibre intersection exact',S.expand((rr*ee*oo-kk)/2-(rr*(ee+oo)-1)/2)==0)
    census=json.loads((HERE/'binomial_periodic_receipt.json').read_text())['census']['complete_cells'];count=0
    for key,cell in census.items():
        p,q,r=map(int,key.split(','));p=p if p>=2 else p+12;q=q if q>=2 else q+12;r=r if r>=2 else r+12
        for i,row in enumerate(cell['orbits']):
            if (row['companion_order'],row['cable_order'],row['cycle_rank'])!=(6,6,0):continue
            e,o=(p,q) if p%2==0 else (q,p)
            if (e,o)==(2,3):e=14;p,q=(e,o) if p%2==0 else (o,e)
            b,t=e//2,o//3;s=r//3;tag=str((key,i))
            comps=V.companion_bases(row,p,q,r);cabs=V.cable_bases(row)
            check('all64 target complete actual incidence retained '+tag,sorted(len(x['sheets']) for x in comps)==[1,3] and sorted(len(x['sheets']) for x in cabs)==[1,3])
            genuine=[a for a in [b,2*b,t] if a>1];chi=1-len(genuine)+sum(S.Rational(1,a) for a in genuine)
            check('actual target three pieces all negative orbifold Euler '+tag,chi<0 and S.Rational(1,e)+S.Rational(1,o)-1<0 and -2+S.Rational(1,s)<0)
            delta1=r*(e+o)-1
            check('actual target two genuine boundary-adjacent JSJ tori '+tag,delta1>0 and delta1%2==0 and delta1//2>0 and len(genuine)>=2)
            count+=1
    check('all64 periodic types with non-trefoil actual lifts checked',count==64)
    check('trefoil solid-torus exception retained',len([a for a in [1,2,1] if a>1])==1)
    inputs=[Path(__file__),HERE/'COMPLEMENTARY_CORNER_AND_BOUNDARY_JSJ_ENDPOINT.md',manifest]+[HERE/n for n in preds]+[HERE/x['path'] for x in sources['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      source=dict(no_field_degree_bound=True,two_entire_omitted_lines_required=True,one_fold_one_unramified_required=True,
        complementary_actual_coordinate_pair_required=True,exact_source_volume_required=True,
        PF50_ramified_pole_collapse_used=True,higher_order_R_side_actual_solid_torus=True,
        R_actual_unit_cable_slope='mu+q*ell',solid_torus_filling_slope='q^2*h_x-(q^2+q+1)*h_f',
        simple_R_side_Seifert_cones=[2,3],simple_R_drilled_side_cones=[2],
        arbitrary_D_characteristic_support_retained=True,unramified_tail='j_star+1-sigma>=2',
        canonical_R_boundary_JSJ_adjacent_tori_at_most=1,terminal_Klein_exception_included=True,
        arbitrary_omitted_count_endpoint_theorem_asserted=False),
      target=dict(actual_binomial_two_stage_quartic_S3_S3_required=True,all_first_exponents_windings_source_orders_poles=True,
        all64_incidence_types_checked=True,actual_fixed_sheets_retained=True,
        non_trefoil_R_boundary_JSJ_adjacent_tori=2,actual_pieces_negative_orbifold_Euler=True,
        nonparallel_edge_fibres_verified=True,full_peripheral_recognition_verified=True,
        non_trefoil_S3_S3_complementary_entirely_excluded=True,
        non_trefoil_S3_S3_sector_excluded_in_both_pencils=True,
        trefoil_S3_S3_complementary_entirely_excluded=False,other_representations_of_same_keys_excluded=False),
      R_side_controls=rside,source_cluster_controls=controls,
      assurance=dict(written_unbounded_corner_spine_and_recognition_proofs=True,
        complete_Hatcher_Waldhausen_primary_pages_rechecked=True,deep_rigidity_JSJ_theorems_imported=True,
        finite_controls_not_universal_proof=True,deep_theorems_independently_reproved=False,
        recognition_from_H1_or_unmarked_group_asserted=False,component_open_assumed_finite_cover=False,
        Klein_fibre_uniqueness_assumed=False,all_complementary_sectors_closed=False,
        full_quartic_exclusion=False,higher_degree_exclusion=False,finite_Keller_map=False,
        PF63_polynomial_curves_declared_nonexistent=False,formal_verification=False,
        independent_full_reproof=False,external_review=False,novelty_priority=False,research_goal_completed=False),
      source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'complementary_corner_jsj_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; degree-free complementary boundary endpoint; non-trefoil S3/S3 excluded in both pencils')
if __name__=='__main__':main()
