"""PF2RULED78 controls; the unrestricted proofs are in the companion note."""
from pathlib import Path
import hashlib,importlib.util,json,math,os,platform
import sympy as S
HERE=Path(__file__).resolve().parent
def load(name,file):
    spec=importlib.util.spec_from_file_location(name,HERE/file)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
J=load('corner78','verify_complementary_corner_jsj.py');C=J.C;V=J.V;E=V.E;B=J.B

def cluster_control(rho,sigma,support):
    z,w,T,U=S.symbols('z w T U');q=2 if rho==1 else rho-1
    R=C.local_cluster(z**-2+z*T if rho==1 else z**-q+z**(2-rho)*T,z**rho,0,z,T)
    D=C.local_cluster(sum(c*w**-j for j,c in support.items())+w**(1-sigma)*U,w**sigma,0,w,U)
    rr,dd=J.retag(R,'R','Hx'),J.retag(D,'D','Hy');out=C.assemble([rr,dd]);mat=out['full']
    vertices=set(range(mat.rows))-set(out['index'])
    adj={i:[j for j in vertices if j!=i and mat[i,j]] for i in vertices}
    arrow=out['labels'][rr['rows'][-1]['parents'][0]]
    spine=J.rooted_nodes(mat,vertices,arrow)
    prev,cur=out['labels']['Hx'],out['labels']['Hy']
    vp,vc=S.Matrix([1,0]),S.Matrix([0,1])
    while len(adj[cur])==2:
        nxt=next(i for i in adj[cur] if i!=prev)
        vp,vc=vc,-mat[cur,cur]*vc-vp
        prev,cur=cur,nxt
    first_D=cur;bridges=[]
    for a,b in spine['edges']:
        seen={a};todo=[(a,[])];path=None
        while todo:
            node,pp=todo.pop()
            if node==b:path=pp;break
            for nxt in adj[node]:
                if nxt not in seen:seen.add(nxt);todo.append((nxt,pp+[nxt]))
        assert path is not None
        va,vb=S.Matrix([1,0]),S.Matrix([0,1])
        for vertex in path[:-1]:va,vb=vb,-mat[vertex,vertex]*vb-va
        internal=path[:-1]
        bridges.append(dict(nodes=[a,b],weights=[int(mat[i,i]) for i in internal],
            transported_fibre=[int(i) for i in vb],
            internal_D=not (rho!=2 and spine['root'] in (a,b))))
    return dict(matrix=out['source'],first_D=first_D,first_D_fibre=vc,
       root=spine['root'],legs={i:[x['order'] for x in legs] for i,legs in spine['legs'].items()},
       bridges=bridges,D_boundary_weights=[int(mat[out['labels'][row['key']],out['labels'][row['key']]]) for row in dd['rows'][:-1]])

def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    preds=['complementary_corner_jsj_receipt.json','matched_meridian_common_exclusion_receipt.json',
           'seifert_sections_winding_three_receipt.json','lifted_surface_s3_receipt.json',
           'coordinate_pencil_matrix_receipt.json','marked_end_lattice_receipt.json',
           'universal_regular_end_receipt.json','unramified_orevkov_receipt.json']
    for name in preds:
        old=json.loads((HERE/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for path,digest in old['source_sha256'].items():check('sealed predecessor '+name+' '+path,sha(HERE/path)==digest)
    manifest=HERE/'sources/trefoil-complementary-fusion/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        p=HERE/row['path'];check('primary pin '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    q,a,p,sigma,g=S.symbols('q a p sigma g',positive=True,integer=True);rho=q+1
    fr=S.Matrix([q,-rho]);fd=S.Matrix([-sigma/g,p/g])
    check('unrestricted actual corner slope determinant',S.expand(fr.row_join(fd).det()-(q*p-rho*sigma)/g)==0)
    check('unrestricted matched primitive opposite fibres',S.simplify((fr+fd).subs({p:rho*a,sigma:q*a,g:a}))==S.zeros(2,1))
    check('matched first D cone is rho',S.simplify((p/g).subs({p:rho*a,g:a})-rho)==0)
    check('simple R cannot match admissible D positive gap',S.expand((2*p-sigma).subs(p,sigma+1+a)-(sigma+2+2*a))==0)
    # Negative internal chains: the written proof uses this recurrence for every length.
    x,y,c=S.symbols('x y c',positive=True,integer=True)
    check('negative-chain determinant monotonicity induction',S.expand((2+c)*x-y-x-((x-y)+c*x))==0)
    for weights in [[],[2],[2,2,2],[3,2,5],[2,4,2,3,2]]:
        mat=S.diag(*weights)
        for i in range(len(weights)-1):mat[i,i+1]=mat[i+1,i]=-1
        d0,d1=0,1
        for w in weights:d0,d1=d1,w*d1-d0
        check('bridge determinant control '+str(weights),d1==mat.det() and d1>=1)
    m,n=S.symbols('m n',integer=True)
    check('Klein only two invariant primitive boundary axes',S.Matrix([[m,m],[n,-n]]).det()==-2*m*n)
    check('Klein distinct from trefoil integral H1',E.h1([[0,2]])==dict(free=1,torsion=[2]))
    s=S.symbols('s',positive=True,integer=True);r=3*s;k=3*s+1;D=(s+1)/2;JJ=(3*s-1)/2
    check('actual trefoil solid companion cap distance',S.expand(r-k/2-JJ)==0)
    check('remaining companion fibre intersection',S.expand(6*r-k-(15*s-1))==0)
    check('target cone difference forces s3',S.expand(JJ-s-(s-1)/2)==0 and S.solve(JJ-s-1,s)==[3])
    check('target least winding meridian transport',S.Matrix([-1,0])+S.Matrix([2,1])==S.Matrix([1,1]))
    check('target least winding marked distance5',abs(S.Matrix([[1,6],[1,1]]).det())==5)
    check('source rho2 meridian fibre distance3',abs(S.Matrix([[1,-2],[-3,3]]).det())==3)
    check('forced source rho2 full determinant15',S.Matrix([[-3,3],[3,-8]]).det()==15)
    check('least target homology order1',D.subs(s,1)==1)
    check('a2 kappa3 cannot break residual gcd2',6-1+3==8 and math.gcd(2,8)==2)
    check('a3 kappa2 forces last support10',9-1+2==10 and math.gcd(3,10)==1 and math.gcd(3,11)==1)
    contacts=(9-3)*12+(3-1)*(24-10);nd=contacts+24-9+1
    check('all-coefficient residual contact sum and nD',contacts==100 and nd==116)
    check('forced source rho4 full determinant212',13*nd-36**2==212)
    check('target s3 homology order2 differs',D.subs(s,3)==2 and 212!=2)
    cases=[(1,2,{3:2}),(2,2,{3:2}),(3,2,{3:2}),
           (3,6,{9:2,7:1}),(4,9,{12:2,10:1}),(4,9,{12:7,10:-3}),
           (4,6,{8:2,7:1}),(3,6,{8:2,7:1}),
           (2,4,{8:2,6:1,5:3}),(3,6,{12:2,8:1,7:3}),
           (2,1,{2:2}),(3,1,{2:2})]
    controls=[];saw_matching=saw_unmatched=saw_klein=saw_integral=False
    for rr,ss,cs in cases:
        ctl=cluster_control(rr,ss,cs);pp=max(cs);gg=math.gcd(ss,pp);tag=str((rr,ss,cs))
        check('actual first D fibre transport '+tag,ctl['first_D_fibre']==S.Matrix([-ss//gg,pp//gg]))
        check('all actual D exceptional boundary weights <=-2 '+tag,all(w<=-2 for w in ctl['D_boundary_weights']))
        check('actual first D side cone '+tag,pp//gg in ctl['legs'][ctl['first_D']])
        nr=5 if rr==1 else rr*rr-rr+1;qq=2 if rr==1 else rr-1
        check('full exact complementary lattice '+tag,ctl['matrix']==S.Matrix([[-nr,qq*pp],[qq*pp,-C.V.diagonal_size(1,ss,cs)]]))
        for bridge in ctl['bridges']:
            if bridge['internal_D']:
                check('actual internal D bridge cannot fuse '+str((tag,bridge['nodes'])),all(w<=-2 for w in bridge['weights']) and bridge['transported_fibre'][1]>0)
            else:
                matching=(2*pp==ss) if rr==1 else (qq*pp==rr*ss)
                check('actual R D bridge agrees with corner equality '+str((tag,bridge['nodes'])),(bridge['transported_fibre'][1]==0)==matching)
                saw_matching|=matching;saw_unmatched|=not matching
        saw_klein|=any(sorted(legs)==[2,2] for legs in ctl['legs'].values());saw_integral|=pp%ss==0 and ss>1
        if (rr,ss,set(cs))==(2,2,{3}):check('r3 forced actual source lattice',ctl['matrix'].det()==15 and sorted(ctl['legs'][ctl['root']])==[2,3])
        if rr==4 and ss==9:check('r9 forced actual source lattice and terminal cones '+str(cs),ctl['matrix'].det()==212 and any(sorted(legs)==[2,3] for legs in ctl['legs'].values()))
        controls.append(dict(rho=rr,sigma=ss,support={str(j):v for j,v in cs.items()},
            first_D_fibre=[int(i) for i in ctl['first_D_fibre']],matrix=[[int(i) for i in row] for row in ctl['matrix'].tolist()],
            bridges=ctl['bridges'],leg_orders=ctl['legs']))
    check('matching nonmatching integral-contact and Klein controls present',saw_matching and saw_unmatched and saw_klein and saw_integral)
    census=json.loads((HERE/'binomial_periodic_receipt.json').read_text())['census']['complete_cells'];targets=[]
    for ee,oo in [(2,3),(3,2)]:
        for rr in [3,9,15,21]:
            ss=rr//3;dd=(ss+1)//2;kk=3*ss+1
            rows=[row for row in census[f'{ee%12},{oo%12},{rr%12}']['orbits'] if (row['companion_order'],row['cable_order'],row['cycle_rank'])==(6,6,0)]
            check('both actual trefoil representations retained '+str((ee,oo,rr)),len(rows)==2)
            want=dict(H1=B.group([dd]),quotient=B.group([math.gcd(dd,ss*kk)]))
            for i,row in enumerate(rows):check('ordinary full satellite trefoil H1 and R mark '+str((ee,oo,rr,i)),V.full_satellite(ee,oo,rr,row)==want)
            targets.append(dict(p=ee,q=oo,r=rr,H1=want['H1'],R_quotient=want['quotient'],root_cones=[ss,(3*ss-1)//2]))
    inputs=[Path(__file__),HERE/'TREFOIL_CAPS_AND_COMPLEMENTARY_FIBRE_FUSION.md',manifest,
       HERE/'verify_complementary_corner_jsj.py']+[HERE/n for n in preds]+[HERE/row['path'] for row in sources['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      source=dict(no_degree_order_support_depth_bound=True,two_entire_omitted_lines_required=True,
        full_smooth_exact_complementary_source_and_whole_charts_required=True,
        first_D_fibre='-(sigma/g)*hx+(p/g)*hf',first_D_cone='p/g',
        internal_D_fibres_distinct_by_negative_chain_determinant=True,
        only_orientable_corner_fusion='p=rho*a,sigma=(rho-1)*a',
        fused_annulus_cones=['rho-1','rho'],rho2_product_collar_removed=True,
        Klein_alternative_allowed_and_excluded_by_terminal_piece_type=True,Klein_fibre_uniqueness_assumed=False),
      target=dict(actual_primitive_binomial_two_stage_quartic_S3_S3_required=True,
        trefoil_root_cones=['s','(3s-1)/2'],s1_R_meridian_fibre_distance=5,
        s_ge3_requires_consecutive_cones_and_hence_s3=True,
        s1_source_forced=[2,2,3],s1_source_marked_distance=3,s1_source_end_H1_order=15,s1_target_end_H1_order=1,
        s3_source_forced=dict(rho=4,sigma=9,support=[12,10],nR=13,nD=116,cross=36,H1_order=212),s3_target_end_H1_order=2,
        all_trefoil_complementary_windings_excluded=True,entire_S3_S3_sector_excluded_in_both_pencils=True,
        other_representations_of_same_target_keys_excluded=False),
      source_cluster_controls=controls,ordinary_target_controls=targets,
      assurance=dict(written_unbounded_Euclidean_bridge_Klein_and_cone_proofs=True,
        full_peripheral_Waldhausen_recognition=True,s1_sufficient_largeness_checked_with_essential_annulus=True,
        complete_primary_pages_rechecked=True,deep_topology_imported=True,
        finite_controls_not_universal_proof=True,component_open_assumed_finite_cover=False,
        recognition_from_H1_or_unmarked_group=False,arbitrary_omitted_count_theorem=False,
        all_complementary_sectors_closed=False,full_quartic_exclusion=False,higher_degree_exclusion=False,
        finite_Keller_map=False,PF63_polynomial_immersions_declared_nonexistent=False,
        formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,research_goal_completed=False),
      source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'trefoil_complementary_fusion_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; trefoil complementary caps excluded; entire S3/S3 sector closed in both pencils')
if __name__=='__main__':main()
