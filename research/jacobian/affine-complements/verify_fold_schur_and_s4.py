"""PF2RULED71: exact controls for the unbounded fold and peripheral proofs."""
from pathlib import Path
import hashlib,importlib.util,json,math,os,platform
import sympy as S
HERE=Path(__file__).resolve().parent
def load(name,file):
    spec=importlib.util.spec_from_file_location(name,HERE/file);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
C=load('common67','verify_common_pencil_mark.py')
V=load('base69','verify_seifert_base_genus.py')
B=load('cone70','verify_peripheral_cone_budget.py');E=V.E

def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    preds=['marked_end_lattice_receipt.json','binomial_periodic_receipt.json',
           'common_pencil_mark_receipt.json','seifert_base_genus_receipt.json','peripheral_cone_budget_receipt.json']
    for name in preds:
        old=json.loads((HERE/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for p,h in old['source_sha256'].items():check('sealed predecessor '+name+' '+p,sha(HERE/p)==h)
    manifest=HERE/'sources/fold-schur-s4/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        p=HERE/row['path'];check('retained primary pin '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    # Local sign invariant: free/crossing rules cannot create an opposite-sign edge.
    for a in range(-4,5):
        check('free edge never crosses sign through zero '+str(a),a*(a+1)>=0)
        for b in range(-4,5):
            if a*b<0:continue
            check('crossing edge sign invariant '+str((a,b)),a*(a+b)>=0 and b*(a+b)>=0)
            if a+b==1:check('b1 crossing has one zero one positive parent '+str((a,b)),sorted([a,b])==[0,1])
    z,w,T,W=S.symbols('z w T W');old=json.loads((HERE/'common_pencil_mark_receipt.json').read_text())
    controls=[]
    for c in old['exact_blowup_controls']:
        if c['index']!=2:continue
        hr={int(k):v for k,v in c['R_support'].items()};hd={int(k):v for k,v in c['D_support'].items()}
        r,s=c['r'],c['s'];lam=c['D_value']
        R=C.local_cluster(sum(v*z**-k for k,v in hr.items())+z**(2-r)*T,z**r,0,z,T)
        D=C.local_cluster(sum(v*w**-k for k,v in hd.items())+w**(1-s)*W,lam+w**s,lam,w,W)
        out=C.assemble([R,D]);G=-out['source'];tag=str((r,hr,s,hd,lam))
        # Check the unnormalized curvette vector recurrence on every created vertex.
        records={row['key']:row for cl in [R,D] for row in cl['rows']};n=len(out['keys'])
        vectors={'Hx':S.zeros(1,n),'Hy':S.zeros(1,n)}
        for j,key in enumerate(out['keys']):
            vector=S.zeros(1,n);vector[0,j]=1
            for parent in records[key]['parents']:vector+=vectors[parent]
            vectors[key]=vector
        check('curvette recurrence equals whole extracted Gram '+tag,
              S.Matrix.vstack(vectors[R['terminal']],vectors[D['terminal']])==out['multiplicities'])
        check('all boundary edges preserve log-order signs '+tag,
              all((out['orders'][i]+1)*(out['orders'][j]+1)>=0 for i in range(out['full'].rows)
                  for j in range(i) if out['full'][i,j]))
        check('fold R free parent is b1 and distinct from D '+tag,
              len(R['rows'][-1]['parents'])==1 and R['rows'][-1]['parents'][0]!=D['terminal'] and
              out['orders'][out['labels'][R['rows'][-1]['parents'][0]]]==0)
        schur=S.cancel(1/G.inv()[0,0]);check('strengthened strict fold Schur '+tag,schur>2 and 0<G.inv()[0,0]<S.Rational(1,2))
        controls.append(dict(source_matrix=out['source'].tolist(),Schur=str(schur)))
    sharp=[]
    for n in [0,1,2,4]:
        R=C.local_cluster(z**-1+T,z*z,0,z,T)
        Ds=[C.local_cluster((i+1)*w**-1+W,w,0,w,W) for i in range(n)]
        out=C.assemble([R]+Ds);G=-out['source'];schur=S.cancel(1/G.inv()[0,0])
        check('sharp family has whole free leaves '+str(n),out['degrees']==[1]*(n+1) and [out['orders'][i] for i in out['index']]==[1]+[0]*n)
        check('sharp family exact all-count Schur '+str(n),schur==2+S.Rational(1,n+1))
        sharp.append(dict(unramified_count=n,Gram=G.tolist(),Schur=str(schur)))
    N,L=S.symbols('N L',positive=True)
    check('unrestricted sibling projection Schur',S.simplify((L+2)-N*L**2/(1+N*L)-(2+L/(1+N*L)))==0)
    # The integer Seifert presentations determine the needed torsion groups.
    for t in [1,3,5,9]:
        alpha=[2,2,t,3*t];mat=S.zeros(4,5)
        for i,a in enumerate(alpha):mat[i,0]=1;mat[i,i+1]=a
        check('companion Seifert torsion group control '+str(t),E.h1(mat.tolist())==B.group([2*t],free=1))
    for s in [1,3,5,9]:
        check('cable filled Seifert torsion group control '+str(s),E.h1([[1,s,0],[1,0,3*s]])==B.group([s],free=1))
    # Exact polynomial relation certificate: 4 times the R mark is zero.
    t,r,eps,s=S.symbols('t r eps s',integer=True);mat=S.Matrix([[0,2*t,0,0],[0,0,0,s],[4,eps,0,0],[-6*t,0,1,0]])
    mark=S.Matrix([[0,0,r,0]])
    check('unrestricted integral order-divides4 certificate',6*r*t*mat.row(2)-3*r*eps*mat.row(0)+4*r*mat.row(3)==4*mark)
    abstract=[]
    for tv,sv,epsilon in [(1,3,1),(3,1,3),(5,3,5),(3,9,1),(5,15,7)]:
        sub={t:tv,s:sv,r:3*sv,eps:epsilon};rows=mat.subs(sub).tolist();m=list(mark.subs(sub))
        check('integral peripheral gluing group '+str((tv,sv,epsilon)),E.h1(rows)==B.group([8*tv,sv]))
        check('integral peripheral gluing quotient mark '+str((tv,sv,epsilon)),E.h1(rows,[m])==B.group([2*tv,sv]))
        abstract.append(dict(t=tv,s=sv,epsilon=epsilon,H1=E.h1(rows),quotient=E.h1(rows,[m])))
    for u in [1,3,5,7]:
        for a in [2,6]:check('order4 element in nonsingular Z8 form has self linking half '+str((u,a)),S.Rational(u*a*a,8)%1==S.Rational(1,2))
    # Complete planar S4/S4 type audit; actual q,r remain unrestricted.
    census=json.loads((HERE/'binomial_periodic_receipt.json').read_text())['census']['complete_cells'];total=0
    for key,cell in census.items():
        p,q,r=map(int,key.split(','));p=p if p>=2 else p+12;q=q if q>=2 else q+12;r=r if r>=2 else r+12
        for i,row in enumerate(cell['orbits']):
            if row['companion_order']!=24 or row['cable_order']!=24 or row['cycle_rank']:continue
            if sum(x['genus'] for x in V.companion_bases(row,p,q,r)):continue
            a,b,beta,tau=map(tuple,row['rho']);even,odd=(a,b) if p%2==0 else (b,a)
            check('all planar S4 S4 cycle types '+str((key,i)),B.lengths(even)==[1,1,2] and B.lengths(odd)==[1,3] and B.lengths(E.mul(E.inv(tau),beta))==[1,3])
            total+=1
    check('all64 planar S4 S4 periodic types audited',total==64)
    full=[]
    for q,r in [(3,9),(9,3),(15,9),(9,27),(15,45)]:
        tv=q//3;sv=r//3;check('actual winding peripheral congruence '+str((q,r)),r*(q+2)%4==1)
        j=(1-r*(q+2))//4;check('cable m torsion killed by unit coefficient '+str((q,r)),math.gcd(j,sv)==1)
        rows=[x for x in census[f'2,{q%12},{r%12}']['orbits'] if x['companion_order']==x['cable_order']==24 and not x['cycle_rank']]
        for i,row in enumerate(rows):
            ans=V.full_satellite(2,q,r,row)
            check('full ordinary satellite group and R quotient '+str((q,r,i)),ans==dict(H1=B.group([8*tv,sv]),quotient=B.group([2*tv,sv])))
            if i==0:check('independent spanning tree full satellite '+str((q,r)),V.full_satellite(2,q,r,row,True)==ans)
            full.append(dict(q=q,r=r,class_index=i,**ans))
    inputs=[Path(__file__),HERE/'FOLD_SCHUR_BOUND_AND_S4_COMMON_EXCLUSION.md',manifest]+[HERE/n for n in preds]+[HERE/x['path'] for x in sources['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
        source=dict(all_finite_field_degrees=True,any_unramified_count=True,sole_volume_order_exactly1_required=True,
          smooth_affine_and_entire_disjoint_A1_required=True,common_regular_f_and_global_P1_x_required=True,
          retained_free_leaf_hypothesis_required=True,no_positive_negative_boundary_edge=True,
          no_zero_or_new_free_b1_descends_positive=True,b1_parent_squared_distance_strictly_greater1=True,
          fold_Schur_strictly_greater2=True,R_self_linking_strictly_between0_and_half=True,
          arbitrary_higher_volume_order_bound_asserted=False,complementary_source_bound_asserted=False),
        target=dict(actual_primitive_binomial_two_stage_quartic_S4_S4_required=True,
          even_first_exponent2_required_by_PF70=True,q='3t',r='3s',s_t_positive_odd=True,
          companion_integral_H1='Z plus Z/(2t)',companion_H='6t*g',companion_m='4g+epsilon*delta, epsilon odd',
          filled_cable_integral_H1='Z*H plus Z/s',filled_cable_m_zero=True,filled_cable_h='rH',
          end_H1='Z/(8t) plus Z/s',marked_quotient='Z/(2t) plus Z/s',R_exact_order=4,
          R_has_no_odd_primary_component=True,two_primary_group_cyclic_order8=True,R_self_linking='1/2',
          entire_S4_S4_common_source_sector_excluded=True,all_first_exponents_windings_source_orders_poles=True,
          complementary_source_sector_excluded=False,remaining_complementary_gate='even first exponent2 and gcd(s,t)=1'),
        inherited_fold_controls=controls,sharp_family=sharp,abstract_integral_gluing=abstract,full_satellite_controls=full,
        assurance=dict(written_unbounded_log_order_ancestry_and_projection_proof=True,
          written_integral_peripheral_homology_proof=True,retained_primary_audits_reused=True,
          finite_controls_not_universal_proof=True,
          arbitrary_noncylclic_pieces_replaced=False,full_linking_square_class_asserted=False,
          affine_source_existence_from_sharp_cluster_asserted=False,full_quartic_exclusion=False,
          higher_degree_exclusion=False,finite_Keller_map=False,PF63_curves_declared_nonexistent=False,
          formal_verification=False,independent_full_reproof=False,external_review=False,
          novelty_priority=False,research_goal_completed=False),
        source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    # Sympy integers in geometry matrices are converted once at the serialization boundary.
    (HERE/'fold_schur_s4_receipt.json').write_text(json.dumps(report,indent=2,default=lambda x:int(x) if isinstance(x,S.Integer) else str(x))+'\n')
    print(f'PASS {len(checks)} controls; strict fold Schur>2; entire S4/S4 common-source sector excluded')
if __name__=='__main__':main()
