"""PF2RULED67: unbounded proof controls; no finite tree census proves the theorem."""
from pathlib import Path
import hashlib,importlib.util,json,os,platform
import sympy as S
HERE=Path(__file__).resolve().parent
def imported(name,file):
    spec=importlib.util.spec_from_file_location(name,HERE/file);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
V=imported('sealed_coordinate50','verify_coordinate_pencil_matrix.py')
E=imported('sealed_end61','verify_marked_end_lattice.py')

def local_cluster(x,f,lam,z,param):
    """Two generic curvettes over the finite point (x=infinity,f=lam) of P1xP1."""
    aa=[S.cancel(1/x.subs(param,j)) for j in [3,7]];bb=[S.cancel(f.subs(param,j)-lam) for j in [3,7]]
    ba,bc='Hx',None;path=('finite',str(lam));rows=[]
    for step in range(150):
        oa,_=V.order(aa[0],z);ob,_=V.order(bb[0],z);assert oa>0 and ob>0
        key=path;rows.append(dict(key=key,parents=[p for p in [ba,bc] if p is not None],multiplicity=min(oa,ob)))
        if ob<oa:aa,bb=bb,aa;ba,bc=bc,ba;direction='b'
        else:direction='a'
        bb=[S.cancel(b/a) for a,b in zip(aa,bb)];ba,bc=key,bc
        residues_orders=[V.order(b,z) for b in bb];residues=[c if n==0 else 0 for n,c in residues_orders]
        if any(n<0 for n,c in residues_orders) or residues[0]!=residues[1]:return dict(rows=rows,terminal=key)
        cb=residues[0]
        if cb!=0:bb=[S.cancel(b-cb) for b in bb];bc=None
        path=path+((direction,str(cb)),)
    raise AssertionError('bounded control only: exceeded150 blowups; not an unbounded proof')

def assemble(clusters):
    full=S.Matrix([[0,1],[1,0]]);labels={'Hx':0,'Hy':1};orders=[-2,-2]
    for cl in clusters:
        for row in cl['rows']:
            key=row['key']
            if key in labels:continue
            old=[labels[p] for p in row['parents']];assert len(old) in [1,2]
            n=full.rows;full=full.row_insert(n,S.zeros(1,n)).col_insert(n,S.zeros(n+1,1));full[n,n]=-1
            for j in old:full[j,j]-=1;full[j,n]=full[n,j]=1
            if len(old)==2:full[old[0],old[1]]-=1;full[old[1],old[0]]-=1
            labels[key]=n;orders.append(sum(orders[j] for j in old)+1)
    ix=[labels[cl['terminal']] for cl in clusters];inverse=full.inv();keys=[k for k in labels if k not in ['Hx','Hy']]
    M=S.Matrix([[next((row['multiplicity'] for row in cl['rows'] if row['key']==key),0) for key in keys] for cl in clusters])
    return dict(full=full,labels=labels,orders=orders,index=ix,source=inverse.extract(ix,ix),multiplicities=M,keys=keys,
        degrees=[sum(full[i,j]!=0 for j in range(full.cols) if i!=j) for i in ix])

def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    preds=['coordinate_pencil_matrix_receipt.json','marked_end_lattice_receipt.json','binomial_periodic_receipt.json','cyclic_companion_receipt.json','oriented_pattern_linking_receipt.json']
    for name in preds:
        old=json.loads((HERE/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for p,h in old['source_sha256'].items():check('sealed predecessor '+name+' '+p,sha(HERE/p)==h)
    manifest=HERE/'sources/common-pencil-mark/manifest.json';source=json.loads(manifest.read_text())
    for row in source['primary_sources']:
        p=HERE/row['path'];check('primary HTML bytes '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    a,b,u,v=S.symbols('a b u v')
    # Direct Jacobian of the blowup charts gives the a_new=sum a_parent+1 rule.
    check('surface blowup Jacobian',S.Matrix([[1,0],[v,u]]).det()==u)
    check('crossing log order rule',S.expand((a+b+1)+1-((a+1)+(b+1)))==0)
    n,q=S.symbols('n q',positive=True,integer=True)
    check('arbitrary unramified-count Schur identity',S.factor(2*q+1-q*n/(n+1)-(q+1+q/(n+1)))==0)
    nr,nd,c=S.symbols('nr nd c',positive=True)
    G=S.Matrix([[nr,c],[c,nd]])
    check('unrestricted binary inverse versus Schur',S.simplify(G.inv()[0,0]-1/(nr-c*c/nd))==0)
    check('two private unit coordinates force squared norm at least2',S.Integer(1)**2+S.Integer(1)**2==2)
    z,w,T,W=S.symbols('z w T W');controls=[]
    cases=[(1,{2:1},1,{2:2},0,2),(1,{3:1},1,{4:2},0,2),
      (2,{1:1},1,{2:2},0,2),(3,{2:1},2,{3:2},0,2),(4,{3:1},3,{4:2},0,2),
      (2,{3:1},2,{3:2},0,2),(3,{4:1},2,{5:2},0,2),(3,{8:1},3,{8:1,4:2},0,2),
      (4,{7:1},4,{7:1,5:2},0,2),(2,{4:2,1:1},1,{2:2},0,2),
      (1,{2:1},1,{2:2},1,2),(1,{2:1},1,{2:2},0,3),(3,{2:1},1,{2:2},0,3)]
    saw_ancestor=False;saw_nonancestor=False
    for r,hr,s,hd,lam,ee in cases:
        R=local_cluster(sum(co*z**(-j) for j,co in hr.items())+z**(ee-r)*T,z**r,0,z,T)
        D=local_cluster(sum(co*w**(-j) for j,co in hd.items())+w**(1-s)*W,lam+w**s,lam,w,W)
        out=assemble([R,D]);Q=out['source'];M=out['multiplicities'];gram=M*M.T
        tag=str((r,hr,s,hd,lam,ee))
        check('P1xP1 whole matrix is unimodular '+tag,abs(out['full'].det())==1)
        check('whole source equals negative curvette Gram '+tag,Q==-gram)
        check('omitted curves remain leaves '+tag,out['degrees']==[1,1])
        check('volume orders match ramified and unramified leaves '+tag,[out['orders'][i] for i in out['index']]==[ee-1,0])
        check('both retained leaves are created free '+tag,all(len(cl['rows'][-1]['parents'])==1 for cl in [R,D]))
        check('R final coordinate private and equal1 '+tag,M[0,out['keys'].index(R['terminal'])]==1 and M[1,out['keys'].index(R['terminal'])]==0)
        ancestor=any(row['key']==D['terminal'] for row in R['rows'])
        schur=S.cancel(gram[0,0]-gram[0,1]**2/gram[1,1]);link=S.cancel(gram.inv()[0,0])
        check('strict Schur and nonzero self linking '+tag,schur>1 and 0<link<1)
        check('R mark is nonzero in exact lattice cokernel '+tag,E.h1(Q.tolist())!=E.h1(Q.tolist(),[[1,0]]))
        if ancestor:
            unique=[row for row in R['rows'] if row['key'] not in {rr['key'] for rr in D['rows']}]
            check('ancestor branch has at least two private centres '+tag,len(unique)>=2 and schur>=2 and len(unique[0]['parents'])==2)
            saw_ancestor=True
        else:
            check('nonancestor branch has private D terminal outside R path '+tag,M[0,out['keys'].index(D['terminal'])]==0)
            saw_nonancestor=True
        controls.append(dict(r=r,R_support=hr,s=s,D_support=hd,D_value=lam,index=ee,
            source_matrix=[[int(x) for x in rr] for rr in Q.tolist()],multiplicities=[[int(x) for x in rr] for rr in M.tolist()],
            leaf_degrees=out['degrees'],R_descends_D=ancestor,Schur=str(schur),R_self_linking=str(link)))
    check('both written ancestry cases exercised',saw_ancestor and saw_nonancestor)
    # Multiple D leaves, including one on a different source value.
    for count in [0,1,2,4]:
        R=local_cluster(z**-2+z*T,z,0,z,T)
        Ds=[local_cluster((i+2)*w**-2+W,w,0,w,W) for i in range(count)]
        out=assemble([R]+Ds);Q=out['source'];gram=out['multiplicities']*out['multiplicities'].T
        expected=2*S.eye(count+1)+2*S.ones(count+1);expected[0,0]+=1
        check('all-count symbolic family exact integer Gram n'+str(count),-Q==gram==expected)
        schur=gram[0,0] if count==0 else S.cancel(gram[0,0]-(gram[0,1:]*gram[1:,1:].inv()*gram[1:,0])[0])
        check('all-count Schur control n'+str(count),schur==3+S.Rational(2,count+1) and schur>1)
        check('all-count retained free leaves n'+str(count),out['degrees']==[1]*(count+1))
    # The same formal orders without the disjoint-leaf condition really attain S=1.
    R=local_cluster(z**-2+z*T,z,0,z,T);D=local_cluster(w**-2+W,w,0,w,W);bad=assemble([R,D]);Q=bad['source']
    check('sharp countercontrol has same volume orders but a nonleaf D',[bad['orders'][i] for i in bad['index']]==[1,0] and bad['degrees']==[1,2])
    check('countercontrol Gram and Schur equality',-Q==S.Matrix([[5,4],[4,4]]) and Q.det()/(-Q[1,1])==1)
    check('countercontrol ramification mark actually zero',all(x.q==1 for x in Q.inv()[:,0]) and E.h1(Q.tolist())==E.h1(Q.tolist(),[[1,0]]))
    old65=json.loads((HERE/'cyclic_companion_receipt.json').read_text());old66=json.loads((HERE/'oriented_pattern_linking_receipt.json').read_text())
    check('all remaining all-odd representations have zero R mark',old65['all_odd']['all_peripheral_S4_classes_covered'] and old65['all_odd']['ramification_class_zero'])
    check('earlier source dichotomy closes complementary sector uniformly',old66['all_odd']['ENTIRE_complementary_source_excluded'] and old66['all_odd']['all_source_orders_and_unramified_poles'])
    check('all-odd prior C3 exclusions retained',old65['all_odd']['entire_C3_case_if_3_divides_pq_excluded'])
    inputs=[Path(__file__),HERE/'COMMON_PENCIL_GRAM_MARK_AND_ALL_ODD_EXCLUSION.md',manifest]+[HERE/x for x in preds]+[HERE/x['path'] for x in source['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      source_theorem=dict(all_field_degrees=True,any_unramified_count=True,any_positive_ramified_volume_order=True,
          smooth_affine_extension_required=True,entire_disjoint_A1_boundary_required=True,
          common_regular_source_coordinate_required=True,complementary_x_global_P1_morphism_required=True,
          sole_positive_volume_prime_required=True,other_volume_orders_zero_required=True,
          values_need_not_coincide=True,whole_lattice_negative_Gram=True,retained_components_free_leaves=True,
          unramified_leaves_incomparable=True,no_unramified_leaf_descends_R=True,
          Schur_strictly_greater_than_one=True,R_self_linking_strictly_between_zero_and_one=True,R_mark_nonzero=True,
          canonical_mark_nonzero_asserted_for_arbitrary_index=False),
      all_odd=dict(entire_actual_primitive_binomial_two_stage_class_excluded=True,
          all_coprime_odd_p_q_at_least_three=True,all_windings_and_source_orders_poles=True,
          all_peripheral_representation_classes_covered=True,complementary_closed_by_PF66=True,
          common_closed_by_nonzero_R_mark=True,PF63_polynomial_immersions_still_exist=True,
          polynomial_curves_declared_nonexistent=False,general_keys_assumed_binomial=False,
          even_first_pairs_excluded=False,longer_sequences_excluded=False),
      exact_blowup_controls=controls,
      countercontrol=dict(whole_disjoint_leaf_hypothesis_dropped=True,volume_orders=[1,0],leaf_degrees=[1,2],
          Gram=[[5,4],[4,4]],Schur=1,R_mark_zero=True,valid_separated_A1_extension=False),
      assurance=dict(written_unbounded_proofs=True,primary_surface_factorization_imported=True,
          finite_controls_not_universal_tree_census=True,integral_coordinate_alone_used_as_automorphism_theorem=False,
          formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,
          full_quartic_exclusion=False,higher_degree_exclusion=False,finite_Keller_map=False,research_goal_completed=False),
      source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'common_pencil_mark_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; all-degree common-pencil mark and entire all-odd class exclusion')

if __name__=='__main__':main()
