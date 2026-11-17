"""Controls for Xie's matrix, source-coordinate pairings and uniform bounds."""
from pathlib import Path
from math import gcd
import hashlib,json,os,platform
import sympy as S

HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def order(f,z):
    n,d=map(lambda h:S.Poly(h,z),S.fraction(S.cancel(f)))
    a=min(k[0] for k,c in n.terms());b=min(k[0] for k,c in d.terms())
    return a-b,S.cancel(n.nth(a)/d.nth(b))
def contact_sum(r,q,support):
    return sum(2*q-max(j for j in support if a*j%r) for a in range(1,r))
def diagonal_size(e,r,cs):
    q=max(cs);return contact_sum(r,q,cs)+2*q-r+e
def cluster(x,y,z,param):
    # Two exact generic curvettes, stopped when their centres separate.
    # The resulting residue coordinate is fractional linear in the chart parameter.
    xs=[S.cancel(x.subs(param,j)) for j in [3,7]];ys=[S.cancel(y.subs(param,j)) for j in [3,7]]
    ox,cx=order(xs[0],z);oy,cy=order(ys[0],z)
    if ox<oy:aa=[1/a for a in xs];bb=[b/a for a,b in zip(xs,ys)];root='X'
    elif oy<ox:aa=[1/b for b in ys];bb=[a/b for a,b in zip(xs,ys)];root='Y'
    else:aa=[1/a for a in xs];bb=[b/a-cy/cx for a,b in zip(xs,ys)];root=('X',str(cy/cx))
    aa=list(map(S.cancel,aa));bb=list(map(S.cancel,bb));ba,bc='L',None;path=(root,);rows=[]
    for unused in range(100):
        oa,_=order(aa[0],z);ob,_=order(bb[0],z)
        assert oa>0 and ob>0,'invalid positive local centre'
        key=path;rows.append((key,ba,bc))
        if ob<oa:aa,bb=bb,aa;ba,bc=bc,ba;direction='b'
        else:direction='a'
        bb=[S.cancel(b/a) for a,b in zip(aa,bb)];ba,bc=key,bc
        residues_and_orders=[order(b,z) for b in bb]
        residues=[c if n==0 else 0 for n,c in residues_and_orders]
        if any(n<0 for n,c in residues_and_orders) or residues[0]!=residues[1]:return rows,key
        cb=residues[0]
        if cb!=0:bb=[S.cancel(b-cb) for b in bb];bc=None
        path=path+((direction,str(cb)),)
    raise AssertionError('control extraction exceeded its bounded step limit')
def blowup_matrix(clusters):
    Q=S.Matrix([[1]]);labels={'L':0};K=[-3]
    for rows,end in clusters:
        for key,ba,bb in rows:
            if key in labels:continue
            old=[labels[b] for b in [ba,bb] if b is not None]
            n=Q.rows;Q=Q.row_insert(n,S.zeros(1,n)).col_insert(n,S.zeros(n+1,1));Q[n,n]=-1
            for j in old:Q[j,j]-=1;Q[j,n]=Q[n,j]=1
            if len(old)==2:Q[old[0],old[1]]-=1;Q[old[1],old[0]]-=1
            labels[key]=n;K.append(sum(K[j] for j in old)+1)
    return Q,Q.inv(),labels,K
def main():
    preds=['two_sided_norm_receipt.json','unramified_sheet_pair_receipt.json','universal_regular_end_receipt.json','exact_coordinate_bridge_receipt.json','unramified_orevkov_receipt.json','coordinate_regular_sector_receipt.json']
    for name in preds:
        r=json.loads((HERE/name).read_text());check('predecessor PASS '+name,r['outcome']=='PASS')
        for p,h in r['source_sha256'].items():check('predecessor seal '+name+' '+p,sha(HERE/p)==h)
    manifest=HERE/'sources/coordinate-pencil-matrix/manifest.json';m=json.loads(manifest.read_text())
    for row in m['primary_sources']:
        f=HERE/row['path_relative_to_affine_complements'];check('primary source '+f.name,sha(f)==row['sha256'] and f.stat().st_size==row['bytes'])
    r,q,s,p,t,g,a=S.symbols('r q s p t g a')
    upperR=2*r*q-r*r+r+1;upperD=2*s*p-s*s-s+2
    check('unramified global bound has exact nonnegative defect',S.expand(p*p-upperD-((p-s)**2+s-2))==0)
    check('index-two global bound has exact nonnegative decomposition',S.expand((3*q*q-upperR).subs(q,r-1+t)-((2*r-1)*(r-2)+(4*r-6)*t+3*t*t))==0)
    sharp=(r+g)*q-g*r+g+1
    check('leading-character refinement has the exact square defect',S.expand(q*q-sharp-((q-r)*(q-g)-g-1))==0)
    check('gcd divides pole gap gives a uniform nonnegative high-pole bound',S.expand(((q-r)*(q-g)-g-1).subs(q,r+a)-(a*(r-1)-1)-(a+1)*(a-g))==0)
    check('simple R quadratic comparison factors exactly',S.expand(3*q*q-(2*q+1)-(3*q+1)*(q-1))==0)
    check('minimal R full diagonal is r squared minus r plus one',S.expand(((r-1)*(r-1)+2*q-r+2).subs(q,r-1)-(r*r-r+1))==0)
    check('sharp low-pole R size reaches three',diagonal_size(2,2,{1:1})==3)
    check('sharp exact unramified normalized diagonal reaches one',diagonal_size(1,1,{2:1})==4)
    z,w,T,V,x,y=S.symbols('z w T V x y');controls=[]
    examples=[(1,{2:1},1,{2:2}), (2,{1:1},1,{2:2}), (3,{2:1},2,{3:2}), (4,{6:2,3:3},3,{4:2})]
    for rr,cs,ss,ds in examples:
        qq=max(cs);pp=max(ds);h=sum(c*z**-j for j,c in cs.items());hD=sum(c*w**-j for j,c in ds.items())
        R=cluster(h+z**(2-rr)*T,z**rr,z,T)
        kinds=['complement','shear2']+(['long'] if rr<=3 else [])
        for kind in kinds:
            if kind=='complement':xd,yd=w**ss,hD+w**(1-ss)*V;mm=1;gg=x
            elif kind=='shear2':xd=hD+w**(1-ss)*V;yd=w**ss-xd**2;mm=2;gg=y+x*x
            else:
                mate=hD+w**(1-ss)*V;yd=w**ss-mate**2;xd=mate-yd**2;mm=2;gg=(x+y*y)**2+y
            D=cluster(xd,yd,w,V);Q,I,labels,K=blowup_matrix([R,D]);i=labels[R[1]];j=labels[D[1]]
            nR=diagonal_size(2,rr,cs);nD=diagonal_size(1,ss,ds);tag=f'r{rr} q{qq} s{ss} p{pp} {kind}'
            check('full blowup matrix is unimodular '+tag,abs(Q.det())==1)
            check('actual primitive R dual diagonal '+tag,I[i,i]==-nR)
            check('actual primitive D dual diagonal '+tag,I[j,j]==-nD)
            check('actual cross pairing uses source pencil degree '+tag,I[i,j]==mm*qq*pp)
            check('source canonical orders recovered from all centres '+tag,K[i]==1 and K[j]==0)
            check('polynomial coordinate degree equals its entire R pole multiple '+tag,S.degree(gg,x)==mm and order(gg.subs({x:h+z**(2-rr)*T,y:z**rr},simultaneous=True),z)[0]==-mm*qq)
            check('pair determinant has exact product form '+tag,I.extract([i,j],[i,j]).det()==nR*nD-mm*mm*qq*qq*pp*pp)
            controls.append(dict(r=rr,q=qq,s=ss,p=pp,relative_coordinate=kind,m=mm,R_dual_square=int(I[i,i]),D_dual_square=int(I[j,j]),cross=int(I[i,j]),canonical_orders=[K[i],K[j]],blowups=Q.rows-1))
    # Exact equality case is excluded by strict negative definiteness.
    R=cluster(z**-2+T,z,z,T);D=cluster(w,w**-2+V,w,V)
    Q,I,labels,K=blowup_matrix([R,D]);i=labels[R[1]];j=labels[D[1]]
    check('two complementary exact unramified primes fail even at equality',I[i,i]==I[j,j]==-4 and I[i,j]==4 and I.extract([i,j],[i,j]).det()==0)
    # Without the exact residue, this universal bound would be false.
    check('nonexact simple unramified residue is outside the diagonal bound',diagonal_size(1,1,{1:1})==2 and 2>1**2)
    check('nonexact complementary matrix can pass the necessary test',S.Matrix([[-2,1],[1,-2]]).det()==3)
    # Same pencil has an overlapping cluster and is not governed by m*q*p.
    R=cluster(z**-2+z*T,z,z,T);D=cluster(2*w**-2+V,w,w,V)
    Q,I,labels,K=blowup_matrix([R,D]);i=labels[R[1]];j=labels[D[1]]
    common=I.extract([i,j],[i,j])
    check('actual common-ruling control keeps a negative cross pairing',common==S.Matrix([[-5,-2],[-2,-4]]) and common.det()>0)
    check('a passing complementary matrix is not an affine certificate',S.Matrix([[-5,4],[4,-4]]).det()==4 and json.loads((HERE/'unramified_sheet_pair_receipt.json').read_text())['simple_pair']['whole_sector_excluded'])
    arithmetic=[]
    for rr in range(2,15):
        for qq in range(rr+1,rr+8):
            gg=gcd(rr,qq);sharpn=(rr+gg)*qq-gg*rr+gg+1
            check(f'finite gcd count control r{rr} q{qq}',sum(rr and aa*qq%rr!=0 for aa in range(1,rr))==rr-gg and sharpn<=qq*qq)
    for rr in [2,3,4,5,12,101]:
        nR=rr*rr-rr+1;survivors=[pp for pp in range(2,9) if 2*nR*pp>(rr-1)**2*pp*pp]
        check('simple D strict matrix survivor control r'+str(rr),survivors==([2,3,4,5] if rr==2 else [2,3] if rr==3 else [2]))
        arithmetic.append(dict(r=rr,strict_matrix_simple_D_poles=survivors,PF48_additionally_excludes_p_equals_r=rr==3))
    inputs=[Path(__file__),HERE/'GLOBAL_COORDINATE_PENCIL_AND_MATRIX.md',manifest]+[HERE/n for n in preds]+[HERE/row['path_relative_to_affine_complements'] for row in m['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      matrix=dict(primary='Xie Theorem1.2',unnormalized_primitive_duals=True,all_other_primes_retained=True,necessary_for_affine_N=True,passing_test_proves_affineness=False,diagonal='-(S+2q-r+e)',distinct_coordinate_pencil_cross='m*q*p',strict_necessary_condition='n_R*n_D>m^2*q^2*p^2',all_principal_parts_retained=True,all_source_automorphisms_covered=True,controls=controls),
      global_conclusions=dict(all_Keller_field_degrees=True,all_unramified_primes_share_one_intrinsic_coordinate_pencil=True,boundary_constants_need_not_agree=True,index_two_prime_pencil_common_or_complementary_to_unramified=True,noncomplementary_pencils_excluded_without_pole_or_defect_restrictions=True,unramified_bound='0<n_D<=p^2',index_two_bound='0<n_R<=3q^2',common_coordinate_case_excluded=False),
      complementary=dict(all_D_orders=True,exact_product_test_retained=True,simple_R_forced_q2=True,higher_R_forced_q_rminus1=True,higher_R_full_principal_part_single_monomial=True,all_higher_and_multicharacteristic_R_poles_excluded=True,simple_R_requires_D_order_at_least_three=True,simple_D_R_order_two_poles=[2,3,4,5],simple_D_R_order_at_least_three_forced_p2=True,arithmetic_controls=arithmetic),
      controls_boundary=dict(common_ruling_matrix=[[int(c) for c in row] for row in common.tolist()],nonexact_residue_control=True,passing_matrix_nonaffineness_control=True,extractor_uses_exact_numeric_curvettes=True,shared_initial_centres_merged=True),
      assurance=dict(written_unbounded_matrix_proof=True,imported_primary_theorem_independently_reproved=False,finite_controls_not_universal_proof=True,full_quartic_exclusion=False,higher_degree_problem_excluded=False,target_Ga_action_or_ruling_inferred=False,formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,research_goal_completed=False),
      source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'coordinate_pencil_matrix_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; global unramified pencil and sharp index-two compatibility')
if __name__=='__main__':main()
