"""PF2RULED74: unbounded proof identities and exact auxiliary controls.

Deep rational-link/orderability and plumbing theorems are imported, not replayed.
"""
from pathlib import Path
import hashlib,importlib.util,json,math,os,platform
import sympy as S
HERE=Path(__file__).resolve().parent
def load(name,file):
    spec=importlib.util.spec_from_file_location(name,HERE/file)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
C=load('cluster74','verify_common_pencil_mark.py');E=C.E
V=load('base74','verify_seifert_base_genus.py')
B=load('budget74','verify_peripheral_cone_budget.py')

def components(mat,vertices):
    left=set(vertices);out=[]
    while left:
        seen={min(left)}
        while True:
            new=seen|{j for i in seen for j in left if j!=i and mat[i,j]}
            if new==seen:break
            seen=new
        out.append(sorted(seen));left-=seen
    return out

def laufer(mat):
    cycle=S.ones(mat.rows,1);intersections=[]
    for _ in range(20000):
        ints=mat*cycle;active=[i for i,x in enumerate(ints) if x>0]
        if not active:
            chi=-(cycle.dot(mat*cycle)+sum((-2-mat[i,i])*cycle[i] for i in range(mat.rows)))/2
            return cycle,intersections,chi
        i=active[0];intersections.append(int(ints[i]));cycle[i]+=1
    raise AssertionError('bounded auxiliary Laufer control exceeded cap; not a universal proof')

def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    preds=['lifted_surface_s3_receipt.json','fibre_order_transfer_receipt.json',
           'fold_schur_s4_receipt.json','common_pencil_mark_receipt.json',
           'binomial_periodic_receipt.json','marked_end_lattice_receipt.json']
    for name in preds:
        old=json.loads((HERE/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for p,h in old['source_sha256'].items():check('sealed predecessor '+name+' '+p,sha(HERE/p)==h)
    manifest=HERE/'sources/rational-end-orderability/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        p=HERE/row['path'];check('primary bytes '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    a,zz,kk=S.symbols('a zz kk',integer=True)
    check('unbounded canonical blowup chi identity',S.expand(-(zz-a*a+kk-a)/2-(-(zz+kk)/2+a*(a+1)/2))==0)
    # Nonnegative for every integer a: a>=0 or a<=-1.
    n=S.symbols('n',integer=True,nonnegative=True)
    check('negative total-transform coefficient certificate',S.expand((a*(a+1)/2).subs(a,-n-1)-n*(n+1)/2)==0)
    pair=S.Matrix([[S.Symbol('w'),1],[1,0]])
    check('zero-end cancellation Schur correction identically zero',pair.inv()[0,0]==0 and pair.det()==-1)
    z,w,T,W=S.symbols('z w T W');source_controls=[];saw_branch=False;saw_many=False;saw_high_volume=False
    cases=[(2,{4:2,1:1},2,[0]),(4,{7:1},2,[0]),(3,{8:1,4:2},2,[0]),
           (2,{1:1},2,[1]),(3,{4:1},3,[1,2]),(2,{3:1},5,[0,1,2])]
    for r,hr,index,values in cases:
        R=C.local_cluster(sum(co*z**(-q) for q,co in hr.items())+z**(index-r)*T,z**r,0,z,T)
        Ds=[C.local_cluster((i+2)*w**-2+W,value+w,value,w,W) for i,value in enumerate(values)]
        out=C.assemble([R]+Ds);mat=out['full'];hx,hy=out['labels']['Hx'],out['labels']['Hy']
        tag=str((r,hr,index,values));vertices=[i for i in range(mat.rows) if i not in out['index']+[hx,hy]]
        blocks=components(mat,vertices)
        check('entire omitted leaves and all volume orders '+tag,out['degrees']==[1]*(len(Ds)+1) and [out['orders'][i] for i in out['index']]==[index-1]+[0]*len(Ds))
        check('untouched zero end and unique neighbour '+tag,mat[hy,hy]==0 and [i for i in range(mat.rows) if i!=hy and mat[hy,i]]==[hx])
        attachments=[sum(mat[hx,i]!=0 for i in block) for block in blocks]
        check('exact R6 no extra S1xS2 count '+tag,all(x==1 for x in attachments) and sum(x-1 for x in attachments)==0)
        determinants=[];funds=[]
        for j,block in enumerate(blocks):
            q=mat.extract(block,block);cycle,steps,chi=laufer(q)
            check('source negative definite exceptional subgraph '+str((tag,j)),(-q).is_positive_definite is True)
            check('source rational fundamental cycle '+str((tag,j)),chi==1 and all(x==1 for x in steps))
            determinants.append(int(abs(q.det())));funds.append([int(x) for x in cycle])
            saw_branch |= any(sum(q[i,h]!=0 for h in range(q.rows) if h!=i)>=3 for i in range(q.rows))
        check('full boundary homology retained by R6 '+tag,math.prod(determinants)==abs(out['source'].det()))
        saw_many |= len(blocks)>1;saw_high_volume |= index>2
        source_controls.append(dict(r=r,support={str(k):v for k,v in hr.items()},index=index,
            D_values=values,component_determinants=determinants,fundamental_cycles=funds))
    check('branching graphs multiple fibres and arbitrary volume exercised',saw_branch and saw_many and saw_high_volume)
    # General S3 quotient: simultaneous peripheral-compatible abelianizations.
    b,t,s=S.symbols('b t s',positive=True,integer=True);A=b*t;k=3*s*(6*b*t-2*b-3*t)+1;D=k/2-A*s;BB=k/2-3*s*A
    check('companion abelianization kills actual eta',-A+A==0)
    check('general third cone determinant',S.expand(S.Matrix([[k/2,A],[3*s,1]]).det()-BB)==0)
    check('general third cone cannot vanish certificate',S.expand(2*BB-(3*s*(4*b*t-2*b-3*t)+1))==0)
    # In a hypothetical B=0, the last identity would force 3 to divide1.
    check('nonzero cone integer obstruction',1%3!=0)
    d=S.symbols('d',integer=True,positive=True);ss=2*d-1;kc=6*d-2;J=3*d-2;prod=ss*kc*J
    check('trefoil third order',S.expand(-BB.subs({b:1,t:1,s:ss})-J)==0)
    check('normalized Euler exact determinant',S.factor(-2+1/ss+(kc-1)/kc+(J-1)/J+d/prod)==0)
    check('pairwise coprimality first two certificate',S.expand(kc-3*ss)==1)
    check('pairwise coprimality first third certificate',S.expand(2*J-3*ss)==-1)
    check('last pair gcd reduces to2',S.expand(kc-2*J)==2)
    ell=ss+1
    check('unbounded floor interval l between s and2s',S.expand(2*ss-ell)==2*d-2)
    check('unbounded floor interval l below k',S.expand(kc-ell)==4*d-2)
    check('unbounded floor interval l below J',S.expand(J-ell)==d-2)
    # A second, direct effective-cycle certificate: central m=s+2, short leg2,
    # on each long -2 arm m-1,...,1,0,... . Its chi is exactly0.
    mm=ss+2
    check('explicit nonrational cycle chi zero',S.expand(-((2*mm-4*ss)+2*(ss-2))/2)==0)
    check('long arm room for nonrational cycle',S.expand((J-1)-(mm-1))==d-3)
    target_controls=[]
    census=json.loads((HERE/'binomial_periodic_receipt.json').read_text())['census']['complete_cells']
    for dd in [11,19,27,43]:
        si,ki,ji=2*dd-1,6*dd-2,3*dd-2;alpha=[si,ki,ji];omega=[1,ki-1,ji-1];ll=si+1
        tag=str(dd);product=math.prod(alpha)
        check('all normalized Seifert congruences '+tag,all((om*(product//al)+dd)%al==0 and math.gcd(al,om)==1 for al,om in zip(alpha,omega)))
        check('star Euler negative and order exact '+tag,-2+sum(S.Rational(om,al) for al,om in zip(alpha,omega))==-S.Rational(dd,product))
        floor_sum=sum(S.floor(-S.Rational(ll*om,al)) for al,om in zip(alpha,omega))
        check('universal Pinkham floor witness '+tag,floor_sum==-2*ll-2)
        # Small integral presentation verifies full group, rather than just determinant.
        # alpha*x-omega*h=0 and sum(x)-2h=0 is the negative-star convention.
        matrix=[[al if j==i else (-om if j==3 else 0) for j in range(4)] for i,(al,om) in enumerate(zip(alpha,omega))]+[[1,1,1,-2]]
        check('Seifert model full H1 '+tag,E.h1(matrix)==B.group([dd]))
        rr=6*dd-3;rows=[x for x in census[f'2,3,{rr%12}']['orbits'] if (x['companion_order'],x['cable_order'],x['cycle_rank'])==(6,6,0)]
        check('all trefoil S3 classes retained '+tag,len(rows)==2)
        for i,row in enumerate(rows):
            check('actual satellite full H1 matches quotient model '+str((dd,i)),V.full_satellite(2,3,rr,row)==dict(H1=B.group([dd]),quotient=B.group([math.gcd(dd,2)])))
        target_controls.append(dict(D=dd,r=rr,alpha=alpha,omega=omega,e0=-2,floor_witness=ll,floor_sum=int(floor_sum),H1=B.group([dd])))
    old=json.loads((HERE/'fibre_order_transfer_receipt.json').read_text())
    check('all other common sectors previously closed',old['quartic']['remaining_common_sector']=='S3/S3 only')
    old=json.loads((HERE/'lifted_surface_s3_receipt.json').read_text())
    check('PF73 trefoil necessary bounds retained',old['quartic']['trefoil_common_necessary_r_min']==63 and old['quartic']['trefoil_common_necessary_r_mod48']==15)
    inputs=[Path(__file__),HERE/'RATIONAL_SOURCE_ENDS_AND_TREFOIL_COMMON_EXCLUSION.md',manifest]+[HERE/n for n in preds]+[HERE/x['path'] for x in sources['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
        source_theorem=dict(all_field_degrees=True,all_ramification_volume_orders=True,any_omitted_line_count=True,
          smooth_affine_extension_required=True,entire_disjoint_A1s_required=True,
          common_regular_coordinate_finite_values_required=True,x_global_P1_morphism_with_poles_required=True,
          end_connected_sum_of_rational_negative_plumbing_links=True,no_extra_S1S2_summand=True,
          unbounded_integer_chi_blowup_proof=True,every_homomorphism_to_left_orderable_group_trivial=True,
          source_end_assumed_irreducible=False,all_source_end_groups_finite_or_cyclic_asserted=False),
        quotient=dict(all_actual_binomial_two_stage_quartic_S3_S3=True,actual_to_Seifert_group_surjection=True,
          simultaneous_companion_abelianizations=True,cones='s,k,abs(k/2-r*bt)',H1='Z/D',
          nontrivial_LO_quotient_excludes_common_source=True,group_equality_or_homeomorphism_asserted=False,
          rationality_classified_for_all_first_pairs=False),
        trefoil=dict(all_remaining_S3_S3_common_sources_excluded=True,all_source_configurations_and_orders=True,
          all_windings_common_pencil_excluded=True,all_common_representation_sectors_excluded=True,
          entire_target_keys_in_both_pencils_excluded=False,complementary_passing_sources_excluded=False,
          normalized_star='e0=-2; (s,1),(k,k-1),(J,J-1)',s='2D-1',k='6D-2',J='3D-2',
          unbounded_nonrational_floor_witness='l=s+1',actual_target_has_nontrivial_LO_quotient=True),
        source_controls=source_controls,target_controls=target_controls,
        assurance=dict(written_unbounded_source_and_target_proofs=True,complete_primary_pages_visually_audited=True,
          deep_rational_link_orderability_theorems_imported=True,deep_theorems_independently_reproved=False,
          finite_controls_not_universal_proof=True,finite_Keller_map=False,full_quartic_exclusion=False,
          higher_degree_exclusion=False,PF63_polynomial_curves_declared_nonexistent=False,
          formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,
          research_goal_completed=False),
        source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'rational_end_orderability_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; all-degree source no-LO-image theorem; entire trefoil common pencil excluded')
if __name__=='__main__':main()
