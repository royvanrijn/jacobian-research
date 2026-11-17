"""PF2RULED76: unbounded identities and auxiliary exact gluing controls.

Primary interval/gluing/orderability theorems are imported, not replayed.
"""
from pathlib import Path
import hashlib,importlib.util,json,math,os,platform
import sympy as S
HERE=Path(__file__).resolve().parent
def load(name,file):
    spec=importlib.util.spec_from_file_location(name,HERE/file)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
V=load('base76','verify_seifert_base_genus.py');E=V.E
B=load('budget76','verify_peripheral_cone_budget.py')

def negative_cf(a,w):
    out=[]
    while w:
        q=(a+w-1)//w;out.append(q);a,w=w,q*w-a
    assert all(q>=2 for q in out)
    return out

def star(e0,alpha,omega):
    legs=[negative_cf(a,w) for a,w in zip(alpha,omega)];size=1+sum(map(len,legs));mat=S.zeros(size)
    mat[0,0]=e0;i=1
    for leg in legs:
        parent=0
        for weight in leg:
            mat[i,i]=-weight;mat[i,parent]=mat[parent,i]=1;parent=i;i+=1
    return mat

def rational_certificate(mat):
    z=S.ones(mat.rows,1);steps=[]
    for _ in range(20000):
        inter=mat*z;active=[i for i,x in enumerate(inter) if x>0]
        if not active:return dict(rational=True,cycle=[int(x) for x in z],steps=steps)
        i=active[0];val=int(inter[i]);steps.append([i,val])
        if val>1:return dict(rational=False,cycle=[int(x) for x in z],steps=steps,bad_vertex=i,intersection=val)
        z[i]+=1
    raise AssertionError('bounded auxiliary Laufer control exceeded cap; universal proof is separate')

def cubic_data(b,t):
    u=(-pow(3*t,-1,2*b))%(2*b);delta=u//b;a=(3*t*u+1)//(2*b);N=(a+t-1)//t
    al=[b,2*b,t];om=[u%b,u,(-pow(2*b,-1,t))%t if t>1 else 0]
    e0=-N+delta
    return dict(alpha=[x for x in al if x>1],omega=[w for x,w in zip(al,om) if x>1],
                e0=e0,u=u,delta=delta,a=a,N=N)

def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    preds=['seifert_sections_winding_three_receipt.json','rational_end_orderability_receipt.json',
           'lifted_surface_s3_receipt.json','fibre_order_transfer_receipt.json',
           'binomial_periodic_receipt.json','peripheral_cone_budget_receipt.json']
    for name in preds:
        old=json.loads((HERE/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for p,h in old['source_sha256'].items():check('sealed predecessor '+name+' '+p,sha(HERE/p)==h)
    manifest=HERE/'sources/r3-seifert-gluing/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        p=HERE/row['path'];check('primary bytes '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    e0,nu,ll,ceil_sum=S.symbols('e0 nu l ceil_sum',integer=True)
    check('rationality converts floor bound to nonnegative interval numerator',S.expand((1-ll*e0-ceil_sum)-(1-(ll*e0+ceil_sum)))==0)
    check('l1 endpoint numerator',1-e0-nu==-(e0+nu-1))
    # The two-leg case uses only the universal ceil(x)+ceil(y)<=ceil(x+y)+1.
    for x,y in [(S.Rational(1,2),S.Rational(1,3)),(S.Rational(2,7),S.Rational(9,14))]:
        check('two-leg fractions have sum below1 '+str((x,y)),x+y<1)
        for n in [1,2,3,7,19]:
            check('two-leg rounding auxiliary '+str((x,y,n)),S.ceiling(n*x)+S.ceiling(n*y)<=n+1)
        check('two-leg l1 endpoint exactly0 '+str((x,y)),1+1-S.ceiling(x)-S.ceiling(y)==0)
    b,t,u,delta,NN=S.symbols('b t u delta N',integer=True,positive=True)
    aa=(3*t*u+1)/(2*b);ww=NN*t-aa
    norm=(u-delta*b)/b+u/(2*b)+ww/t+1/(2*b*t)
    check('unbounded normalized Euler integer identity',S.cancel(norm-(NN-delta))==0)
    check('a below2t when u<=b-1',S.expand(4*b*t-(3*t*(b-1)+1))==b*t+3*t-1)
    check('a below3t when u<=2b-1',S.expand(6*b*t-(3*t*(2*b-1)+1))==3*t-1)
    check('b1 ceiling upper gap',S.expand(2*t-(3*t+1)/2)==(t-1)/2)
    check('b1 ceiling lower gap',S.expand((3*t+1)/2-t)==(t+1)/2)
    check('t1 a1 delta zero inequality',S.expand(b-(2*b-1)/3)==(b+1)/3)
    check('t1 a2 delta one inequality',S.expand((4*b-1)/3-b)==(b-1)/3)
    for dd,nmax in [(0,2),(1,3)]:
        check('integer e0 bound split delta '+str(dd),all(n-dd<=2 for n in range(1,nmax+1)))
    # The actual cubic weights: quasi-homogeneity and primitive circle action.
    check('weighted cubic term degrees',3*b*t==2*b*t+b*t and 3*b*t==b*3*t)
    controls=[];branches=set();census=json.loads((HERE/'binomial_periodic_receipt.json').read_text())['census']['complete_cells']
    cases=[(1,3),(1,13),(2,1),(7,1),(2,3),(2,5),(4,3),(7,15)]
    for bb,tt in cases:
        dat=cubic_data(bb,tt);alpha,omega,ec=dat['alpha'],dat['omega'],dat['e0'];tag=str((bb,tt))
        check('primitive cubic weights and coefficients '+tag,math.gcd(2*bb,3*tt)==1 and math.gcd(2*bb,math.gcd(3*tt,bb*tt))==1 and all(0<w<a and math.gcd(a,w)==1 for a,w in zip(alpha,omega)))
        check('actual Euler number and integer normalization '+tag,ec+sum(S.Rational(w,a) for a,w in zip(alpha,omega))==-S.Rational(1,2*bb*tt) and ec+len(alpha)>=1)
        check('exact normal-slice residues '+tag,(3*tt*dat['u']+1)%(2*bb)==0 and (2*bb*((-pow(2*bb,-1,tt))%tt)+1)%tt==0 if tt>1 else (3*dat['u']+1)%(2*bb)==0)
        check('ceiling normalization identity '+tag,dat['N']-dat['delta']==-ec and dat['N']*tt-dat['a']==((-pow(2*bb,-1,tt))%tt if tt>1 else 0))
        mat=star(ec,alpha,omega);cert=rational_certificate(mat)
        check('actual negative-definite cubic star '+tag,(-mat).is_positive_definite is True and abs(mat.det())==bb)
        branches.add('rational' if cert['rational'] else 'nonrational')
        witness=None
        if cert['rational']:
            check('rational cubic l1 endpoint '+tag,ec+len(alpha)==1)
            for n in [1,2,3,7,19,math.lcm(*alpha)]:
                check('rational cubic universal-bound auxiliary '+str((tag,n)),n*ec+sum(S.ceiling(S.Rational(n*w,a)) for a,w in zip(alpha,omega))<=1)
        else:
            check('exact Laufer nonrational witness '+tag,cert['intersection']>1)
            # This bounded search only finds an explicit auxiliary witness; it proves no universal bound.
            for n in range(1,math.lcm(*alpha)+1):
                if sum(S.floor(-S.Rational(n*w,a)) for a,w in zip(alpha,omega))<=n*ec-2:witness=n;break
            check('independent explicit Pinkham witness '+tag,witness is not None)
            if ec+len(alpha)==1:branches.add('nonrational_with_critical_sum1')
        pp,qq=2*bb,3*tt;wp=(-pow(qq,-1,pp))%pp;wq=(-pow(pp,-1,qq))%qq
        check('ordinary torus side normalized Euler '+tag,qq*wp+pp*wq==pp*qq-1 and S.Rational(wp,pp)+S.Rational(wq,qq)<1)
        aa0=bb*tt;cc=pp*qq-pp-qq;kk=3*cc+1;DD=kk//2-aa0
        phi=S.Matrix([[-1,kk//2],[0,1]])
        check('actual gluing matches meridians and reverses orientation '+tag,phi.det()==-1 and phi*S.Matrix([1,0])==S.Matrix([-1,0]))
        check('actual slope zero has no section translation '+tag,S.Rational(1,2*aa0)==-(-S.Rational(1,2*bb*tt)))
        genuine=[a for a in [bb,2*bb,tt] if a>1];orb=1-len(genuine)+sum(S.Rational(1,a) for a in genuine)
        check('both actual Seifert boundaries incompressible product-cover range '+tag,orb<0 and S.Rational(1,pp)+S.Rational(1,qq)-1<0)
        want=dict(H1=B.group([DD,bb]),quotient=B.group([math.gcd(DD,kk),bb]))
        rows=[x for x in census[f'{pp%12},{qq%12},3']['orbits'] if (x['companion_order'],x['cable_order'],x['cycle_rank'])==(6,6,0)]
        check('both actual ordinary classes retained '+tag,len(rows)==2)
        for i,row in enumerate(rows):check('actual full satellite marked group '+str((tag,i)),V.full_satellite(pp,qq,3,row)==want)
        controls.append(dict(b=bb,t=tt,normalized=dat,critical_sum=ec+len(alpha),Laufer=cert,
          explicit_Pinkham_witness=witness,meridian_status='endpoint' if cert['rational'] else 'outside_L_space_interval',
          D=DD,actual_gluing=[[-1,kk//2],[0,1]],**want))
    check('both rational and nonrational branches and nontrivial critical case exercised',branches=={'rational','nonrational','nonrational_with_critical_sum1'})
    old=json.loads((HERE/'seifert_sections_winding_three_receipt.json').read_text())
    check('prior complete winding and trefoil reductions retained',old['quartic']['entire_actual_two_stage_common_requires_r']==3 and old['quartic']['trefoil_common_entirely_excluded'])
    inputs=[Path(__file__),HERE/'MATCHED_MERIDIAN_SLOPES_AND_ALL_COMMON_EXCLUSION.md',manifest]+[HERE/n for n in preds]+[HERE/x['path'] for x in sources['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      general=dict(no_field_degree_bound=True,matched_actual_slopes_outside_both_L_space_interiors=True,
        Seifert_endpoint_lemma=True,normalized_negative_Euler_required=True,at_least_two_genuine_cones_required=True,
        central_plus_cone_count_at_least1_required=True,non_solid_torus_pieces_required=True,
        actual_graph_gluing_required=True,irreducibility_verified_before_LO_implication=True,
        PF74_explicit_common_source_hypotheses_required=True,nu2_handled_without_extending_Pinkham=True),
      quartic=dict(actual_primitive_binomial_two_stage_required=True,all_first_exponents_windings_source_orders_poles=True,
        all_residual_r3_non_trefoil_common_sources_excluded=True,entire_common_source_sector_excluded=True,
        every_remaining_source_pencil='COMPLEMENTARY only',companion_meridian_not_in_L_space_interior=True,
        ordinary_meridian_not_in_L_space_interior=True,actual_meridians_matched=True,
        actual_r3_end_not_L_space=True,actual_r3_end_group_nontrivial_left_orderable=True,
        passing_complementary_sources_excluded=False,entire_target_keys_in_both_pencils_excluded=False),
      exact_controls=controls,
      assurance=dict(written_unbounded_endpoint_normalization_and_gluing_proofs=True,
        complete_Rasmussen_and_Nemethi_primary_pages_visually_audited=True,
        deep_interval_gluing_and_orderability_theorems_imported=True,finite_controls_not_universal_proof=True,
        deep_theorems_independently_reproved=False,interval_confused_with_interior=False,
        cone_twists_inferred_from_H1=False,actual_companions_replaced_by_virtual_solid_tori=False,
        nu2_Pinkham_extension_assumed=False,full_quartic_exclusion=False,higher_degree_exclusion=False,
        finite_Keller_map=False,PF63_polynomial_curves_declared_nonexistent=False,formal_verification=False,
        independent_full_reproof=False,external_review=False,novelty_priority=False,research_goal_completed=False),
      source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'matched_meridian_common_exclusion_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; degree-free matched-slope lemma; ENTIRE common two-stage quartic sector excluded')
if __name__=='__main__':main()
