"""PF115 exact controls; universal result rests on the written all-branch proof."""
from pathlib import Path
import hashlib,json,platform
import sympy as S
import residual_milnor as R
import triple_boundary_residue as T
import coprime_leading_forms as C

P=Path(__file__).resolve().parent


def local_gradient_length(poly,variables,N):
    x,y=variables
    G=S.groebner([S.diff(poly,x),S.diff(poly,y)]+[x**i*y**(N-i) for i in range(N+1)],x,y)
    lm=[p.LM(order=G.order).exponents for p in G.polys]
    length=sum(1 for a in range(N) for b in range(N-a) if not any(a>=i and b>=j for i,j in lm))
    return dict(order=N,length=length,leading_monomials=lm)


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def ident(name,value):check(name,S.factor(value)==0)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    oldfiles=['torus_cusp_excess_receipt.json','cyclic_quotient_morse_receipt.json','triple_boundary_residue_receipt.json','coprime_leading_forms_receipt.json']
    for name in oldfiles:
        d=json.loads((P/name).read_text());check('sealed predecessor PASS '+name,d['outcome']=='PASS')
        for f,h in d['source_sha256'].items():check('sealed predecessor '+name+' '+f,sha(P/f)==h)
    manifest=json.loads((P/'sources/residual-milnor/manifest.json').read_text())
    for item in manifest['sources']:
        f=P/item['path'];check('exact primary/context '+f.name,sha(f)==item['sha256'] and f.stat().st_size==item['bytes'])
    for item,pages in zip(manifest['sources'][:4],[[2,3,4,10],[15,17,18,20,21],[9],[3]]):
        check('complete primary pages fresh '+item['title'],item['fresh_complete_PDF_pages_rendered_and_visually_inspected']==pages)
    check('charzero imports and PF114 credit explicit',manifest['characteristic_zero_required'] and manifest['polar_Milnor_and_JSJ_inputs_imported'] and manifest['all_residue_one_through_degree15_credited_to_PF114'])
    n,q,eta,c,k=S.symbols('n q eta c k',integer=True,positive=True)
    I=n*q*(c-1)+1
    source_polar_sum=n*eta*(c-1)+k-I
    ident('UNBOUNDED all residual branches Milnor summation',source_polar_sum-k+1-n*(eta-q)*(c-1))
    t=S.symbols('t',integer=True,nonnegative=True)
    negative=n*((q-1-t)-q)*(1+t)
    check('UNBOUNDED eta<q makes Milnor number negative',negative.is_negative)
    e=S.symbols('e',integer=True,positive=True)
    ident('UNBOUNDED e>=4 minimalexcess contradicts torus cost',((4+t)*q-(5*q-2*q))-q*(1+t))
    # A torus cover is an invertible rational map on peripheral H1. It cannot
    # identify the two distinct algebraic companion fibre lines.
    pp,aa,pnext,nu=S.symbols('pp aa pnext nu',positive=True)
    alpha=pnext*pp*aa+nu
    v1=S.Matrix([pp*aa,1]);v2=S.Matrix([alpha,pnext])
    ident('UNBOUNDED algebraic cabling fibre mismatch',S.det(S.Matrix.hstack(v1,v2))+nu)
    a,b,d,z=S.symbols('a b d z');M=S.Matrix([[a,b],[d,z]])
    ident('UNBOUNDED torus-cover determinant retains mismatch',S.det(S.Matrix.hstack(M.inv()*v1,M.inv()*v2))-(-nu)/(a*z-b*d))
    f=S.symbols('f',integer=True,positive=True)
    ident('UNBOUNDED actual nonprimitive source budget',(e*f+(f+1)+(k-1)+1)-((e+1)*f+k+1))
    ident('UNBOUNDED per-parameter excess charge',k*(f+1)-(k*f+k))
    check('UNBOUNDED low-mass degree bound dominates PF114 residue floor',(18*(2+t)-19).is_positive)
    certified=dict(axis_source_certified=True,sole_critical_branch_certified=True,chosen_singular_image_branch_certified=True)
    equality=[]
    for ee,qq,nn,ell in [(2,1,1,3),(3,2,1,8),(3,1,2,4)]:
        gate=R.residual_gate(ee,qq,ell,nn,2,0,**certified)
        check('complete genuine minimal-excess numerical profile '+str((ee,qq,nn)),gate['necessary_arithmetic_passes'] and gate['minimal_excess'] and gate['residual_smooth'])
        equality.append(dict(e=ee,q=qq,n=nn,ell=ell,eta=qq))
    for ee,qq in [(4,1),(5,2),(101,37)]:
        gate=R.residual_gate(ee,qq,(ee+1)*qq,1,2,0,**certified)
        check('arbitrary singular branch minimal excess impossible at e>=4 '+str((ee,qq)),not gate['necessary_arithmetic_passes'] and gate['hypotheses_contradicted'])
    check('hypotheses must be certified before excluding a model',not R.residual_gate(4,1,5,1,2,0)['hypotheses_contradicted'])
    check('wrong Milnor divisibility rejected',not R.residual_gate(5,2,16,2,6,41,**certified)['necessary_arithmetic_passes'])
    check('genuine higher-residue finite A1 family remains arithmetically allowed',R.residual_gate(5,2,16,2,6,40,**certified)['necessary_arithmetic_passes'])
    x,y,u,v,z0=S.symbols('x y u v z')
    # Genuine degree3 cusp, with reduced smooth residual curve.
    pc=x;qc=y**3+x*y;ac=x+3*y*y
    Fc=4*u**3+27*v**2;Hc=4*x+3*y*y
    ident('actual cubic full pullback factorization',Fc.subs({u:pc,v:qc})-ac**2*Hc)
    ident('actual cubic sole local critical branch',S.det(S.Matrix([[S.diff(pc,x),S.diff(pc,y)],[S.diff(qc,x),S.diff(qc,y)]]))-ac)
    check('actual cubic residual is smooth',S.diff(Hc,x).subs({x:0,y:0})!=0)
    # Genuine degree8 smooth and degree4 A1 quotient profiles, SAME upstairs map.
    model=T.rank_zero_model(3);X,Y=model['X'],model['Y'];U,V=model['u'],model['v'];A=model['critical']
    ident('actual smooth e3 q2 full pullback factorization',9*U**2-16*V**3-A**3*(9*X-Y**3)/4)
    ident('actual e3 map descends through minusI in both coordinates',U.subs({X:-X,Y:-Y})-U)
    ident('actual e3 second coordinate descends through minusI',V.subs({X:-X,Y:-Y})-V)
    check('actual smooth e3 q2 and A1 e3 q1 degree data retained',model['local_surface_degree']==8 and model['upstairs_residue_degree']==2 and T.quotient_counterfamily(1)['quotient_degree']==4)
    # Genuine higher-characteristic finite germ. Compute its exact branch and
    # residual equation, rather than prescribing a Milnor number.
    Pg=21*x*x+42*x*y**3-7*y**6+42*x*x*y+28*x*y**4-6*y**7;Qg=x*y;Ag=x-y**3
    F=S.resultant(z0**4-v,u-56*z0**6-64*z0**7,z0)
    ident('higher-characteristic exact parametrization',F.subs({u:56*z0**6+64*z0**7,v:z0**4}))
    ident('higher-characteristic actual unit Jacobian',S.det(S.Matrix([[S.diff(Pg,x),S.diff(Pg,y)],[S.diff(Qg,x),S.diff(Qg,y)]]))-42*(1+y)*Ag**2)
    quotient,remainder=S.div(S.Poly(S.expand(F.subs({u:Pg,v:Qg})),x,y),S.Poly(Ag**3,x,y))
    check('higher-characteristic actual residual exact quotient',remainder.is_zero)
    H=quotient.as_expr()
    image_lengths=[local_gradient_length(F,(u,v),N) for N in (17,18)]
    residual_lengths=[local_gradient_length(H,(x,y),N) for N in (25,26)]
    check('higher-characteristic image exact local conductor16 by Nakayama',all(d['length']==16 for d in image_lengths))
    check('higher-characteristic actual residual exact local Milnor60 by Nakayama',all(d['length']==60 for d in residual_lengths))
    check('higher-characteristic arbitrary-branch identity passes',R.residual_gate(3,1,8,1,16,60,**certified)['necessary_arithmetic_passes'])
    ident('higher-characteristic first conductor intersection',S.Poly(S.expand(H.subs(x,y**3)),y).terms()[-1][0][0]-16)
    check('higher-characteristic branch is genuinely not a torus cusp',S.gcd(4,6)==2 and S.gcd(2,7)==1)
    # Genuine PF107 finite high-index singular model is not excluded locally.
    old=C.explicit_pair(2);X,Y=old['X'],old['Y'];Hp=-(9*X**3+8*Y**3)/64
    ident('genuine PF107 residual triple point',old['P']**3-old['Q']**2-X**9*Hp)
    check('genuine PF107 ordinary triple point exact Milnor4',local_gradient_length(Hp,(X,Y),5)['length']==4)
    check('genuine PF107 finite A1 cover remains allowed',R.residual_gate(9,1,12,2,2,4,**certified)['necessary_arithmetic_passes'])
    for ff,want in [(2,16),(3,23),(5,35),(101,650)]:
        check('PF114 all-degree nonprimitive residue floor '+str(ff),R.nonprimitive_residue_degree_floor(ff)==want)
    check('every omitted residue degree one through genericdegree15',R.nonprimitive_residue_degree_floor(2)==16)
    check('arbitrary nonimmersive image parameter costs keep all preimages',R.nonimmersion_cost(7,3,2,5)==11 and R.nonimmersion_cost(7,3,2)==8)
    text=' '.join((P/'RESIDUAL_MILNOR_IDENTITY_AND_MINIMAL_EXCESS_CLASSIFICATION.md').read_text().lower().split())
    for phrase in ['arbitrary irreducible singular image branch','mu(h)=n*(eta-q)*(c-1)','not the prime-index summand','one local image branch','ordinary preimage curves','not extra omitted primes','p need not be a generic linear source coordinate','negative orbifold euler','closed horizontal torus','original local cyclic-axis pair','no one-puiseux-pair assumption','not a claim of a unique analytic equivalence class','credited to it','not inversion through degree15','nakayama','actual finite high-index cover remains allowed']:
        check('written precise scope '+phrase,phrase in text)
    files=['RESIDUAL_MILNOR_IDENTITY_AND_MINIMAL_EXCESS_CLASSIFICATION.md','residual_milnor.py','verify_residual_milnor.py','sources/residual-milnor/manifest.json']+oldfiles+[x['path'] for x in manifest['sources']]
    result=dict(claim='PF2RULED115',outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={f:sha(P/f) for f in dict.fromkeys(files)},
        theorem=dict(all_degrees_indices_local_residues=True,any_singular_chosen_plane_branch=True,
          smooth_or_cyclic_axis_source_pair_required=True,sole_critical_branch_required=True,
          exact_residual_Milnor_identity=True,normalized_excess_at_least_local_residue=True,
          minimal_excess_complete_numerical_classification=True,
          index_at_least4_forces_excess_at_least_q_plus1=True,
          arbitrary_global_image_nonimmersion_costs=True,
          nonprimitive_omitted_residue_strict_excess=True,
          all_omitted_residues_one_through_degree15=True,
          degree15_corollary_already_follows_from_PF114=True),
        minimal_excess_profiles=equality,
        higher_characteristic_control=dict(image_characteristic=[4,6,7],e=3,q=1,n=1,ell=8,eta=5,
          exact_image_lengths=image_lengths,exact_residual_lengths=residual_lengths,image_conductor=16,residual_Milnor=60,
          global_Keller_map=False,local_finite_germ_excluded=False),
        genuine_finite_high_index_control=dict(degree=12,index=9,source='A1',residual_Milnor=4,excluded_locally=False),
        assurance=dict(characteristic_zero_required=True,polar_Milnor_JSJ_inputs_imported=True,
          written_unbounded_proof=True,one_Puiseux_pair_assumed=False,scalar_Morse_assumed=False,
          source_P_coordinate_assumed_linear=False,residual_curve_assumed_smooth=False,
          global_image_assumed_contractible=False,ordinary_preimages_counted_as_omitted_primes=False,
          PF114_linear_mass_bound_extended_to_arbitrary_branch=False,
          finite_samples_prove_theorem=False,jet_stabilization_assumed_without_Nakayama=False,
          all_degree_residue_one_claimed=False,degree15_inversion_claimed=False,
          unique_analytic_equivalence_class_claimed=False,arithmetic_gates_realize_covers=False,
          actual_finite_algebra_closed=False,formal_verification=False,external_review=False,
          research_goal_completed=False),environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'residual_milnor_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',len(checks),'PF115 controls; arbitrary-branch residual Milnor identity and minimal excess classified; full goal ACTIVE')


if __name__=='__main__':main()
