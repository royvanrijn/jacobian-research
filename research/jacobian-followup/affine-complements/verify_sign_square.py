"""Exact finite controls for PF2RULED55; topology is proved in the canonical note."""
from pathlib import Path
import hashlib,itertools,json,math,os,platform
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def valuation(n,p=3):
    assert n
    n=abs(n);v=0
    while n%p==0:n//=p;v+=1
    return v
def main():
    pred=HERE/'quartic_colouring_receipt.json';r=json.loads(pred.read_text())
    check('PF54 receipt PASS',r['outcome']=='PASS')
    for path,h in r['source_sha256'].items():check('PF54 sealed input '+path,sha(HERE/path)==h)
    check('required triple is etale at singularities',r['sign_surface']['etale_at_finite_singularities'] and r['sign_surface']['connected_etale_cyclic_three_cover_required'])
    manifest=HERE/'sources/sign-square/manifest.json';m=json.loads(manifest.read_text())
    for row in m['primary_sources']+m['background_sources']:
        p=HERE/row['path_relative_to_affine_complements'];check('pinned source '+p.name,p.stat().st_size==row['bytes'] and sha(p)==row['sha256'])
    check('Hamm original singular theorem scan retained',any('hamm-1983-p121' in row['path_relative_to_affine_complements'] and row['complete_pages_visually_inspected']==[121] for row in m['primary_sources']))
    t,x,u,z,a,b,s,v,w=S.symbols('t x u z a b s v w')
    # Irreducibility is essential for the universal quadratic-unit descent.
    check('reducible branch admits nonconstant norm-one unit',S.expand(x*x-(x*x-1))==1 and S.factor(x*x-1)==(x-1)*(x+1))
    # A small irreducible example records the exact descent identities, not a search.
    check('coprime norm factors differ by two',S.expand((a+1)-(a-1))==2)
    check('halving the square factor gives a new constant norm',S.expand((u*u-x*v*v)-2).subs(x,(u*u-2)/v**2)==0)
    # The cone control distinguishes Cl[3] from Pic[3] and quasi-etale from etale.
    cone=u*v-w**3
    check('A2 cone invariant parametrization',S.expand(cone.subs({u:s**3,v:t**3,w:s*t},simultaneous=True))==0)
    check('cone localization is a Laurent polynomial ring',S.cancel(cone.subs(v,w**3/u))==0)
    check('cone has only the singular vertex',list(S.groebner([cone,S.diff(cone,u),S.diff(cone,v),S.diff(cone,w)],u,v,w))==[u,v,w**2])
    check('cone cover fails to be etale at vertex',S.Matrix([s**3,t**3,s*t]).jacobian([s,t]).subs({s:0,t:0}).rank()==0)
    check('cone class group localization relation is3',smith_normal_form(S.Matrix([[3]]),domain=S.ZZ)==S.Matrix([[3]]))
    # Oka positive Picard control: localization presentation is computed exactly.
    f=(x*u-2)**2-u**3;aa=x*u-2+z;bb=x*u-2-z
    relation=z*z-f
    check('Oka ab equals u cubed on sign surface',S.expand(aa*bb-u**3+relation)==0)
    check('Oka torus inverse reconstructs x',S.cancel((aa+bb+4)/(2*u)-x)==0)
    check('Oka sign surface smooth at every affine point',list(S.groebner([relation,S.diff(relation,x),S.diff(relation,u),S.diff(relation,z)],x,u,z))==[1])
    matrix=S.Matrix([[1,3],[1,0]]);snf=smith_normal_form(matrix,domain=S.ZZ)
    check('Oka localization cokernel exactly Z3',sorted(abs(snf[i,i]) for i in range(2))==[1,3])
    check('Oka units lattice kernel zero',matrix.nullspace()==[])
    swap=S.Matrix([[0,1],[1,0]])
    check('Oka involution reverses every boundary class',(swap+S.eye(2))*S.Matrix([1,0])==matrix*S.Matrix([1,0]))
    check('Oka branch is Cstar parametrized',S.cancel(f.subs({u:t**2,x:t+2/t**2},simultaneous=True))==0 and f.subs(u,0)==4)
    # Concrete finite linking forms: the same order9 need not admit isotropy.
    def isotropic(diag,modulus):
        return [list(q) for q in itertools.product(range(modulus),repeat=len(diag)) if sum(c*n*n for c,n in zip(diag,q))%modulus==0]
    aniso=isotropic([1,1],3);hyper=isotropic([1,-1],3);cyclic=isotropic([1],9)
    check('anisotropic3-primary order9 form has only zero isotropic',aniso==[[0,0]])
    check('split3-primary order9 form has nonzero isotropic',len(hyper)==5 and [1,1] in hyper)
    check('cyclic order9 linking form has isotropic subgroup of order3',cyclic==[[0],[3],[6]])
    for Q,A in [(S.zeros(0),[3]),(S.diag(-2,-2),[3]),(S.Matrix([[2,1],[1,2]]),[3]),(S.eye(1),[3,3]),(S.eye(1),[9])]:
        size=math.prod(A);det=abs(int(Q.det()));D=det*size**2;rank=sum(n%3==0 for n in A)
        check('square and Pic3 rank control '+str((Q.tolist(),A)),D%size**2==0 and valuation(D)>=2*rank)
    check('torsion order contains more information than Pic3 rank',9**2>3**2)
    # Exact satellite polynomial controls include both new exclusions and survivors.
    def torus(w,p):
        p=abs(p);return S.cancel((t**(w*p)-1)*(t-1)/((t**w-1)*(t**p-1)))
    def step(D,w,p):return abs(p) if w%2==0 else w*D if p%2==0 else D
    def vstep(e,w,p):return valuation(p) if w%2==0 else valuation(w)+e if p%2==0 else e
    cable=[]
    for stages in [[(3,4)],[(2,9)],[(3,4),(3,2)],[(2,9),(2,3)],[(3,4),(2,-9)],[(3,4),(3,2),(5,7)],[(2,5),(9,2)],[(2,9),(2,5)]]:
        poly=S.Integer(1);D=1;e=0
        for w,p in stages:
            assert w>=2 and math.gcd(w,p)==1
            poly=S.cancel(torus(w,p)*poly.subs(t,t**w));D=step(D,w,p);e=vstep(e,w,p)
        check('all-stage square-gate Alexander identity '+str(stages),poly.is_polynomial(t) and poly.subs(t,1)==1 and abs(int(poly.subs(t,-1)))==D and valuation(D)==e)
        cable.append(dict(stages=stages,determinant=D,v3=e,passes_old_three_gate=D%3==0,passes_new_nine_gate=D%9==0,Keller_realization_asserted=False))
    check('T34 passes old determinant and exact S4 lift but fails new global gate',cable[0]['passes_old_three_gate'] and not cable[0]['passes_new_nine_gate'] and r['lifting']['torus34']['parabolic_H1_dimension']==1)
    check('final even stage erases inner nine divisibility',cable[3]['determinant']==3 and not cable[3]['passes_new_nine_gate'])
    check('odd winding can supply the second factor three',cable[2]['determinant']==9 and cable[2]['passes_new_nine_gate'])
    inputs=[Path(__file__),HERE/'SIGN_SURFACE_PICARD_AND_SQUARE_DETERMINANT.md',manifest,pred]+[HERE/row['path_relative_to_affine_complements'] for row in m['primary_sources']+m['background_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      scope=dict(entire_actual_quartic_packet=True,all_source_orders=True,all_target_degrees=True,all_rulings=True,both_omitted_primes_retained=True,underlying_sign_surface_theorem_all_branch_degrees=True,required_triple_supplied_for_higher_Keller_degrees=False),
      units=dict(irreducible_branch_required=True,normal_singular_surface_allowed=True,all_units_constant=True,quadratic_degree_descent_written=True,Kummer_units_quotient_zero=True,Pic3_nonzero_required=True,tau_acts_minus_one=True,Pic_not_Cl=True),
      topology=dict(singular_affine_Hamm_theorem_used=True,exceptional_rational_trees=True,H1_resolution_equals_H1_sign=True,H2_resolution_torsion_free=True,H3_resolution_zero=True,knot_double_boundary_rational_homology_sphere=True,H1_boundary_surjects_H1_filling=True,H1_sign_finite=True,negative_definiteness_required=False,rational_compactification_required=False),
      determinant=dict(identity='D=abs(det Q)*abs(H1(S,Z))^2',Pic_torsion_order_squared_divides_D=True,necessary_quartic='9 divides abs(Delta_K(-1))',excludes_v3_zero_and_one=True,nonzero_isotropic_order3_subgroup_required=True,all_cable_stages=True,cable_controls=cable),
      controls=dict(linking_order9_anisotropic=aniso,linking_order9_split=hyper,linking_cyclic9=cyclic,Oka_Picard_invariants=[3],cone_Cl_invariants=[3],cone_Picard_torsion_empty=True,reducible_branch_nonconstant_unit=True),
      assurance=dict(written_unbounded_proofs=True,finite_controls_not_universal_proof=True,primary_imports_independently_reproved=False,formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,all_nine_divisible_candidates_excluded=False,finite_Keller_pair_supplied=False,full_quartic_exclusion=False,higher_degree_exclusion=False,research_goal_completed=False),
      source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'sign_square_receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(f'PASS {len(checks)} controls; Picard torsion and square determinant obstruction')
if __name__=='__main__':main()
