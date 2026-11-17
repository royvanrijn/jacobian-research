"""Exact controls for PF2RULED58; the unbounded Frobenius argument is written."""
from pathlib import Path
import hashlib,json,os,platform
import sympy as S
HERE=Path(__file__).resolve().parent
checks=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(name,value):
    assert value,name
    checks.append(name)
def quad_matrix(f,variables):
    return S.hessian(f,variables)/2
def hilbert(f,variables,degree):
    g=S.groebner(f,*variables)
    monoms=[]
    def exps(n,total):
        if n==1:
            yield (total,);return
        for j in range(total+1):
            for tail in exps(n-1,total-j):yield (j,)+tail
    for powers in exps(len(variables),degree):
        m=S.prod(v**p for v,p in zip(variables,powers))
        monoms.append(S.Poly(g.reduce(m)[1],*variables))
    support=sorted({m for f in monoms for m,c in f.terms() if c})
    return S.Matrix([[f.coeff_monomial(m) for f in monoms] for m in support]).rank() if support else 0
def main():
    preds=['canonical_formats_receipt.json','unramified_orevkov_receipt.json','trace_dual_pencil_receipt.json','power_difference_receipt.json']
    for n in preds:
        r=json.loads((HERE/n).read_text());check('predecessor PASS '+n,r['outcome']=='PASS')
        for p,h in r['source_sha256'].items():check('predecessor seal '+n+' '+p,sha(HERE/p)==h)
    manifest=HERE/'sources/canonical-quadrics/manifest.json';m=json.loads(manifest.read_text())
    for row in m['primary_sources']:
        p=HERE/row['path_relative_to_affine_complements']
        check('primary source pin '+p.name,p.stat().st_size==row['bytes'] and sha(p)==row['sha256'])
    x,y,z,u,v=S.symbols('X Y Z u v');variables=(x,y,z)
    def split(a,b):return S.matrix_multiply_elementwise(a,b)
    def onefold(a,b):
        return S.Matrix([a[0]*b[0],a[0]*b[1]+a[1]*b[0],a[2]*b[2],a[3]*b[3]])
    def twofold(a,b):
        return S.Matrix([a[0]*b[0],a[0]*b[1]+a[1]*b[0],a[2]*b[2],a[2]*b[3]+a[3]*b[2]])
    models=[
      ('four simple',split,[1,1,1,1],[[1,0,0,-1],[0,1,0,-1],[0,0,1,-1]],[x*y-x*z,x*y-y*z],[(0,1),(1,0),(1,-1)],None),
      ('one fold',onefold,[0,1,1,1],[[1,0,0,0],[0,1,-1,0],[0,1,0,-1]],[x*(y-z),y*z],[(0,1),(1,0)],None),
      ('two folds',twofold,[0,1,0,1],[[1,0,0,0],[0,0,1,0],[0,1,0,-1]],[x*y,z*z],[(0,1),(1,0)],(0,1))]
    model_rows=[]
    for name,mul,weights,columns,quadrics,roots,rank1 in models:
        basis=[S.eye(4)[:,i] for i in range(4)];lam=S.Matrix(1,4,weights)
        H=S.Matrix.hstack(*[S.Matrix(col) for col in columns])
        gram=S.Matrix(4,4,lambda i,j:(lam*mul(basis[i],basis[j]))[0])
        check(name+' Frobenius pairing nondegenerate',gram.det()!=0)
        check(name+' annihilator basis has rank3',H.rank()==3 and lam*H==S.zeros(1,3))
        prods=[mul(H[:,i],H[:,j]) for i in range(3) for j in range(i,3)]
        products=S.Matrix.hstack(*prods)
        check(name+' products span full fibre',products.rank()==4)
        for qi,q in enumerate(quadrics):
            coeff=S.Matrix([S.Poly(q,*variables).coeff_monomial(variables[i]*variables[j]) for i in range(3) for j in range(i,3)])
            check(name+' actual relation '+str(qi),products*coeff==S.zeros(4,1))
        check(name+' two independent quadrics are full relation space',S.Matrix([[S.Poly(q,*variables).coeff_monomial(variables[i]*variables[j]) for i in range(3) for j in range(i,3)] for q in quadrics]).rank()==2 and len(products.nullspace())==2)
        check(name+' quadrics no common component',S.gcd(*quadrics)==1)
        hf=[hilbert(quadrics,variables,k) for k in range(5)]
        check(name+' canonical Hilbert function',hf==[1,3,4,4,4])
        Q=[quad_matrix(q,variables) for q in quadrics];det=S.factor((u*Q[0]+v*Q[1]).det())
        check(name+' determinant is nonzero binary cubic',det!=0 and S.Poly(det,u,v).total_degree()==3)
        ranks=[]
        for r in roots:
            M=r[0]*Q[0]+r[1]*Q[1];ranks.append(M.rank())
            check(name+' singular parameter '+str(r),M.det()==0 and M.rank()==(1 if r==rank1 else 2))
        factor=u*v*(u+v) if name=='four simple' else u*u*v
        check(name+' complete determinant factorization',S.cancel(det/factor).is_number and S.cancel(det/factor)!=0)
        model_rows.append(dict(name=name,hilbert=hf,determinant=str(det),singular_ranks=ranks))
    # A nondegenerate embedding can still be noncanonical and fail arithmetic Gorensteinness.
    points=[(1,0,0),(0,1,0),(1,1,0),(0,0,1)]
    monoms=[x*x,x*y,x*z,y*y,y*z,z*z]
    evals=S.Matrix([[q.subs(dict(zip(variables,p))) for q in monoms] for p in points])
    bad_quads=[S.expand(sum(c*q for c,q in zip(col,monoms))) for col in evals.nullspace()]
    check('collinear negative control still has two quadrics',len(bad_quads)==2)
    check('collinear negative control has common line',S.gcd(*bad_quads)==z)
    check('collinear negative control needs an extra cubic',S.groebner(bad_quads,*variables).reduce(x*y*(x-y))[1]!=0)
    check('canonical reduced frame has no collinear triple',all(S.Matrix([p for j,p in enumerate([(1,0,0),(0,1,0),(0,0,1),(1,1,1)]) if j!=i]).det()!=0 for i in range(4)))
    # Finite controls for the written arbitrary-n diagonal-pairing proof, not a search bound.
    for n in [2,3,5]:
        zz=S.symbols('z0:'+str(n));relations=[zz[i]*zz[j] for i in range(n) for j in range(i+1,n)]+[zz[i]**2-zz[0]**2 for i in range(1,n)]
        g=S.groebner(relations,*zz)
        check('Artinian pairing control '+str(n)+' has h-vector',(hilbert(relations,zz,0),hilbert(relations,zz,1),hilbert(relations,zz,2),hilbert(relations,zz,3))==(1,n,1,0))
        check('Artinian pairing control '+str(n)+' all cubes killed',all(g.reduce(q**3)[1]==0 for q in zz))
        check('Artinian pairing control '+str(n)+' correct quadratic count',len(relations)==n*(n+1)//2-1)
    t=S.symbols('t')
    check('degree3 needs cubic control',hilbert([t**3],(t,),2)==1 and hilbert([t**3],(t,),3)==0)
    degree=S.symbols('degree',integer=True,positive=True)
    check('all-degree kernel rank identity',S.expand((degree-1)*degree/2-degree-degree*(degree-3)/2)==0)
    check('higher-degree canonical quadrics exceed codimension',S.factor(degree*(degree-3)/2-(degree-2))==(degree-4)*(degree-1)/2)
    a,b,c,d,e,f=S.symbols('a b c d e f')
    q1=x*y+a*y*y+b*y*z+c*z*z;q2=x*z+d*y*y+e*y*z+f*z*z
    residual=d*y**3+(e-a)*y*y*z+(f-b)*y*z*z-c*z**3
    check('universal residual elimination identity',S.expand(y*q2-z*q1-residual)==0)
    Q1=quad_matrix(q1,variables);Q2=quad_matrix(q2,variables)
    check('universal determinant equals full residual cubic',S.expand(-4*(u*Q1+v*Q2).det()-residual.subs({y:v,z:-u},simultaneous=True))==0)
    check('source diagonal tangent matrix exactly unit',S.Matrix([[S.diff(q,y).subs({x:1,y:0,z:0}),S.diff(q,z).subs({x:1,y:0,z:0})] for q in [q1,q2]])==S.eye(2))
    inputs=[Path(__file__),HERE/'CANONICAL_QUADRICS_AND_SOURCE_RESIDUAL.md',manifest]+[HERE/n for n in preds]+[HERE/r['path_relative_to_affine_complements'] for r in m['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      general=dict(all_degrees_at_least_three=True,Gorenstein_hypothesis_explicit=True,all_normalizations_assumed_Gorenstein=False,canonical_uses_dual_annihilator_of_one=True,canonical_embedding_proved=True,arithmetic_Gorensteinness_proved=True,hilbert_function='1,d-1,d,d,...',Artinian_h_vector='1,d-2,1',quadratic_generation_for_all_d_at_least_four=True,quadratic_rank='d(d-3)/2',modules_free_by_Quillen_Suslin=True,non_Gorenstein_cases_excluded=False),
      quartic=dict(all_surviving_rulings=True,all_source_pole_orders=True,whole_normalization_is_global_two_quadric_intersection=True,full_fibre_types=['1+1+1+1','2+1+1','2+2'],rank_one_pencil_members_exactly_double_branch_values=True,determinant_cubic_is_actual_normal_C3=True,conditional_PF57_pencil_gate_now_unconditional=True,every_nonzero_polynomial_member_determinant_nonzero_nonconstant=True,finite_Keller_pair_constructed=False),
      source=dict(full_normalization_fibre_product_retained=True,diagonal_component_split=True,residual_smooth_finite_flat_degree=3,residual_is_actual_normal_cubic_pullback=True,residual_to_N_etale_quasifinite=True,ordinary_source_offdiagonal_automatically_finite=False,source_side_projective_and_pencil_changes=True,changes_asserted_to_descend_to_target=False,universal_residual_identity=True,pullback_branch_smooth=True,pullback_branch_normalizations_assumed_A1=False,residual_assumed_rational_or_plane_or_monogenic=False),
      primary_correction=dict(year=2026,arbitrary_nondegenerate_embedding_not_sufficient=True,canonical_fibre_property_checked_directly=True,original_unrestricted_uniqueness_imported=False),
      finite_model_controls=model_rows,
      assurance=dict(written_unbounded_proofs=True,finite_controls_not_universal_proof=True,primary_imports_independently_reproved=False,formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,full_quartic_exclusion=False,higher_degree_exclusion=False,research_goal_completed=False),
      source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'canonical_quadrics_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; canonical quadrics and full source residual')
if __name__=='__main__':main()
