"""Exact transfer and example controls for the companion written proof."""
from pathlib import Path
import hashlib,json
import sympy as s

z,u,alpha,beta=s.symbols('z u alpha beta',nonzero=True)
e=s.symbols('e',integer=True,positive=True)
j=s.symbols('j',integer=True,nonnegative=True)
m=z**j
A=alpha+beta*z
delta=lambda f:s.diff(f,z)+A*f/z-e*z**(e-1)*u*f
R=lambda f:z**e*f-s.diff(f,u)
Es=lambda f:z*s.diff(f,z)+A*f
assert s.simplify(e*R(u*m)-Es(m)-(1-e)*m+delta(z*m))==0
rows=[]
for ee in [2,3]:
    A0=s.Matrix([[s.Rational(-1,2),1],[0,s.Rational(-1,2)]])
    A1=s.Matrix([[0,1],[1,0]])
    T=s.zeros(2*ee)
    for jj in range(ee):
        T[2*jj:2*jj+2,2*jj:2*jj+2]=(A0+(jj+1-ee)*s.eye(2))/ee
        if jj+1<ee:T[2*(jj+1):2*(jj+1)+2,2*jj:2*jj+2]=A1/ee
    q=s.symbols('q')
    expected=s.prod((q-(s.Rational(-1,2)+jj+1)/ee+1)**2 for jj in range(ee))
    assert s.factor(T.charpoly(q).as_expr()-expected)==0
    vals=[(s.Rational(-1,2)+jj+1)/ee-1 for jj in range(ee)]
    assert all(-1<v<=0 for v in vals)
    rows.append({'ramification':ee,'nonsemisimple_source_residues':['-1/2','-1/2'],'target_residues':[str(v) for v in vals],'all_strictly_above_minus_one':True})
# The genuine trace-square moments follow from the original Haar central term.
rr=s.symbols('r',positive=True)
C=(1-s.sqrt(1-4/rr))/2
ser=s.series(C,rr,s.oo,8).removeO().expand()
mom=[]
for n in range(7):
    cat=s.binomial(2*n,n)/(n+1)
    assert ser.coeff(rr,-n-1)==cat
    mom.append(str(cat))
# Raw quotient versus original-form image: the displayed top primitive is exact.
x,y,v=s.symbols('x y v')
for Q in [s.Integer(1),y,x*y**2]:
    primitive=s.integrate(Q,y)
    phase=x**2
    bx,by=-primitive,s.Integer(0)
    d_top=s.diff(by,x)-s.diff(bx,y)
    wedge_top=s.diff(phase,x)*by-s.diff(phase,y)*bx
    assert s.expand(d_top-v*wedge_top-Q)==0
out={
 'completed':True,'sympy_version':s.__version__,
 'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'symbolic_finite_transfer_Euler_identity_exact':True,
 'nonsemisimple_residue_and_higher_Taylor_block_controls':rows,
 'trace_square_Catalan_moments':mom,
 'raw_Brieskorn_quotient_distinguished_from_original_image':True,
 'universal_direct_image_and_V_lattice_proofs_machine_verified':False,
 'finite_original_image_identified_with_raw_Brieskorn_module':False,
 'primitive_purity_reflection_claimed':False,
 'all_coefficient_ring_topological_or_lattice_promotion_claimed':False,
 'full_exceptional_MZ_cases':'OPEN on BOTH quotients',
 'independent_replay':False,'formal_verification':False,'external_review':False,
 'new_degree_support_or_prefix_campaign':False}
Path(__file__).with_name('haar_finite_recomposition_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
