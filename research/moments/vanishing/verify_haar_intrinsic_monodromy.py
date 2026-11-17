"""Exact extension, canonical-residue and ramification controls."""
from pathlib import Path
import hashlib,json
import sympy as s

a,b,x,y,z=s.symbols('a b x y z')
T=s.eye(3);T[0,2]=a;T[1,2]=b
lift=s.Matrix([x,y,1])
assert (T-s.eye(3))*lift==s.Matrix([a,b,0])
assert (T-s.eye(3))**2==s.zeros(3)

# Normalize a resonant logarithmic system before using its residue.
Araw=s.Matrix([[0,0],[-1,-1/z]])
P=s.diag(1,z)
Acan=s.simplify(P.inv()*Araw*P+P.inv()*P.diff(z))
Rraw=s.simplify(z*Araw).subs(z,0)
Rcan=s.simplify(z*Acan)
assert Rraw==s.diag(0,-1)
assert Rcan==s.Matrix([[0,0],[-1,0]]) and Rcan**2==s.zeros(2)
monodromy=s.eye(2)-2*s.pi*s.I*Rcan
assert monodromy!=s.eye(2)

# Arbitrary fibre nilpotent: the connecting obstruction is ell on ker R_F.
RF=s.Matrix([[0,1],[0,0]])
RU=s.Matrix([[0,a,b],[0,0,1],[0,0,0]])
kernel_lift=s.Matrix([0,1,0])
assert RU*kernel_lift==s.Matrix([a,0,0])
actual=s.Matrix([1,x,y])
assert RU.T*actual==s.Matrix([0,a,b+x])
assert RU.subs(a,0).T*s.Matrix([1,-b,0])==s.zeros(3,1)
assert RU.subs(a,0).T*s.Matrix([1,0,0])==s.Matrix([0,0,b])

def cyclic(M,e):
    n=M.rows;A=s.zeros(n*e)
    for j in range(e-1):A[(j+1)*n:(j+2)*n,j*n:(j+1)*n]=s.eye(n)
    A[:n,(e-1)*n:e*n]=M
    assert A**e==s.diag(*([M]*e))
    return A

def eig1(A):
    q=s.symbols('q');cp=A.charpoly(q).as_expr();m=0
    while s.rem(cp,q-1,q)==0:
        cp=s.div(cp,q-1,q)[0];m+=1
    if not m:return s.zeros(A.rows)
    pol=s.rem(cp*s.invert(cp,(q-1)**m,q),cp*(q-1)**m,q)
    E=s.zeros(A.rows)
    for (j,),c in s.Poly(pol,q).terms():E+=c*A**j
    assert E**2==E and E*A==A*E
    return E

def powers(A,k):return sum((A**v for v in range(k)),s.zeros(A.rows))

controls=[]
J=s.eye(3);J[0,1]=1;J[1,2]=1;N=J-s.eye(3)
for e in (2,3):
    B=s.eye(3)+N/e+s.binomial(s.Rational(1,e),2)*N**2
    assert B**e==J
    assert [((B-s.eye(3))**k).rank() for k in (1,2,3)]==[2,1,0]
    # Every Jordan size survives the unipotent e-th root.
    M=s.diag(1,-1);A=cyclic(M,e);EP=eig1(A);EF=s.diag(1,0)
    d=s.Matrix(range(1,2*e+1));c=EP*d;returns=powers(A,e)*d
    for j in range(e):assert c[2*j:2*j+2,:]==EF*returns[2*j:2*j+2,:]/e
    assert (A-s.eye(2*e))*EP==s.zeros(2*e)
    controls.append({'root_multiplicity':e,'ramified_class_factor':'1/'+str(e),
                     'unipotent_root_Jordan_ranks':[2,1,0]})

# The factor rule is not promoted beyond semisimple fibre1.
J2=s.Matrix([[1,1],[0,1]]);A=cyclic(J2,2)
d=s.Matrix([0,1,0,0]);c=eig1(A)*d;returns=(s.eye(4)+A)*d
assert any(c[2*j:2*j+2,:]!=returns[2*j:2*j+2,:]/2 for j in range(2))

# Genuine SL2 square control: primitive -torus class is divided by2.
A=cyclic(s.ones(1),2);d=s.Matrix([-1,0]);c=eig1(A)*d
assert c==s.Matrix([-s.Rational(1,2),-s.Rational(1,2)])
assert (s.eye(2)+A)*d==s.Matrix([-1,-1])
v=s.symbols('v');marker=s.Matrix([[-s.Rational(1,2),-s.Rational(1,2)]])
assert marker*c==s.Matrix([[s.Rational(1,2)]])

# Distinct roots cannot cancel as homology summands even if scalar sums do.
A=s.diag(cyclic(s.ones(1),2),cyclic(s.ones(1),3))
d=s.Matrix([1,0,-1,0,0]);c=eig1(A)*d
assert list(c)==[s.Rational(1,2),s.Rational(1,2),-s.Rational(1,3),-s.Rational(1,3),-s.Rational(1,3)]
assert sum(c)==0 and c!=s.zeros(5,1)

out={
 'completed':True,'sympy_version':s.__version__,
 'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'semisimple_fibre_lift_independence_exact':True,
 'canonical_resonance_countercontrol_exact':True,
 'raw_residue_eigenvalues':['0','-1'],'canonical_residue':[[str(v) for v in row] for row in Rcan.tolist()],
 'nonsemisimple_connecting_obstruction_and_split_lift_countercontrol_exact':True,
 'ramification_and_Jordan_controls':controls,
 'ramified_factor_requires_semisimple_fibre1':True,
 'actual_SL2_square_projected_class':['-1/2','-1/2'],
 'separate_root_scalar_sum_zero_but_class_nonzero':True,
 'universal_extension_ramification_and_canonical_lattice_proofs_machine_verified':False,
 'general_full_connection_or_canonical_lattice_producer_implemented':False,
 'raw_logarithmic_residue_test_claimed_sufficient':False,
 'zero_connecting_map_claimed_general_positive_theorem':False,
 'all_exceptional_fibres_semisimple_at_eigenvalue1':False,
 'all_exceptional_intrinsic_extensions_proved_zero':False,
 'primitive_purity_reflection_claimed':False,
 'all_coefficient_ring_topological_promotion_claimed':False,
 'abstract_matrix_controls_claimed_as_exceptional_phases':False,
 'full_exceptional_MZ_cases':'OPEN on BOTH quotients',
 'independent_replay':False,'formal_verification':False,'external_review':False,
 'new_degree_support_or_prefix_campaign':False}
Path(__file__).with_name('haar_intrinsic_monodromy_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
