"""Selected exact filtration/extension controls; no universal MHS verifier."""
from pathlib import Path
import hashlib,json
import sympy as s

def relative_coordinate_filtration(N,W,M):
    n=N.rows
    if any(N[i,j]!=0 and M[i]>M[j]-2 for i in range(n) for j in range(n)):
        return False
    if any(N[i,j]!=0 and W[i]>W[j] for i in range(n) for j in range(n)):
        return False
    for k in sorted(set(W)):
        ids=[i for i in range(n) if W[i]==k]
        B=N.extract(ids,ids);weights=[M[i] for i in ids]
        for ell in range(n+1):
            src=[i for i,w in enumerate(weights) if w==k+ell]
            dst=[i for i,w in enumerate(weights) if w==k-ell]
            if len(src)!=len(dst) or (B**ell).extract(dst,src).rank()!=len(src):
                return False
    return True

# Nonsemisimple fibre plus a genuine two-weight-drop Tate extension model.
N=s.zeros(4);N[1,2]=1;N[0,3]=1
W=[4,5,5,6];M=[4,4,6,6]
assert relative_coordinate_filtration(N,W,M)
NK=N[1:,1:]
invariants=NK.nullspace()
assert invariants==[s.Matrix([1,0,0]),s.Matrix([0,0,1])]
connecting=[]
for v in invariants:
    value=N*s.Matrix([0]+list(v))
    assert value[1:]==[0,0,0]
    connecting.append(value[0])
assert connecting==[0,1] # Weight5 invariant killed; top Tate invariant detected.
# A lower-weight connecting value would violate the relative drop here.
invalid=N.copy();invalid[0,1]=1
assert not relative_coordinate_filtration(invalid,W,M)

# The image J of invariants can be smaller than the entire top Tate grade.
K=s.Matrix([[0,1],[0,0]])
assert relative_coordinate_filtration(K,[4,6],[4,6])
assert K.nullspace()==[s.Matrix([1,0])]
assert all(v[1]==0 for v in K.nullspace())

# Independent-block Jordan controls for the invariant-vector bound.
kernel_bound=[]
for centre in (4,5,6):
    for size in (1,2,3,4):
        B=s.zeros(size)
        for i in range(size-1):B[i,i+1]=1
        weights=[centre-(size-1)+2*i for i in range(size)]
        assert relative_coordinate_filtration(B,[centre]*size,weights)
        assert weights[0]<=centre
        kernel_bound.append({'centre':centre,'Jordan_size':size,
                             'kernel_M_weight':weights[0]})

# A split intrinsic extension can leave nonzero actual fibre motion.
A=s.eye(3);A[0,1]=1;A[1,2]=1
T=s.diag(A,s.ones(1));fixed=s.Matrix([0,0,0,1])
actual=s.Matrix([0,0,1,1]);motion=(T-s.eye(4))*actual
assert T*fixed==fixed and motion==s.Matrix([0,1,0,0])
assert (A-s.eye(3)).row_join(motion[:3,:]).rank()==(A-s.eye(3)).rank()
assert s.Matrix([0,0,1]) not in (A-s.eye(3)).nullspace()
t=s.symbols('t');changed_fixed=fixed+s.Matrix([t,0,0,0])
assert (T-s.eye(4))*(actual-changed_fixed)==motion

def cyclic(B,e):
    n=B.rows;A=s.zeros(n*e)
    for j in range(e-1):A[(j+1)*n:(j+2)*n,j*n:(j+1)*n]=s.eye(n)
    A[:n,(e-1)*n:e*n]=B
    return A

def powers(A,k):return sum((A**j for j in range(k)),s.zeros(A.rows))

# Nonsemisimple original eigenvalue1 AND another eigenvalue are retained.
B=s.Matrix([[1,1,0],[0,1,0],[0,0,-1]])
qB=s.Matrix([[0,1,0]])
transfer=[]
for e in (2,3):
    C=cyclic(B,e);sum_map=s.Matrix.hstack(*([s.eye(3)]*e))
    expected=s.zeros(3,3*e);expected[:,3*(e-1):]=B-s.eye(3)
    assert sum_map*(C-s.eye(3*e))==expected
    qC=qB*sum_map
    assert qC*(C-s.eye(3*e))==s.zeros(1,3*e)
    assert 3*e-(C-s.eye(3*e)).rank()==3-(B-s.eye(3)).rank()==1
    d=s.Matrix(range(1,3*e+1));returns=powers(C,e)*d
    for j in range(e):assert qB*returns[3*j:3*j+3,:]==qC*d
    transfer.append({'multiplicity':e,'fibre1_nonsemisimple':True,
                     'coinvariant_factor':'1','quotient_isomorphism_exact':True})

# Distinct roots remain distinct quotient summands even if a scalar sum cancels.
C=s.diag(cyclic(s.ones(1),2),cyclic(s.ones(1),3))
q=s.Matrix([[1,1,0,0,0],[0,0,1,1,1]])
d=s.Matrix([1,0,-1,0,0])
assert q*(C-s.eye(5))==s.zeros(2,5)
assert q*d==s.Matrix([1,-1]) and sum(q*d)==0 and q*d!=s.zeros(2,1)

# Fixed representative factor1/e and quotient factor1 agree in semisimple case.
for e in (2,3):
    fixed_class=s.ones(e,1)/e
    assert s.ones(1,e)*fixed_class==s.ones(1)

out={
 'completed':True,'sympy_version':s.__version__,
 'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'relative_filtration_nonsemisimple_Tate_connecting_control_exact':True,
 'lower_weight_connecting_candidate_rejected_by_relative_drop':True,
 'invariant_image_J_can_be_proper_in_top_Tate_grade':True,
 'invariant_vector_kernel_weight_controls':kernel_bound,
 'zero_intrinsic_class_nonzero_actual_motion_control_exact':True,
 'general_coinvariant_transfer_controls':transfer,
 'different_roots_do_not_cancel_in_coinvariants':True,
 'fixed_factor_one_over_e_sums_to_quotient_factor_one':True,
 'universal_limit_exactness_Tate_factor_and_transfer_proofs_machine_verified':False,
 'Tate_factor_identified_with_entire_top_Tate_space':False,
 'general_actual_Haar_motion_identified_with_Tate_subspace':False,
 'zero_intrinsic_class_claimed_general_positive_MZ_theorem':False,
 'all_exceptional_intrinsic_or_residual_obstructions_proved_zero':False,
 'general_MHS_canonical_lattice_or_cycle_producer_implemented':False,
 'abstract_controls_claimed_as_exceptional_phase_realizations':False,
 'all_coefficient_ring_topological_promotion_claimed':False,
 'full_exceptional_MZ_cases':'OPEN on BOTH quotients',
 'independent_replay':False,'formal_verification':False,'external_review':False,
 'new_degree_support_or_prefix_campaign':False}
Path(__file__).with_name('haar_unipotent_tate_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
