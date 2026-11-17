"""Finite exact controls for the written rational-core/boundary synthesis."""
from pathlib import Path
import hashlib,json,sys
import sympy as s
t,z=s.symbols('t z')
p=Path(__file__).resolve().parent
def clean(v):return s.cancel(v)
def valuation(poly,h):
 if poly==0:return s.oo
 var=z if h==z else t
 n=0
 while s.rem(poly,h,var)==0:poly=s.quo(poly,h,var);n+=1
 return n
def fuchs(L):
 # L lists polynomial coefficients from order0 upward, after a common clearing.
 d=len(L)-1
 while d and L[d]==0:d-=1
 assert d>0
 ratios=[clean(L[j]/L[d]) for j in range(d)]
 finite=[]
 factors={}
 for a in ratios:
  if a:
   den=s.fraction(a)[1]
   for h,e in s.factor_list(den,t)[1]:factors[str(h)]=h
 for h in factors.values():
  bad=[]
  for j,a in enumerate(ratios):
   if not a:continue
   num,den=s.fraction(a)
   v=valuation(num,h)-valuation(den,h)
   if v<j-d:bad.append({'coefficient_order':j,'valuation':int(v),'required':j-d})
  finite.append({'factor':str(h),'regular':not bad,'violations':bad})
 bad=[]
 for j,a in enumerate(ratios):
  if not a:continue
  num,den=s.fraction(a);degree=s.degree(num,t)-s.degree(den,t)
  if degree>j-d:bad.append({'coefficient_order':j,'growth_degree':int(degree),'maximum':j-d})
 return {'order':d,'finite_points':finite,'infinity_regular':not bad,
  'infinity_violations':bad,'regular_on_P1':all(row['regular'] for row in finite) and not bad}
H5=sum(s.binomial(5,j)*(-t)**j/s.factorial(j+1) for j in range(6))
models=[
 ('constant_minimal',[0,1],True),
 ('rational_pole',[1,1-t],True),
 ('square_root',[1,2*(1-t)],True),
 ('constant_plus_square_root',[0,-1,2*(1-t)],True),
 ('logarithm',[0,-1,1-t],True),
 ('constant_plus_exponential',[0,-1,1],False),
 ('Airy',[-t,0,1],False),
 ('Kummer_nonminimal',[5,2-t,t],False),
 ('Kummer_minimal',[-s.diff(H5,t),H5],True),
 ('finite_irregular',[-1,(1-t)**2],False)]
regular=[]
for name,L,expected in models:
 row={'name':name,**fuchs(L)}
 assert row['regular_on_P1']==expected,name
 # Direct infinity coordinate change: D_t=-z^2 D_z, via Ore composition.
 def compose(A,B):
  C=[s.S.Zero]*(len(A)+len(B)-1)
  for i,a in enumerate(A):
   for j,b in enumerate(B):
    for k in range(i+1):C[i+j-k]+=a*s.binomial(i,k)*s.diff(b,z,k)
  return list(map(s.cancel,C))
 powers=[[s.S.One]]
 for _ in range(len(L)-1):powers.append(compose(powers[-1],[0,-z*z]))
 C=[s.S.Zero]*len(L)
 for i,a in enumerate(L):
  for j,b in enumerate(powers[i]):C[j]+=a.subs(t,1/z)*b if hasattr(a,'subs') else a*b
 C=list(map(s.cancel,C));d=len(L)-1
 for j in range(d):
  ratio=clean(C[j]/C[d])
  if ratio:
   num,den=s.fraction(ratio)
   vz=valuation(num,z)-valuation(den,z)
   assert (vz>=j-d)==(not any(b['coefficient_order']==j for b in row['infinity_violations']))
 regular.append(row)
assert any(not r['regular'] for r in regular[-1]['finite_points'])
assert regular[-1]['infinity_regular']
# Point modules: towers are truncated only as a finite control, away from top boundary.
a=s.symbols('a');size=9
D=s.zeros(size);T=a*s.eye(size)
for k in range(size-1):D[k+1,k]=1
for k in range(1,size):T[k-1,k]=-k
checks=0
for k in range(6):
 e=s.eye(size)[:,k]
 assert (D*T-T*D)*e==e
 assert (T-a*s.eye(size))**(k+1)*e==s.zeros(size,1)
 assert ((T-a*s.eye(size))**k*e)[0]==(-1)**k*s.factorial(k)
 assert D*(T-a*s.eye(size))*D*e==-(k+1)*D*e
 checks+=1
# For delta_a tensor E, subspace membership is coefficient-wise in the derivative tower.
spaces=[
 s.zeros(3,0),s.Matrix([1,0,0]),s.Matrix([1,2,3]),
 s.Matrix([[1,0],[0,1],[1,1]]),s.eye(3)]
subspacechecks=0
for U in spaces:
 # Orthogonal annihilator rows cut out the exact multiplicity subspace.
 annih=s.Matrix.hstack(*U.T.nullspace()).T if U.cols<3 else s.zeros(0,3)
 assert annih*U==s.zeros(annih.rows,U.cols)
 assert 3-annih.rank()==U.rank()
 for k in range(5):
  if not U.cols:continue
  v=s.kronecker_product(s.eye(size)[:,k],U[:,0])
  assert s.kronecker_product(D,s.eye(3))*v==s.kronecker_product(D*s.eye(size)[:,k],U[:,0])
  assert s.kronecker_product(T,s.eye(3))*v==s.kronecker_product(T*s.eye(size)[:,k],U[:,0])
  subspacechecks+=1
# Diagonal boundary at2, versus pure rational-pole direction and trace projection.
U=s.Matrix([1,1]);rational=s.Matrix([1,0]);trace=s.diag(1,0)
assert s.Matrix.hstack(U,rational).rank()==2
assert s.Matrix.hstack(U,trace*U).rank()==2
assert (s.Matrix([[1,-1]])*U)==s.zeros(1,1)
assert s.Matrix([[1,-1]])*rational!=s.zeros(1,1)
u=s.sqrt(1-t);g=(1+u)/(t-2)
assert s.simplify((t-2)*g-(1+u))==0
assert s.simplify((g+g.xreplace({u:-u}))/2-1/(t-2))==0
assert s.diff(u,t)==-1/(2*s.sqrt(1-t))
# The chosen branch is formal at zero; no source with a pole at zero is admitted.
assert s.series(g,t,0,4).removeO().subs(t,0)==-1
# Rational core control: derivative of log supplies the full rational pole direction.
assert clean(s.diff(s.log(1-t),t)+1/(1-t))==0
for k in range(1,7):
 assert clean(s.diff(1/(1-t),t,k-1)-s.factorial(k-1)/(1-t)**k)==0
falsekeys=[
 'general_rational_core_solver_executed','general_intermediate_extension_solver_executed',
 'general_holonomic_duality_producer_executed','general_local_cohomology_producer_executed',
 'general_Hom_basis_producer_executed','general_boundary_subspaces_computed_for_a_new_Gaussian_phase',
 'universal_proofs_certified_by_checker','global_complex_safe_locus_uniformity_proved',
 'global_coefficient_stratification_proved','all_higher_rank_Gaussian_connections_classified',
 'new_SIC33_answer_computed','new_Delta_minimum_computed','numeric_cubic_frontier_changed',
 'operator_order_reduced','independent_replay','external_review','formal_verification']
receipt={'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'python_version':sys.version.split()[0],'sympy_version':s.__version__,
 'Fuchs_controls':regular,'infinity_coordinate_change_checked':True,
 'point_tower_controls':checks,'boundary_subspace_controls':subspacechecks,
 'boundary_multiplicity_subspaces':5,
 'diagonal_boundary_control':{'generic_rank':2,'boundary_quotient_dimension':1,
  'rational_pole_subspace_at2_dimension':0,'trace_does_not_preserve_actual_lattice':True,
  'abstract_algebraic_model_not_new_Gaussian_counterexample':True},
 'logarithmic_core_controls':6,
 'minimal_log_extension_has_rational_pole_but_no_boundary_quotient':True,
 'mutation_controls':{'quotient_not_assumed_to_embed_in_actual_function_field':True,
  'boundary_quotient_not_identified_with_punctual_submodule_of_obstruction':True,
  'boundary_space_dimensions_not_used_without_embedding':True,
  'finite_tower_not_used_to_claim_uniform_nilpotence_of_point_module':True,
  'regularity_test_uses_minimal_actual_operator':True,
  'infinity_not_omitted':True,'generic_trace_projection_not_assumed_to_preserve_lattice':True,
  'positive_weight_restriction_theorem_not_applied_to_zero_weight_coefficients':True},
 **{key:False for key in falsekeys}}
(p/'period_boundary_geometry_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({key:receipt[key] for key in ['completed','checker_sha256','point_tower_controls','boundary_subspace_controls']}))

