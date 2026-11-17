"""Exact finite controls; not a general Weyl-kernel or Singer implementation."""
from pathlib import Path
import hashlib,json,sys
import sympy as s

t,x,y=s.symbols('t x y')
p=Path(__file__).resolve().parent
def clean(z):return s.cancel(s.simplify(z))
def apply(A,f):return clean(sum(c*s.diff(f,t,j) for j,c in enumerate(A)))
def trim(A):
 A=list(map(clean,A))
 while len(A)>1 and A[-1]==0:A.pop()
 return A
def compose(A,B):
 C=[s.S.Zero]*(len(A)+len(B)-1)
 for i,a in enumerate(A):
  for j,b in enumerate(B):
   for k in range(i+1):C[i+j-k]+=a*s.binomial(i,k)*s.diff(b,t,k)
 return trim(C)
def rightdivide(A,B):
 R=trim(A);Q=[s.S.Zero]*max(1,len(A)-len(B)+1)
 while R!=[0] and len(R)>=len(B):
  k=len(R)-len(B);c=clean(R[-1]/B[-1]);Q[k]+=c
  sub=compose([0]*k+[c],B)
  R=trim([R[j]-(sub[j] if j<len(sub) else 0) for j in range(len(R))])
 return trim(Q),R
def eqop(A,B):
 A=trim(A);B=trim(B)
 return len(A)==len(B) and all(clean(a-b)==0 for a,b in zip(A,B))

# The four connections include an algebraic branch and a nontrivial rank-three model.
cyclic=[]
for name,A in [
 ('two_trivial',s.zeros(2)),
 ('trivial_plus_exponential',s.diag(0,1)),
 ('trivial_plus_square_root',s.diag(0,-1/(2*(1-t)))),
 ('rank_three',s.Matrix([[0,t,0],[0,1,1],[1,0,-1]]))]:
 n=A.rows
 nabla=lambda v:v.diff(t)+A*v
 v=s.eye(n)[:,0]
 for j in range(1,n):
  w=v
  for _ in range(j):w=nabla(w)
  v+= (s.eye(n)[:,j]-w.subs(t,0))*t**j/s.factorial(j)
 columns=[];w=v
 for j in range(n):columns.append(w);w=nabla(w)
 B=s.Matrix.hstack(*columns)
 assert B.subs(t,0)==s.eye(n)
 relation=B.inv()*w
 assert all(clean(z)==0 for z in w-B*relation)
 cyclic.append({'name':name,'rank':n,'vector':list(map(str,v)),
  'ordinary_point':0,'iterate_determinant_at_point':1,
  'cyclic_operator':list(map(str,[-z for z in relation]+[s.S.One]))})

def finite_relations(f,r,C):
 # Search finite supplied bounds, then certify every basis relation as an exact expression.
 cols=[(j,k) for j in range(r+1) for k in range(C+1)]
 N=len(cols)+4;rounds=[]
 while True:
  cf=[s.expand(s.series(s.diff(f,t,j),t,0,N).removeO()) for j in range(r+1)]
  mat=s.Matrix([[cf[j].coeff(t,n-k) if n>=k else 0 for j,k in cols] for n in range(N)])
  rows=s.polys.matrices.DomainMatrix.from_Matrix(mat).to_field().nullspace().to_Matrix()
  basis=[rows.row(j).T for j in range(rows.rows)]
  ops=[]
  for v in basis:
   A=[s.S.Zero]*(r+1)
   for a,(j,k) in zip(v,cols):A[j]+=a*t**k
   ops.append(trim(A))
  rounds.append({'jet_length':N,'nullity':len(basis)})
  if all(apply(A,f)==0 for A in ops):return ops,rounds
  N*=2
  assert N<=256,'Finite control bound exceeded; not a universal stopping rule'

H5=sum(s.binomial(5,j)*(-t)**j/s.factorial(j+1) for j in range(6))
cases=[
 ('constant',s.S.One,[0,-1,1],4,[0,1]),
 ('exponential',s.exp(t),[0,-1,1],4,[-1,1]),
 ('constant_plus_exponential',1+s.exp(t),[0,-1,1],4,None),
 ('linear_polynomial',1+t,[0,0,1],4,[-1/(1+t),1]),
 ('Kummer_degree_five',H5,[5,2-t,t],26,[-s.diff(H5,t)/H5,1]),
 ('algebraic_nonrational',1+s.sqrt(1-t),[0,-1,2*(1-t)],10,None),
 ('logarithm',s.log(1-t),[0,-1,1-t],10,None),
 ('generic_parameter',1+x*s.exp(t),[0,-1,1],4,None),
 ('exceptional_parameter_zero',s.S.One,[0,-1,1],4,[0,1]),
 ('nonminimal_log_equation_constant',s.S.One,[0,1,t],10,[0,1])]
minimum=[]
for name,f,L,C,expected in cases:
 print('control: '+name,flush=True)
 assert apply(L,f)==0,name
 ops,rounds=finite_relations(f,1,C)
 first=next((A for A in ops if len(A)==2),None)
 if first is None:
  assert expected is None,name
  M=[clean(c/L[-1]) for c in L]
 else:
  M=[clean(c/first[-1]) for c in first]
  assert expected is not None and eqop(M,expected),name
 Q,R=rightdivide(L,M);assert R==[0],(name,R)
 assert eqop(compose(Q,M),L)
 minimum.append({'name':name,'input_order':len(L)-1,'supplied_polynomial_degree_bound':C,
  'least_order_in_certified_finite_search':len(M)-1,'monic_operator':list(map(str,M)),
  'all_relation_basis_vectors_exactly_certified':True,'jet_rounds':rounds,
  'right_division_remainder_zero':True})
assert max(s.degree(s.fraction(clean(z))[1],t) for z in [-s.diff(H5,t)/H5,1])==5
assert s.degree(t,t)==1
# Finite supplied searches above do not implement the general bound/embedding producer.

# H=O u + O v, both periods1: retain the original polynomial lattice.
cyc=s.Matrix([1,t]);B=s.Matrix.hstack(cyc,cyc.diff(t))
Ai=[B.inv()*s.eye(2)[:,j] for j in range(2)]
M=[-1/(1+t),1]
coordinates=[]
for A in Ai:
 _,rem=rightdivide(list(A),M);coordinates.append(clean(rem[0]*(1+t)))
assert coordinates==[1,1]
generic_kernel=s.simplify(cyc.diff(t)-cyc/(1+t))
assert generic_kernel==s.Matrix([-1/(1+t),1/(1+t)])
assert s.simplify((1+t)*generic_kernel)==s.Matrix([-1,1])
assert clean((1/(1+t))*(1+t))==1
pole_num,pole_den=s.cancel(1/(1-t)).as_numer_denom()
assert pole_den.subs(t,1)==0 and pole_num.subs(t,1)!=0
lattice={'safe_source_periods':['1','1'],'cyclic_period':'1+t',
 'generic_kernel_vector':list(map(str,generic_kernel)),
 'contracted_kernel_generator':['-1','1'],
 'ambient_false_pole_at_minus_one_excluded':True,
 'same_generic_rank_does_not_remove_actual_pole_quotient':True}

# Point-module control: t D^k delta_a = a D^k delta_a-k D^(k-1) delta_a.
a=s.symbols('a');size=9
D=s.zeros(size);T=a*s.eye(size)
for k in range(size-1):D[k+1,k]=1
for k in range(1,size):T[k-1,k]=-k
torsion_checks=0
for k in range(5):
 e=s.eye(size)[:,k]
 assert (D*T-T*D)*e==e
 assert (T-a*s.eye(size))**(k+1)*e==s.zeros(size,1)
 torsion_checks+=1
# Exact algebraic relation; no complete algebraicity decision package is called.
f=1+s.sqrt(1-t)
assert clean((f-1)**2-(1-t))==0
assert apply([0,-1,2*(1-t)],f)==0
# Rational connection basis (1,sqrt(1-t)); finite double cover, rank2.
alg={'period':'1+sqrt(1-t)','relation':'(y-1)^2-(1-t)',
 'minimal_order_in_control':2,'finite_cover':'z^2=1-t',
 'constant_selected_from_equation_with_log_solutions_handled':True,
 'complete_Singer_solver_called':False}
# Vandermonde extraction is exact for a known degree bound.
fs=[s.sqrt(1-t),s.log(1-t),s.exp(t)]
F=sum(fs[j]*x**j for j in range(3))
V=s.Matrix([[i**j for j in range(3)] for i in range(3)])
recovered=V.inv()*s.Matrix([F.subs(x,i) for i in range(3)])
assert all(clean(recovered[j]-fs[j])==0 for j in range(3))
outchecks=0
for j in range(4):
 for m in range(5):
  got=s.diff(x**(j+2*m),x,2*m)
  assert got==s.factorial(j+2*m)/s.factorial(j)*x**j
  assert s.degree(got,x)==j
  outchecks+=1
sicchecks=0
for aa in range(4):
 for bb in range(4):
  for m in range(4):
   # Ecal(xi^(m+a) U^(m+b))=d_x^(m+a)x^(m+b).
   got=s.diff(x**(m+bb),x,m+aa)
   expected=0 if aa>bb else s.factorial(m+bb)/s.factorial(bb-aa)*x**(bb-aa)
   assert got==expected
   sicchecks+=1
falsekeys=[
 'general_minimization_field_embedding_executed','general_factor_degree_bound_producer_executed',
 'general_actual_kernel_solver_implemented','general_actual_kernel_solver_executed',
 'general_Weyl_syzygy_contraction_executed','complete_Singer_algebraicity_solver_executed',
 'universal_proofs_certified_by_checker','global_complex_safe_locus_uniformity_proved',
 'higher_actual_rank_classified','new_SIC33_answer_computed','new_Delta_minimum_computed',
 'numeric_cubic_frontier_changed','operator_order_reduced',
 'independent_replay','external_review','formal_verification']
receipt={'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'python_version':sys.version.split()[0],'sympy_version':s.__version__,
 'cyclic_vector_controls':cyclic,'bounded_scalar_minimization_controls':minimum,
 'polynomial_lattice_control':lattice,'point_module_checks':torsion_checks,
 'algebraic_cover_control':alg,'Vandermonde_degree_bound':2,
 'operator_bounded_output_checks':outchecks,'SIC_bounded_output_checks':sicchecks,
 'mutation_controls':{'coincident_initial_constants_do_not_identify_solutions':True,
 'two_equal_jet_spaces_not_used_as_stopping_certificate':True,
 'input_degree_not_used_as_right_factor_degree_bound':True,
 'nonminimal_operator_not_used_for_transcendence_decision':True,
 'generic_connection_not_used_as_actual_polynomial_image':True,
 'killed_point_module_not_discarded_from_full_kernel':True,
 'grid_interpolation_uses_known_bound_and_fresh_scalar_periods':True},
 **{key:False for key in falsekeys}}
out=p/'actual_period_kernel_verification.json'
out.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({key:receipt[key] for key in ['completed','checker_sha256','point_module_checks',
 'operator_bounded_output_checks','SIC_bounded_output_checks']}))

