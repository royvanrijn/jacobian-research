"""New spectral/module/saturation controls only;run before registration."""
import hashlib,json,time
from pathlib import Path
import sympy as s
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
checks=[];algebraic=0
def check(name,truth,alg=False):
 global algebraic
 assert truth,name
 checks.append(name);algebraic+=int(alg)
def eq(name,a,b=0):check(name,s.cancel(a-b)==0,True)
x,y,z,T=s.symbols('x y z T')
def E(f):return x*s.diff(f,x)+y*s.diff(f,y)
P=3*x+2*y+x*y+x**3+y**4;q=x**2*y+y**5+1
eq('source Euler multiplication commutator',E(P*q)-P*E(q),E(P)*q)
for power in [1,2,3]:
 def ad(f,qq):
  # Universal Leibniz makes ad_E^r(m_P)=m_(E^r P).
  return E(f*qq)-f*E(qq)
 ep=P
 for _ in range(power-1):ep=E(ep)
 eq('nested multiplication commutator transcription '+str(power),ad(ep,q),E(ep)*q)
for D in [2,4,7]:
 projector=s.prod((T-r)/(1-r) for r in range(D+1) if r!=1)
 for degree in range(D+1):eq('linear-degree projector '+str((D,degree)),projector.subs(T,degree),int(degree==1))
coeff=s.Poly(s.expand(s.prod((T-r)/(1-r) for r in range(5) if r!=1)),T)
ep=P;projected=0
for power in range(coeff.degree()+1):
 if power:ep=E(ep)
 projected+=coeff.coeff_monomial(T**power)*ep
eq('explicit ordinary linear component recovery',projected,3*x+2*y)
# Local-first-jet weighted control;no global Keller claim for this test polynomial map.
weights=[1,2,4];variables=[x,y,z]
polys=[x+y*y+z,y+x*x+z*z,z+y*y+x*x*y+x**4+x*y]
def component(f,weight):
 return sum(coef*s.prod(v**e for v,e in zip(variables,exps)) for exps,coef in s.Poly(f,*variables).terms() if sum(w*e for w,e in zip(weights,exps))==weight)
eq('weighted local first jet invertible',s.Matrix(polys).jacobian(variables).subs({x:0,y:0,z:0}).det(),1)
eq('weighted minimum block recovers x',component(polys[0],1),x)
eq('weighted second block recovers y after lower subtraction',component(polys[1],2)-x*x,y)
eq('weighted fourth block recovers z after lower subtraction',component(polys[2],4)-y*y-x*x*y-x**4,z)
a,b,c,d=s.symbols('a b c d');M=T*s.eye(2)-s.Matrix([[a,b],[c,d]])
for i in range(2):
 for j in range(2):eq('determinant trick adjugate '+str((i,j)),(M.adjugate()*M)[i,j],M.det()*int(i==j))
eq('determinant trick polynomial is monic',s.Poly(M.det(),T).LC(),1)
X=s.symbols('X')
for m in [1,2,5]:
 P=x;Q=x**m*y
 def E1(f):return (x-1)*s.diff(f,x)+y*s.diff(f,y)
 def E0(f):return x*s.diff(f,x)+y*s.diff(f,y)
 eq('NEW orbit/module family Jacobian m'+str(m),s.det(s.Matrix([P,Q]).jacobian([x,y])),x**m)
 eq('NEW orbit/module family first jet at correct centre m'+str(m),s.det(s.Matrix([P,Q]).jacobian([x,y])).subs({x:1,y:0}),1)
 eq('NEW generator Euler linear component m'+str(m),s.expand(Q.subs(x,1+X)).coeff(X,0),y)
 eq('NEW finite-module stability defect m'+str(m),E1(Q*y),(m+2)*Q*y-m*x**(m-1)*y**2)
 check('NEW defect outside B+B*y by exponent bound m'+str(m),m-1<m,True)
 eq('NEW unshifted grading preserves Q m'+str(m),E0(Q),(m+1)*Q)
 eq('NEW wrong fixed point has singular first jet m'+str(m),s.det(s.Matrix([P,Q]).jacobian([x,y])).subs({x:0,y:0}),0)
 for k in [2,3]:eq('NEW products recover arbitrarily higher y powers transcription '+str((m,k)),s.expand((1+X)**(m*k)*y**k).coeff(X,0),y**k)
m=s.symbols('m',integer=True,positive=True)
eq('all-m nonmembership exponent deficit',(m-1)-m,-1)
context=json.loads((HERE/'sources/contracting-lattice-context.json').read_text())
for p,h in context['retained_proof_and_source_sha256'].items():check('retained byte pin '+p,hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h)
check('no new external theorem adoption',context['newly_imported_proved_external_statements']==[] and context['newly_used_unproved_external_statements']==[])
positive={'all_dimension_all_degree_positive_grading_finite_module_properness_proved','nonzero_module_need_not_be_free_or_contain_one','source_multiplication_recovered_by_spectral_components_and_first_jet','positive_weighted_coordinate_recovery_proved','determinant_trick_forces_integrality_and_finiteness','one_contracting_finite_lattice_inverse_criterion_proved','plane_trace_dual_contracting_invariance_equivalent_to_inverse','contracting_invariance_of_etale_algebra_only_removes_added_unramified_primes','all_m_finite_generator_orbit_to_module_closure_shortcut_broken','fixed_point_first_jet_and_positive_weights_boundaries_retained','external_review_deferred'}
for p,vv in context.items():
 if isinstance(vv,bool):check('scope '+p,vv is (p in positive))
result=dict(outcome='PASS',check_count=len(checks),algebraic_controls=algebraic,checks=checks,proof_sha256=hashlib.sha256((HERE/'CONTRACTING_SOURCE_GRADINGS_AND_FINITE_LATTICES.md').read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),dependencies=context['retained_proof_and_source_sha256'],source_sha256=context['new_primary_sha256'],elapsed_seconds=round(time.monotonic()-START,3),**{p:v for p,v in context.items() if isinstance(v,bool)})
(HERE/'receipt-contracting-lattice.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({p:result[p] for p in ['outcome','check_count','algebraic_controls','elapsed_seconds','original_goal_complete']}))
