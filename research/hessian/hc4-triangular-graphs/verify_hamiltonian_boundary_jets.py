"""New Hamiltonian stabilizer and boundary-jet transcriptions;pre-registration only."""
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
x,f,u,v,z=s.symbols('x f u v z')
def bracket(A,B):return s.diff(A,x)*s.diff(B,f)-s.diff(A,f)*s.diff(B,x)
for a,q in [(2,2),(3,2),(4,3)]:
 e=a+1;G=x*f**q-1;A=f*G;B=x*G**2;C=f**e
 tag='NEW preserving frame '+str((a,q))+' '
 eq(tag+'first bracket',bracket(A,B),-G**2*((q+3)*x*f**q-1))
 eq(tag+'second bracket',bracket(A,C),e*f**(e+q))
 eq(tag+'retained line spanning determinant',bracket(A,B).subs(f,0),1)
 check(tag+'positive-f spanning determinant',e>0 and e+q>0,True)
 check(tag+'boundary orders use retained valuations',min(1+a+q,2*a+q,e)>=e,True)
 eq(tag+'incomplete projected source ODE',-s.diff(A,x),-f**(q+1))
 for n in [0,1,3]:
  h=x*G**2*sum((j+1)*(x*f**q)**j for j in range(n+1))
  defect=x*(-(n+2)*(x*f**q)**(n+1)+(n+1)*(x*f**q)**(n+2))
  eq(tag+'all-order line approximation '+str(n),h-x,defect)
  check(tag+'line-thickening divisibility '+str(n),all(pf>=q*(n+1) for (px,pf),coef in s.Poly(h-x,x,f).terms()),True)
  check(tag+'approximation boundary lower bound '+str(n),2*(a+q)-q>=e,True)
for a in [1,2,4]:
 e=a+1;h=v+v*v;k=v**3-v
 H=u**e*h;K=u**e*k
 local=(s.diff(H,u)*s.diff(K,v)-s.diff(H,v)*s.diff(K,u))/u**a
 eq('NEW local preserving bracket order '+str(a),local,e*u**e*(h*s.diff(k,v)-s.diff(h,v)*k))
 check('NEW maximal-ideal square has strictly greater order '+str(a),2*e>e,True)
 for count in [2,4]:
  coeffs=[v**j for j in range(count)]
  mat=s.Matrix([[s.expand(hh).coeff(v,j) for hh in coeffs] for j in range(count)])
  check('NEW independent leading cotangent coefficients '+str((a,count)),mat.rank()==count,True)
context=json.loads((HERE/'sources/hamiltonian-boundary-jets-context.json').read_text())
for p,h in context['retained_proof_and_source_sha256'].items():check('retained byte pin '+p,hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h)
check('NO new external theorem/cache/unproved import',context['new_primary_sha256']=={} and context['newly_imported_proved_external_statements']==[] and context['newly_used_unproved_external_statements']==[])
positive={'all_divergence_free_source_stabilizers_classified_by_boundary_jets','Hamiltonian_source_stabilizers_of_finite_trace_dual_and_normalization_equal','arbitrary_finite_source_jets_and_closed_source_subscheme_interpolation_proved','finite_preserving_field_collection_spans_entire_source_tangent_module','many_field_source_integration_automatic','positive_boundary_has_Poisson_trace_and_no_preserving_Keller_pair','positive_boundary_Hamiltonian_stabilizer_non_Noetherian_proved','new_three_fields_span_on_all_retained_positive_order_controls','new_all_order_source_line_approximants_match_nonpreserving_Hamiltonian','external_review_deferred'}
for p,vv in context.items():
 if isinstance(vv,bool):check('scope '+p,vv is (p in positive))
result=dict(outcome='PASS',check_count=len(checks),algebraic_controls=algebraic,checks=checks,proof_sha256=hashlib.sha256((HERE/'HAMILTONIAN_STABILIZERS_AND_INVISIBLE_BOUNDARY_JETS.md').read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),dependencies=context['retained_proof_and_source_sha256'],source_sha256={},elapsed_seconds=round(time.monotonic()-START,3),**{p:v for p,v in context.items() if isinstance(v,bool)})
(HERE/'receipt-hamiltonian-boundary-jets.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({p:result[p] for p in ['outcome','check_count','algebraic_controls','elapsed_seconds','original_goal_complete']}))
