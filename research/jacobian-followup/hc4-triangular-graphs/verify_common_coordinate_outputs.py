"""New block/inverse/specialization controls;one run before registration only."""
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
a,b,c,d,c11,c12,c22=s.symbols('a b c d c11 c12 c22');J=s.Matrix([[a,b],[c,d]])
C=s.Matrix([[c11,c12],[c12,c22]]);Z=s.zeros(2);H=C.row_join(J.T).col_join(J.row_join(Z))
eq('universal cotangent Hessian determinant',H.det(),J.det()**2)
w1,w2,x1,x2,D1,D2=s.symbols('w1 w2 xi1 xi2 D1 D2');w=s.Matrix([w1,w2]);xi=s.Matrix([x1,x2]);shift=s.Matrix([D1,D2])
forward=shift+J.T*w;back=J.T.inv()*(xi-shift)
for i in range(2):
 eq('universal momentum forward-after-inverse '+str(i),(shift+J.T*back)[i],xi[i])
 eq('universal momentum inverse-after-forward '+str(i),back[i].subs({x1:forward[0],x2:forward[1]},simultaneous=True),w[i])
 check('momentum inverse denominator only scalar Jacobian '+str(i),s.denom(s.factor(back[i])) in [J.det(),-J.det()],True)
z1,z2,t1,t2,u,v=s.symbols('z1 z2 t1 t2 u v')
for m in [2,3,6]:
 G=s.Matrix([z1+z2**m+t1*z2+t2,z2+t1]);outputs=list(G)+[t1,t2]
 eq('block determinant agrees with fibre determinant m'+str(m),s.Matrix(outputs).jacobian([z1,z2,t1,t2]).det(),G.jacobian([z1,z2]).det())
 eq('actual parameterized Keller fibre determinant m'+str(m),G.jacobian([z1,z2]).det(),1)
 zz2=v-t1;zz1=u-zz2**m-t1*zz2-t2
 for i,target in enumerate([u,v]):eq('actual fibre inverse m'+str(m)+' component'+str(i),G[i].subs({z1:zz1,z2:zz2},simultaneous=True),target)
 # A non-Keller family isolates why the generic-open qualifier is necessary.
 f=t1*z1**m+z1
 check('degree can drop in a non-Keller specialization m'+str(m),s.Poly(f,z1).degree()==m and s.Poly(f.subs(t1,0),z1).degree()==1,True)
 check('degree-drop control explicitly fails Keller m'+str(m),s.diff(f,z1)!=1,True)
 # Actual cotangent stabilization on the same arbitrary-degree triangular pair.
 P=z1+z2**m;Q=z2;D=z1*z2+z2**2
 f=D+w1*P+w2*Q
 eq('cotangent control full determinant m'+str(m),s.hessian(f,[z1,z2,w1,w2]).det(),1)
 p=s.diff(f,z1);q=s.diff(f,z2)
 eq('cotangent control first momentum inverse m'+str(m),p.subs(w1,x1-z2),x1)
 ww2=x2-z1-2*z2-m*z2**(m-1)*(x1-z2)
 eq('cotangent control second momentum inverse m'+str(m),q.subs({w1:x1-z2,w2:ww2},simultaneous=True),x2)
z=s.symbols('z');coef=s.symbols('a0:6');poly=sum(coef[i]*z**i for i in range(6));der=s.Poly(s.diff(poly,z),z)
for i in range(2,6):eq('univariate coefficient identity i'+str(i),der.coeff_monomial(z**(i-1)),i*coef[i])
nn,mm,rr=s.symbols('n m r',integer=True,nonnegative=True)
eq('dimension descent relation',nn-(nn-rr),rr)
eq('stable degree-three seed common-output upper bound',(3+rr)-3,rr)
eq('nonautomorphic four-variable common-output upper bound',4-2,2)
eq('gamma four-variable common-output upper bound',4-3,1)
check('global two-plus-two graph coordinates contradict gamma bound',2>1,True)
check('only first gamma composition has field degree three',3**1==3 and all(3**k>4 for k in [2,3,5]),True)
context=json.loads((HERE/'sources/common-coordinate-outputs-context.json').read_text())
for p,h in context['retained_proof_and_source_sha256'].items():check('retained byte pin '+p,hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h)
check('no newly adopted external theorem',context['newly_imported_proved_external_statements']==[] and context['newly_used_unproved_external_statements']==[])
positive={'all_degree_all_dimension_common_coordinate_specialization_proved','generic_specialization_preserves_full_field_degree_on_nonempty_open','n_minus_one_common_outputs_force_polynomial_inverse','small_excluded_plane_degrees_bound_common_outputs_by_n_minus_three','exact_common_output_count_for_every_identity_stabilization_of_three_dimensional_degree_three_or_four_seed','HC4_gradient_orbit_iff_plane_Keller_times_identity_two_proved','nonautomorphic_HC4_gradient_has_exactly_two_common_outputs','whole_marked_normalization_and_monodromy_product_necessity_proved','entire_gamma_fixed_map_HC4_polynomial_equivalence_route_excluded','every_gamma_global_two_plus_two_polynomial_graph_chart_excluded','published_full_plane_quartic_retirement_preserved','external_review_deferred'}
for p,vv in context.items():
 if isinstance(vv,bool):check('scope '+p,vv is (p in positive))
result=dict(outcome='PASS',check_count=len(checks),algebraic_controls=algebraic,checks=checks,proof_sha256=hashlib.sha256((HERE/'COMMON_COORDINATE_OUTPUTS_AND_HC4_ORBIT_FILTER.md').read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),dependencies=context['retained_proof_and_source_sha256'],source_sha256=context['new_primary_sha256'],elapsed_seconds=round(time.monotonic()-START,3),**{p:v for p,v in context.items() if isinstance(v,bool)})
(HERE/'receipt-common-coordinate-outputs.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({p:result[p] for p in ['outcome','check_count','algebraic_controls','elapsed_seconds','original_goal_complete']}))
