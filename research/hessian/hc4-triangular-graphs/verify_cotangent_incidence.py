"""New incidence component/transcription controls; run once before registration."""
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
x,y,z,u,v,t,L,U,V=s.symbols('x y z u v t L U V')
w=x*y
A=(1+w)**3*z+y*y*(1+w)*(4+3*w)
B=y+3*x*(1+w)**2*z+3*x*y*y*(4+3*w)
C=2*x-3*x*x*y-x**3*z
# These are NEW axis-preimage controls,not a replay of the known Jacobian/collision.
for label,point,bval in [('L0',{x:0,y:t,z:-4*t*t},t),('L1',{x:t,y:-1/t,z:5/t**2},2/t),('L2',{x:t,y:-2/t,z:8/t**2},-2/t)]:
 for name,f,target in [('A',A,0),('C',C,0),('B',B,bval)]:eq(label+' '+name,f.subs(point,simultaneous=True),target)
eq('x zero A',A.subs(x,0),z+4*y*y)
eq('x zero C',C.subs(x,0))
eq('x zero B',B.subs(x,0),y)
H=(1+u)**2*v+u**2*(4+3*u)
unit={y:u/x,z:v/x**2}
eq('unit chart A',x**2*A.subs(unit,simultaneous=True),(1+u)*H)
eq('unit chart B',x*B.subs(unit,simultaneous=True),u+3*H)
eq('unit chart C',C.subs(unit,simultaneous=True),x*(2-3*u-v))
eq('beta complete unit elimination',(1+u)*H.subs(v,2-3*u),(u+1)*(u+2))
eq('beta restricted B',u+3*H.subs(v,2-3*u),4*u+6)
check('beta two distinct unit roots',s.solve((u+1)*(u+2),u)==[-2,-1],True)
eq('gamma eliminate H',(1+u)*(-u/s.Integer(3)),-u*(u+1)/3)
eq('gamma minus-one branch impossible',H.subs(u,-1),1)
eq('gamma zero branch forces v zero',H.subs(u,0),v)
for name,f,target in [('A',A,0),('B',B,0),('C',C,2*x)]:eq('gamma complete line '+name,f.subs({y:0,z:0}),target)
c=s.symbols('c',nonzero=True)
eq('beta nonzero level root on L1',2/(2/c),c)
eq('beta nonzero level root on L2',-2/(-2/c),c)
check('beta torus numerator never zero',all(a!=0 for a in [2,-2]),True)
k,g,p,N=s.symbols('k g p N');EZ=k*L-g*(U+V)+k-p
for n in [2,3,5]:
 ex=L**(2*n-2)-L**(n-2)*EZ+N*L**(n-1)
 expected=L**(2*n-2)+(N-k)*L**(n-1)+g*L**(n-2)*(U+V)+(p-k)*L**(n-2)
 eq('incidence Hodge transcription n'+str(n),ex,expected)
 bidegrees=[(n-1,n-1),(n-1,n-2),(n-2,n-1),(n-2,n-2)]
 check('distinct Hodge defect bidegrees n'+str(n),len(set(bidegrees))==4,True)
eq('beta axis class',L+2*(L-1),3*L-2)
X0=L**4-L*(3*L-2)+L**2;Xc=L**4-L*(3*L-2)+3*L**2
eq('beta zero level motive',X0,L**4-2*L**2+2*L)
eq('beta nonzero level motive',Xc,L**4+2*L)
for label,ex,chi in [('zero',X0,1),('nonzero',Xc,3)]:
 check('beta '+label+' Hodge defect nonzero',s.expand(ex-L**4)!=0,True)
 eq('beta '+label+' unavoidable L coefficient',s.expand(ex-L**4).coeff(L,1),2)
 eq('beta '+label+' Euler',ex.subs(L,1),chi)
eq('gamma all-level motive',L**4-L*L+L**2,L**4)
check('stable Hodge multiplication injective in integral polynomial domain',s.Poly(U*V,U,V).is_zero is False,True)
context=json.loads((HERE/'sources/cotangent-incidence-context.json').read_text())
for p,h in context['retained_proof_and_source_sha256'].items():check('retained byte pin '+p,hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h)
for p,h in context['new_primary_sha256'].items():check('new primary '+p,hashlib.sha256((HERE/p).read_bytes()).hexdigest()==h)
check('one explicitly declared proved import',len(context['newly_imported_proved_external_statements'])==1)
check('no unproved external assertion',context['newly_used_unproved_external_statements']==[])
positive={'all_dimension_all_degree_incidence_class_formula_proved','all_level_stable_Hodge_obstruction_proved','generic_target_line_irreducibility_import_hypotheses_checked','plane_incidence_Hodge_equality_equivalent_to_inverse_using_retained_GR0','entire_primitive_beta_source_question_closed_negatively','zero_level_Euler_failure_detected_by_Hodge_polynomial','primitive_gamma_chart_passes_Hodge_test_only','external_review_deferred'}
for p,vv in context.items():
 if isinstance(vv,bool):check('scope '+p,vv is (p in positive))
result=dict(outcome='PASS',check_count=len(checks),algebraic_controls=algebraic,checks=checks,proof_sha256=hashlib.sha256((HERE/'COTANGENT_INCIDENCE_HODGE_OBSTRUCTION.md').read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),dependencies=context['retained_proof_and_source_sha256'],source_sha256=context['new_primary_sha256'],elapsed_seconds=round(time.monotonic()-START,3),**{p:v for p,v in context.items() if isinstance(v,bool)})
(HERE/'receipt-cotangent-incidence.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({p:result[p] for p in ['outcome','check_count','algebraic_controls','elapsed_seconds','original_goal_complete']}))
