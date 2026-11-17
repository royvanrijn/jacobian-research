"""New balanced exact-volume controls;pre-registration only."""
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
x,f,u,v,q,p=s.symbols('x f u v q p')
def delta(h):return s.expand(x*s.diff(h,x)-f*s.diff(h,f))
for a in [2,3,5]:
 G=x*f-1;H=x*G;T=x**(a+1)*G;xx=u**a*v+u**-1;Gu=u**(a+1)*v
 tag='NEW balanced q1 '+str(a)+' '
 eq(tag+'G chart',G.subs({x:xx,f:u},simultaneous=True),Gu)
 eq(tag+'H chart',H.subs({x:xx,f:u},simultaneous=True),u**a*v*(1+Gu))
 eq(tag+'T chart',T.subs({x:xx,f:u},simultaneous=True),v*(1+Gu)**(a+1))
 eq(tag+'source x on D G',H/G,x)
 eq(tag+'source x on D f',(G+1)/f,x)
 eq(tag+'second v recovery',v*(1+Gu)**(a+1)/(1+Gu)**(a+1),v)
 eq(tag+'two principal opens cover',(G+1)-G,1)
 eq(tag+'source opens cover',x*f-G,1)
 eq(tag+'volume order',s.diff(xx,v),u**a)
 for name,h,w in [('f',f,-1),('G',G,0),('H',H,1),('T',T,a+1)]:
  eq(tag+'stable generator weight '+name,delta(h),w*h)
 eq(tag+'moment source fixed value',(G+1).subs({x:0,f:0}),0)
 eq(tag+'moment added fixed value',(1+Gu).subs({u:0,v:0}),1)
 eq(tag+'primitive mismatch',xx-u**a*v,1/u)
 eq(tag+'source primitive moment constant',x*(-f)+x*f,0)
 eq(tag+'boundary primitive moment constant',u**a*v*(-u)+(1+Gu),1)
 eq(tag+'balanced volume weight',-a+(a+1)-1,0)
# The nonaffine sharpness control:two cotangent charts,not a source reconstruction.
vv=1/q;ww=-q*q*p
eq('NEW cotangent chart volume',s.det(s.Matrix([[s.diff(vv,q),s.diff(vv,p)],[s.diff(ww,q),s.diff(ww,p)]])),1)
eq('NEW cotangent primitives glue dq coefficient',-ww*s.diff(vv,q),-p)
eq('NEW cotangent moment glues',-vv*ww,q*p)
eq('NEW cotangent v weight',q*s.diff(vv,q)-p*s.diff(vv,p),-vv)
eq('NEW cotangent w weight',q*s.diff(ww,q)-p*s.diff(ww,p),ww)
context=json.loads((HERE/'sources/balanced-torus-context.json').read_text())
for pp,h in context['retained_proof_and_source_sha256'].items():check('retained byte pin '+pp,hashlib.sha256((ROOT/pp).read_bytes()).hexdigest()==h)
check('NO new external theorem/cache/unproved import',context['new_primary_sha256']=={} and context['newly_imported_proved_external_statements']==[] and context['newly_used_unproved_external_statements']==[])
positive={'exact_affine_balanced_source_torus_enlargement_rigidity_proved','one_balanced_source_field_preserving_full_normalization_forces_inverse','one_balanced_source_field_preserving_finite_trace_dual_forces_inverse','no_target_field_descent_assumed','global_primitive_forces_all_fixed_moment_values_zero','new_q1_balanced_whole_plane_control_is_nonexact','nonaffine_exact_cotangent_control_marks_affineness_boundary','external_review_deferred'}
for pp,vv in context.items():
 if isinstance(vv,bool):check('scope '+pp,vv is (pp in positive))
result=dict(outcome='PASS',check_count=len(checks),algebraic_controls=algebraic,checks=checks,proof_sha256=hashlib.sha256((HERE/'EXACT_VOLUME_AND_BALANCED_SOURCE_TORUS_RIGIDITY.md').read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),dependencies=context['retained_proof_and_source_sha256'],source_sha256={},elapsed_seconds=round(time.monotonic()-START,3),**{pp:vv for pp,vv in context.items() if isinstance(vv,bool)})
(HERE/'receipt-balanced-torus.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({pp:result[pp] for pp in ['outcome','check_count','algebraic_controls','elapsed_seconds','original_goal_complete']}))
