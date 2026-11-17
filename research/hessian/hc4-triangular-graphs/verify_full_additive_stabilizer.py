"""New positive-boundary affine-source controls only;run before registration."""
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
x,f,u,v=s.symbols('x f u v')
for a,q in [(2,2),(3,2),(4,3),(7,2)]:
 tag='NEW affine positive-boundary chart '+str((a,q))+' '
 N=(q+a+q-1)//q;M=q*N-q-a
 check(tag+'nonnegative source exponent',0<=M<q,True)
 G=x*f**q-1;H=x*G;T=x**N*f**M*G
 sw={x:u**a*v+u**(-q),f:u}
 g=u**(q+a)*v;h=u**a*v*(1+g);t=v*(1+g)**N
 for name,source,other in [('G',G,g),('H',H,h),('T',T,t)]:eq(tag+'ring generator '+name,source.subs(sw),other)
 eq(tag+'two opens cover',(G+1)-G,1)
 eq(tag+'source reconstruction H/G',H/G,x)
 eq(tag+'second chart reconstruction T/(G+1)^N',T/(G+1)**N,(x*f**q-1)/f**(q+a))
 eq(tag+'second chart f reconstruction',f.subs(sw),u)
 eq(tag+'retained source line disjoint from added line',G.subs(f,0),-1)
 eq(tag+'added line free coordinate',t.subs(u,0),v)
 eq(tag+'other boundary equations G',g.subs(u,0),0)
 eq(tag+'other boundary equations H',h.subs(u,0),0)
 # Omega=dx wedge df has dv wedge du coefficient u^a.
 eq(tag+'exact volume order coefficient',s.diff(sw[x],v),u**a)
 KU=x*f/s.Integer(1-q);KW=u**(a+1)*v/s.Integer(1-q)
 Udx=-s.diff(KU,x);Udf=x-s.diff(KU,f)
 Wu=u**a*v-s.diff(KW,u);Wv=-s.diff(KW,v)
 eq(tag+'primitive du compatibility',Udx.subs(sw)*s.diff(sw[x],u)+Udf.subs(sw),Wu)
 eq(tag+'primitive dv compatibility',Udx.subs(sw)*s.diff(sw[x],v),Wv)
 eq(tag+'source primitive differential',s.diff(Udf,x)-s.diff(Udx,f),1)
 eq(tag+'second primitive differential',s.diff(Wu,v)-s.diff(Wv,u),u**a)
 for name,source,other in [('f',f,u),('G',G,g),('H',H,h),('T',T,t)]:
  eq(tag+'minimal additive field compatibility '+name,(f**a*s.diff(source,x)).subs(sw),s.diff(other,v))
  eq(tag+'source finite nilpotence '+name,s.diff(source,x,N+2),0)
 eq(tag+'minimal field translates added line',s.diff(t,v).subs(u,0),1)
 eq(tag+'one-below-threshold pole residue',s.cancel(u*(f**(a-1)*s.diff(T,x)).subs(sw)).subs(u,0),1)
 eq(tag+'minimal field fixes retained source line',(f**a).subs(f,0),0)
 # New ideal's multiplication/commutation,not old KX numerical threshold tests.
 p=f**a*(1+f);r=f**a*(2+f**2);test=x**3+f*x**2
 eq(tag+'all source replicas commute transcription',p*s.diff(r*s.diff(test,x),x)-r*s.diff(p*s.diff(test,x),x),0)
context=json.loads((HERE/'sources/full-additive-stabilizer-context.json').read_text())
for p,h in context['retained_proof_and_source_sha256'].items():check('retained byte pin '+p,hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h)
check('no new external theorem adoption',context['newly_imported_proved_external_statements']==[] and context['newly_used_unproved_external_statements']==[])
positive={'proper_affine_plane_enlargement_has_at_most_one_regular_coordinate_pencil','every_extending_source_LND_produces_integral_source_coordinate','all_extending_source_LNDs_classified_by_one_replica_ideal','finite_trace_dual_and_full_normalization_source_stabilizers_equal','two_distinct_extending_source_additive_fields_force_full_Keller_inverse','nonzero_full_normalization_source_additive_stabilizer_commutative','nonzero_additive_branch_source_automorphism_stabilizer_triangular','whole_plane_exact_volume_positive_canonical_boundary_controls_all_a_q','external_review_deferred'}
for p,vv in context.items():
 if isinstance(vv,bool):check('scope '+p,vv is (p in positive))
result=dict(outcome='PASS',check_count=len(checks),algebraic_controls=algebraic,checks=checks,proof_sha256=hashlib.sha256((HERE/'FULL_NORMALIZATION_ADDITIVE_PENCIL_AND_TRACE_DUAL_STABILIZER.md').read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),dependencies=context['retained_proof_and_source_sha256'],source_sha256=context['new_primary_sha256'],elapsed_seconds=round(time.monotonic()-START,3),**{p:v for p,v in context.items() if isinstance(v,bool)})
(HERE/'receipt-full-additive-stabilizer.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({p:result[p] for p in ['outcome','check_count','algebraic_controls','elapsed_seconds','original_goal_complete']}))
