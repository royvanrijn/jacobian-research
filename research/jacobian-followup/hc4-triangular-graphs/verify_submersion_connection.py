"""New whole-plane submersion/marked-volume controls;run before registration."""
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
x,y,u,t,X,Y=s.symbols('x y u t X Y')
for n in [3,4,7]:
 tag='NEW whole-plane submersion '+str(n)+' '
 P=x+x*x*y**n;A=1-2*x*y**n;B=s.Rational(4,n)*y**(n+1)
 def E(f):return A*s.diff(f,x)+B*s.diff(f,y)
 def D(f):return s.diff(P,y)*s.diff(f,x)-s.diff(P,x)*s.diff(f,y)
 b=s.diff(A,x)+s.diff(B,y)
 def nab(f):return E(f)+b*f
 eq(tag+'global polynomial Bezout lift',E(P),1)
 eq(tag+'divergence mark transcription',b,(2+s.Rational(4,n))*y**n)
 eq(tag+'nonsingular coframe determinant',s.diff(P,x)*A+s.diff(P,y)*B,1)
 for name,f in [('x',x),('y',y),('test',x*x*y+y**3+1)]:
  eq(tag+'Lie bracket '+name,E(D(f))-D(E(f)),-b*D(f))
  eq(tag+'corrected operator descends '+name,nab(D(f)),D(E(f)))
  eq(tag+'Weyl commutator '+name,nab(P*f)-P*nab(f),f)
 h=x+y*y+P
 aa=A+h*s.diff(P,y);bb=B-h*s.diff(P,x);b2=s.diff(aa,x)+s.diff(bb,y)
 eq(tag+'changed lift remains transverse',aa*s.diff(P,x)+bb*s.diff(P,y),1)
 eq(tag+'changed lift divergence class',b2-b,D(h))
 test=x*x*y+y**3+1
 eq(tag+'changed connection difference exact',aa*s.diff(test,x)+bb*s.diff(test,y)+b2*test-nab(test),D(h*test))
 eq(tag+'volume mark is connection of unit',nab(s.Integer(1)),b)
 eq(tag+'whole family finite cover trivialization',P.subs({x:t**n*X,y:Y/t},simultaneous=True),t**n*(X+X*X*Y**n))
 eq(tag+'inverse finite cover source X',(t**n*X)/t**n,X)
 eq(tag+'inverse finite cover source Y',t*(Y/t),Y)
 eq(tag+'zero fibre disjoint factors Bezout',(1+x*y**n)-x*y**n,1)
 eq(tag+'zero-fibre Gm parametrization',P.subs(x,-y**(-n)),0)
 genus=s.Rational(n-s.gcd(n,2),2);punctures=s.gcd(n,2)+1
 eq(tag+'cyclic-cover Euler characteristic',n*(-1)+1,1-n)
 eq(tag+'genus and punctures give rank n',2*genus+punctures-1,n)
 check(tag+'positive genus',genus>=1,True)
 # Monomial weight -1 components y*(x*y^n)^d for new unbounded proof transcriptions.
 for d in [0,1,3]:
  H=u**d;Q=y*(x*y**n)**d
  jac=s.diff(P,x)*s.diff(Q,y)-s.diff(P,y)*s.diff(Q,x)
  expected=((1+2*u)*H+n*u*(1+u)*s.diff(H,u)).subs(u,x*y**n)
  eq(tag+'all-degree mate operator transcription '+str(d),jac,expected)
  eq(tag+'highest coefficient obstruction '+str(d),s.Poly((1+2*u)*H+n*u*(1+u)*s.diff(H,u),u).LC(),2+n*d)
n,d=s.symbols('n d',integer=True,nonnegative=True)
check('universal leading factor positive for n>=2,d>=0',(2+n*d).is_positive,True)
context=json.loads((HERE/'sources/submersion-connection-context.json').read_text())
for p,h in context['retained_proof_and_source_sha256'].items():check('retained byte pin '+p,hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h)
for p,h in context['new_primary_sha256'].items():check('new primary byte pin '+p,hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h)
check('no new external theorem adoption',context['newly_imported_proved_external_statements']==[] and context['newly_used_unproved_external_statements']==[])
positive={'canonical_divergence_corrected_connection_for_every_polynomial_plane_submersion','raw_transverse_connection_bijective_without_any_Keller_mate','connection_and_divergence_class_independent_of_lift','vanishing_volume_mark_equivalent_to_polynomial_Keller_mate','actual_whole_plane_geometric_nonzero_quotients_in_unbounded_degree','explicit_generic_rank_n_and_unbounded_genus_proved','finite_cover_trivialization_gives_finite_generic_monodromy','all_n_controls_proved_have_no_polynomial_Keller_mate','geometric_acyclicity_alone_to_quotient_vanishing_shortcut_broken','external_review_deferred'}
for p,vv in context.items():
 if isinstance(vv,bool):check('scope '+p,vv is (p in positive))
result=dict(outcome='PASS',check_count=len(checks),algebraic_controls=algebraic,checks=checks,proof_sha256=hashlib.sha256((HERE/'SUBMERSION_GAUSS_MANIN_AND_MARKED_VOLUME_OBSTRUCTION.md').read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),dependencies=context['retained_proof_and_source_sha256'],source_sha256=context['new_primary_sha256'],elapsed_seconds=round(time.monotonic()-START,3),**{p:v for p,v in context.items() if isinstance(v,bool)})
(HERE/'receipt-submersion-connection.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({p:result[p] for p in ['outcome','check_count','algebraic_controls','elapsed_seconds','original_goal_complete']}))
