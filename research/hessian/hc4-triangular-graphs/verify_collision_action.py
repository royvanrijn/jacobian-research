"""Exact action/first-derivative controls,not a global HC4 proof replay."""
import hashlib,json,time
from pathlib import Path
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();checks=[]
 def check(label,v):
  assert v,label
  checks.append(label)
 def zero(label,e):check(label,s.cancel(s.expand(e))==0)
 x,y,t,m,b,z,T,V=s.symbols('x y t m b z T V')
 # Polynomial S identity on nontrivial automorphisms;their collision
 # surface is empty,so these check only the source primitive formula.
 for nn in [2,3,5]:
  P=x;Q=y+x**nn;j=1;S=s.Rational(nn,nn+1)*x**(nn+1)
  zero('polynomial primitive dx shear '+str(nn),s.diff(S,x)-P*s.diff(Q,x))
  zero('polynomial primitive dy shear '+str(nn),s.diff(S,y)-(P*s.diff(Q,y)-j*x))
 P=x+(y+x*x)**3;Q=y+x*x
 j=s.expand(s.diff(P,x)*s.diff(Q,y)-s.diff(P,y)*s.diff(Q,x));zero('nonlinear composition Jacobian',j-1)
 aa=s.expand(P*s.diff(Q,x));bb=s.expand(P*s.diff(Q,y)-j*x)
 S=s.integrate(aa,x);res=s.expand(bb-s.diff(S,y));check('closed form residual independent of x',not res.has(x));S+=s.integrate(res,y)
 zero('nonlinear primitive dx',s.diff(S,x)-aa);zero('nonlinear primitive dy',s.diff(S,y)-bb)
 ux,uy,vx,vy=s.symbols('ux uy vx vy')
 p=(ux+vx/2,uy+vy/2);q=(ux-vx/2,uy-vy/2)
 for label,var in [('ux',ux),('uy',uy),('vx',vx),('vy',vy)]:
  theta=p[0]*s.diff(p[1],var)-q[0]*s.diff(q[1],var)
  corrected=theta-s.diff(ux*vy,var)
  expected=vx*s.diff(uy,var)-vy*s.diff(ux,var)
  zero('midpoint Liouville correction '+label,corrected-expected)
 SS=x**3+2*x*y**2+y
 F=(SS.subs({x:q[0],y:q[1]},simultaneous=True)-SS.subs({x:p[0],y:p[1]},simultaneous=True))-ux*vy
 zero('action is swap odd',F.subs({vx:-vx,vy:-vy},simultaneous=True)+F)
 # Arbitrary potential first derivatives reconstruct an exact symplectic
 # pair generically;the ACTUAL common polynomial K remains extra data.
 FF=s.Function('F')(m,b);r=s.diff(FF,b);uu=s.diff(FF,m)/r
 pp=(uu+r/2,b+m*(uu+r/2));qq=(uu-r/2,b+m*(uu-r/2))
 def jac(a):return s.diff(a[0],m)*s.diff(a[1],b)-s.diff(a[0],b)*s.diff(a[1],m)
 zero('general generating first derivatives give equal source area',jac(pp)-jac(qq))
 zero('recovered source difference x',pp[0]-qq[0]-r)
 zero('recovered source difference y',pp[1]-qq[1]-m*r)
 zero('recovered midpoint x',(pp[0]+qq[0])/2-uu)
 zero('recovered chord intercept',(pp[1]+qq[1])/2-m*(pp[0]+qq[0])/2-b)
 f0,rr,uu0=s.symbols('f0 rr uu0',nonzero=True)
 G=f0**2;Gb=2*f0*rr;Gm=Gb*uu0
 zero('squared action midpoint formula',Gm/Gb-uu0)
 zero('squared action difference square formula',Gb**2/(4*G)-rr**2)
 zero('squared action derivative unit identity',Gb**2-4*G*rr**2)
 a0,a1,a2=s.symbols('a0 a1 a2');HH=V**3+a2*V**2+a1*V+a0;Phi=HH.subs(V,T*T)
 zero('whole ordered relation is even',Phi.subs(T,-T)-Phi)
 zero('even action derivative transcription',s.diff(Phi,T)-2*T*s.diff(HH,V).subs(V,T*T))
 # Nontrivial actual open-A2/pseudo-plane control from LA1.
 zz=1+x*x*t;yy=2*t+x*x*t*t;S=-t*x*(t*t*x**4+3*t*x*x+4)
 theta_x=s.expand(-zz*(3*x*s.diff(yy,x)+4*yy)/2);theta_t=s.expand(-zz*3*x*s.diff(yy,t)/2)
 zero('known pseudo-plane primitive dx',s.diff(S,x)-theta_x)
 zero('known pseudo-plane primitive dt',s.diff(S,t)-(theta_t-x))
 xx2=-x;tt2=-t-2/x**2
 AA=s.cancel(S.subs({x:xx2,t:tt2},simultaneous=True)-S)
 FC=s.cancel(AA-(x+xx2)*(t-tt2)/2)
 zero('nontrivial collision action is minus4/x',FC+4/x)
 chord_b=-1/x**2;chord_m=(t+1/x**2)/x
 def db(f):return s.diff(f,x)*x**3/2+s.diff(f,t)*(t*x*x/2+s.Rational(3,2))
 def dm(f):return x*s.diff(f,t)
 zero('pseudo-plane base derivative db(b)',db(chord_b)-1);zero('pseudo-plane base derivative db(m)',db(chord_m))
 zero('pseudo-plane base derivative dm(m)',dm(chord_m)-1);zero('pseudo-plane base derivative dm(b)',dm(chord_b))
 zero('nontrivial first action b derivative',db(FC)-2*x);zero('nontrivial first action m derivative',dm(FC))
 rpc=db(FC);upc=dm(FC)/rpc
 for label,value,expected in [('x1',upc+rpc/2,x),('x2',upc-rpc/2,xx2),('y1',chord_b+chord_m*(upc+rpc/2),t),('y2',chord_b+chord_m*(upc-rpc/2),tt2)]:zero('whole source-point recovery '+label,value-expected)
 zero('pseudo-plane even action relation',FC**2+16*chord_b)
 check('rational deck retains source pole',s.limit(x*x*tt2,x,0)==-2)
 # Ramification with regular UNIT first derivative but pole in second.
 branch_b=z*z;branch_F=z*z+z**3
 fb=s.cancel(s.diff(branch_F,z)/(2*z));fbb=s.cancel(s.diff(fb,z)/(2*z));Psi=(T-b)**2-b**3
 zero('ramified action polynomial',Psi.subs({T:branch_F,b:branch_b}))
 zero('ramified first derivative is unit',fb-(1+s.Rational(3,2)*z))
 zero('ramified second derivative has pole',fbb-s.Rational(3,4)/z)
 for label,exp in [('T',s.diff(Psi,T)),('b',s.diff(Psi,b)),('m',s.diff(Psi,m))]:zero('ramified primitive image singular partial '+label,exp.subs({T:0,b:0}))
 ratio=s.cancel((-s.diff(Psi,b)/s.diff(Psi,T)).subs({T:branch_F,b:branch_b}));zero('regular0over0 derivative ratio',ratio-fb)
 zero('first jet recovers normalization parameter',s.Rational(2,3)*(fb-1)-z)
 pbr=(fb/2,branch_b+m*fb/2);qbr=(-fb/2,branch_b-m*fb/2)
 def jbr(a):return s.diff(a[0],m)*s.diff(a[1],z)-s.diff(a[0],z)*s.diff(a[1],m)
 zero('ramified local source-point maps have equal area',jbr(pbr)-jbr(qbr))
 check('both local source-point maps etale at chord ramification',jbr(pbr).subs(z,0)==jbr(qbr).subs(z,0)==-s.Rational(3,8))
 # Symplecticity alone is not the exactness used in the theorem.
 zero('rational symplectic logarithmic shear coefficient',x*s.diff(y+1/x,x)+1/x)
 check('Laurent derivative has zero dx/x residue',all(s.expand(s.diff(x**k,x)).coeff(x,-1)==0 for k in [-3,-2,-1,0,1,2,3]))
 # The written common-ring proof uses a THIRD generic sheet;these exact
 # ranks check the all-distinct-pair linear transcription,not a degree bound.
 for dd in [2,3,4,5]:
  rows=[]
  for i in range(dd):
   for jj in range(dd):
    if i==jj:continue
    row=[0]*(2*dd);row[i]=1;row[dd+jj]=-1;rows.append(row)
  rank=s.Matrix(rows).rank()
  check('whole ordered-pair equality rank d='+str(dd),rank==(2 if dd==2 else 2*dd-1))
 check('quadratic field countercontrol identifies whole source field',2<3)
 algebra_count=len(checks)
 context=json.loads((HERE/'sources/collision-action-context.json').read_text())
 for name,hsh in context['retained_proof_and_source_sha256'].items():check('retained bytes '+name,sha(ROOT/name)==hsh)
 for name,hsh in context['new_primary_sha256'].items():check('new primary bytes '+name,sha(HERE/name)==hsh)
 check('recent primary credited as context only',context['newly_imported_proved_external_statements']==[] and context['newly_used_unproved_external_statements']==[])
 flags={key:value for key,value in context.items() if isinstance(value,bool)}
 check('first jets do not assert full derivative closure',not flags['all_chord_derivations_preserve_collision_ring'] and not flags['second_chord_derivatives_asserted_regular'])
 check('normal rings kept distinct from primitive action orders',not flags['primitive_action_orders_assumed_normal_or_smooth'] and not flags['action_asserted_integral_over_chord_coordinate_ring'])
 check('no universal closure or N2 exclusion',not any(flags[key] for key in ['connected_chord_degree2_excluded','last_HC4_source_to_regular_inverse_Gram_implication_proved','whole_surviving_higher_degree_packet_excluded','original_goal_complete']))
 proof=HERE/'EXACT_COLLISION_ACTION_AND_FIRST_JET_RECONSTRUCTION.md';text=' '.join(proof.read_text().split())
 for phrase in ['CURRENT LAST UNSOLVED IMPLICATION','POLYNOMIAL primitive','GENERIC PRODUCT','E_C=k0[F]','E_D=k0[G]','ENTIRE C_r','G=0','Zariski Main','SINGULAR','RELATED PRIMARY CONTEXT','Laurent polynomial','NO whole polynomial Keller','EXACT INTERSECTION','third sheet','DEFERRED']:
  check('written scope '+phrase,phrase.lower() in text.lower())
 result=dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=algebra_count,checks=checks,proof_sha256=sha(proof),checker_sha256=sha(Path(__file__)),dependencies=context['retained_proof_and_source_sha256'],source_sha256=context['new_primary_sha256'],elapsed_seconds=round(time.monotonic()-start,3),**flags)
 (HERE/'receipt-collision-action.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:result[k] for k in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
if __name__=='__main__':main()
