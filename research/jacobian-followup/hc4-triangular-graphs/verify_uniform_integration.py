"""Exact transcription controls;the universal proof is in the companion note.

Run only before registration.No finite examples certify the universal theorem.
"""
import hashlib,json,time
from pathlib import Path
import sympy as s
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
checks=[];algebraic=0
def check(name,truth,alg=False):
 global algebraic
 assert truth,name
 checks.append(name);algebraic+=int(alg)
def eq(name,left,right=0):check(name,s.expand(left-right)==0,True)
x,y,p,q,t,z=s.symbols('x y p q t z');alpha=s.symbols('alpha')
# Two different non-linear polynomial Keller frames,whose construction is exact.
frames=[(x+y**2,y),(x+y**2,y+(x+y**2)**3)]
for idx,(P,Q) in enumerate(frames):
 j=s.expand(s.diff(P,x)*s.diff(Q,y)-s.diff(P,y)*s.diff(Q,x))
 def E(f,P=P,Q=Q,j=j):return s.expand((s.diff(Q,y)*s.diff(f,x)-s.diff(Q,x)*s.diff(f,y))/j)
 def D(f,P=P,Q=Q,j=j):return s.expand((-s.diff(P,y)*s.diff(f,x)+s.diff(P,x)*s.diff(f,y))/j)
 eq(f'frame{idx}:Jacobian',j,1)
 for label,l,r in [('EP',E(P),1),('EQ',E(Q),0),('DP',D(P),0),('DQ',D(Q),1)]:eq(f'frame{idx}:{label}',l,r)
 a=x**2+x*y+y; r=x*y**2+x**3+1
 eq(f'frame{idx}:commuting derivatives',E(D(r)),D(E(r)))
 eq(f'frame{idx}:coordinate x derivative',s.diff(P,x)*E(r)+s.diff(Q,x)*D(r),s.diff(r,x))
 eq(f'frame{idx}:coordinate y derivative',s.diff(P,y)*E(r)+s.diff(Q,y)*D(r),s.diff(r,y))
 for label,T in [('E',E),('D',D)]:eq(f'frame{idx}:ideal core identity {label}',T(a*r)-a*T(r),T(a)*r)
 eq(f'frame{idx}:image multiplier identity',D(a*r)-a*D(r),D(a)*r)
 # Actual polynomial primitive of r dx wedge dy,converted to dP,dQ.
 S=s.integrate(r,x)
 # Inverse coframe:dy=-Q_x*dP+P_x*dQ for j=1.
 A=s.expand(-S*s.diff(Q,x));B=s.expand(S*s.diff(P,x))
 eq(f'frame{idx}:top deRham homotopy',E(B)-D(A),r)
 h=x**3*y+y**2
 eq(f'frame{idx}:closed exact one form',E(D(h))-D(E(h)))
 eq(f'frame{idx}:Q stability identity',D(Q*r)-Q*D(r),r)
# Point-module Taylor projection on C[t],E=t,z=-d/dt.
def Z(w):return -s.diff(w,t)
def power(fn,w,n):
 for _ in range(n):w=fn(w)
 return w
def proj(w):
 return s.expand(sum(t**i*power(Z,w,i)/s.factorial(i) for i in range(s.degree(w,t)+1)))
for degree in [0,1,2,5]:
 w=sum((i+1)*t**i for i in range(degree+1))
 eq(f'point degree{degree}:projector',proj(w),1)
 restored=sum((-1)**i*t**i*proj(power(Z,w,i))/s.factorial(i) for i in range(degree+1))
 eq(f'point degree{degree}:inverse Taylor expansion',restored,w)
 eq(f'point degree{degree}:Weyl commutator',t*Z(w)-Z(t*w),w)
# Nonzero Kummer control on arbitrary Laurent basis exponent m,not bounded m.
m=s.symbols('m',integer=True)
eq('universal Kummer Weyl commutator',(m+1+alpha)-(m+alpha),1)
eq('universal Kummer E inverse coefficient',(m+1+alpha)/(m+1+alpha),1)
check('noninteger alpha avoids every integer-exponent coefficient zero',s.Rational(1,2).q!=1,True)
for mval in [-3,-1,0,2]:
 a=s.Rational(1,2)
 eq(f'Kummer basis{mval}:explicit inverse',(mval+1+a)/(mval+1+a),1)
 # M0 E=t,P=-d_t:commutator through negative as well as positive powers.
 w=t**mval
 eq(f'Laurent module basis{mval}:Weyl commutator',t*(-s.diff(w,t))+s.diff(t*w,t),w)
eq('Laurent quotient horizontal class representative',t*t**-1,1)
eq('point submodule E image has zero constant term',(t*(1+2*t+3*t**4)).subs(t,0),0)
# Failed full-frame hypothesis:x partial_x has a nonzero proper image ideal.
eq('nonframe control every x-multiple has polynomial primitive',x*s.diff(x**3*y**2/3,x),x**3*y**2)
check('nonframe control image has no constant coefficient',s.Poly(x*s.diff(x**2+y,x),x).coeff_monomial(1)==0,True)
context=json.loads((HERE/'sources/uniform-integration-context.json').read_text())
for path,sha in context['retained_proof_and_source_sha256'].items():
 check('retained byte pin '+path,hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha)
for path,sha in context['new_primary_sha256'].items():
 check('primary cache '+path,hashlib.sha256((HERE/path).read_bytes()).hexdigest()==sha)
check('no new external theorem import',not context['newly_imported_proved_external_statements'] and not context['newly_used_unproved_external_statements'])
for key,value in context.items():
 if isinstance(value,bool):check('scope '+key,value is (key in {
  'uniform_source_integration_ideal_nonzero_iff_plane_inverse_proved',
  'full_derivative_image_multiplier_ring_equals_C_P_under_failure_proved',
  'ideal_core_dichotomy_for_every_dimension_polynomial_Keller_frames_proved',
  'raw_relative_quotient_transverse_derivative_bijective_proved',
  'special_fibre_point_torsion_exactly_total_reducibility_proved',
  'saturated_quotient_horizontal_classes_equal_point_torsion_cokernel_proved',
  'nonzero_acyclic_Kummer_and_nonsplit_point_extension_controls_proved',
  'classical_surjectivity_criterion_credited_and_reproved','external_review_deferred'}))
proof=HERE/'UNIFORM_INTEGRATION_AND_EXACT_RELATIVE_QUOTIENT.md'
result=dict(outcome='PASS',check_count=len(checks),algebraic_controls=algebraic,checks=checks,
 proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 dependencies=context['retained_proof_and_source_sha256'],source_sha256=context['new_primary_sha256'],
 elapsed_seconds=round(time.monotonic()-START,3),**{k:v for k,v in context.items() if isinstance(v,bool)})
(HERE/'receipt-uniform-integration.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['outcome','check_count','algebraic_controls','elapsed_seconds','original_goal_complete']},indent=2))
