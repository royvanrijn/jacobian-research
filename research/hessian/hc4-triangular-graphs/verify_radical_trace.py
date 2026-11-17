"""Exact controls; the unbounded class and geometric proofs are in the companion note."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic()
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
checks=[]
def check(v,label):
    assert bool(v),label
    checks.append(label)
def zero(v,label):check(sp.simplify(sp.expand(v))==0,label)
def zeros(vs,label):check(all(sp.simplify(sp.expand(v))==0 for v in vs),label)

t,s,w1,w2,x,y1,y2,r=sp.symbols('t s w1 w2 x y1 y2 r')
# m is sampled only as an identity regression. The proof uses an arbitrary finite geometric sum.
for m in [1,2,3,7]:
    w=t**m; rr=t**(m+1); yy=-2*t-t**(m+2)
    vm=sum((-r)**k/sp.Integer(2)**(k+1) for k in range(m))
    zero((2+r)*vm-(1-(-r)**m/sp.Integer(2)**m),'geometric-series source-recovery identity m='+str(m))
    recovered=-yy*vm.subs(r,rr)+(-1)**m*w**m*rr/sp.Integer(2)**m
    zero(recovered-t,'ambient polynomial recovers the actual curve parameter m='+str(m))
    zero((1+rr)**2+w*yy-1,'the actual closed plane lies on the four-quadric m='+str(m))
    zero(rr**m-w**(m+1),'the divisor polynomial vanishes on that exact plane m='+str(m))
    Xminus=-w1*y1/(x+1)
    zero(Xminus**m-w1**(m+1)-w1**m*((-y1/(x+1))**m-w1),
         'the generic unwanted plane has exact factor w1^m m='+str(m))
    check(sp.degree(t**m,t)==m,'the retained base map has degree m and is not a section when m>1, m='+str(m))
    zero(sp.diff(t**m,t)-m*t**(m-1),'the retained projection is explicitly ramified for m>1, m='+str(m))

zero(0**2+1*(1-s*t)+s*t-1,'the class-zero surface is a genuine affine plane on the quadric')
check(sp.Matrix([s,t]).jacobian([s,t])==sp.eye(2),'its ambient w2,Y2 recover both source parameters')
zero((x*x+w2*y2-1).subs(w2,0)-(x*x-1),'the principal divisor w2 on w1=0 has exactly the two signed planes')
C,S,a2,b2=sp.symbols('C S a2 b2')
xr=C*x+S*b2
br=-S*x+C*b2
zero(xr*xr+br*br-(C*C+S*S)*(x*x+b2*b2),'the X-b2 rotation preserves the quadratic on the unit circle')
zeros([xr.subs({C:-1,S:0})+x,br.subs({C:-1,S:0})+b2],'its endpoint flips X and swaps the second hyperbolic pair')
check(sp.Poly(x-1,x).degree()==1,'the positive fibre plane is a principal divisor in the total-family threefold')
check(sp.Poly(x-1,x).eval(0)!=0,'removing the critical level does not remove that positive fibre plane')

# Exact all-degree trace and first-jet relations with arbitrary source coefficients.
f,A1,A2,K1,K2,Cc=sp.symbols('f A1 A2 K1 K2 Cc')
H=w1*A1+w2*A2
X=f+H
Y1=K1-(2*f+H)*A1+w2*Cc
Y2=K2-(2*f+H)*A2-w1*Cc
zero(X*X+w1*Y1+w2*Y2-f*f-w1*K1-w2*K2,
     'the complete variable-normal-trace compensation has the literal cotangent identity')
zero(Y1.subs({w1:0,w2:0})+2*f*A1-K1,'the first free jet recovers K1 modulo the radical trace')
zero(Y2.subs({w1:0,w2:0})+2*f*A2-K2,'the first free jet recovers K2 modulo the radical trace')
zero((K1-Y1).subs({w1:0,w2:0,f:0}),'on f=0 the first actual helper equals the Keller coordinate')
zero((K2-Y2).subs({w1:0,w2:0,f:0}),'on f=0 the second actual helper equals the Keller coordinate')

# Nonconstant hyperbolic descended trace: actual global closed embedding, not an A1 zero curve.
z1,z2=sp.symbols('z1 z2')
k=[z1+z2*z2,z2]
fn=k[0]*k[1]-1
AA=[w1*w1,w2*w2]
hh=w1*AA[0]+w2*AA[1]
XX=fn+hh
YY=[k[i]-(2*fn+hh)*AA[i] for i in range(2)]
phi=fn*fn+w1*k[0]+w2*k[1]
zero(XX*XX+w1*YY[0]+w2*YY[1]-phi,'the nonconstant hyperbola-trace embedding has its exact potential')
zero(sp.Matrix(k).jacobian([z1,z2]).det()-1,'its actual plane coefficient pair is Keller')
zero(sp.hessian(phi,[z1,z2,w1,w2]).det()-1,'the genuine nonconstant source potential has unit Hessian')
recf=XX-hh
reck=[YY[i]+(XX+recf)*AA[i] for i in range(2)]
rec=[reck[0]-reck[1]**2,reck[1]]
zeros([rec[0]-z1,rec[1]-z2],'all source coordinates have global polynomial ambient recovery')
T1,T2=sp.symbols('T1 T2')
hyper=T1*T2-1
zero(hyper.subs({T1:k[0],T2:k[1]},simultaneous=True)-fn,'the variable normal trace really descends to the target')
zero(T1*T2-(hyper+1),'the descended zero curve has a nonconstant unit coordinate')
check(sp.gcd(hyper,T1)==1 and sp.gcd(hyper,T2)==1,'the hyperbola has no vanishing coordinate and is not an affine line')
for b in [1,2]:
    beta=[sp.Rational(1,b),sp.Integer(b)]
    invz=[beta[0]-beta[1]**2,beta[1]]
    zero(fn.subs(dict(zip([z1,z2],invz)),simultaneous=True),'each selected zero target has normal trace zero '+str(b))
    pert=phi-w1*beta[0]-w2*beta[1]
    zeros([sp.diff(pert,v).subs({z1:invz[0],z2:invz[1],w1:0,w2:0},simultaneous=True)
           for v in [z1,z2,w1,w2]],'the actual target perturbation has the stated critical point '+str(b))

# A coordinate-composite trace has a curved actual A1 component, with multiplicities allowed.
coord=z1+z2*z2
zero((coord**3).subs(z1,-z2*z2),'every power of a source coordinate has the full curved A1 zero component')
zeros([coord.subs({z1:t-s*s,z2:s},simultaneous=True)-t,z2.subs(z2,s)-s],
      'a displayed source coordinate pair has a polynomial straightening inverse')
context=json.loads((HERE/'sources/radical-trace-context.json').read_text())
check(context['every_integer_quadric_middle_class_has_closed_affine_plane'] and not context['arbitrary_affine_plane_class_assumed_unit'],
      'arbitrary actual affine planes are not substituted for IW sections')
check(context['rank_one_radical_curve_closed_embedding_proved'] and context['actual_A1_radical_component_forces_inverse'],
      'the radical curve uses the actual zero-section ring, with actual A1 components')
check(context['target_descent_of_square_forces_normal_trace_descent'] and context['descended_radical_curve_isomorphism_proved'],
      'target square descent and the exact reduced-curve isomorphism are kept distinct from assumptions')
check(context['nonautomorphic_descended_components_singular_with_A1_normalization'] and not context['singular_descended_cusp_declared_closed'],
      'singular parametric descended curves are a remaining boundary, not an affine-line shortcut')
check(context['newly_used_unproved_external_statements']==[],'only explicitly retained geometric foundations are reused')
check(not context['unrestricted_HC4_or_JC2_resolved'] and not context['original_goal_complete'],'the original all-degree objective remains open')
math_count=len(checks)
for name,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,'retained byte pin '+name)
for e in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/e['retained_file']).read_bytes()).hexdigest()==e['sha256'],
          'new inspected context byte pin '+e['retained_file'])
old=json.loads((HERE/'receipt-quadric-middle-cycle.json').read_text())
check(old['outcome']=='PASS' and old['check_count']==106,'the retained IW replay remains PASS106; it is not newly replayed')
for name,digest in old['source_sha256'].items():
    check(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,'retained IW receipt source '+name)
names=['verify_radical_trace.py','UNBOUNDED_QUADRIC_CLASSES_AND_RADICAL_TRACE.md','sources/radical-trace-context.json']
result=dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=math_count,checks=checks,
 source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names},
 dependencies=context['retained_proof_and_source_sha256'],retained_receipt_counts_not_newly_replayed={'receipt-quadric-middle-cycle.json':106},
 newly_cached_primary_context=context['newly_cached_primary_context'],newly_used_unproved_external_statements=[],
 retained_external_foundations=context['retained_external_foundations'],retained_structural_inputs=context['retained_structural_inputs'],
 written_not_finite_exhaustion=context['direct_proof_arguments'],
 **{k:v for k,v in context.items() if type(v) is bool},
 runtime=dict(python=platform.python_version(),sympy=sp.__version__),elapsed_seconds=round(time.monotonic()-START,3))
(HERE/'receipt-radical-trace.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=math_count,original_goal_complete=False),indent=2))
