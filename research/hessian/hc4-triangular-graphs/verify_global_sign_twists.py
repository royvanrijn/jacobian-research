"""Exact identity controls for the global all-degree sign-twist proof."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):check(sp.cancel(value)==0,label)
def zeros(values,label):check(all(sp.cancel(v)==0 for v in values),label)
p,X,Y,z,W,x,y,w,P,Q=sp.symbols('p X Y z W x y w P Q')
E12=sp.Matrix([[1,p],[0,1]]);E21=sp.Matrix([[1,0],[-p,1]])
M=sp.Matrix([[1-p*p,2*p-p**3],[-p,1-p*p]])
zeros(list(M-E12*E21*E12),'universal matrix is a product of elementary polynomial automorphisms')
zero(M.det()-1,'universal determinant is literally one')
check(M.subs(p,0)==sp.eye(2),'untwisted root retains the ordered hyperbolic pair')
check(M.subs(p,1)==sp.Matrix([[0,1],[-1,0]]),'twisted root reverses the hyperbolic product sign')
zeros(list(M*M.adjugate()-sp.eye(2)),'universal inverse is a polynomial matrix')
old=M*sp.Matrix([X,Y]);q=sp.expand(2*old[0]*old[1]-2*(1-2*p)*X*Y)
quotient=2*((p+1)*X*X+(2*p*p+2*p-2)*X*Y+(p**3+p*p-2*p-2)*Y*Y)
zero(q-p*(p-1)*quotient,'universal product discrepancy has the exact p(p-1) factor')
zero(sp.rem(q,p*(p-1),p),'divisibility holds identically before any interpolation or degree choice')
zeros([q.subs(p,0),q.subs(p,1)],'both root choices vanish with arbitrary hyperbolic variables')
C=z*z-1;pp=(1-z)/2;MM=M.subs(p,pp)
A=sp.cancel(q.subs(p,pp)/C)
check(sp.Poly(A,X,Y,z).is_multivariate,'the two-root correction is an actual polynomial')
old=MM*sp.Matrix([X,Y]);psi=[old[0],old[1],z,W-A]
check(sp.Matrix(psi).jacobian([X,Y,z,W]).det().expand()==1,'full four-coordinate source Jacobian is one')
recover=MM.adjugate()*sp.Matrix([x,y])
inverse=[recover[0],recover[1],z,w+A.subs({X:recover[0],Y:recover[1]},simultaneous=True)]
zeros([v.subs(dict(zip([x,y,z,w],psi)),simultaneous=True)-r for v,r in zip(inverse,[X,Y,z,W])],
      'the complete polynomial inverse cancels the forward map')
zeros([v.subs(dict(zip([X,Y,z,W],inverse)),simultaneous=True)-r for v,r in zip(psi,[x,y,z,w])],
      'the forward map cancels the complete polynomial inverse')
d=sp.Function('D')(z)
h=2*x*y+C*w+d;target=2*z*X*Y+C*W+d
zero(h.subs(dict(zip([x,y,z,w],psi)),simultaneous=True)-target,
     'the actual source pullback gives C_D for an arbitrary symbolic D(z)')
zero(2*z*X*(z/2)-X*C-X,'the uniform coefficient pair generates the full two-point ideal')
f=sp.Matrix([2*z*X,C]);J=f.jacobian([X,z]).det()
zero(J-4*z*z,'the twisted coefficient Jacobian is globally nonconstant')
zero(sp.hessian(h,[x,y,z,w]).det()-16*z*z,'native h_D Hessian determinant is independent of arbitrary D')
zero(sp.hessian(target,[X,Y,z,W]).det()-16*z**4,'twisted native Hessian is not mistaken for a unit')
for sign in [-1,1]:
    point={x:0,y:0,z:sign,w:-sp.diff(d,z).subs(z,sign)/(2*sign)}
    zeros([sp.diff(h,r).subs(point,simultaneous=True).doit() for r in [x,y,z,w]],
          'the complete h_D critical equations vanish at the stated point '+str(sign))
    zero(sp.hessian(h,[x,y,z,w]).det().subs(point,simultaneous=True)-16,
         'the full critical point is reduced Morse with determinant 16 '+str(sign))
    zero(J.subs({X:0,z:sign})-4,'both unsquared coefficient orientations are literally four '+str(sign))
    zero(h.subs(point,simultaneous=True)-d.subs(z,sign),'critical values retain the arbitrary D data '+str(sign))
K2=w*(2*x*y+z*z-1)+z
for sign in [-1,1]:
    point={x:0,y:0,z:sign,w:-sp.Rational(sign,2)}
    zeros([sp.diff(K2,r).subs(point) for r in [x,y,z,w]],'K2 has the stated critical point '+str(sign))
    zero(sp.hessian(K2,[x,y,z,w]).det().subs(point)-4,'K2 two-point critical determinant '+str(sign))
    zero(K2.subs(point)-sign,'K2 critical value '+str(sign))
# Universal determinant repair is independent of degrees and number of ideal generators' terms.
a,b,c,e,r,s,f1,f2=sp.symbols('a b c e r s f1 f2')
N=sp.Matrix([[a,b],[c,e]]);ff=sp.Matrix([f1,f2]);kk=N*ff
repair=N+sp.Matrix([r,s])*sp.Matrix([[-f2,f1]])
zeros(list(repair*ff-kk),'determinant repair preserves the entire coefficient pair')
zero(repair.det()-N.det()-s*kk[0]+r*kk[1],
     'determinant repair changes det by exactly s k1-r k2')
unit=sp.symbols('unit',nonzero=True)
zero((repair.det()-unit).subs(unit,N.det()+s*kk[0]-r*kk[1]),
     'a common residue orientation gives the exact target unit determinant')
pair=sp.Matrix([P,Q]);trans=N.T*pair
zero((ff.T*trans)[0]-(kk.T*pair)[0],'the transpose fibre change gives the Keller coefficient pair with correct ordering')
zeros(list(N.T.adjugate()*N.T*pair-N.det()*pair),'the fibre map has its actual adjugate-over-unit inverse')
# Generic symmetric cotangent Hessian. D affects only the top-left block.
t11,t12,t22,j11,j12,j21,j22=sp.symbols('t11 t12 t22 j11 j12 j21 j22')
H=sp.Matrix([[t11,t12,j11,j21],[t12,t22,j12,j22],[j11,j12,0,0],[j21,j22,0,0]])
zero(H.det()-(j11*j22-j12*j21)**2,'arbitrary cotangent base block has determinant exactly J(K)^2')
check(sp.hessian(2*x*y+z*z+w*w,[x,y,z,w]).det()!=0,
      'the nondegenerate quadratic endpoint has native unit Hessian and cannot be replaced by a no-unit requirement')
algebra=len(checks);context=json.loads((HERE/'sources/global-sign-twists-context.json').read_text())
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained immutable input '+path)
check(context['newly_imported_proved_external_statements']==[] and context['newly_used_unproved_external_statements']==[],
      'no new external theorem or unproved classification is imported')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
check(sum(flags.values())==5 and not flags['fixed_K2_full_source_orbit_resolved']
      and not flags['two_point_Keller_fibre_constructed'] and not flags['two_point_Keller_fibre_excluded']
      and not flags['h0_identified_with_h1_as_functions'] and not flags['original_goal_complete'],
      'exact existence equivalences and universal redundancy do not claim existence,exclusion or source equivalence')
proof=HERE/'GLOBAL_SIGN_TWISTS_AND_TWO_POINT_KELLER_EQUIVALENCE.md'
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,
 elapsed_seconds=round(time.monotonic()-START,3),python=platform.python_version(),sympy=sp.__version__,
 proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 dependencies=context['retained_proof_and_source_sha256'],
 source_sha256={'sources/global-sign-twists-context.json':hashlib.sha256((HERE/'sources/global-sign-twists-context.json').read_bytes()).hexdigest()},
 independent_replay=False,formal_verification=False,external_review=False,**flags)
(HERE/'receipt-global-sign-twists.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in {'checks','dependencies','source_sha256'}},indent=2))
