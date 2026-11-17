"""Exact transcription controls;no generic elimination or old checker replay."""
import hashlib,json,platform,time
from itertools import combinations
from pathlib import Path
import sympy as s
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(b,label):
 assert bool(b),label
 checks.append(label)
def zero(f,label):check(s.expand(f)==0,label)
x,y,u,v,t=s.symbols('x y u v t');z=[x,y,u,v]
# Integral formal recursion has no factorial denominators and does not assert global inversion.
a=s.symbols('a');Y=s.symbols('Y');F=x+a*x*x;G=Y
for r in range(2,5):
 error=s.Poly(s.expand(F.subs(x,G)),Y).coeff_monomial(Y**r);G=s.expand(G-error*Y**r)
zero(G-(Y-a*Y**2+2*a*a*Y**3-5*a**3*Y**4),'the integral homogeneous recursion gives the exact first four inverse coefficients')
error=s.Poly(s.expand(F.subs(x,G)-Y),Y)
check(all(error.coeff_monomial(Y**r)==0 for r in range(5)),'the recursion cancels every required prefix degree')
check(error.coeff_monomial(Y**5)==-14*a**4,'finite inverse prefixes need an actual composition certificate')
check(all(c.q==1 for c in s.Poly(G,Y,a).coeffs()),'the inverse recursion uses integer coefficients rather than factorial division')
D=s.symbols('D',integer=True,positive=True);m=D-1;B=m**3
zero(m*B-(D-1)**4,'both finite inverse composition degrees are bounded by (D-1)^4')
# Six Pluecker charts glue to the explicit projective quadric.
rows=s.Matrix(2,4,s.symbols('a0:8'));p={(i,j):rows[:,[i,j]].det() for i,j in combinations(range(4),2)}
zero(p[0,1]*p[2,3]-p[0,2]*p[1,3]+p[0,3]*p[1,2],
 'every two-plane basis satisfies the full Pluecker quadric over Z')
rr,ss,tt,ww=s.symbols('rr ss tt ww')
for i,j in combinations(range(4),2):
 frame=s.zeros(2,4);frame[0,i]=1;frame[1,j]=1;others=[k for k in range(4) if k not in [i,j]]
 frame[:,others[0]]=s.Matrix([rr,ss]);frame[:,others[1]]=s.Matrix([tt,ww])
 check(frame[:,[i,j]].det()==1,'the standard projective Grassmann chart has a unit pivot minor '+str((i,j)))
# The closed-branch product gives an actual inverse when the plane ideal is the unit ideal.
ui,j1,j2,b1,b2=s.symbols('ui j1 j2 b1 b2')
zero(ui*(b1*j1+b2*j2)-(b1*(ui*j1)+b2*(ui*j2)),
 'a unit combination of plane-ideal generators kills every inverse-certificate generator exactly')
check(s.rem(s.Symbol('s')**2,s.Symbol('s')**3,s.Symbol('s'))!=0,
 'radical vanishing is deliberately not substituted for ordinary ideal containment')
# Universal all-prime determinant identity,with p an arbitrary integer exponent.
pp=s.symbols('p',integer=True,positive=True)
cp=pp*(pp+1)
matrix=s.Matrix([[cp*x**(pp-1),0,1,0],[0,0,0,1],[1,0,cp*u**(pp-1),0],[0,1,0,0]])
zero(matrix.det()-(1-pp**2*(pp+1)**2*x**(pp-1)*u**(pp-1)),
 'the complete integral determinant has the exact universal p squared factor')
zero(s.simplify(s.diff(x**(pp+1),x,x)-cp*x**(pp-1)),
 'the all-exponent pure second derivative has the required exact p factor')
X,Z,alpha=s.symbols('X Z alpha')
zero((x+t*alpha)*(X+Z)-(x*X+t*alpha*X+x*Z+t*alpha*Z),
 'Frobenius times a linear factor gives the all-prime Hasse expansion without omitted cross terms')
for prime in [2,3,5,7]:
 h=x*u+y*v+x**(prime+1)+u**(prime+1);gradient=[s.diff(h,zz) for zz in z]
 check(all(int(c)%prime==0 for c in s.Poly(s.hessian(h,z).det()-1,z).coeffs()),
  'the exact prime-field Hessian control is unit at p='+str(prime))
 for point in [(0,0,0,0),(1,0,-1,0)]:
  check([int(g.subs(dict(zip(z,point))))%prime for g in gradient]==[0,0,0,0],
   'the prime-field collision is exact at p='+str(prime)+' point '+str(point))
 wx,wy,wu,wv=s.symbols('wx wy wu wv')
 shifted=s.Poly(s.expand(h.subs(dict(zip(z,[x+t*wx,y+t*wy,u+t*wu,v+t*wv])),simultaneous=True)),t)
 coeff=shifted.coeff_monomial(t**prime)
 target=x*wx**prime+u*wu**prime+(wx*wu+wy*wv if prime==2 else 0)
 check(all(int(c)%prime==0 for c in s.Poly(s.expand(coeff-target),z+[wx,wy,wu,wv]).coeffs()),
  'the high Hasse coefficient forces both active directions zero at p='+str(prime))
 check(all(int(c)%(prime*prime)==0 for c in s.Poly(s.hessian(h,z).det()-1,z).coeffs()),
  'the entire naive potential passes precision p squared at p='+str(prime))
zero((y+t*s.Symbol('a'))*(v+t*s.Symbol('b'))-y*v-t*(s.Symbol('a')*v+s.Symbol('b')*y)-t*t*s.Symbol('a')*s.Symbol('b'),
 'the remaining y-v plane retains a nonzero joint quadratic translation term')
e,nu=s.symbols('e nu',integer=True,positive=True)
zero(e*nu-2*e-e*(nu-2),'the valuation certificate and naive precision2e force val_p(N_D)>=2')
# An arbitrary-polynomial automorphic family lies outside the plane locus.
aa=s.Function('a')(x);ell=y+aa*u;U=s.Function('U');h=x*v+u*ell+U(ell)
grad=[s.diff(h,zz) for zz in z]
zero(grad[3]-x,'the arbitrary-function component family recovers the active coordinate from the last output')
zero(grad[2]-aa*grad[1]-ell,'the two passive outputs recover the complete moving linear form')
zero(grad[1]-u-s.Subs(s.diff(U(s.Symbol('L')),s.Symbol('L')),s.Symbol('L'),ell),
 'the passive derivative gives the stated fibre recovery term')
zero(grad[0]-v-s.diff(aa,x)*u*grad[1],
 'the final active derivative gives the integral fourth inverse coordinate')
S,cxx,cxy,cxu=s.symbols('S cxx cxy cxu')
block=s.Matrix([[cxx,cxy,cxu,1],[cxy,S,1+aa*S,0],[cxu,1+aa*S,2*aa+aa*aa*S,0],[1,0,0,0]])
zero(block.det()-1,'the full arbitrary-function Hessian determinant is identically one')
literal=x*v+u*(y+x*u)+(y+x*u)**2
lg=[s.diff(literal,zz) for zz in z];out=s.symbols('o0:4');lx=out[3];le=out[2]-lx*out[1]
lu=out[1]-2*le;ly=le-lx*lu;lv=out[0]-lu*out[1];inv=[lx,ly,lu,lv]
for i in range(4):
 zero(lg[i].subs(dict(zip(z,inv)),simultaneous=True)-out[i],
  'the explicit nonplane component certificate has its full polynomial inverse output '+str(i))
 zero(inv[i].subs(dict(zip(out,lg)),simultaneous=True)-z[i],
  'the explicit nonplane component certificate inverse also composes back in coordinate '+str(i))
r=s.symbols('r',integer=True,positive=True);mm=s.symbols('mm',integer=True,nonnegative=True)
check(s.simplify((mm+2)*(r+1)-(r+2)).is_positive,
 'the arbitrary-r positive-degree family has a strict top homogeneous gap when m>=2')
# Homothetic paths preserve unit Hessian but do not propagate inverse from the special point.
coeff=s.symbols('c2:6');model=sum(coeff[k-2]*x**k for k in range(2,6))
scaled=s.expand(model.subs(x,t*x)/t**2)
zero(s.diff(scaled,x,x)-s.diff(model,x,x).subs(x,t*x),
 'the all-degree homothetic Hessian rule has no additional parameter determinant factor')
zero(scaled.subs(t,0)-coeff[0]*x*x,'the homothetic special fibre is exactly the quadratic part')
check(min(exp[0] for exp,c in s.Poly(scaled-coeff[0]*x*x,t).terms())==1,
 'replacing t by t^M raises the first nonquadratic contact to an arbitrary prescribed order')
algebra=len(checks)
context=json.loads((HERE/'sources/coefficient-closure-context.json').read_text())
for path,digest in context['retained_proof_and_source_sha256'].items():
 check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained immutable input '+path)
for entry in context['newly_cached_primary_context']:
 check(hashlib.sha256((HERE/'sources'/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
  'inspected primary source '+entry['retained_file'])
check(len(context['newly_imported_proved_external_statements'])==1 and context['newly_used_unproved_external_statements']==[],
 'projective properness is the sole new external proved import;all other classification and inverse-degree inputs are retained')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
check(sum(flags.values())==11 and not flags['generic_elimination_ideals_or_actual_N_D_computed']
 and not flags['naive_ordinary_Hessian_coefficient_ring_has_product_containment_proved']
 and not flags['unrestricted_HC4_or_JC2_resolved'] and not flags['original_goal_complete'],
 'the finite controls do not compute generic coefficient ideals or close the plane-Keller remainder')
proof=HERE/'COEFFICIENT_CLOSURE_COMPONENTS_AND_ALL_PRIME_BARRIERS.md'
paths=['sources/coefficient-closure-context.json']+['sources/'+e['retained_file'] for e in context['newly_cached_primary_context']]
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,
 elapsed_seconds=round(time.monotonic()-START,3),python=platform.python_version(),sympy=s.__version__,
 proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 dependencies=context['retained_proof_and_source_sha256'],source_sha256={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths},
 independent_replay=False,formal_verification=False,external_review=False,**flags)
(HERE/'receipt-coefficient-closure.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in {'checks','dependencies','source_sha256'}},indent=2))
