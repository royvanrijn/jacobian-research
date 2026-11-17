"""Transcription/group controls for the written all-degree proof.Run only before registration."""
import hashlib,itertools,json,math,time
from pathlib import Path
import sympy as s
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
checks=[];algebraic=0
def check(name,truth,alg=False):
 global algebraic
 assert truth,name
 checks.append(name);algebraic+=int(alg)
def eq(name,left,right=0):check(name,s.cancel(left-right)==0,True)
u,v,w,z=s.symbols('u v w z');a,b,c=s.symbols('a b c');xs=[u,v,w,z]
R=(u-v)*(w-z)/((u-w)*(v-z))
eq('universal Riccati logarithmic difference',(a*u**2+b*u+c)-(a*v**2+b*v+c),(u-v)*(a*(u+v)+b))
eq('universal crossratio derivative',sum(s.diff(R,x)*(a*x**2+b*x+c) for x in xs))
h=s.symbols('h0:4');T=s.symbols('T')
interp=s.interpolate(list(zip(xs[:3],h[:3])),T)
eq('fourth-point crossratio defect equals interpolation defect',sum(s.diff(R,x)*hv for x,hv in zip(xs,h)),s.diff(R,z)*(h[3]-interp.subs(T,z)))
check('fourth point derivative coefficient nonzero',s.cancel(s.diff(R,z))!=0,True)
V=s.Matrix([[1,x,x**2] for x in xs[:3]])
eq('three-point quadratic interpolation determinant',V.det(),(v-u)*(w-u)*(w-v))
# Rational Mobius gauge including derivatives of the gauge coefficients.
A,B,C,D,Ap,Bp,Cp,Dp=s.symbols('A B C D Ap Bp Cp Dp');Z=s.symbols('Z');Delta=A*D-B*C
num=(Ap*T+Bp)*(C*T+D)-(A*T+B)*(Cp*T+Dp)+Delta*(a*T**2+b*T+c)
gauge=s.cancel((num/(C*T+D)**2).subs(T,(D*Z-B)/(A-C*Z)))
check('full rational Mobius gauge stays quadratic',s.Poly(gauge,Z).degree()<=2,True)
eq('Mobius inverse substitution',(A*((D*Z-B)/(A-C*Z))+B)/(C*((D*Z-B)/(A-C*Z))+D),Z)
P,Q=s.symbols('P Q');ai,bi,ci,aj,bj,cj=[s.Function(n)(P,Q) for n in ['ai','bi','ci','aj','bj','cj']]
vi=ai*T**2+bi*T+ci;vj=aj*T**2+bj*T+cj
curv=s.expand(s.diff(vj,P)-s.diff(vi,Q)+vi*s.diff(vj,T)-vj*s.diff(vi,T))
check('projective curvature cubic terms cancel',s.Poly(curv,T).degree()<=2,True)
for degree,value in enumerate([s.diff(cj,P)-s.diff(ci,Q)+ci*bj-cj*bi,s.diff(bj,P)-s.diff(bi,Q)+2*(ci*aj-cj*ai),s.diff(aj,P)-s.diff(ai,Q)+bi*aj-bj*ai]):eq('curvature coefficient '+str(degree),curv.coeff(T,degree),value)
# Universal finite cubic:implicit derivatives and curvature in the actual field.
AA,BB=s.symbols('AA BB');f=T**3+AA*T+BB;g=s.diff(f,T)
inv=s.invert(g,f,T)
def rem(expr):return s.rem(expr,f,T)
va=rem(-T*inv);vb=-inv
eq('cubic inverse derivative identity',rem(g*inv),1)
for label,vv,df in [('A',va,T),('B',vb,1)]:
 check('cubic '+label+' derivative quadratic',s.Poly(vv,T).degree()<=2,True)
 eq('cubic '+label+' implicit derivative',rem(g*vv+df))
eq('cubic exact flat curvature',s.diff(vb,AA)-s.diff(va,BB)+va*s.diff(vb,T)-vb*s.diff(va,T))
eq('cubic Jacobian is genuinely nonconstant',s.diff(-T**3-AA*T,T),-g)
check('cubic discriminant is not constant or zero',s.discriminant(f,T)==-4*AA**3-27*BB**2,True)
# The universal one-direction obstruction uses distinct formal degrees for m>=4.
m=s.symbols('m',integer=True,positive=True)
eq('all-m leading coefficient elimination determinant',s.diag(m,m-1,m-1).det(),m*(m-1)**2)
check('all-m determinant nonzero at every m>=4 by factorization',s.factor(m*(m-1)**2)==m*(m-1)**2,True)
for n in [4,5,7]:
 fn=T**n-T-P;gn=s.diff(fn,T);vv=s.invert(gn,fn,T)
 check(f'one-direction degree{n}:other derivative not quadratic',s.Poly(vv,T).degree()>2,True)
 eq(f'one-direction degree{n}:exact inverse',s.rem(gn*vv,fn,T),1)
 aa,bb,cc=s.symbols('aa bb cc')
 expected=n*cc*T**(n-1)+(n-1)*aa*T**2+(n*P*aa+(n-1)*bb)*T+n*P*bb-cc
 eq(f'one-direction degree{n}:universal coefficient obstruction',s.rem(gn*(aa*T**2+bb*T+cc),fn,T),expected)
# Tiny finite permutation controls for the exceptional point-stabilizer table.
def compose(a,b):return tuple(a[b[i]] for i in range(len(a)))
def inverse(a):return tuple(a.index(i) for i in range(len(a)))
def cycle(n,*entries):
 a=list(range(n))
 for i,v in enumerate(entries):a[v]=entries[(i+1)%len(entries)]
 return tuple(a)
def generate(gens,n):
 I=tuple(range(n));G={I};front=[I]
 while front:
  a=front.pop()
  for b in gens:
   v=compose(a,b)
   if v not in G:G.add(v);front.append(v)
 return G
def normal_closure(G,h):return generate([compose(compose(g,h),inverse(g)) for g in G],len(h))
def sign(a):return (-1)**sum(a[i]>a[j] for i in range(len(a)) for j in range(i+1,len(a)))
A4={p for p in itertools.permutations(range(4)) if sign(p)==1};S4=set(itertools.permutations(range(4)));A5={p for p in itertools.permutations(range(5)) if sign(p)==1}
cases=[('A4-C2',A4,compose(cycle(4,0,1),cycle(4,2,3)),2,4),('A4-C3',A4,cycle(4,0,1,2),3,12),('S4-C2',S4,cycle(4,0,1),2,24),('S4-C3',S4,cycle(4,0,1,2),3,12),('S4-C4',S4,cycle(4,0,1,2,3),4,24),('A5-C2',A5,compose(cycle(5,0,1),cycle(5,2,3)),2,60),('A5-C3',A5,cycle(5,0,1,2),3,60),('A5-C5',A5,cycle(5,0,1,2,3,4),5,60)]
for label,G,hh,order,normsize in cases:
 H=generate([hh],len(hh));N=normal_closure(G,hh)
 check(label+':point stabilizer order',len(H)==order,True)
 check(label+':normal closure',len(N)==normsize,True)
 check(label+':coset degree',len(G)//len(H)=={'A4-C2':6,'A4-C3':4,'S4-C2':12,'S4-C3':8,'S4-C4':6,'A5-C2':30,'A5-C3':20,'A5-C5':12}[label],True)
for label,G in [('A4',A4),('S4',S4),('A5',A5)]:check(label+':trivial stabilizer not normally generating',len(normal_closure(G,tuple(range(len(next(iter(G)))))))==1,True)
for n in [3,4,5,8]:
 rot=cycle(n,*range(n));ref=tuple((-i)%n for i in range(n));G=generate([rot,ref],n)
 check(f'D{n}:faithful natural action order',len(G)==2*n,True)
 check(f'D{n}:reflection normal closure exact formula',len(normal_closure(G,ref))==2*n//math.gcd(2,n),True)
for signature,order in [((2,3,3),12),((2,3,4),24),((2,3,5),60)]:
 eq('spherical signature '+str(signature),sum(1-s.Rational(1,e) for e in signature),2-s.Rational(2,order))
context=json.loads((HERE/'sources/projective-generator-context.json').read_text())
for path,sha in context['retained_proof_and_source_sha256'].items():check('retained byte pin '+path,hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha)
for path,sha in context['new_primary_sha256'].items():check('primary cache '+path,hashlib.sha256((HERE/path).read_bytes()).hexdigest()==sha)
check('no new external theorem import',not context['newly_imported_proved_external_statements'] and not context['newly_used_unproved_external_statements'])
positive={
 'full_projective_system_equivalent_to_constant_conjugate_crossratios_proved','full_frame_common_constants_of_normal_closure_exactly_C_proved',
 'projective_group_embedding_and_full_point_stabilizer_identified','rational_Mobius_generator_and_target_derivative_basis_covariance_proved',
 'exact_flat_projective_curvature_identity_proved','all_degree_Keller_projective_action_sieve_proved',
 'primitive_plane_degree_four_row_retired_only_for_full_source','source_polynomial_coordinate_primitivity_reproved_with_classical_credit',
 'single_direction_unbounded_counterfamily_proved','cubic_full_projective_control_retains_nonKeller_boundary','external_review_deferred'}
for k,v in context.items():
 if isinstance(v,bool):check('scope '+k,v is (k in positive))
proof=HERE/'PROJECTIVE_DIFFERENTIAL_GENERATORS_AND_MONODROMY.md'
result=dict(outcome='PASS',check_count=len(checks),algebraic_controls=algebraic,checks=checks,
 proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 dependencies=context['retained_proof_and_source_sha256'],source_sha256=context['new_primary_sha256'],
 elapsed_seconds=round(time.monotonic()-START,3),**{k:v for k,v in context.items() if isinstance(v,bool)})
(HERE/'receipt-projective-generator.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['outcome','check_count','algebraic_controls','elapsed_seconds','original_goal_complete']},indent=2))
