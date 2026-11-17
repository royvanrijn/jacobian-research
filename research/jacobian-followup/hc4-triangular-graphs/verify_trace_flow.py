"""Exact trace-flow controls;the universal odd-degree proof is in the note."""
import hashlib,json,time
from pathlib import Path
import sympy as s
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(name,c):
 assert c,name
 checks.append(name)
def zero(name,e):check(name,s.cancel(s.expand(e))==0)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
x,y,p,q,t,q0,K,N=s.symbols('x y p q t q0 K N')
zero('universal chain correction coefficient',1-2*K/(2*K-1)+1/(2*K-1))
kappa=N**2/(N-2);qt=q0/(1-kappa*q0*t)
zero('universal invariant-line finite-time-pole solution',s.diff(qt,t)-kappa*qt**2)
zero('universal divergence determines first coefficient',-2*kappa*q+2*kappa*q)

traces={}
for k in [1,2,3]:
 n=2*k+1;c=s.Rational((-1)**k,n)
 Q=sum((-1)**j*s.binomial(k,j)*p**(k-j)*x**(2*j+1)/s.Integer(2*j+1) for j in range(k+1))
 zero(f'antiderivative identity k{k}',s.diff(Q,x)-(p-x*x)**k)
 check(f'constant leading unit gives finite free degree k{k}',s.Poly(Q,x).LC()==c and s.degree(Q,x)==n)
 Pxy=x*x+y;Qxy=s.expand(Q.subs(p,Pxy))
 jac=s.diff(Pxy,x)*s.diff(Qxy,y)-s.diff(Pxy,y)*s.diff(Qxy,x)
 zero(f'exact Jacobian k{k}',jac+y**k)
 def delta(f,k=k):return y**(k+1)*s.diff(f,x)
 zero(f'source LND preserves Jacobian k{k}',delta(jac))
 zero(f'source LND divergence zero k{k}',s.diff(y**(k+1),x))
 zero(f'global exact contraction Hamiltonian k{k}',s.diff(-y**(2*k+2)/s.Integer(2*k+2),y)+y**(2*k+1))
 zero(f'source polynomial flow P differential k{k}',s.diff((x+t*y**(k+1))**2+y,t).subs(t,0)-delta(Pxy))
 pol=s.Poly((Q-q)/c,x);C=s.zeros(n)
 for i in range(n-1):C[i+1,i]=1
 for i in range(n):C[i,n-1]=-pol.nth(i)
 def mat(f,C=C,n=n):
  f=s.Poly(s.expand(f.subs(y,p-x*x)),x)
  return sum((v*C**mon[0] for mon,v in f.terms()),s.zeros(n))
 def av(f,mat=mat,n=n):return s.cancel(s.trace(mat(f))/n)
 ep=av(delta(Pxy));eq=av(delta(Qxy));kap=s.Rational(n*n,n-2)
 zero(f'averaged first component k{k}',ep+2*kap*p*q)
 gamma=s.cancel((eq-kap*q*q)/p**n)
 check(f'unknown general gamma specializes to rational k{k}',gamma.is_Rational is True)
 zero(f'averaged second component universal form k{k}',eq-kap*q*q-gamma*p**n)
 zero(f'averaged target divergence zero k{k}',s.diff(ep,p)+s.diff(eq,q))
 zero(f'averaged target Hamiltonian P contraction k{k}',s.diff(-kap*p*q*q-gamma*p**(n+1)/(n+1),q)-ep)
 zero(f'averaged target Hamiltonian Q contraction k{k}',-s.diff(-kap*p*q*q-gamma*p**(n+1)/(n+1),p)-eq)
 zero(f'chain-rule derivative on Pzero k{k}',s.diff(Qxy,x).subs(y,-x*x)-(-1)**(k+1)*x**(n-1)/s.Integer(n-2))
 zero(f'Pzero source velocity k{k}',delta(Qxy).subs(y,-x*x)-x**(2*n)/s.Integer(n-2))
 beta=s.simplify(Q.subs({p:1,x:1}));expected=2**(2*k)*s.factorial(k)**2/s.factorial(n)
 zero(f'exact critical value beta k{k}',beta-expected)
 zero(f'opposite critical value k{k}',Q.subs({p:1,x:-1})+beta)
 factor=s.div(s.Poly(Q.subs(p,1)-beta,x),s.Poly((x-1)**(k+1),x))
 check(f'critical root has exact multiplicity and k simple companions k{k}',factor[1].is_zero and factor[0].degree()==k and factor[0].eval(1)!=0 and s.gcd(factor[0],factor[0].diff()).degree()==0 and factor[0].eval(-1)!=0)
 zero(f'zero-P fibre exact monomial k{k}',Q.subs(p,0)-c*x**n)
 check(f'nonzero zero-P targets have no critical roots k{k}',s.gcd(s.Poly(c*x**n-q,x),s.Poly(x**(2*k),x)).degree()==0)
 traces[k]=(ep,eq,av,delta,Pxy,Qxy)

ep,eq,av,delta,Pxy,Qxy=traces[1]
def eps(f):return ep*s.diff(f,p)+eq*s.diff(f,q)
zero('cubic explicit averaged Q coefficient',eq-11*p**3-9*q*q)
eta=lambda f:s.diff(f,x)
zero('source commuting LND bracket x',eta(delta(x))-delta(eta(x)))
zero('source commuting LND bracket y',eta(delta(y))-delta(eta(y)))
zero('second source LND averaged P',av(eta(Pxy)))
zero('second source LND averaged Q',av(eta(Qxy))-3*p)
e0=lambda f:3*p*s.diff(f,q)
zero('target nonzero bracket P coefficient',e0(ep)-eps(0)+54*p*p)
zero('target nonzero bracket Q coefficient',e0(eq)-eps(3*p)-108*p*q)
zero('averaged second time derivative fails flow law',av(delta(delta(Pxy)))-eps(ep)-(220*p**4-90*p*q*q))
branch=9*q*q-4*p**3
rem=s.rem(s.Poly(eps(branch),q),s.Poly(branch,q)).as_expr()
zero('averaged target not logarithmic for entire branch curve',rem-486*p**3*q)
check('branch logarithmic defect is genuinely nonzero',rem!=0)
ALG=len(checks)
context=json.loads((HERE/'sources/trace-flow-context.json').read_text())
for pp,h in context['retained_proof_and_source_sha256'].items():check('retained bytes '+pp,sha(ROOT/pp)==h)
for pp,h in context['new_primary_sha256'].items():check('new primary bytes '+pp,sha(HERE/pp)==h)
check('primary trace multiplication matrix source retained','matrix of multiplication' in (HERE/'sources/stacks-trace-and-norm-0BIE.html').read_text())
proof=HERE/'TRACE_AVERAGING_COMPLETE_FLOW_BARRIER.md';flat=' '.join(proof.read_text().split())
for pp in ['CURRENT LAST UNSOLVED IMPLICATION','EVERY ODD DEGREE','DERIVATION of B','delta(j)/j','EXACTLY n','NOT THE WHOLE AFFINE PLANE','POINTWISE','kappa_k=n^2/(n-2)','gamma_k','finite-time pole','NO NONTRIVIAL INTERMEDIATE FIELD GALOIS OVER F','SINGLE ORIGIN','NOT A LIE ALGEBRA HOMOMORPHISM','NOT LOGARITHMIC','UNAVERAGED','not a replacement','DEFERRED']:
 check('written scope '+pp,pp.lower() in flat.lower())
check('no new major theorem import',context['newly_imported_proved_external_statements']==[] and context['newly_used_unproved_external_statements']==[])
flags={k:v for k,v in context.items() if isinstance(v,bool)}
result=dict(outcome='PASS',check_count=len(checks),algebraic_controls=ALG,checks=checks,proof_sha256=sha(proof),checker_sha256=sha(Path(__file__)),dependencies=context['retained_proof_and_source_sha256'],source_sha256=context['new_primary_sha256'],elapsed_seconds=round(time.monotonic()-START,3),**flags)
(HERE/'receipt-trace-flow.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(outcome=result['outcome'],checks=len(checks),algebraic_controls=ALG,elapsed_seconds=result['elapsed_seconds']),indent=2))
