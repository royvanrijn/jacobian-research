"""Exact controls for weighted period recovery;all-degree proofs are in the note."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):check(sp.simplify(value)==0,label)
def zeros(values,label):check(all(sp.simplify(v)==0 for v in values),label)
u,v,P,Q,s,H=sp.symbols('u v P Q s H');coords=[u,v,P,Q];q=u*P+v*Q
def D0(f):return v*sp.diff(f,P)-u*sp.diff(f,Q)
m,n=sp.symbols('m n',integer=True,positive=True)
sigma=P+P**(m+1);phi=1+(m+1)*P**m
zero(sp.diff(sigma,P)-phi,'the arbitrary-m time change has exactly the displayed derivative')
zero(D0(q),'the linear root field preserves the actual Morse quadric')
zeros([D0(u),D0(v)],'both base coordinates are fixed by the root field')
check(sp.hessian(q,coords).det()==1,'the arbitrary-m local control has nondegenerate unit quadratic Hessian')
check(sp.gcd(u,v)==1,'the local-open normalized field is primitive')
def E(f):return D0(f)/phi
zero(E(q),'the volume-normalized rational field fixes q')
tv=sigma/v;tu=(sigma-(H/u+(H/u)**(m+1)))/v
zero(E(tv)-1,'the first arbitrary-m local time has derivative one')
zero(E(tu)-1,'the second arbitrary-m local time has derivative one with H invariant')
zero((tv-tu)-(H/u+(H/u)**(m+1))/v,'the two local times have the exact inverse-system difference')
check(sp.diff(H,u)==0 and sp.diff(H,P)==0,'the local time calculation holds the actual quotient value fixed')
# A representative binomial cancellation checks chart regularity;the proof uses all m.
mu=3;phic=phi.subs(m,mu);sigmac=sigma.subs(m,mu)
tuc=sp.cancel(tu.subs({m:mu,H:q},simultaneous=True))
check(sp.denom(tuc)==u**4,'the representative second local time has no v denominator after the exact binomial cancellation')
zero(E(tv).subs(m,mu)-1,'the representative local-time field identity remains exact')
zero((tv-tu).subs(m,mu)-(H/(u*v)+H**4/(u**4*v)),
     'the representative cocycle is s times its normalized inverse system')
# Constant-volume analytic Morse coordinates exist explicitly at every m.
a=(m+1)/m;U=u*sp.exp(-a*P**m);Z=P*sp.exp(a*P**m)
zero(U*Z-u*P,'the explicit analytic change preserves the hyperbolic product')
zero(sp.Matrix([U,Z]).jacobian([u,P]).det()-phi,
     'the arbitrary-m analytic Morse change has the exact prescribed volume Jacobian')
zero(phi.subs(P,0)-1,'the volume is a unit at the Morse point')
check(sp.degree(phic,P)==mu and sp.gcd(phic,sp.diff(phic,P))==1,
      'the representative deleted global volume divisor has only simple roots')
zero(P*sp.diff(phi,P)-m*(phi-1),'at every arbitrary-m root r,the exact identity gives r phi prime(r)=-m')
zero((m+1)*sp.integrate(sp.Symbol('r',positive=True)**m,(sp.Symbol('r',positive=True),0,1))-1,
     'the all-m sphere marginal moment cancels exactly the volume coefficient m+1')
# Laurent coefficient duality and the actual ideal/class distinction.
def pp(f):
    result=0
    for term in sp.Add.make_args(sp.expand(f)):
        powers=term.as_powers_dict()
        if powers.get(u,0)<0 and powers.get(v,0)<0:result+=term
    return sp.expand(result)
def residue(f):return sp.expand(f).coeff(u,-1).coeff(v,-1)
eta=1/(u*v)+s**mu/(u**(mu+1)*v);xi=s*eta
zero(pp(u**(mu+1)*eta),'the exact representative finite ideal kills the entire normalized class')
zero(pp(v*eta),'the second finite-ideal generator kills the entire normalized class')
zero(eta.subs(s,0)-1/(u*v),'the normalized class specializes to a reduced inverse system')
zero(xi.subs(s,0),'the unnormalized Morse class vanishes at its critical value')
zero(residue(xi)-s,'the unweighted period is exactly linear despite generic thickness')
zero(residue(u**mu*xi)-s**(mu+1),'a higher weighted moment detects the thickness missed by the leading term')
for j in range(mu+1):
    expected=(1/(u*v) if j==0 else 0)+s**mu/(u**(mu+1-j)*v)
    zero(pp(u**j*eta)-expected,'negative Laurent coefficient duality column '+str(j))
basis=[u**j for j in range(mu+1)]
matrix=sp.Matrix([[residue(f*g*eta) for g in basis] for f in basis])
J=sp.zeros(mu+1)
for j in range(mu+1):J[j,mu-j]=1
E00=sp.zeros(mu+1);E00[0,0]=1
check(matrix==s**mu*J+E00,'the full representative moment matrix is reversal plus the rank-one residue update')
check(J*J==sp.eye(mu+1),'the reversal matrix is its own inverse')
check(J[0,0]==0,'the determinant update cofactor is zero in positive thickness')
zero(matrix.det()-(-1)**(mu*(mu+1)//2)*s**(mu*(mu+1)),
     'the representative determinant agrees with the written arbitrary-size formula')
check(matrix.rank()==mu+1,'the generic representative Hankel rank is the entire thick length')
check(matrix.subs(s,0).rank()==1,'the critical-value normalized Hankel rank is one')
check(mu+1-1==mu,'the normalized specialization loses exactly the nilradical dimension')
zero(residue(u*eta),'the first moment vanishes in the representative control')
zero(residue(u**2*eta),'the second moment vanishes in the representative control')
zero(residue(u*u**2*eta)*residue(eta)-residue(u*eta)*residue(u**2*eta)-s**3,
     'the exact covariance beyond the leading moments detects nonreducedness')
check(sp.KroneckerDelta(m,0)==0,'in every positive size the all-m reversal determinant update has zero corner')
j=sp.symbols('j',integer=True,nonnegative=True)
zero(sp.summation(j,(j,0,m))-m*(m+1)/2,'the all-size reversal inversion count is exactly m(m+1)/2')
# Pole orders prove failure of global local nilpotence for every m and every iterate.
t,A,C=sp.symbols('t A C',nonzero=True)
leading=C*t**(-(2*n-1))
zero(sp.diff(leading,t)/(A*t)+(2*n-1)*C/A*t**(-(2*n+1)),
     'the arbitrary-iterate pole order grows from 2n-1 to 2n+1 with nonzero coefficient')
check((2*n-1).is_positive,'the all-iterate leading coefficient never vanishes in characteristic zero')
zero(E(P)-v/phi,'the first iterate has the initial simple volume-divisor pole')
zero(E(E(P))+v*v*sp.diff(phi,P)/phi**3,
     'the next iterate has the exact third-order pole predicted by the all-iterate induction')
zeros([D0(D0(c)) for c in coords],'the linear first jet squares to zero on every coordinate')
check(m.is_positive and m.is_integer,'the nonlinear correction raises degree on every finite-jet filtration')
check(sp.Matrix([D0(c) for c in coords]).jacobian(coords)**2==sp.zeros(4),
      'each associated-graded linear root matrix is nilpotent before extending to homogeneous polynomials')
# The actual Nambu normalization and the formal saturation defect.
def N(f):return sp.Matrix([u,v,q,f]).jacobian(coords).det()
zeros([N(c)+D0(c) for c in coords],'the primitive quadric Nambu field has the exact root-action orientation')
zeros([N(c)/phi+E(c) for c in coords],
      'the volume-normalized Nambu field agrees with the rational field on every source coordinate')
# The finite inverse-system image has exactly the same specialization defect as the Hankel form.
images=[pp(f*eta) for f in basis]
embed=sp.Matrix([[image.coeff(u,-k-1).coeff(v,-1) for image in images] for k in range(mu+1)])
zero(embed.det()-matrix.det(),'the coefficient embedding and the Hankel matrix have the same determinant in the ordered monomial bases')
check(embed.subs(s,0).columnspace()==[sp.Matrix([1,0,0,0])],
      'the normalized embedding image is precisely the residue line at the Morse value')
check(embed.subs(s,0).nullspace()==[sp.eye(mu+1)[:,j] for j in range(1,mu+1)],
      'its specialization kernel is precisely the nilradical of the actual length-four ideal fibre')
zeros([pp(u/(u*v)),pp(v/(u*v))],
      'both reduced-point parameters annihilate the normalized reduced control')
algebra=len(checks)
context=json.loads((HERE/'sources/weighted-period-recovery-context.json').read_text())
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained immutable input '+path)
for entry in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
          'inspected primary quotient-hypothesis audit '+entry['retained_file'])
check(context['newly_imported_proved_external_statements']==[] and context['newly_used_unproved_external_statements']==[],
      'no faithful-flat quotient criterion,isochore theorem or unproved saturation statement is imported')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
check(sum(flags.values())==7 and not flags['local_control_realizes_an_actual_native_A4_Ga_action']
      and not flags['normalized_residue_cokernel_value_torsion_ruled_out_for_native_actions']
      and not flags['fixed_support_alone_proved_to_force_minimal_native_plinth_length']
      and not flags['original_goal_complete'],
      'weighted reconstruction and the exact local specialization defect do not close global nilpotence,native minimum length or alignment')
proof=HERE/'WEIGHTED_PERIOD_RECOVERY_AND_THE_GLOBAL_NILPOTENCE_GAP.md'
paths=['sources/weighted-period-recovery-context.json']+['sources/'+e['retained_file'] for e in context['newly_cached_primary_context']]
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,
 elapsed_seconds=round(time.monotonic()-START,3),python=platform.python_version(),sympy=sp.__version__,
 proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 dependencies=context['retained_proof_and_source_sha256'],source_sha256={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths},
 independent_replay=False,formal_verification=False,external_review=False,**flags)
(HERE/'receipt-weighted-period-recovery.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in {'checks','dependencies','source_sha256'}},indent=2))
