"""MX0 new supporting controls, not an independent whole-proof verification."""
import hashlib,json,time
from pathlib import Path
import sympy as S
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
receipt=HERE/'receipt-logarithmic-cyclic-source-images.json'
assert not receipt.exists(),'Sealed receipt exists: never replay.'
assert not any(e['id']=='HC4MYMX0' for e in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries'])
started=time.perf_counter();checks=[];algebraic=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(name,value,alg=True):
 assert bool(value),name
 checks.append(name)
 if alg:algebraic.append(name)
t,u,x,z,g1,g2,g3=S.symbols('t u x z g1 g2 g3')
check('unit relation exact multiplicative cocycle',S.cancel((g3/g2)*(g2/g1)-g3/g1)==0)
a,b,c=S.symbols('a b c')
check('dlog relation exact additive cocycle',S.expand((c-b)-(c-a)+(b-a))==0)
for order in [1,4]:
 check(f'nonzero valuation dlog has simple pole order {order}',S.cancel(S.diff(t**order,t)/t**order-order/t)==0)
 check(f'exact rational primitive has higher pole order {order}',S.limit(t**(order+1)*S.diff(t**(-order),t),t,0)==-order)

# Two fixed transcriptions. The proof is uniform in all g,k,B.
for degree,k,Bu in [(2,2,u-1),(3,2,u**2+u+1)]:
 tag=f'g{degree}k{k}'
 sigmau=-Bu+degree*(1-k)*u*S.diff(Bu,u)
 Bt=Bu.subs(u,t**degree);sigma=sigmau.subs(u,t**degree)
 at=t**2*Bt**k;G=t**-1*Bt**(1-k)
 check(tag+' base squarefree and nonzero at zero',S.gcd(Bu,S.diff(Bu,u))==1 and Bu.subs(u,0)!=0)
 check(tag+' moving label differential identity',S.cancel(S.diff(G,t)-sigma/at)==0)
 check(tag+' moving labels separate all bad fibres',S.gcd(sigmau,u*Bu)==1)
 check(tag+' a has cyclic weight two',all(int(mon[0])%degree==2%degree for mon in S.Poly(at,t).monoms()))
 check(tag+' sigma is cyclic invariant',all(int(mon[0])%degree==0 for mon in S.Poly(sigma,t).monoms()))
 vv=sigma+at*x
 zz=S.cancel((vv**degree-sigma**degree)/at)
 check(tag+' chart numerator clears polynomially',S.denom(zz)==1)
 check(tag+' complete binomial chart formula',S.expand(zz-sum(S.binomial(degree,j)*sigma**(degree-j)*at**(j-1)*x**j for j in range(1,degree+1)))==0)
 check(tag+' chart pulls back volume to dxdt',S.cancel(S.Matrix([vv,t]).jacobian([x,t]).det()/at-1)==0)
 modulus=u*Bu**k
 C0=S.rem(-Bu*S.invert(sigmau,modulus,u),modulus,u)
 D0=S.cancel((Bu+C0*sigmau)/modulus)
 check(tag+' quotient Bezout certificate',S.degree(S.denom(D0),u)==0)
 hp=degree*t**(degree-1);beta=S.diff(G,t)/hp
 sl=x/hp+beta
 F0=t*C0.subs(u,t**degree)*x+t**(degree-1)*D0.subs(u,t**degree)
 check(tag+' rational slice has prescribed primitive',S.cancel(beta*hp-S.diff(G,t))==0)
 check(tag+' quotient invariant slice formula',S.cancel(sl-vv/(degree*t**(degree+1)*Bt**k))==0)
 check(tag+' EXPLICIT polynomial source representative',S.cancel(F0-G-degree*t**degree*C0.subs(u,t**degree)*sl)==0)
 check(tag+' primitive trace zero character',(-1)%degree!=0)
 r=S.degree(Bu,u)
 check(tag+' quotient Euler and free Picard ranks agree',1+r*(degree-1)==(1+(degree-1)*(1+degree*r))/degree)

# In the exact g=2,B=u-1 model, t-1 represents a free component class.
# It cannot have closed-source logarithmic splitting, since its mth power
# is not invariant under the ACTUAL cyclic field action on t.
check('exact cyclic model has source polynomial free logarithmic defect',S.expand((t-1)**2-(-t-1)**2)==-4*t)

cp=HERE/'sources/logarithmic-cyclic-source-images-context.json';c=json.loads(cp.read_text())
for path,digest in c['retained_proof_and_source_sha256'].items():check('retained byte pin '+path,sha(ROOT/path)==digest,False)
for path,digest in c['retained_primary_sha256'].items():check('retained primary pin '+path,sha(ROOT/path)==digest,False)
for key in ['logarithmic_quotient_equals_Picard','closed_source_part_equals_torsion','cyclic_exact_surjective_plane_source_all_degrees','arbitrary_free_Picard_rank','explicit_polynomial_representative','PF142_minimal_mixed_retirement_retained','full_trace_goal_open']:
 check('scope '+key,c['scope_flags'][key],False)
for key in ['logarithmic_Picard_identified_with_trace_defect','closed_source_splitting_automatic','plane_target_constructed','source_shear_asserted_polynomial','unrestricted_HC4_or_JC2_resolved']:
 check('boundary '+key,not c['scope_flags'][key],False)
check('no new cache download import or registered replay',c['new_primary_caches']==c['new_primary_downloads']==c['new_external_theorem_imports']==0 and c['registered_mathematical_checkers_replayed']==[],False)
out=dict(outcome='PASS',claim='HC4MYMX0',check_count=len(checks),algebraic_controls=len(algebraic),checks=checks,elapsed_seconds=time.perf_counter()-started,
 proof_sha256=sha(HERE/'LOGARITHMIC_DESCENT_AND_CYCLIC_EXACT_SOURCE_IMAGES.md'),checker_sha256=sha(Path(__file__)),context_sha256=sha(cp),registered_mathematical_checkers_replayed=[],independent_replay=False,formal_verification=False,external_review=False)
receipt.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
