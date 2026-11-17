"""One-use exact finite controls; never infer universal claims from these jets.

Writes a new certificate, runtime binding and receipt only when absent.
"""
from pathlib import Path
import hashlib,json,platform,sys
import sympy as s
import gaussian_saddle_arithmetic as A

p=Path(__file__).resolve().parent
family='gaussian_saddle_arithmetic'
out=[p/(family+tail+'.json') for tail in ['_certificate','_runtime','_verification']]
assert all(not f.exists() for f in out),'Sealed output exists; do not rerun.'
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
dump=lambda o:json.dumps(o,indent=2,ensure_ascii=False)+'\n'
x,y=s.symbols('x y');t=A.t
controls=[];matrix=[];count_integral=0
fixtures=[
 ('linear',x,1+x,(x,),6),
 ('quadratic',x**2/2,1+x**4,(x,),6),
 ('critical_quadratic',x+x**2,x**2,(x,),6),
 ('critical_cubic',x+x**3,x+1,(x,),6),
 ('quartic_divergent',x**4,1+x**2,(x,),5),
 ('mixed_no_linear',x**2*y+x*y**2,x**2*y**2+1,(x,y),4),
 ('rational_coefficients',x/3+x**3/2,1+x/5,(x,),5),
 ('constant',s.Rational(2,3),1+x**2,(x,),4),
]
for name,P,Q,X,n in fixtures:
 g,psi,B,L=A.critical(P,X,n)
 f=A.period(P,Q,X,n)
 H=A.jet(s.exp(-psi)*f,n)
 saddle=A.normalized_via_saddle(P,Q,X,n)
 assert A.jet(H-saddle,n)==0,name
 D=int(s.ilcm(*[c.q for c in s.Poly(P,*X).coeffs()],1))
 DQ=int(s.ilcm(*[c.q for c in s.Poly(Q,*X).coeffs()],1))
 assert A.scaled_integral(H,4*D,DQ),name
 count_integral+=1
 for z in (L.T*B*L-s.eye(len(X))):assert A.jet(z,n)==0,name
 assert all(A.jet(g[i]-t*s.diff(P,X[i]).subs(dict(zip(X,g)),simultaneous=True),n)==0 for i in range(len(X))),name
 controls.append(dict(name=name,phase=str(P),marker=str(Q),variables=list(map(str,X)),jet_order=n,critical_branch=list(map(str,g)),critical_action=str(psi),original_period=str(f),normalized_period=str(H),saddle_period=str(saddle),phase_denominator=D,marker_denominator=DQ,scale=4*D))

# A coupled matrix root, with t-dependent coefficients that need not commute
# coefficientwise. The full power series matrix identity is what is used.
P=x+y+x**2*y+x*y**2+x**2
g,psi,B,L=A.critical(P,(x,y),3)
assert all(A.jet(z,3)==0 for z in L.T*B*L-s.eye(2))
matrix.append(dict(phase=str(P),jet_order=3,branch=list(map(str,g)),action=str(psi),B=[[str(B[i,j]) for j in range(2)] for i in range(2)],L=[[str(L[i,j]) for j in range(2)] for i in range(2)]))

occupations=[]
for alpha,counts,gamma in [
 ([(2,),(3,)],[3,2],(0,)),
 ([(2,0),(1,1),(0,3)],[5,4,2],(0,0)),
 ([(2,0),(1,1),(0,2)],[9,4,3],(2,0)),
 ([(2,1),(1,2)],[3,5],(1,1)),
 ([(4,)],[7],(2,)),
]:
 value=A.occupation(alpha,counts,gamma)
 assert (value*2**sum(counts)).is_Integer
 occupations.append(dict(vertices=alpha,counts=counts,marker=gamma,value=str(value),two_power=sum(counts)))

no_linear=[]
for P,Q,X in [(x**2+x**3,1+x,(x,)),(x*y+x**2*y**2,1+x**2,(x,y)),(x**4,1+x**4,(x,))]:
 v=A.period(P,Q,X,6)
 assert A.scaled_integral(v,2)
 no_linear.append(dict(phase=str(P),marker=str(Q),variables=list(map(str,X)),period=str(v),jet_order=6,scale=2))

gvar=s.Symbol('g')
tt=gvar/(1+3*gvar**2)
action=gvar**2*(1-gvar**2)/(2*(1+3*gvar**2))
conjugate=s.cancel(action.subs(gvar,1/(3*gvar)))
assert s.cancel(tt.subs(gvar,1/(3*gvar))-tt)==0
assert s.cancel(action-conjugate)!=0
assert s.cancel(tt*(gvar+gvar**3)-gvar**2/2-action)==0
obstruction=dict(phase='x**3 + x',critical_equation='3*t*g**2-g+t',parameter_time=str(tt),action=str(action),other_root='1/(3*g)',conjugate_action=str(conjugate),difference=str(s.factor(action-conjugate)),conclusion='selected critical action nonrational;unit logarithmic derivative not rational;actual rank>=2 by written theorem,not finite jets')

models=[]
for R,q,N in [
 (0,1/(1-t),2),
 (t**2,1-t,1),
 (t**3/(1-t),((1-t)**2/(1+t))**2,6),
 (0,(1-t**2)**2,3),
]:
 R=s.sympify(R);q=s.sympify(q)
 r=s.cancel(s.diff(R,t)+s.diff(q,t)/(N*q))
 cert=A.kummer_certificate(R,r)
 assert s.cancel(s.diff(s.sympify(cert['q']),t)/s.sympify(cert['q'])-cert['N']*(r-s.diff(R,t)))==0
 models.append(cert)

curves=[]
for prime in [3,5,7]:
 for a in [s.Rational(1,2)/(1-t),2/(t-1)+1/(t+1),1/(1+t**2)]:
  z=s.Integer(1)
  for _ in range(prime):z=A.mod_rational(s.diff(z,t)+a*z,prime)
  c=A.curvature(a,prime)
  assert A.mod_rational(c-z,prime)==0
  curves.append(dict(log_derivative=str(a),prime=prime,closed_formula=str(c),iterated_derivative=str(z)))
assert A.curvature(1/(1+t**2),3)!=0
assert A.curvature(1/(1+t**2),5)==0
assert A.curvature(1/(1+t**2),7)!=0
R=t**2/(1-t);rr=s.diff(R,t)+s.Rational(1,2)/(1-t)
assert A.mod_rational(A.curvature(rr,5)-s.diff(R,t)**5,5)==0

rejects=[]
for desc,fn in [
 ('linear phase vertices in degree>=2 denominator lemma',lambda:A.occupation([(1,)],[3],(1,))),
 ('irrational normalized residues exp(arctan t)',lambda:A.kummer_certificate(0,1/(1+t**2))),
 ('remaining polynomial logarithmic derivative',lambda:A.kummer_certificate(0,t)),
 ('remaining higher-order pole',lambda:A.kummer_certificate(0,1/(1-t)**2)),
 ('action not normalized at origin',lambda:A.kummer_certificate(1,0)),
]:
 try:fn()
 except ValueError:rejects.append(desc)
 else:raise AssertionError(desc)
# The raw linear exponential is deliberately outside the denominator lemma.
assert not A.scaled_integral(A.period(x,1,(x,),10),4)
assert '1/(1 - t)' not in str(A.period(x**4,1,(x,),5))

packet=dict(date='2026-11-17',universal_claims_use_written_proofs=True,actual_rank_inferred_from_finite_jets=False,Gaussian_realization_inferred_from_rational_models=False,higher_rank_pcurvature_conjecture_assumed=False,field_safety_implies_purity_claimed=False,regular_transcendental_rank_one_excluded_by_written_proof=True,critical_normalization_controls=controls,matrix_controls=matrix,occupation_controls=occupations,no_linear_controls=no_linear,critical_conjugate_obstruction=obstruction,Kummer_models=models,pcurvature_controls=curves,rejections=rejects,raw_linear_normalization_boundary=True)
runtime=dict(python=platform.python_version(),sympy=s.__version__,implementation='finite exact rational polynomial/matrix arithmetic',sealed_imports_sha256={},new_primary_sources_sha256={f.name:sha(f) for f in (p/'sources').iterdir() if f.name.startswith(('bostan_caruso_roques_2304_05061v2','stacks_0559','stacks_0578','stacks_054K','furnsinn_pannier_2510_00892v3'))})
out[0].write_text(dump(packet));out[1].write_text(dump(runtime))
receipt=dict(completed=True,kind='finite exact controls;universal field/arithmetic descent is written mathematics',canonical_proof_sha256=sha(p/'CRITICAL_NORMALIZATION_AND_RANK_ONE_ARITHMETIC.md'),producer_sha256=sha(p/'gaussian_saddle_arithmetic.py'),checker_sha256=sha(Path(__file__)),certificate_sha256=sha(out[0]),runtime_sha256=sha(out[1]),exact_controls=dict(critical_normalizations=len(controls),matrix_controls=len(matrix),occupation_controls=len(occupations),no_linear_denominator_controls=len(no_linear),critical_conjugate_obstructions=1,Kummer_models=len(models),pcurvature_controls=len(curves),rejections=len(rejects)),independent_replay=False,formal_verification=False,external_review=False)
out[2].write_text(dump(receipt));print(dump(receipt))
