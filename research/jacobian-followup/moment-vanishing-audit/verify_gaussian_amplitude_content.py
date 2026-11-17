"""One-use finite Taylor/source/Wick/content controls; outputs must be absent."""
from pathlib import Path
from functools import lru_cache
import hashlib,json,platform
import sympy as s
import gaussian_amplitude_content as G

p=Path(__file__).resolve().parent;family='gaussian_amplitude_content'
outputs=[p/(family+tail+'.json') for tail in ['_certificate','_runtime','_verification']]
assert all(not f.exists() for f in outputs),'Sealed outputs exist;do not rerun.'
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
dump=lambda x:json.dumps(x,indent=2,ensure_ascii=False)+'\n'
x,y,z,w,t,u,a,b,c,e=s.symbols('x y z w t u a b c e')

def wick(poly,X,C):
 C=s.Matrix(C)
 @lru_cache(None)
 def moment(alpha):
  if not any(alpha):return s.Integer(1)
  if sum(alpha)%2:return s.Integer(0)
  i=next(j for j,k in enumerate(alpha) if k)
  v=list(alpha);v[i]-=1;answer=0
  for j,k in enumerate(v):
   if k:
    vv=v.copy();vv[j]-=1;answer+=k*C[i,j]*moment(tuple(vv))
  return answer
 return s.expand(sum(cc*moment(aa) for aa,cc in s.Poly(s.expand(poly),*X).terms()))

def period(amp,phases,X,times,C,N):
 answer=0
 for counts in G.multiindices(len(times),N):
  poly=amp*s.prod(P**k for P,k in zip(phases,counts))
  answer+=wick(poly,X,C)*s.prod(v**k/s.factorial(k) for v,k in zip(times,counts))
 pp=s.Poly(s.expand(answer),*times)
 return s.expand(sum(cc*s.prod(v**k for v,k in zip(times,aa)) for aa,cc in pp.terms() if sum(aa)<=N))

fixtures=[
 (2+3*x+x**2,1+x**2,(x**3+x,),(x,),(t,),[[1]],[(0,),(1,),(2,)]),
 (1+z+w+z*w,w**2+z,(z+w**2,),(z,w),(t,),[[0,1],[1,0]],[(0,0),(1,0),(0,1),(1,1)]),
 (a+b*x+c*x**2,1+x,(x**2,),(x,),(t,),[[2]],[(0,),(1,),(2,)]),
 (1+x+x**2,x**3,(x**2,x**3),(x,),(t,u),[[1]],[(0,),(2,)])]
witnesses=[]
for weight,amp,phases,X,times,C,indices in fixtures:
 for alpha in indices:
  cert=G.coefficient_witness(weight,alpha,amp,phases,X,times,C)
  recovery=s.sympify(cert['recovery_amplitude'])
  lhs=period(s.sympify(cert['coefficient'])*amp,phases,X,times,C,4)
  rhs=period(weight*recovery,phases,X,times,C,4)
  assert s.expand(lhs-rhs)==0,cert
  cert.update(Wick_jet_order=4,original_Wick=str(lhs),weighted_recovery_Wick=str(rhs))
  witnesses.append(cert)

tops=[]
for weight,alpha,amp,phases,X,times,C in [(2+x**2,(2,),1+x,(x**3,),(x,),(t,),[[1]]),((1+z)*(1+w),(1,1),w**2,(z,),(z,w),(t,),[[0,1],[1,0]])]:
 cert=G.coefficient_witness(weight,alpha,amp,phases,X,times,C,top_only=True)
 lhs=period(s.sympify(cert['coefficient'])*amp,phases,X,times,C,4)
 rhs=period(weight*s.sympify(cert['recovery_amplitude']),phases,X,times,C,4)
 assert s.expand(lhs-rhs)==0
 tops.append(cert)

rho=a+(1-a)*x;amp=1+x**2;P=x**2
recovery=s.expand(amp+G.adjoint((1-x)*amp,0,(P,),(x,),(t,),[[1]]))
assert s.expand(rho+(1-x)*s.diff(rho,x)-1)==0
assert s.expand(period(amp,(P,),(x,),(t,),[[1]],5)-period(rho*recovery,(P,),(x,),(t,),[[1]],5))==0
primitive=dict(weight=str(rho),coefficient_content_witness='a+(1-a)=1',amplitude=str(amp),recovery=str(recovery),coefficient_units_assumed=False)

def basic_witt(A,B):
 return s.factorial(B)*t**(B-A)/s.factorial(B-A) if B>=A else s.Integer(0)

positive=[]
rho=(1+z)*(1+w)
for aa,bb in [(0,0),(0,1),(1,0),(0,2),(1,1),(2,0),(1,2),(0,3)]:
 exact=s.expand((basic_witt(aa,bb)+basic_witt(aa+1,bb)+basic_witt(aa,bb+1)+basic_witt(aa+1,bb+1))/2)
 actual=period(rho*z**aa*w**bb,(z,),(z,w),(t,),[[0,1],[1,0]],aa+bb+3)/2
 assert s.expand(exact-actual)==0
 assert s.degree(exact,t)<=aa+bb+1
 positive.append(dict(exponents=[aa,bb],exact_period=str(exact),degree_bound=aa+bb+1))
seed=positive[0]['exact_period'];assert s.expand(s.sympify(seed)-1-t/2)==0
assert s.cancel((1-s.sympify(seed))/t)==-s.Rational(1,2)
for degree in range(5):assert s.expand((basic_witt(0,degree)+basic_witt(1,degree)+basic_witt(0,degree+1)+basic_witt(1,degree+1))/2).coeff(t,degree+1)*s.factorial(degree+1)==s.factorial(degree+1)/2
assert period(z,(z,),(z,w),(t,),[[0,1],[1,0]],4)==0
assert period(z*w,(z,),(z,w),(t,),[[0,1],[1,0]],4)==1

relations=s.groebner([a*e,e**2],a,e,order='lex')
mask=[]
for power in range(1,6):
 remainder=relations.reduce(s.expand(e*a**power))[1];assert remainder==0
 mask.append(dict(power=power,remainder=str(remainder)))

transfers=[G.cutoff_transfer('1','0',2,1,'original_to_weighted'),G.cutoff_transfer('1','1',2,1,'weighted_to_original'),G.cutoff_transfer('3/2','2',3,4,'weighted_to_original')]
assert [v['E'] for v in transfers]==['2','5','37/2']

rejects=[]
for desc,fn in [
 ('zero weight',lambda:G.coefficient_witness(0,(0,),1,(x,),(x,),(t,),[[1]])),
 ('invalid covariance',lambda:G.coefficient_witness(x,(1,),1,(x,),(x,),(t,),[[0]])),
 ('invalid top index',lambda:G.coefficient_witness(1+x**2,(1,),1,(x,),(x,),(t,),[[1]],top_only=True)),
 ('noninteger multiindex',lambda:G.coefficient_witness(x,(s.Rational(1,2),),1,(x,),(x,),(t,),[[1]])),
 ('negative cutoff',lambda:G.cutoff_transfer('-1','0',1,1,'original_to_weighted'))]:
 try:fn()
 except ValueError:rejects.append(desc)
 else:raise AssertionError(desc)

assert wick(s.diff(x,x),(x,),[[1]])==1  # source HTML m-1 boundary at n0,m1
packet=dict(date='2026-11-17',universal_claims_use_written_proofs=True,actual_image_kernel_computed=False,unweighted_nonpure_safe_counterexample_claimed=False,field_safety_implies_purity_resolved=False,source_identities= witnesses,top_coefficient_identities=tops,primitive_nonunit_content=primitive,positive_tilt=dict(phase='z',weight='(1+z)*(1+w)',mass='2',normalized_seed=seed,exact_controls=positive,sharp_C='1',sharp_E='1',finite_seed_shift_eigenvalue='-1/2',unweighted_seed='1',unweighted_E='0',weighted_measure_is_unweighted_Gaussian=False),zero_scalar_period_nonzero_full_ideal_image=True,ring_mask=dict(ring='QQ[a,e]/(a*e,e**2)',phase='a*x**2',weight='e',normalization_ideal='(e)',content='(e)',weighted_normalization='R',single_structural_annihilation_relation='a*e=0',finite_controls=mask),cutoff_transfers=transfers,rejections=rejects,source_HTML_monomial_boundary=dict(n=0,m=1,correct_factor=1,displayed_m_minus_one=0,PDF_display_newly_checked=False))
imports=['STEIN_ACTION_AND_RANK_ONE_LATTICES.md','GAUSSIAN_G_FUNCTIONS_AND_GLOBAL_REGULAR_CORE.md','gaussian_gfunction_regular_core_verification.json','BERNSTEIN_CLASSICAL_CREDIT_AND_FORMAL_EXTENSIONS.md']
runtime=dict(python=platform.python_version(),sympy=s.__version__,sealed_imports_sha256={nm:sha(p/nm) for nm in imports})
outputs[0].write_text(dump(packet));outputs[1].write_text(dump(runtime))
receipt=dict(completed=True,kind='finite original Wick and source-content controls;universal identities use written proofs',canonical_proof_sha256=sha(p/'GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md'),producer_sha256=sha(p/'gaussian_amplitude_content.py'),checker_sha256=sha(Path(__file__)),certificate_sha256=sha(outputs[0]),runtime_sha256=sha(outputs[1]),exact_controls=dict(coefficient_recoveries=len(witnesses),top_recoveries=len(tops),primitive_nonunit_content=1,positive_tilt_markers=len(positive),ring_masks=1,cutoff_transfers=len(transfers),zero_scalar_seed_control=1,source_HTML_boundary=1,rejections=len(rejects)),independent_replay=False,formal_verification=False,external_review=False)
outputs[2].write_text(dump(receipt));print(dump(receipt))
