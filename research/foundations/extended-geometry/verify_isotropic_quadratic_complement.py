#!/usr/bin/env python3
"""Exact source authentication and BOUNDED rational support for the written theorem."""
from pathlib import Path,PurePosixPath
from fractions import Fraction as Q
from math import factorial,prod,comb
from functools import lru_cache
import argparse,copy,hashlib,json,re,sys
HERE=Path(__file__).resolve().parent
PN='ISOTROPIC_QUADRATIC_COMPLEMENT_COUNTEREXAMPLES.md'
CN='verify_isotropic_quadratic_complement.py'
IN='isotropic_quadratic_complement_inputs.json'
LN='isotropic_quadratic_complement_source_ledger.json'
RN='isotropic_quadratic_complement_verification.json'
SCOPE={'schema': 1, 'reviewed_proof_sha256': 'fd84d668b215a80b5f4ac949e9d307e5bbcd6055e150454528f513518186bb19', 'base_field': 'EVERY characteristic-zero K; every extension', 'quadratic_scope': 'EVERY nondegenerate ISOTROPIC q, dimension N>=3; choose q=2XY+sum a_i Z_i^2, all a_i nonzero', 'M': 'N-2>=1', 'functional': 'unique normalized same-field SO(q) x G_m invariant Reynolds trace', 'seed': 'p=(1+X)^3*Y/q-(X+1)*(X+2)/2', 'diagonal_moments': 'mu_k=k!/[2^k*(N/2)_k], all k>=0', 'mixed_moments': 'ALL k>=1,j>=0: mu_k*binom(k-M,k-j) for j<=k; zero j>k; generalized upper binomial', 'pure_seed_moments': 'ALL k>=1: (-1)^k*binom(M-1,k)*mu_k', 'pure_power_iff': 'EVERY integer a>=1: p^a is pure IFF a>=M', 'minimal_monomial_marker': 'j<M vanishes for all k>=M; j>=M survives all k>=max(M,j); least detecting exponent M', 'power_phase': 'F0=p^M', 'nonpower_phase': 'Fstar=p^M*(1+p)', 'fixed_marker': 'H=X^M, independent of m', 'marked_power': 'ALL m>=1: Lambda(H*F0^m)=mu_(M*m)!=0', 'marked_nonpower': 'ALL m>=1: sum_(r=0)^m binom(m,r)*mu_(M*m+r), positive RATIONAL hence nonzero in K', 'absolute_nonpower': 'over Kbar, Fstar != unit*g^e for EVERY e>=2; ord_C1(Fstar)=1 in Kbar[V,q^-1]', 'raw_Milne_required': True, 'old_Hopf_context_only': True, 'prior_cover_or_fulltrace_or_even_packet_required': False, 'ordering_or_whole_field_embedding_required': False, 'norm_one_K_point_required': False, 'anisotropic_odd_forms_claimed': False, 'dimension_two_claimed': False, 'radial_primitivity_or_Image_degree_claimed': False, 'all_phases_or_all_markers_classified': False, 'finite_models_prove_universal_quantifiers': False, 'controls': {'dimensions': [3, 4, 5, 6, 7, 8, 9, 11], 'seed_k_bound': '2*(N-2)+2', 'nonpower_m_max': 3, 'Lie_M': [1, 2, 4], 'Lie_k_max': 4}, 'assurance_flags': {'external_review': False, 'formal_verification': False, 'independent_replay': False, 'worldwide_priority': False}}
FLAGS=SCOPE['assurance_flags']
ROLE={
 'written_proof':(PN,'Complete_Sections1--9_all_field_Lie_endpoint_sharp_iff_marker_and_absolute_nonpower_proof',SCOPE['reviewed_proof_sha256']),
 'exact_checker':(CN,'Source_authentication_and_BOUNDED_exact_sparse_Lie_factor_models_NOT_universal_proof',None),
 'typed_inputs':(IN,'Exact_isotropic_all_char0_mixed_moments_sharp_seed_threshold_nonpower_witness_and_boundaries',None),
 'milne_complete_reducibility_reused':('sources/milne-algebraic-groups-2017.pdf','NECESSARY_REUSED_SRC_RHT_MILNE_COMPLETE_REDUCIBILITY_22_42_22_43_printed477_478_and_12l_printed248_249_same_field_only','503f1a7d511632af05541c52fb2bf6fc36f5a75b0323f2f9110f03dcce11fd67'),
 'hopf_endpoint_context_reused':('HOPF_LIFT_CLASSIFICATION.md','CONTEXT_ONLY_REUSED_SRC_HOPF_LIFT_CLASSIFICATION_Theorem7_1_equations7_1--7_7_M_equals1_endpoint_reproduced_locally_NOT_core_premise','1688857a80e9c09db5727bf66d8e59f85126d19a410db2dbef5c7929bd4f8a60')}
def need(v,msg):
 if not v:raise ValueError(msg)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def typed_equal(a,b):
 if type(a)is not type(b):return False
 if isinstance(a,dict):return a.keys()==b.keys() and all(typed_equal(a[k],b[k])for k in a)
 if isinstance(a,list):return len(a)==len(b)and all(typed_equal(x,y)for x,y in zip(a,b))
 return a==b
def pairs(xs):
 d={}
 for k,v in xs:
  need(k not in d,'Duplicate JSON key: '+k);d[k]=v
 return d
def parse(text):return json.loads(text,object_pairs_hook=pairs,parse_constant=lambda s:(_ for _ in ()).throw(ValueError('Nonfinite JSON: '+s)))
def source_path(name):
 need(type(name)is str and '\\'not in name,'POSIX source path')
 p=PurePosixPath(name)
 need(not p.is_absolute()and '..'not in p.parts and str(p)==name and name not in ('','.'),'Canonical confined relative path')
 q=(HERE/name).resolve()
 need(q.is_relative_to(HERE)and q.is_file(),'Source inside package: '+name)
 return q
def scope_guard(z):need(typed_equal(z,SCOPE),'Exact typed isotropic theorem scope and unchanged assurances')
def ledger_guard(z):
 need(type(z)is dict and set(z)=={'schema','assurance_flags','sources'},'Exact ledger fields')
 need(type(z['schema'])is int and z['schema']==1 and typed_equal(z['assurance_flags'],FLAGS),'Typed ledger schema and assurances')
 need(type(z['sources'])is list,'Ledger sources array')
 seen=set();out=[]
 for row in z['sources']:
  need(type(row)is dict and set(row)=={'role','path','relation','sha256'},'Exact source-row fields')
  role=row['role'];need(type(role)is str and role in ROLE and role not in seen,'Unique known role');seen.add(role)
  path,relation,pin=ROLE[role]
  q=source_path(row['path'])
  need(row['path']==path and row['relation']==relation,'Exact source role path and scope')
  need(type(row['sha256'])is str and re.fullmatch('[0-9a-f]{64}',row['sha256'])and sha(q)==row['sha256'],'Missing or stale source: '+role)
  need(pin is None or row['sha256']==pin,'Exact reviewed source body: '+role)
  out.append(row)
 need(seen==set(ROLE),'All mandatory authentication roles')
 return sorted(out,key=lambda r:r['role'])
# Sparse polynomials, exact rational coefficients, arbitrary finite variable count.
def const(n,c=1):return {(0,)*n:Q(c)}if c else {}
def var(n,i):
 e=[0]*n;e[i]=1;return {tuple(e):Q(1)}
def add(*ps):
 z={}
 for p in ps:
  for e,c in p.items():z[e]=z.get(e,Q(0))+c
 return {e:c for e,c in z.items()if c}
def scale(p,c):return {e:v*c for e,v in p.items()if v*c}
def mul(p,q):
 z={}
 for e,c in p.items():
  for f,d in q.items():
   g=tuple(a+b for a,b in zip(e,f));z[g]=z.get(g,Q(0))+c*d
 return {e:c for e,c in z.items()if c}
def power(p,k,n):
 z=const(n)
 for _ in range(k):z=mul(z,p)
 return z
def diff(p,i):
 z={}
 for e,c in p.items():
  if e[i]:
   f=list(e);f[i]-=1;z[tuple(f)]=c*e[i]
 return z
@lru_cache(None)
def mu(N,k):return Q(factorial(k),prod(N+2*i for i in range(k)))
def binomial(n,r):
 if r<0:return 0
 if n>=0:return comb(n,r)if r<=n else 0
 return (-1)**r*comb(r-n-1,r)
def trace(p,N,j=0):return sum((c*mu(N,b)for(a,b),c in p.items()if a+j==b),Q(0))
def models():
 seed={(3,1):Q(1),(2,1):Q(3),(1,1):Q(3),(0,1):Q(1),(2,0):Q(-1,2),(1,0):Q(-3,2),(0,0):Q(-1)}
 controls=SCOPE['controls'];dims=controls['dimensions'];limit=max(max(2*(N-2)+2,(N-1)*controls['nonpower_m_max'])for N in dims)
 ps=[const(2)]
 for _ in range(limit):ps.append(mul(ps[-1],seed))
 rows=[];comparisons=0
 for N in dims:
  M=N-2;kmax=2*M+2;count=0
  for k in range(1,kmax+1):
   for j in range(k+3):
    want=mu(N,k)*binomial(k-M,k-j)if j<=k else Q(0)
    assert trace(ps[k],N,j)==want;count+=1
   assert (trace(ps[k],N)==0)==(k>=M)
  assert trace(ps[1],N)==Q(3,N)-1
  phase=mul(ps[M],add(const(2),seed));phase_power=const(2);np=[]
  for m in range(1,controls['nonpower_m_max']+1):
   phase_power=mul(phase_power,phase)
   marked=trace(phase_power,N,M);want=sum((comb(m,r)*mu(N,M*m+r)for r in range(m+1)),Q(0))
   assert trace(phase_power,N)==0 and marked==want and marked>0
   assert trace(ps[M*m],N)==0 and trace(ps[M*m],N,M)==mu(N,M*m)
   np.append(dict(m=m,pure='0',marked=str(marked),power_marked=str(mu(N,M*m))))
  rows.append(dict(N=N,M=M,k_max=kmax,marker_comparisons=count,first_seed_moment=str(trace(ps[1],N)),least_pure_seed_power=M,nonpower_phase_moments=np));comparisons+=count
 lie=[]
 for M in controls['Lie_M']:
  n=M+2;X=var(n,0);Y=var(n,1);Z=[var(n,i+2)for i in range(M)];aa=[(-1)**i*(i+1)for i in range(M)]
  q=add(scale(mul(X,Y),2),*(scale(mul(z,z),a)for a,z in zip(aa,Z)))
  for k in range(1,controls['Lie_k_max']+1):
   summed={}
   for i,(a,z)in enumerate(zip(aa,Z)):
    D=lambda f:add(scale(mul(z,diff(f,1)),a),scale(mul(X,diff(f,i+2)),-1))
    assert not D(q)
    f=mul(mul(power(X,k-1,n),power(Y,k,n)),z)
    rhs=add(scale(mul(power(mul(X,Y),k-1,n),mul(z,z)),k*a),scale(power(mul(X,Y),k,n),-1))
    assert D(f)==rhs;summed=add(summed,D(f))
   rhs=add(scale(mul(power(mul(X,Y),k-1,n),q),k),scale(power(mul(X,Y),k,n),-(2*k+M)))
   assert summed==rhs;lie.append(dict(M=M,k=k,remainder='0'))
 n=3;X=var(n,0);Y=var(n,1);R=var(n,2);one=const(n);q=add(scale(mul(X,Y),2),R)
 C0=mul(add(X,one),add(scale(Y,2),scale(mul(add(X,scale(one,2)),R),-1)))
 C1=add(scale(mul(add(scale(X,3),one),Y),2),scale(mul(mul(X,add(X,scale(one,3))),R),-1))
 original=add(scale(mul(power(add(X,one),3,n),Y),2),scale(mul(mul(add(X,one),add(X,scale(one,2))),q),-1))
 assert original==C0 and add(C1,scale(C0,-1))==scale(q,2)
 assert -Q(-1,3)*(Q(-1,3)+3)==Q(8,9)
 return dict(sparse_seed_models=rows,actual_mixed_coefficient_comparisons=comparisons,Lie_numerator_models=lie,factor_identities_for_formal_nonzero_q0=True,primitive_linear_coefficient_at_X_minus_one_third='(8/9)*q0 != 0',negative_control='For N>3, Lambda(p)=3/N-1 !=0; omitting the threshold power fails purity',meaning='BOUNDED rational algebra support only; not an irreducibility algorithm or proof of unrestricted field/dimension/order statements')
def guard_tests(ledger):
 labels=[]
 def reject(label,fn):
  try:fn()
  except (ValueError,TypeError,KeyError):labels.append(label)
  else:raise AssertionError('Bad source/scope accepted: '+label)
 for key in SCOPE:
  z=copy.deepcopy(SCOPE);z.pop(key);reject('missing_scope:'+key,lambda z=z:scope_guard(z))
 for i,row in enumerate(ledger['sources']):
  for change in ['missing','bad_hash','scope']:
   z=copy.deepcopy(ledger)
   if change=='missing':z['sources'].pop(i)
   elif change=='bad_hash':z['sources'][i]['sha256']='0'*64
   else:z['sources'][i]['relation']='context or novelty inflation'
   reject(change+':'+row['role'],lambda z=z:ledger_guard(z))
 for raw in ['{"a":1,"a":2}','{"a":{"b":1,"b":2}}','{"a":NaN}','{"a":Infinity}','{"a":-Infinity}']:reject('json:'+raw,lambda raw=raw:parse(raw))
 z=copy.deepcopy(SCOPE);z['controls']['dimensions'][0]=True;reject('boolean_dimension',lambda:scope_guard(z))
 z=copy.deepcopy(SCOPE);z['controls']['nonpower_m_max']=3.0;reject('floating_integer',lambda:scope_guard(z))
 for name in ['../outside','/tmp/outside','./bad','sources//bad','sources\\bad']:reject('path:'+name,lambda name=name:source_path(name))
 return labels
def main():
 if not __debug__:raise SystemExit('Optimized execution is not accepted.')
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--repo-root',type=Path,default=HERE);ap.add_argument('--output',type=Path,default=HERE/RN);a=ap.parse_args()
 need(a.repo_root.is_dir(),'Existing repo-root directory')
 inputs=parse((HERE/IN).read_text());scope_guard(inputs)
 ledger=parse((HERE/LN).read_text());sources=ledger_guard(ledger)
 protected={source_path(v[0])for v in ROLE.values()}|{(HERE/LN).resolve()}
 need(a.output.resolve()not in protected and not(a.output.exists()and any(a.output.samefile(p)for p in protected)),'Receipt cannot overwrite proof, checker, inputs, ledger or retained source')
 report=dict(schema=1,status='PASS',proof_sha256=sha(HERE/PN),checker_sha256=sha(HERE/CN),inputs_sha256=sha(HERE/IN),source_ledger_sha256=sha(HERE/LN),sources=sources,scope=inputs,models=models(),source_scope_negative_controls=guard_tests(ledger),python_version=sys.version,meaning='Source authentication and declared BOUNDED exact algebra PASS. The written Lie recurrence, endpoint divisibility, rational nonvanishing and valuation proof establish the universal theorem. Finite models do not prove arbitrary fields/dimensions/orders or absolute irreducibility.',assurance_flags=FLAGS)
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps(dict(status='PASS',output=str(a.output),assurance_flags=FLAGS),sort_keys=True))
if __name__=='__main__':main()
