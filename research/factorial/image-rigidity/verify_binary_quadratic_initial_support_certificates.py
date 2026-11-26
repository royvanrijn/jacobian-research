"""Verify exact initial-window certificates for every binary quadratic support.
All computations are over finite prime fields. No Singular or third-party module
is imported. The written argument, not this checker alone, transfers projective
emptiness to characteristic zero. This covers START1 ONLY.
"""
from pathlib import Path
from math import factorial,isqrt
from functools import lru_cache
import argparse,copy,hashlib,json,time
MONS=((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))
LABELS=('1','X','Y','X2','XY','Y2')
def need(ok,msg):
 if not ok:raise ValueError(msg)
def sha(b):return hashlib.sha256(b).hexdigest()
def compositions(n,r):
 if r==0:
  if n==0:yield ()
  return
 if r==1:yield(n,);return
 for i in range(n+1):
  for tail in compositions(n-i,r-1):yield(i,)+tail
def multi(n,a):
 v=factorial(n)
 for j in a:v//=factorial(j)
 return v
@lru_cache(None)
def native_forms(mask):
 """Direct original-support multinomial expansion followed by M1 elimination.
 The eliminated coefficient is -sum L_i u_i; every other is L_pivot*u_i.
 This is an integral linear parameterization of ker(M1), invertible when
 L_pivot (1 or2) is invertible.
 """
 ids=tuple(i for i in range(6) if mask>>i&1)
 need(ids and 0<mask<64,'nonempty support mask')
 pivot=ids[0];rest=ids[1:];r=len(rest)
 L=lambda i:factorial(MONS[i][0])*factorial(MONS[i][1])
 forms=[]
 for j in range(2,len(ids)+1):
  f={}
  for alpha in compositions(j,len(ids)):
   t=alpha[0]
   sx=sum(a*MONS[i][0] for a,i in zip(alpha,ids))
   sy=sum(a*MONS[i][1] for a,i in zip(alpha,ids))
   raw=multi(j,alpha)*factorial(sx)*factorial(sy)*L(pivot)**(j-t)*(-1)**t
   for beta in compositions(t,r):
    exponent=tuple(a+b for a,b in zip(alpha[1:],beta))
    value=raw*multi(t,beta)
    for i,b in zip(rest,beta):value*=L(i)**b
    f[exponent]=f.get(exponent,0)+value
  forms.append({a:v for a,v in f.items() if v})
 return ids,pivot,rest,forms
def form_hash(forms):
 d=[[[v,list(a)] for a,v in sorted(f.items())] for f in forms]
 return sha(json.dumps(d,separators=(',',':')).encode())
def prime(p):
 return isinstance(p,int) and p>12 and all(p%d for d in range(2,isqrt(p)+1))
def verify_case(c):
 mask=c['mask'];ids,pivot,rest,forms=native_forms(mask);s=len(ids);r=s-1
 need(c['support']==[LABELS[i] for i in ids],'support inventory')
 need(c['pivot']==LABELS[pivot] and c['variables']==[LABELS[i] for i in rest],'M1 parameterization inventory')
 need(c['native_form_sha256']==form_hash(forms),'native coefficient hash')
 if s==1:
  need(c['singleton_first_moment_multiplier']==factorial(MONS[pivot][0])*factorial(MONS[pivot][1]),'singleton moment')
  return dict(mask=mask,support_size=s,polynomial_identities=0,terms=0)
 p=c['prime'];need(prime(p),'prime >12')
 D=s*(s-1)//2+1;need(c['power_degree']==D,'exact declared power degree')
 U=c['multipliers'];need(len(U)==r and all(len(row)==r for row in U),'multiplier dimensions')
 base=D+1;weights=[base**i for i in range(r)]
 enc=lambda a:sum(x*w for x,w in zip(a,weights))
 F=[[(enc(a),v%p) for a,v in f.items() if v%p] for f in forms]
 total=0
 for col in range(r):
  result={}
  for row in range(r):
   seen=set()
   for term in U[row][col]:
    need(isinstance(term,list) and len(term)==2,'term shape')
    coeff,a=term
    need(isinstance(coeff,int) and 0<coeff<p,'canonical nonzero coefficient')
    need(isinstance(a,list) and len(a)==r and all(isinstance(x,int) and x>=0 for x in a),'exponent vector')
    need(sum(a)==D-(row+2),'homogeneous multiplier degree')
    aa=tuple(a);need(aa not in seen,'duplicate multiplier monomial');seen.add(aa)
    t=enc(a);total+=1
    for u,v in F[row]:
     key=t+u;result[key]=(result.get(key,0)+coeff*v)%p
  result={k:v for k,v in result.items() if v}
  need(result=={D*weights[col]:1},'pure-power identity failed mask '+str(mask)+' column '+str(col))
 return dict(mask=mask,support_size=s,prime=p,power_degree=D,polynomial_identities=r,terms=total)
def verify(data):
 need(data['scope']=='START1_EVERY_NONEMPTY_BINARY_QUADRATIC_MONOMIAL_SUPPORT','scope')
 need(data['monomial_universe']==list(LABELS),'universe')
 cases=data['cases'];need(len(cases)==63 and sorted(c['mask'] for c in cases)==list(range(1,64)),'all63 supports exactly once')
 return [verify_case(c) for c in cases]
def corruptions(data):
 labels=[]
 def reject(label,mutate):
  bad=copy.deepcopy(data);mutate(bad)
  try:verify(bad)
  except (ValueError,KeyError,TypeError):labels.append(label);return
  raise RuntimeError('corruption accepted '+label)
 full=lambda d:next(c for c in d['cases'] if c['mask']==63)
 reject('missing support face',lambda d:d['cases'].pop())
 reject('incorrect first-moment pivot',lambda d:full(d).update(pivot='X'))
 reject('nonprime modulus',lambda d:full(d).update(prime=15))
 reject('wrong native factorial moment',lambda d:full(d).update(native_form_sha256='0'*64))
 reject('wrong power degree',lambda d:full(d).update(power_degree=15))
 def changed(d):
  t=full(d)['multipliers'][0][0][0];t[0]=t[0]%100+1
 reject('one multiplier coefficient corrupted',changed)
 reject('first-support scope overstated as all-start',lambda d:d.update(scope='ALL_STARTS'))
 return labels
def main():
 if not __debug__:raise SystemExit('ordinary Python required; optimized execution disabled')
 a=argparse.ArgumentParser();a.add_argument('--data',type=Path,default=Path(__file__).with_name('binary_quadratic_initial_support_certificates.json'));a.add_argument('--output',type=Path,required=True);args=a.parse_args()
 need(not args.output.exists(),'output already exists')
 b=args.data.read_bytes();data=json.loads(b);start=time.monotonic();rows=verify(data);negative=corruptions(data)
 result=dict(status='PASS',scope=data['scope'],checker_sha256=sha(Path(__file__).read_bytes()),data_sha256=sha(b),cases=rows,support_count=len(rows),nontrivial_support_count=sum(x['support_size']>1 for x in rows),pure_power_identities=sum(x['polynomial_identities'] for x in rows),multiplier_terms=sum(x['terms'] for x in rows),negative_controls=negative,seconds=time.monotonic()-start,all_coefficients_exact=True,finite_certificate_for_universal_coefficients=True,all_start_SFC_claim=False,actual_dependency_engine_evaluations=0,canonical_mutations=False)
 with args.output.open('x') as h:json.dump(result,h,indent=2);h.write('\n')
 print(json.dumps({k:v for k,v in result.items() if k!='cases'},indent=2))
if __name__=='__main__':main()
