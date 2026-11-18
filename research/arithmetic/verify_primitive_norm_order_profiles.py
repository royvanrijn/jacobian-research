#!/usr/bin/env python3
"""Finite conductor norm/profile audit; no infinite analytic assertion."""
if not __debug__:
 raise RuntimeError("This certificate checker requires active Python assertions; run without -O/-OO.")

import argparse
from fractions import Fraction
from itertools import product
from math import gcd,isqrt,prod
from pathlib import Path
import hashlib,json

def vp_unit(a,p):
 e=0
 while a%p==0:a//=p;e+=1
 return e,a

def legendre(a,p):
 r=pow(a%p,(p-1)//2,p)
 return -1 if r==p-1 else r

def form(d,p,c):
 t=1 if d%4==1 else 0
 u=(1-d)//4 if t else -d//4
 return (1,p**c*t,p**(2*c)*u)

def image(d,p,c,h):
 b=(h+1)//2;mod=p**h;base=p**b;lift=mod//base
 a,t,u=form(d,p,c);cosets=set()
 for x,y in product(range(base),repeat=2):
  if x%p==y%p==0:continue
  g=gcd(gcd(2*a*x+t*y,t*x+2*u*y),lift)
  period=base*g
  cosets.add(((a*x*x+t*x*y+u*y*y)%period,period))
 flags=bytearray(mod)
 for r,period in cosets:
  for n in range(r,mod,period):flags[n]=1
 return flags

def predicted(d,p,c,n):
 if not n:return (d%8==1) if p==2 else (d%p!=0 and legendre(d,p)==1)
 e,u=vp_unit(n,p)
 if p!=2:
  if e<2*c:return e%2==0 and legendre(u,p)==1
  if d%p==0:
   return (e==2*c and legendre(u,p)==1) or (e==2*c+1 and legendre(u,p)==legendre(-d//p,p))
  if legendre(d,p)==-1:return e==2*c
  return e>2*c or (e==2*c and (p!=3 or legendre(u,p)==-1))
 u%=8
 if d%2:
  split=d%8==1
  if c==1:return e==0 or (split and e>=3) or (not split and e==2)
  b=c-1
  if e==0:return u in (1,5) if b==1 else u==1
  if e<2*b:
   if e%2:return False
   k=e//2
   return u==(5 if b-k==1 else 1)
  if e==2*b:return u in (3,7)
  return (split and e>=2*b+3) or (not split and e==2*b+2)
 beta,_=vp_unit(d,2)
 if beta==2:
  eps=-d//4
  if e==0:return u in (1,5) if c==1 else u==1
  if e<2*c:
   return e%2==0 and u==(5 if c-e//2==1 else 1)
  if e==2*c:return u in (1,5)
  return e==2*c+1 and u%4==((1+eps)//2)%4
 eps=-d//8
 if e<2*c:return e%2==0 and u==1
 if e==2*c:return u==(1+2*eps)%8
 return e==2*c+1 and u in (eps%8,(eps+2)%8)

def distinct_prime_factors(n):
    result = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            result.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        result.append(n)
    return result

def jacobi(a, n):
    assert n > 0 and n % 2 == 1
    a %= n
    answer = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                answer = -answer
        a, n = n, a
        if a % 4 == n % 4 == 3:
            answer = -answer
        a %= n
    return answer if n == 1 else 0

def valuation_unit(a, p):
    assert a != 0
    e = 0
    while a % p == 0:
        a //= p
        e += 1
    return e, a

def hilbert_integer(a, d, p):
    assert a > 0 and d < 0
    e, u = valuation_unit(a, p)
    beta, v = valuation_unit(d, p)
    if p == 2:
        eps = lambda z: ((z - 1) // 2) % 2
        omega = lambda z: ((z * z - 1) // 8) % 2
        exponent = eps(u) * eps(v) + e * omega(v) + beta * omega(u)
        return -1 if exponent % 2 else 1
    answer = -1 if (e * beta * ((p - 1) // 2)) % 2 else 1
    if beta % 2:
        answer *= jacobi(u, p)
    if e % 2:
        answer *= jacobi(v, p)
    return answer

def splitting_sign(d, p):
    if d % p == 0:
        return 0
    if p == 2:
        return 1 if d % 8 == 1 else -1
    return jacobi(d, p)

def cyclic_selfsum(flags):
    modulus = len(flags)
    bits = sum(1 << n for n, present in enumerate(flags) if present)
    full_mask = (1 << modulus) - 1
    covered = 0
    for a, present in enumerate(flags):
        if present:
            covered |= ((bits << a) | (bits >> (modulus - a))) & full_mask
    return bytearray((covered >> n) & 1 for n in range(modulus))

def histogram(a_values, b_values, a_coeff, b_coeff, modulus):
    counts = [0] * modulus
    for a in a_values:
        for b in b_values:
            counts[(a_coeff * a + b_coeff * b) % modulus] += 1
    return counts

def order_form(d,f):
 t=1 if d%4==1 else 0
 u=(1-d)//4 if t else -d//4
 return (1,f*t,f*f*u)

def direct_image(d,f,p,h):
 b=(h+1)//2;mod=p**h;base=p**b;lift=mod//base
 a,t,u=order_form(d,f);cosets=set()
 for x,y in product(range(base),repeat=2):
  if x%p==y%p==0:continue
  q=a*x*x+t*x*y+u*y*y
  g=gcd(gcd(2*a*x+t*y,t*x+2*u*y),lift)
  period=base*g
  cosets.add((q%period,period))
 flags=bytearray(mod)
 for r,period in cosets:
  for n in range(r,mod,period):flags[n]=1
 if mod<=256:
  brute=bytearray(mod)
  for x,y in product(range(mod),repeat=2):
   if x%p or y%p:brute[(a*x*x+t*x*y+u*y*y)%mod]=1
  assert brute==flags
 return flags

def local_pred(d,f,p,n):
 c=vp_unit(f,p)[0]
 if c:return predicted(d,p,c,n)
 sign=splitting_sign(d,p)
 if sign==1:return True
 if sign==-1:return n%p!=0
 if n==0:return False
 e,_=vp_unit(n,p)
 return e<=1 and hilbert_integer(n,d,p)==1

def global_flags(d,f,bound):
 a,t,u=order_form(d,f);disc=f*f*d
 flags=bytearray(bound+1)
 ymax=isqrt((4*bound)//(-disc))
 for y in range(-ymax,ymax+1):
  h=isqrt(4*bound+disc*y*y)
  for x in range(-((h+t*y)//2),(h-t*y)//2+1):
   n=a*x*x+t*x*y+u*y*y
   if 0<n<=bound and gcd(x,y)==1:flags[n]=1
 return flags

def generic_table_audit():
 rows=[]
 fields=(-3,-4,-7,-8,-11,-15,-20,-23,-24,-40,-56)
 for d in fields:
  for p,c in [(2,c) for c in range(1,6)]+[(3,1),(3,2),(5,1),(7,1)]:
   h=2*c+(4 if p==2 else 2)
   flags=image(d,p,c,h)
   assert flags==bytearray(predicted(d,p,c,n) for n in range(len(flags)))
   rows.append([d,p,c,h,len(flags),sum(flags),hashlib.sha256(flags).hexdigest()])
 return {'cases':len(rows),'complete_residue_comparisons':sum(r[4] for r in rows),'rows':rows}

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--output',type=Path,help='Write the deterministic JSON receipt to this path.')
 args=parser.parse_args()
 orders=[(-4,2),(-4,3),(-4,4),(-4,5),(-4,6),(-4,15),(-4,25),(-3,2),(-3,4),(-3,8),(-3,6),(-3,9),
         (-7,2),(-7,4),(-7,8),(-7,6),(-7,12),(-8,2),(-8,4),(-8,3),(-8,9),(-8,27),
         (-20,2),(-24,2),(-15,2),(-23,2)]
 records=[];residue_checks=0;profile_total=0;allperiod=0;table_cases=0
 controls=[];control_images={};prefix_label_witness=None
 for d,f in orders:
  bad=distinct_prime_factors(2*f*(-d));level=8*prod(p for p in bad if p!=2)
  caps={}
  exponents={}
  local_targets=[];local_rows=[]
  for p in bad:
   c=vp_unit(f,p)[0];sign=splitting_sign(d,p)
   if c:cap=2*c+2 if sign==1 else 2*c+(sign==0)
   else:cap=1 if (sign==0 or p==2 and sign==1) else 0
   coarse_cap=2*(2*c+vp_unit(d,p)[0])+2
   assert cap<=coarse_cap
   caps[p]=cap
   unit_reps=(1,3,5,7) if p==2 else range(1,p)
   exponents[p]=[e for e in range(cap+1) if any(local_pred(d,f,p,p**e*u) for u in unit_reps)]
   h=(2*c+4 if p==2 else 2*c+2) if c else (7 if p==2 else 3)
   flags=direct_image(d,f,p,h)
   assert flags==bytearray(local_pred(d,f,p,n) for n in range(p**h)),(d,f,p)
   residue_checks+=len(flags);table_cases+=1
   sums=cyclic_selfsum(flags)
   coarse=p**min(h,cap+(3 if p==2 else 1))
   assert all(sums[n]==sums[n%coarse] for n in range(len(sums)))
   local_targets.append((coarse,sums[:coarse]))
   unit_image=[u for u in unit_reps if flags[u]]
   assert all(local_pred(d,f,p,a*b) for a,b in product(unit_image,repeat=2))
   local_rows.append({'prime':p,'conductor_exponent':c,'image_modulus':p**h,'image_count':sum(flags),
                      'cap_for_additive_profiles':cap,'coarse_adjugate_Hensel_cap':coarse_cap,'nonempty_valuation_layers':exponents[p],
                      'unit_square_class_representatives':unit_image,
                      'direct_primitive_image_sha256':hashlib.sha256(flags).hexdigest(),
                      'pair_target_coarse_modulus':coarse,'local_pair_target_count':sum(sums[:coarse])})
   if (d,f,p) in [(-8,3,3),(-4,3,3),(-4,2,2),(-7,2,2)]:control_images[(d,f,p)]=flags
  prefixes=sorted(prod(p**e for p,e in zip(bad,es)) for es in product(*(exponents[p] for p in bad)))
  z=prod(p**caps[p] for p in bad);period=z*level
  units=[a for a in range(1,level) if gcd(a,level)==1]
  state=lambda a:(a%8,tuple(legendre(a,p) for p in bad if p!=2))
  ar={}
  for r in prefixes:
   by_state={}
   ar[r]=[]
   for a in units:
    good=all(local_pred(d,f,p,r*a) for p in bad)
    if state(a) in by_state:assert by_state[state(a)]==good
    by_state[state(a)]=good
    if good:
     assert jacobi(d,a)==1
     ar[r].append(a)
  assert all(ar[r] for r in prefixes)
  profile_mask=bytearray(period);positive_profiles=0;profiles=0
  for r in prefixes:
   for s in prefixes:
    g=gcd(r,s);a,b=r//g,s//g
    counts=histogram(ar[r],ar[s],a,b,level)
    if prefix_label_witness is None:
     wrong_ar=[u for u in units if all(local_pred(d,f,p,a*u) for p in bad)]
     wrong_bs=[u for u in units if all(local_pred(d,f,p,b*u) for p in bad)]
     wrong=histogram(wrong_ar,wrong_bs,a,b,level)
     if wrong!=counts:
      t=next(t for t in range(level) if wrong[t]!=counts[t])
      prefix_label_witness=[d,f,r,s,g,t,counts[t],wrong[t]]
    if any(counts):positive_profiles+=1
    for t,count in enumerate(counts):
     if count:
      for n in range(g*t,period,g*level):profile_mask[n]=1
    profiles+=1
  direct_mask=bytearray(all(flags[n%mod] for mod,flags in local_targets) for n in range(period))
  assert profile_mask==direct_mask,(d,f)
  # Exact full-mask equality must reject a deleted representation and a
  # spurious locally excluded target, not just compare scalar counts.
  admitted=next(i for i,x in enumerate(direct_mask) if x)
  excluded=next(i for i,x in enumerate(direct_mask) if not x)
  erased=profile_mask[:];erased[admitted]=0
  injected=profile_mask[:];injected[excluded]=1
  assert erased!=direct_mask and injected!=direct_mask
  profile_total+=profiles;allperiod+=period
  primitive=global_flags(d,f,4096)
  for n,present in enumerate(primitive):
   if not present:continue
   assert all(local_pred(d,f,p,n) for p in bad)
   for q in distinct_prime_factors(n):
    if q not in bad:assert splitting_sign(d,q)==1
  records.append({'fundamental_discriminant':d,'conductor':f,'order_discriminant':f*f*d,
                  'principal_norm_form':order_form(d,f),'unit_state_period_L':level,
                  'complete_target_period':period,'local_layers':local_rows,'prefix_count':len(prefixes),
                  'all_profile_count':profiles,'positive_profile_count':positive_profiles,
                  'locally_admissible_target_residues':sum(profile_mask),
                  'profile_and_direct_target_sha256':hashlib.sha256(profile_mask).hexdigest(),
                  'bounded_global_primitive_values_4096':sum(primitive),
                  'bounded_global_primitive_flags_sha256':hashlib.sha256(primitive).hexdigest()})
 # Purposeful wrong-conductor rules challenged by independent coordinates.
 a=control_images[(-8,3,3)]
 assert not a[9] and a[18] and a[27]
 controls.append('split3 threshold square-unit exception cannot be discarded')
 a=control_images[(-4,3,3)]
 assert a[9]
 controls.append('inert conductor prime factors cannot be excluded as if the order were maximal')
 a=control_images[(-4,2,2)]
 assert a[4] and a[8] and not a[16]
 controls.append('maximal-order dyadic criterion cannot replace the conductor criterion')
 a=control_images[(-7,2,2)]
 assert not a[2] and not a[4] and a[8] and a[16] and a[32]
 controls.append('split conductor image has a high-valuation tail, not a universal valuation bound')
 assert prefix_label_witness is not None
 assert prefix_label_witness[-2]!=prefix_label_witness[-1]
 controls.append('cancelling a common prefix must not change the ORIGINAL local norm labels')
 for row in records:
  assert row['all_profile_count']==row['prefix_count']**2
 # A norm value need not determine a unique genus at conductor primes.
 forms=((1,0,25),(2,2,13));witnesses=((0,1),(2,1))
 assert all(b*b-4*a*c==-100 for a,b,c in forms)
 assert all(gcd(x,y)==1 and a*x*x+b*x*y+c*y*y==25
            for (a,b,c),(x,y) in zip(forms,witnesses))
 unit_images=[]
 for a,b,c in forms:
  unit_images.append(sorted({(a*x*x+b*x*y+c*y*y)%5
                             for x,y in product(range(5),repeat=2)
                             if (a*x*x+b*x*y+c*y*y)%5}))
 assert unit_images==[[1,4],[2,3]]
 controls.append('the conductor norm value25 does not determine a unique genus at discriminant-100')
 multiple_genus={'discriminant':-100,'value':25,'forms':forms,
                 'primitive_witnesses':witnesses,'unit_images_mod5':unit_images,
                 'meaning':'The two forms have different local unit-value images at5, but both primitively represent25.'}
 report={'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'scope':'Finite exact conductor tables, complete profile/direct-local target masks, square-class-state checks, and bounded coordinate inclusion only. Infinite arithmetic and analytic proofs are separate.',
         'generic_conductor_tables':generic_table_audit(),
         'multiple_genus_conductor_example':multiple_genus,
         'orders':records,'order_count':len(records),'local_image_cases':table_cases,
         'direct_image_residue_comparisons':residue_checks,'profile_count':profile_total,
         'complete_period_target_residue_comparisons':allperiod,'negative_controls':controls,
         'full_mask_corruptions_rejected':2*len(records),'wrong_prefix_label_witness':prefix_label_witness}
 encoded=json.dumps(report,indent=2,sort_keys=True)+'\n'
 if args.output:
  args.output.write_text(encoded,encoding='utf-8')
  summary={k:report[k] for k in ['order_count','local_image_cases','direct_image_residue_comparisons','profile_count','complete_period_target_residue_comparisons']}
  summary['checker_sha256']=report['checker_sha256']
  summary['receipt_sha256']=hashlib.sha256(encoded.encode()).hexdigest()
  print(json.dumps(summary,sort_keys=True))
 else:
  print(encoded,end='')

if __name__=='__main__':main()
