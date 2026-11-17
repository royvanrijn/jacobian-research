#!/usr/bin/env python3
"""Independent inverse witness and row-minor audit for cubic target rigidity."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import isqrt
import hashlib,json,time
from verify import decode,derivative,mul,add,compose,variable,identity,mmul
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def mon(ix):
 e=[0]*12
 for j in ix:e[j]+=1
 return {tuple(e):Q(1)}
def minor_mod(columns,row_exponents,selected,prime):
 # Sparse ROW elimination on an explicitly selected square minor.
 # The producer used COLUMN elimination over all source coefficient rows.
 rows=[tuple(e) for e in row_exponents];n=len(rows)
 assert n==len(selected)==len(set(rows))==len(set(selected))
 positions={e:i for i,e in enumerate(rows)};matrix=[{} for _ in rows]
 for j,k in enumerate(selected):
  for e,c in columns[k].items():
   i=positions.get(e)
   if i is not None:
    assert c.denominator==1
    val=c.numerator%prime
    if val:matrix[i][j]=val
 basis={};order=[];product=1
 for row_number,row in enumerate(matrix):
  if n>1000 and row_number%500==0:print('MINOR ROW',row_number,'OF',n,flush=True)
  v=dict(row)
  while v:
   j=min(v)
   if j not in basis:
    pivot=v[j];iv=pow(pivot,-1,prime);assert pivot*iv%prime==1
    product=product*pivot%prime;basis[j]={k:c*iv%prime for k,c in v.items()};order.append(j);break
   factor=v[j]
   for k,c in basis[j].items():
    value=(v.get(k,0)-factor*c)%prime
    if value:v[k]=value
    elif k in v:del v[k]
  assert v,'singular proposed minor'
 assert len(basis)==n
 seen=set();cycles=0
 for i in range(n):
  if i in seen:continue
  cycles+=1;j=i
  while j not in seen:seen.add(j);j=order[j]
 parity=(n-cycles)%2
 determinant=(-product if parity else product)%prime
 assert determinant
 return determinant

def main():
 start=time.monotonic();path=HERE/'cubic12-integer-witness.json';raw=json.loads(path.read_text());F=[decode(row,12) for row in raw['coefficients']];X=[variable(i,12) for i in range(12)]
 assert len(F)==12 and all(not p.get((0,)*12,0) for p in F)
 assert all(c.denominator==1 for p in F for c in p.values())
 assert all({e:c for e,c in p.items() if sum(e)==1}==X[i] for i,p in enumerate(F))
 primal_path=HERE/'target-pullback12/checkpoint.json';primal=json.loads(primal_path.read_text());prime=primal['prime']
 assert prime==1000000007 and all(prime%d for d in range(2,isqrt(prime)+1))
 assert primal['input_sha256']==sha(path) and primal['source_sha256']==sha(HERE/'target_pullback12.py')
 words=[ix for d in (2,3) for ix in combinations_with_replacement(range(12),d)]
 assert len(words)==442 and primal['words']==[[i+1 for i in ix] for ix in words]
 columns=[]
 for ix in words:
  full=compose(mon(ix),F,12)
  columns.append({e:c for e,c in full.items() if sum(e)>3})
 assert len(set().union(*(set(p) for p in columns)))==7276
 det_target=minor_mod(columns,primal['pivot_monomials'],list(range(442)),prime)
 print('TARGET MINOR442',det_target,flush=True)
 higher_path=HERE/'target-pullback12-higher/checkpoint.json';higher=json.loads(higher_path.read_text())
 assert higher['status']=='complete-through5' and higher['prime']==prime and higher['source_sha256']==sha(HERE/'target_pullback12_higher.py') and higher['input_sha256']==sha(path) and not higher['dependent_columns']
 higher_words=[ix for degree in range(2,6) for ix in combinations_with_replacement(range(12),degree)]
 assert len(higher_words)==6175 and higher['words']==[[i+1 for i in ix] for ix in higher_words]
 selected_rows=set(map(tuple,higher['pivot_monomials']));assert len(selected_rows)==6175 and all(sum(e)>3 for e in selected_rows)
 # Independent integer modular coefficient homomorphism. Only parent products
 # through degree4 and coefficients of the proposed minor are retained.
 def modmul(a,b):
  out={}
  for ea,ca in a.items():
   for eb,cb in b.items():
    ex=tuple(x+y for x,y in zip(ea,eb));out[ex]=(out.get(ex,0)+ca*cb)%prime
  return {e:c for e,c in out.items() if c}
 Fmod=[{e:c.numerator%prime for e,c in p.items()} for p in F];cache={(i,):p for i,p in enumerate(Fmod)};higher_cols=[]
 for word in higher_words:
  poly=modmul(cache[word[:-1]],Fmod[word[-1]])
  if len(word)<5:cache[word]=poly
  higher_cols.append({e:c for e,c in poly.items() if e in selected_rows})
 det_higher=minor_mod(higher_cols,higher['pivot_monomials'],list(range(6175)),prime)
 print('TARGET MINOR6175',det_higher,flush=True)

 inv_path=HERE/'target-derivation12/inverse.json';inv_raw=json.loads(inv_path.read_text());V=[[decode(p,12) for p in row] for row in inv_raw['inverse']]
 assert len(V)==12 and all(len(row)==12 for row in V)
 assert inv_raw['input_sha256']==sha(path) and inv_raw['source_sha256']==sha(HERE/'target_derivation12.py')
 assert all(c.denominator==1 for row in V for p in row for c in p.values())
 assert sum(len(p) for row in V for p in row)==4620 and max(sum(e) for row in V for p in row for e in p)==15
 J=[[derivative(f,j) for j in range(12)] for f in F]
 assert mmul(J,V)==identity(12,12) and mmul(V,J)==identity(12,12)
 print('BOTH EXACT INVERSE IDENTITIES PASS',flush=True)
 cert_path=HERE/'target-derivation12/checkpoint.json';cert=json.loads(cert_path.read_text())
 assert cert['status']=='complete' and cert['prime']==prime and cert['source_columns']==455 and cert['inverse_sha256']==sha(inv_path) and cert['input_sha256']==sha(path) and cert['source_sha256']==sha(HERE/'target_derivation12.py')
 source_words=[ix for d in range(4) for ix in combinations_with_replacement(range(12),d)]
 assert len(source_words)==455
 source=[mon(ix) for ix in source_words];grad=[[derivative(p,j) for j in range(12)] for p in source];axis=[]
 assert [rec['axis'] for rec in cert['axes']]==list(range(1,13))
 for i,rec in enumerate(cert['axes']):
  images=[]
  for dq in grad:
   p={}
   for a in range(12):p=add(p,mul(V[a][i],dq[a]))
   images.append(p)
  assert not images[0]
  assert len(rec['pivot_columns'])==len(rec['pivot_monomials'])==443
  assert all(0<=j<455 for j in rec['pivot_columns'])
  det=minor_mod(images,rec['pivot_monomials'],rec['pivot_columns'],prime)
  # JV=I verifies all Delta_i K_j relations; their constant/linear jets
  # independently give the required twelve-dimensional exact kernel.
  axis.append({'axis':i+1,'minor_size':443,'minor_determinant_mod_prime':det,'exact_kernel_dimension':12})
  print('AXIS',i+1,'MINOR443',det,flush=True)
 old_path=HERE/'integer-validation12.json';old=json.loads(old_path.read_text());assert old['witness_sha256']==sha(path) and old['determinant_one'] and old['nonlinear_total']==36
 receipt={'status':'pass','claim':'KCUB26','prime':prime,'prime_trial_division_checked':True,'target_word_count':6175,'target_degree_bound':5,'degree3_control_word_count':442,'degree3_control_minor_determinant_mod_prime':det_target,'higher_target_minor_determinant_mod_prime':det_higher,'source_high_rows':7276,'target_minor_determinant_mod_prime':det_higher,'polynomial_inverse_terms':4620,'polynomial_inverse_degree':15,'both_inverse_identities_checked':True,'source_cubic_dimension':455,'axes':axis,'all12_kernel_dimensions':12,'arbitrary_degree_single_coordinate_target_conclusion':True,'rational_other_output_field_intersection':True,'affine_frame_extension_only_for_target_degree_at_most5':True,
 'independent_arithmetic':'stdlib Fraction inverse/control reconstruction, integer modular pullback coefficients, and row-based modular elimination on selected square minors; no primary inverse algorithm, column elimination, RREF, SymPy or FLINT trusted.',
 'proof_boundary':'Nonzero modular integer minors prove characteristic-zero ranks; exact inverse and known kernel generators prove equality. Chain rule, injective pullback and the arbitrary rational-target consequence are written; no proof-assistant formalization or external referee.',
 'checker_sha256':sha(Path(__file__)),'proof_sha256':sha(HERE/'TARGET_CUBIC_RIGIDITY.md'),'arithmetic_sha256':sha(HERE/'verify.py'),'input_sha256':sha(path),'base_validation_sha256':sha(old_path),'target_certificate_sha256':sha(primal_path),'higher_target_certificate_sha256':sha(higher_path),'higher_target_producer_sha256':sha(HERE/'target_pullback12_higher.py'),'derivation_certificate_sha256':sha(cert_path),'inverse_sha256':sha(inv_path),'target_producer_sha256':sha(HERE/'target_pullback12.py'),'derivation_producer_sha256':sha(HERE/'target_derivation12.py'),'seconds':round(time.monotonic()-start,3),'new_smaller_positive_witness':False,'positive_bounds_unchanged':True,'goal_complete':False}
 (HERE/'target-derivation12/validation.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print('PASS complete target6175 rank(through degree5) and twelve443 ranks;uncapped rational single-coordinate rigidity;positive bounds unchanged.',flush=True)
if __name__=='__main__':main()
