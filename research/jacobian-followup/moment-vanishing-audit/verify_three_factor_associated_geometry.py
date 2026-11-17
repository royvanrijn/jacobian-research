#!/usr/bin/env python3
"""Audit or replay exact isolated-component and associated-prime evidence."""
from pathlib import Path
import argparse,hashlib,json,re,subprocess,sys,time
import sympy as S
P=Path(__file__).resolve().parent
EXPECTED_WICK='7254eff7dad198b2e15fa727da0703bab5f8e3b8704697b27a0ed212527821d6'
def sha(name):return hashlib.sha256((P/name).read_bytes()).hexdigest()
assert sha('verify_nodal_product_reduction.py')==EXPECTED_WICK
from verify_nodal_product_reduction import moments

def blocks_of(log):
 ms=list(re.finditer(r'(?m)^([A-Z][A-Z0-9_]*)$',log))
 blocks={}
 for i,m in enumerate(ms):
  assert m[1] not in blocks
  blocks[m[1]]=log[m.end():ms[i+1].start() if i+1<len(ms) else len(log)].strip()
 return blocks

def main():
 start=time.monotonic();ap=argparse.ArgumentParser();ap.add_argument('--audit-existing',action='store_true');args=ap.parse_args()
 if not args.audit_existing:
  for runner in ['verify_three_factor_row_family_v2.py','run_three_factor_minimal_primary.py','run_three_factor_associated_primes_v2.py']:
   x=subprocess.run([sys.executable,str(P/runner)],text=True,capture_output=True);assert x.returncode==0,(runner,x.stdout,x.stderr)
 runs={};logs={};sources={};provenance={}
 for stem in ['three_factor_minimal_primary','three_factor_associated_primes_v2']:
  rec=json.loads((P/(stem+'_run.json')).read_text());assert rec['terminal'] and rec['completed'] and not rec['timed_out'] and rec['exit_code']==0
  assert rec['input_sha256']==sha(stem+'.sing') and rec['log_sha256']==sha(stem+'_stdout.txt')
  log=(P/(stem+'_stdout.txt')).read_text();assert '?' not in log and log.rstrip().endswith('END')
  runs[stem]=rec;logs[stem]=blocks_of(log);sources[stem]=(P/(stem+'.sing')).read_text()
  for name in [stem+'.sing',stem+'_stdout.txt',stem+'_run.json']:provenance[name]=sha(name)
 row=json.loads((P/'three_factor_row_family_v2_certificate.json').read_text())
 assert row['completed'] and row['terminal'] and row['checker_sha256']==sha('verify_three_factor_row_family_v2.py') and row['input_checker_sha256']==EXPECTED_WICK
 assert [x['free_rank'] for x in row['data']]==[1,2,3,4,4,4,3,2,1,0]
 z=S.symbols('z');bez=row['degree10_Bezout'];assert S.expand(S.sympify(bez['u'])*S.sympify(bez['M2_a_coefficient'])+S.sympify(bez['v'])*S.sympify(bez['M4_a2_coefficient']))==1
 for item in row['data']:assert all(S.degree(S.sympify(v),z)==0 for v in item['Smith_nonzero_diagonal'])
 # Star central-variance-zero fiber: a=0,d=e=1 and eliminate f by M2.
 a,b,c,d,e,f=S.symbols('a b c d e f');raw=moments([[a,d,e],[d,b,f],[e,f,c]],10)
 fiber=[S.expand(v.subs({a:0,d:1,e:1}).subs(f,-(b+c)/4)) for v in raw[1::2]]
 assert fiber[0]==0
 g2,g3=fiber[1:3];res=S.resultant(g2.subs(c,1),g3.subs(c,1),b);assert res!=0 and S.Poly(g2,b,c).coeff_monomial(b*b)!=0
 gb=S.groebner([g2,g3],b,c,domain=S.QQ);assert all(gb.reduce(b**i*c**(4-i))[1]==0 for i in range(5))
 assert all(gb.reduce(v)[1]==0 for v in fiber[3:])
 for v in raw[1::2]:
  m=S.Poly(v,a,b,c,d,e,f).total_degree()//3
  normalized=S.expand(v.subs({a:1,d:1,e:1}).subs({b:a*b,c:a*c,f:a*f},simultaneous=True))
  assert S.expand(normalized-a**m*v.subs({d:1,e:1}))==0
 # Uniform finiteness uses the retained, already proved J^13 containment.
 old=json.loads((P/'three_factor_colon_geometry_verification.json').read_text());assert old['exact_nilradical_indices']['I10_equals_I_infinity']==13
 geo=logs['three_factor_minimal_primary'];soc=logs['three_factor_associated_primes_v2'];t=S.symbols('t')
 def vector(data,name):
  v=[int(x) for x in re.findall(r'-?\d+',data[name])]
  while v and v[-1]==0:v.pop()
  return v
 def h(v):return sum(c*t**k for k,c in enumerate(v))
 assert vector(geo,'ROW_EXPONENT')==[13] and vector(geo,'STAR_EXPONENT')==[11]
 assert vector(geo,'FULL_IN_UNMIXED')==[] # zero is trimmed
 assert vector(geo,'UNMIXED_IN_FULL')==[103]
 full=vector(geo,'FULL_HILBERT');unmixed=vector(geo,'UNMIXED_HILBERT')
 hd=S.cancel((h(full)-h(unmixed))/(1-t)**5);assert S.denom(hd)==1 and hd.subs(t,1)==1617
 dh=[int(S.expand(hd).coeff(t,k)) for k in range(S.degree(hd,t)+1)]
 socle={}
 for label,n,count,dimension in [('PLANE_AD',4,0,0),('AXIS_A',5,13,7),('AXIS_D',5,0,0),('ORIGIN',6,8,8)]:
  assert soc[label+'_SOCLE_GENERATORS']==str(count)
  assert soc[label+'_FORWARD']==soc[label+'_REVERSE']=='0'
  hf=S.cancel((h(vector(soc,label+'_I_HILBERT'))-h(vector(soc,label+'_K_HILBERT')))/(1-t)**n)
  assert S.denom(hf)==1 and hf.subs(t,1)==dimension
  socle[label]={'nonzero_normal_form_generator_count':count,'socle_dimension':dimension,'filtered_Hilbert_polynomial':str(S.expand(hf)),'witness':soc[label+'_WITNESS']}
 assert S.expand(S.sympify(socle['ORIGIN']['filtered_Hilbert_polynomial'])-(5*t**21+3*t**30))==0
 # Check the moment input and the actual quotient/saturation calls being audited.
 minimal_source=sources['three_factor_minimal_primary']
 assert 'ideal I='+','.join(str(v).replace('**','^') for v in raw[1:8:2])+';' in minimal_source
 assert 'poly P='+str(raw[9]).replace('**','^')+';' in minimal_source
 expected=','.join(str(v).replace('**','^') for v in raw[1::2])
 assert sources['three_factor_associated_primes_v2'].count('ideal I='+expected+';')==4
 assert sources['three_factor_associated_primes_v2'].count('ideal K=std(quotient(G,m));')==4
 for label,field,variables in [('PLANE_AD','(0,a,d)','b,c,e,f'),('AXIS_A','(0,a)','b,c,d,e,f'),('AXIS_D','(0,d)','a,b,c,e,f'),('ORIGIN','0','a,b,c,d,e,f')]:
  assert 'ring r'+label+'='+field+',('+variables+'),dp;' in sources['three_factor_associated_primes_v2']
 assert 'sat_with_exp(G,ideal(b*c))' in sources['three_factor_minimal_primary'] and 'sat_with_exp(G,ideal(a*d*e))' in sources['three_factor_minimal_primary']
 # Enumerate candidate coordinate faces after the six primary open sets.
 facets=[set(x) for x in ['bcf','ace','abd','ade','bdf','cef']]
 from itertools import combinations
 faces={''.join(sorted(v)) for F in facets for k in range(4) for v in combinations(sorted(F),k)}
 remaining=sorted(F for F in faces if not any(set(pair)<=set(F) for pair in ['bc','ac','ab','de','df','ef']))
 assert remaining==['','a','ad','ae','b','bd','bf','c','ce','cf','d','e','f']
 for name in ['verify_three_factor_row_family_v2.py','three_factor_row_family_v2_certificate.json','run_three_factor_minimal_primary.py','run_three_factor_associated_primes_v2.py','three_factor_colon_geometry_verification.json']:
  provenance[name]=sha(name)
 cert={'claim':'G3PRODUCTASS20261114','field':'QQ, characteristic zero','associated_primes':{
 'minimal':['(a,d,e)','(b,d,f)','(c,e,f)','(b,c,f)','(a,c,e)','(a,b,d)'],
 'embedded':['(b,c,d,e,f)','(a,c,d,e,f)','(a,b,d,e,f)','(a,b,c,d,e,f)']},
 'row_flat_family':{'base':'QQ[z], b=c=1,f=z; no excluded z','weights':{'a':2,'d':1,'e':1},'ranks':[1,2,3,4,4,4,3,2,1],'rank':24},
 'star_flat_family':{'base':'QQ[a], d=e=1; no excluded a','rank':6,'a_zero_quadratic':str(g2),'a_zero_cubic':str(g3),'binary_resultant':str(res),'uniform_nilpotence_bound_from_J13':13},
 'isolated_primary_components':{'row':'I10:(b*c)^infinity and permutations','star':'I10:(d*e)^infinity = I10:(a*d*e)^infinity and permutations','computed_row_saturation_power':13,'computed_star_ade_saturation_power':11,'row_basis_generators':len(geo['ROW_BASIS'].split(',')),'star_basis_generators':len(geo['STAR_BASIS'].split(',')),'unmixed_basis_generators':len(geo['UNMIXED_BASIS'].split(','))},
 'socle_tests':socle,'unmixed_defect':{'reduced_Hilbert_numerator_over_1_minus_t':dh,'multiplicity':1617,'generic_diagonal_axis_torsion_length':539,'interpretation':'Intrinsic finite-length torsion after localizing at a diagonal-axis prime, not a chosen primary-component length'},
 'remaining_coordinate_faces':remaining,'provenance':provenance,
 'boundary':'Full associated-prime classification and six isolated primary components; explicit axis/origin primary components of I10 remain uncomputed. Do not identify the earlier I8:M10 colon components with these full-ideal components. Internal exact CAS and written proof, no independent replay, external review or formalization.'}
 fpath=P/'three_factor_associated_geometry_certificate.json';fpath.write_text(json.dumps(cert,indent=2)+'\n')
 receipt={'claim':cert['claim'],'completed':True,'associated_prime_count':10,'minimal_primary_components':6,'embedded_primes':4,'generic_diagonal_axis_torsion_length':539,'origin_socle_dimension':8,'certificate_sha256':sha(fpath.name),'checker_sha256':sha(Path(__file__).name),'elapsed_seconds':time.monotonic()-start,'audit_existing_only':args.audit_existing,'assurance':cert['boundary']}
 (P/'three_factor_associated_geometry_verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
