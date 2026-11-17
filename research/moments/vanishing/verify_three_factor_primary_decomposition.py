#!/usr/bin/env python3
"""Exact Hilbert certificates for a constructive ten-component decomposition."""
from pathlib import Path
import argparse,hashlib,json,re,subprocess,sys,time
import sympy as S
P=Path(__file__).resolve().parent
EXPECTED_ASSOC='a6fd844bfdb11c64b64cfab60d68cbfbd8a12d62ec70d1dbf8cc59549c9a8f6c'
EXPECTED_WICK='7254eff7dad198b2e15fa727da0703bab5f8e3b8704697b27a0ed212527821d6'
def sha(n):return hashlib.sha256((P/n).read_bytes()).hexdigest()
assert sha('three_factor_associated_geometry_certificate.json')==EXPECTED_ASSOC
assert sha('verify_nodal_product_reduction.py')==EXPECTED_WICK
from verify_nodal_product_reduction import moments

def blocks_of(log):
 ms=list(re.finditer(r'(?m)^([A-Z][A-Z0-9_]*)$',log));out={}
 for i,m in enumerate(ms):
  assert m[1] not in out
  out[m[1]]=log[m.end():ms[i+1].start() if i+1<len(ms) else len(log)].strip()
 return out

def homogeneous_degree(raw,weights):
 degrees=set()
 for token in re.findall(r'[+-]?[^+-]+',raw.strip().strip(',').replace(' ','')):
  token=token.lstrip('+-');num=re.match(r'\d+(?:/\d+)?',token);rest=token[len(num[0]):] if num else token;deg=0
  while rest:
   m=re.match(r'([a-f])(\d*)',rest);assert m,(raw,rest)
   deg+=weights[m[1]]*(int(m[2]) if m[2] else 1);rest=rest[len(m[0]):]
  degrees.add(deg)
 assert len(degrees)==1,(raw,degrees)
 return next(iter(degrees))

def main():
 start=time.monotonic();ap=argparse.ArgumentParser();ap.add_argument('--audit-existing',action='store_true');args=ap.parse_args()
 if not args.audit_existing:
  for runner in ['run_three_factor_axis_weighted.py','run_three_factor_origin_coordinate.py']:
   x=subprocess.run([sys.executable,str(P/runner)],text=True,capture_output=True);assert x.returncode==0,(runner,x.stdout,x.stderr)
 logs={};sources={};provenance={}
 for stem in ['three_factor_axis_weighted','three_factor_origin_coordinate']:
  rec=json.loads((P/(stem+'_run.json')).read_text());assert rec['terminal'] and rec['completed'] and rec['exit_code']==0 and not rec['timed_out']
  assert rec['input_sha256']==sha(stem+'.sing') and rec['log_sha256']==sha(stem+'_stdout.txt')
  log=(P/(stem+'_stdout.txt')).read_text();assert '?' not in log and log.rstrip().endswith('END');logs[stem]=blocks_of(log);sources[stem]=(P/(stem+'.sing')).read_text()
  for n in [stem+'.sing',stem+'_stdout.txt',stem+'_run.json','run_'+stem+'.py']:provenance[n]=sha(n)
 oldlog=blocks_of((P/'three_factor_minimal_primary_stdout.txt').read_text())
 assert 'ideal U='+oldlog['UNMIXED_BASIS']+';' in sources['three_factor_axis_weighted']
 assert 'map deh=r,1,b,c,d,e,f;' in sources['three_factor_axis_weighted']
 assert 'ring r1=0,(b,c,d,e,f),wp(2,2,1,1,2);' in sources['three_factor_axis_weighted']
 assert 'intvec w=2,2,1,1,2;' in sources['three_factor_axis_weighted']
 origin_source=sources['three_factor_origin_coordinate']
 for v in ['a','b','c']:assert 'std(sat(G,ideal('+v+')))' in origin_source
 assert 'std(intersect(intersect(Ga,Gb),Gc))' in origin_source
 a,b,c,d,e,f=S.symbols('a b c d e f');raw=moments([[a,d,e],[d,b,f],[e,f,c]],10)
 expectedI='ideal I='+','.join(str(raw[k-1]).replace('**','^') for k in [2,4,6,8])+';'
 expectedP='poly P='+str(raw[9]).replace('**','^')+';'
 for source in sources.values():assert expectedI in source and expectedP in source
 w={'b':2,'c':2,'d':1,'e':1,'f':2}
 for k in range(1,6):
  q=S.Poly(raw[2*k-1].subs(a,1),b,c,d,e,f)
  assert all(sum(n*w[v.name] for n,v in zip(exp,q.gens))==4*k for exp,coef in q.terms())
 for label in ['DEHOMOG_FULL_BASIS','DEHOMOG_UNMIXED_BASIS']:
  for poly in logs['three_factor_axis_weighted'][label].split(','):homogeneous_degree(poly,w)
 for poly in logs['three_factor_origin_coordinate']['ORIGIN_SAT_BASIS'].split(','):homogeneous_degree(poly,dict.fromkeys('abcdef',1))
 t=S.symbols('t')
 def vector(log,name):
  v=[int(x) for x in re.findall(r'-?\d+',log[name])]
  while v and v[-1]==0:v.pop()
  return v
 def h(v):return sum(c*t**k for k,c in enumerate(v))
 def coeffs(poly):
  assert S.denom(poly)==1
  out=[int(S.expand(poly).coeff(t,k)) for k in range(S.degree(poly,t)+1)];assert all(v>=0 for v in out)
  return out
 axis=logs['three_factor_axis_weighted'];org=logs['three_factor_origin_coordinate']
 ai=vector(axis,'DEHOMOG_FULL_WEIGHTED_HILBERT');au=vector(axis,'DEHOMOG_UNMIXED_WEIGHTED_HILBERT')
 hd=S.cancel((h(ai)-h(au))/((1-t)**2*(1-t*t)**3));dh=coeffs(hd)
 assert sum(dh)==539 and len(dh)==34 and next(k for k,v in enumerate(dh) if v)==16
 oi=vector(org,'FULL_HILBERT');os=vector(org,'ORIGIN_SAT_HILBERT');ht=S.cancel((h(oi)-h(os))/(1-t)**6);th=coeffs(ht)
 assert sum(th)==77 and len(th)==31 and next(k for k,v in enumerate(th) if v)==18
 assert th[18:]==[9,15,18,8]+[3]*9
 assert oi==vector(oldlog,'FULL_HILBERT')
 assoc=json.loads((P/'three_factor_associated_geometry_certificate.json').read_text());assert assoc['unmixed_defect']['generic_diagonal_axis_torsion_length']==sum(dh)
 axis_powers={v:33//wt+1 for v,wt in w.items()};assert axis_powers=={'b':17,'c':17,'d':34,'e':34,'f':17}
 assert all(axis_powers[v]*wt>33 for v,wt in w.items())
 components=[{'name':'row_A','prime':'(a,d,e)','ideal':'I:(b*c)^infinity'}, {'name':'row_B','prime':'(b,d,f)','ideal':'I:(a*c)^infinity'}, {'name':'row_C','prime':'(c,e,f)','ideal':'I:(a*b)^infinity'}, {'name':'star_A','prime':'(b,c,f)','ideal':'I:(d*e)^infinity'}, {'name':'star_B','prime':'(a,c,e)','ideal':'I:(d*f)^infinity'}, {'name':'star_C','prime':'(a,b,d)','ideal':'I:(e*f)^infinity'}, {'name':'axis_A','prime':'(b,c,d,e,f)','ideal':'(I+(b^17,c^17,d^34,e^34,f^17)):a^infinity'}, {'name':'axis_B','prime':'(a,c,d,e,f)','ideal':'(I+(a^17,c^17,d^34,e^17,f^34)):b^infinity'}, {'name':'axis_C','prime':'(a,b,d,e,f)','ideal':'(I+(a^17,b^17,d^17,e^34,f^34)):c^infinity'}, {'name':'origin','prime':'(a,b,c,d,e,f)','ideal':'I+(a^31,b^31,c^31,d^31,e^31,f^31)'}]
 assert {x['prime'] for x in components}==set(assoc['associated_primes']['minimal']+assoc['associated_primes']['embedded'])
 retained_failures={}
 for stem in ['three_factor_origin_weighted','three_factor_embedded_primary']:
  rr=json.loads((P/(stem+'_run.json')).read_text());assert rr['terminal'] and rr['timed_out'] and not rr['completed']
  assert rr['input_sha256']==sha(stem+'.sing') and rr['log_sha256']==sha(stem+'_stdout.txt')
  retained_failures[stem]='Terminal50-second timeout; not used as a primary-ideal or intersection certificate.'
 for n in ['verify_nodal_product_reduction.py','three_factor_associated_geometry_certificate.json','three_factor_minimal_primary_stdout.txt']:provenance[n]=sha(n)
 cert={'claim':'G3PRODUCTPRIMARY20261114','field':'QQ, characteristic zero','ideal':'I=I_infinity=I10=(M2,M4,M6,M8,M10)',
 'primary_decomposition':components,'intersection':'I is the irredundant intersection of these ten explicitly defined ideals; saturation notation is retained, expanded embedded standard bases are not asserted.',
 'origin_torsion':{'definition':'H^0_m(R/I)=S/I, S=intersect(I:a^infinity,I:b^infinity,I:c^infinity)','Hilbert_function':th,'length':77,'first_degree':18,'last_degree':30},
 'diagonal_axis_torsion':{'normalization':'a=1 is an exact homogeneous localization reparametrization, not a sampled coefficient specialization','weights':w,'Hilbert_function':dh,'length':539,'first_weight':16,'last_weight':33,'truncation_powers':axis_powers},
 'Hilbert_numerators':{'axis_full':ai,'axis_unmixed':au,'origin_full':oi,'origin_saturated':os},
 'proof_mechanism':'All relation ideals are homogeneous for the asserted grading. A kernel supported in degrees<=B meets any ideal generated in degrees>B trivially. Apply this to the localized axis defect and to intrinsic origin torsion. The ten-associated-prime theorem and six known isolated primary components exclude any other support; this gives the written full intersection identity and irredundancy.',
 'provenance':provenance,'retained_failures':retained_failures,'boundary':'Constructive written full primary decomposition with exact Hilbert/saturation inputs; not a claim of expanded embedded Groebner bases or a raw full-intersection replay. Intrinsic torsion lengths77 and539 are distinct from noncanonical chosen primary-component lengths. No constrained-pullback or positive-characteristic conclusion, independent replay, external review or formal verification.'}
 fpath=P/'three_factor_primary_decomposition_certificate.json';fpath.write_text(json.dumps(cert,indent=2)+'\n')
 receipt={'claim':cert['claim'],'completed':True,'primary_components':10,'origin_torsion_length':77,'generic_axis_torsion_length':539,'highest_degrees':{'origin':30,'weighted_axis':33},'checker_sha256':sha(Path(__file__).name),'certificate_sha256':sha(fpath.name),'audit_existing_only':args.audit_existing,'elapsed_seconds':time.monotonic()-start,'assurance':cert['boundary']}
 (P/'three_factor_primary_decomposition_verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
