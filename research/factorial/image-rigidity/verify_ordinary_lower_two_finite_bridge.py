#!/usr/bin/env python3
"""Complete exact actual lower-two ordinary consecutive finite bridge."""
if not __debug__:raise SystemExit('Refusing -O/-OO: exact verification is required')
from pathlib import Path
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor,as_completed
from hashlib import sha256
from math import comb,factorial,isqrt
import copy,json,struct,sys
import flint
from flint import nmod_poly
HERE=Path(__file__).resolve().parent
PROOF='ORDINARY_LOWER_TWO_STRONG_FACTORIAL.md'
INPUTS='ordinary_lower_two_finite_bridge_inputs.json'
LEDGER='ordinary_lower_two_finite_bridge_source_ledger.json'
RECEIPT='ordinary_lower_two_finite_bridge_verification.json'
def need(q,m):
 if not q:raise ValueError(m)
def digest(p):return sha256(p.read_bytes()).hexdigest()
def domain(d):
 need(set(d)=={'schema_version','lower_exponent','prime','minimum_gap','maximum_gap','minimum_terminal','tail_minimum','tail_linear_factor','expected_cases','gaps'},'wrong input schema')
 need([d[k]for k in ['schema_version','lower_exponent','prime','minimum_gap','maximum_gap','minimum_terminal','tail_minimum','tail_linear_factor','expected_cases']]==[1,2,1000003,6,45,2,4000,128,176816],'wrong mathematical domain')
 p=d['prime'];need(all(p%j for j in range(2,isqrt(p)+1)),'nonprime modulus')
 need(len(d['gaps'])==40,'missing complete gap family')
 rows=[]
 for g,r in zip(range(6,46),d['gaps']):
  last=max(4000,128*(g+2))-1
  need(set(r)=={'gap','first_terminal','last_terminal','tested','case_stream_sha256'},'bad gap schema')
  need([r[k]for k in ['gap','first_terminal','last_terminal','tested']]==[g,2,last,last-1],'wrong actual endpoints')
  need(isinstance(r['case_stream_sha256'],str)and len(r['case_stream_sha256'])==64,'bad case digest')
  need((g+2)*last<p,'factorial specialization loses degree')
  rows.append((g,last,p,r['case_stream_sha256']))
 need(sum(n-1 for g,n,p,h in rows)==176816 and max(n for g,n,p,h in rows)==6015,'incomplete box')
 need(max((g+2)*n for g,n,p,h in rows)==282705,'wrong maximal factorial')
 return rows
def tables(p,n):
 need(n<p,'nonunit factorial')
 f=[1]*(n+1)
 for j in range(1,n+1):f[j]=f[j-1]*j%p
 v=[1]*(n+1);v[n]=pow(f[n],-1,p)
 for j in range(n,0,-1):v[j-1]=v[j]*j%p
 return f,v
def coefficients(n,g,p,tab):
 f,v=tab;pref=f[n]*v[2*n]%p
 return[pref*f[2*n+g*k]*v[k]*v[n-k]%p for k in range(n+1)]
def coprime(a,b,p,N):
 f=nmod_poly(a,p);g=nmod_poly(b,p)
 need(f.degree()==N-1 and g.degree()==N,'actual pair degrees were lost')
 h=f.gcd(g);need(h.degree()==0 and int(h[0])==1,'actual pair has no gcd-one certificate')
def worker(row):
 g,last,p,expected=row;tab=tables(p,(g+2)*last);h=sha256()
 for N in range(2,last+1):
  lo=coefficients(N-1,g,p,tab);hi=coefficients(N,g,p,tab)
  need(lo[0]==hi[0]==1 and lo[-1]!=0 and hi[-1]!=0,'actual endpoint normalization lost')
  coprime(lo,hi,p,N);h.update(struct.pack('<IIII',g,N,p,0))
 need(h.hexdigest()==expected,'complete case stream differs')
 return{'gap':g,'first_terminal':2,'last_terminal':last,'tested':last-1,'prime':p,'all_full_degree_actual_pair_gcds_equal_one':True,'case_stream_sha256':h.hexdigest()}
def controls(inp):
 p=inp['prime'];literal=0
 for g in[1,2,6,45]:
  tab=tables(p,(g+2)*9)
  for n in range(1,9):
   c=coefficients(n,g,p,tab)
   need(c==[comb(n,k)*(factorial(2*n+g*k)//factorial(2*n))%p for k in range(n+1)],'wrong physical Gamma shape')
   literal+=1
 tab=tables(p,64);n=4;g=6
 actual=coefficients(n,g,p,tab)
 wrong=[comb(n,k)*(factorial(2*n-1+g*k)//factorial(2*n-1))%p for k in range(n+1)]
 need(wrong!=actual and wrong[1]!=actual[1],'incorrect Gamma shift was undetected')
 a=coefficients(3,g,p,tab);b=coefficients(4,g,p,tab)
 a[1]=-(a[0]+sum(a[2:]))%p;b[1]=-(b[0]+sum(b[2:]))%p
 try:coprime(a,b,p,4)
 except ValueError:pass
 else:raise ValueError('actual full-degree common-root corruption accepted')
 bads=[]
 x=copy.deepcopy(inp);x['gaps'].pop();bads.append(x)
 x=copy.deepcopy(inp);x['gaps'][-1]['last_terminal']-=1;bads.append(x)
 x=copy.deepcopy(inp);x['prime']=1000005;bads.append(x)
 for x in bads:
  try:domain(x)
  except ValueError:pass
  else:raise ValueError('malformed finite domain accepted')
 return{'literal_actual_factorial_controls':literal,'negative_controls':['incorrect_actual_Gamma_shift','full_degree_common_root_corruption','three_malformed_domains']}
def root_path(explicit):
 if explicit:return explicit.resolve()
 for p in[HERE,*HERE.parents]:
  if(p/'research/MATH_STATUS.json').is_file():return p
 raise SystemExit('Use --repo-root for relocated replay')
def main():
 ap=ArgumentParser(description=__doc__);ap.add_argument('--repo-root',type=Path);ap.add_argument('--jobs',type=int,default=4);ap.add_argument('--output','--write-receipt',dest='output',type=Path);a=ap.parse_args()
 need(1<=a.jobs<=32,'bad worker count');root=root_path(a.repo_root);ledger=json.loads((HERE/LEDGER).read_text())
 required={('packet',p)for p in[PROOF,Path(__file__).name,INPUTS]}|{('repository','research/factorial/image-rigidity/'+p)for p in['LOWER_TWO_PRIME_GAP_SFC.tex','LOWER_TWO_FOURTH_POWER_GAP_SFC.tex','UNIFORM_PADE_PARABOLA_SFC.tex','lower_exponent_initial_and_boundary.tex']}
 keys={(x['location'],x['path'])for x in ledger['files']};need(required<=keys and len(keys)==len(ledger['files']),'missing indispensable source pins')
 pins=[]
 for v in ledger['files']:
  need(v['location']in['packet','repository'],'unknown source location')
  p=(HERE if v['location']=='packet'else root)/v['path'];need(p.is_file()and digest(p)==v['sha256'],'source pin mismatch:'+str(p));pins.append({'path':v['path'],'sha256':v['sha256']})
 inp=json.loads((HERE/INPUTS).read_text());rows=domain(inp);cc=controls(inp);completed=[]
 with ProcessPoolExecutor(max_workers=a.jobs)as pool:
  futures=[pool.submit(worker,row)for row in rows]
  for future in as_completed(futures):
   row=future.result();completed.append(row);print('PASS gap',row['gap'],'all',row['tested'],'actual consecutive pairs',flush=True)
 completed.sort(key=lambda x:x['gap'])
 need(len(completed)==40 and sum(x['tested']for x in completed)==176816,'complete run not finished')
 result={'status':'PASS','scope':'Complete actual lower-two ordinary consecutive finite bridge d6..45,2<=N<max(4000,128(d+2)); the full all-terminal theorem additionally requires the written analytic tail/overlap and retained small-gap proofs.',
 'lower_exponent':2,'actual_physical_shapes':['2*N-1 at degree N-1','2*N+1 at degree N'],'prime':1000003,
 'completed_cases':176816,'complete_gap_families':40,'maximum_terminal':6015,'maximum_factorial_argument':282705,'gaps':completed,'controls':cc,
 'versions':{'python':sys.version.split()[0],'python_flint':flint.__version__,'FLINT':flint.__FLINT_VERSION__},
 'source_sha256':{'proof':digest(HERE/PROOF),'checker':digest(Path(__file__)),'inputs':digest(HERE/INPUTS),'ledger':digest(HERE/LEDGER)},'verified_source_pins':pins,'external_review_or_formal_verification':False}
 dest=a.output or HERE/RECEIPT;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print('PASS complete finite bridge:176816 actual pairs,40 gap families')
if __name__=='__main__':main()
