#!/usr/bin/env python3
"""Bounded complete native coefficient-power identities for unresolved mixed rays."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import json,hashlib,time,subprocess
from verify import variable,add,mul,derivative,determinant3
HERE=Path(__file__).resolve().parent;OUT=HERE/'weight-resonance7';CASES=OUT/'mixed-cases';CASES.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def enc(p):return [{'exponent':list(e),'coefficient':str(c)} for e,c in sorted(p.items(),reverse=True)]
def render(p,names):
 return '+'.join('*'.join(['('+str(c)+')']+[z+('^'+str(k) if k!=1 else '') for z,k in zip(names,e) if k]) for e,c in sorted(p.items(),reverse=True)) or '0'
def main():
 start=time.monotonic();census=OUT/'census.json';data=json.loads(census.read_text());records=[r for r in data['records'] if r['gate']=='unresolved' and r['mixed_positive_negative_magnitudes'] and r['mixed_positive_negative_magnitudes']!=[1,1,1]];cp=OUT/'mixed-lift-checkpoint.json';old=json.loads(cp.read_text()) if cp.exists() else None;assert old is None or (old['source_sha256']==sha(Path(__file__)) and old['census_sha256']==sha(census));results=old['cases'] if old else [];done={tuple(r['magnitudes']) for r in results if r['status'] in ['complete-coefficient-powers','partial-coefficient-powers','timeout']}
 def checkpoint(status):cp.write_text(json.dumps({'status':status,'source_sha256':sha(Path(__file__)),'census_sha256':sha(census),'case_count':len(records),'cases':results,'this_invocation_seconds':time.monotonic()-start},indent=2)+'\n')
 for rec in records:
  p,q,r=rec['mixed_positive_negative_magnitudes'];key=(p,q,r)
  if key in done:continue
  if time.monotonic()-start>55:checkpoint('invocation-time-cap');print('TIME CAP',len(done),'of',len(records),flush=True);return
  name='p%dq%dr%d'%key;weight=[p,-q,-r];basis=[(i,e) for i in range(3) for e in product(range(8),repeat=3) if 2<=sum(e)<=7 and sum(e[j]*weight[j] for j in range(3))==weight[i]];assert len(basis)==rec['nonlinear_slots'];N=len(basis);n=3+N;F=[variable(i,n) for i in range(3)]
  for j,(i,e) in enumerate(basis):F[i]=add(F[i],mul({e+(0,)*N:Q(1)},variable(3+j,n)))
  J=determinant3([[derivative(z,i) for i in range(3)] for z in F],n);J=add(J,{(0,)*n:Q(1)},Q(-1));eq=[{e[3:]:c for e,c in J.items() if e[:3]==mon} for mon in sorted({e[:3] for e in J},reverse=True)];first=[j for j,(i,e) in enumerate(basis) if i==0];names=['c'+str(j) for j in range(N)];inp={'status':'complete-native-model','source_sha256':sha(Path(__file__)),'magnitudes':list(key),'weights':weight,'canonical_census_weight':rec['weight'],'basis':[{'output':i,'exponent':list(e)} for i,e in basis],'parameter_names':names,'generators':[enc(z) for z in eq],'positive_output_parameters':first};ip=CASES/(name+'.json');ip.write_text(json.dumps(inp,indent=2)+'\n');script=CASES/(name+'.sing');lines=['ring rr=0,('+','.join(names)+'),dp;','short=0;','ideal I='+','.join(render(z,names) for z in eq)+';','ideal G=std(I);','int i;int k;poly target;matrix M;poly pp;']
  for j in first:
   lines+=['k=1;target='+names[j]+';','while(k<=32 && reduce(target,G)!=0){k=2*k;target='+names[j]+'^k;}','if(k>32){print("MISS:'+str(j)+'");}else{','print("TARGET:'+str(j)+':"+string(k));','M=lift(I,ideal(target));','if((matrix(I)*M)[1,1]!=target){print("ERROR");quit;}','for(i=1;i<=nrows(M);i++){print("GEN:"+string(i));pp=M[i,1];while(pp!=0){print("COEF:"+string(leadcoef(pp)));print("EXP:"+string(leadexp(pp)));pp=pp-lead(pp);}print("END");}','}']
  lines+=['quit;'];script.write_text('\n'.join(lines)+'\n');entry={'magnitudes':list(key),'name':name,'input_sha256':sha(ip),'script_sha256':sha(script),'status':'running'};results.append(entry);checkpoint('running');tick=time.monotonic()
  try:z=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=5);output=z.stdout+z.stderr
  except subprocess.TimeoutExpired as e:
   output=e.stdout or b'';output=output.decode() if isinstance(output,bytes) else output;(CASES/(name+'.log')).write_text(output);entry.update(status='timeout',seconds=time.monotonic()-tick);done.add(key);checkpoint('running');print('TIMEOUT',name,flush=True);continue
  (CASES/(name+'.log')).write_text(output);entry['seconds']=time.monotonic()-tick
  if z.returncode!=0 or '?' in output or 'ERROR' in output:entry['status']='producer-error';checkpoint('producer-error');raise RuntimeError(name+':'+output[:2000])
  raw=output.splitlines();pos=0;targets=[];miss=[]
  for j in first:
   if raw[pos]=='MISS:'+str(j):miss.append(j);pos+=1;continue
   prefix='TARGET:'+str(j)+':';assert raw[pos].startswith(prefix);power=int(raw[pos][len(prefix):]);pos+=1;multipliers=[]
   for i in range(len(eq)):
    assert raw[pos]=='GEN:'+str(i+1);pos+=1;row=[]
    while raw[pos]!='END':
     assert raw[pos].startswith('COEF:') and raw[pos+1].startswith('EXP:');c=Q(raw[pos][5:]);e=list(map(int,raw[pos+1][4:].split(',')));assert len(e)==N;row.append({'exponent':e,'coefficient':str(c)});pos+=2
    pos+=1;multipliers.append(row)
   targets.append({'parameter':j,'power':power,'multipliers':multipliers})
  assert pos==len(raw);cert={'status':'exact-primary-coefficient-powers','source_sha256':sha(Path(__file__)),'input_sha256':sha(ip),'script_sha256':sha(script),'targets':targets,'unresolved_parameters':miss};cp0=CASES/(name+'-certificate.json');cp0.write_text(json.dumps(cert,indent=2)+'\n');entry.update(status='partial-coefficient-powers' if miss else 'complete-coefficient-powers',certificate_sha256=sha(cp0),targets=len(targets),unresolved_parameters=miss,multiplier_terms=sum(len(z) for t in targets for z in t['multipliers']),max_power=max([0]+[t['power'] for t in targets]));done.add(key);checkpoint('running')
  if miss or len(done)%50==0:print('PROGRESS',len(done),'/',len(records),name,entry['status'],entry['multiplier_terms'],'terms',flush=True)
 checkpoint('complete');counts={}
 for z in results:counts[z['status']]=counts.get(z['status'],0)+1
 print('COMPLETE',counts,'multiplier terms',sum(z.get('multiplier_terms',0) for z in results),flush=True)
if __name__=='__main__':main()
