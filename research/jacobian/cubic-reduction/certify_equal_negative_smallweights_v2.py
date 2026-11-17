#!/usr/bin/env python3
"""Exact bounded Singular lifts, saved for independent rational multiplication."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,time,subprocess
import sympy as sp
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def render(p,V):
 terms=[]
 for e,c in sp.Poly(p,*V).terms():
  c=sp.Rational(c);atoms=['(%s/%s)'%(c.p,c.q)]
  atoms.extend(str(v)+('^'+str(k) if k!=1 else '') for v,k in zip(V,e) if k);terms.append('*'.join(atoms))
 return '+'.join(terms) or '0'
def main():
 folder=HERE/'equal-negative-smallweights-v2';check=json.loads((folder/'checkpoint.json').read_text());assert check['status']=='complete' and len(check['cases'])==14;start=time.monotonic();records=[]
 def checkpoint(status):(folder/'lift-checkpoint.json').write_text(json.dumps({'status':status,'source_sha256':sha(Path(__file__)),'producer_checkpoint_sha256':sha(folder/'checkpoint.json'),'cases':records,'seconds':round(time.monotonic()-start,3)},indent=2)+'\n')
 for case in check['cases']:
  if time.monotonic()-start>60:checkpoint('total-time-cap');return
  path=folder/case['model'];d=json.loads(path.read_text());V=sp.symbols(' '.join(d['parameter_names']),seq=True);eqs=[sum(sp.Rational(t['coefficient'])*sp.prod(v**k for v,k in zip(V,t['exponent'])) for t in row) for row in d['generators']];script=folder/(case['name']+'.sing')
  script.write_text('\n'.join(['ring rr=0,('+','.join(map(str,V))+'),dp;','short=0;','ideal I='+','.join(render(p,V) for p in eqs)+';','matrix M=lift(I,ideal(1));','matrix CC=matrix(I)*M;','if(CC[1,1]!=1){print("ERROR");quit;}','int i;poly pp;','for(i=1;i<=nrows(M);i++){print("GEN:"+string(i));pp=M[i,1];while(pp!=0){print("COEF:"+string(leadcoef(pp)));print("EXP:"+string(leadexp(pp)));pp=pp-lead(pp);}print("END");}','quit;'])+'\n');rec={'name':case['name'],'model_sha256':sha(path),'script_sha256':sha(script),'status':'running'};records.append(rec);checkpoint('running');tick=time.monotonic()
  try:r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=10);output=r.stdout+r.stderr
  except subprocess.TimeoutExpired as e:
   output=e.stdout or b'';output=output.decode() if isinstance(output,bytes) else output;(folder/(case['name']+'.log')).write_text(output);rec['status']='timeout';checkpoint('running');print('TIMEOUT',case['name'],flush=True);continue
  (folder/(case['name']+'.log')).write_text(output);assert r.returncode==0 and '?' not in output and 'ERROR' not in output,output[:1000];lines=output.splitlines();pos=0;multipliers=[]
  for i in range(len(eqs)):
   assert lines[pos]=='GEN:%d'%(i+1);pos+=1;row=[]
   while lines[pos]!='END':
    assert lines[pos].startswith('COEF:') and lines[pos+1].startswith('EXP:');coefficient=Q(lines[pos][5:]);exponent=list(map(int,lines[pos+1][4:].split(',')));assert len(exponent)==len(V);row.append({'exponent':exponent,'coefficient':str(coefficient)});pos+=2
   pos+=1;multipliers.append(row)
  assert pos==len(lines);payload={'status':'exact-primary-unit-identity','source_sha256':sha(Path(__file__)),'model_sha256':sha(path),'singular_source_sha256':sha(script),'multipliers':multipliers};cp=folder/(case['name']+'-unit.json');cp.write_text(json.dumps(payload,indent=2)+'\n');rec.update(status='exact-primary-unit-identity',certificate_sha256=sha(cp),multiplier_terms=sum(map(len,multipliers)),max_multiplier_degree=max((sum(t['exponent']) for row in multipliers for t in row),default=0),seconds=round(time.monotonic()-tick,3));checkpoint('running');print(case['name'],'TERMS',rec['multiplier_terms'],'MAXDEG',rec['max_multiplier_degree'],'SECONDS',rec['seconds'],flush=True)
 checkpoint('complete')
if __name__=='__main__':main()
