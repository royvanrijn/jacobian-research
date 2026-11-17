#!/usr/bin/env python3
"""Bounded exact primary units for the complete odd nonlinear-first strata."""
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
import json,hashlib,subprocess,time
HERE=Path(__file__).resolve().parent;OUT=HERE/'mixed-weight2'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def enc(p,V):return [{'exponent':list(e),'coefficient':str(c)} for e,c in s.Poly(s.expand(p),*V).terms() if c]
def render(p,V):
 return '+'.join('*'.join(['('+str(c)+')']+[str(v)+('^'+str(k) if k!=1 else '') for v,k in zip(V,e) if k]) for e,c in s.Poly(p,*V).terms() if c) or '0'
def main():
 records=[];started=time.monotonic();cp=OUT/'odd-lift-checkpoint.json'
 def checkpoint(status):cp.write_text(json.dumps({'status':status,'source_sha256':sha(Path(__file__)),'cases':records,'elapsed_seconds':time.monotonic()-started},indent=2)+'\n')
 for r in [5,7]:
  src=OUT/('r'+str(r)+'.json');d=json.loads(src.read_text());full=s.symbols(' '.join(d['parameter_names']),seq=True);target=[(r+3)//2,1,1];nu_index=d['basis'][0].index(target);nu=full[nu_index];V=[x for x in full if x!=nu];gens=[s.expand(sum(s.Rational(t['coefficient'])*s.prod(v**e for v,e in zip(full,t['exponent'])) for t in row).subs(nu,1)) for row in d['native_generators']];gens=[z for z in gens if z!=0];model={'status':'complete-primary-nu1-stratum','source_sha256':sha(Path(__file__)),'full_model_sha256':sha(src),'r':r,'normalized_parameter':str(nu),'normalization':'nonzero coefficient nu of first-output a^((r+3)/2)*b*c is made1 by diagonal source/output conjugacy in c','parameter_names':list(map(str,V)),'generators':[enc(z,V) for z in gens]};mp=OUT/('r'+str(r)+'-nu1.json');mp.write_text(json.dumps(model,indent=2)+'\n');script=OUT/('r'+str(r)+'-nu1.sing');script.write_text('\n'.join(['ring rr=0,('+','.join(map(str,V))+'),dp;','short=0;','ideal I='+','.join(render(z,V) for z in gens)+';','ideal G=std(I);','if(size(G)!=1||G[1]!=1){print("NONUNIT");print(G);quit;}','matrix M=lift(I,ideal(1));','matrix CC=matrix(I)*M;','if(CC[1,1]!=1){print("ERROR");quit;}','int i;poly pp;','for(i=1;i<=nrows(M);i++){print("GEN:"+string(i));pp=M[i,1];while(pp!=0){print("COEF:"+string(leadcoef(pp)));print("EXP:"+string(leadexp(pp)));pp=pp-lead(pp);}print("END");}','quit;'])+'\n');rec={'r':r,'input':mp.name,'input_sha256':sha(mp),'script_sha256':sha(script),'status':'running'};records.append(rec);checkpoint('running');tick=time.monotonic()
  try:z=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);output=z.stdout+z.stderr
  except subprocess.TimeoutExpired as e:
   output=e.stdout or b'';output=output.decode() if isinstance(output,bytes) else output;(OUT/('r'+str(r)+'-nu1.log')).write_text(output);rec.update(status='timeout',seconds=time.monotonic()-tick);checkpoint('running');print('TIMEOUT',r,flush=True);continue
  (OUT/('r'+str(r)+'-nu1.log')).write_text(output);rec['seconds']=time.monotonic()-tick
  if 'NONUNIT' in output:rec['status']='primary-nonunit';checkpoint('running');print('NONUNIT',r,output[:200],flush=True);continue
  assert z.returncode==0 and '?' not in output and 'ERROR' not in output,output[:1000];lines=output.splitlines();pos=0;multipliers=[]
  for i in range(len(gens)):
   assert lines[pos]=='GEN:%d'%(i+1);pos+=1;row=[]
   while lines[pos]!='END':
    assert lines[pos].startswith('COEF:') and lines[pos+1].startswith('EXP:');c=Q(lines[pos][5:]);e=list(map(int,lines[pos+1][4:].split(',')));assert len(e)==len(V);row.append({'exponent':e,'coefficient':str(c)});pos+=2
   pos+=1;multipliers.append(row)
  assert pos==len(lines);cert={'status':'exact-primary-unit','source_sha256':sha(Path(__file__)),'input_sha256':sha(mp),'singular_source_sha256':sha(script),'multipliers':multipliers};certp=OUT/('r'+str(r)+'-nu1-unit.json');certp.write_text(json.dumps(cert,indent=2)+'\n');rec.update(status='exact-primary-unit',certificate_sha256=sha(certp),multiplier_terms=sum(map(len,multipliers)),max_multiplier_degree=max((sum(t['exponent']) for row in multipliers for t in row),default=0));checkpoint('running');print('UNIT',r,rec['multiplier_terms'],rec['max_multiplier_degree'],rec['seconds'],flush=True)
 checkpoint('complete')
if __name__=='__main__':main()
