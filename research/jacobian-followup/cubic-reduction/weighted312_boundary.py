#!/usr/bin/env python3
"""Changed mu1 queries using the complete b=0 determinant boundary split."""
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
import json,hashlib,time,subprocess
HERE=Path(__file__).resolve().parent;OUT=HERE/'weighted-312-412';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def enc(p,V):return [{'exponent':list(e),'coefficient':str(c)} for e,c in s.Poly(s.expand(p),*V).terms() if c]
def render(p,V):return '+'.join('*'.join(['('+str(c)+')']+[str(v)+('^'+str(k) if k!=1 else '') for v,k in zip(V,e) if k]) for e,c in s.Poly(p,*V).terms() if c) or '0'
def main():
 parent=OUT/'p3.json';d=json.loads(parent.read_text());V=list(s.symbols(' '.join(d['parameter_names']),seq=True));original=[sum(s.Rational(t['coefficient'])*s.prod(x**k for x,k in zip(V,t['exponent'])) for t in row) for row in d['native_generators']];records=[];start=time.monotonic();cp=OUT/'boundary-checkpoint.json'
 def checkpoint(status):cp.write_text(json.dumps({'status':status,'source_sha256':sha(Path(__file__)),'parent_model_sha256':sha(parent),'normalization_cover':'L=mu=1. b=0 coefficients give T*A*W=0 and A*(8W+3T)=3T*(3T-1) after eliminating H. Complete field cover: T0/W0; T0/A0; W0/A=3T-1; A0/T=1/3. No unchanged mu1 ideal rerun.','cases':records,'seconds':time.monotonic()-start},indent=2)+'\n')
 cases=[('mu1-T0W0',{4:s.Integer(1),0:s.Integer(1),2:s.Integer(0),13:s.Integer(0)}),('mu1-T0A0',{4:s.Integer(1),0:s.Integer(1),2:s.Integer(0),7:s.Integer(0)}),('mu1-W0A3Tminus1',{4:s.Integer(1),0:s.Integer(1),13:s.Integer(0),7:3*V[2]-1}),('mu1-A0Tthird',{4:s.Integer(1),0:s.Integer(1),7:s.Integer(0),2:s.Rational(1,3)})]
 for name,fix in cases:
  fixed={V[j]:c for j,c in fix.items()};eqs=[s.expand(z.subs(fixed)) for z in original];sub={x:fixed.get(x,x) for x in V};active=[x for x in V if x not in fixed];steps=[];priority=[V[j] for j in [14,11,10,7,8,6,5,12,13,9]]
  while True:
   found=None
   for var in priority:
    if var not in active:continue
    for j,eq in enumerate(eqs):
     co=s.diff(eq,var)
     if co.is_Rational and co!=0:found=(var,j,co,s.expand(var-eq/co));break
    if found:break
   if not found:break
   var,j,pivot,expr=found;steps.append({'parameter':V.index(var),'native_equation_index':j,'constant_pivot':str(pivot),'expression':enc(expr,V)});eqs=[s.expand(z.subs(var,expr)) for z in eqs];sub={x:s.expand(z.subs(var,expr)) for x,z in sub.items()};active.remove(var)
  indices=[j for j,z in enumerate(eqs) if z!=0];gens=[eqs[j] for j in indices];model={'status':'complete-primary-normalized-reduction','source_sha256':sha(Path(__file__)),'parent_model_sha256':sha(parent),'case':name,'fixed_parameters':{str(j):enc(c,V) for j,c in fix.items()},'pivot_steps':steps,'free_parameter_names':list(map(str,active)),'substitutions':[enc(sub[x],active) for x in V],'native_generator_indices':indices,'generators':[enc(z,active) for z in gens]};ip=OUT/(name+'.json');ip.write_text(json.dumps(model,indent=2)+'\n');script=OUT/(name+'.sing');script.write_text('\n'.join(['ring rr=0,('+','.join(map(str,active))+'),dp;','short=0;','ideal I='+','.join(render(z,active) for z in gens)+';','ideal G=std(I);','if(size(G)!=1||G[1]!=1){print("NONUNIT");print(G);quit;}','matrix M=lift(I,ideal(1));','if((matrix(I)*M)[1,1]!=1){print("ERROR");quit;}','int i;poly pp;','for(i=1;i<=nrows(M);i++){print("GEN:"+string(i));pp=M[i,1];while(pp!=0){print("COEF:"+string(leadcoef(pp)));print("EXP:"+string(leadexp(pp)));pp=pp-lead(pp);}print("END");}','quit;'])+'\n');rec={'case':name,'input_sha256':sha(ip),'script_sha256':sha(script),'fixed':{str(j):str(c) for j,c in fix.items()},'free_variables':len(active),'constant_pivots':len(steps),'generators':len(gens),'status':'running'};records.append(rec);checkpoint('running');tick=time.monotonic();print('MODEL',name,len(active),'variables',len(steps),'constant pivots',len(gens),'generators',flush=True)
  try:z=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=10);output=z.stdout+z.stderr
  except subprocess.TimeoutExpired as e:
   output=e.stdout or b'';output=output.decode() if isinstance(output,bytes) else output;(OUT/(name+'.log')).write_text(output);rec.update(status='changed-reduced-timeout10',seconds=time.monotonic()-tick);checkpoint('running');print('TIMEOUT',name,flush=True);continue
  (OUT/(name+'.log')).write_text(output);rec['seconds']=time.monotonic()-tick
  if 'NONUNIT' in output:rec['status']='primary-nonunit';checkpoint('running');print('NONUNIT',name,output[:300],flush=True);continue
  assert z.returncode==0 and '?' not in output and 'ERROR' not in output,output[:2000];lines=output.splitlines();pos=0;mult=[]
  for i in range(len(gens)):
   assert lines[pos]=='GEN:'+str(i+1);pos+=1;row=[]
   while lines[pos]!='END':
    assert lines[pos].startswith('COEF:') and lines[pos+1].startswith('EXP:');co=Q(lines[pos][5:]);e=list(map(int,lines[pos+1][4:].split(',')));assert len(e)==len(active);row.append({'exponent':e,'coefficient':str(co)});pos+=2
   pos+=1;mult.append(row)
  assert pos==len(lines);cert={'status':'exact-primary-unit','source_sha256':sha(Path(__file__)),'input_sha256':sha(ip),'script_sha256':sha(script),'multipliers':mult};op=OUT/(name+'-unit.json');op.write_text(json.dumps(cert,indent=2)+'\n');rec.update(status='exact-primary-unit',certificate_sha256=sha(op),multiplier_terms=sum(map(len,mult)),max_multiplier_degree=max((sum(t['exponent']) for row in mult for t in row),default=0));checkpoint('running');print('UNIT',name,rec['multiplier_terms'],'terms',rec['seconds'],'seconds',flush=True)
 checkpoint('complete');print('COMPLETE',[(r['case'],r['status']) for r in records],flush=True)
if __name__=='__main__':main()
