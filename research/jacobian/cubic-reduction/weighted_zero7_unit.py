#!/usr/bin/env python3
"""Structured exact unit lifts for reduced full-ring unmarked branches."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,subprocess,time,sys
import sympy as sp
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 mode=sys.argv[1];out=HERE/'weighted-zero7'/mode;source=out/'input.json';d=json.loads(source.read_text());V=sp.symbols(' '.join(d['variables']),seq=True);eqs=[sp.Add(*(sp.Rational(c)*sp.prod(v**a for v,a in zip(V,e)) for e,c in terms)) for terms in d['equations']];script=out/'lift.sing';script.write_text('\n'.join(['ring r=0,('+','.join(map(str,V))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs)+';','matrix M=lift(I,ideal(1));','matrix Check=matrix(I)*M;','if(Check[1,1]!=1){print("ERROR");quit;}','int i; poly q;','for(i=1;i<=nrows(M);i++){print("GEN:"+string(i));q=M[i,1];while(q!=0){print("COEF:"+string(leadcoef(q)));print("EXP:"+string(leadexp(q)));q=q-lead(q);}print("END");}','quit;'])+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=30);text=r.stdout+r.stderr
 except subprocess.TimeoutExpired as e:
  text=(e.stdout or b'').decode() if isinstance(e.stdout,bytes) else str(e.stdout or '');(out/'lift.log').write_text(text);(out/'lift-checkpoint.json').write_text(json.dumps(dict(status='timeout',source_sha256=sha(Path(__file__)),seconds=time.monotonic()-started),indent=2)+'\n');print('TIMEOUT',mode,flush=True);return
 (out/'lift.log').write_text(text);assert r.returncode==0 and '?' not in text and 'ERROR' not in text,text[:1000];raw=text.splitlines();pos=0;multipliers=[]
 for i in range(len(eqs)):
  assert raw[pos]=='GEN:%d'%(i+1);pos+=1;terms=[]
  while raw[pos]!='END':
   assert raw[pos].startswith('COEF:') and raw[pos+1].startswith('EXP:');c=Q(raw[pos][5:]);e=list(map(int,raw[pos+1][4:].split(',')));assert len(e)==len(V);terms.append([e,str(c)]);pos+=2
  pos+=1;multipliers.append(terms)
 assert pos==len(raw);payload=dict(status='exact-primary-unit-identity',source_sha256=sha(Path(__file__)),input_sha256=sha(source),singular_source_sha256=sha(script),multipliers=multipliers,seconds=round(time.monotonic()-started,3));(out/'certificate.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print(mode,'TERMS',sum(map(len,multipliers)),'SECONDS',payload['seconds'],flush=True)
if __name__=='__main__':main()
