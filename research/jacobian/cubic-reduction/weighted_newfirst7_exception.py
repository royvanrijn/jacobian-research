#!/usr/bin/env python3
"""Exact rational exceptional parameter replay for the parabolic boundary line."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,subprocess,time
import sympy as s
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 out=HERE/'weighted-newfirst7/parameter-field/parabolic-exception';out.mkdir(exist_ok=True);source=out.parent/'parabolic-boundary/input.json';d=json.loads(source.read_text());V=s.symbols(' '.join(d['variables']),seq=True);theta=s.Symbol('f4');eqs=[s.expand(sum(s.Rational(c)*theta**e[0]*s.prod(v**a for v,a in zip(V,e[1:])) for e,c in terms).subs(theta,s.Rational(7,4))) for terms in d['equations']];eqs=[f for f in eqs if f];ip=out/'input.json';ip.write_text(json.dumps(dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),variables=list(map(str,V)),specialization={'f4':'7/4'},generators=[[[list(e),str(c)] for e,c in s.Poly(f,*V).terms() if c] for f in eqs]),indent=2,sort_keys=True)+'\n')
 lines=['ring r=0,('+','.join(map(str,V))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs)+';','matrix M=lift(I,ideal(1));','matrix Check=matrix(I)*M;','if(Check[1,1]!=1){print("ERROR");quit;}','int i; poly q;','for(i=1;i<=nrows(M);i++){print("GEN:"+string(i));q=M[i,1];while(q!=0){print("COEF:"+string(leadcoef(q)));print("EXP:"+string(leadexp(q)));q=q-lead(q);}print("END");}','quit;'];script=out/'lift.sing';script.write_text('\n'.join(lines)+'\n');started=time.monotonic()
 try:r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr
 except subprocess.TimeoutExpired as e:
  (out/'checkpoint.json').write_text(json.dumps(dict(status='timeout',seconds=time.monotonic()-started),indent=2)+'\n');return
 (out/'lift.log').write_text(text);assert r.returncode==0 and '?' not in text and 'ERROR' not in text,text[:1000];raw=text.splitlines();pos=0;multipliers=[]
 for i in range(len(eqs)):
  assert raw[pos]=='GEN:%d'%(i+1);pos+=1;terms=[]
  while raw[pos]!='END':
   assert raw[pos].startswith('COEF:') and raw[pos+1].startswith('EXP:');c=Q(raw[pos][5:]);e=list(map(int,raw[pos+1][4:].split(',')));assert len(e)==len(V);terms.append([e,str(c)]);pos+=2
  pos+=1;multipliers.append(terms)
 assert pos==len(raw);payload=dict(status='exact-primary-unit-identity',source_sha256=sha(Path(__file__)),input_sha256=sha(ip),singular_source_sha256=sha(script),multipliers=multipliers,seconds=round(time.monotonic()-started,3));(out/'certificate.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print(json.dumps({k:v for k,v in payload.items() if k!='multipliers'}),flush=True)
if __name__=='__main__':main()
