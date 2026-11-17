#!/usr/bin/env python3
"""Export compact structured unit certificates for weighted coefficient branches."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,subprocess,time
import sympy as s
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def expression(terms,V):return s.Add(*(s.Rational(c)*s.prod(v**a for v,a in zip(V,e)) for e,c in terms))
def main():
 out=HERE/'weighted-unit-certificates';out.mkdir(exist_ok=True);cases=[('small-h',HERE/'weighted-tail7/small-h/input.json'),('high-w',HERE/'weighted-tail7/high-w/generic/input.json'),('mid-w',HERE/'weighted-tail7/mid-w/input.json'),('theta-zero',HERE/'weighted-newfirst6/input.json')];logs=[];started=time.monotonic()
 for tag,source in cases:
  d=json.loads(source.read_text());V=s.symbols(' '.join(d['variables']));eqs=[expression(f,V) for f in d['equations']]
  if tag=='theta-zero':
   theta=next(v for v in V if str(v)=='f2');eqs=[s.expand(f.subs(theta,0)) for f in eqs];V=tuple(v for v in V if v!=theta);eqs=[f for f in eqs if f]
  inputdata=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),variables=list(map(str,V)),generators=[[[list(e),str(c)] for e,c in s.Poly(f,*V).terms() if c] for f in eqs],specialization={'f2':'0'} if tag=='theta-zero' else {})
  ip=out/(tag+'-input.json');ip.write_text(json.dumps(inputdata,indent=2,sort_keys=True)+'\n')
  lines=['ring r=0,('+','.join(map(str,V))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs)+';','matrix M=lift(I,ideal(1));','matrix Check=matrix(I)*M;','if(Check[1,1]!=1){print("ERROR");quit;}','int i; poly q;',
   'for(i=1;i<=nrows(M);i++){print("GEN:"+string(i));q=M[i,1];while(q!=0){print("COEF:"+string(leadcoef(q)));print("EXP:"+string(leadexp(q)));q=q-lead(q);}print("END");}','quit;']
  script=out/(tag+'.sing');script.write_text('\n'.join(lines)+'\n');r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr;(out/(tag+'.log')).write_text(text);assert r.returncode==0 and '?' not in text and 'ERROR' not in text,text[:1000]
  raw=text.splitlines();pos=0;multipliers=[]
  for i in range(len(eqs)):
   assert raw[pos]=='GEN:%d'%(i+1);pos+=1;terms=[]
   while raw[pos]!='END':
    assert raw[pos].startswith('COEF:') and raw[pos+1].startswith('EXP:');c=Q(raw[pos][5:]);e=list(map(int,raw[pos+1][4:].split(',')));assert len(e)==len(V);terms.append([e,str(c)]);pos+=2
   pos+=1;multipliers.append(terms)
  assert pos==len(raw);payload=dict(tag=tag,source_sha256=sha(Path(__file__)),input_sha256=sha(ip),singular_source_sha256=sha(script),multipliers=multipliers)
  cp=out/(tag+'.json');cp.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');logs.append(dict(tag=tag,variables=len(V),generators=len(eqs),input_sha256=sha(ip),certificate_sha256=sha(cp)));print('CERTIFIED',tag,flush=True)
 receipt=dict(status='exact-primary-unit-identities',source_sha256=sha(Path(__file__)),cases=logs,seconds=round(time.monotonic()-started,3),scope='Independent universal-equation reconstruction, branch coverage proofs and polynomial replay remain required.')
 (out/'manifest.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
