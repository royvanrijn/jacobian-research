#!/usr/bin/env python3
"""Export exact membership of complete135-parameter necessary relations."""
from pathlib import Path
import hashlib,json,subprocess,time
import sympy as s
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 folder=HERE/'carrier-reparam-full';ip=folder/'input.json';d=json.loads(ip.read_text());names=d['variables'];V=s.symbols(' '.join(names));by=dict(zip(names,V));eqs=[sum(s.Rational(c)*s.prod(v**k for v,k in zip(V,e)) for e,c in terms) for terms in d['high_equations']];active={1,2,3,4,5,6,7,8,9,10,11,12,13,17,18,19,22,26,27,28,30,34,35,38,39,42,45,79,80,83,96,117,118};texts=['z'+str(j) for j in range(1,136) if j not in active]+['z3-1','z1-z2','z7-z9','z8-z9','z10-z12','z2*z9-1','z5*z12-z1','z4-z9*z12','z11-z4-3*z22','z5*z42+z22','z2*z17-2*z13','z2*z18-z5*z22','z5*z38-3*z5*z96-z17','z12*z30-z96','z12*z34-z117','z12*z35-z118','z2*z80-z27','z9*z28-z118','z9*z26-z79','z12*z79-z12*z117+3*z28','z9*z19+3*z34','z6*z9-z9-z19','z5*z39+z96','z5*z45-z80','z83-z118'];targets=[s.sympify(t,locals=by) for t in texts];assert len(targets)==127
 script=folder/'classification-lift.sing';script.write_text('\n'.join(['ring r=0,('+','.join(names)+'),dp;','short=0;','ideal I='+','.join(map(sing,eqs))+';','ideal T='+','.join(map(sing,targets))+';','matrix M=lift(I,T);','matrix Check=matrix(I)*M-matrix(T);','int j;int i;poly q;','for(j=1;j<=ncols(Check);j++){if(Check[1,j]!=0){print("ERROR");quit;}}','for(j=1;j<=ncols(M);j++){for(i=1;i<=nrows(M);i++){q=M[i,j];if(q!=0){print("ENTRY:"+string(j)+","+string(i));while(q!=0){print("COEF:"+string(leadcoef(q)));print("EXP:"+string(leadexp(q)));q=q-lead(q);}print("END");}}}','quit;'])+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=30);text=r.stdout+r.stderr;assert r.returncode==0 and '?' not in text and 'ERROR' not in text,text[:1000]
 except subprocess.TimeoutExpired as e:
  text=(e.stdout or b'').decode() if isinstance(e.stdout,bytes) else str(e.stdout or '');(folder/'classification-lift.log').write_text(text);(folder/'classification-lift-checkpoint.json').write_text(json.dumps(dict(status='timeout',source_sha256=sha(Path(__file__)),seconds=round(time.monotonic()-started,3)),indent=2)+'\n');print('TIMEOUT');return
 (folder/'classification-lift.log').write_text(text);raw=text.splitlines();pos=0;entries=[]
 while pos<len(raw):
  assert raw[pos].startswith('ENTRY:'),raw[pos];j,i=map(int,raw[pos][6:].split(','));pos+=1;terms=[]
  while raw[pos]!='END':
   assert raw[pos].startswith('COEF:') and raw[pos+1].startswith('EXP:');c=str(s.Rational(raw[pos][5:]));e=list(map(int,raw[pos+1][4:].split(',')));assert len(e)==135;terms.append([e,c]);pos+=2
  pos+=1;entries.append(dict(target=j-1,generator=i-1,multiplier=terms))
 enc=lambda f:[[list(e),str(c)] for e,c in s.Poly(f,*V).terms() if c];payload=dict(status='exact-primary-memberships',source_sha256=sha(Path(__file__)),input_sha256=sha(ip),singular_source_sha256=sha(script),variables=names,target_texts=texts,targets=[enc(f) for f in targets],entries=entries,seconds=round(time.monotonic()-started,3));(folder/'classification-certificate.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print('TARGETS',len(targets),'ENTRIES',len(entries),'TERMS',sum(len(e['multiplier']) for e in entries),'SECONDS',payload['seconds'])
if __name__=='__main__':main()
