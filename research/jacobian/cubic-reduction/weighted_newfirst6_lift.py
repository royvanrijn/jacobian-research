#!/usr/bin/env python3
"""Lift the generic parameter-field unit and expose every exceptional theta denominator."""
from pathlib import Path
import hashlib,json,subprocess,time
import sympy as s
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 out=HERE/'weighted-newfirst6/parameter-field';source=HERE/'weighted-newfirst6/input.json';d=json.loads(source.read_text());V=s.symbols(' '.join(d['variables']));theta=s.Symbol('theta');byname=dict(zip(d['variables'],V));W=tuple(v for v in V if str(v)!='f2');P=(theta,)+W
 eqs=[s.expand(sum(s.Rational(c)*s.prod(v**a for v,a in zip(V,e)) for e,c in terms).subs(byname['f2'],theta)) for terms in d['equations']]
 lines=['ring r=(0,theta),('+','.join(map(str,W))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs)+';','matrix M=lift(I,ideal(1));','matrix Check=matrix(I)*M;','if(Check[1,1]!=1){print("ERROR");quit;}','int i; poly q;',
 'for(i=1;i<=nrows(M);i++){print("GEN:"+string(i));q=M[i,1];while(q!=0){print("COEF:"+string(leadcoef(q)));print("EXP:"+string(leadexp(q)));q=q-lead(q);}print("END");}','quit;']
 script=out/'lift.sing';script.write_text('\n'.join(lines)+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=30);text=r.stdout+r.stderr;assert r.returncode==0 and '?' not in text and 'ERROR' not in text,text[:1000]
 except subprocess.TimeoutExpired as e:
  text=str(e.stdout or '')+'\n'+str(e.stderr or '');(out/'lift.log').write_text(text);(out/'lift-checkpoint.json').write_text(json.dumps(dict(status='timeout',source_sha256=sha(Path(__file__)),seconds=time.monotonic()-started),indent=2)+'\n');return
 (out/'lift.log').write_text(text);raw=text.splitlines();pos=0;data=[];D=s.Poly(1,theta,domain=s.QQ);cache={}
 for i in range(len(eqs)):
  assert raw[pos]=='GEN:%d'%(i+1);pos+=1;terms=[]
  while raw[pos]!='END':
   assert raw[pos].startswith('COEF:') and raw[pos+1].startswith('EXP:');string=raw[pos][5:]
   if string not in cache:
    value=s.cancel(s.sympify(string.replace('^','**'),locals={'theta':theta}));cache[string]=value;den=s.Poly(s.denom(value),theta,domain=s.QQ);D=D.lcm(den).monic()
   e=tuple(map(int,raw[pos+1][4:].split(',')));assert len(e)==len(W);terms.append((e,cache[string]));pos+=2
  pos+=1;data.append(terms)
 assert pos==len(raw);multipliers=[]
 for terms in data:
  poly={}
  for e,c in terms:
   value=s.Poly(s.cancel(D.as_expr()*c),theta,domain=s.QQ)
   for (power,),co in value.terms():
    if co:poly[(power,)+e]=poly.get((power,)+e,0)+co
  multipliers.append([[list(e),str(c)] for e,c in sorted(poly.items(),reverse=True) if c])
 target=[[[power]+[0]*len(W),str(c)] for (power,),c in D.terms() if c];generators=[[[list(e),str(c)] for e,c in s.Poly(f,*P).terms() if c] for f in eqs]
 factors=[dict(polynomial=str(f.as_expr()),multiplicity=k) for f,k in s.factor_list(D)[1]]
 payload=dict(status='exact-primary-parameter-certificate',source_sha256=sha(Path(__file__)),input_sha256=sha(source),singular_source_sha256=sha(script),variables=list(map(str,P)),generators=generators,multipliers=multipliers,target=target,denominator=str(D.as_expr()),factors=factors,seconds=round(time.monotonic()-started,3),scope='Cleared generic unit identity D(theta)=sum M_i I_i. Only roots of D can support a coefficient solution; independent reconstruction and coefficient replay required before an exclusion is accepted.')
 (out/'certificate.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print(json.dumps({k:v for k,v in payload.items() if k not in ['generators','multipliers','target']}),flush=True)
if __name__=='__main__':main()
