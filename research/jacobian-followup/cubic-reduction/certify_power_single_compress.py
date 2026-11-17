#!/usr/bin/env python3
"""Export exact unit-ideal identities for the arbitrary single-cube obstruction."""
from pathlib import Path
import hashlib,json,subprocess,time
import sympy as sp
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def poly(terms,variables):return sp.Add(*(sp.Rational(c)*sp.prod(v**a for v,a in zip(variables,e)) for e,c in terms))
def main():
 out=HERE/"power-single-compress";primary=json.loads((out/"checkpoint.json").read_text());assert primary["status"]=="primary-all-closed"
 assert primary["source_sha256"]==sha(HERE/"power_single_compress.py")
 logs=[];started=time.monotonic()
 for entry in primary["cases"]:
  case=entry["case"];source=out/("case-%02d.json"%case);d=json.loads(source.read_text());variables=sp.symbols(" ".join(d["variables"]))
  if not isinstance(variables,tuple):variables=(variables,)
  constraints=[poly(f,variables) for f in d["constraints"]];coords=[poly(f,variables) for f in d["coordinates"]]
  lines=['ring r=0,('+','.join(d["variables"])+'),dp;','short=0;','ideal J; matrix M; matrix P; int h;']
  for j,k,ratio in d["pairs"]:
   generators=constraints+[coords[j]-1,coords[k]-sp.Rational(ratio)]
   lines+=['J='+','.join(sing(f) for f in generators)+';',
    'M=lift(J,ideal(1));','P=matrix(J)*M;','if(P[1,1]!=1){print("ERROR");quit;}',
    'print("CERT:%d:%d:%d");'%(j,k,len(generators)),
    'for(h=1;h<=nrows(M);h++){print(string(M[h,1]));}']
  lines+=['quit;'];script=out/("cert-%02d.sing"%case);script.write_text('\n'.join(lines)+'\n')
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=30);text=r.stdout+r.stderr
  (out/("cert-%02d.log"%case)).write_text(text);assert r.returncode==0 and '?' not in text and 'ERROR' not in text,text[:1500]
  raw=text.splitlines();position=0;certificates=[]
  for j,k,ratio in d["pairs"]:
   count=len(constraints)+2;assert raw[position]=='CERT:%d:%d:%d'%(j,k,count),(case,raw[position]);position+=1
   multipliers=[]
   for line in raw[position:position+count]:
    f=sp.Poly(sp.sympify(line.replace('^','**'),locals={str(v):v for v in variables}),*variables)
    multipliers.append([[list(e),str(c)] for e,c in f.terms() if c])
   position+=count;certificates.append(dict(pair=[j,k],ratio=ratio,multipliers=multipliers))
  assert position==len(raw)
  payload=dict(case=case,source_sha256=sha(Path(__file__)),case_sha256=sha(source),singular_source_sha256=sha(script),certificates=certificates)
  cp=out/("cert-%02d.json"%case);cp.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
  logs.append(dict(case=case,identities=len(certificates),certificate_sha256=sha(cp)));print('CERTIFIED',case,len(certificates),flush=True)
 receipt=dict(status='exact-primary-identities',source_sha256=sha(Path(__file__)),primary_checkpoint_sha256=sha(out/'checkpoint.json'),
   cases=logs,identities=sum(e['identities'] for e in logs),seconds=round(time.monotonic()-started,3),scope='Unit-ideal certificates require independent reconstruction and polynomial replay before a theorem is registered.')
 (out/'certification.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
