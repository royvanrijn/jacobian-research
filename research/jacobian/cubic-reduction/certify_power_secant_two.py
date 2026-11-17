#!/usr/bin/env python3
"""Export exact rational unit identities for all607 unrestricted two-new-cube systems."""
from pathlib import Path
import hashlib,json,subprocess,time
import sympy as sp
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def poly(terms,V):return sp.Add(*(sp.Rational(c)*sp.prod(x**a for x,a in zip(V,e)) for e,c in terms))
def main():
 out=HERE/'power-secant-two';primary=json.loads((out/'checkpoint.json').read_text());assert primary['status']=='primary-all-unit' and primary['source_sha256']==sha(HERE/'power_secant_two.py')
 logs=[];started=time.monotonic()
 for entry in primary['cases']:
  if entry['status']=='coordinate-zero-exclusion':continue
  i,j=entry['pair'];tag='pair-%02d-%02d'%(i,j);source=out/(tag+'.json');d=json.loads(source.read_text());V=sp.symbols(' '.join(d['variables']));constraints=[poly(f,V) for f in d['membership']];coords=[poly(f,V) for f in d['coordinates']]
  lines=['ring r=0,('+','.join(d['variables'])+'),dp;','short=0;','ideal J; matrix M; matrix P; int h;']
  for triple in d['triples']:
   ids=[a for a,c in triple];generators=constraints+[coords[a]-sp.Rational(c) for a,c in triple]
   lines+=['J='+','.join(sing(f) for f in generators)+';','M=lift(J,ideal(1));','P=matrix(J)*M;','if(P[1,1]!=1){print("ERROR");quit;}',
     'print("CERT:%d:%d:%d:%d");'%(*ids,len(generators)),'for(h=1;h<=nrows(M);h++){print(string(M[h,1]));}']
  lines+=['quit;'];script=out/(tag+'-cert.sing');script.write_text('\n'.join(lines)+'\n')
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=15);text=r.stdout+r.stderr;(out/(tag+'-cert.log')).write_text(text)
  assert r.returncode==0 and '?' not in text and 'ERROR' not in text,text[:1000]
  raw=text.splitlines();position=0;certs=[]
  for triple in d['triples']:
   ids=[a for a,c in triple];count=len(constraints)+3;assert raw[position]=='CERT:%d:%d:%d:%d'%(*ids,count);position+=1;multipliers=[]
   for line in raw[position:position+count]:
    f=sp.Poly(sp.sympify(line.replace('^','**'),locals={str(v):v for v in V}),*V);multipliers.append([[list(e),str(c)] for e,c in f.terms() if c])
   position+=count;certs.append(dict(indices=ids,multipliers=multipliers))
  assert position==len(raw);payload=dict(pair=[i,j],source_sha256=sha(Path(__file__)),input_sha256=sha(source),singular_sha256=sha(script),certificates=certs)
  cp=out/(tag+'-cert.json');cp.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');logs.append(dict(pair=[i,j],identities=len(certs),certificate_sha256=sha(cp)))
  if len(logs)%10==0:print('CERTIFIED',len(logs),'PAIRS',sum(x['identities'] for x in logs),'IDENTITIES',flush=True)
 receipt=dict(status='exact-primary-identities',source_sha256=sha(Path(__file__)),primary_checkpoint_sha256=sha(out/'checkpoint.json'),cases=logs,identities=sum(x['identities'] for x in logs),seconds=round(time.monotonic()-started,3),scope='All607 identities need independent coefficient replay and full support/membership coverage before a theorem is accepted.')
 (out/'certification.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k!='cases'}),flush=True)
if __name__=='__main__':main()
