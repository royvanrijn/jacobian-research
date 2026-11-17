#!/usr/bin/env python3
"""Exact exploratory components of linear forms whose cubes lie in the saved span."""
from pathlib import Path
import hashlib,json,subprocess,time
from power_single_compress import sing
import sympy as sp
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 out=HERE/"power-span-components";out.mkdir(exist_ok=True);logs=[];started=time.monotonic()
 for case in range(15):
  cp=HERE/"power-single-compress"/("case-%02d.json"%case);d=json.loads(cp.read_text());variables=sp.symbols(" ".join(d["variables"]));variables=variables if isinstance(variables,tuple) else (variables,)
  polys=[sp.Add(*(sp.Rational(c)*sp.prod(x**a for x,a in zip(variables,e)) for e,c in terms)) for terms in d["constraints"]]
  lines=['LIB "primdec.lib";','ring r=0,('+','.join(d["variables"])+'),dp;','short=0;',
   'ideal I='+(','.join(sing(f) for f in polys) or '0')+';','list L=minAssGTZ(I);','int a; int b; ideal P;',
   'print("COUNT:"+string(size(L)));',
   'for(a=1;a<=size(L);a++){P=L[a];print("PRIME:"+string(a)+":"+string(dim(std(P)))+":"+string(size(P)));for(b=1;b<=size(P);b++){print(string(P[b]));}}','quit;']
  spath=out/("case-%02d.sing"%case);spath.write_text('\n'.join(lines)+'\n')
  try:
   r=subprocess.run(['Singular','-q',str(spath)],capture_output=True,text=True,timeout=15);text=r.stdout+r.stderr
   status='primary-components' if r.returncode==0 and '?' not in text and 'COUNT:' in text else 'error'
  except subprocess.TimeoutExpired as e:text=str(e.stdout or '')+'\n'+str(e.stderr or '');status='timeout'
  (out/("case-%02d.log"%case)).write_text(text)
  components=[]
  if status=='primary-components':
   raw=text.splitlines();assert raw[0].startswith('COUNT:');total=int(raw[0].split(':')[1]);pos=1
   for k in range(total):
    tag,num,dimension,count=raw[pos].split(':');assert tag=='PRIME';pos+=1;count=int(count)
    equations=[]
    for line in raw[pos:pos+count]:
     f=sp.Poly(sp.sympify(line.replace('^','**'),locals={str(v):v for v in variables}),*variables)
     equations.append([[list(e),str(c)] for e,c in f.terms() if c])
    pos+=count;components.append(dict(dimension=int(dimension),equations=equations))
   assert pos==len(raw)
   payload=dict(case=case,source_sha256=sha(Path(__file__)),case_sha256=sha(cp),components=components,
      scope='Primary minimal primes over Q; no independent radical decomposition certificate yet.')
   (out/("case-%02d.json"%case)).write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
  row=dict(case=case,status=status,component_count=len(components),dimensions=[c['dimension'] for c in components],
    max_equation_degree=max((sum(e) for c in components for f in c['equations'] for e,a in f),default=0),seconds=round(time.monotonic()-started,3));logs.append(row);print(json.dumps(row),flush=True)
  receipt=dict(status='running',source_sha256=sha(Path(__file__)),cases=logs)
  (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
 receipt['status']='primary-components' if all(r['status']=='primary-components' for r in logs) else 'inconclusive'
 receipt['scope']='Discovery decomposition of the cubed-form membership ideals; no minimality or smaller Keller witness claimed.'
 (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
