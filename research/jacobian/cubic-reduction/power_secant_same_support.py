#!/usr/bin/env python3
"""Seek two new cubes outside the old span, initially on common source cliques."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
from collections import defaultdict
import hashlib,json,subprocess,time
import sympy as sp
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def expression(terms,variables):return sp.Add(*(sp.Rational(c)*sp.prod(v**e for v,e in zip(variables,p)) for p,c in terms))
def main():
 out=HERE/"power-secant-same-support";out.mkdir(exist_ok=True);source=HERE/"power-basis-anneal/candidate.json";d=json.loads(source.read_text());B=[[Q(c) for c in row] for row in d["B"]];groups=defaultdict(list)
 for j in range(89):
  v=tuple(row[j] for row in B);c=next(a for a in v if a);groups[tuple(a/c for a in v)].append((j,c))
 triples=[triple for group in groups.values() for triple in combinations(group,3)];logs=[];started=time.monotonic()
 for case in range(15):
  casepath=HERE/"power-single-compress"/("case-%02d.json"%case);data=json.loads(casepath.read_text());m=len(data["support"])
  U=sp.symbols("u1:"+str(m+1));V=sp.symbols("v1:"+str(m+1));vars=U+V
  constraints=[expression(f,U)+expression(f,V) for f in data["constraints"]]
  coords=[expression(f,U)+expression(f,V) for f in data["coordinates"]]
  for triple in triples:
   ids=[j for j,c in triple]
   if any(not data['coordinates'][j] for j in ids):continue
   tag='case-%02d-%02d-%02d-%02d'%(case,*ids);eqs=constraints+[coords[j]-sp.Rational(str(c)) for j,c in triple]
   payload=dict(case=case,support=data['support'],indices=ids,scalars=[str(c) for j,c in triple],variables=list(map(str,vars)),
     equations=[[[list(e),str(c)] for e,c in sp.Poly(f,*vars).terms() if c] for f in eqs],
     case_sha256=sha(casepath),input_sha256=sha(source),source_sha256=sha(Path(__file__)))
   jp=out/(tag+'.json');jp.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
   lines=['ring r=0,('+','.join(map(str,vars))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs)+';',
     'ideal G=std(I);','if(size(G)==1 && G[1]==1){print("UNIT");}else{print("PROPER");print(G);}','quit;']
   script=out/(tag+'.sing');script.write_text('\n'.join(lines)+'\n');t=time.monotonic()
   try:
    r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=3);text=r.stdout+r.stderr
    status='unit' if text.strip()=='UNIT' and r.returncode==0 else ('proper' if text.startswith('PROPER') and '?' not in text and r.returncode==0 else 'error')
   except subprocess.TimeoutExpired as exc:text=str(exc.stdout or '')+'\n'+str(exc.stderr or '');status='timeout'
   (out/(tag+'.log')).write_text(text);row=dict(case=case,indices=ids,status=status,equations=len(eqs),variables=len(vars),seconds=round(time.monotonic()-t,3),input_sha256=sha(jp));logs.append(row);print(json.dumps(row),flush=True)
   receipt=dict(status='running',source_sha256=sha(Path(__file__)),input_sha256=sha(source),systems=logs,timeout_per_system=3)
   (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
 receipt['status']='proper-systems-require-extraction' if any(r['status']=='proper' for r in logs) else ('primary-all-unit' if all(r['status']=='unit' for r in logs) else 'inconclusive')
 receipt['seconds']=round(time.monotonic()-started,3);receipt['scope']='Two unrestricted forms on one of the15 saved source cliques, canceling at least three proportional old columns while re-solving all retained columns. Proper ideals need point extraction and witness validation; primary unit status needs independent certificate replay. Does not cover mixed source cliques or arbitrary supports.'
 (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
