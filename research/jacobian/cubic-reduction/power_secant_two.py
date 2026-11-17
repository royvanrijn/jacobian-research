#!/usr/bin/env python3
"""Exact unrestricted two-new-cube replacement systems across all source-clique pairs."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations,combinations_with_replacement
from collections import defaultdict
import hashlib,json,subprocess,time
import sympy as sp
from flint import fmpq_mat
from power_single_compress import sing
from power_line_compress import cubic
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def encode(f,V):return [[list(e),str(c)] for e,c in sp.Poly(f,*V).terms() if c]
def main():
 out=HERE/"power-secant-two";out.mkdir(exist_ok=True);source=HERE/"power-basis-anneal/candidate.json";d=json.loads(source.read_text());B=[[Q(c) for c in row] for row in d['B']];forms=d['forms'];n=89
 cubes=[cubic(f) for f in forms];mons=sorted(set().union(*(set(f) for f in cubes)));mi={e:i for i,e in enumerate(mons)}
 R,_=fmpq_mat([[str(f.get(e,0)) for f in cubes]+[str(int(i==j)) for j in range(len(mons))] for i,e in enumerate(mons)]).rref();assert all(R[i,i]==1 for i in range(n))
 groups=defaultdict(list)
 for j in range(n):
  v=tuple(row[j] for row in B);c=next(x for x in v if x);groups[tuple(x/c for x in v)].append((j,c))
 triples=[t for ids in groups.values() for t in combinations(ids,3)];assert len(triples)==214
 supports=[json.loads((HERE/'power-single-compress'/('case-%02d.json'%i)).read_text())['support'] for i in range(15)]
 cached=[]
 for case,S in enumerate(supports):
  V=sp.symbols('q1:'+str(len(S)+1));coords=[sp.Integer(0)]*89;residual={}
  for ids in combinations_with_replacement(range(len(S)),3):
   e=[0]*18
   for a in ids:e[S[a]]+=1
   mult=1 if ids[0]==ids[2] else (3 if ids[0]==ids[1] or ids[1]==ids[2] else 6);term=mult*sp.prod(V[a] for a in ids);e=tuple(e)
   if e not in mi:residual[('o',e)]=term
   else:
    for i in range(len(mons)):
     value=R[i,n+mi[e]]
     if not value:continue
     if i<89:coords[i]+=sp.Rational(str(value))*term
     else:residual[('i',i)]=residual.get(('i',i),0)+sp.Rational(str(value))*term
  cached.append((V,[sp.expand(f) for f in coords],{e:sp.expand(f) for e,f in residual.items()}))
 logs=[];started=time.monotonic();systems=0
 for i,j in combinations_with_replacement(range(15),2):
  S,T=supports[i],supports[j];U=sp.symbols('u1:'+str(len(S)+1));V=sp.symbols('v1:'+str(len(T)+1));variables=U+V
  vi,ci,ri=cached[i];vj,cj,rj=cached[j];si=dict(zip(vi,U));sj=dict(zip(vj,V))
  coords=[sp.expand(f.xreplace(si)+g.xreplace(sj)) for f,g in zip(ci,cj)]
  selected=[t for t in triples if all(coords[a] for a,c in t)]
  if not selected:
   logs.append(dict(pair=[i,j],status='coordinate-zero-exclusion',systems=0));continue
  equations=[];seen=set()
  for label in sorted(set(ri)|set(rj)):
   f=sp.expand(ri.get(label,sp.Integer(0)).xreplace(si)+rj.get(label,sp.Integer(0)).xreplace(sj))
   if not f:continue
   f=sp.Poly(f,*variables).monic().as_expr();key=str(f)
   if key not in seen:seen.add(key);equations.append(f)
  tag='pair-%02d-%02d'%(i,j);inputdata=dict(pair=[i,j],supports=[S,T],variables=list(map(str,variables)),
    membership=[encode(f,variables) for f in equations],coordinates=[encode(f,variables) for f in coords],
    triples=[[[a,str(c)] for a,c in triple] for triple in selected],source_sha256=sha(Path(__file__)),input_sha256=sha(source))
  jp=out/(tag+'.json');jp.write_text(json.dumps(inputdata,indent=2,sort_keys=True)+'\n')
  lines=['ring r=0,('+','.join(map(str,variables))+'),dp;','short=0;','ideal I='+(','.join(sing(f) for f in equations) or '0')+';','ideal J; ideal G;']
  for triple in selected:
   ids=[a for a,c in triple];lines+=['J=I,'+','.join(sing(coords[a]-sp.Rational(str(c))) for a,c in triple)+';',
      'G=std(J);','if(size(G)==1 && G[1]==1){print("UNIT:%d:%d:%d");}else{print("PROPER:%d:%d:%d");print(G);}'%(*ids,*ids)]
  lines+=['quit;'];script=out/(tag+'.sing');script.write_text('\n'.join(lines)+'\n');t=time.monotonic()
  try:
   r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=3);text=r.stdout+r.stderr
   status='unit' if r.returncode==0 and text.count('UNIT:')==len(selected) and '?' not in text and 'PROPER:' not in text else ('proper' if 'PROPER:' in text and '?' not in text else 'error')
  except subprocess.TimeoutExpired as exc:text=str(exc.stdout or '')+'\n'+str(exc.stderr or '');status='timeout'
  (out/(tag+'.log')).write_text(text);systems+=len(selected);row=dict(pair=[i,j],status=status,systems=len(selected),membership_equations=len(equations),variables=len(variables),seconds=round(time.monotonic()-t,3),input_sha256=sha(jp));logs.append(row);print(json.dumps(row),flush=True)
  receipt=dict(status='running',source_sha256=sha(Path(__file__)),input_sha256=sha(source),cases=logs,systems=systems,timeout_per_case=3)
  (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
 receipt.update(status='primary-all-unit' if all(row['status'] in ['unit','coordinate-zero-exclusion'] for row in logs) else 'proper-or-inconclusive',
   total_support_pairs=len(logs),seconds=round(time.monotonic()-started,3),scope='Unrestricted two-new-form systems using the clique support lemma and fixed-span optimality. Unit statuses require independent membership reconstruction and polynomial certificate replay before a theorem is accepted.')
 (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k!='cases'}),flush=True)
if __name__=='__main__':main()
