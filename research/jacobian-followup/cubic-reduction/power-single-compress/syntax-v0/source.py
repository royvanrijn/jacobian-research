#!/usr/bin/env python3
"""Search arbitrary one-cube compressions of the saved89-cube basis by exact ideals."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations,combinations_with_replacement
from collections import defaultdict
import hashlib,json,subprocess,time
import sympy as sp
from flint import fmpq_mat
from power_line_compress import cubic
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def cliques(forms):
 adj={i:set() for i in range(18)}
 for f in forms:
  for i,j in combinations([a for a,c in enumerate(f) if c],2):adj[i].add(j);adj[j].add(i)
 out=[]
 def bk(R,P,X):
  if not P and not X:out.append(sorted(R));return
  u=max(P|X,key=lambda u:len(P&adj[u]))
  for v in sorted(P-adj[u]):bk(R|{v},P&adj[v],X&adj[v]);P.remove(v);X.add(v)
 bk(set(),set(adj),set());return sorted(out)
def encode(poly,variables):return [[list(e),str(c)] for e,c in sp.Poly(poly,*variables).terms() if c]
def sing(poly):return str(poly).replace("**","^")
def main():
 out=HERE/"power-single-compress";out.mkdir(exist_ok=True);source=HERE/"power-basis-anneal/candidate.json";d=json.loads(source.read_text())
 forms=d["forms"];n=len(forms);assert n==89;B=[[Q(c) for c in row] for row in d["B"]]
 cubes=[cubic(f) for f in forms];mons=sorted(set().union(*(set(p) for p in cubes)));mi={e:i for i,e in enumerate(mons)}
 joint=fmpq_mat([[str(p.get(e,0)) for p in cubes]+[str(int(i==j)) for j in range(len(mons))] for i,e in enumerate(mons)])
 rr,rank=joint.rref();assert rank==len(mons);assert all(rr[i,i]==1 for i in range(n))
 groups=defaultdict(list)
 for j in range(n):
  v=tuple(row[j] for row in B);c=next(c for c in v if c);groups[tuple(a/c for a in v)].append(j)
 pairs=[(j,k) for key,ids in groups.items() for j,k in combinations(ids,2)]
 logs=[];started=time.monotonic()
 for case,S in enumerate(cliques(forms)):
  variables=sp.symbols("u1:"+str(len(S)+1));coords=[sp.Integer(0)]*len(mons);outside=[]
  for ids in combinations_with_replacement(range(len(S)),3):
   e=[0]*18
   for i in ids:e[S[i]]+=1
   mult=1 if ids[0]==ids[2] else (3 if ids[0]==ids[1] or ids[1]==ids[2] else 6)
   term=mult*sp.prod(variables[i] for i in ids);e=tuple(e)
   if e not in mi:outside.append(term)
   else:
    for i in range(len(mons)):
     c=rr[i,n+mi[e]]
     if c:coords[i]+=sp.Rational(str(c))*term
  coords=[sp.expand(p) for p in coords];constraints=outside+[p for p in coords[n:] if p]
  # Remove duplicate homogeneous equations up to nonzero scalar.
  unique={}
  for f in constraints:
   f=sp.Poly(f,*variables).monic().as_expr();unique[str(f)]=f
  constraints=list(unique.values());selected=[]
  for j,k in pairs:
   if not coords[j] or not coords[k]:continue
   v=next(o for o in range(18) if B[o][j]);ratio=B[v][k]/B[v][j]
   selected.append((j,k,ratio))
  data=dict(case=case,support=S,variables=[str(v) for v in variables],constraints=[encode(f,variables) for f in constraints],
    coordinates=[encode(f,variables) for f in coords[:n]],pairs=[[j,k,str(r)] for j,k,r in selected])
  jp=out/("case-%02d.json"%case);jp.write_text(json.dumps(data,indent=2)+"\n")
  lines=['ring r=0,('+','.join(map(str,variables))+'),dp;','option(redSB);','ideal I='+(','.join(sing(f) for f in constraints) or '0')+';','ideal J; ideal G;']
  for j,k,ratio in selected:
   lines+=['J=I,'+sing(coords[j]-1)+','+sing(coords[k]-sp.Rational(str(ratio)))+';',
    'G=std(J);','if(size(G)==1 && G[1]==1){print("PASS:%d:%d");}else{print("OPEN:%d:%d");print(G);}'%(j,k,j,k)]
  lines+=['quit;'];script=out/("case-%02d.sing"%case);script.write_text('\n'.join(lines)+'\n')
  try:
   process=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=15)
   text=process.stdout+process.stderr;status='closed' if process.returncode==0 and text.count('PASS:')==len(selected) and 'OPEN:' not in text and '?' not in text else 'open-or-error'
  except subprocess.TimeoutExpired as exc:text=(exc.stdout or b'').decode() if isinstance(exc.stdout,bytes) else (exc.stdout or '');status='timeout'
  (out/("case-%02d.log"%case)).write_text(text)
  row=dict(case=case,support=S,constraints=len(constraints),pairs=len(selected),status=status,seconds=round(time.monotonic()-started,3));logs.append(row);print(json.dumps(row),flush=True)
  receipt=dict(status='running',source_sha256=sha(Path(__file__)),input_sha256=sha(source),cases=logs)
  (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
 receipt['status']='primary-all-closed' if all(r['status']=='closed' for r in logs) else 'inconclusive-or-candidate'
 receipt['scope']='All coordinate-support cliques for a single new cube in the saved cube span; proper ideals require extraction, unit ideals require independent certificates. No global optimum claim.'
 (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
