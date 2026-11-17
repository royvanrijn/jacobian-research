#!/usr/bin/env python3
"""Coefficient specialization search on the established raw13-term three-variable map."""
from pathlib import Path
import hashlib,json,subprocess,time
import sympy as sp
from power_search import cubic12,X,DB
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 out=HERE/'raw13-specialization';out.mkdir(exist_ok=True);base,_,_=cubic12(2);Y=X[:3];assert all(sp.Poly(f,*Y).total_degree()<=7 for f in base)
 supports=[[(e,c) for e,c in sp.Poly(f,*Y).terms() if sum(e)>=2] for f in base];assert list(map(len,supports))==[2,5,6]
 assert sp.expand(sp.det(sp.Matrix(base).jacobian(Y)))==1
 A=sp.symbols('a1:14');d=sp.Symbol('d');F=[];cursor=0;atoms=[]
 for row,terms in enumerate(supports):
  f=Y[row]
  for e,c in terms:f+=A[cursor]*sp.prod(x**k for x,k in zip(Y,e));atoms.append(dict(row=row,exponent=list(e),coefficient=str(c)));cursor+=1
  F.append(f)
 determinant=sp.Poly(sp.expand(sp.det(sp.Matrix(F).jacobian(Y))-1),*Y);eqs=[sp.expand(c) for e,c in determinant.terms() if c]
 # The first point is (0,1,c). Its first two images are always0 and1;
 # c is recovered uniquely from the third collision equation after solving.
 point=dict(zip(Y,[-1,1,d]));eqs += [sp.expand(F[0].subs(point)),sp.expand(F[1].subs(point)-1)]
 variables=tuple(A)+(d,);logs=[];started=time.monotonic();payload=dict(source_sha256=sha(Path(__file__)),base_source_sha256=sha(DB/'verify.py'),atoms=atoms,variables=list(map(str,variables)),equations=[[[list(e),str(c)] for e,c in sp.Poly(f,*variables).terms() if c] for f in eqs],scope='Raw J(0)=I coordinates, not collision-axis normalized. Variable d permits the second point(-1,1,d), with first point(0,1,c) recovered linearly.')
 (out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
 for i,a in enumerate(A):
  lines=['ring r=0,('+','.join(map(str,variables))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs+[a])+';','ideal G=std(I);','if(size(G)==1 && G[1]==1){print("UNIT");}else{print("PROPER");print(G);}','quit;']
  script=out/('coefficient-%02d.sing'%i);script.write_text('\n'.join(lines)+'\n');t=time.monotonic()
  try:
   r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=5);text=r.stdout+r.stderr;status='unit' if text.strip()=='UNIT' else ('proper' if text.startswith('PROPER') and '?' not in text else 'error')
  except subprocess.TimeoutExpired as exc:text=str(exc.stdout or '')+'\n'+str(exc.stderr or '');status='timeout'
  (out/('coefficient-%02d.log'%i)).write_text(text);row=dict(coefficient=i,status=status,seconds=round(time.monotonic()-t,3));logs.append(row);print(json.dumps(row),flush=True)
 receipt=dict(status='proper-systems-require-extraction' if any(e['status']=='proper' for e in logs) else ('primary-all-unit' if all(e['status']=='unit' for e in logs) else 'inconclusive'),source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),systems=logs,seconds=round(time.monotonic()-started,3),scope='Coefficient deletions in the exact raw13-atom pattern with the specified two-point collision family. Unit statuses are primary only. No general12-term exclusion, normalized sparsity change or new counterexample is claimed.')
 (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
