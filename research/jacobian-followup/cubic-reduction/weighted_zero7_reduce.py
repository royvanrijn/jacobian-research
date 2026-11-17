#!/usr/bin/env python3
"""Reduce necessary boundary equations before full-ring exceptional-stratum tests."""
from pathlib import Path
import hashlib,json,subprocess,time,sys
import sympy as sp
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(f,V):return [[list(e),str(c)] for e,c in sp.Poly(f,*V).terms() if c]
def dec(terms,V):return sp.Add(*(sp.Rational(c)*sp.prod(v**a for v,a in zip(V,e)) for e,c in terms))
def main():
 case,mode=sys.argv[1:];assert mode in ['purecubic','parabola'];base='cubicgraph' if mode=='purecubic' else mode;source=HERE/'weighted-zero7'/case/base/'input.json';d=json.loads(source.read_text());P=sp.symbols(' '.join(d['parameters']),seq=True) if d['parameters'] else ();W=sp.symbols(' '.join(d['variables']),seq=True);V=P+W;sub={}
 if mode=='purecubic':sub[P[1]]=0
 eqs=[sp.expand(dec(f,V).subs(sub)) for f in d['equations']];eqs=[f for f in eqs if f];linear=[sp.expand(dec(f,V).subs(sub)) for f in d['boundary_equations']];linear=[f for f in linear if f];P=tuple(p for p in P if p not in sub);M,b=sp.linear_eq_to_matrix(linear,W);red,pivots=M.row_join(b).rref();assert all(i<len(W) for i in pivots);free=[i for i in range(len(W)) if i not in pivots];Z=tuple(W[i] for i in free);solve={W[p]:sp.factor(red[k,len(W)]-sum(red[k,j]*W[j] for j in free)) for k,p in enumerate(pivots)};D=sp.Poly(1,*P,domain=sp.QQ) if P else None
 for value in solve.values():
  den=sp.denom(sp.cancel(value))
  if P:D=D.lcm(sp.Poly(den,*P,domain=sp.QQ)).monic()
  else:assert not den.free_symbols
 D=D.as_expr() if P else sp.Integer(1);reduced=[]
 for f in eqs:
  num,den=sp.fraction(sp.cancel(f.subs(solve,simultaneous=True)))
  if num:reduced.append(sp.expand(num))
 nonzero={'cubicgraph':sp.Symbol('f1'),'purecubic':sp.Symbol('f1'),'quadraticgraph':sp.Symbol('f2'),'parabola':sp.Symbol('f4'),'parallel':sp.Symbol('f4')*(4*sp.Symbol('f4')-1),'linear':sp.Integer(1)}[mode];sat=sp.Symbol('branch_inv');saturation=sp.expand(D*nonzero);ZV=Z+P+((sat,) if saturation!=1 else ());all_eqs=reduced+([1-sat*saturation] if saturation!=1 else []);out=HERE/'weighted-zero7'/case/'reduced'/mode;out.mkdir(parents=True,exist_ok=True);payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),mode=mode,variables=list(map(str,ZV)),linear_rank=len(pivots),free_tail_variables=list(map(str,Z)),parameters=list(map(str,P)),specialization={str(k):str(v) for k,v in sub.items()},linear_solution={str(k):str(v) for k,v in solve.items()},linear_denominator=str(D),required_nonzero=str(nonzero),saturation=str(saturation),equations=[enc(f,ZV) for f in all_eqs],scope='Full polynomial-ring test on the stated nonzero branch and away from recorded linear-solution denominator roots. Every excluded denominator stratum must be checked separately; no generic inference.');(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print('REDUCED',mode,len(pivots),len(Z),'D',sp.factor(D),flush=True);script=out/'model.sing';script.write_text('\n'.join(['ring r=0,('+','.join(map(str,ZV))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in all_eqs)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr;status='primary-full-ring-basis' if r.returncode==0 and text.startswith('DIM:') and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=(e.stdout or b'').decode() if isinstance(e.stdout,bytes) else str(e.stdout or '');status='timeout'
 (out/'model.log').write_text(text);receipt=dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),mode=mode,seconds=round(time.monotonic()-started,3));(out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True);print(text[:8000],flush=True)
if __name__=='__main__':main()
