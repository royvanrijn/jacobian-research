#!/usr/bin/env python3
"""Complete boundary and affine-tail reductions before changed full-ring solves."""
from pathlib import Path
import hashlib,json,subprocess,time,sys
import sympy as sp
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(f,V):return [[list(e),str(c)] for e,c in sp.Poly(f,*V).terms() if c]
def main():
 case,degree=sys.argv[1:];assert case in ['both','mu-zero'] and degree in ['cubic','quad','linear'];source=HERE/'weighted-113/input.json';d=json.loads(source.read_text());V=sp.symbols(' '.join(d['variables']),seq=True);by=dict(zip(d['variables'],V));mu=-1 if case=='both' else 0;sub={by['c1']:mu,by['c3']:-1,by['c5']:0};P=by['c2'];R=by['c4'];parameters=[P,R]
 if degree!='cubic':sub[R]=0;parameters=[P]
 if degree=='linear':sub[P]=0;parameters=[]
 W=V[5:];u,w,t=sp.symbols('u w t');g,h=[sp.sympify(s,locals=by|{'u':u,'w':w}) for s in d['quotients'][1:]];B=1+mu*t+(P*t**2 if degree!='linear' else 0)+(R*t**3 if degree=='cubic' else 0);GH=[sp.expand(f.subs({u:t,w:B},simultaneous=True)) for f in [g,h]];limits=[2,6] if degree=='cubic' else [1,3];boundary=[]
 for f,limit in zip(GH,limits):boundary += [c for (k,),c in sp.Poly(f,t).terms() if k>limit and c]
 eqs=[sp.expand(sum(sp.Rational(c)*sp.prod(v**a for v,a in zip(V,e)) for e,c in terms).subs(sub,simultaneous=True)) for terms in d['equations']];eqs=[f for f in eqs if f];linear=boundary+[f for f in eqs if sp.Poly(f,*W).total_degree()<=1];M,b=sp.linear_eq_to_matrix(linear,W);rr,pivots=M.row_join(b).rref();out=HERE/'weighted-113/boundary'/(case+'-'+degree);out.mkdir(parents=True,exist_ok=True);payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),case=case,degree=degree,parameters=list(map(str,parameters)),variables=list(map(str,W)),substitution={str(k):str(v) for k,v in sub.items()},equations=[enc(f,tuple(parameters)+W) for f in eqs+boundary],boundary_equations=[enc(f,tuple(parameters)+W) for f in boundary],linear_equations=[enc(f,tuple(parameters)+W) for f in linear],scope='Boundary G Hprime-3Gprime H=1 with complete source degree bounds. Cubic graph degrees<=6,7 force degG<=2,degH<=6. Quadratic/linear graphs degrees<=4,5 force degG<=1,degH<=3. All boundary and affine-in-tail equations are retained before row reduction. Every denominator/nonzero branch remains explicit.');(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
 if len(W) in pivots:
  (out/'checkpoint.json').write_text(json.dumps(dict(status='generic-linear-inconsistent',source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),parameters=list(map(str,parameters)),rank=len(pivots)),indent=2,sort_keys=True)+'\n');print(case,degree,'GENERIC LINEAR INCONSISTENT',flush=True);return
 Z=tuple(W[i] for i in range(len(W)) if i not in pivots);solve={W[p]:sp.factor(rr[k,len(W)]-sum(rr[k,j]*W[j] for j in range(len(W)) if j not in pivots)) for k,p in enumerate(pivots)};D=sp.Poly(1,*parameters,domain=sp.QQ) if parameters else None
 for f in solve.values():
  if parameters:D=D.lcm(sp.Poly(sp.denom(sp.cancel(f)),*parameters,domain=sp.QQ)).monic()
 D=D.as_expr() if parameters else sp.Integer(1);nonzero=R if degree=='cubic' else (P if degree=='quad' else sp.Integer(1));sat=sp.expand(D*nonzero);inv=sp.Symbol('branch_inv');N=Z+tuple(parameters)+((inv,) if sat!=1 else ());reduced=[]
 for f in eqs:
  num=sp.fraction(sp.cancel(f.subs(solve,simultaneous=True)))[0]
  if num:reduced.append(sp.expand(num))
 if sat!=1:reduced.append(1-inv*sat)
 reducedout=out/'reduced';reducedout.mkdir(exist_ok=True);rd=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),variables=list(map(str,N)),parameters=list(map(str,parameters)),free_tail_variables=list(map(str,Z)),solution={str(k):str(v) for k,v in solve.items()},linear_rank=len(pivots),linear_denominator=str(D),saturation=str(sat),equations=[enc(f,N) for f in reduced],scope='Full polynomial-ring test only on the declared saturated branch. Every linear denominator root outside the intended nonzero branch needs separate treatment.');(reducedout/'input.json').write_text(json.dumps(rd,indent=2,sort_keys=True)+'\n');print(case,degree,'REDUCED',len(pivots),len(Z),'D',sp.factor(D),flush=True);script=reducedout/'model.sing';script.write_text('\n'.join(['ring r=0,('+','.join(map(str,N))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in reduced)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr;status='primary-reduced-basis' if r.returncode==0 and text.startswith('DIM:') and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=(e.stdout or b'').decode() if isinstance(e.stdout,bytes) else str(e.stdout or '');status='timeout'
 (reducedout/'model.log').write_text(text);receipt=dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(reducedout/'input.json'),seconds=round(time.monotonic()-started,3));(reducedout/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True);print(text[:2000],flush=True)
if __name__=='__main__':main()
