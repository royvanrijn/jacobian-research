#!/usr/bin/env python3
from pathlib import Path
import sympy as sp
root=Path(__file__).resolve().parent
t=sp.symbols('t'); av=sp.symbols('a2:8'); dv=sp.symbols('d2:13')
A=t+sum(av[i-2]*t**i for i in range(2,8))+t**8
D=sum(dv[i-2]*t**i for i in range(2,13))
R=sp.Poly(sp.expand(2*A*sp.diff(D,t)-3*sp.diff(A,t)*D-t**2),t)
subs={}; residual=[]
for k in range(R.degree()+1):
 e=sp.cancel(R.nth(k).subs(subs))
 if e==0:continue
 solved=False
 for v in dv:
  if e.has(v) and sp.Poly(e,v).degree()==1 and not any(e.has(w) for w in dv if w!=v):
   pv=sp.Poly(e,v);c=pv.nth(1);subs[v]=sp.cancel(-(e-c*v)/c);solved=True;break
 if not solved:residual.append(sp.factor(e))
assert (len(subs),len(residual))==(11,6)
polys=[]
for e in residual:
 num=sp.together(e).as_numer_denom()[0]
 polys.append(sp.Poly(num,*av,domain=sp.QQ).primitive()[1].as_expr())
ss=lambda q:str(q).replace('**','^')
s=root/'firstblock_Q_exact.sing'
s.write_text('LIB "modstd.lib";\nLIB "resources.lib";\nlist sc=Resources::setcores_subtree(1);\nring rd=0,('+','.join(map(str,av))+'),dp;\nideal I=\n'+',\n'.join(ss(q) for q in polys)+';\noption(redSB); timer=1;\nideal G=modStd(I);\nprint("DP_SIZE="+string(size(G)));\nring rl=0,('+','.join(map(str,av))+'),lp;\nideal G=fetch(rd,G);\nideal L=fglm(rd,G);\nprint("LEX_SIZE="+string(size(L))); L;\n')
(root/'firstblock_data.py').write_text('DSUBS='+repr({str(k):str(v) for k,v in subs.items()})+'\n')
print(s)
