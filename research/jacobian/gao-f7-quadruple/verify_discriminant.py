"""Exact source, discriminant, primitive remainder and norm-factor replay.
No symbolic discriminant/resultant, factorizer or Groebner search is called.
The separate norm-grid receipt must be fresh for the inputs checked here.
"""
from pathlib import Path
import hashlib,json,math,time,sympy as S
from fractions import Fraction
from flint import fmpz_mpoly_ctx,fmpz_poly,fmpz_mat,nmod_poly
HERE=Path(__file__).resolve().parent;START=time.monotonic();checks=0
if not __debug__:raise RuntimeError('Run without -O')
def yes(c):
 global checks
 assert c;checks+=1
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def emit(stage,**kw):print(json.dumps({'stage':stage,'seconds':time.monotonic()-START,**kw}),flush=True)
def load(name):return json.loads((HERE/name).read_text())
def bareiss(M):
 W=[r[:] for r in M];prev=1;sign=1
 for k in range(len(W)-1):
  j=next((i for i in range(k,len(W)) if W[i][k]),None)
  if j is None:return 0
  if j!=k:W[k],W[j]=W[j],W[k];sign=-sign
  v=W[k][k]
  for i in range(k+1,len(W)):
   u=W[i][k]
   for j in range(k+1,len(W)):
    q,r=divmod(W[i][j]*v-u*W[k][j],prev);assert not r;W[i][j]=q
   W[i][k]=0
  prev=v
 return sign*W[-1][-1]
def sylv(f,g):
 m,n=len(f)-1,len(g)-1;M=[[0]*(m+n) for _ in range(m+n)]
 for i in range(n):
  for j,c in enumerate(reversed(f)):M[i][i+j]=c
 for i in range(m):
  for j,c in enumerate(reversed(g)):M[n+i][i+j]=c
 return M
# Fresh source relation on the already certified quadruple incidence.
src=load('quadruple-input.json');old=json.loads((HERE.parent/'gao-f7-linear-cube/high-root-bivariate-input.json').read_text());can=load('generic-triple-input.json');base=load('receipt.json')
yes(base['result']=='PASS' and base['checker_sha256']==sha(HERE/'verify.py') and base['source_sha256']==sha(HERE/'quadruple-input.json') and base['generic_curve_source_sha256']==sha(HERE.parent/'gao-f7-linear-cube/high-root-bivariate-input.json'))
d,a,s,b=S.symbols('d a s b');ctx=fmpz_mpoly_ctx.get(['d','a','s','b']);dn,an,sn,bn=ctx.gens()
def poly(e):return ctx.from_dict({ex:int(c) for ex,c in S.Poly(e,d,a,s,b).terms()})
def coeff(F,k,j):return ctx.from_dict({ex[:k]+(0,)+ex[k+1:]:int(c) for ex,c in F.to_dict().items() if ex[k]==j})
f2=poly(S.sympify(src['equations'][0]['polynomial']));L=S.Poly(S.sympify(old['linear_equations'][0]),b);P1,Q1=poly(L.nth(1)),poly(L.nth(0));A=coeff(f2,3,2);D=sn*(3*sn**3-2)
Q=sum((poly(S.sympify(r['polynomial']))*(dn+sn)**r['coefficient'] for r in src['quotient_coefficients']),ctx.constant(0));LG=A*Q-coeff(Q,3,2)*f2;yes(LG.degrees()[3]<=1);CR=P1*coeff(LG,3,0)-Q1*coeff(LG,3,1);C0=poly(S.sympify(can['equations'][1]['polynomial']));yes(CR==8*D*C0);yes(2*coeff(C0,0,8)==D*P1);yes(poly(S.sympify(can['equations'][0]['polynomial']))==poly(S.sympify(old['equations'][3]['polynomial'])))
yes(list(map(int,C0.degrees()))==[8,7,27,0])
# Entire bivariate discriminant identity by a complete INTEGER grid.
data=load('residual-discriminant-candidate.json');yes(data['source_sha256']==sha(HERE/'generic-triple-input.json'));delta=ctx.from_dict({(ex[0],ex[1],ex[2],0):int(c) for ex,c in data['terms']});yes(all(ex[0]==ex[3]==0 for ex in delta.to_dict()))
Crows=[{} for _ in range(9)]
for (i,j,k,z),c in C0.to_dict().items():Crows[int(i)].setdefault(int(j),{})[int(k)]=int(c)
Drows={}
for (i,j,k,z),c in delta.to_dict().items():Drows.setdefault(int(j),{})[int(k)]=int(c)
def rowpolys(row):return [fmpz_poly([row.get(j,{}).get(k,0) for k in range(max(row.get(j,{0:0}))+1)]) for j in range(max(row)+1)]
Crows=list(map(rowpolys,Crows));Drows=rowpolys(Drows);yes(delta.degrees()[1]<=105-3 and delta.degrees()[2]<=405-18)
# Res_d(C0,C0') has degree <=15*7 in a, <=15*27 in s.
# LC_d(C0)*delta has the same verified upper bounds.
yes((coeff(C0,0,8)*delta).degrees()[1]<=105 and (coeff(C0,0,8)*delta).degrees()[2]<=405)
for si,sv in enumerate(range(-203,203)):
 ca=[fmpz_poly([int(p(sv)) for p in row]) for row in Crows];da=fmpz_poly([int(p(sv)) for p in Drows])
 for av in range(-53,53):
  f=[int(p(av)) for p in ca];g=[(j+1)*f[j+1] for j in range(8)];M=sylv(f,g);yes(len(M)==15);expect=f[-1]*int(da(av));yes(int(fmpz_mat(M).det())==expect);yes(bareiss(M)==expect)
 if (si+1)%64==0:emit('discriminant_grid',s_values=si+1)
emit('discriminant_verified',points=106*406)
# Literal fraction-free identities, with all deleted contents retained.
rc=load('primitive-remainder-candidate.json');yes(rc['input_sha256']==sha(HERE/'residual-discriminant-candidate.json') and rc['curve_sha256']==sha(HERE/'generic-triple-input.json'));G=poly(S.sympify(can['equations'][0]['polynomial']));ell=coeff(G,1,8);R=delta
allowed=[fmpz_poly([0,1]),fmpz_poly([-2,0,0,3]),fmpz_poly([-1,0,0,3]),fmpz_poly([1,0,0,3])]
for j,row in enumerate(rc['removed_contents']):
 cp=fmpz_poly(list(map(int,row['content'])));C=ctx.from_dict({(0,0,k,0):int(c) for k,c in enumerate(cp.coeffs()) if c});yes(bool(C))
 for v in allowed:
  while True:
   u,r=divmod(cp,v)
   if r:break
   cp=u
 yes(cp.degree()==0 and cp[0]!=0)
 if j:
  deg=int(R.degrees()[1]);yes(deg==row['cancelled_a_degree'] and deg>=8);R=ell*R-coeff(R,1,deg)*an**(deg-8)*G
 R,rem=divmod(R,C);yes(not rem)
yes(len(rc['removed_contents'])==88 and R.degrees()[1]==7);stored=ctx.from_dict({(0,ex[0],ex[1],0):int(c) for ex,c in rc['terms']});yes(R==stored)
# The separately replayed norm is bound to precisely this R and G.
nr=load('norm-grid-receipt.json');yes(nr['result']=='PASS' and nr['checker_sha256']==sha(HERE/'replay_norm.py') and nr['c_checker_sha256']==sha(HERE/'verify_norm_modular.c') and nr['remainder_sha256']==sha(HERE/'primitive-remainder-candidate.json') and nr['curve_input_sha256']==sha(HERE/'generic-triple-input.json') and nr['norm_sha256']==sha(HERE/'crt-norm-candidate.json'))
nc=load('crt-norm-candidate.json');fac=load('crt-norm-factors-candidate.json');yes(fac['norm_sha256']==sha(HERE/'crt-norm-candidate.json'));N=fmpz_poly(list(map(int,nc['coefficients'])));H=fmpz_poly(list(map(int,fac['residual_coefficients'])));B=fmpz_poly(list(map(int,fac['coordinate_pole_factor']['coefficients'])));prod=fmpz_poly([int(fac['content'])])*H*B**fac['coordinate_pole_factor']['power']
yes([v['power'] for v in fac['pole_factors']]==[1645,805,805,448] and fac['coordinate_pole_factor']['power']==28)
for row in fac['pole_factors']:prod*=fmpz_poly(row['coefficients'])**row['power']
yes(prod==N and N.degree()==11167 and H.degree()==912 and B.degree()==87)
# Literal inverses in the actual rank912 finite-free presentation.
prime=1000000007;yes(all(prime%j for j in range(2,math.isqrt(prime)+1)));hp=nmod_poly([int(c)%prime for c in H.coeffs()],prime);yes(hp.degree()==H.degree() and hp[hp.degree()]!=0);hm=[int(c)*pow(int(hp[hp.degree()]),-1,prime)%prime for c in hp.coeffs()]
def trim(v):
 while v and not v[-1]:v.pop()
 return v
def mulrem(f,g):
 out=[0]*max(0,len(f)+len(g)-1)
 for i,c in enumerate(f):
  for j,d in enumerate(g):out[i+j]=(out[i+j]+c*d)%prime
 out=trim(out)
 while len(out)>=len(hm):
  k=len(out)-len(hm);v=out[-1]
  for j,c in enumerate(hm):out[k+j]=(out[k+j]-v*c)%prime
  out=trim(out)
 return out
units=[]
for name,q in [('squarefree',H.derivative()),('coordinate_pole',B)]+[('leading_pole_'+str(j),f) for j,f in enumerate(allowed)]:
 qp=nmod_poly([int(c)%prime for c in q.coeffs()],prime);g,u,v=qp.xgcd(hp);yes(g.degree()==0 and int(g[0])!=0);u=u*pow(int(g[0]),-1,prime)%hp;yes(qp*u%hp==1);yes(mulrem(list(map(int,qp.coeffs())),list(map(int,u.coeffs())))==[1]);units.append(name)
receipt={'result':'PASS','scope':'Exact C0=(D*P1/2)*Q8 on quadruple incidence; complete integer discriminant grid;87 literal primitive reductions; coupled fully replayed integer norm; exact factorization with squarefree degree912 residual away from all presentation poles. Chart closure is a separate component.','checks':checks,'discriminant_grid_points':106*406,'integer_determinant_engines':['FLINT','plain Python Bareiss'],'remainder_steps':87,'norm_degree':11167,'residual_degree':912,'coordinate_pole_degree':87,'coordinate_pole_power':28,'finite_free_units':units,'checker_sha256':sha(Path(__file__)),'hashes':{f:sha(HERE/f) for f in ['quadruple-input.json','generic-triple-input.json','residual-discriminant-candidate.json','primitive-remainder-candidate.json','crt-norm-candidate.json','crt-norm-factors-candidate.json','norm-grid-receipt.json']},'seconds':time.monotonic()-START,'assurance':'Exact arithmetic and written determinant/discriminant/norm arguments; no independent end-to-end, formal or external review.'}
(HERE/'discriminant-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
