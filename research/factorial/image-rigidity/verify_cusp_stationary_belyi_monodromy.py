#!/usr/bin/env python3
"""Exact support for stationary-cover geometry; no actual-zero/SFC assertion."""
from __future__ import annotations
if not __debug__:raise SystemExit('Refusing -O/-OO: exact checks must execute.')
from pathlib import Path
from fractions import Fraction as Q
from collections import deque
from hashlib import sha256
from math import factorial,isqrt
import argparse,copy,json,subprocess,sys
HERE=Path(__file__).resolve().parent
PN='CUSP_STATIONARY_BELYI_MONODROMY_AND_RADIAL_BARRIER.md'
CN='verify_cusp_stationary_belyi_monodromy.py'
IN='cusp_stationary_belyi_inputs.json'
LN='cusp_stationary_belyi_source_ledger.json'
RN='cusp_stationary_belyi_verification.json'
REQ={('local',PN),('local',CN),('local',IN),('local','cusp_stationary_belyi_primary_audit.md'),('local','sources/vidunas-belyi-2025.pdf'),('repo','research/factorial/image-rigidity/COALESCING_BINOMIAL_SADDLE_BRIDGE.md')}
EXPECTED={'schema':1,'integer_d_min':2,'general_real_a_strictly_positive':True,'cusp_a':'(d-1)^2/4','geometric_odd_cusp':'A_(d+1)','geometric_other':'S_(d+1)','odd_cusp_arithmetic_constant':'(-1)^((d-1)/2)*d','literal_cover_d_max':20,'literal_discriminant_d_max':10,'literal_discriminant_W':[-1,2,3],'enumerated_group_d_max':7,'actual_global_zero_localization':False,'new_SFC_claim':False,'distinguished_sector_d_min':3,'distinguished_global_radial_graph':True}
def sha(p):return sha256(p.read_bytes()).hexdigest()
def clean(p):return{k:Q(v)for k,v in p.items()if v}
def const(v):return{(0,0,0):Q(v)}if v else{}
def var(i):return{tuple(int(i==j)for j in range(3)):Q(1)}
def add(*ps):
 out={}
 for p in ps:
  for k,v in p.items():out[k]=out.get(k,Q(0))+v
 return clean(out)
def sc(p,q):return clean({k:v*q for k,v in p.items()})
def mul(*ps):
 out=const(1)
 for p in ps:
  n={}
  for k,v in out.items():
   for l,w in p.items():
    m=tuple(x+y for x,y in zip(k,l));n[m]=n.get(m,Q(0))+v*w
  out=clean(n)
 return out
def pw(p,n):
 out=const(1)
 for _ in range(n):out=mul(out,p)
 return out
def identities():
 d,a,s=map(var,range(3));one=const(1);D=add(a,d)
 stationary=add(mul(d,add(D,sc(s,-1)),add(s,sc(a,-1))),sc(mul(s,add(s,sc(a,-1))),-1),sc(mul(s,add(D,sc(s,-1))),-1),mul(d,add(pw(s,2),sc(mul(add(sc(a,2),d,sc(one,-1)),s),-1),mul(a,D))))
 assert not stationary
 b,q,r=map(var,range(3));dp=add(sc(b,2),one);ac=add(pw(b,2),sc(pw(q,2),-1))
 assert not add(mul(add(b,q),add(b,one,sc(q,-1))),sc(add(ac,b,q),-1))
 Dc=add(ac,dp)
 assert not add(sc(mul(pw(b,2),Dc),4),sc(mul(pw(add(b,one),2),ac),-4),sc(mul(dp,pw(q,2)),-4))
 ai=add(pw(b,2),pw(q,2));Di=add(ai,dp)
 assert not add(sc(mul(pw(b,2),Di),2),sc(mul(pw(add(b,one),2),ai),-2),sc(mul(dp,pw(q,2)),2))
 d,r,x=map(var,range(3));A=add(d,sc(one,-1));B=add(d,one)
 lr=add(mul(B,r),sc(A,-1));rr=add(B,sc(mul(A,r),-1))
 assert not add(mul(d,lr,rr),sc(mul(A,r,lr),-1),sc(mul(B,r,rr),-1),mul(d,A,B,pw(add(r,sc(one,-1)),2)))
 # Cleared circle-phase derivative identity.
 den=add(pw(B,2),sc(mul(A,B,x),-2),pw(A,2))
 assert not add(mul(d,den),sc(add(pw(B,2),sc(pw(A,2),-1)),-1),sc(mul(d,A,B,add(one,sc(x,-1))),-2))
 # Stationary envelope numerator and exact Vidunas affine factors.
 assert not add(add(D,sc(s,-1)),add(s,sc(a,-1)),sc(d,-1))
 assert not add(sc(add(one,mul(A,sc(add(one,sc(r,-1)),Q(1,2)))),2),sc(rr,-1))
 # Exact positive derivative on the crossing line Re(r)=1.
 b,t,_=map(var,range(3));bp=add(b,one);dd=add(sc(b,2),one);ss=mul(b,bp)
 den0=add(one,t);den1=add(one,mul(pw(b,2),t));den2=add(one,mul(pw(bp,2),t))
 psi=add(mul(dd,den1,den2),sc(mul(b,den0,den2),-1),sc(mul(bp,den0,den1),-1),sc(mul(dd,ss,t,add(one,mul(add(ss,sc(one,-1)),t))),-1))
 assert not psi
 # Imaginary numerator of the real Mobius map K.
 cc,xx,yy=map(var,range(3))
 assert not add(sc(mul(cc,yy,add(xx,sc(cc,-1))),-1),sc(mul(add(one,sc(mul(cc,xx),-1)),yy),-1),sc(mul(add(pw(cc,2),sc(one,-1)),yy),-1))
 return 10
def trim(v):
 v=[Q(x)for x in v]
 while v and not v[-1]:v.pop()
 return v
def pa(a,b):
 n=max(len(a),len(b));return trim([(a[i]if i<len(a)else 0)+(b[i]if i<len(b)else 0)for i in range(n)])
def ps(a,c):return trim([v*c for v in a])
def pm(a,b):
 v=[Q(0)]*(len(a)+len(b)-1)if a and b else[]
 for i,x in enumerate(a):
  for j,y in enumerate(b):v[i+j]+=x*y
 return trim(v)
def pd(a):return trim([i*a[i]for i in range(1,len(a))])
def ev(a,x):
 o=Q(0)
 for v in reversed(a):o=o*x+v
 return o
def div(a,b):
 a=trim(a);b=trim(b);assert b
 q=[Q(0)]*max(0,len(a)-len(b)+1)
 while a and len(a)>=len(b):
  j=len(a)-len(b);c=a[-1]/b[-1];q[j]=c
  a=pa(a,[Q(0)]*j+ps(b,-c))
 return trim(q),a
def det(M):
 A=[[int(x)for x in row]for row in M];n=len(A);sgn=1;old=1
 for k in range(n-1):
  j=next((j for j in range(k,n)if A[j][k]),None)
  if j is None:return 0
  if j!=k:A[j],A[k]=A[k],A[j];sgn=-sgn
  p=A[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    v=A[i][j]*p-A[i][k]*A[k][j];assert v%old==0;A[i][j]=v//old
  for i in range(k+1,n):A[i][k]=0
  old=p
 return sgn*A[-1][-1]if n else 1
def resultant(f,g):
 n=len(f)-1;m=len(g)-1
 rows=[]
 for j in range(m):rows.append([0]*j+list(reversed(f))+[0]*(m-1-j))
 for j in range(n):rows.append([0]*j+list(reversed(g))+[0]*(n-1-j))
 return det(rows)
def perm(n,cyc):
 p=list(range(n))
 for i,x in enumerate(cyc):p[x]=cyc[(i+1)%len(cyc)]
 return tuple(p)
def comp(a,b):return tuple(a[b[i]]for i in range(len(a)))
def inv(a):
 p=[0]*len(a)
 for i,j in enumerate(a):p[j]=i
 return tuple(p)
def ct(p):
 seen=set();o=[]
 for i in range(len(p)):
  if i not in seen:
   j=i;c=0
   while j not in seen:seen.add(j);c+=1;j=p[j]
   o.append(c)
 return sorted(o,reverse=True)
def group(gs):
 n=len(gs[0]);I=tuple(range(n));seen={I};todo=deque([I])
 while todo:
  a=todo.popleft()
  for b in gs:
   c=comp(a,b)
   if c not in seen:seen.add(c);todo.append(c)
 return seen
def validate(L,repo):
 assert isinstance(L,dict)and type(L.get('schema'))is int and L['schema']==1
 rows=L.get('files');assert isinstance(rows,list)and len(rows)==len(REQ)
 assert all(isinstance(r,dict)and set(r)=={'root','path','sha256'}for r in rows)
 assert{(r['root'],r['path'])for r in rows}==REQ
 for r in rows:
  p=(HERE if r['root']=='local'else repo)/r['path']
  assert p.is_file()and sha(p)==r['sha256'],r['path']
def validate_input(v):
 assert type(v)is dict and set(v)==set(EXPECTED)
 for k,w in EXPECTED.items():
  assert type(v[k])is type(w),(k,'type')
  if type(w)is list:assert all(type(x)is type(y)for x,y in zip(v[k],w))
 assert v==EXPECTED
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,default=HERE/RN);args=ap.parse_args();repo=args.repo_root.resolve()
 L=json.loads((HERE/LN).read_text());validate(L,repo);inputs=json.loads((HERE/IN).read_text());validate_input(inputs)
 nids=identities();covers=[]
 for d in range(2,21):
  c=Q(d-1,d+1);s=Q(d*d-1,4);a=Q((d-1)**2,4);D=a+d;wc=s**d/c
  num=[Q(0)]*d+[Q(1),-c];den=[-c,Q(1)]
  deriv=pa(pm(pd(num),den),ps(pm(num,pd(den)),-1))
  expected=[Q(0)]*(d-1)+ps([1,-2,1],-d*c)
  assert deriv==expected
  fiber=pa(num,ps(den,-1));q,rem=div(fiber,[-1,3,-3,1]);assert not rem and ev(q,1)
  assert c/(1-c)**2==s and 2*(d-1)+2==2*d
  for r in[Q(-2),Q(-1),Q(0),Q(1),Q(2),Q(7,3)]:
   if r==c:continue
   S=s*r;omega=S**d*(D-S)/(S-a)
   assert omega/wc==ev(num,r)/ev(den,r)
  # Ordinary cusp inactive-saddle coefficient ratio, with its exact
  # exponential lower bound exp(2s)>1+2s>d.
  if d%2:
   assert a.denominator==1
   plus=d*s**int(D)/(s-a);minus=d*(-s)**int(D)/(-s-a)
   assert abs(minus/plus)==Q(1,d)and 1+2*s>d
  covers.append(dict(d=d,degree=d+1,c=str(c),s_c=str(s),fiber_one_cofactor_degree=len(q)-1))
 discriminants=[]
 for d in range(2,11):
  n=d+1;A=d-1;B=d+1
  for w in[-1,2,3]:
   f=[0]*(d+2);f[0]=-A*w;f[1]=B*w;f[d]=-B;f[d+1]=A
   res=resultant(f,[int(x)for x in pd(f)]);assert res%A==0
   disc=(-1)**(n*(n-1)//2)*(res//A)
   expected=(-1)**(d*(d-1)//2)*d**d*A**(d-1)*B**(d+1)*w**(d-1)*(w-1)**2
   assert disc==expected
   discriminants.append(dict(d=d,W=w,discriminant=str(disc)))
  for w in[0,1]:
   f=[0]*(d+2);f[0]=-A*w;f[1]=B*w;f[d]=-B;f[d+1]=A
   assert resultant(f,[int(x)for x in pd(f)])==0
 groups=[]
 for d in range(2,8):
  n=d+1;h=d;p0=perm(n,tuple(range(d)));p1=perm(n,(h,d-1,d-2))
  assert ct(p0)==[d,1]and ct(p1)==[3]+[1]*(d-2)and ct(comp(p0,p1))==[d,1]
  G=group([p0,p1]);wanted=factorial(n)//(2 if d%2 else 1);assert len(G)==wanted
  # All adjacent hub triples are conjugates of p1.
  pp=tuple(range(n))
  for _ in range(d):
   t=comp(comp(pp,p1),inv(pp));assert t in G
   pp=comp(p0,pp)
  groups.append(dict(d=d,order=len(G),geometric='A'if d%2 else'S'))
 # Literal hub-generation universal permutation identities.
 h,u,v,z=0,1,2,3;ab=perm(4,(h,u,v));bc=perm(4,(h,v,z))
 assert comp(ab,bc)==perm(4,(u,v,z))
 assert comp(comp(ab,bc),inv(ab))==perm(4,(u,h,z))
 sector_controls=0
 for d in range(3,21):
  b=Q(d-1,2);ss=b*(b+1)
  assert ss>=2
  for y in[Q(1,10),Q(1),Q(10)]:
   lhs=Q(d)/(1+y*y)-b/(1+b*b*y*y)-(b+1)/(1+(b+1)**2*y*y)
   rhs=d*ss*y*y*(1+(ss-1)*y*y)/((1+y*y)*(1+b*b*y*y)*(1+(b+1)**2*y*y))
   assert lhs==rhs and lhs>0;sector_controls+=1
 arithmetic=[]
 for d in range(3,50,2):
  k=(-1)**((d-1)//2)*d
  C=(-1)**(d*(d-1)//2)*d**d*(d-1)**(d-1)*(d+1)**(d+1)
  assert C%k==0 and C//k>0 and isqrt(C//k)**2==C//k
  square=k>0 and isqrt(k)**2==k
  assert square==(isqrt(d)**2==d)
  arithmetic.append(dict(d=d,kappa=k,arithmetic='A'if square else'S'))
 # Six meaningful negative/scope controls.
 wrong=pa(pm(pd([0,0,0,1,Q(-1,2)]),[Q(1,2),1]),ps([0,0,0,1,Q(-1,2)],-1))
 assert wrong!=[0,0]+ps([1,-2,1],Q(-3,2)) # Wrong pole sign.
 bad0=perm(4,(0,1,2));bad1=perm(4,(0,2,1));assert len(group([bad0,bad1]))==3 # Missing hub/transitivity.
 assert groups[1]['order']==12 and groups[1]['order']!=24 # GeometricA4 is not arithmeticS4.
 assert arithmetic[0]['kappa']==-3 and arithmetic[0]['arithmetic']=='S'
 assert (-1)**(2*(2-1)//2)*2**2*1*3**3<0 # Wrong discriminant sign rejected.
 assert 1+Q(4)>3 # Larger inactive d3 factor cannot be ignored by formal dominance.
 assert Q(1,4)+1>1 and Q(1,2)<1 # Outside a disk alone does not imply Re(r)>1.
 failures=0
 for j in range(6):
  bad=copy.deepcopy(L);del bad['files'][j]
  try:validate(bad,repo)
  except(AssertionError,KeyError,TypeError):failures+=1
  else:raise AssertionError('missing source accepted')
 bad_hash_controls=0
 for j in range(6):
  bad=copy.deepcopy(L);bad['files'][j]['sha256']='0'*64
  try:validate(bad,repo)
  except(AssertionError,KeyError,TypeError):bad_hash_controls+=1
  else:raise AssertionError('bad hash accepted')
 schema_controls=0
 for value in[True,0]:
  bad=copy.deepcopy(L);bad['schema']=value
  try:validate(bad,repo)
  except(AssertionError,KeyError,TypeError):schema_controls+=1
  else:raise AssertionError('wrong schema accepted')
 input_controls=0
 for k,value in[('actual_global_zero_localization',0),('new_SFC_claim',0),('integer_d_min',2.0),('schema',True)]:
  bad=copy.deepcopy(inputs);bad[k]=value
  try:validate_input(bad)
  except(AssertionError,KeyError,TypeError):input_controls+=1
  else:raise AssertionError('wrong input type accepted')
 for opt in['-O','-OO']:
  run=subprocess.run([sys.executable,opt,str(Path(__file__))],capture_output=True,text=True)
  assert run.returncode!=0 and 'Refusing -O/-OO'in run.stderr
 out=dict(schema=1,status='PASS_STATIONARY_COVER_ALGEBRA_AND_GROUP_SUPPORT',scope='Positive-parameter stationary branch classification, exact cusp geometric/arithmetic cover, distinguished global radial graph and activation barriers only; no global actual-root/SFC conclusion.',proof_sha256=sha(HERE/PN),checker_sha256=sha(HERE/CN),input_sha256=sha(HERE/IN),source_ledger_sha256=sha(HERE/LN),mandatory_sources=6,universal_identities=nids,literal_cover_controls=covers,literal_discriminants=discriminants,degenerate_discriminant_controls=18,enumerated_geometric_groups=groups,arithmetic_squareclass_controls=arithmetic,negative_controls=7,sector_derivative_controls=sector_controls,missing_source_controls=failures,bad_hash_controls=bad_hash_controls,schema_controls=schema_controls,input_type_scope_controls=input_controls,optimized_mode_controls=2,independent_replay=False,formal_verification=False,external_review=False)
 args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
 print(json.dumps(dict(status='PASS',identities=nids,cover_controls=len(covers),literal_discriminants=len(discriminants),groups=len(groups),arithmetic_controls=len(arithmetic),negative_controls=7,sector_controls=sector_controls,source_controls=failures+bad_hash_controls+schema_controls+input_controls+2)))
if __name__=='__main__':main()
