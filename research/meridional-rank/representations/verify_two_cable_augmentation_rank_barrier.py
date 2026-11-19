#!/usr/bin/env python3
"""Exact controls for the all-parameter two-cable rank lower bound.
The written proof supplies the unrestricted parameter induction and cycle
formula; finite b/n checks below do not replace that proof.
"""
from pathlib import Path
from math import gcd
import argparse,hashlib,json,platform
HERE=Path(__file__).resolve().parent
if not __debug__:
 raise SystemExit('Run without -O: exact assertion checks must be enabled.')
U=[1,2,-3,-3,2,-3,-3,2,1,1,2,2,2,-3]
def mul(a,b):return tuple(a[b[i]]for i in range(len(a)))
def inv(a):
 q=[0]*len(a)
 for i,j in enumerate(a):q[j]=i
 return tuple(q)
def power(p,n):
 if n<0:return power(inv(p),-n)
 q=tuple(range(len(p)))
 while n:
  if n&1:q=mul(q,p)
  p=mul(p,p);n//=2
 return q
def conj(a,p):return mul(mul(p,a),inv(p))
def tr(n,a,b):
 p=list(range(n));p[a]=b;p[b]=a;return tuple(p)
def hz(t,word):
 t=list(t)
 for l in word:
  i=abs(l)-1;a,b=t[i:i+2]
  t[i:i+2]=[b,conj(a,inv(b))]if l>0 else[conj(b,a),a]
 return tuple(t)
def cable(word):
 answer=[]
 for l in word:
  j=abs(l);part=[2*j,2*j-1,2*j+1,2*j]
  answer.extend(part if l>0 else[-i for i in part[::-1]])
 return answer
def packet_conj(P,a):return tuple(conj(p,a)for p in P)
def packet_hz(t,word):
 t=list(t)
 for l in word:
  i=abs(l)-1;a,b=t[i:i+2];A=mul(*a);B=mul(*b)
  t[i:i+2]=[b,packet_conj(a,inv(B))]if l>0 else[packet_conj(b,A),a]
 return tuple(t)
def packets(b):
 n=2*b;pairs=[[(0,1),(2,3)],[(2,3),(1,4)],[(0,5),(1,6)],[(1,2),(6,7)]]
 if b>=5:pairs.append([(0,8),(7,9)])
 for j in range(5,b):pairs.append([(2*j-2,2*j),(2*j-1,2*j+1)])
 return tuple(tuple(tr(n,a,c)for a,c in pair)for pair in pairs)
def flat(P):return tuple(p for pair in P for p in pair)
def edge(p):return next((i,j)for i,j in enumerate(p)if i<j and p[j]==i)
def graph_connected(t):
 n=len(t[0]);seen={0};q=[0]
 while q:
  x=q.pop()
  for p in t:
   if p[x]not in seen:seen.add(p[x]);q.append(p[x])
 return len(seen)==n
# Verify one positive and negative 2-parallel block in a FREE group.
def red(w):
 q=[]
 for i in w:
  if q and q[-1]==-i:q.pop()
  else:q.append(i)
 return tuple(q)
def wi(w):return tuple(-i for i in w[::-1])
def wm(*ws):return red(i for w in ws for i in w)
def whz(t,word):
 t=list(t)
 for l in word:
  i=abs(l)-1;a,b=t[i:i+2]
  t[i:i+2]=[b,wm(wi(b),a,b)]if l>0 else[wm(a,b,wi(a)),a]
 return tuple(t)
f=tuple((i,)for i in range(1,5));A=wm(f[0],f[1]);B=wm(f[2],f[3])
assert whz(f,[2,1,3,2])==(f[2],f[3],wm(wi(B),f[0],B),wm(wi(B),f[1],B))
assert whz(f,[-2,-3,-1,-2])==(wm(A,f[2],wi(A)),wm(A,f[3],wi(A)),f[0],f[1])
# The four-bridge induction uses two exact finite identities.
T=flat(packets(4));E=hz(T,cable(U+[2]*3))
assert E==tuple([T[1],T[0],*T[2:]])
assert hz(E,cable([2]*6))==E
assert hz(E,[1])==T and hz(E,[-1])==T
controls=0
for n in range(8):
 for q in [-5,-1,1,3,7]:
  beta=U+[2]*(3+6*n)
  assert hz(T,cable(beta)+([1]*q if q>0 else[-1]*(-q)))==T
  assert packet_hz(packets(4),beta)==tuple(tuple(E[2*j:2*j+2])for j in range(4))
  controls+=1
# Uniform b,k controls; the note proves the formula for all b by explicit p cycles.
for b in range(5,13):
 P=packets(b);top=flat(P);d=b-2;delta=[2,3,4,-3]+list(range(5,b))
 after=packet_hz(P,U)
 assert after[0]==P[0][::-1]and after[3:]==P[3:]
 C=mul(*P[3]);framed=(after[1],after[2])+tuple(packet_conj(a,C)for a in P[4:])
 p=tuple(range(2*b))
 for a in framed:p=mul(p,mul(*a))
 assert packet_conj(P[1],p)==after[2]
 assert packet_conj(P[2],power(p,d))==after[1]
 for i,z in enumerate(framed[2:],1):assert z==packet_conj(P[2],power(p,i))
 assert power(p,d*(d+1))==tuple(range(2*b))
 for k in range(3):
  beta=U+delta*((d+1)*(1+d*d*k))
  result=hz(top,cable(beta))
  assert result==tuple([top[1],top[0],*top[2:]])
  assert result==flat(packet_hz(P,beta))
  assert hz(result,[1])==top and hz(result,[-1])==top
  assert graph_connected(top)
  controls+=1
# Independent ORIGINAL companion transport plus torus-pattern group relations.
P=tuple(mul(*a)for a in packets(4));xs=list(P);src=list(range(4));js=[tuple(range(8))]*4
beta=U+[2]*3
for l in beta:
 i=abs(l)-1;a,b=xs[i:i+2];sa,sb=src[i:i+2]
 if l>0:js[sa]=mul(js[sa],b);xs[i:i+2]=[b,conj(a,inv(b))]
 else:js[sb]=mul(js[sb],inv(a));xs[i:i+2]=[conj(b,a),a]
 src[i:i+2]=[sb,sa]
assert tuple(xs)==P
j=0;J=tuple(range(8));seen=[]
while j not in seen:seen.append(j);J=mul(J,js[j]);j=src.index(j)
assert len(seen)==4
assert J==(3,2,0,1,5,7,4,6)
w=sum(1 if a>0 else-1 for a in beta);M=P[0];L=mul(J,power(M,-w))
assert L==(2,3,1,0,5,7,4,6)
assert mul(L,M)==mul(M,L)
x,y=packets(4)[0]
assert mul(x,y)==M
assert conj(x,inv(L))==y and conj(y,inv(L))==x
for Q in [-13,-1,1,15,27]:
 action=hz((x,y),[1]*Q if Q>0 else[-1]*(-Q))
 assert action==(conj(x,inv(L)),conj(y,inv(L)))
# Odd-p central fiber acts nonscalarly on im(M-I): e0-e1 goes to +/-e2-+e3.
v=(1,-1,0,0,0,0,0,0)
def act(p,v):return tuple(v[inv(p)[i]]for i in range(8))
for p in [3,5,7,9]:
 for Q in range(1,p):
  if gcd(p,Q)!=1:continue
  Z=mul(power(M,Q),power(L,-p));zv=act(Z,v)
  assert zv not in[v,tuple(-a for a in v)]
# Wrong crossing, parity and packet order must be rejected.
assert hz(T,cable(beta))!=T
assert hz(T,cable(beta)+[1,1])!=T
rejected=0
for i in range(len(cable(beta))):
 bad=cable(beta);bad[i]*=-1
 rejected+=hz(T,bad+[1])!=T
assert rejected>0
sources=[Path(__file__),HERE/'TWO_CABLE_AUGMENTATION_RANK_BARRIER.md',HERE/'two_cable_augmentation_rank_barrier_source_ledger.json']
ledger=json.loads((HERE/'two_cable_augmentation_rank_barrier_source_ledger.json').read_text())
for entry in ledger['sources']+ledger['repository_inputs']:
 p=HERE/entry['file']
 assert hashlib.sha256(p.read_bytes()).hexdigest()==entry['sha256']
 sources.append(p)
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--write-receipt',type=Path,default=None)
args=parser.parse_args()
receipt={'status':'passed','python':platform.python_version(),'arithmetic':'exact finite permutations and freely reduced words',
 'source_sha256':{p.relative_to(HERE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()for p in sources},
 'four_bridge_period_identity':True,'uniform_all_parameter_scope':'written ordered-packet and explicit cycle proof; finite b4..12/n or k controls supplement it',
 'parameter_controls':controls,'independent_cabling_methods':['literal4crossing replacement','ordered packet Hurwitz action','original companion peripheral transport + mappingtorus amalgam'],
 'independent_group_case':{'b':4,'n':0,'blackboard_q':1,'Seifert_Q':15,'meridian':M,'longitude':L,'transport':J},
 'conclusion':'ar_C(2-cable)>=2b-1; mr<=2b; actual gap<=1 for every released family member and every odd slope',
 'odd_p_scope':'same natural companion representation cannot extend to a KCH representation; no bound for other companion representations',
 'external_review':False,'formal_verification':False}
if args.write_receipt:args.write_receipt.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':'passed','parameter_controls':controls,'group_case_Q':15,'odd_p_negative_control':True}))

