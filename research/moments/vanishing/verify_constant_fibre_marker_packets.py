"""One-use exact polynomial transcription check. Never replay after success.

The universal theorem is the written proof, not this bounded check.
Every retained identity is checked in the ambient polynomial ring,
with explicit spatial IBP witnesses and normal-ordered parameter operators.
"""
from pathlib import Path
from functools import lru_cache
import hashlib, json, time
import sympy as S

root=Path(__file__).resolve().parent
output=root/'constant_fibre_marker_packet_verification.json'
assert not output.exists(), 'Successful packet is sealed;do not replay.'
started=time.monotonic()
t,a,x,z,y=S.symbols('t a x z y')
U=(x,z)
B=x*x+x*z+1
g=x*x+z**3
P=a*y*y+B*y+g
s=1-2*a*t
u=1-a*t
zero=(S.Integer(0),)*3
ex=S.expand

def add(*packets):
    rep={};w=list(zero)
    for rr,ww in packets:
        for key,v in rr.items():rep[key]=ex(rep.get(key,0)+v)
        w=[ex(v+q) for v,q in zip(w,ww)]
    return {k:v for k,v in rep.items() if v!=0},tuple(w)

def scale(packet,c):
    rr,ww=packet
    assert not S.sympify(c).has(x,z,y)
    return {k:ex(c*v) for k,v in rr.items()},tuple(ex(c*v) for v in ww)

def nabla(packet):
    rr,ww=packet;rep={}
    for (powers,k),v in rr.items():
        rep[(powers,k)]=ex(rep.get((powers,k),0)+S.diff(v,t))
        rep[(powers,k+1)]=ex(rep.get((powers,k+1),0)+v)
    return {k:v for k,v in rep.items() if v!=0},tuple(ex(S.diff(v,t)+P*v) for v in ww)

def base(Q):
    Q=ex(Q);assert not Q.has(y)
    return {((i,j),0):c for (i,j),c in S.Poly(Q,x,z).terms()},zero

def witness(i,v):
    w=list(zero);w[i]=ex(v);return {},tuple(w)

@lru_cache(None)
def conditional(m,Q):
    if m==0:return base(Q)
    if m==1:return add(base(t*B*Q),witness(2,-Q))
    return add(scale(conditional(m-1,ex(B*Q)),t),
               scale(conditional(m-2,Q),(m-1)*s),
               witness(2,-s**(m-1)*Q*y**(m-1)))

@lru_cache(None)
def connection(j,Q):
    if j==1:
        return add(scale(nabla(base(Q)),t),base(Q-t*(g*Q+S.diff(Q,t))),witness(2,-Q*y))
    k=j-1;p=2*k-1
    return add(scale(nabla(connection(k,Q)),t*u),
               scale(connection(k,Q),t*p*a+(2*k+1)*u),
               scale(connection(k,ex(g*Q+S.diff(Q,t))),-t*u),
               witness(2,-u**(p+1)*Q*y**(2*k+1)))

@lru_cache(None)
def even(j,Q):
    if j==0:return base(Q)
    p=2*j-1;q=2*j;K=p+q-1
    aa=sum(S.binomial(K,v)*2**v*(-1)**(K-v)*u**(v-p)*s**(K-v) for v in range(p,K+1))
    bb=sum(S.binomial(K,v)*2**v*(-1)**(K-v)*u**v*s**(K-v-q) for v in range(p))
    assert ex(aa*u**p+bb*s**q-1)==0
    return add(scale(connection(j,Q),ex(aa)),scale(conditional(2*j,Q),ex(bb)))

def derivative(alpha):
    return S.diff(B,x,alpha[0],z,alpha[1])

@lru_cache(None)
def source(j,alpha,Q):
    if alpha==(0,0):
        return add(scale(even(j+1,Q),s),scale(even(j,Q),-(2*j+1)),witness(2,Q*y**(2*j+1)))
    i=next(i for i,n in enumerate(alpha) if n)
    prev=list(alpha);prev[i]-=1;prev=tuple(prev)
    F=derivative(prev);Bi=S.diff(B,U[i])
    LQ=ex((U[i]-t*S.diff(g,U[i]))*Q-S.diff(Q,U[i]))
    return add(source(j,prev,LQ),scale(even(j+1,ex(Bi*F*Q)),-t*t),witness(i,t*F*Q*y**(2*j+1)))

@lru_cache(None)
def coefficient(j,alpha,Q):
    packets=[]
    for gamma in [(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]:
        beta=tuple(v+w for v,w in zip(alpha,gamma))
        if derivative(beta)==0:continue
        mult=(-x)**gamma[0]*(-z)**gamma[1]
        mult/=S.factorial(alpha[0])*S.factorial(alpha[1])*S.factorial(gamma[0])*S.factorial(gamma[1])
        packets.append(source(j,beta,ex(mult*Q)))
    return add(*packets)

def odd(j,alpha,Q):
    c=S.Poly(B,x,z).coeff_monomial(x**alpha[0]*z**alpha[1])
    return add(nabla(coefficient(j,alpha,Q)),
               scale(coefficient(j,alpha,ex(S.diff(Q,t)+g*Q)),-1),
               scale(coefficient(j+1,alpha,Q),-a),
               scale(even(j+1,ex(B*Q)),-c*t))

def evaluate(packet):
    return ex(sum(c*x**powers[0]*z**powers[1]*P**k for (powers,k),c in packet[0].items()))

def spatial(ww):
    return ex(S.diff(ww[0],x)+(t*S.diff(P,x)-x)*ww[0]
             +S.diff(ww[1],z)+(t*S.diff(P,z)-z)*ww[1]
             +S.diff(ww[2],y)+(-s*y+t*B)*ww[2])

def residual(target,packet):
    return ex(target-evaluate(packet)-spatial(packet[1]))

checks=[];saved=[]
def check(name,target,packet,retain=False):
    assert residual(target,packet)==0,name
    checks.append(name)
    if retain:
        saved.append(dict(name=name,target=str(target),normal_ordered=[dict(base_powers=list(powers),parameter_order=k,time_coefficient=str(c)) for (powers,k),c in sorted(packet[0].items())],spatial_witnesses={str(v):str(w) for v,w in zip((x,z,y),packet[1])}))
    print('PASS',name,flush=True)

Q=x+2*z
check('conditional-degree4',s**4*Q*y**4,conditional(4,Q))
check('even-degree2',Q*y**2,even(1,Q),True)
check('even-degree4',Q*y**4,even(2,Q))
check('coefficient-constant-Taylor',t*Q*y,coefficient(0,(0,0),Q))
check('coefficient-quadratic-Taylor',t*Q*y,coefficient(0,(2,0),Q))
packet=odd(0,(0,0),Q)
check('odd-degree1-global-no-time-denominator',Q*y,packet,True)
bad=add(packet,scale(even(1,ex(B*Q)),t))
assert residual(Q*y,bad)!=0,'Mutation deleting the final even term must fail.'
badw=list(packet[1]);badw[0]=0
assert residual(Q*y,(packet[0],tuple(badw)))!=0,'Mutation deleting a base spatial witness must fail.'
assert (B.subs({x:0,z:0}))==1
data=dict(date='2026-11-17',status='pass',purpose='Bounded exact ambient polynomial transcription of universal written packet theorem;not an independent proof of universal quantifiers.',example=dict(base_variables=['x','z'],fibre='y',quadratic_coefficient='a (symbolic,including a=0)',source=str(B),base_phase=str(g),marker=str(Q)),checks=checks,rejected_mutations=2,retained_packets=saved,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),elapsed_seconds=time.monotonic()-started,independent_replay=False,formal_verification=False,external_review=False,rerun='Successful packet sealed;never rerun this checker.')
output.write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps(dict(status='pass',checks=len(checks),rejected_mutations=2,seconds=data['elapsed_seconds'],artifact=str(output))),flush=True)
