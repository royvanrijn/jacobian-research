#!/usr/bin/env python3
"""Exact symbolic replay for the two-sided protected-U proof.
The unrestricted proof is in TWO_SIDED_PROTECTED_OBSTRUCTIONS.md. This program
checks its propagation and scalar identities, and an actual sharp U coloring;
it is not a search replacing that proof.
"""
from __future__ import annotations
import argparse, hashlib, json, platform
from pathlib import Path
import sympy as s

mu,x,X,y,Y,z,Z=s.symbols('mu x X y Y z Z')
eta=1-mu
G=s.Matrix([[eta,x,y],[X,eta,z],[Y,Z,eta]])
e=[s.eye(3)[:,i] for i in range(3)]
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
checks=[]
def zero(a):
    if isinstance(a,s.MatrixBase): return all(zero(b) for b in a)
    return s.cancel(a)==0
def require(condition,name):
    if not condition: raise AssertionError(name)
    checks.append(name)
def pairing(f,u): return (f.T*G*u)[0]
def conjugate(o,a,inverse=False):
    u,f=o; v,g=a
    if inverse: return (v+u*pairing(f,v)/mu, g-f*pairing(g,u))
    return (v-u*pairing(f,v), g+f*pairing(g,u)/mu)
def same(a,b): return zero(a[0]-b[0]) and zero(a[1]-b[1])
A=(e[0],e[0]); B=(e[1],e[1]); C=(e[2],e[2])
p15=conjugate(C,B)
p2=conjugate(p15,C)
p3=conjugate(A,p2,True)
p4=conjugate(B,p3)
p13=conjugate(p3,B,True)
h=1+z*Z/mu
t=-z*x/mu+h*y
ell=Z*X+h*Y
kappa=t*X/mu+z*(h-eta/mu)
Lambda=-ell*x+(eta+h)*Z
expected={
 2:(s.Matrix([0,-z/mu,h]),s.Matrix([0,Z,h])),
 3:(s.Matrix([t/mu,-z/mu,h]),s.Matrix([-ell,Z,h])),
 4:(s.Matrix([t/mu,-z/mu-kappa,h]),s.Matrix([-ell,Z+Lambda/mu,h])),
 13:(s.Matrix([Lambda*t/mu**2,1-Lambda*z/mu**2,Lambda*h/mu]),
     s.Matrix([kappa*ell,1-kappa*Z,-kappa*h]))}
for i,pair in [(2,p2),(3,p3),(4,p4),(13,p13)]:
    require(same(pair,expected[i]),f'independent formal propagation row {i}')
    require(zero(pairing(pair[1],pair[0])-eta),f'row {i} displacement normalization')
# Mutation: reversing the M3 crossing must change the propagated pair.
require(not same(conjugate(A,p2),expected[3]),'reject reversed M3 crossing')
# These evaluations use only the coefficient consequences kappa*Z=1,
# Lambda=-mu*Z, not independence of the formal root vectors.
f13_reduced=s.Matrix([kappa*ell,0,-kappa*h])
S=pairing(f13_reduced,expected[2][0])
require(zero(S-kappa*(ell*t-h+mu*h*h)),'S identity after kappa Z=1')
u4_reduced=s.Matrix([t/mu,-z/mu-1/Z,h])
Q=pairing(e[2],u4_reduced)
require(zero(Q-(t*Y/mu-mu*h)),'Q identity after kappa Z=1')
require(zero(mu*(kappa*Z-1)-(t*X*Z+mu**2*h**2-mu*h-mu**2)),
        'kappa definition gives scalar equation C')
require(zero(mu*(Q*h+mu)-(t*Y*h-mu**2*(h*h-1))),
        'outer covector equation gives scalar equation B')
require(zero((t*ell-mu*h)-((t*X*Z+mu**2*h**2-mu*h-mu**2)
                 +(t*Y*h-mu**2*(h*h-1)))), 'B plus C gives t ell=mu h')
# From equations A,D: (h-mu h²)-mu h=h(eta-mu h).
H,M,P=s.symbols('h M P')
require(zero((H-M*H*H)-M*H-H*(1-M-M*H)),
        'A minus D forces h=eta/mu')
# Independent short last step, using f13=-kappa f4, hence f4(u2)=0.
f4_reduced=s.Matrix([-ell,0,h])
require(zero(pairing(expected[4][1],e[1])+Z-(Lambda+mu*Z)/mu),
        'f4(v)=-Z modulo Lambda=-mu Z')
require(zero((-z/mu)*(-Z)+h*P-(h-1+h*P)),
        'S=0 gives h(P+1)=1')
forced_h=(1-M)/M
forced_P=forced_h-1-1/M
require(zero(forced_P+2),'first root relation evaluated by g forces P=-2')
require(zero((1+forced_h)-1/M),'1+h=1/mu under mu h=eta')
# The contradiction uses only integer coefficients: characteristic two is included.
require(zero((H*(P+2)-(1+H))-(H*(P+1)-1)),
        'last contradiction identity with integer coefficients')
require(not zero((H*(P+2)-(1-H))-(H*(P+1)-1)),
        'reject corrupted sign in final identity')

# Replay all three preliminary second-fold conjugations in dimension two.
old_G=G
x1,y1=s.symbols('x1 y1')
G=s.Matrix([[eta,y1],[x1,eta]])
d0=s.Matrix([1,0]); p0=s.Matrix([0,1])
Dpair=(d0,d0); Ppair=(p0,p0)
Npair=conjugate(Ppair,Dpair,True)
Qpair=conjugate(Npair,Ppair,True)
Cpair=conjugate(Qpair,Npair,True)
h1=1+x1*y1/mu
for name,actual,expected_fold in [
 ('N',Npair,(d0+(x1/mu)*p0,d0-y1*p0)),
 ('Q',Qpair,(y1*d0+h1*p0,-(x1/mu)*d0+h1*p0)),
 ('C',Cpair,(h1*d0+(x1/mu)*(1+h1)*p0,h1*d0-y1*(1+h1)*p0))]:
    require(same(actual,expected_fold),'second-fold conjugation '+name)
G=old_G
# The first fold uses only its covector comparison, then applies f to a
# root equality; these are rational identities before the final comparison.
c,lam=s.symbols('c lam')
require(zero(lam*(mu/lam+(mu/(lam*c))*c*eta/mu)-1),
        'first-fold application of f gives mu+eta=1')
require(zero(eta-(-lam)*mu/lam-1),'first-fold nu(a)=-lambda')

# A genuine U representation attains both local span dimensions three.
# Roots in Q(sqrt(5)); all reflections have determinant -1.
phi=(1+s.sqrt(5))/2
roots=[s.Matrix([-2,0,0]),s.Matrix([-1,phi-1,-phi]),
       s.Matrix([phi-1,phi,-1]),s.Matrix([phi-1,phi,-1])]
matrices=[s.eye(3)-v*v.T/2 for v in roots]
U=[1,2,-3,-3,2,-3,-3,2,1,1,2,2,2,-3]
def clean(a): return a.applyfunc(s.simplify)
def hurwitz(values,braid):
    values=list(values)
    for letter in braid:
        i=abs(letter)-1; a,b=values[i:i+2]
        values[i:i+2]=[b,clean(b.inv()*a*b)] if letter>0 else [clean(a*b*a.inv()),a]
    return values
after=hurwitz(matrices,U)
require(zero(after[0]-matrices[0]) and zero(after[3]-matrices[3]),
        'actual rank-three reflection U tuple has both outer returns')
require(s.Matrix.hstack(*roots).rank()==3,'actual U roots span three')
require(s.Matrix.vstack(*(v.T/2 for v in roots)).rank()==3,'actual U covectors span three')
require(all(zero(m.det()+1) and (s.eye(3)-m).rank()==1 for m in matrices),
        'actual U meridians have determinant -1 and displacement rank one')
# Root/covector transpose-inverse interchange is an exact identity.
u=s.Matrix(s.symbols('u0:3')); f=s.Matrix(s.symbols('f0:3'))
D=s.eye(3)-u*f.T
require(zero((s.eye(3)+f*u.T/(1-(f.T*u)[0]))*D.T-s.eye(3)),
        'contragredient interchanges root and covector lines')
sources=[Path(__file__),HERE/'TWO_SIDED_PROTECTED_OBSTRUCTIONS.md',
         HERE/'MULTIPLE_PROTECTED_OBSTRUCTIONS.md',
         ROOT/'results/prime-knot-augmentation-gaps/representations.tex']
receipt={'status':'passed',
 'scope':'exact symbolic support for unrestricted written two-sided proof, plus actual sharp rank-three U control',
 'python':platform.python_version(),'sympy':s.__version__,
 'arithmetic':'rational functions over Q; exact Q(sqrt(5)) reflection control',
 'checks':checks,'check_count':len(checks),
 'source_sha256':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
 'assurance':{'internal_replay':True,'external_review':False,'formal_verification':False},
 'does_not_establish':['prime realization','hyperbolicity','compatible transitive global coloring','augmentation gap two']}
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--write-receipt',type=Path)
args=parser.parse_args()
if args.write_receipt:
    args.write_receipt.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':'passed','check_count':len(checks)}))
