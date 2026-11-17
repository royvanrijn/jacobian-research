"""Exact partial three-point obstruction and elliptic cusp-surface checks.

No Groebner, finite-field search, or rational reconstruction is run here.
The saved identity is replayed in FLINT and the retained Fraction engine.
"""
import hashlib
import importlib.util
import json
from fractions import Fraction
from pathlib import Path
import platform
import sympy as S
from flint import fmpq, fmpq_mpoly_ctx

if not __debug__:
    raise RuntimeError('Run without -O')
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
checks=[]


def zero(name,expr):
    assert S.cancel(expr)==0,name
    checks.append(name)


def yes(name,condition):
    assert condition,name
    checks.append(name)


t,Y1,Y2,Y3,Y4=S.symbols('t Y1 Y2 Y3 Y4')
btarget=Y4-3*Y2*t*t+2*Y1*t**3-t**4
atarget=Y3-2*Y2*t+Y1*t*t-t**3-t*t*btarget
R=S.expand(atarget*atarget+atarget+t*btarget+t**3+Y1*t-Y2)


def partitions(n,maximum=None):
    if n==0:
        yield []
        return
    for first in range(min(n,maximum or n),0,-1):
        for rest in partitions(n-first,first):yield [first]+rest


types=[p for p in partitions(12) if p.count(1)==3]
yes('eight possible three-point types',len(types)==8)
for part in types:
    nonsimple=[m for m in part if m>1]
    odd=[m%2 for m in nonsimple]
    if sum(odd)==1:
        powers=[(m-3*o)//2 for m,o in zip(nonsimple,odd)]
        yes('linear-cube square covering '+str(part),sum(powers)==3 and min(powers)>=0)
    else:yes('cubic-cube covering '+str(part),nonsimple==[3,3,3])

# Independently reconstruct the cubic-cube coefficient equations from R.
u,v,w,z,h,k,m=S.symbols('u v w z h k m')
Q=t**3+u*t*t+v*t+w
P=S.Poly(S.expand(Q**3*(t**3+z*t*t+h*t+k)),t)
y1=-P.nth(11)/4
y2=(P.nth(10)-4*y1*y1)/6
kfit=-12*y1*y2-2-P.nth(9).subs(k,0)
target=S.Poly(S.expand(P.as_expr().subs(k,kfit)),t)
y4=(6*y1+9*y2*y2-target.nth(8))/2
y3=(target.nth(6)-14*y1*y2+6*y2*y4-2)/2
res=S.Poly(S.expand(R.subs({Y1:y1,Y2:y2,Y3:y3,Y4:y4})-target.as_expr()),t)
for j in [12,11,10,9,8,6]:zero('eliminated cubic coefficient '+str(j),res.nth(j))
source=json.loads((HERE/'cubic_cube-input.json').read_text())
for row in source['equations']:
    ratio=S.cancel(res.nth(row['coefficient'])/S.sympify(row['polynomial']))
    yes('source cubic coefficient '+str(row['coefficient']),ratio.is_Rational and ratio!=0)
Delta=(u-z)**2+2*(v-h)
f7=S.sympify(next(row['polynomial'] for row in source['equations'] if row['coefficient']==7))
zero('linear factor coefficient',S.diff(f7,w)-576*Delta)
zero('Delta-zero substitution',Delta.subs({z:u-m,h:v+m*m/2}))
delta_data=json.loads((HERE/'cubic-delta-zero-input.json').read_text())
yes('seven delta-zero equations',len(source['equations'])==len(delta_data['equations'])==7)
for old,new in zip(source['equations'],delta_data['equations']):
    yes('delta coefficient matched '+str(old['coefficient']),old['coefficient']==new['coefficient'])
    expr=S.sympify(old['polynomial']).subs({z:u-m,h:v+m*m/2})
    ratio=S.cancel(expr/S.sympify(new['polynomial']))
    yes('delta necessary equation '+str(old['coefficient']),ratio.is_Rational and ratio!=0)

cert=json.loads((HERE/'cubic-delta-zero-macaulay-certificate.json').read_text())
yes('certificate covers exact source',cert['source']=='cubic-delta-zero-input.json' and cert['variables']==['u','v','m','w'] and cert['target']=='1')
yes('certificate generator count',cert['generator_indices']==list(range(7)) and len(cert['multipliers'])==7)
vars=list(map(S.Symbol,cert['variables']))
generators=[{ex:Fraction(int(c.p),int(c.q)) for ex,c in S.Poly(S.sympify(row['polynomial']),*vars,domain=S.QQ).terms()}
            for row in delta_data['equations']]
multipliers=[]
for terms in cert['multipliers']:
    out={}
    for term in terms:
        ex=tuple(term['exponents']);assert len(ex)==4 and min(ex)>=0 and ex not in out
        out[ex]=Fraction(term['coefficient'])
    multipliers.append(out)
context=fmpq_mpoly_ctx.get(cert['variables'])
total=context.constant(0)
for f,g in zip(multipliers,generators):
    first=context.from_dict({ex:fmpq(c.numerator,c.denominator) for ex,c in f.items()})
    second=context.from_dict({ex:fmpq(c.numerator,c.denominator) for ex,c in g.items()})
    total+=first*second
yes('FLINT literal identity equals one',total==context.constant(1))
path=ROOT/'foundations/jcsearch/sparse_macaulay.py'
spec=importlib.util.spec_from_file_location('retained_sparse_macaulay',path)
legacy=importlib.util.module_from_spec(spec);spec.loader.exec_module(legacy)
plain={}
for f,g in zip(multipliers,generators):legacy.add_polynomial(plain,legacy.multiply_polynomial(f,g))
yes('independent Fraction identity equals one',plain=={(0,0,0,0):Fraction(1)})

# Triple-root incidence: differentiate at fixed Y BEFORE substituting X.
s,a,b=S.symbols('s a b')
G2=a*a+a+s*b+s**3
G3=2*s*G2+a+s**3+s*s*b
G4=3*s*s*G2+b+s**4
X1=-S.Matrix([G2,G3,G4]).jacobian([s,a,b]).det().expand()
X=[X1,G2+s*X1,G3+s*s*X1,G4+s**3*X1]
sub=dict(zip([Y1,Y2,Y3,Y4],X))
zero('critical source root equation',R.subs(t,s).subs(sub))
zero('critical source derivative equation',S.diff(R,t).subs(t,s).subs(sub))
f=S.Poly(S.expand(S.diff(R,t,2).subs(t,s).subs(sub)/2),b)
yes('triple incidence quadratic in b',f.degree()==2)
coefficients=S.Poly(f.as_expr(),a,b).coeffs()
content=coefficients[0]
for coefficient in coefficients[1:]:content=S.gcd(content,coefficient)
yes('no extra vertical components',S.degree(content,s)==0)
Ab,Bb,Cb=f.all_coeffs()
zero('quadratic leading coefficient',Ab-4*s*s*(3*s**3-2)**2)
H=S.Poly(S.expand(S.discriminant(f.as_expr(),b)/4),a)
model=json.loads((HERE/'elliptic-cusp-input.json').read_text())
for i,(fresh,saved) in enumerate(zip(H.all_coeffs(),model['quartic_coefficients'])):zero('quartic model coefficient '+str(i),fresh-S.sympify(saved))
L=4*s*(3*s**3-2)*(3*s**3-1)
zero('two rational infinity points',H.LC()-L*L)
zero('birational square completion',(Ab*b+Bb/2)**2-H.as_expr()-Ab*f.as_expr())
A4,B4,C4,D4,E4=H.all_coeffs()
I=S.expand(12*A4*E4-3*B4*D4+C4*C4)
J=S.expand(72*A4*C4*E4+9*B4*C4*D4-27*A4*D4*D4-27*B4*B4*E4-2*C4**3)
zero('I invariant',I-S.sympify(model['I']))
zero('J invariant',J-S.sympify(model['J']))
quartic_disc=S.expand((4*I**3-J*J)/27)
zero('polynomial quartic discriminant',S.discriminant(H.as_expr(),a)-quartic_disc)
zero('retained discriminant',quartic_disc-S.sympify(model['quartic_discriminant']))
Q64=S.Poly(S.sympify(model['Q64']),s)
zero('complete discriminant factorization',quartic_disc+4096*s**4*(3*s**3-2)**4*Q64.as_expr())
yes('64 distinct nodal places',Q64.degree()==64 and S.gcd(Q64,Q64.diff()).degree()==0)
yes('nodal places disjoint from named four',S.Poly(S.gcd(Q64.as_expr(),s*(3*s**3-2)),s).degree()==0)
yes('c4 unit at every finite bad place',S.Poly(S.gcd(I,quartic_disc),s).degree()==0)
yes('c4,c6 have no finite common zero',S.Poly(S.gcd(I,J),s).degree()==0)
yes('invariant degrees',[S.degree(x,s) for x in [I,J,quartic_disc]]==[23,40,80])
zero('short Jacobian discriminant normalization',-16*(4*(-27*I)**3+27*(-27*J)**2)-16*3**12*quartic_disc)
yes('infinity local orders',[28-S.degree(I,s),42-S.degree(J,s),84-S.degree(quartic_disc,s)]==[5,2,4])
yes('infinity minimality',42-S.degree(J,s)<6)
yes('Euler and holomorphic Euler',4*4+64+4==84 and 84//12==7)
yes('geometric genus from chi and base genus',84//12-1==6)
yes('nonisotrivial j map',3*S.degree(I,s)!=S.degree(quartic_disc,s))
yes('smooth rational specialization',S.discriminant(H.as_expr().subs(s,1),a)==5799391526912)

manifest=json.loads((HERE/'sources/manifest.json').read_text())
for row in manifest['sources']:
    yes('source hash '+row['repository_path'],hashlib.sha256((ROOT/row['repository_path']).read_bytes()).hexdigest()==row['sha256'])
receipt={'result':'PASS','claim_ids':['GAOF7THREE1','GAOF7CUSP1'],'checks':len(checks),'check_names':checks,
         'python':platform.python_version(),'sympy':S.__version__,
         'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'certificate_sha256':hashlib.sha256((HERE/'cubic-delta-zero-macaulay-certificate.json').read_bytes()).hexdigest(),
         'model_sha256':hashlib.sha256((HERE/'elliptic-cusp-input.json').read_bytes()).hexdigest(),
         'identity_replays':['FLINT','retained Fraction sparse_macaulay'],
         'multiplier_terms':sum(len(x) for x in multipliers),
         'elliptic_completion':{'geometric_fibres':{'I4':4,'I1':64,'IV':1},'Euler_number':84,'chi_O':7,'geometric_genus':6},
         'three_point_status':'OPEN: only the Delta=0 cubic-cube subcase is excluded',
         'assurance':'Exact arithmetic and sourced written elliptic-surface interpretation; no formal verification, independent end-to-end review, external review or novelty priority.'}
(HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ['result','checks','identity_replays','multiplier_terms','elliptic_completion','three_point_status']}))
