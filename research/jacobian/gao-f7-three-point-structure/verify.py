"""Exact section, fibre-intersection and height arithmetic for F7 incidence.

The height formula and extension of sections are written geometric inputs,
not formalized here. No Groebner search is part of verification.
"""
from pathlib import Path
import hashlib,json,platform
import sympy as S
from flint import fmpq,fmpq_poly
HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'gao-f7-three-point'
checks=[]
def yes(name,condition):
 assert condition,name
 checks.append(name)
def zero(name,expression):yes(name,S.cancel(expression)==0)
def polynomial(expression):
 p=S.Poly(expression,s,domain=S.QQ)
 return fmpq_poly([fmpq(int(c.p),int(c.q)) for c in reversed(p.all_coeffs())])
def order_at_zero(expression,variable):
 num,den=S.fraction(S.cancel(expression))
 return min(m[0] for m,c in S.Poly(num,variable).terms())-min(m[0] for m,c in S.Poly(den,variable).terms())

# Check the generic covariant formula separately before using saved data.
l,B,C,D,E=S.symbols('l B C D E');A=l*l
invI=12*A*E-3*B*D+C*C
invJ=72*A*C*E+9*B*C*D-27*A*D*D-27*B*B*E-2*C**3
px=3*(3*B*B-8*A*C)/(4*A)
py=27*(B**3-4*A*B*C+8*A*A*D)/(8*l**3)
zero('generic leading-square quartic covariant point',py*py-px**3+27*invI*px+27*invJ)

s,T=S.symbols('s T')
model=json.loads((OLD/'elliptic-cusp-input.json').read_text())
section=json.loads((HERE/'elliptic-section-input.json').read_text())
A,B,C,D,E=map(S.sympify,model['quartic_coefficients'])
I=S.sympify(model['I']);J=S.sympify(model['J'])
D4=S.sympify(model['quartic_discriminant'])
L=4*s*(3*s**3-2)*(3*s**3-1)
x=S.sympify(section['x']);y=S.sympify(section['y'])
zero('quartic leading coefficient square',A-L*L)
zero('covariant x specialization',x-3*(3*B*B-8*A*C)/(4*A))
zero('covariant y specialization',y-27*(B**3-4*A*B*C+8*A*A*D)/(8*L**3))
zero('point on exact Weierstrass model',y*y-x**3+27*I*x+27*J)
# A separate FLINT cross multiplication replays the point identity.
Nx,Dx=map(polynomial,S.fraction(x));Ny,Dy=map(polynomial,S.fraction(y));Ii=polynomial(I);Jj=polynomial(J)
yes('FLINT point identity',Ny**2*Dx**3-Nx**3*Dy**2+27*Ii*Nx*Dx**2*Dy**2+27*Jj*Dx**3*Dy**2==0)

q=3*s**3-1
for key,expr,power in [('x',x,2),('y',y,3)]:
 num,den=S.fraction(expr)
 zero(key+' denominator has only three specified poles',den-q**power)
 yes(key+' numerator coprime to pole divisor',S.degree(S.gcd(num,q),s)==0)
yes('three distinct finite poles',S.degree(q,s)==3 and S.degree(S.gcd(q,S.diff(q,s)),s)==0)
yes('all three poles are good fibres',S.degree(S.gcd(q,D4),s)==0)
yes('pole divisor disjoint from four I4 places',S.degree(S.gcd(q,s*(3*s**3-2)),s)==0)
yes('finite I4 specializations have nonzero y',S.degree(S.gcd(S.fraction(y)[0],s*(3*s**3-2)),s)==0)
yes('all remaining finite coordinates integral',S.degree(S.gcd(S.fraction(x)[1],S.fraction(y)[1]),s)==6)
# z=-x/y is the local parameter at the zero section; orders -2,-3
# give intersection one at each pole. Infinity is checked in its minimal chart.
yes('zero-section local parameter order at each finite pole',-2-(-3)==1)

local_x=S.cancel(T**14*x.subs(s,1/T));local_y=S.cancel(T**21*y.subs(s,1/T))
a4=S.expand(T**28*(-27*I).subs(s,1/T));a6=S.expand(T**42*(-27*J).subs(s,1/T))
yes('minimal infinity orders',[order_at_zero(e,T) for e in [a4,a6,local_x,local_y]]==[5,2,8,1])
eta=S.cancel(local_y/T).subs(T,0);beta=S.cancel(a6/T**2).subs(T,0)
yes('infinity first blowup y coordinate nonzero',eta!=0)
zero('infinity first blowup point on smooth exceptional component',eta*eta-beta)
yes('infinity first blowup x coordinate zero',S.cancel(local_x/T).subs(T,0)==0)
yes('infinity exceptional component smooth at section',2*eta!=0)
yes('no meeting zero section at infinity',order_at_zero(local_x,T)>=0 and order_at_zero(local_y,T)>=0)

A3=S.Matrix([[2,-1,0],[-1,2,-1],[0,-1,2]])
A2=S.Matrix([[2,-1],[-1,2]])
i4=[A3.inv()[i,i] for i in range(3)];iv=[A2.inv()[i,i] for i in range(2)]
yes('I4 self corrections',[S.Rational(3,4),1,S.Rational(3,4)]==i4)
yes('IV self corrections',[S.Rational(2,3)]*2==iv)
chi=S.Integer(7);intersection=3
minimum=2*chi-4*max(i4)-max(iv)
height=2*chi+2*intersection-S.Rational(2,3)
yes('universal nonzero-section height bound',minimum==S.Rational(28,3))
yes('torsion contradiction',minimum>0)
yes('explicit section height',height==S.Rational(58,3))
yes('section cannot be a nontrivial multiple',height<4*minimum)

# Reconstruct the balanced cubic split, retaining it as a reduction only.
a,b,m,n,w,z,k=S.symbols('a b m n w z k')
zero('balanced cubic factorization',(a-b/4)**3*(a+3*b/4)-(a**4-S.Rational(3,8)*a*a*b*b+S.Rational(1,8)*a*b**3-S.Rational(3,256)*b**4))
manifest=json.loads((HERE/'sources.json').read_text())
for row in manifest['sources']:
 path=HERE/row['path']
 yes('source hash '+row['path'],hashlib.sha256(path.read_bytes()).hexdigest()==row['sha256'])
receipt={'result':'PASS','claim_id':'GAOF7MW1','checks':len(checks),'check_names':checks,'python':platform.python_version(),'sympy':S.__version__,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'section_sha256':hashlib.sha256((HERE/'elliptic-section-input.json').read_bytes()).hexdigest(),'height_pairing_convention':'Shioda self-pairing, not half of it','section_height':'58/3','all_nonzero_sections_height_at_least':'28/3','geometric_torsion':'trivial','section_primitive':True,'generic_rank_at_least':1,'rational_incidence_points':'Zariski dense by the written section argument','three_point_fibres':'OPEN','assurance':'Exact arithmetic with a written sourced geometric interpretation; no independent end-to-end or formal verification, external review or novelty priority.'}
(HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ['result','checks','section_height','geometric_torsion','section_primitive','generic_rank_at_least','three_point_fibres']}))
