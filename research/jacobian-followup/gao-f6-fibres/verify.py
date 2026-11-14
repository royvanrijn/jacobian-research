"""Exact bounded checks for Gao F6 fibres. See PROOF.md for geometry."""
import hashlib
import json
from pathlib import Path
import platform
import sympy as S

if not __debug__:
    raise RuntimeError("Assertions require Python without -O")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
checks = []


def zero(name, expr):
    assert S.cancel(expr) == 0, name
    checks.append(name)


def yes(name, condition):
    assert condition, name
    checks.append(name)


q = S.Rational
s,a,b,g,u,v,w = S.symbols('s a b g u v w')
Y1,Y2,Y3,Y4 = S.symbols('Y1 Y2 Y3 Y4')
G2=a+s*b+s*s-s**3
G3=2*s*G2+a
psi=2*s**4-2*s**5
G4=b-G2**2+psi
X1=-(S.Matrix([G2,G3,G4]).jacobian([s,a,b])).det().expand()
zero('sweep first coordinate', X1-(2*a+2*b*s-b-10*s**5+8*s**4-2*s**3+5*s*s-2*s))
SW=S.Matrix([X1,G2+s*X1,G3+s*s*X1,G4+a*X1])+g*S.Matrix([1,s,s*s,a])
E=[S.cancel(f.subs({s:g*u,a:g**3*v,b:g**4*w})/g**i)
   for i,f in enumerate(SW,1)]
for i,f in enumerate(E,1):
    yes(f'E{i} polynomial', S.denom(f)==1)
boundaryE=[1-2*u,u-u*u,u*u+v,u**4-2*u*v+v+w]
for i,(f,be) in enumerate(zip(E,boundaryE),1):
    zero(f'gamma-zero stage {i}',f.subs(g,0)-be)

# Weighted jets prove all missing low-weight coefficients vanish. Since
# the unexpanded expressions are polynomials in v1..v4, this certifies
# global polynomiality without expanding the high-degree component list.
e,y,z1,z2,z3=S.symbols('e y z1 z2 z3')
stage={g:1-29*e*y+999*e**2*y*y+355*e**3*y*z1-41553*e**3*y**3+e**4*z3,
       u:1+27*e*y-5*e*e*z1,
       v:e*y-12*e*e*z1+q(2128,5)*e*e*y*y+e**3*z2,
       w:-4*e*y-10*e*e*z1}
H=[]
expected=[y,29*z1-q(4636,5)*y*y,
          (6039*y**3-5205*y*z1+25*z2)/5,
          (74000249*y**4-463585*y*y*z1+125*y*z2-26025*z1*z1-50*z3)/25]
for i,f in enumerate(E,1):
    jet=S.series(f.subs(stage),e,0,i+1).removeO().expand()
    for j in range(i):zero(f'no source weight {j} below output weight {i}',jet.coeff(e,j))
    H.append(jet.coeff(e,i))
    zero(f'boundary x-zero output {i}',H[-1]-expected[i-1])
yes('boundary triangular diagonal',[S.diff(f,t) for f,t in zip(H,[y,z1,z2,z3])]==[1,29,5,-2])

v1,v2,v3,v4=S.symbols('v1 v2 v3 v4')
stage_v=S.Matrix([1-29*v1+999*v1*v1+355*v1*v2-41553*v1**3+v4,
                  1+27*v1-5*v2,v1-12*v2+q(2128,5)*v1*v1+v3,-4*v1-10*v2])
zero('stage automorphism determinant',stage_v.jacobian([v1,v2,v3,v4]).det()+290)
inv1=(u-1-w/2)/29
inv2=-w/10-2*inv1/5
inv3=v-inv1+12*inv2-q(2128,5)*inv1*inv1
inv4=g-1+29*inv1-999*inv1*inv1-355*inv1*inv2+41553*inv1**3
for i,(f,want) in enumerate(zip(stage_v,[g,u,v,w])):
    zero(f'stage inverse {i}',f.subs({v1:inv1,v2:inv2,v3:inv3,v4:inv4})-want)
# In row order (gamma,s,a,b) the sweep has determinant gamma.
zero('sweep Jacobian',SW.jacobian([g,s,a,b]).det(method='domain-ge')-g)

x,A,B,D,Z=S.symbols('x A B D Z')
u0=(1-x*A)/2
v0=x**3*D-u0*u0
w0=x**4*Z-u0**4+(2*u0-1)*v0
delta=A*A+4*B
for i,be in enumerate(boundaryE):
    residual=be.subs({u:u0,v:v0,w:w0})/x**(i+1)-[A,B,D,Z][i]
    zero(f'gamma-zero target reconstruction {i}',residual-((1-delta*x*x)/(4*x*x) if i==1 else 0))

astar=Y3-2*s*Y2+s*s*Y1
bstar=Y4+Y2*Y2-Y1*Y3-psi
gstar=Y1-X1.subs({a:astar,b:bstar})
R=2*s**6-2*s**5-s**3+(Y1+1)*s*s+(Y4+Y2*Y2-Y1*Y3-2*Y2+Y1)*s+Y3-Y2
zero('Gao lifting derivative identity',gstar-S.diff(R,s)+2*R)
for i,f in enumerate(SW):
    residual=f.subs({a:astar,b:bstar,g:gstar},simultaneous=True)-[Y1,Y2,Y3,Y4][i]
    zero(f'sweep lifting residual {i}',S.rem(S.together(residual).expand(),R,s))
A0,B0,C0,h=S.symbols('A0 B0 C0 h')
P=s**6-s**5-s**3/2+A0*s*s+B0*s+C0
coeff={A0:(Y1+1)/2,B0:(Y4+Y2*Y2-Y1*Y3-2*Y2+Y1)/2,C0:(Y3-Y2)/2}
zero('monic coefficient family',P.subs(coeff)-R/2)
invY={Y1:2*A0-1,Y2:h,Y3:h+2*C0,Y4:2*B0-h*h+(2*A0-1)*(h+2*C0)+2*h-(2*A0-1)}
for name,expr,want in [('A0',coeff[A0],A0),('B0',coeff[B0],B0),('C0',coeff[C0],C0)]:
    zero('coefficient coordinate inverse '+name,expr.subs(invY)-want)

# The only partitions with no simple root have already been considered
# by Gao. The forced square is retained and its three roots are distinct.
Q=s**3-s*s/2-s/8-q(5,16)
square=S.expand(Q*Q)
for j,want in [(5,-1),(4,0),(3,-q(1,2))]:zero(f'forced square coefficient {j}',square.coeff(s,j)-want)
yes('forced cubic has three distinct roots',S.discriminant(Q,s)!=0)
zero('two-triple obstruction',S.expand((s*s-s/3-q(1,9))**3).coeff(s,3)-q(5,27))
r=S.symbols('r')
five_one=S.expand((s-r)**5*(s-(1-5*r)))
yes('five-one first constraints',S.factor(five_one.coeff(s,4))== -5*r*(3*r-1))
for root in [0,q(1,3)]:yes('five-one third constraint '+str(root),five_one.coeff(s,3).subs(r,root)!=-q(1,2))
yes('sixfold root impossible',S.expand((s-q(1,6))**6).coeff(s,4)!=0)

K4=40*r**3-20*r*r-1
K32=4320*r**6-4320*r**5+540*r**4-448*r**3+384*r*r+36*r+35
T3=s**3+(3*r-1)*s*s+(6*r*r-3*r)*s+10*r**3-6*r*r-q(1,2)
triple=S.expand((s-r)**3*T3)
for j,want in [(5,-1),(4,0),(3,-q(1,2))]:zero(f'triple family coefficient {j}',triple.coeff(s,j)-want)
zero('quadruple root condition',T3.subs(s,r)-K4/2)
zero('triple-double condition',S.discriminant(T3,s)+K32/4)
yes('three quadruple points',S.discriminant(K4,r)==-75200)
yes('six triple-double points',S.discriminant(K32,r)==-21557457167069587580544614400)
yes('triple degenerations disjoint',S.gcd(K4,K32)==1)

sig=S.symbols('sig')
pi=(8*sig**3-6*sig*sig-1)/(12*sig-4)
Q2=s*s-sig*s+pi
T2=s*s+(2*sig-1)*s+3*sig*sig-2*sig-2*pi
two_double=S.Poly(S.expand(Q2**2*T2),s)
for j,want in [(5,-1),(4,0),(3,-q(1,2))]:zero(f'two-double family coefficient {j}',two_double.coeff_monomial(s**j)-want)
J4=5*sig**3-5*sig*sig-1
J222=8*sig**3-8*sig*sig+sig+3
J321=1440*sig**6-2880*sig**5+2220*sig**4-728*sig**3+48*sig*sig+24*sig+13
zero('double roots collide',S.discriminant(Q2,s)+J4/(3*sig-1))
zero('residual roots collide',S.discriminant(T2,s)+J222/(3*sig-1))
zero('double meets simple root',S.resultant(Q2,T2,s)-J321/(16*(3*sig-1)**2))
for name,f in [('quadruple',J4),('three-double',J222),('triple-double',J321)]:
    yes('two-double '+name+' divisor squarefree',S.discriminant(f,sig)!=0)
divisors=[3*sig-1,J4,J222,J321]
for i,f in enumerate(divisors):
    for j,g0 in enumerate(divisors[:i]):yes(f'two-double punctures disjoint {i},{j}',S.gcd(f,g0)==1)

parts=[[1,1,1,1,1,1],[2,1,1,1,1],[2,2,1,1],[2,2,2],
       [3,1,1,1],[3,2,1],[3,3],[4,1,1],[4,2],[5,1],[6]]
excluded=[[3,3],[4,2],[5,1],[6]]
table=[]
for part in parts:
    test=S.prod((s-i)**m for i,m in enumerate(part))
    dg=S.degree(S.gcd(test,S.diff(test,s)),s)
    dgg=S.degree(S.gcd(S.gcd(test,S.diff(test,s)),S.diff(test,s,2)),s)
    yes('partition simple-root count '+str(part),6-2*dg+dgg==part.count(1))
    table.append({'partition':part,'gcd_degrees':[int(dg),int(dgg)],'fibre_points':part.count(1),'occurs':part not in excluded})
samples=[('six',2*s*(s+1)*(s-q(1,2))*(s-q(3,2))*(s*s+q(5,4)),6),
         ('four',s*s*(s-1)*(2*s**3-1),4),
         ('three',2*triple.subs(r,0),3),
         ('two',2*two_double.as_expr().subs(sig,0),2),
         ('zero',2*square,0)]
for name,f,count in samples:
    monic=S.Poly(f,s).monic().as_expr()
    for j,want in [(5,-1),(4,0),(3,-q(1,2))]:zero(f'witness {name} coefficient {j}',monic.coeff(s,j)-want)
    dg=S.degree(S.gcd(f,S.diff(f,s)),s)
    dgg=S.degree(S.gcd(S.gcd(f,S.diff(f,s)),S.diff(f,s,2)),s)
    yes('witness count '+name,6-2*dg+dgg==count)

# Discriminant normalization and saturated pullback. Target coordinates
# (c,A,k,l,q) are Gao's polynomial target change; q is an identity factor.
c,k,l=S.symbols('c k l')
pulled=P.subs({A0:(1+c*A)/2,B0:(c*A+2*c*c*k+c**3*l)/2,C0:c*c*k/2})
disc=S.Poly(S.discriminant(pulled,s),c,A,k,l)
yes('discriminant c-valuation two',min(m[0] for m in disc.monoms())==2)
sat=S.cancel(disc.as_expr()/c**2)
print_boundary=S.factor(sat.subs(c,0))
yes('boundary discriminant only quadratic',S.cancel(print_boundary/(A*A-4*k)).is_Rational)
yes('boundary discriminant nonzero',S.cancel(print_boundary/(A*A-4*k))!=0)
base=s**6-s**5-s**3/2
BB=-S.diff(base,s).subs(s,r)-2*A0*r
CC=r*S.diff(base,s).subs(s,r)-base.subs(s,r)+A0*r*r
zero('double incidence equation',P.subs({s:r,B0:BB,C0:CC}))
zero('double incidence derivative',S.diff(P,s).subs({s:r,B0:BB,C0:CC}))
discall=S.discriminant(P,s)
yes('discriminant reduced',S.diff(discall,C0).subs({A0:q(1,2),B0:0,C0:0})!=0)
Pzero=s*s*(s-1)*(s**3-q(1,2))
yes('boundary escape roots simple',S.gcd(S.diff(Pzero,s),2*s**3-1)==1)
yes('boundary escape roots not one',S.gcd(s-1,2*s**3-1)==1)

manifest=json.loads((HERE/'sources/manifest.json').read_text())
for row in manifest['sources']:
    yes('source hash '+row['repository_path'],hashlib.sha256((ROOT/row['repository_path']).read_bytes()).hexdigest()==row['sha256'])
receipt={'result':'PASS','claim_ids':['GAOF6FIB1','GAOF6NP1'],
         'checks':len(checks),'check_names':checks,
         'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'python':platform.python_version(),'sympy':S.__version__,
         'partition_table':table,'cardinalities':[0,1,2,3,4,6],
         'boundary_original_discriminant':'A^2+4*B',
         'boundary_transformed_discriminant':'A^2-4*k',
         'saturated_discriminant_boundary':str(print_boundary),
         'source_manifest_sha256':hashlib.sha256((HERE/'sources/manifest.json').read_bytes()).hexdigest(),
         'assurance':'Exact identities plus written geometry, not formal verification, external review or novelty priority.'}
(HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ['result','checks','cardinalities','saturated_discriminant_boundary']}))
