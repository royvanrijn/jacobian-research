"""Bounded exact F7 certificates. No finite-field search or Groebner run."""
import hashlib
import json
from pathlib import Path
import platform
import sympy as S
from sympy.polys.domains import QQ
from sympy.polys.rings import ring

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


q=S.Rational
t,a,b,g,u,v,w=S.symbols('t a b g u v w')
Y1,Y2,Y3,Y4=S.symbols('Y1 Y2 Y3 Y4')
YY=[Y1,Y2,Y3,Y4]
G2=a*a+a+t*b+t**3
G3=2*t*G2+a+t**3+t*t*b
G4=3*t*t*G2+b+t**4
X1=-S.Matrix([G2,G3,G4]).jacobian([t,a,b]).det().expand()
X=S.Matrix([X1,G2+t*X1,G3+t*t*X1,G4+t**3*X1])
SW=X+g*S.Matrix([1,t,t*t,t**3])
E=[S.cancel(f.subs({t:g*u,a:g**3*v,b:g**4*w})/g**i)
   for i,f in enumerate(SW,1)]
be=[S.Integer(1),u,u**3+u*u+v,u**4+u**3+w]
for i,f in enumerate(E):
    yes(f'scaled sweep coordinate {i} polynomial',S.denom(f)==1)
    zero(f'gamma-zero scaled coordinate {i}',f.subs(g,0)-be[i])

e,y,z1,z2,z3=S.symbols('e y z1 z2 z3')
stage={g:1-61*e*y+9012*e*e*y*y+q(238,19)*e**3*y*z1-q(35823670,19)*e**3*y**3+8*e**3*z2+e**4*z3,
       u:1+59*e*y+e*e*z1-5178*e*e*y*y+19*e**3*z2,
       v:-7*e*y-27*e*e*z1,w:-1-234*e*y-5*e*e*z1}
H=[y,-722*y*y+61*z1,
   (72851307*y**3-24938*y*z1-1957*z2)/19,
   (130524671288*y**4+1724330*y*y*z1-111055*y*z2-22572*z1*z1-361*z3)/19]
for i,f in enumerate(E,1):
    jet=S.series(f.subs(stage),e,0,i+1).removeO().expand()
    for j in range(i):zero(f'weight {j} below output {i}',jet.coeff(e,j))
    zero(f'x-zero output {i}',jet.coeff(e,i)-H[i-1])
yes('triangular boundary diagonal',[S.diff(f,z) for f,z in zip(H,[y,z1,z2,z3])]==[1,61,-103,-19])
yes('boundary determinant',61*103*19==119377)
v1,v2,v3,v4=S.symbols('v1 v2 v3 v4')
stageV=S.Matrix([1-61*v1+9012*v1*v1+q(238,19)*v1*v2-q(35823670,19)*v1**3+8*v3+v4,
                1+59*v1+v2-5178*v1*v1+19*v3,-7*v1-27*v2,-1-234*v1-5*v2])
zero('stage determinant',stageV.jacobian([v1,v2,v3,v4]).det()-119377)
iv1=(5*v-27*(w+1))/6283
iv2=(7*(w+1)-234*v)/6283
iv3=(u-1-59*iv1-iv2+5178*iv1*iv1)/19
iv4=g-1+61*iv1-9012*iv1*iv1-q(238,19)*iv1*iv2+q(35823670,19)*iv1**3-8*iv3
for i,(f,expected) in enumerate(zip(stageV,[g,u,v,w])):
    zero(f'stage inverse {i}',f.subs({v1:iv1,v2:iv2,v3:iv3,v4:iv4})-expected)
zero('sweep determinant',SW.jacobian([g,t,a,b]).det(method='domain-ge')-g)
x,A,B,D,Z=S.symbols('x A B D Z')
u0=B*x*x;v0=x**3*D-u0**3-u0*u0;w0=x**4*Z-u0**4-u0**3
for i,f in enumerate(be,1):
    residual=f.subs({u:u0,v:v0,w:w0})/x**i-[A,B,D,Z][i-1]
    zero(f'gamma-zero boundary reconstruction {i}',residual-(1/x-A if i==1 else 0))

bt=Y4-3*Y2*t*t+2*Y1*t**3-t**4
at=Y3-2*Y2*t+Y1*t*t-t**3-t*t*bt
R=S.expand(at*at+at+t*bt+t**3+Y1*t-Y2)
gt=S.expand(Y1-X1.subs({a:at,b:bt}))
quotient,remainder=S.div(gt-S.diff(R,t),R,t)
zero('Gao simple-root congruence',remainder)
for i,f in enumerate(SW):
    residual=f.subs({a:at,b:bt,g:gt},simultaneous=True)-YY[i]
    zero(f'full sweep reconstruction {i}',S.rem(S.expand(residual),R,t))
L=t*bt+t**3+Y1*t-Y2-q(1,4)
zero('square completion',R-(at+q(1,2))**2-L)
yes('square obstruction degree and leading coefficient',S.degree(L,t)==5 and S.expand(L).coeff(t,5)==-1)

H7=(t-1)*(t**6+t**5+t**4-t**3-t*t-t+1)
zero('known seven-point polynomial',R.subs(dict(zip(YY,[0,0,0,0])))-t**5*H7)
yes('H7 squarefree discriminant',S.discriminant(H7,t)==713)
yes('H7 constant term',H7.subs(t,0)==-1)
samples=[('twelve',[1,2,-1,3],12),('ten',[-1,0,0,1],10),
         ('nine',[0,0,-1,0],9),('seven',[0,0,0,0],7)]
for name,target,count in samples:
    f=S.Poly(R.subs(dict(zip(YY,target))),t)
    dg=S.degree(S.gcd(f.as_expr(),f.diff().as_expr()),t)
    dd=S.degree(S.gcd(S.gcd(f.as_expr(),f.diff().as_expr()),f.diff().diff().as_expr()),t)
    yes('rational target '+name,12-2*dg+dd==count)
f10=R.subs(dict(zip(YY,[-1,0,0,1])))
yes('simple-double discriminant witness',S.gcd(f10,S.diff(f10,t))==t)
yes('transverse discriminant witness',S.diff(R,Y2).subs(dict(zip(YY,[-1,0,0,1]))).subs(t,0)==-1)

# A separate exact number-field witness for a quadruple root and eight
# simple roots, avoiding numerical approximations to the quintic root.
z=S.symbols('z')
f5=16*z**5+48*z**4+48*z**3+10*z*z-13*z-5
yes('eight-point witness field degree five',S.Poly(f5,z).is_irreducible)
alpha=S.CRootOf(f5,0);K=QQ.algebraic_field(alpha);TR,tt=ring('t',K);aa=K.unit
bb=aa*(6*aa**3+12*aa*aa+7*aa+1)/(2*aa+1)
y8=[4*aa**3+6*aa*aa+2*aa-bb,aa*aa+aa,aa,bb]
b8=y8[3]-3*y8[1]*tt**2+2*y8[0]*tt**3-tt**4
a8=y8[2]-2*y8[1]*tt+y8[0]*tt**2-tt**3-tt**2*b8
r8=a8*a8+a8+tt*b8+tt**3+y8[0]*tt-y8[1]
yes('eight-point witness exact gcd',r8.gcd(r8.diff(tt))==tt**3)

# Independently differentiate the coefficient-factorization equations
# in SymPy. Preflight used FLINT column formulas; no search is repeated.
liftdata=json.loads((HERE/'lifting-points.json').read_text())
lift_receipts={}
for count,row in sorted(liftdata['points'].items()):
    p=row['modulus'];sub=dict(zip(YY,row['Y']));unknowns=list(YY);product=S.Integer(1);factors=[]
    for m,coeffs in sorted(row['multiplicity_factors'].items(),key=lambda z:int(z[0])):
        m=int(m);degree=len(coeffs)-1;variables=S.symbols(f'c{m}_0:{degree}')
        unknowns.extend(variables);sub.update(dict(zip(variables,coeffs[:-1])))
        factor=t**degree+sum(vv*t**j for j,vv in enumerate(variables))
        product*=factor**m
        fixed=S.Poly(sum(int(vv)*t**j for j,vv in enumerate(coeffs)),t,modulus=p)
        factors.append((m,fixed))
        yes(f'lift {count} factor {m} squarefree',S.gcd(fixed,fixed.diff()).degree()==0)
    for i,(_,f) in enumerate(factors):
        for j,(_,h) in enumerate(factors[:i]):yes(f'lift {count} factors {i},{j} disjoint',S.gcd(f,h).degree()==0)
    polys=S.Poly(S.expand(R-product),t)
    equations=S.Matrix([polys.coeff_monomial(t**j) for j in range(12)])
    yes(f'lift {count} exact factorization mod p',all(int(vv.subs(sub))%p==0 for vv in equations))
    matrix=equations.jacobian(unknowns).subs(sub)
    columns=row['pivot_columns'];minor=matrix[:,columns]
    det=int(minor.det())%p
    yes(f'lift {count} Jacobian minor',det!=0)
    simple_degree=sum(f.degree() for m,f in factors if m==1)
    yes(f'lift {count} exact simple-root degree',simple_degree==int(count))
    lift_receipts[count]={'modulus':p,'variables':len(unknowns),'equations':12,
                          'pivot_columns':columns,'minor_determinant_mod_p':det,
                          'fixed_columns':{str(i):int(sub[unknowns[i]]) for i in range(len(unknowns)) if i not in columns},
                          'multiplicity_factor_degrees':{str(m):f.degree() for m,f in factors}}

# Weighted Hensel calculation determines the saturated discriminant's
# boundary term without expanding the full four-parameter discriminant.
rho,h=S.symbols('rho h')
pulled=R.subs({t:rho*h,Y1:rho**4*A,Y2:rho**8*B,Y3:rho**12*D,Y4:rho**16*Z})
weighted=S.cancel(pulled/rho**5)
yes('weighted degeneration polynomial',S.denom(weighted)==1)
zero('five-root leading polynomial',weighted.subs(rho,0)+h**5-A*h)
zero('five-root leading discriminant',S.discriminant(h**5-A*h,h)+256*A**5)
yes('saturated boundary constant',-256*713==-182528)
yes('escaping roots different from one',S.gcd(H7/(t-1),t-1)==1)

def partitions(n,maximum=None):
    if n==0:
        yield []
        return
    for first in range(min(n,maximum or n),0,-1):
        for rest in partitions(n-first,first):yield [first]+rest

allparts=list(partitions(12))
emptytypes=[part for part in allparts if 1 not in part and any(m%2 for m in part)]
threetypes=[part for part in allparts if part.count(1)==3]
yes('all degree-twelve partitions counted',len(allparts)==77)
yes('ten remaining omitted-fibre types',len(emptytypes)==10)
yes('eight three-point types',len(threetypes)==8)
yes('omitted-fibre gcd threshold',min(sum(m-1 for m in part) for part in emptytypes)==7)
yes('three-point gcd threshold',min(sum(m-1 for m in part) for part in threetypes)==5)

manifest=json.loads((HERE/'sources/manifest.json').read_text())
for row in manifest['sources']:
    yes('source hash '+row['repository_path'],hashlib.sha256((ROOT/row['repository_path']).read_bytes()).hexdigest()==row['sha256'])
receipt={'result':'PASS','checks':len(checks),'check_names':checks,
         'claim_ids':['GAOF7BND1','GAOF7CARD1','GAOF7NP1'],
         'python':platform.python_version(),'sympy':S.__version__,
         'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'lifting_input_sha256':hashlib.sha256((HERE/'lifting-points.json').read_bytes()).hexdigest(),
         'lifting_certificates':lift_receipts,
         'proved_occurring_cardinalities':[1,2,4,5,6,7,8,9,10,12],
         'excluded_cardinality':[0,11],'unresolved_cardinalities':[3],
         'imported_exclusion':{'0':'GAOF7SURJ1; separate power-factorization proof and checker, not reverified in this checker'},
         'saturated_discriminant_c_valuation':5,'saturated_discriminant_boundary':'-182528*A^5',
         'hypothetical_empty_fibre_partitions_excluded_separately':emptytypes,'remaining_three_point_partitions':threetypes,
         'assurance':'Exact algebra and nonsingular finite-field lifting inputs plus written geometry; no independent end-to-end review, formal verification, external review or novelty priority.'}
(HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ['result','checks','proved_occurring_cardinalities','unresolved_cardinalities','lifting_certificates']}))
