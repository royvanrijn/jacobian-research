"""Bounded exact checks for the complete F4 fibre classification; no searches.

The geometric proof is PROOF.md. This checker verifies its polynomial
identities and nonvanishing tests, not all geometric implications.
"""
import hashlib
import json
from pathlib import Path
import platform
import sympy as S

if not __debug__:
    raise RuntimeError("Run without -O: assertions are required")

HERE = Path(__file__).resolve().parent
checks = []


def zero(name, expression):
    assert S.cancel(expression) == 0, name
    checks.append(name)


def yes(name, condition):
    assert condition, name
    checks.append(name)


x, y, z, t = S.symbols("x y z t")
s, a, g, u, v, c = S.symbols("s a g u v c")
Y1, Y2, Y3 = S.symbols("Y1 Y2 Y3")
q = S.Rational
stage_g = 1 + x*y + x*x*t
stage_u = 1 + x*z
stage_v = -q(7, 9) + q(22, 9)*x*y + q(8, 3)*x*z
A = S.Matrix([s/3-4*s*s/3, 2*s*s/3-8*s**3/9,
              2*s+7*s**3/9-2*s**4/3])
B = S.Matrix([s*s-1, s**3, s**4+s*s+2])
direction = S.Matrix([1, s, s*s])
sweep = A+a*B+g*direction
numerators = [1+u/3-v-q(4, 3)*g*u*u+g*g*v*u*u,
              u+q(2, 3)*u*u-q(8, 9)*g*u**3+g*g*v*u**3,
              2*u+2*v+g*g*(u*u+q(7, 9)*u**3+v*u*u)
              -q(2, 3)*g**3*u**4+g**4*v*u**4]
for i, (component, num, power) in enumerate(zip(sweep, numerators, [1, 2, 1])):
    zero(f"numerator agrees with Gao sweep {i}",
         component.subs({s:g*u,a:g*v})-g**power*num)
F = [stage_g*x]
for i, (num, power) in enumerate(zip(numerators, [1, 2, 1]), 2):
    poly = S.Poly(S.expand(num.subs({g: stage_g, u: stage_u, v: stage_v})),
                  x, y, z, t, domain=S.QQ)
    yes(f"F{i} numerator divisible by x^{power}",
        all(m[0] >= power for m in poly.monoms()))
    terms = {(m[0]-power, *m[1:]): coef for m, coef in poly.terms()}
    F.append(S.Poly.from_dict(terms, (x, y, z, t)).as_expr())
yes("component degrees", [S.Poly(f, x, y, z, t).total_degree() for f in F]
    == [4, 11, 12, 21])
stage_det = S.Matrix([x, stage_g, stage_u, stage_v]).jacobian([x,y,z,t]).det()
zero("source stage determinant", stage_det+q(22,9)*x**4)
zero("rescaling determinant", S.Matrix([x,g,g*u,g*v]).jacobian([x,g,u,v]).det()-g*g)
padded = S.Matrix([g*x, *sweep])
zero("padded sweep determinant", padded.jacobian([x,g,s,a]).det()-2*g*g)
C, T1, T2, T3 = S.symbols("C T1 T2 T3")
zero("target twist determinant", S.Matrix([C,T1/C,T2/C**2,T3/C]).jacobian([C,T1,T2,T3]).det()-C**-4)
zero("chain Jacobian constant", stage_det*g*g*2*g*g/(g*x)**4+q(44,9))

# Boundary x=0: the polynomial extension, independently expanded above.
boundary0 = [S.Integer(0), -q(26,9)*y-q(35,9)*z,
             q(37,9)*y*y+q(16,3)*y*z+q(11,3)*z*z-q(22,9)*t,
             q(20,3)*y+q(29,3)*z]
for i, f in enumerate(F):
    zero(f"x-zero restriction {i}", f.subs(x,0)-boundary0[i])
av, bv, dv = S.symbols("av bv dv")
y0 = -(87*av+35*dv)/18
z0 = (30*av+13*dv)/9
t0 = q(9,22)*(q(37,9)*y0*y0+q(16,3)*y0*z0+q(11,3)*z0*z0-bv)
for i, f in enumerate(boundary0):
    zero(f"x-zero inverse {i}", f.subs({y:y0,z:z0,t:t0})-[0,av,bv,dv][i])

# The other boundary piece has gamma=0 and x nonzero. Work in exact
# rational functions of x, then reduce its one quadratic relation.
K = 2*av+dv
boundary_delta = 3*K*K-32*bv
u0 = q(3,8)*K*x-q(3,4)
v0 = dv*x/2-u0
zero("gamma-zero second target", numerators[0].subs({g:0,u:u0,v:v0})/x-av)
zero("gamma-zero fourth target", numerators[2].subs({g:0,u:u0,v:v0})/x-dv)
zero("gamma-zero third target residual",
     numerators[1].subs({g:0,u:u0,v:v0})/x**2-bv
     -(boundary_delta*x*x-12)/(32*x*x))
zero("quadratic nonzero roots product", S.Poly(boundary_delta*x*x-12,x).TC()+12)

# Universal stage inversion, including gamma=0 when x is retained.
xy = q(9,22)*(v+q(7,9)-q(8,3)*(u-1))
inverse_stage = {y:xy/x,z:(u-1)/x,t:(g-1-xy)/x**2}
for name, expression, expected in [("gamma",stage_g,g),("u",stage_u,u),("v",stage_v,v)]:
    zero("stage inverse "+name, expression.subs(inverse_stage)-expected)

phi = 2*s-q(2,9)*s**3-q(2,9)*s**4
astar = (Y3-2*s*Y2+s*s*Y1-phi)/2
gstar = Y1-A[0]-astar*B[0]
R = q(2,9)*s**5+q(2,9)*s**4+(q(8,9)+Y1)*s**3 \
    -(q(4,3)+2*Y2)*s*s+(2*Y1+Y3)*s-2*Y2
zero("Gao simple-root derivative identity", S.diff(R,s)-2*gstar-s*R)
residuals = [S.expand(f.subs({a:astar,g:gstar})-Y)
             for f,Y in zip(sweep,[Y1,Y2,Y3])]
zero("first sweep equation", residuals[0])
zero("second sweep equation", 2*residuals[1]-R)
zero("third sweep equation", residuals[2]-s*R)
A0, B0, D0 = S.symbols("A0 B0 D0")
P = s**5+s**4+A0*s**3+B0*s*s+D0*s+B0+6
zero("monic family", P.subs({A0:4+q(9,2)*Y1,B0:-6-9*Y2,D0:9*Y1+q(9,2)*Y3})-q(9,2)*R)
zero("quintuple-root impossibility", ((s+q(1,5))**5).expand().coeff(s,0)
     -((s+q(1,5))**5).expand().coeff(s,2)-6+q(18999,3125))

r, b, sig, pi, W = S.symbols("r b sig pi W")
d, w = S.symbols("d w")
single_double = S.expand((s-r)**2*(s**3+(1+2*r)*s*s+b*s+d))
H2 = (r*r-1)*d+2*r*b-r*r*(1+2*r)-6
zero("double incidence coefficient constraint",single_double.coeff(s,0)-single_double.coeff(s,2)-6-H2)
Tr = r*r*(1+2*r)+6
bd = {b:r*Tr/2+(r*r-1)*w,d:-Tr-2*r*w}
zero("double incidence affine plane",H2.subs(bd))
zero("double incidence inverse parameter", -(bd[b]-r*Tr/2)
     -r*(bd[d]+Tr)/2-w)
den = r*(r*r-3)
H4 = 4*r**5+r**4-20*r**3-6*r*r-6
h = 3*r**5+2*r**4+2*r**3-2*r*r-r+8
triple = S.expand((s-r)**3*(s*s+(1+3*r)*s+b))
H3 = -b*r**3+3*b*r-8*r**3-3*r*r-6
zero("triple family coefficient constraint", triple.coeff(s,0)-triple.coeff(s,2)-6-H3)
br = -(8*r**3+3*r*r+6)/den
zero("triple family solved", H3.subs(b,br))
yes("triple excluded denominators", S.gcd(8*r**3+3*r*r+6,den)==1)
Qr = s*s+(1+3*r)*s+br
zero("quadruple-root divisor", Qr.subs(s,r)-H4/den)
zero("triple-double divisor", S.discriminant(Qr,s)-3*h/den)
yes("five distinct quadruple points", S.discriminant(H4,r)==-43861243392)
yes("Gao five triple-double points", S.discriminant(h,r)==1292134656)
for name, f1, f2 in [("quadruple versus triple-double",H4,h),
                     ("quadruple versus denominator",H4,den),
                     ("triple-double versus denominator",h,den)]:
    yes(name+" disjoint", S.gcd(f1,f2)==1)

double = S.expand((s*s-sig*s+pi)**2*(s+1+2*sig))
H22 = (1+2*sig)*(pi*pi-sig*sig)-2*(1+sig)*pi-6
zero("two-double constraint", double.coeff(s,0)-double.coeff(s,2)-6-H22)
coeff_double = [-3*sig*sig-2*sig+2*pi,
                sig*sig*(1+2*sig)+2*pi*(1+sig),
                pi*pi-2*sig*pi*(1+2*sig)]
for i, (power, coefficient) in enumerate(zip([3,2,1],coeff_double)):
    zero(f"two-double coefficient {i}", double.coeff(s,power)-coefficient)
quartic = 4*sig**4+4*sig**3+2*sig*sig+14*sig+7
zero("elliptic square completion", ((1+2*sig)*pi-(1+sig))**2-quartic-(1+2*sig)*H22)
yes("elliptic quartic squarefree", S.discriminant(quartic,sig)==-12854016)
zero("exceptional affine point", H22.subs({sig:-q(1,2),pi:-6}))
zero("exceptional affine point smooth", S.diff(H22,pi).subs({sig:-q(1,2),pi:-6})+1)
T = S.symbols("T")
projective = (T+2*sig)*(pi*pi-sig*sig)-2*(T+sig)*pi*T-6*T**3
zero("cubic dehomogenization", projective.subs(T,1)-H22)
for point in [(0,1,0),(1,1,0),(1,-1,0)]:
    sub = dict(zip([sig,pi,T],point))
    zero("point at infinity "+str(point),projective.subs(sub))
    yes("smooth at infinity "+str(point),any(S.diff(projective,j).subs(sub)!=0 for j in [sig,pi,T]))
zero("all infinity points",projective.subs(T,0)-2*sig*(pi-sig)*(pi+sig))
zero("quadruple punctures on cubic", H22.subs({sig:2*r,pi:r*r})-H4)
other = (-1-3*r)/2
zero("triple-double punctures on cubic", H22.subs({sig:r+other,pi:r*other})+q(3,4)*h)

partitions = [[1,1,1,1,1],[2,1,1,1],[2,2,1],[3,1,1],[3,2],[4,1],[5]]
partition_table = []
for parts in partitions:
    sample = S.prod((s-i)**m for i,m in enumerate(parts))
    ga = S.degree(S.gcd(sample,S.diff(sample,s)),s)
    gb = S.degree(S.gcd(S.gcd(sample,S.diff(sample,s)),S.diff(sample,s,2)),s)
    simple = parts.count(1)
    yes("gcd cardinality "+str(parts),5-2*ga+gb==simple)
    partition_table.append({"partition":parts,"gcd_degrees":[int(ga),int(gb)],
                            "fibre_points":simple,"occurs":parts!=[5]})

# Rational targets witnessing all positive cardinalities, and explicit
# source checks at the new boundary and two-double stratum.
witnesses = {
    "five":s**5+s**4+4*s**3-15*s*s-9,
    "three":s*s*(s-1)*(s*s+2*s+6),
    "two":(s-1)**3*(s*s+4*s+q(17,2)),
    "one":s*(s*s+s/2-6)**2,
}
for name, expected in [("five",5),("three",3),("two",2),("one",1)]:
    poly = S.expand(witnesses[name])
    zero("witness family "+name,poly.coeff(s,0)-poly.coeff(s,2)-6)
    ga = S.degree(S.gcd(poly,S.diff(poly,s)),s)
    gb = S.degree(S.gcd(S.gcd(poly,S.diff(poly,s)),S.diff(poly,s,2)),s)
    yes("witness count "+name,5-2*ga+gb==expected)
source_witnesses = [
    ((0,0,0,0),(0,0,0,0)),
    ((2,q(10,11),-q(1,2),-q(31,44)),(0,0,0,1)),
    ((q(1,4),q(383,44),-4,q(145,11)),(1,-q(7,2),0,15)),
]
for j,(point,target) in enumerate(source_witnesses):
    for i,f in enumerate(F):
        zero(f"expanded source witness {j} component {i}",f.subs(dict(zip([x,y,z,t],point)))-target[i])

# Discriminant pullback and escape at every point of the target boundary.
pulled_P = S.expand(P.subs({A0:4+q(9,2)*c*av,B0:-6-9*c*c*bv,
                            D0:9*c*av+q(9,2)*c*dv}))
pulled_disc = S.Poly(S.discriminant(pulled_P,s),c,av,bv,dv)
yes("discriminant c valuation exactly two",min(m[0] for m in pulled_disc.monoms())==2)
E = S.cancel(pulled_disc.as_expr()/c**2)
yes("remaining discriminant polynomial",S.denom(E)==1)
zero("boundary intersection of discriminant",E.subs(c,0)+393660*boundary_delta)
discY = S.discriminant(q(9,2)*R,s)
yes("discriminant reduced at the single-double witness",
    S.diff(discY,Y2).subs({Y1:0,Y2:0,Y3:0})!=0)
Pzero = witnesses['three']
qzero = s*s+2*s+6
yes("boundary escape roots simple",S.gcd(S.diff(Pzero,s),qzero)==1)
zero("boundary escape numerator",S.rem(s-S.diff(Pzero,s)/9,qzero,s)+3*s+4)
yes("boundary escape numerator nonzero",S.gcd(3*s+4,qzero)==1)

# Bind all retained primary-source inputs as well as the checker itself.
manifest_path = HERE/'sources/manifest.json'
manifest = json.loads(manifest_path.read_text())
for row in manifest['sources']:
    yes("source SHA256 "+row['file'],
        hashlib.sha256((HERE/'sources'/row['file']).read_bytes()).hexdigest()==row['sha256'])
receipt = {
    "claim_ids":["GAOF4FIB1","GAOF4ELL1","GAOF4NP1"],
    "result":"PASS", "checks":len(checks), "check_names":checks,
    "python":platform.python_version(),"sympy":S.__version__,
    "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "source_files":{str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in sorted((HERE/'sources').iterdir()) if p.is_file()},
    "partition_table":partition_table,
    "possible_fibre_cardinalities":[0,1,2,3,5],
    "boundary_discriminant":"3*(2*v2+v4)^2-32*v3",
    "elliptic_completion_punctures":{"infinity":3,"quadruple":5,"triple_double":5},
    "assurance":"Exact identities and nonvanishing tests plus a separate written geometric proof; no formal verification, independent end-to-end review or novelty priority."
}
(HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({"result":"PASS","checks":len(checks),"cardinalities":[0,1,2,3,5]}))
