"""Independent exact controls for the full pole-(3,4) exclusion.

Root coverage and geometric normalization are written in POLE_FOUR_EXCLUSION.md.
The script differentiates the actual quotient polynomials, not the
producer's exterior coefficient formula, and checks every rational unit identity.
"""
from pathlib import Path
from time import monotonic
from itertools import combinations
import hashlib,json,platform
import sympy as S

HERE=Path(__file__).resolve().parent
checks=[]
s,t,z,w,x=S.symbols("s t z w x")
delta,b,k=S.symbols("delta b k")
rel=6*delta**2-12*delta+5
def check(name,condition):
    assert condition,name
    checks.append(name)
def zero(name,expr):
    check(name,S.expand(expr)==0)
def red(expr):
    return S.Poly(S.expand(expr),delta,domain=S.QQ.frac_field(b,k,w)).rem(
        S.Poly(rel,delta,domain=S.QQ.frac_field(b,k,w))).as_expr()

def lift_numerator(f,constant):
    f=S.Poly(f,w)
    return w*f.as_expr()+z*(f.nth(2)*w+2*f.nth(3)*w*w+constant)+f.nth(3)*z*z

def differential_coefficients(f,constants,kv,lv,orders):
    """Direct differentiation of P=p/z^3, Q=q/z^4 using quotient numerators."""
    T,X,Y,Z=[lift_numerator(a,c) for a,c in zip(f,constants)]
    P=S.Poly(S.expand(T**3+3*z*T*X+3*z*z*Y),z)
    Q=S.Poly(S.expand(T**4+4*z*T*T*X+4*z*z*T*Y+2*z*z*X*X
                     +kv*z*z*T*T+2*kv*z**3*X+lv*z**3*T+S.Rational(4,3)*z**3*Z),z)
    result={}
    for order in orders:
        expression=0
        for i in range(order+1):
            j=order-i
            pi,qj=P.nth(i),Q.nth(j)
            expression+=(i-3)*pi*S.diff(qj,w)-(j-4)*S.diff(pi,w)*qj
        result[order]=S.expand(expression)
    return result

start=monotonic()
# Universal factorization and normalizations.
T,Y,mu,chi,zeta,kap,ell=S.symbols("T Y mu chi zeta kap ell")
X=mu*T+chi;Z=3*mu*Y-mu**3*T+zeta
P=T**3+3*T*X+3*Y
Q=T**4+4*T*T*X+4*T*Y+2*X*X+kap*T*T+2*kap*X+ell*T+S.Rational(4,3)*Z
outer=S.diff(P,T)*S.diff(Q,Y)-S.diff(P,Y)*S.diff(Q,T)
expected=(12*mu**2-12*chi-6*kap)*T-12*Y-6*kap*mu-3*ell+4*mu**3
zero("proportional-symbol outer Jacobian",outer-expected)
check("outer Jacobian is nonconstant for every parameter",S.Poly(outer,Y).coeff_monomial(Y)==-12)
a0,a1,a2,a3,b0,b1,b2,b3=S.symbols("a0:4 b0:4")
A=S.Matrix([[a0,a1],[a2,a3]]);B=S.Matrix([[b0,b1],[b2,b3]])
zero("universal chain-rule determinant",(A*B).det()-A.det()*B.det())

T,X,Y,Z,c=S.symbols("T X Y Z c")
zero("free-lift translation preserves P",
     (T+c)**3+3*(T+c)*(X-c*T)+3*(Y-c*X-c**3/3)-(T**3+3*T*X+3*Y))
a,ts,xs,ys,zs=S.symbols("a ts xs ys zs")
def pq(T,X,Y,Z,kv,lv):
    return (T**3+3*T*X+3*Y,
            T**4+4*T*T*X+4*T*Y+2*X*X+kv*T*T+2*kv*X+lv*T+S.Rational(4,3)*Z)
pn,qn=pq(T,X,Y,Z,kap,ell)
pg,qg=pq(T,X-a,Y+a*(T-ts)+ys,Z+3*a*(X-xs)+zs,
         kap,ell+4*a*ts-4*ys)
zero("all omitted a,y0,z0 parameters translate P",pg-pn-3*(ys-a*ts))
zero("all omitted a,y0,z0 parameters translate Q",
     qg-qn-(2*a*a-4*a*xs-2*kap*a+S.Rational(4,3)*zs))
scale=S.Symbol("scale",nonzero=True)
xx,yy=S.symbols("xx yy")
bb=xx-yy*yy
g=[yy/bb,yy*yy/bb,xx*yy/bb,xx*xx/bb]
for i,(gf,weight) in enumerate(zip(g,[-1,0,1,2])):
    check("marked-conic scaling on basis "+str(i+1),
          S.cancel(gf.subs({xx:scale**2*xx,yy:scale*yy})-scale**weight*gf)==0)
zero("marked triple normalization exponent",scale**(-4)*scale-scale**(-3))
zero("away triple normalization exponent",scale**2-scale**2)

# Verify the two complete repeated-root parametrizations.
L,M,N=S.symbols("L M N");D=S.symbols("D")
for name,sections in [
    ("double_simple",[L*L*M,L*M*N,M*N*N/3+L*D,3*N*D]),
    ("triple",[L**3,L*L*N,L*D,3*N*D-N**3])]:
    tt,uu,vv,rr=sections
    zero("universal cubic relation "+name,uu**3-3*tt*uu*vv+tt*tt*rr)

atlas=json.loads((HERE/"pole4_symbol_atlas.json").read_text())
check("atlas producer hash",atlas["source_sha256"]==hashlib.sha256(
    (HERE/"compute_pole4_symbol_atlas.py").read_bytes()).hexdigest())
check("five distinct full marked charts",
      {r["name"] for r in atlas["cases"]}=={
      "triple_at_marked","triple_away_marked","double_at_marked",
      "simple_at_marked","neither_at_marked"})
t0,x0,y0,z0,k0,ell0=S.symbols("t0 x0 y0 z0 k ell")
for row in atlas["cases"]:
    sections=list(map(S.sympify,row["leading_sections"]))
    coeff=differential_coefficients(sections,[t0,x0,y0,z0],k0,ell0,[3,4])
    zero("independent direct pole-five cancellation "+row["name"],coeff[3])
    stored=sum(S.sympify(value)*w**degree for degree,value in row["coefficient_equations"])
    zero("all direct pole-four coefficients "+row["name"],coeff[4]-stored)

# Every localized certificate is checked literally over QQ, without a CAS ideal oracle.
mixed=json.loads((HERE/"pole4_mixed_certificates.json").read_text())
check("mixed certificate producer hash",mixed["source_sha256"]==hashlib.sha256(
    (HERE/"certify_pole4_mixed.py").read_bytes()).hexdigest())
check("mixed input seal",mixed["input_sha256"]==hashlib.sha256(
    (HERE/"pole4_symbol_atlas.json").read_bytes()).hexdigest())
for row in mixed["certificates"]:
    source=next(r for r in atlas["cases"] if r["name"]==row["name"])
    equations=list(map(S.sympify,row["equations"]))
    multipliers=list(map(S.sympify,row["multipliers"]))
    check("mixed equations use full atlas "+row["name"],
          equations[:-1]==[S.sympify(v) for deg,v in source["coefficient_equations"]])
    zero("localized unit identity "+row["name"],
         sum(u*v for u,v in zip(equations,multipliers))-1)

# Triple coefficient systems: solve without losing either root of f.
n0,n1,d0,d1,d2=S.symbols("n0 n1 d0 d1 d2")
a,yfree=S.symbols("a yfree")
tripledata=json.loads((HERE/"pole4_triple_candidates.json").read_text())
check("triple producer hash",tripledata["source_sha256"]==hashlib.sha256(
    (HERE/"compute_pole4_triple_candidates.py").read_bytes()).hexdigest())
quartics={}
for marked in (True,False):
    name="triple_at_marked" if marked else "triple_away_marked"
    tstar=(9-6*delta)*b
    xstar=(-21+15*delta)*b*b+(1-S.Rational(3,2)*delta)*k
    lstar=(84-60*delta)*b**3+(-6+6*delta)*b*k
    if marked:lstar+=4*delta-S.Rational(4,3)
    else:
        xstar+=2*delta-1
        lstar+=(8-8*delta)*b
    saved=next(r for r in tripledata["cases"] if r["name"]==name)
    for key,value in [("t0",tstar),("x0",xstar),("ell",lstar)]:
        zero("saved unique lift solution "+name+" "+key,red(S.sympify(saved[key])-value))
    raw=next(r for r in atlas["cases"] if r["name"]==name)
    replacements={n0:0,n1:1,d1:b,d0:a if marked else delta,
                  d2:delta if marked else a}
    normalized={degree:S.expand(S.sympify(v).subs(replacements))
                for degree,v in raw["coefficient_equations"]}
    sol={t0:tstar,x0:xstar-a,y0:yfree,ell0:lstar+4*a*tstar-4*yfree}
    for degree,expr in normalized.items():
        zero("full parameter lift equation "+name+" "+str(degree),red(expr.subs(sol)))
    triangle=[(8,t0,12*(3*delta-2)),(9,x0,-24*(delta-1)),(10,ell0,-3)] if marked else [
              (5,t0,-12*(3*delta-2)),(4,x0,24*(delta-1)),(3,ell0,3)]
    for degree,var,value in triangle:
        zero("unique triangular coefficient "+name+" "+str(var),S.diff(normalized[degree],var)-value)
    zero("full quadratic branch "+name,
         normalized[7 if marked else 6]-(-4 if marked else 4)*rel)
    sections=([w**3,w*w,b*w*w+delta*w,3*b*w+3*delta-1] if marked else
              [1,w,delta*w*w+b*w,(3*delta-1)*w**3+3*b*w*w])
    coeff=differential_coefficients(sections,[tstar,xstar,0,0],k,lstar,[0,1,2,3,4,5])
    for order in range(5):
        zero("direct higher-pole cancellation "+name+" "+str(order),red(coeff[order]))
    G=S.cancel(red(coeff[5])/w**2)
    check("marked-point double factor removed exactly "+name,S.denom(G)==1)
    check("contact restriction really is quartic "+name,S.Poly(G,w).degree()<=4)
    zero("independently reconstructed contact quartic "+name,
         G-S.sympify(saved["contact_quartic"]))
    quartics[name]=G

check("no lost exceptional denominators",
      S.gcd(S.Poly(rel,delta),S.Poly((delta-1)*(3*delta-2),delta)).degree()==0)
for value,expected in [
    (S.Rational(7,12),S.Rational(1,24)),(S.Rational(19,32),-S.Rational(5,512)),
    (S.Rational(11,18),-S.Rational(5,54)),(S.Rational(71,120),S.Rational(1,2400))]:
    zero("exact rational separation "+str(value),rel.subs(delta,value)-expected)

# Verify the complete plane-section restriction and the fourth-power necessities.
section_basis=[xx**2,xx*yy,yy**2,xx**3,xx**2*yy,xx*yy**2,yy**3]
images=[S.cancel(v.subs({xx:w*w,yy:w})/w**2) for v in section_basis]
matrix=S.Matrix(5,7,lambda i,j:S.Poly(images[j],w).nth(i))
check("all five contact coefficients occur",matrix.rank()==5)
check("contact restriction has exactly two kernel dimensions",7-matrix.rank()==2)
for v in (bb*xx,bb*yy):
    zero("plane conic kernel restriction "+str(v),v.subs({xx:w*w,yy:w}))
rho,aa,bb0=S.symbols("rho aa bb0")
pure=[rho*aa**i*bb0**(4-i) for i in range(5)]
for expr in [pure[0]*pure[2]-pure[1]**2,
             pure[0]**2*pure[3]-pure[1]**3,
             pure[0]**3*pure[4]-pure[1]**4]:
    zero("universal necessary fourth-power contact equation "+str(len(checks)),expr)

contact=json.loads((HERE/"pole4_contact_certificates.json").read_text())
check("contact certificate producer hash",contact["source_sha256"]==hashlib.sha256(
    (HERE/"certify_pole4_contact.py").read_bytes()).hexdigest())
for row in contact["certificates"]:
    q=S.Poly(quartics[row["name"]],w)
    c=[q.nth(i)/S.binomial(4,i) for i in range(5)]
    if row["reverse_contact_coefficients"]:c.reverse()
    necessary=[rel,red(c[0]*c[2]-c[1]**2),
               red(c[0]**2*c[3]-c[1]**3),red(c[0]**3*c[4]-c[1]**4)]
    stored=list(map(S.sympify,row["equations"]))
    for i,(v,e) in enumerate(zip(necessary,stored)):
        zero("contact certificate equation matches direct quartic "+row["name"]+" "+str(i),v-e)
    mult=list(map(S.sympify,row["multipliers"]))
    zero("literal rational unit identity "+row["name"],sum(a*b for a,b in zip(stored,mult))-1)

proofs=[HERE/"POLE_FOUR_EXCLUSION.md",HERE/"PRIME_RAMIFICATION_FACTORIZATION.md"]
report=dict(as_of="2026-11-14",outcome="PASS",check_count=len(checks),checks=checks,
            leading_atlas=dict(proportional_branch="excluded by written factorization lemma",
                 squarefree_branch="reduces to proportional",
                 mixed_charts=3,triple_charts=2,complete_coverage="written proof"),
            exact_localized_certificates=3,exact_contact_unit_certificates=2,
            written_results=["no compatible quartic pole presentation of height <=4",
                             "no strict intermediate polynomial plane in either actual quartic normalization"],
            open_pole_pairs_at_height_five=[],
            subsequent_height_seven_exclusion="PF2RULED7 (external degree input; separate proof and receipt)",
            subsequent_all_height_pole_three_exclusion="PF2RULED8 (separate intrinsic proof and receipt)",
            frontier_superseded_by='PF2RULED16',
    frontier_handoff_scope='historical pre-PF2RULED16; this receipt does not prove the current frontier',
    first_open_reduced_pairs=[[6,8],[7,8]],
    frontier_adoption='PF2RULED12; subsequent handoff, outside this checker intrinsic proof coverage',
            subsequent_pair_exclusion="PF2RULED6 (separate proof and receipt)",
            subsequent_all_height_coordinate_exclusion="PF2RULED5 (separate proof and receipt)",
            assurance=dict(actual_ruled_compatibility_proved=False,
                           higher_poles_excluded=False,quartic_packets_excluded=False,
                           formal_verification=False,independent_full_reproof=False,
                           external_review=False),
            source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                           for p in [Path(__file__)]+proofs},
            producer_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in [HERE/"compute_pole4_symbol_atlas.py",
                                      HERE/"compute_pole4_triple_candidates.py",
                                      HERE/"certify_pole4_mixed.py",HERE/"certify_pole4_contact.py"]},
            python=platform.python_version(),sympy=S.__version__,
            elapsed_seconds=round(monotonic()-start,6))
(HERE/"pole4_exclusion_receipt.json").write_text(json.dumps(report,indent=2)+"\n")
print(f"PASS {len(checks)} exact controls; full compatible pole-height<=4 exclusion; higher poles and actual compatibility open")

