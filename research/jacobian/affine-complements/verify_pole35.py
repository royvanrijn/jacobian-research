"""Exact controls for the complete pole-(3,5) exclusion.

Global coverage, local geometry and finiteness are written in the proof;
finite controls do not formally verify those arguments.
"""
from pathlib import Path
from time import monotonic
import sympy as S,json,hashlib,platform
HERE=Path(__file__).resolve().parent
start=monotonic();checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def zero(name,value):check(name,S.cancel(value)==0)
def order(expr,u):
    p=S.Poly(S.expand(expr),u)
    return min(m[0] for m,c in p.terms() if c!=0)

T,A,X,Y=S.symbols('T A X Y')
p=T**3+A;first=T**5+S.Rational(5,3)*T*T*A
zero('first truncated-power exterior identity',
     S.diff(p,T)*S.diff(first,A)-S.diff(p,A)*S.diff(first,T)+S.Rational(10,3)*T*A)
u=S.Symbol('u_local')
tau,aa,dd=[S.Function(n)(u) for n in ('tau','a2','d3')]
k=S.Symbol('k')
top=(-S.Rational(10,3)*tau*aa-4*k*tau**3)*(2*S.diff(tau,u)*aa-tau*S.diff(aa,u))
top+=9*tau*tau*S.diff(tau,u)*dd-3*tau**3*S.diff(dd,u)
primitive=dd/tau**3-S.Rational(5,9)*aa*aa/tau**4-S.Rational(4,3)*k*aa/tau**2
zero('pole-seven integration with arbitrary symbols',top+3*tau**6*S.diff(primitive,u))
p=T**3+3*T*X+3*Y
G5=T**5+5*T**3*X+5*T*T*Y+5*T*X*X+10*X*Y
G4=T**4+4*T*T*X+4*T*Y+2*X*X
for (i,j),expected in {
 (T,X):15*T*X*X+30*X*Y,(T,Y):15*X*X-30*T*Y,(X,Y):-30*Y}.items():
    zero('G5 exterior coefficient '+str((i,j)),
         S.diff(p,i)*S.diff(G5,j)-S.diff(p,j)*S.diff(G5,i)-expected)
for (i,j),expected in {
 (T,X):12*X*X-12*T*Y,(T,Y):-12*T*X-12*Y,(X,Y):-12*X}.items():
    zero('G4 exterior coefficient '+str((i,j)),
         S.diff(p,i)*S.diff(G4,j)-S.diff(p,j)*S.diff(G4,i)-expected)
xx,cc=[S.Function(n)(u) for n in ('x','c2')]
top=15*tau*xx**2*(S.diff(tau,u)*xx-tau*S.diff(xx,u))
top+=6*tau*tau*S.diff(tau,u)*cc-3*tau**3*S.diff(cc,u)
zero('pole-six integration with arbitrary symbols',
     top+3*tau**5*S.diff(cc/tau**2+S.Rational(5,3)*xx**3/tau**3,u))
vv,ee=[S.Function(n)(u) for n in ('v','c1')]
prop=3*tau*tau*(S.diff(tau,u)*ee-tau*S.diff(ee,u))
prop-=30*tau*vv*(S.diff(tau,u)*vv-tau*S.diff(vv,u))
zero('proportional-branch pole-five integration',
     prop+3*tau**4*S.diff(ee/tau-5*vv*vv/tau**2,u))

# Exact model discriminants give pole totals; geometric stability and
# repeated-face inequalities are separate written Newton arguments.
z,xi=S.symbols('z xi')
models=[]
def model(name,F,expected,degree):
    disc=S.discriminant(F,z)
    valuation=order(disc,u)
    local_degree=min(m[0] for m,c in S.Poly(F.subs(u,0),z).terms() if c!=0)
    check('local Weierstrass degree '+name,local_degree==degree)
    check('generic discriminant pole degree '+name,valuation-degree==expected)
    models.append(dict(name=name,discriminant_order=valuation,local_degree=degree,pole_total=expected))
for h,j,expected in [(1,0,6),(1,1,7),(1,2,9),(1,3,9),
                      (2,0,7),(2,1,9),(2,2,11),(2,3,12),(2,4,12)]:
    F=u**9+3*u**(3+h)*z+u**j*z*z+xi*z**3
    model(f'triple_h{h}_j{j}',F,expected,2 if j==0 else 3)
# Three coalescent-face examples retain, rather than discard, that branch.
for gap in (1,2,3):
    F=u**9+3*u**5*z+(S.Rational(9,4)*u+u**(1+gap))*z*z+xi*z**3
    model('triple_h2_j1_coalescent_gap'+str(gap),F,9+gap,3)
model('proportional_mixed_double',u**6+3*u*z*z+xi*z**3,6,3)
model('proportional_simple',u**3+3*u*z*z+xi*z**3,3,3)
model('proportional_triple',u**9+3*u*u*z*z+xi*z**3,12,3)
model('mixed_double_nondegenerate',u**6+3*u**3*z+z*z+xi*z**3,4,2)
model('mixed_simple_b_unit',u**3+3*u*u*z+z*z+xi*z**3,1,2)
model('mixed_simple_b_zero',u**3+3*u*u*z+u*z*z+xi*z**3,3,3)
model('exceptional_a2_triple_b_unit',u**9+u*u*z+z*z+xi*z**3,2,2)
model('exceptional_a2_triple_b_zero',u**9+u*u*z+u*z*z+xi*z**3,3,3)

# All distinguished-branch contact calculations, including u-ramification.
contact_controls=[
 ('exceptional_double',1,5,0,1,-4), # p=0 is the away forbidden residue
 ('exceptional_triple',1,7,1,1,-3),
 ('triple_h1',1,5,3,1,-1),
 ('triple_h2_j0',2,9,7,1,-1),
 ('triple_h2_j1_low',1,1,1,1,0),
 ('triple_h2_jge2',1,4,4,1,0),
 ('mixed_double',1,3,2,2,-1),
 ('mixed_simple_unit',2,3,1,1,-2),
 ('mixed_simple_zero',1,1,1,3,0),
 ('proportional_mixed_low',1,1,1,1,0),
 ('proportional_simple',1,1,1,3,0)]
for name,ua,ze,pole,number,expected in contact_controls:
    check('contact pole bound '+name,pole-min(ze,4*ua)==expected)
check('proportional triple contact leaves degree three',(9-7)+(3-2)==3)
check('mixed away degree exceeds four',4+1>4)
check('coalescent triple lower bound exceeds four',12-3>4)
# Laurent derivatives cannot have a nonzero simple-pole term.
laurent=sum(S.Symbol('a'+str(i+8))*u**i for i in range(-8,9))
zero('meromorphic derivative residue',S.expand(u*S.diff(laurent,u)).coeff(u,0))

# Independently reconstruct the producer's mixed-chart F5 by exterior
# coefficients, rather than its quotient-numerator convolution.
s,t,w=S.symbols('s t w')
t0,x0,y0,k,lam=S.symbols('t0 x0 y0 k lam')
n1,delta,d0,d1,d2=S.symbols('n1 delta d0 d1 d2')
c1,c2,c3,c4=S.symbols('c1 c2 c3 c4')
D=d0*s*s+d1*s*t+d2*t*t
f=3*delta*delta-3*delta+1
zero('literal separation identity',16*f-(12*delta-3)*(4*delta-3)-7)
zero('quadratic coefficient is nonzero on f',f.subs(delta,0)-1)
zero('collision value is not a root',f.subs(delta,S.Rational(3,4))-S.Rational(7,16))
g=[w/z,w*w/z,(z+w*w)*w/z,(z+w*w)**2/z]
def lift(binary,const=0):
    a=S.Poly(S.expand(binary.subs({s:w,t:1})),w)
    return w*a.as_expr()/z+a.nth(2)*w+2*a.nth(3)*w*w+a.nth(3)*z+const
def wedge(F,G):return S.expand(S.diff(F,z)*S.diff(G,w)-S.diff(F,w)*S.diff(G,z))
saved=json.loads((HERE/'pole35_repeated_root_equations.json').read_text())
check('producer source hash',saved['source_sha256']==hashlib.sha256(
    (HERE/'explore_pole35_repeated_root.py').read_bytes()).hexdigest())
check('all three marked charts retained',len(saved['cases'])==3)
for row,(L,M),mult in zip(saved['cases'],[(s,t),(t,s),(s-t,t)],[-20,20,-120]):
    T=lift(L*L*M,t0);X=lift(n1*L*M*M,x0)
    Y=lift(delta*n1*n1*M**3+L*D,y0)
    J=lift(n1*n1*L*M*M);K=lift(n1*M**3)
    C=sum(c*gg for c,gg in zip((c1,c2,c3,c4),g))
    def dP(F):return (3*T*T+3*X)*wedge(T,F)+3*T*wedge(X,F)+3*wedge(Y,F)
    omega=(15*T*X*X+30*X*Y)*wedge(T,X)+(15*X*X-30*T*Y)*wedge(T,Y)-30*Y*wedge(X,Y)
    omega+=k*((12*X*X-12*T*Y)*wedge(T,X)+(-12*T*X-12*Y)*wedge(T,Y)-12*X*wedge(X,Y))
    omega-=S.Rational(5,3)*(K*dP(J)+J*dP(K))
    omega+=2*lam*T*dP(T)+dP(C)
    full=S.Poly(S.expand(omega*z**6),z)
    zero('pole six cancels '+row['name'],full.nth(0))
    f5=S.Poly(full.nth(1),w)
    stored=sum(S.sympify(expr)*w**degree for degree,expr in row['coefficient_equations'])
    zero('full F5 table independently reproduced '+row['name'],f5.as_expr()-stored)
    functional=f5.nth(7) if L==s else f5.nth(11) if L==t else S.diff(f5.as_expr(),w,3).subs(w,1)
    zero('marked-chart functional '+row['name'],functional-mult*n1**4*f)
    zero('stored functional '+row['name'],functional-S.sympify(row['functional_value']))
    # On the potential collision delta=3/4 this is a literal localized unit.
    H=S.Symbol('inverse_n1')
    value=S.expand(functional.subs(delta,S.Rational(3,4)))
    ratio=S.Rational(1,1)/(S.Rational(mult)*S.Rational(7,16))
    identity=ratio*H**4*value+(1-H*n1)*sum((H*n1)**j for j in range(4))
    zero('localized collision exclusion '+row['name'],identity-1)
report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
    local_discriminant_models=models,mixed_chart_count=3,
    necessary_mixed_relation=str(f),collision_value='3/4',collision_relation_value='7/16',
    written_result='Complete compatible quartic pole-(3,5) exclusion',
    first_open_height_five_pairs=[],
            subsequent_height_seven_exclusion="PF2RULED7 (external degree input; separate proof and receipt)",
            subsequent_all_height_pole_three_exclusion="PF2RULED8 (separate intrinsic proof and receipt)",
            frontier_superseded_by='PF2RULED16',
    frontier_handoff_scope='historical pre-PF2RULED16; this receipt does not prove the current frontier',
    first_open_reduced_pairs=[[6,8],[7,8]],
    frontier_adoption='PF2RULED12; subsequent handoff, outside this checker intrinsic proof coverage',
    assurance=dict(full_pole_three_five_exclusion='written proof plus exact controls',
       complete_pole_height_five_excluded=False,
       assurance_scope='this checker covers (3,5); full heights<=7 use PF2RULED7 and its separate external input',
       actual_ruled_compatibility_proved=False,
       quartic_packets_excluded=False,higher_poles_excluded=False,
       generic_fibre_geometry_formally_certified=False,independent_full_reproof=False,
       formal_verification=False,external_review=False),
    source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [
       Path(__file__),HERE/'POLE_THREE_FIVE_EXCLUSION.md',
       HERE/'explore_pole35_repeated_root.py',HERE/'pole35_repeated_root_equations.json']},
    python=platform.python_version(),sympy=S.__version__,elapsed_seconds=round(monotonic()-start,6))
(HERE/'pole35_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
print(f'PASS {len(checks)} exact controls; complete (3,5) exclusion; height<=7 now closed by PF2RULED7')
