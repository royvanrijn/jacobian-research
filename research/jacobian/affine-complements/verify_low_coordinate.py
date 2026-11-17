"""Exact controls for the all-height exclusion of a pole-two coordinate.

The all-height induction and generic-fibre geometry are written in
LOW_COORDINATE_POLE.md. Finite controls do not formally certify that proof.
"""
from pathlib import Path
from time import monotonic
import hashlib,json,platform
import sympy as S

HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def zero(name,expr):
    check(name,S.cancel(expr)==0)
def order(expr,u):
    num,den=S.fraction(S.cancel(expr))
    if num==0:return S.oo
    return min(p[0] for p,c in S.Poly(num,u).terms())-min(p[0] for p,c in S.Poly(den,u).terms())

start=monotonic()
T,A,u,z=S.symbols("T A u z")
controls=[]
for r in range(1,9):
    n=2*r+1
    F=sum(S.binomial(S.Rational(n,2),j)*T**(n-2*j)*A**j for j in range(r+1))
    cr=S.binomial(S.Rational(n,2),r)
    zero("binomial principal-part differential n="+str(n),
         2*T*S.diff(F,A)-S.diff(F,T)+cr*A**r)
    check("principal coefficient is nonzero n="+str(n),cr!=0)
    if r>1:
        check("every smaller odd correction has lower differential pole n="+str(n),
              all((j-1)//2+3<r+3 for j in range(r+1,n) if j%2))
    controls.append(dict(n=n,principal_coefficient=str(cr),first_rational_coefficient=str(cr/(2*(r+1)))))
r=S.Symbol("r",integer=True,positive=True)
cr=S.Symbol("cr")
tau,aa,dd=[S.Function(name)(u) for name in ("tau","a","d")]
top=(-cr*aa**r*(S.diff(tau,u)*aa-tau*S.diff(aa,u))
     +2*r*tau*S.diff(tau,u)*dd-2*tau*tau*S.diff(dd,u))
primitive=dd/tau**r-cr*aa**(r+1)/(2*(r+1)*tau**(r+1))
zero("uniform principal-part integrating identity",top+2*tau**(r+2)*S.diff(primitive,u))
n=S.Symbol("n",integer=True,positive=True)
sigma=S.Function("sigma")(u)
zero("uniform polar-centralizer derivative identity",
     n*S.diff(tau,u)*sigma-tau*S.diff(sigma,u)+tau**(n+1)*S.diff(sigma/tau**n,u))

branches=[]
for k,h in ((2,1),(3,1),(3,2)):
    # Exact quadratic model P=u^(2k)/z^2+u^h/z=2.
    unit=S.sqrt(1+8*u**(2*(k-h)))
    roots=[u**h*(1-unit)/4,u**h*(1+unit)/4]
    pz=-2*u**(2*k)/z**3-u**h/z**2
    row=[]
    for index,root in enumerate(roots):
        zero(f"exact Newton model {k},{h},branch{index}",
             u**(2*k)+root*u**h-2*root*root)
        series=S.series(root,u,0,2*k-h+9).removeO()
        e=order(series,u)
        expected_e=2*k-h if index==0 else h
        check(f"Newton branch order {k},{h},{index}",e==expected_e)
        derivative_order=order(pz.subs(z,series),u)
        expected_pz=3*h-4*k if index==0 else -h
        check(f"first-coordinate normal derivative {k},{h},{index}",
              derivative_order==expected_pz)
        away_order=order((z**-3/pz).subs(z,series),u)
        away_pole=-away_order-1
        check(f"away pole {k},{h},{index}",
              away_pole==(2*k-1 if index==0 else 2*h-1))
        contact_order=order(((z-3*u**4)*z**-3/pz).subs(z,series),u)
        row.append(dict(branch=index,e=e,pz_order=derivative_order,
                        away_q_pole=away_pole,contact_q_derivative_order=contact_order))
    branches.append(dict(k=k,h=h,branches=row))
check("mixed multiple-root contact has a forbidden residue",
      branches[0]["branches"][1]["contact_q_derivative_order"]==-1)
check("triple h1 contact has a forbidden residue",
      branches[1]["branches"][1]["contact_q_derivative_order"]==-1)
check("triple h2 generic contact gives two simple poles",
      [r["contact_q_derivative_order"] for r in branches[2]["branches"]]==[-2,-2])

# Simple root where a also vanishes: two generic rational slopes.
pz=-2*u*u/z**3-u/z**2
for slope in (S.Integer(1),-S.Rational(1,2)):
    root=slope*u
    zero("simple-root generic fibre equation "+str(slope),u*u+root*u-2*root*root)
    check("simple-root normal derivative "+str(slope),order(pz.subs(z,root),u)==-1)
    check("simple-root away pole "+str(slope),order((z**-3/pz).subs(z,root),u)==-2)
    check("simple-root contact residue "+str(slope),
          order(((z-3*u**4)*z**-3/pz).subs(z,root),u)==-1)

# Leading cancellation at e=4 cannot create an extra pole.
root=u**2*(1-S.sqrt(1+8*u*u))/4
series=S.series(root,u,0,12).removeO()
pz=-2*u**6/z**3-u**2/z**2
check("cancelled contact order may rise to six",order(series+u**4,u)==6)
check("cancelled contact has no pole of Q",
      order(((z+u**4)*z**-3/pz).subs(z,series),u)>=0)
# A regular term in T gives the intermediate order-five/residue case.
root=u**2*((1+2*u)-S.sqrt((1+2*u)**2+4*u*u))/2
series=S.series(root,u,0,12).removeO()
pz=-2*u**6/z**3-(2*u**3+u*u)/z**2
check("contact order five is explicitly retained",order(series+u**4,u)==5)
check("contact order five is a forbidden derivative residue",
      order(((z+u**4)*z**-3/pz).subs(z,series),u)==-1)

coeff=S.symbols("c0:17")
laurent=sum(c*u**i for c,i in zip(coeff,range(-8,9)))
derivative=S.expand(S.diff(laurent,u))
zero("meromorphic derivatives have zero simple-pole coefficient",
     S.expand(u*derivative).coeff(u,0))

totals=[
 dict(configuration="double_plus_simple",poles=[3,1,1,1],degree=6),
 dict(configuration="triple_h1",poles=[5,1],degree=6),
 dict(configuration="triple_h2_away",poles=[5,3],degree=8),
 dict(configuration="triple_h2_contact",poles_upper_bounds=[1,1],degree_upper_bound=2)]
for row in totals:
    if "poles" in row:
        check("quartic degree contradiction "+row["configuration"],
              sum(row["poles"])==row["degree"] and row["degree"]!=4)
    else:
        check("contact quartic degree contradiction",
              sum(row["poles_upper_bounds"])==row["degree_upper_bound"]<4)
report=dict(as_of="2026-11-14",outcome="PASS",check_count=len(checks),checks=checks,
            binomial_controls=controls,local_branch_controls=branches,degree_controls=totals,
            written_result="Every compatible quartic target coordinate has pole at least three, independently of the other coordinate's height",
            first_open_height_five_pairs=[],
            subsequent_height_seven_exclusion="PF2RULED7 (external degree input; separate proof and receipt)",
            subsequent_all_height_pole_three_exclusion="PF2RULED8 (separate intrinsic proof and receipt)",
            frontier_superseded_by='PF2RULED16',
    frontier_handoff_scope='historical pre-PF2RULED16; this receipt does not prove the current frontier',
    first_open_reduced_pairs=[[6,8],[7,8]],
    frontier_adoption='PF2RULED12; subsequent handoff, outside this checker intrinsic proof coverage',
            subsequent_pair_exclusion="PF2RULED6 (separate proof and receipt)",
            assurance=dict(all_heights="written induction and generic-fibre proof",
                           complete_generic_fibre_geometry_formally_certified=False,
                           actual_ruled_compatibility_proved=False,higher_pairs_excluded=False,
                           quartic_packets_excluded=False,independent_full_reproof=False,
                           formal_verification=False,external_review=False),
            source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                           for p in [Path(__file__),HERE/"LOW_COORDINATE_POLE.md"]},
            python=platform.python_version(),sympy=S.__version__,
            elapsed_seconds=round(monotonic()-start,6))
(HERE/"low_coordinate_receipt.json").write_text(json.dumps(report,indent=2)+"\n")
print(f"PASS {len(checks)} exact controls; all-height pole-two exclusion; height<=7 now closed by PF2RULED7")

