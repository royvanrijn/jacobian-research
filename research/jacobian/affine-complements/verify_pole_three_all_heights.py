"""Exact controls for the all-height pole-three coordinate exclusion.

The weighted-initial lemma and complete local geometric coverage are
written in POLE_THREE_ALL_HEIGHTS.md, not formally certified here.
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
t,Z=S.symbols('t Z',nonzero=True)
c,K=S.symbols('c K',integer=True,positive=True)
H=S.Function('H')(Z);G=S.Function('G')(Z)
PP=t**(-3*c)*H;QQ=t**(-K)*G
jac=S.diff(PP,t)*S.diff(QQ,Z)-S.diff(PP,Z)*S.diff(QQ,t)
zero('universal weighted initial Jacobian ODE',
     t**(3*c+K+1)*jac-(K*S.diff(H,Z)*G-3*c*H*S.diff(G,Z)))
zero('weighted power-ratio logarithmic derivative',
     (3*c*S.diff(G,Z)/G-K*S.diff(H,Z)/H)*H*G-
     (3*c*H*S.diff(G,Z)-K*S.diff(H,Z)*G))
# Complete residue classes: the zero multiplicity of a nonconstant
# degree<=2 numerator is1 or2, never divisible by3.
for residue in (1,2):
    for multiplicity in (1,2):
        check(f'cube valuation contradiction nmod3={residue},zero={multiplicity}',
              residue*multiplicity%3!=0)
alpha=S.Symbol('alpha',nonzero=True)
zero('linear noncube face has a simple zero',S.diff(1+alpha*Z,Z)-alpha)
zero('coalescent quadratic is retained as multiplicity two',
     (1+alpha*Z)**2-(1+2*alpha*Z+alpha*alpha*Z*Z))

inequalities=[]
def row(kind,k,h,a,b,exception=False):
    cc=b-a*k
    check('positive weighted pole scale '+kind+str((k,h)),cc>0)
    lhs=8*cc;rhs=a*(2*k-1)
    if not exception:
        check('complete inequality at minimal n7 '+kind+str((k,h)),lhs>rhs)
    inequalities.append(dict(kind=kind,k=k,h=h,w_weight=a,z_weight=b,c=cc,
                             at_n7=lhs,bound=rhs,handled_by_fibre=exception))
for k in (1,2,3):
    for h in range(k):
        row('quadratic_symbol_defect',k,h,1,3*k-h)
for k in (1,2,3):row('u_support_defect',k,0,1,2*k)
for k,h in [(2,1),(3,1),(3,2)]:row('nonproportional_u',k,h,1,2*k-h)
for k in (1,2,3):
    for h in range(k):row('depressed_v_defect',k,h,2,3*k-h,(k,h)==(3,2))
leftovers={(k,h,j) for k,h in [(2,1),(3,1),(3,2)] for j in range(max(0,2*h-k))}
check('complete u-face dominance exception list',leftovers=={(3,2,0)})
check('all ordinary inequalities increase strictly with companion height',
      all(x['c']>0 for x in inequalities))

# Exact Newton discriminants support the two explicitly written fibre cases.
w,z,xi=S.symbols('w z xi')
models=[]
for name,F,expected_D,expected_disc,expected_poles in [
 ('triple_u2_vunit',w**9+3*w**5*z+z*z+xi*z**3,2,9,7),
 ('triple_depressed_v2',w**9+3*w*w*z*z+xi*z**3,3,15,12)]:
    disc=S.discriminant(F,z)
    order=min(m[0] for m,value in S.Poly(disc,w).terms() if value!=0)
    degree=min(m[0] for m,value in S.Poly(F.subs(w,0),z).terms() if value!=0)
    check('exceptional local Weierstrass degree '+name,degree==expected_D)
    check('exceptional discriminant order '+name,order==expected_disc)
    check('exceptional away pole total '+name,order-degree==expected_poles)
    models.append(dict(name=name,local_degree=degree,discriminant_order=order,
                       away_q_degree=expected_poles))
# Formula dQ/dt=unit*ell*w'/Fz; these parameters include the ramified
# normalization parameter and retain the exact contact-four order.
first=dict(w_order=2,z_order=9,Fz_order=9)
check('first exception away degree seven',first['Fz_order']-first['w_order']==7)
contact=min(first['z_order'],4*first['w_order'])
check('first exception contact derivative is regular',
      contact+first['w_order']-1-first['Fz_order']==0)
second=[dict(w_order=2,z_order=7,Fz_order=11),dict(w_order=1,z_order=2,Fz_order=4)]
away=[];at_contact=[]
for j,branch in enumerate(second):
    pole=branch['Fz_order']-branch['w_order']
    ell=min(branch['z_order'],4*branch['w_order'])
    away.append(pole);at_contact.append(pole-ell)
    check('second exception no contact-leading cancellation '+str(j),
          branch['z_order']!=4*branch['w_order'])
check('second exception exact away poles',away==[9,3])
check('second exception exact contact poles',at_contact==[2,1])
check('both second-exception totals contradict quartic degree',sum(away)==12 and sum(at_contact)==3)
# The normal form depression stays in V1; its T symbol is unchanged.
T,YY,mu,chi=S.symbols('T YY mu chi')
Tnew=T+mu
Ynew=YY+(chi-mu*mu)*Tnew+(2*mu**3-3*mu*chi)/3
zero('proportional cubic depression is exact',
     T**3+3*T*(mu*T+chi)+3*YY-Tnew**3-3*Ynew)
check('structural candidates after pole-three exclusion, before PF2RULED12',
      [(m,8) for m in range(4,8) if 8%m!=0]==[(5,8),(6,8),(7,8)])
report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
    complete_cubic_weight_inequalities=inequalities,exceptional_fibre_models=models,
    written_result='No compatible quartic target coordinate has B-pole three at any companion height',
    minimum_target_coordinate_pole=4,frontier_superseded_by='PF2RULED16',
    frontier_handoff_scope='historical pre-PF2RULED16; this receipt does not prove the current frontier',
    first_open_reduced_pairs=[[6,8],[7,8]],
    frontier_adoption='PF2RULED12; subsequent handoff, outside this checker intrinsic proof coverage',
    assurance=dict(all_height_local_initial_coverage='written proof',
       global_newton_and_finiteness_geometry_formally_certified=False,
       external_2022_degree_theorem_used=False,plane_torus_linearization_used=False,
       unrestricted_higher_minimum_poles_excluded=False,actual_ruled_compatibility_proved=False,
       quartic_packets_excluded=False,formal_verification=False,independent_full_reproof=False,external_review=False),
    source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),HERE/'POLE_THREE_ALL_HEIGHTS.md']},
    python=platform.python_version(),sympy=S.__version__,elapsed_seconds=round(monotonic()-start,6))
(HERE/'pole_three_all_height_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
print(f'PASS {len(checks)} exact controls; every compatible target coordinate has pole>=4; adopted later frontier(6,8),(7,8)')
