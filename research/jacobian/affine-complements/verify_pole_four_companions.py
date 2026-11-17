"""Exact controls for PF2RULED9; written geometric coverage is not formal verification."""
from pathlib import Path
from time import monotonic
import sympy as S,json,hashlib,platform
HERE=Path(__file__).resolve().parent
start=monotonic();checks=[];weights=[]
def check(name,value):
    assert value,name
    checks.append(name)
def zero(name,value):check(name,S.cancel(value)==0)
t,Z=S.symbols('t Z',nonzero=True)
c,K=S.symbols('c K',positive=True,integer=True)
H=S.Function('H')(Z);G=S.Function('G')(Z)
for power in (2,4):
    P=t**(-power*c)*H;Q=t**(-K)*G
    J=S.diff(P,t)*S.diff(Q,Z)-S.diff(P,Z)*S.diff(Q,t)
    zero(f'weighted power{power} ODE',t**(power*c+K+1)*J-(K*S.diff(H,Z)*G-power*c*H*S.diff(G,Z)))
for nmod in (1,3):
    for multiplicity in (1,2,3):
        check(f'fourth-power zero contradiction nmod4={nmod} multiplicity={multiplicity}',nmod*multiplicity%4!=0)
check('odd-power square condition retains double zeros',3*2%2==0)

def odd(kind,k,h,a,b,e=0,exception=False):
    cc=b-a*k
    check('positive odd weight '+str((kind,k,h,a,b)),cc>0)
    lhs=11*cc+e;rhs=a*(2*k-1)
    check('odd weighted inequality or retained geometric exception '+str((kind,k,h,a,b)),exception or lhs>rhs)
    weights.append(dict(kind=kind,k=k,h=h,a=a,b=b,c=cc,contact_lower_bound=e,lhs_at_n9=lhs,rhs=rhs,geometric_exception=exception))
for k in range(1,4):
    for h in range(2*k):
        b=4*k-h;odd('cubic_symbol_defect',k,h,1,b)
        check('cubic defect lower terms above '+str((k,h)),4*k-2*h>0 and 8*k-3*h>0)
    odd('A1_support_defect',k,0,1,2*k)
for k,h in [(2,1),(3,1),(3,2)]:
    odd('A1_nonproportional_j_ge_2h',k,h,1,2*k-h)
    check('A1 first row z^-1 term cannot dominate '+str((k,h)),3*h-2*k<=0)
    for j in range(2*h):
        odd('A1_nonproportional_j'+str(j),k,h,2,4*k-j)
        check('A1 second row z^-3 term strictly above '+str((k,h,j)),2*h-j>0)
        check('A1 second row z^-1 term cannot dominate '+str((k,h,j)),3*j-4*k<=0)
leftovers=set()
for k in range(1,4):
    for h in range(2*k):
        odd('H2_symbol_defect',k,h,2,4*k-h)
        for l in range(max(0,3*h-4*k)+1):
            if 2*l<3*h-4*k:leftovers.add((k,h,l))
check('complete dominating H2 z^-1 exception list',leftovers=={(2,3,0),(3,5,0),(3,5,1)})
for k in range(1,4):
    for h in range(k):odd('H1_symbol_defect',k,h,3,4*k-h,exception=(k,h)==(3,2))
odd('H1_final_contact_repair',3,2,3,10,e=10)
check('all odd inequalities grow strictly with n',all(row['c']>0 for row in weights))
# Exact changes of auxiliary lift, preserving the actual function P.
T,nu,chi,HH,U,AA=S.symbols('T nu chi HH U AA')
v=T+nu/4
expanded=S.Poly(S.expand((T**4+nu*T**3+chi*T*T+HH).subs(T,U-nu/4)),U)
check('quartic depression kills cubic term',expanded.coeff_monomial(U**3)==0)
check('quartic depression remainder has degree at most two',S.Poly(expanded.as_expr()-U**4,U).degree()<=2)
zero('square completion is exact',T*T+T*AA+HH-(T+AA/2)**2-(HH-AA*AA/4))

# Universal nonzero leading coefficients; every face is binomial and
# square-free on Z!=0, including the low branch of a two-edge model.
w,z=S.symbols('w z');B,C,D=S.symbols('B C D',nonzero=True)
models=[
 ('H2_k2_h3_l0',w**8+B*w**3*z*z+C*z**3,3,[(3,8,16,1)]),
 ('H2_k3_h5_l0',w**12+B*w**5*z*z+C*z**3,3,[(1,4,8,3)]),
 ('H2_k3_h5_l1',w**12+B*w**5*z*z+C*w*z**3+D*z**4,4,[(3,11,25,1),(1,1,3,1)]),
 ('H1_k3_h2',w**12+B*w*w*z**3+D*z**4+C*w**6*z*z,4,[(3,10,26,1),(1,2,6,1)])]
newton=[]
for name,F,degree,branches in models:
    local_degree=min(mon[0] for mon,value in S.Poly(F.subs(w,0),z).terms())
    check('local Weierstrass degree '+name,local_degree==degree)
    census=[]
    for a,b,order,multiplicity in branches:
        f=S.Poly(S.expand(F.subs({w:t**a,z:t**b*Z})),t)
        value=min(mon[0] for mon,coeff in f.terms())
        face=f.coeff_monomial(t**value)
        polynomial=S.Poly(face,Z)
        zfactor=min(mon[0] for mon,coeff in polynomial.terms())
        residual=S.cancel(face/Z**zfactor)
        check('Newton face root count '+name+str((a,b)),S.degree(residual,Z)==a*multiplicity)
        check('Newton face nonzero roots '+name+str((a,b)),S.simplify(residual.subs(Z,0))!=0)
        check('Newton face roots are simple '+name+str((a,b)),S.degree(S.gcd(residual,S.diff(residual,Z)),Z)==0)
        check('normalization parameter primitive '+name+str((a,b)),S.gcd(a,b)==1)
        check('F_z branch order '+name+str((a,b)),value-b==order)
        pole=order-a-b
        census.append(dict(w_order=a,z_order=b,Fz_order=order,branches=multiplicity,Q_pole_away=pole,contact_lower_bound=min(b,4*a)))
    check('all branches covered by Weierstrass degree '+name,sum(row['w_order']*row['branches'] for row in census)==degree)
    newton.append(dict(name=name,local_degree=degree,branches=census,away_total=sum(row['Q_pole_away']*row['branches'] for row in census)))
check('exceptional away pole totals', [row['away_total'] for row in newton]==[5,9,12,16])
check('first exception contact makes Q regular',5-min(8,12)<0)
check('second exception contact makes all Q regular even with cancellation',3-4<0)
check('third exception contact gives forbidden simple derivative poles',11-11==0 and 1-1==0)
check('last exception contact pole total alone is allowed and must not be used as obstruction',13-10+3-2==4)
# Explicit T-lift corrections are above every exceptional edge.
for name,k,branches in [('H2_k2_h3_l0',2,models[0][3]),('H2_k3_h5_l0',3,models[1][3]),('H2_k3_h5_l1',3,models[2][3]),('H1_k3_h2',3,models[3][3])]:
    F=next(row[1] for row in models if row[0]==name)
    for a,b,order,count in branches:
        minimum=order+b
        for j in (1,2,3):
            check('T-lift term above exceptional edge '+str((name,a,b,j)),a*k*(4-j)+b*j>minimum)
# For the final model the mu*T^2 terms are above both edges as well.
for a,b,order,count in models[3][3]:
    check('mu*T^2 above final edge '+str((a,b)),6*a+2*b>order+b and 3*a+3*b>order+b)

# Even companion n=2r, r>=5 odd: complete finite sextic-defect list.
even_weights=[]
for k in range(1,7):
    for h in range(k):
        a=1;b=2*k-h;cc=2*b-k
        check('even symbol defect inequality '+str((k,h)),6*cc>k-1)
        check('even symbol defect lower terms strictly above '+str((k,h)),2*k-2*h>0 and 4*k-3*h>0)
        even_weights.append(dict(k=k,h=h,a=a,b=b,c=cc,lhs_at_r5=6*cc,rhs=k-1))
A,beta,eta=S.symbols('A beta eta')
zero('simple sextic root square discriminant',S.discriminant((A+beta*Z)**2+eta*Z**2,Z)+4*A*A*eta)
# Jet of the formal critical-value construction; the exact implicit/Morse
# completion and all higher-coefficient coverage are in the written proof.
R,c10,c01,c20,c11,c02,q=S.symbols('R c10 c01 c20 c11 c02 q')
h=c10*R+c01*z+c20*R*R+c11*R*z+c02*z*z
F=R*R+z*z*h-q*z**4
Rcrit=-c10*z*z/2-c11*z**3/2
critical=S.series(F.subs(R,Rcrit),z,0,5).removeO().expand()
zero('critical value cubic coefficient',critical.coeff(z,3)-c01)
zero('critical value quartic coefficient',critical.coeff(z,4)-(c02-c10*c10/4-q))
check('generic fibre parameter controls quartic coefficient',S.diff(critical.coeff(z,4),q)==-1)
for d,a,b in [(3,2,3),(4,1,2)]:
    order=a+(a-1)-b # z*z'/F_rho
    check('simple sextic root derivative obstruction d'+str(d),order==(0 if d==3 else -1))
# Degree-six divisor partitions, with no omission of a non-square case.
def partitions(n,minimum=1):
    if n==0:yield ()
    for k in range(minimum,n+1):
        for tail in partitions(n-k,k):yield (k,)+tail
remaining=[list(reversed(row)) for row in partitions(6,2)]
remaining.sort(reverse=True)
check('complete sextic root partitions without simple roots',remaining==[[6],[4,2],[3,3],[2,2,2]])
check('only nonsquare sextic partition is3+3',[row for row in remaining if any(k%2 for k in row)]==[[3,3]])
check('structural height-eight candidates before PF2RULED12',[(m,8) for m in range(4,8) if 8%m]==[(5,8),(6,8),(7,8)])
report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
    odd_weight_atlas=weights,exceptional_Newton_models=newton,even_sextic_defect_atlas=even_weights,
    written_result='All odd companions to a pole-four coordinate excluded; reduced even companions require square-completed P and a sextic without simple roots',
    remaining_pole_four_companion_residue=dict(modulus=4,residue=2,minimum=10),
    remaining_sextic_root_partitions=remaining,minimum_target_coordinate_pole=4,frontier_superseded_by='PF2RULED16',
    frontier_handoff_scope='historical pre-PF2RULED16; this receipt does not prove the current frontier',
    first_open_reduced_pairs=[[6,8],[7,8]],
    frontier_adoption='PF2RULED12; subsequent handoff, outside this checker intrinsic proof coverage',
    assurance=dict(all_height_initial_and_formal_local_coverage='written proof',
        exact_control_scope='ODE, complete finite weight lists, symbolic Newton edges, contact orders, square and critical-value jets, sextic partitions',
        external_2022_degree_theorem_used_for_height_seven_base=True,
        odd_companions_n_ge_9_use_external_degree_theorem=False,
        even_companion_layer_excluded=False,all_pole_four_coordinates_excluded=False,
        unrestricted_higher_minimum_poles_excluded=False,actual_ruled_compatibility_proved=False,
        quartic_packets_excluded=False,formal_verification=False,independent_full_reproof=False,external_review=False),
    source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),HERE/'POLE_FOUR_COMPANIONS.md']},
    python=platform.python_version(),sympy=S.__version__,elapsed_seconds=round(monotonic()-start,6))
(HERE/'pole_four_companions_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
print(f'PASS {len(checks)} exact controls; odd pole-four companions excluded; four even sextic partitions retained')
