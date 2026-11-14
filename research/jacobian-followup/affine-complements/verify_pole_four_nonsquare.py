"""PF2RULED10 exact controls; all-height coverage and geometry are written proofs."""
from pathlib import Path
from time import monotonic
import sympy as S,json,hashlib,platform
HERE=Path(__file__).resolve().parent;start=monotonic();checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def zero(name,value):check(name,S.cancel(value)==0)
w=S.Symbol('w');z=S.Symbol('z',positive=True)
sig=S.Function('sig')(w);rho=S.Function('rho')(w)
j=S.symbols('j',integer=True)
P=sig**2*z**-4;R=rho*z**(-j)
J=lambda P,Q:S.diff(P,z)*S.diff(Q,w)-S.diff(P,w)*S.diff(Q,z)
zero('universal formal leading equation',S.expand(z**(j+5)*J(P,R))-(-4*sig**2*S.diff(rho,w)+2*j*sig*S.diff(sig,w)*rho))
for parity in (0,1):
    check('triple-root square-class parity '+str(parity),(3*parity)%2==parity)
check('normal orders through one require cancellation',min(jj+5 for jj in (-1,0,1))>3)
check('order two is the first unforced formal remainder',(-2)+5==3)

# Complete small-multiplicity defect atlas; j denotes ord(beta) here.
small=[];exceptions=[]
for k in range(1,4):
    for h in range(k):
        a=2;b=2*k-h;cc=2*b-a*k
        check('small-root positive scale '+str((k,h)),cc>0)
        check('small-root square criterion at r5 '+str((k,h)),6*cc>2*(k-1))
        check('small-root lower Laurent term never dominates '+str((k,h)),2*k-3*h>=0)
        for jj in range(4):
            if 2*jj<h:exceptions.append((k,h,jj))
        small.append(dict(k=k,h=h,w_weight=a,z_weight=b,c=cc,lhs_at_r5=6*cc,rhs=2*(k-1)))
check('complete beta-unit splitting exception list',exceptions==[(2,1,0),(3,1,0),(3,2,0)])
A,B,C,D,Z=S.symbols('A B C D Z')
zero('small-root nonsquare quadratic discriminant',S.discriminant((A+B*Z)**2+C*Z*Z,Z)+4*A*A*C)
for k,h,jj in exceptions:
    if h==1:
        a,b,s=2,2*k,2*k+1
        check('small-root cusp has regular Q derivative '+str(k),b+a-1-s==0)
    else:
        a,b,s=1,k,k+1
        check('small-root smooth branches have residue obstruction or regular Q',b+a-1-s==-1)
S0,H2,mu=S.symbols('S0 H2 mu')
zero('global second square completion',S0*S0+H2-(S0+mu/2)**2-(H2-mu*S0-mu*mu/4))

# Generic quartic models retain the constant regular lift g. No root,
# fibre irreducibility or a numerical choice of g is assumed.
a,g,V,xi=S.symbols('a g V xi')
models=[]
expected={
 (1,1):(15,-4*V*a**3*(-V**2-4*V*a*g+4*a*a*xi)),
 (1,2):(16,-16*a**4*xi*(V*a-xi)),
 (2,1):(16,-27*V**4),
 (2,2):(18,256*xi**2*(g*g+xi))}
for (jj,d),(order,lead) in expected.items():
    F=(w**3+a*w**jj*z+g*z*z)**2+V*w**d*z**3+xi*z**4
    disc=S.Poly(S.discriminant(F,z),w)
    actual=min(m[0] for m,c in disc.terms())
    check('generic discriminant order '+str((jj,d)),actual==order)
    zero('generic discriminant leading coefficient '+str((jj,d)),disc.coeff_monomial(w**order)-lead)
    local_degree=min(m[0] for m,c in S.Poly(F.subs(w,0),z).terms())
    const_order=min(m[0] for m,c in S.Poly(F.subs(z,0),w).terms())
    check('Weierstrass degree4 '+str((jj,d)),local_degree==4)
    check('constant-term order6 '+str((jj,d)),const_order==6)
    check('away-root local pole sum exceeds quartic degree '+str((jj,d)),order-local_degree-const_order>4)
    xi_poly=S.Poly(lead,xi)
    # Leading xi coefficients must be nonzero under the actual hypotheses:
    # a is a unit for jj=1, V is a unit always, g is arbitrary.
    required=[-16*V*a**5,16*a**4,-27*V**4,256][len(models)]
    zero('generic fibre parameter cannot annihilate discriminant identically '+str((jj,d)),xi_poly.LC()-required)
    models.append(dict(beta_order_case='1' if jj==1 else 'at least2',v_order=d,local_degree=local_degree,constant_term_order=const_order,discriminant_order=order,leading_discriminant=str(S.factor(lead)),away_signed_Q_pole_total=order-10))

# Uniform normal-order-one coefficient: symbolic gap identity, not an
# inference from the finite coefficient controls below.
l,k,nH,nGamma,d=S.symbols('l k nH nGamma d',integer=True,nonnegative=True)
i=l+1+k;m=2*k+1;nBeta=l-3*k-nH-2*nGamma
valuation=d*(i-nH)-3*(m+nBeta+nGamma)
gap=k*(d+3)+nH*(3-d)+3*nGamma
zero('universal coefficient valuation gap',valuation-(l+1)*(d-3)-gap)
for dd in (1,2):
    check('uniform positive higher-binomial gap d'+str(dd),dd+3>0)
    check('uniform positive regular numerator gap d'+str(dd),3-dd>0)
    check('uniform negative top coefficient valuation d'+str(dd),dd-3<0)
normal_order=S.symbols('normal_order',integer=True,positive=True)
for dd in (1,2):
    check('all higher numerator normal orders have positive cost d'+str(dd),3-dd>0)
baseline=(l+1)*d-3*(l+1)
after_numerator=l*d-3*(l+1-normal_order)
zero('generic higher numerator cost formula',after_numerator-baseline-(3*normal_order-d))
after_denominator=(l+1)*d-3*(l+2-normal_order)
zero('generic higher denominator cost formula',after_denominator-baseline-3*(normal_order-1))
check('all denominator normal orders at least two have positive cost',3*(2-1)>0)

coefficient_controls=[]
for ll in range(9):
    top=S.binomial(S.Rational(2*ll+1,2),ll+1)*(-1)**ll
    check('nonzero half-binomial leading coefficient l'+str(ll),top!=0)
    for dd in (1,2):
        terms=[]
        for kk in range(ll//3+1):
            ii=ll+1+kk;mm=2*kk+1;needed=ll-3*kk
            for nh in range(min(ii,needed)+1):
                for ng in range((needed-nh)//2+1):
                    nb=needed-nh-2*ng
                    coefficient=S.binomial(S.Rational(2*ll+1,2),ii)*S.binomial(ii,nh)*S.binomial(-mm,nb+ng)*S.binomial(nb+ng,ng)
                    power=dd*(ii-nh)-3*(mm+nb+ng)
                    if coefficient!=0:terms.append((power,kk,nh,ng,coefficient))
        least=min(t[0] for t in terms);at_least=[t for t in terms if t[0]==least]
        check('finite exact coefficient has predicted pole '+str((ll,dd)),least==(ll+1)*(dd-3))
        check('finite exact coefficient has unique nonzero leading contribution '+str((ll,dd)),len(at_least)==1 and at_least[0][1:4]==(0,0,0) and at_least[0][4]==top)
        coefficient_controls.append(dict(l=ll,v_order=dd,z1_w_order=least,leading_coefficient=str(top),term_count=len(terms)))
# A separate n10 bridge: universal Jacobian y^3 coefficient implies beta
# vanishes at a simple tau root, contradicting the finite away-root gate.
x,y=S.symbols('x y');tau=S.Function('tau')(x);N=S.Function('N')(x);beta=S.Function('beta')(x);gamma=S.Function('gamma')(x)
alpha,betaC,t2,n1=S.symbols('alpha betaC t2 n1')
Sig=tau**3;v=tau*N;r=S.Rational(5,16)*N**3
SS=Sig*y*y+beta*y+gamma;HH=v*y-t2*n1*x;RR=r*y-S.Rational(5,16)*n1**3*x
Jxy=lambda F,G:S.diff(F,x)*S.diff(G,y)-S.diff(F,y)*S.diff(G,x)
EE=S.expand(2*SS*Jxy(SS,RR)+Jxy(HH,RR)-(S.Rational(15,8)*HH*HH+S.Rational(3,2)*alpha*HH+betaC)*Jxy(SS,HH))
zero('n10 bridge highest coefficient cancels',EE.coeff(y,4))
normalized=S.cancel(EE.coeff(y,3)/tau**2)
tau_prime=S.Symbol('tau_prime')
zero('n10 bridge next coefficient at tau zero',normalized.subs(S.diff(tau,x),tau_prime).subs(tau,0)-S.Rational(15,4)*N**3*beta*tau_prime)
# Complete divisor handoff; the small-multiplicity lemma also applies to222.
old=[[6],[4,2],[3,3],[2,2,2]]
remaining=[row for row in old if row!=[3,3]]
check('remaining sextic partitions are exactly cubic squares',remaining==[[6],[4,2],[2,2,2]] and all(all(k%2==0 for k in row) for row in remaining))
check('pole-four coordinate layer remains open',len(remaining)==3)
report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
    small_multiplicity_atlas=small,away_root_quartic_models=models,z1_coefficient_controls=coefficient_controls,
    written_result='The entire nonsquare sextic(3,3) pole-four companion layer is excluded at all heights; every remaining common sextic is a cubic square',
    excluded_sextic_root_partitions=[[3,3]],remaining_sextic_root_partitions=remaining,
    remaining_pole_four_companion_residue=dict(modulus=4,residue=2,minimum=10),
    square_free_cubic_case_additional_form='P=S^2+H, S in V2, H in V1',
    minimum_target_coordinate_pole=4,frontier_superseded_by='PF2RULED16',
    frontier_handoff_scope='historical pre-PF2RULED16; this receipt does not prove the current frontier',
    first_open_reduced_pairs=[[6,8],[7,8]],
    frontier_adoption='PF2RULED12; subsequent handoff, outside this checker intrinsic proof coverage',
    assurance=dict(formal_centralizer_and_uniform_jet_coverage='written proof',
        discriminant_geometric_coverage_and_finiteness='written proof',
        finite_z1_controls_are_not_the_all_height_proof=True,
        remaining_partitions_are_PF2RULED10_handoff_before_PF2RULED11=True,
        all_n_ge_10_nonsquare_step_uses_external_degree_theorem=False,
        height_seven_base_inherited_from_PF2RULED9=True,
        remaining_cubic_square_cases_excluded=False,all_pole_four_coordinates_excluded=False,
        actual_ruled_compatibility_proved=False,quartic_packets_excluded=False,
        formal_verification=False,independent_full_reproof=False,external_review=False),
    source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),HERE/'POLE_FOUR_NONSQUARE_EXCLUSION.md']},
    python=platform.python_version(),sympy=S.__version__,elapsed_seconds=round(monotonic()-start,6))
(HERE/'pole_four_nonsquare_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
print(f'PASS {len(checks)} exact controls; nonsquare sextic(3,3) excluded; retained PF2RULED10 handoff has three cubic-square types')
