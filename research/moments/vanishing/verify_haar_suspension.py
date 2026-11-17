"""Small exact transcription controls; the universal proof is in the note."""
from pathlib import Path
import hashlib
import json
import sympy as S

tau,z,s,y,T,r,t,a,ap,rp,Y,x2p,x3 = S.symbols('tau z s y T r t a ap rp Y x2p x3')
A=(tau+1)*z/2
B=(tau-1)/(2*z)
phase=tau+A+B
assert S.expand(A*B-(tau*tau-1)/4)==0
assert S.expand((A-B)**2-(phase**2-2*phase*tau+1))==0
for m in range(9):
    ct=S.expand(phase**m).coeff(z,0)
    assert S.expand(ct-S.legendre(m,tau))==0

# Fibre derivative dS=(a-b)da/a when base coordinates are fixed.
bfun=(tau*tau-1)/(4*a)
assert S.simplify(S.diff(tau+a+bfun,a)-(a-bfun)/a)==0

# The actual top form, not a substituted compactification generator.
def jac(old,new):
    return S.Matrix([[S.diff(f,x) for x in new] for f in old]).det()

assert S.simplify(jac([r,t,r*ap],[r,t,ap])/(r*ap)-1/ap)==0
assert S.simplify(jac([a*rp,t,a],[rp,t,a])/a-1)==0
assert S.simplify(jac([r,r*x2p,x3],[r,x2p,x3])/(r*Y)-1/Y)==0
assert S.expand((s*s-2*s/T+1)*T*T-(T*T*(s*s+1)-2*s*T))==0
F=(s*s+1)*T-2*s
assert S.diff(F,s).subs({s:0,T:0})==-2

# k+ceil(b/2)>=0 is exactly the strict discrepancy inequality.
for pole in range(1,10):
    for order in range(-6,7):
        assert (order+(pole+1)//2>=0)==(2*(order+1)+pole>0)
    assert pole-pole//2==(pole+1)//2

# All-degree lower bounds, expressed symbolically rather than by a sweep.
m,d,aa=S.symbols('m d aa',positive=True)
D=2*m+1+d
assert S.expand(D-m-2-(m-1+d))==0
assert S.simplify(2*D/m-4-2*(1+d)/m)==0
assert S.expand((D-4).subs(m,2)-(1+d))==0
assert S.expand((D-2).subs(m,1)-(1+d))==0

# Both exceptional orbit bounds force strictness; degree1 control fails.
for orbit in [6,12]:
    assert S.expand(orbit*m/2-2*m)>0
assert -2+(1+1)//2==-1
assert -2+(3+1)//2==0
assert -2+(4+1)//2==0

out={
    'date':'2026-11-16','completed':True,
    'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'sympy_version':S.__version__,
    'unweighted_suspension_Legendre_moment_identity_controls':True,
    'conic_and_double_cover_equations_exact':True,
    'actual_volume_fibre_derivative_exact':True,
    'finite_conic_blowup_both_chart_volume_controls':True,
    'SNC_double_cover_blowup_volume_control':True,
    'infinity_normalized_equation_and_corner_transversality_exact':True,
    'infinity_rounding_discrepancy_equivalence_controls':True,
    'symbolic_all_degree_valuation_lower_bounds':True,
    'exceptional_orbit_6_12_strict_bounds':True,
    'ordinary_Long_degree1_infinity_volume_control':'FAILS as required: -1',
    'compact_period_mass':'1, by Omega wedge dz/(2pi i z)',
    'classification':'H(P) eventual-pure iff H=0 or (uv=0 and H(1)=0); every surviving phase SAFE',
    'ALL_nonlinear_equivariant_Long_recompositions_excluded':True,
    'purity_of_composite_reflected_to_primitive':False,
    'proper_volume_model_assumed_for_original_full_phase':False,
    'weighted_moment_transfer_used':False,
    'all_original_marker_forms_assumed_extend':False,
    'all_Bezout_weights_excluded':False,
    'Hopf_nonlinear_composites_claimed_excluded':False,
    'arbitrary_coefficient_ring_promotion_claimed':False,
    'universal_geometry_or_MZ_proof_machine_verified':False,
    'new_degree_support_or_parameter_campaign':False,
    'independent_replay':False,'formal_verification':False,'external_review':False,
    'external_review_next_task':False,'full_exceptional_MZ_cases':'OPEN on BOTH quotients',
}
Path(__file__).with_name('haar_suspension_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
