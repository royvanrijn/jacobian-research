"""Exact selected controls for written real-seed density and trace proofs.

No moment-prefix, degree, seed or density search; no numerical conformal map.
The connected-tube, curve-period and Haar conclusions are not machine proofs.
"""
from pathlib import Path
import hashlib,json,platform
import sympy as s
here=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
t,u,S,e=s.symbols('t u S e',positive=True)

# A closed positive source loses the interior lower bound if C1 regularity
# is omitted: the transformed sine law has one-sided mass o(epsilon).
density=2*u/(s.pi*s.sqrt(1-u**4))
cdf=s.asin(e**2)/s.pi
assert s.simplify(s.diff(cdf,e)-density.subs(u,e))==0
assert s.limit(cdf/e,e,0,dir='+')==0
assert cdf.subs(e,1)-cdf.subs(e,0)==s.Rational(1,2)

# Divided-difference reconstruction coefficients, including nonmonic p.
divided_controls=[]
for p in [2*t**3+3*t-1,t**4-s.Rational(3,2)*t**2+s.Rational(1,4)+s.I*(t-t**3)]:
 d=int(s.degree(p,t));coeff=[s.expand(p).coeff(t,k) for k in range(d+1)]
 bs=[sum(coeff[k]*t**(k-1-j) for k in range(j+1,d+1)) for j in range(d)]
 assert s.expand((t-u)*sum(bs[j]*u**j for j in range(d))-(p-p.subs(t,u)))==0
 assert s.expand(sum(bs[j]*t**j for j in range(d))-s.diff(p,t))==0
 divided_controls.append({'degree':d,'leading_coefficient':str(coeff[d]),'kernel_and_diagonal_exact':True})

# A two-atom measure verifies both trace directions exactly without treating
# finite controls as a full-interval nonrationality proof.
p=2*t**2+3*t-1
G=1/(3*(t+1))+2/(3*(t-2))
F0=s.Rational(1,3)/(S+2)+s.Rational(2,3)/(S-13)
F1=-s.Rational(1,3)/(S+2)+s.Rational(4,3)/(S-13)
roots=[(-3+s.sqrt(8*S+17))/4,(-3-s.sqrt(8*S+17))/4]
for j,F in [(0,F0),(1,F1)]:
 trace=sum(xi**j*G.subs(t,xi)/(4*xi+3) for xi in roots)
 assert s.simplify(trace-F)==0
assert s.cancel((2*t+3)*F0.subs(S,p)+2*F1.subs(S,p)-G)==0

# Nonzero total mass and the differential pole hypothesis are essential.
# On the circle f=z,omega=dz/z has mass2*pi*i and all positive moments0,
# but its pole at0 does not lie over f=infinity.
z=s.Symbol('z')
assert s.residue(1/z,z,0)==1
assert s.residue(z**2/z,z,0)==s.residue(z**5/z,z,0)==0
# A zero-mass polynomial density can have every moment0 with nonzero phase.
m=s.Symbol('m',integer=True,nonnegative=True)
pzero=t*(t-1);qzero=s.diff(pzero,t)
primitive=pzero**(m+1)/(m+1)
assert s.simplify(s.diff(primitive,t)-pzero**m*qzero)==0
assert pzero.subs(t,0)==pzero.subs(t,1)==0
assert s.integrate(qzero,(t,0,1))==0

# Finite branch-point pullbacks of polynomial/hyperelliptic densities have
# holomorphic differentials. These are selected local controls, not all genus.
F=t-t**3
assert s.discriminant(F,t)!=0
assert s.expand(F.subs(t,u**2))==u**2-u**6
holomorphic_pullback=2/s.sqrt(1-u**4)
assert s.limit(holomorphic_pullback,u,0)==2
assert holomorphic_pullback.series(u,0,5).removeO()==u**4+2
# Every holomorphic local form has worst Cauchy exponent1/k-1>-1.
local_controls=[]
for k in [2,3,5]:
 exponent=s.Rational(1,k)-1
 assert exponent>-1
 assert s.limit(e**(exponent+1)*s.log(e),e,0,dir='+')==0
 local_controls.append({'ramification':k,'worst_exponent':str(exponent),'subpolar':True})

receipt={
 'date':'2026-11-15','completed':True,'checker_sha256':sha(Path(__file__)),
 'python_version':platform.python_version(),'sympy_version':s.__version__,
 'non_C1_circle_seed_density':str(density),
 'non_C1_one_sided_mass':str(cdf),'non_C1_mass_over_epsilon_limit':'0',
 'divided_difference_controls':divided_controls,
 'two_atom_nonmonic_trace_and_reconstruction_exact':True,
 'two_atom_control_claimed_to_have_full_interval_support':False,
 'pole_hypothesis_boundary_circle_zeroth_residue':1,
 'zero_mass_polynomial_density':str(qzero),
 'hyperelliptic_branch_pullback_limit':'2','local_Cauchy_controls':local_controls,
 'universal_interior_seed_bound':'h_*measure >=c_J Lebesgue on every compact J inside h(M),including critical values',
 'quartic_preceding_specific_law_Haar_realization':'IMPOSSIBLE for a C1 real seed on a connected positive smooth source,in ANY positive dimension',
 'quartic_polynomial_graph_or_other_laws_excluded':False,
 'seed_trace_formula':'F_j(s)=sum_(p(xi)=s) xi^j G(xi)/p_prime(xi),0<=j<deg p',
 'seed_inverse_formula':'G(t)=sum_j b_j(t) F_j(p(t))',
 'all_marked_basis_resolvents_may_be_finite_Laurent_under_nonzero_purity':False,
 'compact_curve_pole_criterion':'poles(omega) subset f^-1(infinity);rational Cauchy period implies zero,including finite endpoint chains and finite sums',
 'published_compact_genus_scope':'Gavrilov-Pakovich Section2 treats closed curves;finite endpoint subpolar argument is written here',
 'all_algebraic_or_all_Haar_seed_laws_assumed_pole_compatible':False,
 'new_degree_density_seed_or_moment_prefix_campaign':False,
 'universal_tube_analytic_and_curve_period_proofs_machine_verified':False,
 'either_exceptional_MZ_decision':'OPEN',
 'independent_replay':False,'formal_verification':False,'external_review':False}
(here/'seed_density_trace_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
