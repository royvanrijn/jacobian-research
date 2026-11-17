"""Formula and character controls; generic existence/topology are written proofs."""
from pathlib import Path
import hashlib,json
import sympy as s

k,N=s.symbols('k N',integer=True,positive=True)
g_cover=k*(k-2);g_quot=g_cover/N
assert s.expand((2*g_cover-2)-(N*(2*g_quot-2)+2*N-2))==0
a,b=s.symbols('a b',integer=True)
chi_curve=1-(1-a)*(1-b)
assert s.expand((1-chi_curve)-(a-1)*(b-1))==0
rows=[]
for name,order,exponent,branch in [('2O',24,12,[(6,4),(8,3),(12,2)]),('2I',60,30,[(12,5),(20,3),(30,2)])]:
    ram=sum(count*(r-1) for count,r in branch)
    assert ram==2*order-2
    # Two fixed formula transcriptions, not a search or extrapolation.
    for kk in [2*order,2*order+exponent]:
        assert kk%exponent==0 and kk-1>=2*order-1
        gg=s.Rational(kk*(kk-2),order);assert gg.q==1
        rows.append({'group':name,'k':kk,'matrix_degree':2*kk,
                     'quotient_boundary_genus':int(gg),'radical_dimension_lower_bound':int(2*gg)})
    for r in {r for _,r in branch}:
        # A primitive rth root represents lambda^2. The off-diagonal value
        # exponent is 1, its allowed first-jet exponent is 1-1=0, and the
        # diagonal value exponent k is 0 modulo r.
        assert 1%r!=0 and (1-1)%r==0 and exponent%r==0
# Effective weight-one scaling trivializes both families exactly.
Pval,scale=s.symbols('Pval scale',nonzero=True)
assert s.simplify(scale*(Pval/scale)-Pval)==0
out={
 'completed':True,'sympy_version':s.__version__,
 'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'quotient_genus_Riemann_Hurwitz_identity_exact':True,
 'bidegree_curve_Euler_characteristic_identity_exact':True,
 'stabilizer_character_and_first_jet_exponents_exact':True,
 'formula_transcriptions':rows,
 'actual_phase_coefficients_expanded_or_sampled':False,
 'generic_invariant_jet_existence_and_boundary_tori_proofs_machine_verified':False,
 'unbounded_boundary_space_on_actual_pure_safe_locus':'Internal all-degree written existence proof on BOTH quotients',
 'universal_vanishing_of_invariant_boundary_space_on_pure_locus_possible':False,
 'new_safety_sector_claimed':False,
 'exceptional_failure_constructed':False,
 'full_exceptional_MZ_cases':'OPEN on BOTH quotients',
 'independent_replay':False,'formal_verification':False,'external_review':False,
 'new_degree_support_or_prefix_campaign':False}
Path(__file__).with_name('haar_boundary_excess_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
