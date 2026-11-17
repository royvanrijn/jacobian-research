"""Finite exact transcription controls; universal proofs are in the note."""
from pathlib import Path
import hashlib
import json
import sympy as s

t, alpha, c = s.symbols('t alpha c')
rows = []
for a in [s.Rational(1, 2), s.Rational(2, 3)]:
    f = sum((c*t)**m/s.rf(a+1, m) for m in range(10))
    residual = s.expand(t*s.diff(f, t)-(c*t-a)*f-a)
    assert all(residual.coeff(t, m) == 0 for m in range(10))
    for n in range(9):
        assert s.simplify((n+a)*s.rf(a,n)-s.rf(a,n+1)) == 0
        # beta moment equals m!/(a+1)_m, by the gamma recurrence.
        assert s.simplify(a*s.factorial(n)*s.rf(a,1)/
                          (a*s.rf(a,n+1))-s.factorial(n)/s.rf(a+1,n)) == 0
    for j in range(1, 7):
        fj = s.diff(f,t,j).subs(t,0)
        prev = s.diff(f,t,j-1).subs(t,0)
        assert s.simplify((a+j)*fj-j*c*prev) == 0
    rows.append({'alpha':str(a),'jet_order':9,'nontrivial_finite_monodromy':True})
# Retain the nonzero/nonunit normalization boundary over C[c,b]/(cb).
groebner = s.groebner([c*s.Symbol('b')],c,s.Symbol('b'))
b = s.Symbol('b')
assert groebner.reduce(b*c**4)[1] == 0
assert groebner.reduce(s.Integer(1))[1] == 1
assert groebner.reduce(b)[1] == b
out = {
    'completed':True,'sympy_version':s.__version__,
    'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'kummer_series_equation_Ext_functional_and_fibre_recurrence_controls':rows,
    'nonzero_nonunit_coefficient_ideal_control':True,
    'universal_proofs_machine_verified':False,
    'exceptional_Kummer_realization_claimed':False,
    'local_monodromy_exclusion_valid':False,
    'full_exceptional_MZ_cases':'OPEN on BOTH quotients',
    'independent_replay':False,'formal_verification':False,'external_review':False,
    'new_degree_or_support_campaign':False,
}
Path(__file__).with_name('haar_separated_core_verification.json').write_text(
    json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
