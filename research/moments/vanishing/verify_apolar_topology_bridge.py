"""Small exact transcription checks for the fixed-ideal topology adoption audit.

The all-index obstruction is the written Hermite translation proof, not this
finite replay. No search, recurrence fitting, or general gap decision occurs.
"""
from pathlib import Path
from math import factorial
import hashlib
import json
import sympy as s


def apply_symbol(symbol, variables, dual, polynomial):
    result = s.S.Zero
    for powers, coefficient in s.Poly(symbol, *dual).terms():
        term = polynomial
        for variable, degree in zip(variables, powers):
            if degree:
                term = s.diff(term, variable, degree)
        result += coefficient * term
    return s.expand(result)


def verify():
    x, y, X, E, L, a, b = s.symbols('x y X E L a b')
    v, u, V, U = s.symbols('v u V U')
    f = x*x + y
    full = x*x + v + u
    hermite, contractions, actual = [], [], []
    for n in range(1, 7):
        H = s.expand(factorial(n)*sum((-E)**j * X**(n-2*j)
            / (s.Integer(factorial(n-2*j))*factorial(j))
            for j in range(n//2+1)))
        killed = apply_symbol(H, (x,y), (X,E), f**(n-1))
        square = apply_symbol(H*H, (x,y), (X,E), f**n)
        assert killed == 0
        assert square == 2**n*factorial(n)**2
        # Simultaneous translation cancels the individual squares exactly.
        translated = s.expand((f + 2*x*(a+b) + 2*a*b)**n)
        coefficient = translated.coeff(a,n).coeff(b,n)
        assert coefficient == 2**n
        hermite.append(dict(n=n, H=str(H), killed_power=n-1,
            square_surviving_power=n, square_output=str(square)))
    for n in range(1, 5):
        R = s.expand(factorial(2*n)*sum((-1)**j * L**(n-j)*E**(2*j)
            / (s.Integer(factorial(2*n-2*j))*factorial(j))
            for j in range(n+1)))
        H = s.expand(factorial(2*n)*sum((-E)**j*X**(2*n-2*j)
            / (s.Integer(factorial(2*n-2*j))*factorial(j))
            for j in range(n+1)))
        assert s.expand(R.subs(L,X*X*E)-E**n*H) == 0
        symbol = s.expand(R.subs({L:X*X*E}))
        killed = apply_symbol(symbol, (x,y), (X,E), f**(3*n-1))
        square = apply_symbol(symbol*symbol, (x,y), (X,E), f**(4*n))
        target = 2**(2*n)*factorial(4*n)*factorial(2*n)
        assert killed == 0 and square == target
        contractions.append(dict(n=n, R=str(R), killed_power=3*n-1,
            square_surviving_power=4*n, square_output=str(square)))
        if n <= 2:
            full_symbol = s.expand(R.subs({L:X*X*V,E:U}))
            assert apply_symbol(full_symbol,(x,v,u),(X,V,U),full**(3*n-1)) == 0
            full_square = apply_symbol(full_symbol*full_symbol,(x,v,u),(X,V,U),full**(4*n))
            assert full_square == target
            actual.append(dict(n=n, full_passive_suspension=True, square_output=str(full_square)))
    profiles=[]
    for m in range(1, 9):
        outputs=[apply_symbol((X*X*V)**k,(x,v,u),(X,V,U),full**m)
                 for k in range(m+1)]
        height=max(k for k,output in enumerate(outputs) if output != 0)
        assert height == m//2 and outputs[m] == 0
        boundary=apply_symbol((X*X*V)**m,(x,v,u),(X,V,U),v**m*full**m)
        assert boundary == factorial(2*m)*factorial(m)
        profiles.append(dict(m=m,height=height,deficit=m-height,boundary_output=str(boundary)))
    guards=0
    for invalid in [dict(fixed_symbol_height_subadditivity_refuted=True),
                    dict(general_actual_gap_decided=True),
                    dict(all_swanson_representations_excluded=True)]:
        try:
            assert all(value is False for value in invalid.values())
        except AssertionError:
            guards+=1
    assert guards==3
    here=Path(__file__).resolve()
    return dict(completed=True,checker_sha256=hashlib.sha256(here.read_bytes()).hexdigest(),
        hermite_controls=hermite,contracted_controls=contractions,
        full_passive_suspension_controls=actual,actual_profiles=profiles,
        cutoff='1',height_gap_fraction='1/2',guard_refusals=guards,
        all_index_proof='Written simultaneous-translation identity and exact factorial formula,not bounded checks.',
        general_actual_gap_decided=False,actual_eventual_but_nonlinear_constructed=False,
        fixed_symbol_height_subadditivity_refuted=False,
        all_swanson_representations_excluded=False,
        finite_symbolic_rees_generation_required_by_swanson=False,
        external_review=False,independent_replay=False,formal_verification=False)


if __name__=='__main__':
    print(json.dumps(verify(),indent=2))
