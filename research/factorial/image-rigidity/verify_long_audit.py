"""Exact transcription checks for the written EIT/FC2 audit.

No finite control here certifies Beukers' theorem, analytic continuation,
semisimplicity of every module, or universal factorial nonvanishing.
"""
from pathlib import Path
import hashlib
import json
import sympy as s

HERE = Path(__file__).resolve().parent
x, z = s.symbols('x z')
counts = {}


def check(group, condition):
    if not condition:
        raise AssertionError(group)
    counts[group] = counts.get(group, 0)+1


def coordinate(poly, dim):
    p = s.Poly(poly, x)
    return s.Matrix([p.nth(j) for j in range(dim)])


def moment(poly, cycle):
    primitive = s.Poly(poly, x).integrate()
    return sum(c*primitive.eval(a) for a, c in cycle)


def adjoint(A):
    # Column-stacked vector convention: vec(AP-PA).
    eye = s.eye(A.rows)
    return s.kronecker_product(eye, A)-s.kronecker_product(A.T, eye)


cases = [x**2+x, x**3-3*x, x**2*(x-1)**3, (x*x-1)**3+1, x**7]
for q in cases:
    n = s.degree(q, x)
    dim = n-1
    qp = s.diff(q, x)
    quotients, remainders = zip(*(s.div(s.expand(q*x**j), qp, x) for j in range(dim)))
    C0 = s.Matrix.hstack(*(coordinate(r, dim) for r in remainders))
    D0 = s.Matrix.hstack(*(coordinate(-s.diff(h, x), dim) for h in quotients))
    C, D = C0.T, D0.T
    check('global_residue_spectrum', D.charpoly(z).as_expr() == s.prod(z+s.Rational(j, n) for j in range(1, n)).expand())
    labels, columns = [], []
    critical = s.roots(qp, x)
    check('critical_points_complete', sum(critical.values()) == dim)
    for tau, multiplicity in critical.items():
        value = s.expand(q).subs(x, tau)
        for k in range(1, multiplicity+1):
            g = s.cancel(qp/(x-tau)**k).expand()
            column = coordinate(g, dim)
            check('multiplication_critical_eigenvector', C0*column == value*column)
            columns.append(column)
            labels.append((tau, multiplicity, k, value))
    G = s.Matrix.hstack(*columns)
    check('critical_basis_invertible', G.det() != 0)
    compressed = G.inv()*D0*G
    for j, (tau, multiplicity, k, value) in enumerate(labels):
        check('local_flag_diagonal', compressed[j,j] == -s.Rational(multiplicity+1-k, multiplicity+1))
        for i, (other, _, h, other_value) in enumerate(labels):
            if value == other_value and (other != tau or h > k):
                check('collision_flag_offdiagonal', compressed[i,j] == 0)

    for cycle in [[(0,-1),(1,1)], [(0,-2),(1,1),(2,1)]]:
        forcing = {}
        for a, coefficient in cycle:
            value = s.expand(q).subs(x, a)
            vector = s.Matrix([h.subs(x,a) for h in quotients])
            forcing[value] = forcing.get(value, s.zeros(dim,1))+coefficient*vector
        previous = None
        q_power = s.Poly(1, x)
        for m in range(9):
            actual = s.Matrix([moment(q_power.as_expr()*x**j, cycle) for j in range(dim)])
            expected = sum((vector*value**m for value, vector in forcing.items()), s.zeros(dim,1))
            if m:
                expected += m*C*previous
            check('actual_moment_ibp_system', (m*s.eye(dim)-D)*actual == expected)
            previous = actual
            q_power *= s.Poly(q, x)

    lo, hi, width = -1, 2, dim*dim
    adC, adD = adjoint(C), adjoint(D)
    system = s.zeros((hi-lo+2)*width, (hi-lo+1)*width)
    for row, power in enumerate(range(lo-1,hi+1)):
        if lo <= power <= hi:
            system[row*width:(row+1)*width,(power-lo)*width:(power-lo+1)*width] = -adC
        if lo <= power+1 <= hi:
            system[row*width:(row+1)*width,(power+1-lo)*width:(power+2-lo)*width] = (power+1)*s.eye(width)-adD
    solutions = system.nullspace()
    constants = adC.col_join(adD).nullspace()
    check('bounded_laurent_endomorphism_dimension', len(solutions) == len(constants))
    for solution in solutions:
        check('bounded_laurent_endomorphisms_constant', all(solution[(j-lo)*width+i] == 0
              for j in range(lo,hi+1) if j != 0 for i in range(width)))

w = s.Symbol('omega')
for degree in range(1,7):
    cyclotomic = s.Poly(s.cyclotomic_poly(degree,w), w)
    for power in range(19):
        filtered = s.rem(s.Poly(sum(w**(k*power) for k in range(degree)),w),cyclotomic).as_expr()/degree
        check('mittag_leffler_root_filter', filtered == (1 if power % degree == 0 else 0))
    negative = sum(w**((-k) % degree) for k in range(degree))/degree
    negative = s.rem(s.Poly(negative,w), cyclotomic).as_expr()
    check('degree_one_root_filter_correction', negative == (1 if degree == 1 else 0))

for degree in (4,6,8):
    cyclotomic = s.Poly(s.cyclotomic_poly(degree,w),w)
    for divisor in s.divisors(degree):
        if divisor == degree:
            continue
        for residue in range(divisor):
            trace = s.Poly(sum(w**j for j in range(residue,degree,divisor)),w)
            check('simultaneous_proper_block_trace', s.rem(trace,cyclotomic).is_zero)

volterra = -s.sqrt(x)/(z-x)+s.log(z-x)/(2*s.sqrt(x))
check('volterra_sign_negative_control', s.simplify(s.diff(s.sqrt(x)*s.log(z-x),x)-volterra) == 0)
check('volterra_corrected_density_zero', s.simplify(-s.sqrt(z)+s.integrate(1/(2*s.sqrt(x)),(x,0,z))) == 0)
omega = (-1+s.sqrt(3)*s.I)/2
beta = [s.Integer(1)/3,omega/3,omega**2/3]
check('bilinear_self_pairing_can_vanish', s.expand(sum(c*c for c in beta)) == 0)
check('conjugate_pairing_positive', s.simplify(sum(c*s.conjugate(c) for c in beta)) == s.Rational(1,3))

q = x*(x-1)
r = 2*x*x-x+1
check('polynomial_amplitude_scope_boundary', s.expand(s.diff(x,x)+x*s.diff(q,x)-r) == 0)
check('polynomial_amplitude_nonzero_initial_value', moment(r,[(0,-1),(1,1)]) == s.Rational(7,6))
check('transcendental_flat_limit_can_vanish', s.integrate(s.exp(2*s.pi*s.I*x),(x,0,1)) == 0)
check('transcendental_flat_first_moment_survives', 6+4*s.pi*s.I != 0)

u,v,eps,eta,y1,y2 = s.symbols('s t eps eta y1 y2')
ordinary = s.groebner([u*u*v*v],u,v)
for m in range(1,10):
    check('ordinary_normalization_t_squared', ordinary.reduce(v*v*u**m*s.factorial(m))[1] == 0 if m >= 2
          else ordinary.reduce(v*v*u**m)[1] != 0)
    check('reduced_normalization_t_is_wrong', ordinary.reduce(v*u**m*s.factorial(m))[1] != 0)
nil = s.groebner([eps**2,eta**3],eps,eta,y1,y2)
phase = eps*y1+eta*y2
check('finite_nilpotent_cutoff', nil.reduce(s.expand(phase**4))[1] == 0)
check('nilpotent_cutoff_sharp_control', nil.reduce(s.expand(phase**3))[1] != 0)
dual = s.groebner([eps**2],eps,y1,y2)
pure_nil = eps*(y1-y2)
check('nonzero_pure_nilpotent_phase', dual.reduce(pure_nil)[1] != 0)
check('pure_nilpotent_all_higher_powers', dual.reduce(s.expand(pure_nil**2))[1] == 0)
check('pure_nilpotent_first_moment', s.expand(eps*(s.factorial(1)-s.factorial(1))) == 0)

receipt = {
    'date':'2026-11-16', 'status':'pass', 'checks':counts, 'total':sum(counts.values()),
    'scope':'Finite exact transcriptions and negative controls for the written EIT/FC2 proof audit; no universal analytic theorem is inferred from the finite grids.',
    'proof_basis':'LONG_EIT_AND_FC2_AUDIT.tex plus pinned source proofs and identified published inputs',
    'files':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
             for name in ['verify_long_audit.py','LONG_EIT_AND_FC2_AUDIT.tex']},
    'independent_replay':False,'formal_verification':False,'external_review':False,
}
(HERE/'long_audit_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
