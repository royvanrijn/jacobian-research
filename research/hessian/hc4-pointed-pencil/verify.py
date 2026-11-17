"""Exact pencil-control identities; global geometry is proved in PROOF.md.

The g=0..8 representatives regress the written arbitrary-genus proof.
They are not an enumeration of Keller maps or a genus-bound search.
"""
from pathlib import Path
import hashlib
import json
import platform
import time
import sympy as S

HERE = Path(__file__).resolve().parent
checks = []


def check(name, condition):
    assert condition, name
    checks.append(name)


def zero(name, expr):
    values = list(expr) if isinstance(expr, S.MatrixBase) else [expr]
    check(name, all(S.cancel(v) == 0 for v in values))


def main():
    start = time.monotonic()
    x, y, t, W, v, a, b = S.symbols('x y t W v a b')
    # Universal blowup determinant formula, with no chosen G.
    aa, bb = S.Function('A')(a, b), S.Function('B')(a, b)
    original = S.Matrix([a*aa, a*bb]).jacobian((a, b)).det()
    lifted = S.Matrix([a*aa, bb/aa]).jacobian((a, b)).det()
    zero('universal blowup determinant equals J_original/A', lifted-original/(a*aa))
    # A nonlinear Keller example, not a nonproperness control.
    g1 = y+x*x
    g2 = 2*x+2*g1*g1
    row = S.Matrix([[g1, g2]])
    zero('nonlinear automorphism orientation', row.jacobian((x, y)).det()+2)
    xx = g2/2-g1*g1
    yy = g1-xx*xx
    zero('nonlinear automorphism inverse x', xx-x)
    zero('nonlinear automorphism inverse y', yy-y)
    blow = row.subs({x: a, y: a*b}, simultaneous=True)
    aa_example = S.cancel(blow[0]/a)
    bb_example = S.cancel(blow[1]/a)
    zero('nonlinear example lifted etale determinant',
         S.Matrix([a*aa_example, bb_example/aa_example]).jacobian((a, b)).det()
         +2/aa_example)
    zero('nonlinear exceptional map is projectivized linear map',
         (bb_example/aa_example).subs(a, 0)-2/b)
    controls = []
    ell = S.Matrix([[y, 2*x]])
    for genus in range(9):
        h1 = y+x*y*y+x**(2*genus+1)
        h2 = 2*x
        hrow = S.Matrix([[h1, h2]])
        D = S.Matrix([[1, 0], [(y*y+x**(2*genus))/2, 1]])
        tag = f'g={genus}'
        zero(tag+' elementary matrix determinant', D.det()-1)
        zero(tag+' row identity', ell*D-hrow)
        zero(tag+' literal x ideal certificate', h2/2-x)
        zero(tag+' literal y ideal certificate', h1-(y*y+x**(2*genus))*h2/2-y)
        zero(tag+' same entire nonconstant Jacobian',
             hrow.jacobian((x, y)).det()+2*(1+2*x*y))
        F = 1-4*x**(2*genus+2)-8*t*x*x
        zero(tag+' hyperelliptic pencil identity',
             (2*x*y+1)**2-F-4*x*(h1+t*h2))
        poly = S.Poly(F, x, domain=S.QQ.frac_field(t))
        check(tag+' generic squarefree double cover', poly.gcd(poly.diff()).degree() == 0)
        check(tag+' degree and Riemann-Hurwitz genus',
              poly.degree() == 2*genus+2 and (poly.degree()-2)//2 == genus)
        check(tag+' t=0 squarefree specialization',
              S.gcd(F.subs(t, 0), S.diff(F, x).subs(t, 0), x) == 1)
        zero(tag+' origin branch is regular', F.subs(x, 0)-1)
        zero(tag+' regular branch has W-1 order at least two', S.diff(F, x).subs(x, 0))
        # The remaining three punctures and the universal genus statement
        # are written geometric deductions, not separate CAS genus output.
        controls.append(dict(parameter=genus, branch_degree=int(poly.degree()),
                             genus_from_Riemann_Hurwitz=genus,
                             punctures_from_written_inverse=3,
                             first_Betti_number_from_written_completion=2*genus+2))
    # Exact inherited-row control.
    P = 18*v**5+81*v**4+120*v**3+60*v*v-1
    Q = (v+1)*(v+2)
    gamma = S.Matrix([[-y*P.subs(v, x*y), x*Q.subs(v, x*y)]])
    C0 = S.Matrix([
        [(3*v*v+6*v+2)/2, x*x*(2*v+3)],
        [3*y*y*(18*v**4+63*v**3+69*v*v+21*v-1)/2,
         36*v**5+108*v**4+87*v**3+3*v*v-3*v+1]
    ]).subs(v, x*y)
    zero('retained normalizer determinant', C0.det()-1)
    zero('retained normalizer changes row to linear pencil', gamma*C0-ell)
    branch = v*P*Q
    check('inherited pencil branch degree eight', S.degree(branch, v) == 8)
    check('inherited pencil squarefree', S.gcd(branch, S.diff(branch, v), v) == 1)
    zero('inherited row Jacobian is branch-polynomial derivative',
         gamma.jacobian((x, y)).det()-S.diff(branch, v).subs(v, x*y))
    # W=y*P and gamma1+t*gamma2=0 imply W^2=t*v*P*Q.
    zero('inherited hyperelliptic identity on pencil',
         (y*P.subs(v, x*y))**2-t*branch.subs(v, x*y)
         +y*P.subs(v, x*y)*(gamma[0]+t*gamma[1]))
    # Euler stratification and completion arithmetic, using symbolic sizes.
    d, k, total_deficit, gg, rr = S.symbols('d k total_deficit gg rr')
    # k nonzero exceptional values, sum of their deficits=total_deficit.
    euler = d*(1-(k+1))+1+k*d-total_deficit
    zero('singleton Euler cancellation is independent of degree', euler-(1-total_deficit))
    zero('Euler completion rearranges to deficit=2g+r-1',
         (1-(2*gg+rr-1))-(2-2*gg-rr))
    # Independent nonproper calibration on the retained F1 source-A3 map.
    z = S.symbols('z')
    D = z**3+z-2
    q = z**3-z+4
    radial = 2*z*z/D
    zero('F1 radial slice derivative', S.diff(radial, z)+2*z*q/D**2)
    check('F1 three distinct poles', S.gcd(D, S.diff(D, z), z) == 1)
    check('F1 three additional distinct critical points', S.gcd(q, S.diff(q, z), z) == 1)
    check('F1 seven punctures are disjoint', S.gcd(D, z*q, z) == 1 and q.subs(z, 0) != 0)
    R = 28*t**3-18*t*t-t+16
    zero('F1 nonzero critical value resultant', S.resultant(q, t*D-2*z*z, z)+8*R)
    check('F1 three distinct nonzero critical values',
          R.subs(t, 0) != 0 and S.gcd(R, S.diff(R, t), t) == 1)
    X = D/(z*q)
    Y = 2*z*(z-3)/D
    Z = S.cancel((2*X-3*X**2*Y-radial)/X**3)
    u = S.cancel(1+X*Y)
    values = [u**3*Z+Y**2*u*(4+3*X*Y),
              Y+3*X*u**2*Z+3*X*Y**2*(4+3*X*Y),
              2*X-3*X**2*Y-X**3*Z]
    zero('F1 complete rational reconstruction on radial curve', S.Matrix(values)-S.ones(3, 1)*radial)
    zero('F1 rational root coordinate inverse', Y+1/X-z)
    check('F1 infinity is retained source origin', all(S.limit(e, z, S.oo) == 0 for e in [X, Y, Z]))
    check('F1 zero critical point is a pole of X', S.limit(z*X, z, 0) != 0)
    check('F1 pole roots are poles of Y', S.gcd(D, 2*z*(z-3), z) == 1)
    check('F1 extra critical roots are poles of X', S.gcd(q, D, z) == 1)
    check('F1 Euler and nonzero-deficit budgets agree', 3-(2+3*2) == 2-7 and 3*2 == 2*0+7-1)
    manifest_path = HERE/'sources/manifest.json'
    manifest = json.loads(manifest_path.read_text())
    for item in manifest['sources']:
        path = HERE/'sources'/item['file']
        check('primary pin '+item['file'], hashlib.sha256(path.read_bytes()).hexdigest() == item['sha256'])
    result = dict(
        as_of='2026-11-13', outcome='PASS', checks=checks, check_count=len(checks),
        runtime_seconds=round(time.monotonic()-start, 6),
        checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        python=platform.python_version(), sympy=S.__version__,
        controls=controls,
        F1_radial_calibration=dict(source_dimension=3, cover_degree=3, genus=0,
                                  punctures=7, zero_deficit=2,
                                  nonzero_critical_value_count=3,
                                  each_nonzero_deficit=2, Euler_characteristic=-5),
        finite_representatives_are_not_the_universal_proof=True,
        global_geometry='Written blowup, Stein-factor, normalization, DVR and Euler arguments; not certified by this checker',
        Keller_counterexample_found=False, unrestricted_pointed_Keller_gate='OPEN',
        formal_verification=False, external_review=False,
    )
    (HERE/'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ('outcome', 'check_count', 'runtime_seconds')}))


if __name__ == '__main__':
    main()
