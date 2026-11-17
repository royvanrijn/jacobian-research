"""Exact local jet, denominator and rational-boundary certificates.

This does not certify Tang's complete HC4 classification or the JC2 conjecture.
Run from the repository root with .venv/bin/python, with a 60-second worker cap.
"""
from __future__ import annotations

import hashlib
from itertools import combinations_with_replacement, permutations
import json
from pathlib import Path
import subprocess
import sys
import sympy as S

HERE = Path(__file__).resolve().parent


def zero(value):
    entries = list(value) if isinstance(value, S.MatrixBase) else [value]
    assert all(S.cancel(v) == 0 for v in entries)


def jet_certificate():
    pairs = list(combinations_with_replacement(range(3), 2))
    K, T, variables = {}, {}, []
    for i in range(3):
        for j, k in pairs:
            symbol = S.Symbol(f'K{i+1}{j+1}{k+1}')
            K[i, j, k] = K[i, k, j] = symbol
            variables.append(symbol)
    for i, j, k in combinations_with_replacement(range(3), 3):
        symbol = S.Symbol(f'T{i+1}{j+1}{k+1}')
        for triple in set(permutations((i, j, k))):
            T[triple] = symbol
        variables.append(symbol)
    J = S.Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
    R = S.Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    equations, labels = [], []
    for ell in range(3):
        dN = S.Matrix(3, 3, lambda i, j: K[i, j, ell])
        dH = S.Matrix(3, 3, lambda i, j: T[i, j, ell])
        for name, expression in [('trace', S.trace(dN)),
                                  ('trace-square/2', S.trace(J*dN)),
                                  ('trace-cube/3', S.trace(J**2*dN))]:
            labels.append(f'd{ell+1} {name}')
            equations.append(expression)
        M = dH*J + R*dN
        for i, j in combinations_with_replacement(range(3), 2):
            if i != j:
                labels.append(f'd{ell+1} skew(HN) {i+1}{j+1}')
                equations.append(M[i, j] - M[j, i])
    A, constant = S.linear_eq_to_matrix(equations, variables)
    zero(constant)
    target = K[2, 1, 1] - K[2, 2, 0]
    target_vector = S.Matrix([[S.diff(target, v) for v in variables]])
    family = next(iter(S.linsolve((A.T, target_vector.T))))
    free = set().union(*(entry.free_symbols for entry in family))
    coefficients = S.Matrix([[entry.subs(dict.fromkeys(free, 0)) for entry in family]])
    zero(coefficients*A-target_vector)
    def fraction(q):
        return [int(S.numer(q)), int(S.denom(q))]
    certificate = dict(variables=list(map(str, variables)), labels=labels,
                       rows=[[fraction(q) for q in A.row(i)] for i in range(A.rows)],
                       weights=[fraction(q) for q in coefficients],
                       target=[fraction(q) for q in target_vector],
                       interpretation='K322-K331 is the normalized Frobenius coefficient of im(DV)',
                       scope='Universal first-order identities at regular nilpotent J3 and anti-diagonal Hessian metric; no degree truncation')
    (HERE/'jet-certificate.json').write_text(json.dumps(certificate, indent=2)+'\n')
    return dict(equations=A.rows, variables=A.cols, rank=int(A.rank()),
                certificate_nonzero_weights={labels[i]: str(q) for i, q in enumerate(coefficients) if q})


def rational_family():
    # k0,k1,k2 stand for k(t), k'(t), k''(t). These identities impose no bound on k.
    t, y, s, k0, k1, k2 = S.symbols('t y s k0 k1 k2')
    a = S.Matrix([1, -t, 0])
    dt = a/y
    Ht = S.Matrix([[0, -1/y**2, 0], [-1/y**2, 2*t/y**2, 0], [0, 0, 0]])
    ds = S.Matrix([-k1, -k0+t*k1, 1])
    H = Ht + ds*ds.T - s*k2/y*(a*a.T)
    direction = S.Matrix([t, 1, k0])
    N = direction*dt.T
    zero(N*N)
    zero(H*direction+dt/y)
    zero(H*N+(dt*dt.T)/y)
    zero(H.det()+1/y**4)

    x, y, z = S.symbols('x y z')
    t = x/y
    V = S.Matrix([t**2/2, t, t**3/3])
    f = t+(z-x*x/y)**2/2
    N = V.jacobian((x, y, z))
    H = S.hessian(f, (x, y, z))
    zero(H*N-(H*N).T)
    zero(N*N)
    zero(H.det()+1/y**4)
    F = S.Matrix([x, y, z])+V
    zero(F.jacobian((x, y, z)).det()-1)
    p = {x: 0, y: 1, z: 0}
    q = {x: -2, y: -1, z: -S.Rational(8, 3)}
    zero(F.subs(p)-F.subs(q))
    Y1, Y2, Y3, r = S.symbols('X Y Z r')
    inverse = {x: r*(Y2-r), y: r, z: Y3-(Y2-r)**3/3}
    result = F.subs(inverse, simultaneous=True)
    zero(result-S.Matrix([(Y2**2-r**2)/2, Y2, Y3]))
    # The forward map identifies r with y and Delta with y^2.
    zero(F[1]**2-2*F[0]-y**2)
    theta = S.Symbol('theta')
    relation_matrix = S.Matrix([[S.expand(v).coeff(theta, j) for j in range(4)]
                               for v in (1, theta**2/2, theta, theta**3/3)])
    assert relation_matrix.det() != 0
    return dict(metric_determinant='-1/y^4', nilpotence='N^2=0',
                generic_degree=2, source='Spec Q[x,y,z,y^-1]',
                target='Spec Q[X,Y,Z,(Y^2-2X)^-1]',
                collision=[[0, 1, '0'], [-2, -1, '-8/3']],
                universal_family='V=(t^2/2,t,K(t)), f=t+(z-y*K_prime(t))^2/2, t=x/y; deg K>=3',
                nonplanarity='Independence of 1,t,t^2,K(t); concrete coefficient determinant checked',
                scope='Rational metric and punctured affine source; not a polynomial HC4 or JC2 counterexample')


def denominator_identity():
    q, p = S.symbols('q p')
    M = S.Matrix(3, 3, lambda i, j: S.Symbol(f'm{i}{j}'))
    g = S.Matrix(S.symbols('g0:3'))
    # Matrix determinant lemma with no assumption that M is invertible.
    zero((q*M+2*p*g*g.T).det()-q**2*(q*M.det()+2*p*(g.T*M.adjugate()*g)[0]))
    x, y, z, w = S.symbols('x y z w')
    h = (x*x+y*y+z*z)/(2*(1+w))
    full = S.hessian(h, (x, y, z, w))
    assert full.rank() == 3
    f = h.subs(w, 0)
    A = S.hessian(f, (x, y, z))
    v = A.inv()*S.Matrix([S.diff(h, u, w).subs(w, 0) for u in (x, y, z)])
    zero(v+S.Matrix([x, y, z]))
    zero(S.Matrix([S.diff(h, u) for u in (x, y, z)]).subs(
        {x: x-w*v[0], y: y-w*v[1], z: z-w*v[2]}, simultaneous=True)-S.Matrix([x, y, z]))
    D = (S.eye(3)-w*v.jacobian((x, y, z))).det()
    zero(D-(1+w)**3)
    zero(S.hessian(h, (x, y, z)).det()*D-1)
    grad = S.Matrix([S.diff(h, u) for u in (x, y, z, w)])
    zero(sum(u*u for u in grad[:3])+2*grad[3])
    return dict(transport='det(A)=det(Hess_X h(X-w*v,w))*det(I-w*Dv)',
                denominator_bound='det(I-w*Dv) divides Q(X-w*v,w)^3 over K(X)[w] for h=P/Q in four variables',
                general_bound='In n transverse variables the exponent is n, and is sharp',
                intermediate_hessian_denominator_bound='2*n+1; support plus degree(D)<=n improves this to exponent n',
                rational_failure=dict(potential='(x^2+y^2+z^2)/(2*(1+w))',
                                      slice_metric='I3', normal_field='v=-X',
                                      characteristic_determinant='(1+w)^3', nilpotent=False,
                                      gradient_relation='U1^2+U2^2+U3^2+2*U4=0; no nonzero translation direction'),
                proof_boundary='The universal divisibility implication is proved in README; exact matrix identity and rational example are replayed here')


def phase_minor_identities():
    c, alpha, FL, FLL = S.symbols('c alpha FL FLL')
    ell = S.Matrix(S.symbols('l0:3'))
    ell_prime = S.Matrix(S.symbols('lp0:3'))
    b = alpha*ell+FL*ell_prime
    H = S.Matrix([[c, *list(b)]]).col_join(b.row_join(FLL*ell*ell.T))
    for i, j in ((0, 1), (0, 2), (1, 2)):
        zero(H.extract([0, i+1, j+1], [0, i+1, j+1]).det()
             +FLL*FL**2*(ell_prime[i]*ell[j]-ell_prime[j]*ell[i])**2)
    r, rp, Bp = S.symbols('r rp Bp')
    ell = S.Matrix([r, 0, 1])
    b = alpha*ell+FL*S.Matrix([rp, 0, 0])+Bp*S.Matrix([0, 1, 0])
    H = S.Matrix([[c, *list(b)]]).col_join(b.row_join(FLL*ell*ell.T))
    normal = S.Matrix([0, -Bp, FL*rp, r*Bp])
    zero(H*normal)
    zero(H.adjugate()+FLL*normal*normal.T)
    return dict(three_phase_cofactors='PASS for symbolic phase and first derivative',
                two_phase_affine_partner_kernel_and_adjugate='PASS',
                wall_exhaustion_and_termination_verified=False)


def nonintegrable_control():
    x, y, z = S.symbols('x y z')
    u = y-x*x
    V = S.Matrix([u, z+2*x*u, -u*u])
    N = V.jacobian((x, y, z))
    omega = S.Matrix([[2*u, 0, 1]])
    zero(N**3)
    assert N.rank() == 2
    zero(omega*N)
    curl = S.Matrix([S.diff(omega[2], y)-S.diff(omega[1], z),
                     S.diff(omega[0], z)-S.diff(omega[2], x),
                     S.diff(omega[1], x)-S.diff(omega[0], y)])
    coefficient = S.expand((omega*curl)[0])
    assert coefficient == -2
    return dict(V=list(map(str, V)), nilpotence='N^3=0', image_annihilator=list(map(str, omega)),
                frobenius_coefficient=str(coefficient),
                consequence='No generically nonsingular Hessian metric, even rational, can symmetrize DV')


def main():
    source = json.loads((HERE/'source-receipt.json').read_text())
    assert hashlib.sha256((HERE/'source-v1.pdf').read_bytes()).hexdigest() == source['sha256']
    result = dict(jet=jet_certificate(), rational_family=rational_family(),
                  denominator=denominator_identity(), phase_minors=phase_minor_identities(),
                  control=nonintegrable_control(),
                  full_Tang_normal_form_independently_verified=False,
                  HC4_or_JC2_resolved=False, novelty_priority_claimed=False)
    independent = subprocess.run([sys.executable, str(HERE/'check_jet.py')],
                                  capture_output=True, text=True, timeout=15, check=True)
    result['independent_jet_check'] = independent.stdout.strip()
    (HERE/'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
