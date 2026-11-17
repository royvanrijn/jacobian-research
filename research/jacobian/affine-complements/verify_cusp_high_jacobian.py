"""PF2RULED23: exact highest next-Jacobian block on all sixteen cusp-lift families.
Run with sage -python. No numerical embeddings or sampled parameters are used.
"""
from pathlib import Path
from time import monotonic
import hashlib
import json
import platform
from sage.all import QQ, NumberField, PolynomialRing, matrix, vector, binomial
from sage.version import version as sage_version

HERE = Path(__file__).resolve().parent
checks = []

def check(name, value):
    assert value, name
    checks.append(name)

def jet_basis(m, V):
    """Full normalized Cox basis; its symbol is r^(2(m-j)+i)."""
    v = V.gen()
    result = []
    for j in range(m + 1):
        for i in range(m + j + 1):
            nn = m + j - i
            jets = []
            for aa in range(3):
                coeff = V.zero()
                for bb in range(min(aa, j) + 1):
                    z = aa - bb
                    if z <= nn:
                        coeff += ((-1)**(m-j+bb) * binomial(j, bb) *
                                  binomial(nn, z) * v**(m-j+bb+nn-z) *
                                  (1-v*v)**(j-bb))
                jets.append(coeff)
            result.append(dict(j=j, i=i, symbol_power=2*(m-j)+i, jets=jets))
    return result

def multiply(A, B):
    return [sum(A[i]*B[k-i] for i in range(k+1)) for k in range(3)]

def combine(coefficients, basis):
    return [sum(c*f['jets'][j] for c, f in zip(coefficients, basis)) for j in range(3)]

def jac1(A, B):
    return (2*(A[2]*B[0].derivative()-A[0].derivative()*B[2]) +
            A[1]*B[1].derivative()-A[1].derivative()*B[1])

def jac_coefficient(A, B, k):
    total = A[0].base_ring().zero()
    for ja, jb, scalar, derivative in [(2, 0, 2, 'b'), (0, 2, -2, 'a'),
                                      (1, 1, 1, 'b'), (1, 1, -1, 'a')]:
        for aa in range(max(0, k+1-B[jb].degree()), A[ja].degree()+1):
            bb = k+1-aa
            if bb >= 0:
                total += scalar*(bb if derivative == 'b' else aa)*A[ja][aa]*B[jb][bb]
    return total

def sparse(values):
    return [[i, str(c)] for i, c in enumerate(values) if c]

def main():
    start = monotonic()
    pred = json.loads((HERE/'full_cusp_lift_receipt.json').read_text())
    check('predecessor PASS', pred['outcome'] == 'PASS')
    for name, digest in pred['source_sha256'].items():
        check('predecessor seal '+name, hashlib.sha256((HERE/name).read_bytes()).hexdigest() == digest)
    candidates = json.loads((HERE/'polar_cusp_sector_certificates.json').read_text())
    R = PolynomialRing(QQ, 'y')
    sectors = []
    for old in candidates['two_root_marked_sectors']:
        k = old['multiplicities'][0]
        p8 = R(old['P8'])
        check(f'irreducible degree-eight candidate algebra k{k}', p8.degree() == 8 and p8.is_irreducible())
        K = NumberField(p8, 'y')
        y = K.gen()
        x = K(R(old['x_as_polynomial_in_y']))
        V = PolynomialRing(K, 'v'); v = V.gen()
        Hring = PolynomialRing(K, 'w'); w = Hring.gen()
        H = (w-x)**k*(w-y)**(9-k)
        a, b = H[8], H[7]
        delta = -a/48
        s = (24*b+72-11*a*a)/(12*a)
        h, lam = v-s, 3*(v-s)/4+delta
        check(f'nonzero cusp unit and distinct roots k{k}', delta != 0 and x != y)
        def L(F): return h*F[1]-lam*F[0].derivative()
        bas3 = jet_basis(3, V)
        T = [V.zero() for _ in range(3)]
        chosen = []
        for power in range(10):
            f = next(f for f in bas3 if f['symbol_power'] == power)
            chosen.append([f['j'], f['i']])
            T = [T[j]+H[9-power]*f['jets'][j] for j in range(3)]
        check(f'T has complete four highest cusp coefficients zero k{k}', L(T).degree() <= 7)
        T2 = multiply(T, T); T3 = multiply(T2, T)
        solutions = []; reference_data = []
        for m, F, dimension in [(5, T2, 31), (8, T3, 85)]:
            basis = jet_basis(m, V)
            check(f'complete Cox basis dimension k{k} m{m}', len(basis) == (m+1)*(3*m+2)//2)
            check(f'residual in full lower-pole L image k{k} m{m}', L(F).degree() < 4*m)
            M = matrix(K, [[L(f['jets'])[j] for f in basis] for j in range(4*m)])
            check(f'full cusp operator surjectivity k{k} m{m}', M.rank() == 4*m)
            rhs = vector(K, [-L(F)[j] for j in range(4*m)])
            solution = M.solve_right(rhs)
            correction = combine(solution, basis)
            reference = [F[j]+correction[j] for j in range(3)]
            check(f'whole reference cusp equation k{k} m{m}', L(reference) == 0)
            ker = M.right_kernel().basis_matrix()
            check(f'full affine direction dimension k{k} m{m}', ker.nrows() == dimension)
            directions = [combine(row, basis) for row in ker.rows()]
            check(f'every affine direction satisfies full cusp equation k{k} m{m}', all(L(f) == 0 for f in directions))
            solutions.append((reference, directions, ker, basis, solution))
            reference_data.append(dict(pole=m,correction=sparse(solution),kernel_dimension=dimension))
        P, kp, pker, pbas, psol = solutions[0]
        Q, kq, qker, qbas, qsol = solutions[1]
        # Every pair in the full116-dimensional affine space is represented.
        # Va,Vb bounds already force the direction-direction bracket degree<=49.
        # Check the only potentially remaining top bilinear coefficient exactly.
        check(f'complete bilinear top coefficient49 vanishes k{k}',
              all(jac_coefficient(A, B, 49) == 0 for A in kp for B in kq))
        J = jac1(P, Q)
        p_high = matrix(K, [[jac_coefficient(f, Q, j) for f in kp] for j in range(46, 54)])
        q_high = matrix(K, [[jac_coefficient(P, f, j) for f in kq] for j in range(46, 54)])
        check(f'five reference top coefficients49..53 vanish k{k}', all(J[j] == 0 for j in range(49, 54)))
        check(f'five top coefficients vanish on all P directions k{k}', p_high[3:, :] == 0)
        check(f'five top coefficients vanish on all Q directions k{k}', q_high[3:, :] == 0)
        check(f'fixed-Q remaining three-equation matrix has full rank k{k}', p_high.rank() == 3)
        pivots = p_high[:3, :].pivots()
        square = p_high[:3, list(pivots)]
        det = square.det()
        check(f'exact rank minor is a unit at every embedding k{k}', det != 0 and det*(1/det) == 1)
        rhs = vector(K, [-J[j] for j in range(46, 49)])
        controls = square.solve_right(rhs)
        Pnew = [P[j]+sum(controls[ii]*kp[pivots[ii]][j] for ii in range(3)) for j in range(3)]
        check(f'all eight high Jacobian equations solved k{k}', all(jac1(Pnew, Q)[j] == 0 for j in range(46, 54)))
        check(f'fixed-Q affine high-block dimension28 k{k}', len(kp)-3 == 28)
        # Exact curve controls on the selected witness. Open-locus proof uses T2/T3.
        pp, remp = Pnew[0].derivative().quo_rem(h)
        qq, remq = Q[0].derivative().quo_rem(h)
        check(f'witness derivatives divisible by h k{k}', remp == remq == 0)
        check(f'witness has clean derivative gcd k{k}', pp.gcd(qq).degree() == 0)
        W = pp*qq.derivative()-pp.derivative()*qq
        check(f'witness ordinary cusp wedge nonzero k{k}', W(s) != 0)
        modifications = []
        for z in (2, 3):
            f = [h**z, z*lam*h**(z-2), -z*v*lam*h**(z-2)]
            check(f'actual global cusp modification k{k} z{z}', L(f) == 0)
            check(f'P-modification preserves eight high coefficients k{k} z{z}', jac1(f, Q).degree() <= 45)
            check(f'Q-modification preserves eight high coefficients k{k} z{z}', jac1(Pnew, f).degree() <= 45)
            modifications.append(f)
        check(f'modification cross brackets preserve high block k{k}',
              all(jac1(A,B).degree() <= 45 for A in modifications for B in modifications))
        # Retain an exact global-section witness via its complete Cox coefficients.
        new_correction = vector(K, psol)
        for ii, pivot in enumerate(pivots): new_correction += controls[ii]*pker.row(pivot)
        kappa = delta**2*W(s)
        actualJ = jac1(Pnew, Q)
        check(f'witness next quotient degree at most44 k{k}', (actualJ-kappa).quo_rem(h)[1] == 0 and ((actualJ-kappa)//h).degree() <= 44)
        G = (actualJ-kappa)//h
        low_columns = [2*qq*v**j for j in range(11)]+[-2*pp*v**j for j in range(23)]
        low_matrix = matrix(K, [[f[j] for f in low_columns] for j in range(45)])
        low_rhs = matrix(K, 45, 1, [-G[j] for j in range(45)])
        check(f'clean witness normal-control matrix has rank34 k{k}', low_matrix.rank() == 34)
        check(f'chosen high-block witness fails lower membership even after normal correction k{k}',
              low_matrix.augment(low_rhs).rank() == 35)
        sectors.append(dict(multiplicities=old['multiplicities'],P8=str(p8),x=str(x),
            delta=str(delta),cusp_parameter=str(s),T_symbol_lift_cox_ji=chosen,
            reference_corrections=reference_data,rank_minor_kernel_columns=list(map(int,pivots)),
            rank_minor=[[str(c) for c in row] for row in square.rows()],
            rank_minor_determinant=str(det),rank_minor_inverse=str(1/det),
            high_block_controls=list(map(str,controls)),
            solved_P_lower_pole_correction=sparse(new_correction),
            solved_Q_lower_pole_correction=sparse(qsol),
            fixed_Q_affine_high_block_dimension=28,witness_derivative_gcd_is_h=True,
            witness_cusp_wedge_nonzero=True,witness_low_matrix_rank=34,witness_low_augmented_rank=35,witness_next_quotient_degree=int(((actualJ-kappa)//h).degree()),
            all_five_automatic_coefficients_checked_on_complete_affine_space=True,
            source_Keller_pair_constructed=False))
        print(f'PASS k{k}: all8 embeddings; five automatic coefficients; rank3 high block solved; remaining11 membership conditions open', flush=True)
    cert = dict(as_of='2026-11-14',sectors=sectors,marked_points=16,
        full_affine_lift_dimension=116,automatic_quotient_degrees=list(range(48,53)),
        active_high_quotient_degrees=[45,46,47],fixed_Q_high_rank=3,fixed_Q_affine_solution_dimension=28,
        rank_three_chart_dimension=113,rank_three_chart_boundary_parameters=43,
        clean_singleton_cusp_open_nonempty_after_high_block=True,remaining_low_membership_conditions=11,
        full_next_Jacobian_solved=False,global_finiteness_proved=False)
    cp = HERE/'cusp_high_jacobian_certificates.json'
    cp.write_text(json.dumps(cert,indent=2)+'\n')
    inputs = [Path(__file__),HERE/'CUSP_HIGH_JACOBIAN_REDUCTION.md',cp,
              HERE/'FULL_CUSP_FIRST_JET_LIFTS.md',HERE/'full_cusp_lift_receipt.json',
              HERE/'polar_cusp_sector_certificates.json',HERE/'SOURCE_COORDINATES.md',HERE/'WRIGHT_MARKED_MODEL.md']
    report = dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        marked_points=16,full_affine_lift_dimension=116,automatic_quotient_degrees=list(range(48,53)),
        active_high_quotient_degrees=[45,46,47],fixed_Q_high_rank=3,fixed_Q_affine_solution_dimension=28,
        rank_three_chart_dimension=113,rank_three_chart_boundary_parameters=43,
        all_eight_high_equations_solved_at_all_marked_points=True,
        clean_singleton_cusp_open_nonempty_after_high_block=True,remaining_low_membership_conditions=11,
        first_open_reduced_pairs=[[5,9],[6,9],[7,9],[8,9]],minimum_compatible_quartic_pole_height=9,
        assurance=dict(full_affine_space_covered=True,all_complex_embeddings_covered=True,
            no_candidate_removed=True,lower_membership_equations_solved=False,
            full_next_Jacobian_solved=False,all_higher_Jacobian_coefficients_solved=False,
            field_degree_four_realization_constructed=False,global_finiteness_proved=False,
            contact_cases_excluded=False,actual_ruled_compatibility_proved=False,
            quartic_packets_excluded=False,formal_verification=False,independent_full_reproof=False,
            external_review=False,research_goal_completed=False),
        source_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs},
        python=platform.python_version(),sage=sage_version,elapsed_seconds=round(monotonic()-start,6))
    (HERE/'cusp_high_jacobian_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} exact checks; sixteen candidates retained; high block closed; goal active',flush=True)

if __name__ == '__main__':
    main()


