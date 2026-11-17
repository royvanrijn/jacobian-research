"""Focused exact transcriptions; the universal proofs are in the linked note.

Writes this packet only. Does not rerun old Gaussian/Oaku registrations.
"""
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
from math import factorial, prod, comb
from pathlib import Path
import json
import time
import gaussian_seed_filtration as G

s = G.s
HERE = Path(__file__).resolve().parent


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def wick(e):
    p, q, az, bw, c, d, h = e
    if p != q or az != bw or c != d or h % 2:
        return 0
    return factorial(p)*factorial(az)*factorial(c)*prod(range(1, h, 2))


def direct_moments(phase, markers, N):
    variables = (G.u, G.a, G.z, G.w, G.y, G.v, G.T)
    terms = {e: Fraction(c) for e, c in s.Poly(phase, variables).terms()}
    power = {(0,)*7: Fraction(1)}
    result = [[] for _ in markers]
    for n in range(N+1):
        for i, alpha in enumerate(markers):
            result[i].append(sum(c*wick(tuple(e[j]+alpha[j] for j in range(7)))
                                 for e, c in power.items()))
        if n < N:
            next_power = defaultdict(Fraction)
            for e, c in power.items():
                for f, b in terms.items():
                    next_power[tuple(e[j]+f[j] for j in range(7))] += c*b
            power = {e:c for e,c in next_power.items() if c}
    return result


def coeffs(expr, N):
    expression = s.series(expr, G.t, 0, N+1).removeO().expand()
    return [expression.coeff(G.t, n) for n in range(N+1)]


def check():
    start = time.monotonic()
    markers = [(0,)*7, (1,1,0,0,0,0,0), (0,2,0,0,1,0,0),
               (0,0,1,1,0,0,0), (0,0,0,0,0,1,0), (0,0,0,0,2,1,0),
               (0,0,0,0,1,2,2), (0,0,0,0,0,0,2), (0,0,0,0,0,0,1),
               (0,0,1,0,0,0,0), (2,1,0,0,2,1,0), (0,0,2,2,1,0,0)]
    packet = dict(schema_version=1, universal_claims_use_written_proofs=True,
                  independent_replay=False, formal_verification=False, external_review=False,
                  computed_Galois_groups=False, actual_image_kernel_computed=False,
                  original_marker_controls=[], filtration_controls=[], power_controls=[],
                  transfer_controls=[], recurrence_controls=[], airy_controls=[], rejections=[])
    for f, M, five, N in [(1-G.t, 0, True, 4), (1-G.t, 2, True, 6),
                           (1-G.t-G.t**2, 1, False, 4)]:
        phase = G.phase(f, M, five_variables=five)
        values = direct_moments(phase, markers, N)
        for alpha, moments in zip(markers, values):
            if five and (alpha[0] or alpha[1]):
                # U,A are absent independent unused Gaussians in this control.
                continue
            period = G.marker_period(f, M, alpha)
            assert coeffs(period, N) == [s.Rational(m.numerator,m.denominator)/s.factorial(n)
                                        for n,m in enumerate(moments)], (f,M,alpha)
            quotient = s.cancel(period/f**(M-sum(alpha)))
            assert s.fraction(quotient)[1] == 1
            packet['original_marker_controls'].append(dict(f=str(f), M=M, marker=alpha,
                period=str(period), moments=[str(m) for m in moments], N=N))
        for d in range(M+3):
            period = sum(s.binomial(d,c)*G.marker_period(f,M,(0,0,0,0,c,0,0))
                         for c in range(d+1))
            assert s.cancel(period-f**(M-d)) == 0
            packet['filtration_controls'].append(dict(f=str(f), M=M, d=d,
                                                      generator=str(s.factor(period))))
    # Every degree-at-most-two marker in all seven original variables.
    from itertools import combinations_with_replacement
    all_small = {(0,)*7}
    for d in [1,2]:
        for indices in combinations_with_replacement(range(7),d):
            e=[0]*7
            for i in indices:e[i]+=1
            all_small.add(tuple(e))
    for alpha in sorted(all_small):
        period=G.marker_period(1-G.t,2,alpha)
        assert s.fraction(s.cancel(period))[1] == 1
        assert s.fraction(s.cancel(period/(1-G.t)**(2-sum(alpha))))[1] == 1
    packet['all_degree_at_most_two_markers_checked']=len(all_small)
    # Pure powers use original Wick contractions of P^K, without expanding Q^n twice.
    for M,K in [(0,2),(1,2),(2,3)]:
        N=6
        coords=[(0,0,0,0,c,0,0) for c in range(M+2)]
        moments=direct_moments(G.phase(1-G.t,M,five_variables=True),coords,N)
        assert all(moments[0][n] == 0 for n in range(K,N+1,K))
        bad=[sum(comb(M+1,c)*moments[c][n] for c in range(M+2))
             for n in range(0,N+1,K)]
        assert bad == [factorial(n) for n in range(0,N+1,K)]
        data=G.power_data(M,K)
        assert data['first_bad_marker_degree']==M+1
        packet['power_controls'].append(dict(M=M,K=K,original_power_indices=list(range(0,N+1,K)),
            unit_moments=[str(moments[0][n]) for n in range(0,N+1,K)],
            bad_marker=str(data['bad_marker']),bad_moments=[str(x) for x in bad],
            full_common_image='D_t*h_'+str(K),actual_rank=K))
    # A mixed polynomial/pole packet checks exact finite-coefficient extraction.
    for K in [2,3,4]:
        old=3+2*G.t**K+5*G.t**(K+1)+7/(1-G.t)**2-4/(1-G.t)**3
        polynomial,op=G.power_transfer(old,K)
        hc=G.h_coefficients(K,6)
        actual=coeffs(old,K*6)
        wanted=[actual[K*n]*s.factorial(K*n)/s.factorial(n) for n in range(7)]
        got=[s.expand(polynomial).coeff(G.t,n)+op.subs(G.theta,n)*hc[n] for n in range(7)]
        assert got == wanted
        packet['transfer_controls'].append(dict(K=K,old_period=str(old),polynomial=str(polynomial),
                                                  Euler_operator=str(op),coefficients=[str(x) for x in got]))
        b=K*s.prod(K*G.theta+j for j in range(1,K))
        assert all(hc[n+1]==b.subs(G.theta,n)*hc[n] for n in range(6))
        packet['recurrence_controls'].append(dict(K=K,operator=str(G.quotient_operator(K)),
                                                     coefficients=[str(x) for x in hc]))
        x=s.Symbol('x');F=s.Function('F')(x);value=x*F
        # Substitute K*theta_t=-theta_x into the ORIGINAL A_K, then gauge y=xF.
        for j in range(1,K):value=-x*s.diff(value,x)+j*value
        transformed=x*F-(-1)**(K-1)*x**(-K)*value
        assert s.simplify(transformed+(s.diff(F,x,K-1)-x*F)) == 0
        packet['airy_controls'].append(dict(K=K,cover=str(G.airy_cover(K)['t']),
                                              gauge='y=x*u',identity='A_K[x*u]=-(D_x^(K-1)-x)*u'))
    # Arbitrary polynomial seed and an ordinary nonreduced specialization.
    f=1+2*G.t+3*G.t**2
    assert G.marker_period(f,1,(0,)*7) == f
    packet['arbitrary_seed_control']=dict(f=str(f),pure_moments=[1,2,6,0,0])
    e=s.Symbol('e');f=1-e*G.t
    inverse=sum(s.binomial(2,c)*G.marker_period(f,1,(0,0,0,0,c,0,0)) for c in range(3))
    assert s.expand(s.series(inverse,e,0,2).removeO()-1-e*G.t) == 0
    packet['nilpotent_control']=dict(ring='Q[e]/(e^2)',seed='1-e*t',inverse='1+e*t',
                                     safe=True,pure=False,uses_written_all_marker_unit_criterion=True)
    for call in [lambda:G.phase(2-G.t,0),lambda:G.phase(1-G.t**2,1,five_variables=True),
                 lambda:G.power_data(2,2),lambda:G.marker_period(1-G.t,0,(0,)*6),
                 lambda:G.power_transfer(1/(1+G.t),2),lambda:G.airy_cover(1)]:
        try:call()
        except ValueError:packet['rejections'].append(True)
        else:raise AssertionError('Invalid input accepted')
    imports=['ALGEBRAIC_FIXED_POINT_LIFTS_AND_UNIVERSAL_GAUSSIAN_IMAGES.md',
             'OAKU_GRAPH_CERTIFICATES_AND_ACTUAL_POLYNOMIAL_DECISION.md',
             'BERNSTEIN_CLASSICAL_CREDIT_AND_FORMAL_EXTENSIONS.md']
    runtime=dict(python=__import__('sys').version,sympy=s.__version__,
                 sealed_imports_sha256={n:digest(HERE/n) for n in imports})
    (HERE/'gaussian_seed_filtration_runtime.json').write_text(json.dumps(runtime,indent=2)+'\n')
    certificate=HERE/'gaussian_seed_filtration_certificate.json'
    certificate.write_text(json.dumps(packet,indent=2)+'\n')
    controls={k:len(packet[k]) for k in ['original_marker_controls','filtration_controls','power_controls',
                                        'transfer_controls','recurrence_controls','airy_controls','rejections']}
    receipt=dict(completed=True,elapsed_seconds=round(time.monotonic()-start,3),exact_controls=controls,
                 all_degree_at_most_two_markers_checked=len(all_small),
                 canonical_proof_sha256=digest(HERE/'POLYNOMIAL_SEEDS_AND_PURE_MARKER_DELAY.md'),
                 certificate_sha256=digest(certificate),producer_sha256=digest(HERE/'gaussian_seed_filtration.py'),
                 checker_sha256=digest(Path(__file__)),runtime_sha256=digest(HERE/'gaussian_seed_filtration_runtime.json'),
                 universal_claims_use_written_proofs=True,independent_replay=False,
                 formal_verification=False,external_review=False,computed_Galois_groups=False)
    (HERE/'gaussian_seed_filtration_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__ == '__main__':
    check()
