"""Directed finite-block controls for RADIAL_ANGULAR_MODE_ISOLATION.

This certifies one two-dimensional radial block at rho=99/100. It does
not count every adverse radial mode or prove the all-partition inequality.
The all-rho compactness, global polynomial extraction and finite angular
endpoint reduction use the written proof.
"""
from fractions import Fraction as Q
import json
from flint import arb, ctx
from check_three_cell_coupled_radial import quadrant, first_mode, enclose


def certificate():
    ctx.prec = 160
    rho = arb(99)/100
    d = 1-rho*rho
    exponents = (0, 1, 4, 16, 64, 256, 1024)
    f = (Q(-21, 34), Q(1), Q(0), Q(2), Q(0), Q(0), Q(0))
    gc = tuple(Q(n, 1000) for n in (27, 37, 74, 181, 145, 1000))
    g = (-sum((c/Q(1+a) for a, c in zip(exponents[1:], gc)), Q(0)),)+gc
    assert all(sum((c/Q(1+a) for a, c in zip(exponents, v)), Q(0)) == 0
               for v in (f, g))
    matrices = [[], []]
    for a in exponents:
        rows = [[], []]
        for b in exponents:
            def pairing(z, diagonal=False, angular_derivative=False):
                prod = 1+(a+b)*d if diagonal else (1+a*d)*(1+b*d)
                if angular_derivative:
                    # DG24 first_mode=j1/2 times integral Q*kappa*exp(-kappa).
                    return 2*first_mode(prod, rho, d, arb(1))
                return d*quadrant(prod, -rho*z)
            n1, nh, nm = (pairing(z) for z in (arb(1), -arb(1)/2, -arb(1)))
            d1, dh, dm = (pairing(z, True) for z in (arb(1), -arb(1)/2, -arb(1)))
            na = pairing(0, angular_derivative=True)
            da = pairing(0, True, True)
            old = (n1+2*nh-nm-d1+dm)/2
            exact_arc = (n1+nh-d1+dh)/2
            new_penalty = (10*(dh-nh)-10*(dm-nm)-4*(da-na))/4
            rows[0].append(old)
            rows[1].append(exact_arc-new_penalty)
        for out, row in zip(matrices, rows):
            out.append(row)
    vectors = [tuple(arb(c.numerator)/c.denominator for c in v) for v in (f, g)]
    blocks = []
    for matrix in matrices:
        def bilinear(v, w):
            return sum((v[i]*w[j]*matrix[i][j] for i in range(7) for j in range(7)), arb(0))
        ff = bilinear(vectors[0], vectors[0])
        fg = bilinear(vectors[0], vectors[1])
        gg = bilinear(vectors[1], vectors[1])
        det = ff*gg-fg*fg
        blocks.append((ff, fg, gg, det))
    off, ofg, ogg, odet = blocks[0]
    nff, nfg, ngg, ndet = blocks[1]
    assert off < 0 and ogg > 0 and odet < 0
    shift = arb(1)/30000
    shifted_det = (nff+shift)*(ngg+shift)-nfg*nfg
    assert nff+shift < 0 and shifted_det > 0
    # The older block has exactly one positive direction. The new block
    # is negative on the entire span, with the stated coefficient margin.
    positive_eigenvalue = (off+ogg+((off-ogg)**2+4*ofg*ofg).sqrt())/2
    assert positive_eigenvalue > 0
    eigenvector_f_over_g = ofg/(positive_eigenvalue-off)
    center = -ofg/off
    radius = (-odet).sqrt()/(-off)
    return {'status': 'pass', 'rho': '99/100', 'arb_precision_bits': 160,
            'arithmetic': 'directed Arb closed Gaussian quadrant and derivative formulas',
            'radial_exponents': list(exponents),
            'basis_coefficients': [[str(c) for c in v] for v in (f, g)],
            'exact_radial_means': ['0', '0'],
            'old_block': {'ff': enclose(off), 'fg': enclose(ofg), 'gg': enclose(ogg),
                          'determinant': enclose(odet),
                          'positive_eigenvalue_in_coefficient_norm': enclose(positive_eigenvalue),
                          'positive_eigenvector_f_over_g': enclose(eigenvector_f_over_g),
                          'positive_cone_f_over_g': [enclose(center-radius), enclose(center+radius)]},
            'two_mode_block': {'ff': enclose(nff), 'fg': enclose(nfg), 'gg': enclose(ngg),
                               'determinant': enclose(ndet),
                               'negative_coefficient_margin': '1/30000',
                               'shifted_determinant': enclose(shifted_det)},
            'positive_control_old_bound_rejected': True,
            'scope': 'The older exact-arc/secant quadratic form has one positive eigenvalue in span{f,g}; the new two-angular-mode form is <=-(c_f^2+c_g^2)/30000 on that entire span. Closed formulas have no radial/angular quadrature or sampled-only gaps. Written uniform Taylor control gives a common small-profile Gaussian comparison for these directions with arbitrary angular shapes. This certificate does not cover other radial directions at ordertwo or finite amplitudes; the written finite-higher-order existence theorem is separate.'}


if __name__ == '__main__':
    print(json.dumps(certificate(), indent=2))
