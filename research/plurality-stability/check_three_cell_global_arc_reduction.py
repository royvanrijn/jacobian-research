"""Directed explicit-instance controls for the global arc obstruction DG14.

The all-rho construction, its comparison with every measurable family of
consecutive arcs, and the interior-profile extension use the written proof.
This program certifies its constants at rho=1/2 and exact inner mass 10^-8.
"""
from fractions import Fraction as Q
from itertools import permutations
import json
import math
from flint import arb, ctx


def enclose(x):
    return [math.nextafter(float(x.lower()), -math.inf),
            math.nextafter(float(x.upper()), math.inf)]


def circular_distance_pi(x):
    x = x % 2
    return min(x, 2-x)


def constants(alpha_exact):
    ctx.prec = 100
    rho = arb(1)/2
    d = 1-rho*rho
    alpha = arb(alpha_exact.numerator)/alpha_exact.denominator
    beta = arb(-2).exp() - (arb(-9)/2).exp()
    available = (arb(-1)/2).exp() - (arb(-25)/32).exp()
    qmin = (-169*rho*rho/(32*d)).exp()/d
    delta = qmin*(rho/d).sinh()/72
    coercivity = 8*rho*(-9*rho*rho/d).exp()/(arb.pi()**2*d*d)
    lipschitz = 2*(15*rho/(4*d)).sinh()/(3*arb.pi()*d)
    cross_error = arb(2).sqrt()*arb.pi()*lipschitz
    inner_bound = (25*rho/(16*(1+rho))).exp()/d
    quadratic_cost = inner_bound + cross_error**2/(4*coercivity)
    gap = 2*delta*alpha*beta-quadratic_cost*alpha*alpha
    assert 0 < alpha < available and alpha < arb(1)/12
    assert delta > 0 and coercivity > 0 and gap > 0
    inner_radius_squared = -2*((arb(-1)/2).exp()-alpha).log()
    assert 1 < inner_radius_squared < arb(25)/16
    assert ((arb(-1)/2).exp()-(-inner_radius_squared/2).exp()-alpha).contains(0)
    remaining = [arb(1)/3-alpha/3-beta/2]*2+[arb(1)/3-alpha/3]
    assert all(x > 0 for x in remaining)
    assert (sum(remaining, arb(0))-(1-alpha-beta)).contains(0)
    return {'rho': '1/2', 'inner_mass': str(alpha_exact),
            'outer_annulus': ['2', '3'], 'inner_radius_squared': enclose(inner_radius_squared),
            'outer_mass': enclose(beta), 'inner_mass_available': enclose(available),
            'cross_gap_delta': enclose(delta), 'outer_phase_coercivity': enclose(coercivity),
            'cross_phase_lipschitz': enclose(lipschitz),
            'cross_error_coefficient': enclose(cross_error),
            'inner_self_kernel_bound': enclose(inner_bound),
            'quadratic_cost': enclose(quadratic_cost),
            'stability_advantage_over_every_consecutive_partition': enclose(gap),
            'remaining_pure_radial_masses': [enclose(x) for x in remaining]}


def certificate():
    # Up to a common rotation these exhaust the six labeled orders of
    # three equal consecutive arcs. Offsets from opposite half-circle
    # centers are always separated by pi/3.
    count = 0
    for centers in permutations((Q(0), Q(2, 3), Q(4, 3))):
        assert circular_distance_pi(centers[0]-(centers[1]-1)) == Q(1, 3)
        count += 1
    assert Q(1, 2)*Q(1, 6)**2 == Q(1, 72)
    result = constants(Q(1, 100_000_000))
    interior_tau = Q(1, 100_000_000_000_000)
    tau = arb(interior_tau.numerator)/interior_tau.denominator
    gap_lower = arb(result['stability_advantage_over_every_consecutive_partition'][0])
    assert gap_lower-40*tau > 0
    try:
        constants(Q(1, 100))
    except AssertionError:
        rejected = True
    else:
        raise AssertionError('the larger-annulus control must fail this sufficient margin')
    return {'status': 'pass', 'arb_precision_bits': 100,
            'arithmetic': 'python-flint Arb directed real balls; exact Fraction order controls',
            'labeled_consecutive_orders': count,
            'oversized_inner_mass_control_rejected': rejected,
            'explicit_instance': result,
            'interior_profile_control': {
                'mixing_parameter': str(interior_tau),
                'minimum_circular_mass': str(interior_tau/3),
                'remaining_gap_lower_bound': enclose(gap_lower-20*tau)},
            'scope': 'One explicit balanced two-annulus instance has a positive directed lower gap over every consecutive-arc partition with its circular masses, using the written uniform comparison DG60-DG65. The all-rho theorem and interior-profile extension are analytic, not inferred from this selected correlation. No partition beating the regular simplex is asserted.'}


if __name__ == '__main__':
    print(json.dumps(certificate(), indent=2))
