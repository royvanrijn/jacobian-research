"""Exact finite controls for the reflection comparisons in DG12--DG13.

The Gaussian orbit identity and the balanced small-ball realization use the
written analytic proof. This program exhausts its finite label identities,
orientation choices and rational robustness margin; it is not a Gaussian
partition search or a proof of the all-correlation simplex inequality.
"""
from fractions import Fraction as Q
from itertools import combinations, product
import json


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), Q(0))


def finite_controls():
    unit = [tuple(Q(i == j) for j in range(3)) for i in range(3)]
    count = 0
    for a, b, c, d in product(range(3), repeat=4):
        m = tuple((x+y)/2 for x, y in zip(unit[a], unit[b]))
        n = tuple((x+y)/2 for x, y in zip(unit[c], unit[d]))
        delta = tuple(x-y for x, y in zip(unit[a], unit[b]))
        epsilon = tuple(x-y for x, y in zip(unit[c], unit[d]))
        same = dot(unit[a], unit[c]) + dot(unit[b], unit[d])
        reflected = dot(unit[a], unit[d]) + dot(unit[b], unit[c])
        assert same == 2*dot(m, n) + dot(delta, epsilon)/2
        assert reflected == 2*dot(m, n) - dot(delta, epsilon)/2
        assert abs(dot(delta, epsilon)) == dot(tuple(abs(x) for x in delta),
                                              tuple(abs(x) for x in epsilon))
        count += 1

    incidence = ((1, -1, 0), (0, 1, -1), (-1, 0, 1))
    for size in (1, 2):
        for types in combinations(range(3), size):
            gauges = list(product((-1, 1), repeat=size))
            assert any(all(g[i]*g[j]*dot(incidence[types[i]], incidence[types[j]]) > 0
                           for i in range(size) for j in range(size)) for g in gauges)
    assert not any(all(g[i]*g[j]*dot(incidence[i], incidence[j]) > 0
                       for i in range(3) for j in range(3))
                   for g in product((-1, 1), repeat=3))

    # Every cyclic three-type triangle has a frustrated edge. Integration
    # of this pointwise inequality, not a sample of masses, proves DG56.
    triangle_values = [sum(1+s[i]*s[j] for i in range(3) for j in range(i+1, 3))
                       for s in product((-1, 1), repeat=3)]
    assert min(triangle_values) == 2 and max(triangle_values) == 6

    types = (0, 1, 2, 0)
    weights = {step: Q(1, 2)**(step*step) for step in (1, 2, 3)}

    def energy(signs):
        return sum((weights[j-i]*dot(incidence[types[i]], incidence[types[j]])
                    *signs[i]*signs[j] for i in range(4) for j in range(i+1, 4)), Q(0))

    choices = [(energy(s), s) for s in product((-1, 1), repeat=4)]
    best = max(value for value, signs in choices)
    constant = max(value for value, signs in choices if signs[0] == signs[3])
    split = (1, -1, 1, -1)
    assert energy(split) == best == Q(351, 256)
    assert constant == Q(161, 256)
    assert best-constant == Q(95, 128)
    absolute_mass = sum((weights[j-i]*abs(dot(incidence[types[i]], incidence[types[j]]))
                         for i in range(4) for j in range(i+1, 4)), Q(0))
    assert absolute_mass == Q(417, 256)
    relative_error = Q(1, 50)
    robust_gap = best-constant-2*relative_error*absolute_mass
    assert robust_gap == Q(4333, 6400) > 0
    # e^(1/100)<=1/(1-1/100)<1+1/50; and e^(-1/100)>=1-1/100.
    assert Q(100, 99) < 1+relative_error
    assert 1-Q(1, 100) > 1-relative_error

    # For an unordered interaction, a flip of exactly one endpoint gives
    # minus twice that interaction, as required by the full ordered cut.
    for signs in product((-1, 1), repeat=4):
        for bits in product((0, 1), repeat=4):
            flipped = tuple((-1 if bit else 1)*sign for sign, bit in zip(signs, bits))
            cut = sum((weights[j-i]*dot(incidence[types[i]], incidence[types[j]])
                       *signs[i]*signs[j] for i in range(4) for j in range(i+1, 4)
                       if bits[i] != bits[j]), Q(0))
            assert energy(flipped)-energy(signs) == -2*cut

    return {'status': 'pass', 'arithmetic': 'exact Python Fraction arithmetic',
            'four_label_identity_cases': count,
            'orientation_patterns': len(choices),
            'type_constant_patterns': sum(s[0] == s[3] for v, s in choices),
            'flip_identity_cases': len(choices)*16,
            'two_type_positive_gauges_verified': True,
            'three_type_positive_gauge_excluded': True,
            'scalar_polarization_identity_cases': count,
            'triangle_sign_cases': len(triangle_values),
            'triangle_defect_range': [min(triangle_values), max(triangle_values)],
            'types': ['12', '23', '31', '12'],
            'normalized_neighbor_weights': {str(k): str(v) for k, v in weights.items()},
            'best_split_orientation': list(split), 'split_energy': str(best),
            'best_type_constant_energy': str(constant), 'exact_gap': str(best-constant),
            'absolute_interaction_mass': str(absolute_mass),
            'relative_kernel_error': str(relative_error), 'robust_gap': str(robust_gap),
            'scope': 'Finite identities, all16 four-patch orientations and all8 three-type triangle signs are exhausted. Written DG12 supplies the Gaussian small-ball realization; DG13 supplies the integrated frustration bound, measurable reflection-order rigidity and compactness argument. This is not a finite proof of Gaussian simplex optimality.'}


if __name__ == '__main__':
    print(json.dumps(finite_controls(), indent=2))
