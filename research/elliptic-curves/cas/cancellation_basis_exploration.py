"""Model-independent box exposure for the basis amplification comparison.

A zero fitted radius mass never removes a job. Only an exactly equal completed
box (same actual anchor, projective signed permutation, and height) is skipped.
The small boundary probe proves coordinate exposure when it succeeds; failure
to find such a witness is UNKNOWN and does not suppress a point search.
"""
from fractions import Fraction as F

from pointed_box_equivalence import box_key
from search_observability import primitive


def inverse(matrix):
    a, b, c, d = map(F, matrix)
    if a*d == b*c:
        raise ArithmeticError('singular chart')
    return d, -b, -c, a  # projective inverse


def apply(matrix, pair):
    a, b, c, d = map(F, matrix)
    m, n = pair
    return primitive(a*m+b*n, c*m+d*n)


def height(pair):
    return max(map(abs, pair))


def key(anchor, matrix, bound):
    return tuple(map(F, anchor)), box_key(matrix), bound


def in_box(raw_coordinate, matrix, bound):
    return height(apply(inverse(matrix), raw_coordinate)) <= bound


def tail_witness(mapping, base, bound, completed):
    """Exact finite witness beyond the old fitted H125000 range and prior boxes.

    completed contains only boxes on this same actual anchor. This probes at
    most 256 boundary addresses, with no square test and no point search.
    """
    if bound <= 125000 and box_key(mapping['matrix']) == box_key(base['matrix']):
        return None
    for offset in range(1, 33):
        for pair in ((bound, offset), (offset, bound), (bound, bound-offset),
                     (bound-offset, bound)):
            for sign in (1, -1):
                new = primitive(sign*pair[0], pair[1])
                raw = apply(mapping['matrix'], new)
                old = apply(inverse(base['matrix']), raw)
                if height(old) <= 125000:
                    continue
                if any(in_box(raw, row['matrix'], row['height']) for row in completed):
                    continue
                return {'new_coordinate': list(new), 'old_coordinate': list(old),
                        'new_height': height(new), 'old_height': height(old),
                        'outside_completed_boxes': True,
                        'boundary': 'Coordinate witness only; no quartic square or elliptic point is asserted.'}
    return None


class Exposure:
    """Persistent exact completed boxes; epoch-local statistical state is absent."""
    def __init__(self):
        self.completed = {}

    def seen(self, anchor, mapping, bound):
        return key(anchor, mapping['matrix'], bound) in self.completed

    def inspect(self, anchor, mapping, base, bound):
        previous = [r for k, r in self.completed.items() if k[0] == tuple(map(F, anchor))]
        return tail_witness(mapping, base, bound, previous)

    def observe(self, anchor, mapping, bound, complete):
        if complete:
            self.completed[key(anchor, mapping['matrix'], bound)] = {
                'matrix': mapping['matrix'], 'height': bound}


def prepare(model, basis, word, mapper):
    from finite_cancellation_features import neighbour_maps
    base = mapper.mapping(model, basis, {'representative': word})
    neighbours, attempts = neighbour_maps(base, mapper.pari)
    return {'models': [{'name': 'factor_free', 'mapping': base}] + neighbours,
            'attempts': attempts,
            'policy': 'Factor-free H125000, then each distinct retained prime neighbour H125000; no fitted scores.'}
