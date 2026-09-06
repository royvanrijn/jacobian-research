"""Exact span-preserving compression of the involution P -> C-P on a short curve.

This optional helper does not change existing search/admission workers. A caller
must retain the complete raw point cloud and the returned skip witnesses.
The conclusion is equality of generated integral subgroups after adjoining the
supplied subgroup, not finite-reduction dependence or an upper rank bound.
"""
from fractions import Fraction as F
from alternate_quartic_covers import point_on_short_curve, short_add
from half_lattice_pointed_sieve import linear_combination


def compress(model, subgroup, centre_coefficients, points):
    model = tuple(F(str(x)) for x in model)
    subgroup = tuple(tuple(F(str(x)) for x in p) for p in subgroup)
    points = tuple(tuple(F(str(x)) for x in p) for p in points)
    word = tuple(F(str(c)) for c in centre_coefficients)
    if len(model) != 5 or any(model[:3]):
        raise ValueError('requires a short Weierstrass model')
    if len(word) != len(subgroup) or any(c.denominator != 1 for c in word):
        raise ValueError('centre must be an integral word in the supplied subgroup')
    if any(not point_on_short_curve(model, p) for p in (*subgroup, *points)):
        raise ValueError('input point is off the exact curve')
    centre = linear_combination(model, subgroup, tuple(map(int, word)))
    if centre is None:
        raise ValueError('pointed centre is at infinity')
    kept, skipped, seen = [], [], {}
    for i, point in enumerate(points):
        partner = short_add(model, centre, (point[0], -point[1]))
        if partner in seen:
            j = seen[partner]
            if short_add(model, points[j], point) != centre:
                raise ArithmeticError('involution relation failed')
            skipped.append({'index': i, 'partner_index': j,
                            'relation': 'P[index] + P[partner_index] = centre'})
        else:
            kept.append(i)
            seen[point] = i
    return {'centre': list(map(str, centre)), 'centre_coefficients': list(map(int, word)),
            'kept_indices': kept, 'skipped': skipped,
            'claim_boundary': 'Together with the supplied subgroup, kept points generate the same integral subgroup as all input points. Raw points remain necessary for complete-cloud replay; no exact rank or unrecorded-point claim.'}
