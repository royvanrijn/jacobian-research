"""A scale-aware numerical consistency guard, never an optimality certificate."""
from math import isfinite, ulp


def checked_distance_error(norm, reported_distance, degree):
    if type(norm) is not int or norm < 0 or type(degree) is not int or degree < 1:
        raise ValueError('nonnegative exact integer norm and positive integer degree required')
    distance = float(reported_distance)
    expected = float(norm)
    scaled = degree*degree*distance
    if not isfinite(distance) or distance < 0 or not isfinite(scaled) or not isfinite(expected):
        raise ArithmeticError('nonfinite or negative numerical distance')
    error = abs(scaled-expected)
    tolerance = max(1e-6, 64*ulp(expected), 64*ulp(scaled))
    if error > tolerance:
        raise ArithmeticError(f'inconsistent CVP norm={norm}, error={error}, tolerance={tolerance}')
    return error
