"""Boundary Weil-class vectors derived from finite decomposition orbit data."""
from math import gcd
from functools import reduce

def ramified_boundary_profile(ef):
    """One coordinate per omitted prime, not per geometric sheet/cycle."""
    ramified=[(e,f) for e,f in ef if e>1]
    return {'mandatory_primes':len(ramified),
      'canonical_coefficients':sorted(e-1 for e,f in ramified),
      'residue_degrees':sorted(f for e,f in ramified)}

def canonical_divisibility(coefficients):
    nonzero=[abs(c) for c in coefficients if c]
    return reduce(gcd,nonzero) if nonzero else 0
