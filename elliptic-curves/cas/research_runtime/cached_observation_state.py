"""Optional exact observation-membership cache, preserving MWState records.

The copied validator is bound to the frozen original source. Only repeated
membership evaluation is memoized; every other validation and all admission
semantics are inherited unchanged. This does not cache rank assertions.
"""
from dataclasses import fields
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
from . import mw_state as original
from .arithmetic import rationals
from .mw_state import MWState

ORIGINAL_SHA256 = "e53c24a786afe1096a3bc90eff94396b1b3b6b1fd3ec52c24e51c69622b4596d"
if sha256(Path(original.__file__).read_bytes()).hexdigest() != ORIGINAL_SHA256:
    raise RuntimeError("cached validator requires the pinned original MWState source")

@lru_cache(maxsize=65536)
def _contains(model, point):
    # Both arguments are immutable exact values. False results remain false.
    return model.contains(point)

class CachedObservationMWState(MWState):
    """Same portable schema, keys, observations and proof gates as MWState."""
    def __post_init__(self):
        from mod2_reduction_independence import _is_prime
        if self.arithmetic.model != self.reductions.model:
            raise ValueError("subgroup and arithmetic curve models differ")
        if not self.reductions.independent_images:
            raise ValueError("basis has no certified independent finite images")
        if self.no_two_torsion_prime is None:
            if self.rank:
                raise ValueError("a nonempty independent basis needs a torsion witness")
        else:
            p = self.no_two_torsion_prime
            polynomial = tuple(Fraction(c) for c in self.arithmetic.model.two_division_polynomial)
            if not _is_prime(p) or any(c.denominator%p==0 for c in polynomial):
                raise ValueError("invalid no-2-torsion prime for the chosen rational presentation")
            residues = tuple(c.numerator*pow(c.denominator,-1,p)%p for c in polynomial)
            if any(sum(c*x**i for i,c in enumerate(residues))%p==0 for x in range(p)):
                raise ValueError("independence descent requires a no-rational-2-torsion witness")
        gram=tuple(tuple(None if c is None else str(Fraction(str(c))) for c in row) for row in self.height_gram)
        object.__setattr__(self,"height_gram",gram)
        object.__setattr__(self,"kummer_classes",tuple(rationals(row) for row in self.kummer_classes))
        object.__setattr__(self,"observations",tuple(self.observations))
        if len(gram)!=self.rank or any(len(row)!=self.rank for row in gram):
            raise ValueError("height Gram dimensions disagree with the basis")
        if any(gram[i][j]!=gram[j][i] for i in range(self.rank) for j in range(self.rank)):
            raise ValueError("height Gram is not symmetric")
        if self.height_kind not in ("unknown","approximate","exact"):
            raise ValueError("height data need an explicit exact/approximate/unknown label")
        if len(self.kummer_classes)!=self.rank or any(self.kummer_class(p)!=beta for p,beta in zip(self.basis,self.kummer_classes)):
            raise ValueError("Kummer classes are not bound to the model transport and points")
        if any(not _contains(self.model, row.point) for row in self.observations):
            raise ValueError("observed point is off the exact curve")
        for row in (*self.parity.generic_coordinates,*self.parity.quotient_coordinates):
            if len(row)!=self.rank:raise ValueError("stale parity-lattice presentation")

    def adjoin(self, point, *, cache, extra_primes=(), height_row=None):
        result = super().adjoin(point, cache=cache, extra_primes=extra_primes, height_row=height_row)
        # The original method returns the base class after a new admission.
        # Preserve this optional validator for the next operation as well.
        if type(result) is MWState:
            result = type(self)(**{field.name: getattr(result, field.name) for field in fields(result)})
        return result
