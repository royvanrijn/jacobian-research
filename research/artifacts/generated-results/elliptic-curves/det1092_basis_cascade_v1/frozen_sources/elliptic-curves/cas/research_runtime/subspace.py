"""Descent on a supplied global squareclass subspace, without ambient BNF.

Arithmetic backends provide exact witnesses and separate replay methods.
The linear layer never upgrades a partial local image or a CT radical into
a complete Selmer group or a soluble class.
"""

from dataclasses import dataclass
from typing import Protocol

from .binary import BinaryBasis, combine, kernel_masks, pack, unpack
from .store import digest


@dataclass(frozen=True)
class GlobalSquareclasses:
    algebra_key: str
    representatives: tuple[tuple[str, ...], ...]
    independence_witness: str

    def __post_init__(self):
        object.__setattr__(self, "representatives", tuple(tuple(row) for row in self.representatives))
        if not self.algebra_key or not self.independence_witness:
            raise ValueError("global squareclasses require labelled algebra and independence evidence")

    @property
    def dimension(self):
        return len(self.representatives)

    @property
    def key(self):
        return digest(self)


def local_intersection(dimension, restrictions):
    """One kernel of the map to all local quotients; no class-by-class filter.

    Each restriction has one row per global generator. Its target is the
    local squareclass space modulo the *complete* local Kummer image.
    Completeness and maps must be established by arithmetic replay upstream.
    """
    columns = [[] for _ in range(dimension)]
    for rows in restrictions:
        if len(rows) != dimension:
            raise ValueError("local map is not bound to the whole global basis")
        width = len(rows[0]) if rows else 0
        if any(len(row) != width for row in rows):
            raise ValueError("ragged local map")
        for i, row in enumerate(rows):
            pack(row)
            columns[i].extend(row)
    return kernel_masks(columns)


def restricted_radical(admissible_masks, pairing, *, global_dimension):
    """Return necessary point classes and the exact obstruction count.

    This is the radical of the *restricted* CT pairing. Classes inside it
    retain point-or-Sha status UNKNOWN. It gives no full-curve rank upper bound.
    """
    masks = tuple(admissible_masks)
    basis = BinaryBasis(global_dimension)
    for mask in masks:
        basis, dependency = basis.append(mask)
        if dependency is not None:
            raise ValueError("dependent admissible basis")
    n = len(masks)
    if len(pairing) != n or any(len(row) != n for row in pairing):
        raise ValueError("incomplete restricted pairing matrix")
    for i, row in enumerate(pairing):
        pack(row)
        if row[i] or any(row[j] != pairing[j][i] for j in range(n)):
            raise ValueError("Cassels-Tate matrix must be alternating")
    kernel = kernel_masks(pairing, width=n)
    return {"admissible_dimension": n, "pairing_rank": n-len(kernel),
            "radical_dimension": len(kernel), "radical_coordinates": list(kernel),
            "radical_global_masks": [combine(mask, masks) for mask in kernel],
            "obstructed_class_count": (1 << n) - (1 << len(kernel)),
            "nonzero_compatible_class_count": (1 << len(kernel))-1,
            "full_selmer_computed": False, "full_curve_rank_upper": None,
            "radical_point_or_sha_status": "UNKNOWN"}


class SubspaceBackend(Protocol):
    """Discovery and cheap witness replay are deliberately separate methods."""

    def verify_global(self, context, classes): ...
    def required_places(self, context, classes): ...
    def local_map(self, context, classes, place): ...
    def verify_local(self, context, classes, place, record): ...
    def cover(self, context, classes, mask): ...
    def verify_cover(self, context, classes, mask, record): ...
    def ct_pairing(self, context, classes, masks, covers): ...
    def verify_ct(self, context, classes, masks, covers, record): ...


class SubspaceDescent:
    def __init__(self, context, classes: GlobalSquareclasses, backend: SubspaceBackend):
        if context.two_torsion.key != classes.algebra_key:
            raise ValueError("squareclasses belong to a different labelled algebra")
        self.context, self.classes, self.backend = context, classes, backend

    def run(self, *, retained=None, include_ct=True):
        """Build a witness package, or replay it without any discovery calls.

        The backend verifier returns exactly True only after checking its
        mathematical witness. Missing/incomplete arithmetic fails closed.
        A callers' complete-Selmer requirement belongs in a separate backend.
        """
        context, classes, backend = self.context, self.classes, self.backend
        identity = {"arithmetic_context": context.key, "global_classes": classes.key}
        if retained is not None and retained.get("identity") != identity:
            raise ValueError("subspace witness identity mismatch")
        if backend.verify_global(context, classes) is not True:
            raise ArithmeticError("unverified global squareclass independence")
        places = tuple(backend.required_places(context, classes))
        if len(set(places)) != len(places):
            raise ValueError("duplicate required places")
        # Rational descent must at least include 2, infinity and the curve's
        # bad primes. The backend additionally includes representative support.
        if not {2, "infinity", *context.bad_primes} <= set(places):
            raise ValueError("incomplete mandatory local-place support")
        if retained is not None and retained.get("places") != list(places):
            raise ValueError("retained local-place scope mismatch")
        local_records = []
        for index, place in enumerate(places):
            row = retained["local_maps"][index] if retained is not None else backend.local_map(context, classes, place)
            if backend.verify_local(context, classes, place, row) is not True:
                raise ArithmeticError(f"incomplete/unverified local quotient at {place}")
            local_records.append(row)
        masks = local_intersection(classes.dimension, [row["quotient_rows"] for row in local_records])
        result = {"schema": "elliptic-curves.subspace-descent.v1", "identity": identity,
                  "places": list(places), "local_maps": local_records,
                  "admissible_masks": list(masks), "full_selmer_computed": False,
                  "full_curve_rank_upper": None, "status": "LOCAL_INTERSECTION_ONLY"}
        if retained is not None and retained.get("admissible_masks") != list(masks):
            raise ArithmeticError("tampered admissible subspace")
        if not include_ct:
            return result
        covers = []
        for index, mask in enumerate(masks):
            cover = retained["covers"][index] if retained is not None else backend.cover(context, classes, mask)
            if backend.verify_cover(context, classes, mask, cover) is not True:
                raise ArithmeticError("unverified explicit covering map")
            covers.append(cover)
        ct = retained["ct"] if retained is not None else backend.ct_pairing(context, classes, masks, covers)
        if backend.verify_ct(context, classes, masks, covers, ct) is not True:
            raise ArithmeticError("incomplete/unverified restricted CT pairing")
        radical = restricted_radical(masks, ct["matrix"], global_dimension=classes.dimension)
        if retained is not None and retained.get("radical") != radical:
            raise ArithmeticError("tampered radical witness")
        result.update({"covers": covers, "ct": ct, "radical": radical,
                       "status": "VERIFIED_RESTRICTED_CT_RADICAL"})
        return result

    def search(self, state, witness, masks, *, search, limits, registry=None, pruning_audit=None):
        """Search explicitly requested radical classes, updating MWState.

        Masks are supplied lazily by the common search policy. No exponential
        subset list is constructed here. A bounded miss supplies no theorem.
        """
        from .supervisor import Limits
        from .mw_state import MWState
        from .pruning import PruningRegistry, SearchRequest
        from pathlib import Path
        if not isinstance(limits, Limits):
            raise ValueError("bounded subspace search requires declared worker limits")
        if not isinstance(state,MWState) or state.arithmetic.key != self.context.key:
            raise ValueError("subspace search requires an MWState on this exact arithmetic context")
        replay = self.run(retained=witness)
        radical = BinaryBasis(self.classes.dimension)
        registry=registry or PruningRegistry(Path(__file__).resolve().parents[3])
        scope=(('curve',state.model.key),('squareclass_space',digest(
            {'algebra':self.classes.algebra_key,'basis':self.classes.representatives})))
        for mask in replay["radical"]["radical_global_masks"]:
            radical, _ = radical.append(mask)
        for mask in masks:
            if not mask or radical.reduce(mask)[0]:
                raise ValueError("requested class is outside the restricted radical")
            decision=registry.decision(SearchRequest('kummer_class',scope,class_mask=mask,
                                                     class_dimension=self.classes.dimension))
            if pruning_audit is not None:pruning_audit.append({'mask':mask,**decision})
            if not decision['search_allowed']:continue
            state = search(state, self.context, self.classes, mask, limits)
            if not isinstance(state,MWState) or state.arithmetic.key != self.context.key:
                raise ValueError("point search returned a misbound subgroup state")
        return state
