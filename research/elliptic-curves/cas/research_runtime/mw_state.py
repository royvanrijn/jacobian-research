"""Immutable subgroup state with exact point admission and cache invalidation."""

from dataclasses import dataclass, replace
from fractions import Fraction

from .arithmetic import ArithmeticContext, rationals
from .finite_reduction import IncrementalReductions
from .store import canonical, digest


@dataclass(frozen=True)
class ParityLattice:
    generic_coordinates: tuple[tuple[str, ...], ...] = ()
    quotient_coordinates: tuple[tuple[str, ...], ...] = ()
    metric_id: str | None = None
    cvp_checkpoint: str | None = None
    seen_holes: tuple[int, ...] = ()

    def __post_init__(self):
        object.__setattr__(self,"generic_coordinates",tuple(rationals(row) for row in self.generic_coordinates))
        object.__setattr__(self,"quotient_coordinates",tuple(rationals(row) for row in self.quotient_coordinates))
        object.__setattr__(self,"seen_holes",tuple(self.seen_holes))


@dataclass(frozen=True)
class PointObservation:
    point: tuple[str, ...]
    status: str
    finite_relation_mask: int | None


@dataclass(frozen=True)
class MWState:
    arithmetic: ArithmeticContext
    reductions: IncrementalReductions
    no_two_torsion_prime: int | None
    height_gram: tuple[tuple[str | None, ...], ...]
    height_kind: str
    kummer_classes: tuple[tuple[str, ...], ...]
    parity: ParityLattice = ParityLattice()
    observations: tuple[PointObservation, ...] = ()
    parent_state: str | None = None

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
        if any(not self.model.contains(row.point) for row in self.observations):
            raise ValueError("observed point is off the exact curve")
        for row in (*self.parity.generic_coordinates,*self.parity.quotient_coordinates):
            if len(row)!=self.rank:raise ValueError("stale parity-lattice presentation")

    @property
    def model(self):return self.arithmetic.model

    @property
    def basis(self):return self.reductions.points

    @property
    def rank(self):return len(self.basis)

    @property
    def key(self):return digest(self)

    def kummer_class(self,point):
        x=Fraction(point[0]);u,r,_,_=map(Fraction,self.arithmetic.kummer_to_input)
        z=4*(x-r)/u**2
        alpha=tuple(map(Fraction,self.arithmetic.curve_generator_in_algebra))
        return rationals((z-alpha[0],-alpha[1],-alpha[2]))

    @classmethod
    def empty(cls,arithmetic,*,cache,primes,no_two_torsion_prime):
        reductions=IncrementalReductions.empty(arithmetic.model.coefficients,primes,cache)
        return cls(arithmetic,reductions,no_two_torsion_prime,(),"unknown",())

    def adjoin(self,point,*,cache,extra_primes=(),height_row=None):
        """Return a new state; only verified new directions enlarge its basis.

        Ambiguity is recorded as an observation and never called dependence.
        Extra primes are requested lazily only if the new column is ambiguous.
        Known Gram entries survive; missing new heights remain explicitly null.
        A basis change invalidates every parity/CVP presentation and order.
        """
        point=rationals(point)
        if not self.model.contains(point):raise ValueError("point is off the exact curve")
        a1,_,a3,_,_=map(Fraction,self.model.coefficients)
        negative=rationals((Fraction(point[0]),-Fraction(point[1])-a1*Fraction(point[0])-a3))
        if point in self.basis or negative in self.basis:
            return replace(self,observations=(*self.observations,PointObservation(point,"KNOWN_POINT_UP_TO_SIGN",None)),parent_state=self.key)
        if self.no_two_torsion_prime is None:
            return replace(self,observations=(*self.observations,PointObservation(point,"TORSION_HYPOTHESIS_UNKNOWN",None)),parent_state=self.key)
        trial,relation=self.reductions.append(point,cache)
        if not trial.independent_images:
            for prime in extra_primes:
                if prime in trial.primes:continue
                try:trial=trial.escalate(prime,cache)
                except ValueError:continue  # bad reduction or point denominator
                if trial.independent_images:break
        if not trial.independent_images:
            observation=PointObservation(point,"AMBIGUOUS_FINITE_REDUCTIONS",relation)
            return replace(self,observations=(*self.observations,observation),parent_state=self.key)
        if height_row is None:
            row=(None,)*(self.rank+1)
        else:
            row=tuple(None if c is None else str(Fraction(str(c))) for c in height_row)
            if len(row)!=self.rank+1:raise ValueError("new height row has wrong dimension")
        gram=tuple((*old,row[i]) for i,old in enumerate(self.height_gram))+ (row,)
        observation=PointObservation(point,"CERTIFIED_INDEPENDENT_POINT",None)
        return MWState(self.arithmetic,trial,self.no_two_torsion_prime,gram,self.height_kind,
            (*self.kummer_classes,self.kummer_class(point)),ParityLattice(),
            (*self.observations,observation),self.key)

    def with_arithmetic(self, context):
        """Attach newly prepared arithmetic without reconstructing the subgroup."""
        if context.model != self.model:
            raise ValueError("arithmetic upgrade changed the exact curve")
        u,r,_,_=map(Fraction,context.kummer_to_input)
        alpha=tuple(map(Fraction,context.curve_generator_in_algebra))
        classes=tuple(rationals((4*(Fraction(p[0])-r)/u**2-alpha[0],-alpha[1],-alpha[2])) for p in self.basis)
        return replace(self,arithmetic=context,kummer_classes=classes,parity=ParityLattice(),parent_state=self.key)

    def with_geometry(self,gram,*,height_kind,parity=ParityLattice()):
        return replace(self,height_gram=tuple(map(tuple,gram)),height_kind=height_kind,
                       parity=parity,parent_state=self.key)

    def verify(self, cache):
        """Replay a retained basis using the cached finite quotient witnesses.

        This is a trust-boundary check, not the point-admission operation.
        Cached signatures are O(1) per old point/place on subsequent calls.
        Heights remain labelled data, not a height certificate.
        """
        actual = IncrementalReductions.empty(self.model.coefficients, self.reductions.primes, cache)
        for point in self.basis:
            actual, _ = actual.append(point, cache)
        if actual != self.reductions:
            raise ArithmeticError("retained MWState has incorrect reduction signatures")
        return True

    def record(self):
        return {"schema": "elliptic-curves.mw-state.v1", "key": self.key,
                "state": canonical(self)}

    @classmethod
    def from_record(cls, record, *, cache):
        from .arithmetic import CurveModel, TwoTorsionContext
        from .binary import BinaryBasis
        if record.get("schema") != "elliptic-curves.mw-state.v1":
            raise ValueError("unknown MWState schema")
        row = record["state"]
        context = row["arithmetic"]
        context = ArithmeticContext(CurveModel(**context["model"]), CurveModel(**context["minimal_model"]) if context["minimal_model"] else None,
            context["minimal_to_input"], context["discriminant_factorization"],
            TwoTorsionContext(**context["two_torsion"]), context["curve_generator_in_algebra"])
        reductions = dict(row["reductions"])
        basis = reductions["basis"]
        reductions["basis"] = BinaryBasis(basis["width"], tuple(map(tuple, basis["pivots"])), basis["column_count"])
        reductions["model"] = CurveModel(**reductions["model"])
        state = cls(context, IncrementalReductions(**reductions), row["no_two_torsion_prime"],
            row["height_gram"], row["height_kind"], row["kummer_classes"], ParityLattice(**row["parity"]),
            tuple(PointObservation(tuple(observation["point"]), observation["status"], observation["finite_relation_mask"])
                  for observation in row["observations"]), row["parent_state"])
        if state.key != record["key"]:
            raise ValueError("MWState content identity mismatch")
        state.verify(cache)
        return state
