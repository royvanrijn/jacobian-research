"""Sage/PARI discovery adapter; no BNF or Selmer work occurs by default."""

from __future__ import annotations

import base64
from pathlib import Path
import tempfile
from dataclasses import replace

from .arithmetic import ArithmeticContext, CurveModel, TwoTorsionContext, rationals
from .store import default_store


class SageArithmetic:
    def __init__(self, store=None):
        from sage.all import pari
        from sage.version import version
        import os
        if os.environ.get("EC_PARI_STACK_BYTES"):
            limit = int(os.environ["EC_PARI_STACK_BYTES"])
            pari.allocatemem(min(64_000_000, limit), limit, silent=True)
        self.pari = pari
        self.store = store or default_store()
        self.version = f"sage-{version}/pari-{pari.version()}/arithmetic-1"
        self._live = {}

    def _blob(self, value):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "pari.bin"
            self.pari.writebin(str(path), value)
            return base64.b64encode(path.read_bytes()).decode("ascii")

    def _restore(self, blob):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "pari.bin"
            path.write_bytes(base64.b64decode(blob, validate=True))
            return self.pari.read(str(path))

    def _fact(self, kind, identity, build, discover):
        if discover:
            return self.store.discover(kind, identity, build, version=self.version)
        return self.store.require(kind, identity, version=self.version)

    def prepare(self, coefficients, *, factor_primes=(), discover=False):
        """Normalize/factor once; all subsequent workers require this snapshot.

        Existing prime hints are registered *before* minimalization. A supplied
        support must factor the input discriminant completely, so a typo cannot
        cause a hidden general factorization. With no support, explicit discovery
        may factor; it should run under the common supervisor's declared budget.
        """
        from sage.all import EllipticCurve, QQ, ZZ
        model = CurveModel(tuple(coefficients))

        def build():
            primes = tuple(sorted(set(map(int, factor_primes))))
            if not primes:
                primes=tuple(sorted({p for value in (model.discriminant.numerator,model.discriminant.denominator)
                                     for p,_ in self.factor_integer(value,discover=True)}))
            if primes:
                remaining = abs(QQ(model.discriminant.numerator) / model.discriminant.denominator)
                for p in primes:
                    if not ZZ(p).is_prime(proof=True):
                        raise ValueError("a supplied discriminant divisor is not prime")
                    remaining /= QQ(p)**remaining.valuation(p)
                if remaining != 1:
                    raise ValueError("supplied primes do not cover the model discriminant")
                self.pari.addprimes(list(primes))
            curve = EllipticCurve(QQ, list(model.coefficients))
            minimal = curve.global_minimal_model()
            transport = curve.isomorphism_to(minimal).tuple()
            delta = abs(ZZ(minimal.discriminant()))
            factors = tuple((p, int(delta.valuation(p))) for p in primes if delta % p == 0)
            minimal_model = CurveModel(tuple(minimal.a_invariants()))
            return ArithmeticContext(model, minimal_model, rationals(transport), factors,
                                     TwoTorsionContext(minimal_model.two_division_polynomial),
                                     ("0", "1", "0")).record()

        row = self._fact("arithmetic-context", {"model": model.coefficients}, build, discover)
        return ArithmeticContext.from_record(row)

    def prepare_context(self, context, *, factor_primes=(), discover=False):
        """Upgrade a raw context, retaining explicit labelled generator maps."""
        from fractions import Fraction
        if context.prepared:return context
        generator=tuple(Fraction(c)/4 for c in context.curve_generator_in_algebra)
        return self.prepare_congruent(context.model.coefficients,context.two_torsion,generator,
                                      factor_primes=factor_primes,discover=discover)

    def prepare_congruent(self,coefficients,algebra,generator_in_raw_model,*,factor_primes=(),discover=False):
        """Attach an explicitly 2-congruent curve to a shared labelled algebra.

        generator_in_raw_model represents the raw x-coordinate of a 2-torsion
        point. Minimalization and z=4*x transport it; ArithmeticContext checks
        the resulting cubic identity and invertibility exactly.
        """
        from fractions import Fraction
        from math import lcm
        context=self.prepare(coefficients,factor_primes=factor_primes,discover=discover)
        u,r,_,_=map(Fraction,context.minimal_to_input)
        alpha=list(map(lambda c:Fraction(str(c)),generator_in_raw_model))
        if len(alpha)!=3:raise ValueError("three raw generator coordinates required")
        # PARI returns [nf,map] for a nonintegral defining polynomial. Normalize
        # the generator explicitly, before any order/local data are requested;
        # treating that pair as an nf would corrupt every power-basis coordinate.
        scale=lcm(*(Fraction(c).denominator for c in algebra.polynomial))
        if scale!=1:
            original=algebra
            algebra=TwoTorsionContext(tuple(Fraction(c)*scale**(3-i) for i,c in enumerate(original.polynomial)),
                (*original.labels,f'integral-generator={scale}*labelled-generator',f'source-algebra={original.key}'))
            alpha=[c/scale**i for i,c in enumerate(alpha)]
        alpha[0]-=r
        return replace(context,two_torsion=algebra,curve_generator_in_algebra=rationals(4*c/u**2 for c in alpha))

    def factor_integer(self,value,*,discover=False):
        """Factor cover-specific integers once, separately from field setup."""
        from sage.all import ZZ
        value=abs(ZZ(value))
        if value==0:raise ValueError("zero has no finite prime factorization")
        def build():return {"factors":[[int(p),int(e)] for p,e in value.factor(proof=True)]}
        return self._fact("integer-factorization",str(value),build,discover)["factors"]

    def square_root(self,context,coefficients,*,discover=False):
        """Retain an exact square-root witness in the already prepared field."""
        from sage.all import QQ
        coefficients=rationals(coefficients)
        if len(coefficients)!=3:raise ValueError("three cubic power coordinates required")
        identity={"algebra":context.key,"value":coefficients}
        nf=self.nf(context);theta=self.pari.Mod("y",nf.nf_get_pol())
        value=sum(self.pari(QQ(c))*theta**i for i,c in enumerate(coefficients))
        def build():
            roots=self.pari.nfroots(nf,self.pari("x")**2-value)
            return {"root":None if not len(roots) else [str(self.pari.polcoef(self.pari.lift(roots[0]),i)) for i in range(3)]}
        row=self._fact("two-torsion/square-root",identity,build,discover)
        if row["root"] is None:return None
        root=sum(self.pari(QQ(c))*theta**i for i,c in enumerate(row["root"]))
        if root*root!=value:raise ArithmeticError("retained field square root is invalid")
        return tuple(row["root"])

    def field(self, context: TwoTorsionContext, *, factor_primes=(), discover=False):
        """Maximal etale cubic order and polredabs maps, retained by labels.

        The nf structures remain in the original labelled generators. Reducible
        cubics retain their rational/number-field factors separately.
        """
        from sage.all import PolynomialRing, QQ, ZZ
        if any(QQ(c).denominator()!=1 for c in context.polynomial):
            raise ValueError('order setup needs an integral labelled generator; prepare the curve context first')

        def build():
            ring = PolynomialRing(QQ, "y")
            polynomial = ring(list(context.polynomial))
            primes = tuple(sorted(set(map(int, factor_primes))))
            disc = polynomial.discriminant()
            if not primes:
                primes=tuple(sorted({p for value in (disc.numerator(),disc.denominator())
                                     for p,_ in self.factor_integer(value,discover=True)}))
            remainder = abs(disc)
            for p in primes:
                if not ZZ(p).is_prime(proof=True):
                    raise ValueError("field support contains a composite")
                remainder /= QQ(p)**remainder.valuation(p)
            if remainder != 1:
                raise ValueError("incomplete cubic discriminant support")
            self.pari.addprimes(list(primes))
            components = []
            for factor, multiplicity in polynomial.factor():
                if multiplicity != 1:
                    raise ValueError("non-etale algebra")
                if factor.degree() == 1:
                    components.append({"polynomial": str(factor), "degree": 1,
                                       "rational_root": str(-factor[0]/factor[1])})
                    continue
                nf = self.pari.nfinit([self.pari(factor), list(primes)])
                if list(self.pari.nfcertify(nf)):
                    raise ArithmeticError("maximal order is not certified")
                reduced = self.pari.polredabs(nf, 1)
                # polredabs flag 1 returns the ORIGINAL generator in the
                # reduced field; retain its inverse as well, checking composition.
                inverse = None
                for candidate in self.pari.nfisisom(reduced[0], self.pari(factor)):
                    if self.pari.subst(self.pari.lift(reduced[1]), "y", candidate) == self.pari.Mod("y", self.pari(factor)):
                        inverse = candidate
                        break
                if inverse is None:
                    raise ArithmeticError("polredabs inverse did not compose to the labelled generator")
                components.append({"polynomial": str(factor), "degree": int(factor.degree()),
                                   "nf_binary": self._blob(nf), "maximal_order_basis": str(nf.nf_get_zk()),
                                   "field_discriminant": str(nf[2]),
                                   "reduced_polynomial": str(reduced[0]),
                                   "original_generator_in_reduced": str(reduced[1]),
                                   "reduced_generator_in_original": str(inverse)})
            return {"labelled_algebra": context.key, "components": components,
                    "discriminant_support": list(primes), "maximal_order_certified": True}

        return self._fact("two-torsion/order", context.key, build, discover)

    def nf(self, context, *, factor_primes=(), discover=False, component=0):
        key = (context.key, component)
        if key not in self._live:
            row = self.field(context, factor_primes=factor_primes, discover=discover)
            item = row["components"][component]
            if item["degree"] == 1:
                raise ValueError("a rational algebra component has no nf structure")
            self._live[key] = self._restore(item["nf_binary"])
        return self._live[key]

    def reduced_field(self, context, *, factor_primes=(), discover=False):
        """Transport an already certified maximal basis to the polredabs model.

        PARI accepts [polynomial, integral_basis, ramified_primes]. This avoids
        factoring the discriminant of the new defining order, whose index can
        contain primes absent from the original discriminant support.
        """
        from sage.all import QQ
        source = self.field(context, factor_primes=factor_primes, discover=discover)
        if len(source["components"]) != 1 or source["components"][0]["degree"] != 3:
            raise ValueError("reduced-field adapter requires an irreducible cubic")
        item = source["components"][0]
        polynomial = self.pari(item["reduced_polynomial"])
        target = TwoTorsionContext(tuple(str(self.pari.polcoef(polynomial, i)) for i in range(4)))
        if target.key == context.key:
            return target
        def build():
            nf = self.nf(context)
            image = self.pari(item["original_generator_in_reduced"])
            basis = [self.pari.lift(self.pari.subst(b, "y", image)) for b in nf.nf_get_zk()]
            delta = abs(QQ(item["field_discriminant"]))
            primes = [p for p in source["discriminant_support"] if delta.valuation(p)]
            remaining = delta
            for p in primes:
                remaining /= QQ(p)**remaining.valuation(p)
            if remaining != 1:
                raise ArithmeticError("transported maximal order lacks certified ramification support")
            reduced_nf = self.pari.nfinit([polynomial, basis, primes])
            if str(reduced_nf[2]) != item["field_discriminant"] or list(self.pari.nfcertify(reduced_nf)):
                raise ArithmeticError("maximal basis transport failed")
            identity = str(self.pari.Mod("y", polynomial))
            component = {"polynomial": str(polynomial), "degree": 3,
                         "nf_binary": self._blob(reduced_nf),
                         "maximal_order_basis": str(reduced_nf.nf_get_zk()),
                         "field_discriminant": str(reduced_nf[2]),
                         "reduced_polynomial": str(polynomial),
                         "original_generator_in_reduced": identity,
                         "reduced_generator_in_original": identity}
            return {"labelled_algebra": target.key, "components": [component],
                    "discriminant_support": primes, "maximal_order_certified": True,
                    "transported_from": context.key, "source_generator": str(image),
                    "defining_polynomial_discriminant_factored": False}
        self._fact("two-torsion/order", target.key, build, discover)
        return target

    def prime_ideals(self, context, prime, *, discover=False, component=0):
        from sage.all import ZZ
        prime = int(prime)
        if not ZZ(prime).is_prime(proof=True):
            raise ValueError("prime ideals require a rational prime")
        identity = {"algebra": context.key, "component": component, "prime": prime}
        def build():
            value = self.pari.idealprimedec(self.nf(context, component=component), prime)
            return {"binary": self._blob(value), "ideals": [str(item) for item in value]}
        row = self._fact("two-torsion/prime-ideals", identity, build, discover)
        return self._restore(row["binary"])

    def local_factorization(self, context, prime, precision, *, discover=False):
        identity = {"algebra": context.key, "prime": int(prime), "precision": int(precision)}
        if int(precision) < 1:
            raise ValueError("positive local precision required")
        def build():
            from sage.all import PolynomialRing, QQ, ZZ
            if not ZZ(prime).is_prime(proof=True):
                raise ValueError("a completion requires a prime")
            pol = self.pari(PolynomialRing(QQ, "y")(list(context.polynomial)))
            factors = self.pari.factorpadic(pol, int(prime), int(precision))
            return {"binary": self._blob(factors), "factors": str(factors)}
        row = self._fact("two-torsion/local-factors", identity, build, discover)
        return self._restore(row["binary"])

    def bnf(self, context, *, requirement=None, discover=False, component=0,
            flag=1, tech=(), certify_flag=0):
        """Explicit upper-bound/completeness path; never a scheduling default."""
        if discover and requirement not in ("unconditional-upper-bound", "complete-selmer", "frozen-bnf-experiment", "two-primary-upper-bound"):
            raise ValueError("full BNF discovery requires an explicit completeness purpose")
        if flag not in (0, 1) or certify_flag not in (0, 1):
            raise ValueError("unsupported BNF flag or certification mode")
        if certify_flag == 1 and requirement not in ("two-primary-upper-bound", "frozen-bnf-experiment"):
            raise ValueError("one-sided BNF certification cannot supply a complete Selmer group")
        identity = {"algebra": context.key, "component": component}
        namespace = "two-torsion/bnf"
        if flag != 1 or tech or certify_flag != 0:
            from fractions import Fraction
            identity.update({"flag": flag, "tech": [str(Fraction(str(c))) for c in tech],
                             "certify_flag": certify_flag})
            namespace = "two-torsion/bnf-with-options"
        def build():
            value = self.pari.bnfinit(self.nf(context, component=component), flag, list(tech))
            if int(self.pari.bnfcertify(value, certify_flag)) != 1:
                raise ArithmeticError("BNF certification failed")
            return {"binary": self._blob(value), "class_group": str(value.bnf_get_cyc()),
                    "units": str(value.bnf_get_fu()), "certified": certify_flag == 0,
                    "certification": "full" if certify_flag == 0 else "one-sided-class-group-bound"}
        row = self._fact(namespace, identity, build, discover)
        return self._restore(row["binary"])

    def scheduling_features(self, context):
        """Read known 2-primary/local features; never start BNF to fill gaps."""
        field = self.field(context)
        components = []
        for index, item in enumerate(field["components"]):
            known = self.store.get("two-torsion/bnf", {"algebra": context.key, "component": index}, version=self.version)
            components.append({"degree": item["degree"], "class_group": known["class_group"] if known else None,
                               "class_group_status": "CERTIFIED" if known else "UNKNOWN"})
        return {"algebra": context.key, "components": components,
                "discriminant_support": field["discriminant_support"],
                "full_bnf_requested": False, "used_as_mathematical_exclusion": False}

    def fast_features(self, context:TwoTorsionContext, *, primes=(2,3,5,7)):
        """2-primary/local scheduling path, without order or BNF discovery.

        At a squarefree reduction, factor degrees also describe completions.
        Repeated factors are only residual-polynomial features: their true
        completion data remain UNKNOWN until a separately prepared context
        supplies them. No class-group completeness is inferred from this path.
        """
        import json
        from sage.all import GF,PolynomialRing,QQ,ZZ
        from .store import FiniteFieldFacts
        ring=PolynomialRing(QQ,"y");polynomial=ring(list(context.polynomial))
        discriminant=QQ(polynomial.discriminant())
        local=[];facts=FiniteFieldFacts(self.store)
        for prime in primes:
            prime=int(prime)
            if not ZZ(prime).is_prime(proof=True):raise ValueError("feature places must be prime")
            def build(prime=prime):
                if any(QQ(c).denominator()%prime==0 for c in context.polynomial):
                    return {"prime":prime,"status":"NONINTEGRAL_PRESENTATION","completion_factor_degrees":None}
                reduced=PolynomialRing(GF(prime),"y")(list(map(QQ,context.polynomial)))
                factors=[{"degree":int(f.degree()),"multiplicity":int(e),
                          "coefficients":[str(c) for c in f]} for f,e in reduced.factor()]
                squarefree=bool(reduced.is_squarefree())
                return {"prime":prime,"status":"UNRAMIFIED_COMPLETIONS" if squarefree else "RESIDUAL_FEATURE_ONLY",
                        "polynomial_discriminant_valuation":int(discriminant.valuation(prime)),
                        "factorization":factors,"squarefree_reduction":squarefree,
                        "completion_factor_degrees":[r["degree"] for r in factors] if squarefree else None}
            local.append(facts.query({"labelled_cubic":context.key},prime,"scheduling-local",
                                    build=build,version="2-primary-features-1"))
        order=self.store.get("two-torsion/order",context.key,version=self.version)
        classes=[]
        if order:
            for i,component in enumerate(order["components"]):
                if component["degree"]==1:
                    classes.append({"component":i,"class_group_two_rank":0,"status":"RATIONAL_COMPONENT"})
                    continue
                known=self.store.get("two-torsion/bnf",{"algebra":context.key,"component":i},version=self.version)
                classes.append({"component":i,"class_group_two_rank":
                    sum(int(n)%2==0 for n in json.loads(known["class_group"])) if known else None,
                    "status":"CERTIFIED" if known else "UNKNOWN"})
        return {"schema":"elliptic-curves.local-scheduling-features.v1","algebra":context.key,
                "real_signature":[3,0] if discriminant>0 else [1,1],
                "local_features":local,"known_class_group_two_parts":classes or None,
                "maximal_order_discovery_requested":False,"full_bnf_requested":False,
                "used_as_mathematical_exclusion":False}

    def full_selmer(self, context: ArithmeticContext, *, requirement=None, discover=False):
        """Optional complete descent using the prepared BNF and local contexts.

        This explicitly bypasses monolithic ellrankinit. Simon's global
        S-squareclass/norm computation and simultaneous local intersections
        consume the retained field and per-place ppinit structures. Curve
        discriminant support is supplied rather than factored in the backend.
        """
        context.require_prepared()
        if discover and requirement not in ("unconditional-upper-bound", "complete-selmer"):
            raise ValueError("complete Selmer discovery requires an upper-bound/completeness purpose")
        from sage.all import PolynomialRing, QQ, pari
        from sage.env import SAGE_EXTCODE
        from hashlib import sha256
        from run_elkies_2026_relative_2selmer_checkpointed import SIMON_GP_FUNCTION
        simon=Path(SAGE_EXTCODE)/"pari"/"simon"
        names=("ellQ.gp","ell.gp","qfsolve.gp","resultant3.gp")
        implementation={name:sha256((simon/name).read_bytes()).hexdigest() for name in names}
        implementation["adapted_function"]=sha256(SIMON_GP_FUNCTION.encode()).hexdigest()
        implementation["context_adapter"]="prepared-support-reduced-4"
        identity={"context":context.key,"implementation":implementation}

        def build():
            field=self.field(context.two_torsion)
            if len(field["components"])!=1 or field["components"][0]["degree"]!=3:
                raise ValueError("this complete-descent backend requires an irreducible cubic")
            working_context=self.reduced_field(context.two_torsion,discover=True)
            working_field=self.field(working_context)
            bnf=self.bnf(working_context,requirement=requirement,discover=True)
            from .simon import prepared_simon
            primes=sorted({2,*context.bad_primes,*working_field["discriminant_support"]})
            c,b,a,_=map(QQ,context.minimal_model.two_division_polynomial)
            curve=pari.ellinit([0,a,0,b,c])
            nf=self.nf(context.two_torsion)
            component=field["components"][0]
            theta=pari(component["original_generator_in_reduced"])
            inverse=pari(component["reduced_generator_in_original"])
            alpha=sum(pari(QQ(c))*theta**i for i,c in enumerate(context.curve_generator_in_algebra))
            raw=prepared_simon(self,curve,bnf,alpha,primes=primes,
                class_data_id={"algebra":working_context.key,"certification":"full","bnf_flag":1},discover=True)
            dimension=len(raw[1])
            basis=[]
            for column in raw[1]:
                beta=pari(1)
                for i,bit in enumerate(column):
                    if int(bit):beta*=raw[0][i]
                beta=pari.subst(pari.lift(beta),"y",inverse)
                coefficients=[str(pari.polcoef(pari.lift(beta),i)) for i in range(3)]
                basis.append({"beta_power_coordinates":coefficients,
                              "norm":str(pari.nfeltnorm(nf,beta))})
            return {"schema":"elliptic-curves.cached-complete-selmer.v1",
                    "arithmetic_context":context.key,"labelled_algebra":context.two_torsion.key,
                    "working_labelled_algebra":working_context.key,
                    "working_generator_in_input_algebra":str(inverse),
                    "full_selmer_dimension":dimension,"basis":basis,
                    "local_places":[int(p) for p in raw[8]],
                    "local_audit":[[int(v) for v in row] for row in raw[3]],
                    "local_metadata":[[int(v) for v in row] for row in raw[10]],
                    "raw_binary":self._blob(raw),"class_group_certified":True,
                    "status":"COMPLETE_CERTIFIED_SELMER_BASIS_NO_POINT_CLAIMS"}
        return self._fact("arithmetic/full-selmer",identity,build,discover)
