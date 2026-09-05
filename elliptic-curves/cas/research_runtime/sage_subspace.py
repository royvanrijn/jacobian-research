"""Subspace descent over a labelled irreducible cubic, with exact replay.

Known classes may be certified by an MWState on any explicitly 2-congruent
curve, or an external exact global-squareclass certificate adapter. No BNF
is called. Run discovery under research_runtime.supervisor: rational conic
solving and witness finding have bounded-search failure, never rank meaning.
"""

from hashlib import sha256
from itertools import combinations

from sage.all import AA, GF, PolynomialRing, QQ, ZZ, matrix, pari, prod, vector

from .binary import BinaryBasis, pack, quotient_rows
from .fisher import fisher_gamma, hilbert_symbol, invariants, local_square, require
from .local_kummer import LocalSquareclasses
from .store import digest
from .subspace import GlobalSquareclasses


def multiply(left, right, polynomial):
    """Power-coordinate product, also over a polynomial coefficient ring."""
    coefficients = [sum(left[i]*right[k-i] for i in range(3) if 0 <= k-i < 3)
                    for k in range(5)]
    for k in (4, 3):
        for j in range(3):
            coefficients[k-3+j] -= coefficients[k]*QQ(polynomial[j])
    return tuple(coefficients[:3])


def local_candidates(cubic, p, cap):
    """Deterministic bounded witnesses; reaching the known image rank proves completeness."""
    seen = set()
    def candidates():
        for i in range(4097 if p <= 127 else 257):
            yield QQ(i)
            yield QQ(-i)
        if p != 2:
            ring = PolynomialRing(GF(p), "z")
            reduced = ring(cubic)
            for root in reduced.gcd(reduced.derivative()).roots(multiplicities=False):
                for offset in range(-64, 65):
                    yield QQ(ZZ(root)+offset*p)
            for i in range(512):
                yield QQ(int.from_bytes(sha256(f"{p}:{i}".encode()).digest(), "big") % p)
        for exponent in range(1, 9 if p == 2 else 4):
            for numerator in range(1, 1025 if p == 2 else 129):
                if numerator % p:
                    yield QQ(numerator)/ZZ(p)**exponent
                    yield -QQ(numerator)/ZZ(p)**exponent
    for value in candidates():
        if value in seen:
            continue
        if len(seen) >= cap:
            return
        seen.add(value)
        yield value


def quartic_local_witness(q, gamma, place, *, node_cap=10000):
    """Find a rational point in a local open set; only a witness is retained."""
    from collections import deque
    x = q.parent().gen()
    def valid(xx):
        return bool(q(xx) and gamma(xx) and
                    (q(xx) > 0 if place == "infinity" else local_square(q(xx), place)))
    if place == "infinity":
        for i in range(min(node_cap, 1001)):
            for xx in (QQ(i), QQ(-i)):
                if valid(xx):
                    return xx
        roots = q.roots(AA, multiplicities=False)
        for lo, hi in zip(roots, roots[1:]):
            for fraction in (QQ(1)/2, QQ(1)/3, QQ(2)/3):
                xx = QQ(((1-fraction)*lo+fraction*hi).n(256))
                if valid(xx):
                    return xx
        raise RuntimeError("real witness cap reached; local solubility remains UNKNOWN")
    p = ZZ(place)
    scaled = q.denominator()**2*q
    queue = deque([(scaled, ZZ(0), ZZ(1), False),
                   (q.parent()([scaled[4-i] for i in range(5)]), ZZ(0), ZZ(1), True)])
    for _ in range(node_cap):
        if not queue:
            break
        g, shift, scale, reciprocal = queue.popleft()
        valuation = min(c.valuation(p) for c in g if c)
        gn = g/p**(2*(valuation//2))
        for t in range(2 if p == 2 else min(int(p), 64)):
            argument = shift+scale*t
            if reciprocal and not argument:
                for exponent in range(1, 21):
                    xx = 1/QQ(shift+scale*p**exponent)
                    if valid(xx):
                        return xx
            else:
                xx = 1/QQ(argument) if reciprocal else QQ(argument)
                if valid(xx):
                    return xx
        valuation = min(c.valuation(p) for c in gn if c)
        roots = ([0, 1] if p == 2 else
                 PolynomialRing(GF(p), "z")(gn/p**valuation).roots(multiplicities=False))
        for t in roots:
            t = ZZ(t)
            next_g = gn(p*x+t)
            if next_g[0] and all(not c or c.valuation(p) > next_g[0].valuation(p)+(2 if p == 2 else 0)
                                for c in next_g.list()[1:]):
                if not local_square(next_g[0], p):
                    continue
            queue.append((next_g, shift+scale*t, scale*p, reciprocal))
    raise RuntimeError("local witness cap reached; solubility remains UNKNOWN")


class PointClassWitness:
    """Independent point Kummer classes are unramified away from origin bad primes and 2."""
    def __init__(self, state, cache):
        self.state, self.cache = state, cache

    @property
    def classes(self):
        return GlobalSquareclasses(self.state.arithmetic.two_torsion.key,
                                   self.state.kummer_classes, self.state.key)

    def verify(self, context, classes):
        state = self.state
        require(classes == self.classes and context.two_torsion == state.arithmetic.two_torsion,
                "misbound point-squareclass certificate")
        require(state.verify(self.cache), "invalid source MWState")
        return True

    def support(self):
        return (2, *self.state.arithmetic.bad_primes)


class SageSubspaceBackend:
    """Concrete local maps, explicit quartics and Fisher CT on any known subspace.

    global_witness.verify must prove independence, square norms, and that the
    supplied support contains every ramified prime of these representatives.
    PointClassWitness supplies those facts without computing global units.
    Other kinds of squareclass certificates can implement the same two methods.
    """
    def __init__(self, arithmetic, context, global_witness, *, local_candidate_cap=40000):
        context.require_prepared()
        self.arithmetic, self.context, self.global_witness = arithmetic, context, global_witness
        self.cap = local_candidate_cap
        if type(self.cap) is not int or self.cap <= 0:
            raise ValueError("positive local witness cap required")
        self.R = PolynomialRing(QQ, "x")
        self.polynomial = tuple(map(QQ, context.two_torsion.polynomial))
        require(self.R(self.polynomial).is_irreducible(), "this backend requires a cubic field")
        self.K = self.R.quotient(self.R(self.polynomial), "theta")
        self.theta = self.K.gen()
        self.alpha = self.K(list(map(QQ, context.curve_generator_in_algebra)))
        self.cubic = self.R(list(map(QQ, context.minimal_model.two_division_polynomial)))
        self.a2, self.a4, self.a6 = self.cubic[2], self.cubic[1], self.cubic[0]
        self.I, self.J = self.a2**2-3*self.a4, -2*self.a2**3+9*self.a2*self.a4-27*self.a6
        self.nf = arithmetic.nf(context.two_torsion)
        self.pari_theta = pari.Mod("y", self.nf.nf_get_pol())
        self.locals = {}

    def element(self, row):
        return sum(pari(QQ(c))*self.pari_theta**i for i, c in enumerate(row))

    def field_element(self, row):
        return self.K(list(map(QQ, row)))

    def beta(self, classes, mask):
        require(type(mask) is int and 0 < mask < 1 << classes.dimension, "invalid cover mask")
        return prod(self.field_element(row) for i, row in enumerate(classes.representatives) if mask >> i & 1)

    def coordinates(self, value):
        return [str(value.lift()[i]) for i in range(3)]

    def verify_global(self, context, classes):
        require(context == self.context and self.global_witness.verify(context, classes) is True,
                "unverified global squareclasses")
        for row in classes.representatives:
            require(len(row) == 3, "incorrect cubic coordinate width")
            norm = QQ(pari.nfeltnorm(self.nf, self.element(row)))
            require(norm != 0 and norm.is_square(), "class has nonsquare norm")
        return True

    def required_places(self, context, classes):
        support = sorted({2, *context.bad_primes, *self.global_witness.support()})
        require(all(ZZ(p).is_prime(proof=True) for p in support), "nonprime local support")
        return (*support, "infinity")

    def local(self, p):
        if p not in self.locals:
            self.locals[p] = LocalSquareclasses(self.nf, p, arithmetic=self.arithmetic,
                                               context=self.context.two_torsion)
        return self.locals[p]

    def _real_map(self, classes):
        roots = self.R(self.polynomial).roots(AA, multiplicities=False)
        embeddings = [sum(QQ(c)*root**i for i, c in enumerate(self.context.curve_generator_in_algebra))
                      for root in roots]
        image = []
        if len(roots) == 3:
            smallest = min(range(3), key=lambda i: embeddings[i])
            image = [[int(i != smallest) for i in range(3)]]
        rows = [[int(sum(QQ(c)*root**i for i, c in enumerate(row)) < 0) for root in roots]
                for row in classes.representatives]
        quotient, dimension = quotient_rows(image, rows, width=len(roots))
        return {"place": "infinity", "image_dimension": len(image), "image_rows": image,
                "quotient_dimension": dimension, "quotient_rows": quotient}

    def _finite_map(self, classes, p, xs):
        local = self.local(p)
        alpha = self.element(self.context.curve_generator_in_algebra)
        image = []
        for text in xs:
            xx = QQ(text)
            require(self.cubic(xx) != 0 and local_square(self.cubic(xx), p), "invalid local Kummer point")
            image.append(local.signature(pari(xx)-alpha))
        require(len(image) == local.point_kummer_dimension, "incomplete local point image")
        rows = [local.signature(self.element(row)) for row in classes.representatives]
        width = len(local.signature(pari(1)))
        quotient, dimension = quotient_rows(image, rows, width=width)
        return {"place": p, "basis_x": list(map(str, xs)), "image_dimension": len(image),
                "image_rows": list(map(list, image)), "quotient_dimension": dimension,
                "quotient_rows": quotient}

    def local_map(self, context, classes, place):
        if place == "infinity":
            return self._real_map(classes)
        local = self.local(place)
        alpha = self.element(context.curve_generator_in_algebra)
        def discover():
            basis = BinaryBasis(len(local.signature(pari(1))))
            xs = []
            if basis.rank < local.point_kummer_dimension:
                for xx in local_candidates(self.cubic, place, self.cap):
                    if not self.cubic(xx) or not local_square(self.cubic(xx), place):
                        continue
                    signature = pack(local.signature(pari(xx)-alpha))
                    if basis.reduce(signature)[0]:
                        basis, _ = basis.append(signature)
                        xs.append(str(xx))
                    if basis.rank == local.point_kummer_dimension:
                        break
            require(basis.rank == local.point_kummer_dimension,
                    "local image remains incomplete at the declared witness cap")
            return {"basis_x": xs}
        record = self.arithmetic.store.discover("curve/local-kummer-image",
            {"context": context.key, "prime": place}, discover, version="z-model-witness-1")
        return self._finite_map(classes, place, record["basis_x"])

    def verify_local(self, context, classes, place, record):
        expected = (self._real_map(classes) if place == "infinity" else
                    self._finite_map(classes, place, record["basis_x"]))
        require(record == expected, "local map witness mismatch")
        return True

    def _cover_forms(self, beta):
        ring = PolynomialRing(QQ, "a,b,c")
        gamma = ring.gens()
        return ring, multiply(tuple(QQ(c) for c in self.coordinates(beta)),
                              multiply(gamma, gamma, self.polynomial), self.polynomial)

    def cover(self, context, classes, mask):
        def discover():
            beta = self.beta(classes, mask)
            ring, forms = self._cover_forms(beta)
            alpha = list(map(QQ, context.curve_generator_in_algebra))
            conic = alpha[2]*forms[1]-alpha[1]*forms[2]
            gram = matrix(QQ, 3, 3, lambda i, j: conic.derivative(ring.gen(i)).derivative(ring.gen(j))/2)
            integral = matrix(ZZ, gram/QQ(pari.content(gram)))
            factors = self.arithmetic.factor_integer(integral.det(), discover=True)
            signed = ([[-1, 1]] if integral.det() < 0 else []) + factors
            fmat = pari.matrix(len(signed), 2, [c for row in signed for c in row])
            point = pari.qfsolve([pari(integral), fmat])
            require(point.type() == "t_COL", "no rational conic parametrization witness")
            parametrization = matrix(QQ, pari.qfparam(pari(integral), point, 1))
            parametrization /= QQ(pari.content(parametrization))
            x = self.R.gen()
            gamma = parametrization*vector([x*x, x, 1])
            coefficients = multiply(tuple(QQ(c) for c in self.coordinates(beta)),
                                    multiply(gamma, gamma, self.polynomial), self.polynomial)
            index = 1 if alpha[1] else 2
            raw = self.R(-coefficients[index]/alpha[index])
            transform = matrix(QQ, [[1, 0], [0, 1]])
            if raw.degree() < 4:
                k = next(k for k in range(5) if raw(k))
                transform = matrix(QQ, [[k, 1], [1, 0]])
                s, t = transform*vector([x, 1])
                raw = self.R(sum(raw[i]*s**i*t**(4-i) for i in range(5)))
            ii, jj = invariants(raw)
            scale = (QQ(self.I/ii).nth_root(4) if self.I else QQ(self.J/jj).nth_root(6))
            quartic = raw*scale**2
            require(invariants(quartic) == (self.I, self.J), "cover normalization mismatch")
            phi = -3*self.alpha-self.a2
            cubic_invariant = (4*quartic[4]*phi+3*quartic[3]**2-8*quartic[4]*quartic[2])/3
            class_root = self.arithmetic.square_root(context.two_torsion,
                self.coordinates(cubic_invariant/beta), discover=True)
            require(class_root is not None, "quartic cubic invariant does not represent the requested class")
            return {"mask": mask, "beta": self.coordinates(beta),
                    "parametrization": [[str(c) for c in row] for row in parametrization.rows()],
                    "parameter_transform": [[str(c) for c in row] for row in transform.rows()],
                    "d_over_quartic_y": str(1/scale),
                    "cubic_invariant_over_beta_square_root": list(class_root),
                    "quartic": [str(quartic[i]) for i in range(5)]}
        return self.arithmetic.store.discover("subspace/explicit-cover",
            {"context": context.key, "classes": classes.key, "mask": mask}, discover, version="quadrics-2")

    def verify_cover(self, context, classes, mask, record):
        beta = self.beta(classes, mask)
        require(record["mask"] == mask and record["beta"] == self.coordinates(beta), "cover class mismatch")
        q = self.R(list(map(QQ, record["quartic"])))
        require(q.degree() == 4 and invariants(q) == (self.I, self.J), "cover invariant mismatch")
        parametrization = matrix(QQ, record["parametrization"])
        transform = matrix(QQ, record["parameter_transform"])
        require(parametrization.dimensions() == (3, 3) and parametrization.det() != 0,
                "degenerate conic parametrization")
        require(transform.dimensions() == (2, 2) and transform.det() != 0, "degenerate chart")
        s, t = transform*vector([self.R.gen(), 1])
        gamma = parametrization*vector([s*s, s*t, t*t])
        coefficients = multiply(tuple(QQ(c) for c in self.coordinates(beta)),
                                multiply(gamma, gamma, self.polynomial), self.polynomial)
        alpha = list(map(QQ, context.curve_generator_in_algebra))
        scale = QQ(record["d_over_quartic_y"])
        require(scale != 0 and all(coefficients[i]+alpha[i]*scale**2*q == 0 for i in (1, 2)),
                "quartic map fails the two-quadric identities")
        phi = -3*self.alpha-self.a2
        cubic_invariant = (4*q[4]*phi+3*q[3]**2-8*q[4]*q[2])/3
        root = self.field_element(record["cubic_invariant_over_beta_square_root"])
        require(root != 0 and root**2*beta == cubic_invariant, "incorrect quartic squareclass transport")
        return True

    def point_from_cover(self, context, classes, mask, record, point):
        """Exact quartic point -> input curve, through both retained transports."""
        self.verify_cover(context, classes, mask, record)
        U, V, Y = map(QQ, point)
        q = list(map(QQ, record["quartic"]))
        require((U or V) and Y**2 == sum(q[i]*U**i*V**(4-i) for i in range(5)),
                "invalid quartic point")
        d = QQ(record["d_over_quartic_y"])*Y
        if not d:
            return None  # the elliptic identity, not a new affine direction
        s, t = matrix(QQ, record["parameter_transform"])*vector([U, V])
        gamma = self.K(list(matrix(QQ, record["parametrization"])*vector([s*s, s*t, t*t])))
        beta = self.beta(classes, mask)
        value = beta*gamma**2
        alpha = list(map(QQ, context.curve_generator_in_algebra))
        z = value.lift()[0]/d**2+alpha[0]
        norm_beta = QQ(pari.nfeltnorm(self.nf, self.element(self.coordinates(beta))))
        norm_gamma = QQ(pari.nfeltnorm(self.nf, self.element(self.coordinates(gamma))))
        Yz = norm_beta.sqrt()*norm_gamma/d**3
        require(Yz**2 == self.cubic(z), "invalid norm-cover recovery")
        a1, _, a3, _, _ = map(QQ, context.minimal_model.coefficients)
        xm = z/4
        ym = (Yz-4*a1*xm-4*a3)/8
        u, r, s, t = map(QQ, context.minimal_to_input)
        recovered = (str(u**2*xm+r), str(u**3*ym+s*u**2*xm+t))
        require(context.model.contains(recovered), "cover point fails exact model transport")
        return recovered

    def cover_map(self, context, classes, mask, record, *, y_denominator=1):
        """Rational functions for a retained cover, in the exact input model."""
        self.verify_cover(context, classes, mask, record)
        ring = PolynomialRing(QQ, "x,y")
        x, y = ring.fraction_field().gens()
        s, t = matrix(QQ, record["parameter_transform"])*vector([x, 1])
        gamma = tuple(matrix(QQ, record["parametrization"])*vector([s*s, s*t, t*t]))
        beta = self.beta(classes, mask)
        beta_coordinates = tuple(QQ(c) for c in self.coordinates(beta))
        value = multiply(beta_coordinates, multiply(gamma, gamma, self.polynomial), self.polynomial)
        multiplication = matrix(ring.fraction_field(), 3, 3,
            lambda i, j: multiply(gamma, tuple(int(k == j) for k in range(3)), self.polynomial)[i])
        norm_beta = QQ(pari.nfeltnorm(self.nf, self.element(beta_coordinates)))
        d = QQ(record["d_over_quartic_y"])*y/QQ(y_denominator)
        alpha0 = QQ(context.curve_generator_in_algebra[0])
        z, Yz = value[0]/d**2+alpha0, norm_beta.sqrt()*multiplication.det()/d**3
        a1, _, a3, _, _ = map(QQ, context.minimal_model.coefficients)
        xm = z/4
        ym = (Yz-4*a1*xm-4*a3)/8
        u, r, s, t = map(QQ, context.minimal_to_input)
        return u**2*xm+r, u**3*ym+s*u**2*xm+t

    def search_cover(self, state, context, classes, mask, limits, *, cache, height):
        """Direct enumeration of one allowed cover inside the caller's supervisor.

        Coordinate normalization and chart choice remain explicit separate
        stages. A bounded miss is merely a search observation.
        """
        from .supervisor import Limits
        require(isinstance(limits, Limits) and type(height) is int and height > 0,
                "point enumeration requires declared resource and height bounds")
        cover = self.cover(context, classes, mask)
        self.verify_cover(context, classes, mask, cover)
        q = self.R(list(map(QQ, cover["quartic"])))
        denominator = q.denominator()
        raw_points = pari.hyperellratpoints(q*denominator**2, height)
        points = [(QQ(p[0]), QQ(1), QQ(p[1])/denominator) for p in raw_points]
        if q[4].is_square():
            points.append((QQ(1), QQ(0), q[4].sqrt()))
        for point in points:
            recovered = self.point_from_cover(context, classes, mask, cover, point)
            if recovered is not None:
                state = state.adjoin(recovered, cache=cache)
        return state

    def _pair_support(self, quartics, *, retained=None):
        # Remark 3.3 also requires integrality. Include denominator primes for
        # rational quartics; gamma is normalized primitive integral in fisher_gamma.
        q, other, _ = quartics
        values = [abs(ZZ(q.discriminant().numerator())), ZZ(q.discriminant().denominator()),
                  ZZ(q.denominator()), abs(ZZ(other[4].numerator())), ZZ(other[4].denominator())]
        facts = []
        support = {2, 3, 5, 7}
        if retained is not None:
            require(len(retained) == len(values), "incomplete pairing support witness")
        for i, value in enumerate(values):
            require(value != 0, "degenerate pairing support")
            fac = self.arithmetic.factor_integer(value, discover=True) if retained is None else retained[i]
            require(all(type(e) is int and e > 0 and ZZ(p).is_prime(proof=True) for p, e in fac)
                    and len({p for p, _ in fac}) == len(fac)
                    and prod(ZZ(p)**e for p, e in fac) == value, "invalid support factorization")
            support.update(int(p) for p, _ in fac)
            facts.append(fac)
        return sorted(support)+["infinity"], facts

    def _pair(self, classes, masks, records, *, retained=None):
        quartics = [self.R(list(map(QQ, row["quartic"]))) for row in records]
        if retained is None:
            phi = -3*self.alpha-self.a2
            coordinates = matrix(QQ, 3, 3, lambda i, j: (phi**j).lift()[i])
            z = prod((4*q[4]*phi+3*q[3]**2-8*q[4]*q[2])/3 for q in quartics)
            root = self.arithmetic.square_root(self.context.two_torsion, self.coordinates(z), discover=True)
            require(root is not None, "pairing product is not a square")
            sqrt_phi = list(map(str, coordinates.inverse()*vector(QQ, root)))
        else:
            require(retained["masks"] == list(masks), "misbound pairing")
            sqrt_phi = retained["square_root_phi_coefficients"]
        gamma = self.R(fisher_gamma(quartics, sqrt_phi, self.I, self.J))
        places, factors = self._pair_support(quartics, retained=None if retained is None else retained["support_factors"])
        terms, answer = [], 1
        if retained is not None:
            require([term["place"] for term in retained["local_terms"]] == places, "incomplete local pairing")
        for i, place in enumerate(places):
            xx = (quartic_local_witness(quartics[0], gamma, place, node_cap=self.cap)
                  if retained is None else QQ(retained["local_terms"][i]["x"]))
            qv, gv = quartics[0](xx), gamma(xx)
            require(qv != 0 and gv != 0 and (qv > 0 if place == "infinity" else local_square(qv, place)),
                    "invalid local pairing witness")
            symbol = int(hilbert_symbol(quartics[1][4], gv, place))
            answer *= symbol
            terms.append({"place": place, "x": str(xx), "q_value": str(qv),
                          "gamma_value": str(gv), "hilbert_symbol": symbol})
        result = {"masks": list(masks), "square_root_phi_coefficients": sqrt_phi,
                  "gamma": [str(gamma[i]) for i in range(3)], "support_factors": factors,
                  "local_terms": terms, "value": int(answer == -1)}
        if retained is not None:
            require(result == retained, "pairing witness mismatch")
        return result

    def ct_pairing(self, context, classes, masks, covers):
        lookup = dict(zip(masks, covers))
        pairs = []
        matrix_rows = [[0]*len(masks) for _ in masks]
        for i, j in combinations(range(len(masks)), 2):
            summed = masks[i] ^ masks[j]
            if summed not in lookup:
                lookup[summed] = self.cover(context, classes, summed)
            self.verify_cover(context, classes, summed, lookup[summed])
            trio = (masks[i], masks[j], summed)
            pair = self.arithmetic.store.discover("subspace/fisher-pairing",
                {"context": context.key, "classes": classes.key,
                 "covers": [digest(lookup[m]) for m in trio]},
                lambda: self._pair(classes, trio, [lookup[m] for m in trio]), version="fisher-2022-1")
            pairs.append(pair)
            matrix_rows[i][j] = matrix_rows[j][i] = pair["value"]
        return {"matrix": matrix_rows, "pairs": pairs,
                "sum_covers": [lookup[m] for m in sorted(set(lookup)-set(masks))]}

    def verify_ct(self, context, classes, masks, covers, record):
        lookup = dict(zip(masks, covers))
        for cover in record["sum_covers"]:
            require(cover["mask"] not in lookup, "duplicate sum cover")
            self.verify_cover(context, classes, cover["mask"], cover)
            lookup[cover["mask"]] = cover
        expected_sums = {a ^ b for a, b in combinations(masks, 2)} - set(masks)
        require(set(lookup)-set(masks) == expected_sums, "incomplete or extraneous sum covers")
        indices = list(combinations(range(len(masks)), 2))
        require(len(record["pairs"]) == len(indices), "incomplete CT matrix")
        rows = [[0]*len(masks) for _ in masks]
        for (i, j), pair in zip(indices, record["pairs"]):
            trio = (masks[i], masks[j], masks[i] ^ masks[j])
            value = self._pair(classes, trio, [lookup[m] for m in trio], retained=pair)["value"]
            rows[i][j] = rows[j][i] = value
        require(rows == record["matrix"], "incorrect CT matrix")
        return True
