"""Cached local characters evaluated on a whole squareclass matrix."""

from .binary import BinaryBasis, unpack


class LocalSquareclasses:
    """One fixed local coordinate map for every class and related curve.

    At odd primes, valuation and residue characters suffice. At two, unit
    characters come from O_K/P^(2e+1); the omitted principal units are squares
    by Hensel's lemma. Only even cyclic factors contribute to units modulo
    squares. No class group or global unit basis is used.
    """

    def __init__(self, nf, p, *, arithmetic=None, context=None):
        from sage.all import pari
        self.pari, self.nf, self.p = pari, nf, int(p)
        self._signatures = {}
        self._arithmetic, self._context = arithmetic, context

        def build():
            primes = (arithmetic.prime_ideals(context, self.p, discover=True)
                      if arithmetic else pari.idealprimedec(nf, self.p))
            data = []
            for prime in primes:
                uniformizer = pari.nfbasistoalg(nf, pari.idealappr(nf, prime))
                if int(pari.idealval(nf, uniformizer, prime)) != 1:
                    raise ArithmeticError("invalid local uniformizer")
                if self.p == 2:
                    bid = pari.idealstar(nf, pari.idealpow(nf, prime, 2*int(prime[2])+1), 2)
                    data.append([prime, uniformizer, bid])
                else:
                    data.append([prime, uniformizer, pari.nfmodprinit(nf, prime)])
            return data

        if arithmetic:
            if context is None:
                raise ValueError("cached local arithmetic requires a labelled algebra")
            record = arithmetic.store.discover("two-torsion/local-characters",
                {"algebra": context.key, "prime": self.p},
                lambda: {"binary": arithmetic._blob(pari(build()))}, version=arithmetic.version + "/characters-1")
            self.data = list(arithmetic._restore(record["binary"]))
        else:
            self.data = build()
        self.primes = [row[0] for row in self.data]

    @property
    def point_kummer_dimension(self):
        return len(self.primes) - 1 + (1 if self.p == 2 else 0)

    def signature(self, value):
        pari, nf = self.pari, self.nf
        if value == 0:
            raise ValueError("zero is not a squareclass")
        key = str(pari.nfalgtobasis(nf, value))
        if key not in self._signatures:
            def build():
                bits = []
                for prime, uniformizer, local_map in self.data:
                    valuation = int(pari.idealval(nf, value, prime))
                    bits.append(valuation % 2)
                    unit = value / uniformizer**valuation
                    if self.p == 2:
                        logs = pari.ideallog(nf, unit, local_map)
                        cyclic = local_map.bid_get_cyc()
                        bits.extend(int(logs[i]) % 2 for i in range(len(cyclic)) if int(cyclic[i]) % 2 == 0)
                    else:
                        bits.append(0 if pari.issquare(pari.nfmodpr(nf, unit, local_map)) else 1)
                return {"bits": bits}
            if self._arithmetic:
                row = self._arithmetic.store.discover("two-torsion/local-class",
                    {"algebra": self._context.key, "prime": self.p, "element": key}, build,
                    version=self._arithmetic.version + "/characters-1")
            else:
                row = build()
            self._signatures[key] = tuple(row["bits"])
        return self._signatures[key]

    def is_square(self, value):
        return not any(self.signature(value))

    def coordinates(self, elements):
        from .binary import pack
        if not elements:
            return [], []
        signatures = [self.signature(element) for element in elements]
        basis = BinaryBasis(len(signatures[0]))
        representatives, coordinates = [], []
        for element, signature in zip(elements, signatures):
            residual, combination = basis.reduce(pack(signature))
            if residual:
                # Only independent representatives enter this basis, so its
                # provenance coordinates refer precisely to representatives.
                basis, _ = basis.append(pack(signature))
                combination = 1 << len(representatives)
                representatives.append(element)
            coordinates.append(combination)
        return representatives, [list(unpack(row, len(representatives))) for row in coordinates]
