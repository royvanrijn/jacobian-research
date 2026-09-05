"""Lazy exact best-first lattice enumeration and distinct parity cosets.

A partial LDL assignment supplies an admissible lower bound. Each node keeps
one lazy stream of integer children, so neither integer boxes nor 2^r parity
subsets are built in advance. The first vector of each parity is a closest
representative of that coset. This is an exact ordering, not a polynomial-time
CVP claim; ill-conditioned lattices can still require many visited nodes.
"""
from dataclasses import dataclass
from fractions import Fraction
import heapq

from .store import digest


def rational(value):
    return Fraction(str(value))


def ldl(gram):
    gram = tuple(tuple(map(rational, row)) for row in gram)
    n = len(gram)
    if any(len(row) != n for row in gram) or any(gram[i][j] != gram[j][i] for i in range(n) for j in range(n)):
        raise ValueError("CVP requires a symmetric square metric")
    lower = [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    diagonal = []
    for i in range(n):
        diagonal.append(gram[i][i]-sum(lower[i][k]**2*diagonal[k] for k in range(i)))
        if diagonal[-1] <= 0:
            raise ValueError("CVP metric must be positive definite")
        for j in range(i+1, n):
            lower[j][i] = (gram[j][i]-sum(lower[j][k]*lower[i][k]*diagonal[k] for k in range(i)))/diagonal[i]
    return gram, lower, diagonal


@dataclass(frozen=True)
class Hole:
    mask: int
    doubled_coordinates: tuple[int, ...]
    squared_distance: Fraction

    @property
    def coordinates(self):
        return tuple(Fraction(z, 2) for z in self.doubled_coordinates)

    def record(self):
        return {"mask": self.mask, "doubled_coordinates": list(self.doubled_coordinates),
                "squared_distance": str(self.squared_distance)}

    @classmethod
    def from_record(cls, row):
        return cls(row["mask"], tuple(row["doubled_coordinates"]), rational(row["squared_distance"]))


class VoronoiIterator:
    """Request successive unseen M/2M centres in a declared metric.

    ``target`` and returned centres use the original lattice coordinates.
    Zero is omitted by default. A positive node budget yields a resumable
    TimeoutError, never a claim that a coset has no representative.
    """
    def __init__(self, gram, *, target=None, seen=(), include_zero=False, binding=None):
        self.gram, self.lower, self.diagonal = ldl(gram)
        self.rank = len(self.gram)
        self.target = tuple(map(rational, target if target is not None else [0]*self.rank))
        if len(self.target) != self.rank:
            raise ValueError("target dimension differs from the metric")
        self.binding = binding
        self.seen = set(seen)
        if any(type(mask) is not int or not 0 <= mask < 1 << self.rank for mask in self.seen):
            raise ValueError("invalid seen parity")
        if not include_zero:
            self.seen.add(0)
        self.heap, self.pending = [], []
        self.serial, self.visited_nodes = 0, 0
        if self.rank:
            self._start(())
        elif not self.seen:
            self.pending.append(Hole(0, (), Fraction(0)))

    @property
    def identity(self):
        return {"gram": [[str(x) for x in row] for row in self.gram],
                "target": list(map(str, self.target)), "binding": self.binding}

    def _parameters(self, suffix):
        i = self.rank-1-len(suffix)
        residual = [Fraction(z)-2*self.target[i+1+j] for j, z in enumerate(suffix)]
        centre = 2*self.target[i]-sum(self.lower[i+1+j][i]*z for j, z in enumerate(residual))
        cost = Fraction(0)
        for j, z in enumerate(residual):
            k = i+1+j
            cost += self.diagonal[k]*(z+sum(self.lower[k+1+h][k]*v for h, v in enumerate(residual[j+1:])))**2
        return i, centre, cost

    def _start(self, suffix):
        _, centre, _ = self._parameters(suffix)
        floor = centre.numerator//centre.denominator
        self._push(suffix, floor, floor+1)

    def _push(self, suffix, left, right):
        i, centre, base = self._parameters(suffix)
        z = left if abs(left-centre) <= abs(right-centre) else right
        cost = (base+self.diagonal[i]*(z-centre)**2)/4
        self.serial += 1
        # Depth breaks equal-cost ties without expanding a whole sphere first.
        heapq.heappush(self.heap, (cost, -len(suffix), self.serial, tuple(suffix), left, right))

    def __iter__(self):
        return self

    def __next__(self):
        return self.next_hole()

    def next_hole(self, *, node_budget=None):
        if self.pending:
            return self.pending.pop(0)
        if node_budget is not None and (type(node_budget) is not int or node_budget < 1):
            raise ValueError("positive integer node budget required")
        start = self.visited_nodes
        while len(self.seen) < 1 << self.rank and self.heap:
            if node_budget is not None and self.visited_nodes-start >= node_budget:
                raise TimeoutError("CVP node budget reached; retain checkpoint and resume")
            cost, _, _, suffix, left, right = heapq.heappop(self.heap)
            _, centre, _ = self._parameters(suffix)
            z = left if abs(left-centre) <= abs(right-centre) else right
            self._push(suffix, left-1 if z == left else left, right+1 if z == right else right)
            self.visited_nodes += 1
            vector = (z, *suffix)
            if len(vector) < self.rank:
                self._start(vector)
                continue
            mask = sum((entry % 2) << i for i, entry in enumerate(vector))
            if mask not in self.seen:
                self.seen.add(mask)
                return Hole(mask, vector, cost)
        raise StopIteration

    def next_holes(self, count, *, diversity_window=1, previous=(), node_budget=None, allowed=None):
        """Choose from a bounded nearest-centre window using metric separation.

        Window=1 preserves global distance order. Larger windows greedily favor
        distance from earlier chosen centres inside that window; unselected
        holes remain pending. ``allowed`` is queried again on emission so new
        theorems can invalidate a saved queue. Rejected masks are recorded.
        """
        if type(count) is not int or count < 0 or type(diversity_window) is not int or diversity_window < 1:
            raise ValueError("nonnegative count and positive diversity window required")
        pool, result = [], []
        start = self.visited_nodes
        def separation(a, b):
            delta = [x-y for x, y in zip(a.coordinates, b.coordinates)]
            return sum(delta[i]*self.gram[i][j]*delta[j] for i in range(self.rank) for j in range(self.rank))
        try:
            for _ in range(count*diversity_window):
                remaining = None if node_budget is None else node_budget-(self.visited_nodes-start)
                if remaining is not None and remaining <= 0:
                    raise TimeoutError("CVP batch budget reached")
                try:
                    hole = self.next_hole(node_budget=remaining)
                except StopIteration:
                    break
                pool.append(hole)
            pool = [hole for hole in pool if allowed is None or allowed(hole)]
            for _ in range(min(count, len(pool))):
                anchors = [*previous, *result]
                if diversity_window == 1 or not anchors:
                    index = min(range(len(pool)), key=lambda i: (pool[i].squared_distance, pool[i].mask))
                else:
                    index = max(range(len(pool)), key=lambda i: (min(separation(pool[i], p) for p in anchors),
                                -pool[i].squared_distance, -pool[i].mask))
                result.append(pool.pop(index))
        finally:
            self.pending = sorted([*self.pending, *pool], key=lambda h: (h.squared_distance, h.mask))
        return result

    def checkpoint(self):
        row = {"schema": "elliptic-curves.lazy-cvp.v1", "identity": self.identity,
               "seen": sorted(self.seen), "serial": self.serial, "visited_nodes": self.visited_nodes,
               "heap": [[str(cost), depth, serial, list(suffix), left, right]
                        for cost, depth, serial, suffix, left, right in self.heap],
               "pending": [h.record() for h in self.pending]}
        return {**row, "sha256": digest(row)}

    @classmethod
    def resume(cls, record, *, binding=None):
        row = dict(record); expected = row.pop("sha256")
        if digest(row) != expected or row["schema"] != "elliptic-curves.lazy-cvp.v1":
            raise ValueError("corrupt CVP checkpoint")
        identity = row["identity"]
        if identity["binding"] != binding:
            raise ValueError("CVP checkpoint belongs to a different MWState or policy")
        obj = cls(identity["gram"], target=identity["target"], seen=row["seen"], include_zero=True, binding=binding)
        obj.serial, obj.visited_nodes = row["serial"], row["visited_nodes"]
        obj.heap = [(rational(cost), depth, serial, tuple(suffix), left, right)
                    for cost, depth, serial, suffix, left, right in row["heap"]]
        heapq.heapify(obj.heap)
        obj.pending = [Hole.from_record(h) for h in row["pending"]]
        return obj
