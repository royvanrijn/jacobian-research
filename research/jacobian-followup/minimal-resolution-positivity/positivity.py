"""Necessary minimal-resolution tests; numerical acceptance is not realization.

The geometric inputs (a complete boundary, an actual morphism and its Keller
open) must be established outside this module. Blowdowns use the smooth
surface formula and retain arithmetic genera after the boundary ceases to be
SNC. No singular Stein-surface Schur complement replaces the minimal smooth
resolution here.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import sympy as S


@dataclass(frozen=True)
class Model:
    names: tuple[str, ...]
    intersection: S.Matrix
    canonical_intersections: S.Matrix
    hyperplane_intersections: S.Matrix

    @classmethod
    def from_boundary(cls, names: Sequence[str], matrix, genera, poles):
        m = S.Matrix(matrix)
        n = len(names)
        if m.shape != (n, n) or len(genera) != n or len(poles) != n:
            raise ValueError('boundary dimensions disagree')
        if m != m.T or m.det() == 0:
            raise ValueError('intersection matrix must be symmetric and invertible')
        if any(x.is_Integer is not True for x in m):
            raise ValueError('integral intersection matrix required')
        if any(m[i, j] < 0 for i in range(n) for j in range(n) if i != j):
            raise ValueError('distinct curve intersections must be nonnegative')
        if len(set(names)) != n or any(not isinstance(g, int) or g < 0 for g in genera):
            raise ValueError('distinct names and nonnegative arithmetic genera required')
        p = S.Matrix(poles)
        if any(x.is_Integer is not True or x < 0 for x in p):
            raise ValueError('effective integral target-line divisor required')
        k = S.Matrix([2*g - 2 - m[i, i] for i, g in enumerate(genera)])
        result = cls(tuple(names), m, k, m*p)
        result.validate()
        return result

    @property
    def rank(self):
        return len(self.names)

    @property
    def canonical(self):
        return self.intersection.inv()*self.canonical_intersections

    @property
    def poles(self):
        return self.intersection.inv()*self.hyperplane_intersections

    @property
    def genera(self):
        return tuple(S.cancel((self.intersection[i, i] + self.canonical_intersections[i])/2 + 1)
                     for i in range(self.rank))

    def pairing(self, a, b):
        return S.cancel((a.T*self.intersection.inv()*b)[0])

    @property
    def invariants(self):
        k, h = self.canonical_intersections, self.hyperplane_intersections
        n, kh, kk = self.pairing(h, h), self.pairing(k, h), self.pairing(k, k)
        r, q = k + 3*h, 2*k + 3*h
        return dict(rank=self.rank, degree=n, genus=S.cancel(1+(n+kh)/2),
                    canonical_square=kk, ramification_square=self.pairing(r, r),
                    adjoint_square=self.pairing(q, q),
                    ramification_arithmetic_genus=S.cancel(1+self.pairing(r, k+r)/2))

    def validate(self):
        if self.rank == 0 or self.intersection.shape != (self.rank, self.rank):
            raise ValueError('nonempty boundary basis required')
        if any(x < 0 for x in self.hyperplane_intersections):
            raise ValueError('target pullback is not nef on the boundary')
        if any(x.is_Integer is not True or x < 0 for x in self.genera):
            raise ValueError('invalid arithmetic genus after contraction')
        if abs(self.intersection.det()) != 1:
            raise ValueError('whole affine-plane boundary must be unimodular')
        if any(x.is_Integer is not True for x in self.canonical):
            raise ValueError('canonical boundary class is not integral')
        if any(x.is_Integer is not True or x < 0 for x in self.poles):
            raise ValueError('target-line class must be effective and integral')
        if self.invariants['canonical_square'] + self.rank != 10:
            raise ValueError('rational-surface canonical identity fails')
        if self.invariants['degree'] <= 0:
            raise ValueError('positive generic degree required')

    def blow_down(self, index):
        if (self.intersection[index, index], self.genera[index],
                self.hyperplane_intersections[index]) != (-1, 0, 0):
            raise ValueError('only a smooth rational H-null minus-one curve may be removed')
        keep = [i for i in range(self.rank) if i != index]
        v = self.intersection.extract(keep, [index])
        m = self.intersection.extract(keep, keep) + v*v.T
        k = self.canonical_intersections.extract(keep, [0]) - v
        h = self.hyperplane_intersections.extract(keep, [0])
        result = Model(tuple(self.names[i] for i in keep), m, k, h)
        result.validate()
        before, after = self.invariants, result.invariants
        for key in ('degree', 'genus'):
            assert before[key] == after[key], key
        assert after['canonical_square'] == before['canonical_square'] + 1
        assert after['ramification_square'] == before['ramification_square'] + 1
        assert after['adjoint_square'] == before['adjoint_square'] + 4
        assert after['ramification_arithmetic_genus'] == before['ramification_arithmetic_genus'] + 1
        return result

    def minimalize(self, reverse=False):
        model, history = self, []
        while True:
            choices = [i for i in range(model.rank)
                       if (model.intersection[i, i], model.genera[i],
                           model.hyperplane_intersections[i]) == (-1, 0, 0)]
            if not choices:
                return model, history
            index = choices[-1] if reverse else choices[0]
            following = model.blow_down(index)
            history.append(dict(curve=model.names[index], before=model.invariants,
                                after=following.invariants))
            model = following

    def positivity_obstructions(self, include_source_corollaries=True):
        """Necessary tests on an already relatively minimal smooth model."""
        if any((self.intersection[i, i], self.genera[i], self.hyperplane_intersections[i])
               == (-1, 0, 0) for i in range(self.rank)):
            raise ValueError('minimalize before applying positivity')
        inv = self.invariants
        if inv['degree'] == 1:
            return ['OUTSIDE_SCOPE: geometric degree one']
        k, h = self.canonical_intersections, self.hyperplane_intersections
        r, q = k+3*h, 2*k+3*h
        result = []
        if inv['adjoint_square'] <= 0:
            result.append('Q_SQUARE_NOT_POSITIVE')
        if inv['ramification_square'] <= 0:
            result.append('R_SQUARE_NOT_POSITIVE')
        if inv['ramification_arithmetic_genus'] < 2:
            result.append('RAMIFICATION_ARITHMETIC_GENUS_LT_TWO')
        for i, name in enumerate(self.names):
            null = (h[i] == 0 and self.intersection[i, i] == -2 and self.genera[i] == 0)
            for label, value in (('Q', q[i]), ('R', r[i])):
                if value < 0 or (value == 0) != null:
                    result.append(f'{label}_BOUNDARY_INTERSECTION:{name}:{value}')
        if include_source_corollaries and inv['genus'] < 2:
            result.append('GENERAL_LINE_GENUS_LT_TWO')
        if include_source_corollaries and self.rank < 2:
            result.append('MINIMAL_PICARD_RANK_LT_TWO')
        return result


REQUIRED_GEOMETRY = ('complete_whole_A2_boundary', 'actual_projective_morphism',
                     'etale_on_whole_A2', 'target_line_missed')


def evaluate(data):
    """Fail closed on absent geometric evidence, before reading a partial matrix."""
    missing = [x for x in REQUIRED_GEOMETRY if data.get('geometry', {}).get(x) is not True]
    if missing:
        return dict(status='NOT_APPLICABLE', missing_or_false=missing,
                    reason=data.get('reason', 'geometric hypotheses not certified'))
    model = Model.from_boundary(data['names'], data['matrix'], data['genera'], data['poles'])
    minimal, history = model.minimalize()
    if minimal.invariants['degree'] == 1:
        return dict(status='OUTSIDE_SCOPE', reason='geometric degree one',
                    invariants=minimal.invariants, realization_certified=False)
    obstructions = minimal.positivity_obstructions()
    return dict(status='NUMERICALLY_EXCLUDED' if obstructions else 'NECESSARY_TESTS_PASS',
                minimal_names=minimal.names, invariants=minimal.invariants,
                blowdowns=history, obstructions=obstructions,
                realization_certified=False)
