"""Target-blind, costed union-of-boxes policy for pointed quartics.

The probabilities are scheduling surrogates, never local/global solubility
certificates. Completed misses remove surrogate exposure, not rational points.
No oracle, corpus or endpoint files are read by this module.
"""
from fractions import Fraction as F
from itertools import product
from math import gcd, sqrt

import numpy as np

from finite_cancellation_features import ev, neighbour_maps, xmap
from finite_cancellation_validation_features import q_tree
from search_observability import transform

HEIGHTS = (8000, 32000, 125000)


def root_distribution(q, p, roots):
    """q-only differential mass; unresolved branches remain explicit.

    Each category is a *base-coordinate* residue ball. Keeping one partition
    per prime prevents incompatible same-prime balls being multiplied.
    """
    tree = q_tree(q, p, depth=4, cap=512)
    masses = {r: 0.0 for r in roots}
    masses[None] = 0.0
    for leaf in tree['leaves']:
        if leaf['status'] != 'square':
            continue
        root = leaf['residue'] % p if leaf['kind'] == 'affine' else None
        category = root if root in roots else None
        masses[category] += float(F(leaf['mass'])) * p**(leaf['vq']/2)
    total = sum(masses.values())
    uniform = {r: 1/(p+1) for r in roots}
    uniform[None] = 1-sum(uniform.values())
    probabilities = {r: w/total for r, w in masses.items()} if total else uniform
    return probabilities, tree


def real_grid(q, size=32):
    """Fixed midpoint quadrature for |m dn-n dm|/sqrt(q(m,n)).

    The two square-boundary charts cover RP1 once up to endpoints. This is a
    numerical feature, not a certified integral or a global point measure.
    """
    scale = max(map(abs, q))
    coefficients = [a/scale for a in q]
    rows = []
    for i in range(size):
        t = -1+(2*i+1)/size
        for m, n in ((t, 1.0), (1.0, t)):
            value = ev(coefficients, m, n)
            if value > 0:
                rows.append((m, n, 1/sqrt(value)))
    if not rows:
        # No theorem follows from a quadrature miss; use the full boundary.
        rows = [(0., 1., 1.), (1., 0., 1.)]
    total = sum(r[2] for r in rows)
    return [(m, n, w/total) for m, n, w in rows]


def neighbour_law(n, d, matrix, p):
    """Exact jointly primitive content, with no support factorization."""
    T = tuple(map(F, matrix))
    if any(x.denominator != 1 for x in T) or abs(T[0]*T[3]-T[1]*T[2]) != p:
        raise ArithmeticError('not an integral prime neighbour')
    values = transform(n, T)+transform(d, T)
    content = gcd(*(int(v) for v in values))
    e = 0
    while content % p == 0:
        content //= p
        e += 1
    if content != 1 or not 0 <= e <= 4:
        raise ArithmeticError('prime-neighbour content law failed')
    return e


def prepare(curve, anchor, mapper, local_weight):
    """All candidate creation is charged, including unselected models."""
    old = mapper.mapping(tuple(map(F, curve)), [tuple(map(F, anchor))],
                         {'representative': [1]})
    neighbours, attempts = neighbour_maps(old, mapper.pari)
    n, d, q = xmap(F(curve[3]), F(curve[4]), anchor, old)
    models = [{'name': 'factor_free', 'mapping': old}] + neighbours
    roots = {}
    for model in neighbours:
        p, r = model['prime'], model['residue']
        model['content_exponent'] = neighbour_law(n, d, model['mapping']['second_matrix'], p)
        roots.setdefault(p, set()).add(r)
    partitions = []
    for p, rs in sorted(roots.items()):
        rs = sorted(rs)
        uniform = {r: 1/(p+1) for r in rs}
        uniform[None] = 1-sum(uniform.values())
        if local_weight:
            differential, tree = root_distribution(q, p, rs)
        else:
            differential, tree = uniform, None
        probabilities = {r: (1-local_weight)*uniform[r]+local_weight*differential[r]
                         for r in uniform}
        partitions.append({'prime': p, 'categories': [
            {'residue': r, 'probability': probabilities[r]} for r in rs+[None]],
            'tree': tree})
    # Empty product is one state. Same-prime categories are mutually exclusive.
    states = []
    for cells in product(*(p['categories'] for p in partitions)):
        weight = 1.
        residues = {}
        for partition, cell in zip(partitions, cells):
            weight *= cell['probability']
            residues[partition['prime']] = cell['residue']
        states.append((residues, weight))
    ratios, weights = [], []
    for m, v, rw in real_grid(q):
        for residues, pw in states:
            row = [1.]
            for model in neighbours:
                a, b, c, z = map(F, model['mapping']['second_matrix'])
                # z_new = adj(T) w / gcd(adj(T)w). The gcd is p in the
                # distinguished old ball and 1 outside: exact for primitive w.
                g = model['prime'] if residues[model['prime']] == model['residue'] else 1
                row.append(max(abs(float(z)*m-float(b)*v), abs(-float(c)*m+float(a)*v))/g)
            ratios.append(row)
            weights.append(rw*pw)
    return {'curve': curve, 'anchor': anchor, 'models': models, 'attempts': attempts,
            'partitions': partitions, 'ratios': ratios, 'weights': weights,
            'local_weight': local_weight,
            'boundary': 'Finite-grid product-measure exposure surrogate; unknown local mass is not pruned. Models search complete overlapping boxes, not exclusive residue classes.'}


def prior_for(fit, index):
    bucket = 0 if index < 4 else 1 if index < 16 else 2
    return fit['anchor_priors'][bucket]


class AnchorExposure:
    def __init__(self, prepared, fit, index):
        self.index, self.prepared = index, prepared
        self.prior = prior_for(fit, index)
        self.weights = np.array(prepared['weights'])
        self.covered = np.zeros(len(self.weights))
        radial = np.array(fit['radial_heights'])
        ratios = np.array(prepared['ratios'])
        self.jobs, self.exposures = [], []
        for mi in range(len(prepared['models'])):
            for h in HEIGHTS:
                self.jobs.append((mi, h))
                # Empirical conditional radius distribution from certified
                # development hits, not literal corpus height regressions.
                self.exposures.append(np.searchsorted(radial, h/ratios[:, mi], side='right')/len(radial))
        self.exposures = np.array(self.exposures)
        self.attempted = set()

    def options(self, cost):
        previous = float(self.covered @ self.weights)
        increments = np.maximum(self.exposures-self.covered, 0) @ self.weights
        for j, ((mi, height), mass) in enumerate(zip(self.jobs, increments)):
            if j in self.attempted or mass <= 1e-12:
                continue
            probability = self.prior*float(mass)/(1-self.prior*previous)
            price = cost.estimate(height, probability)
            yield (probability/price, self.index, j, mi, height, probability, price)

    def observe(self, job, complete):
        self.attempted.add(job)
        if complete:
            self.covered = np.maximum(self.covered, self.exposures[job])
        # A timeout has no completed prefix; its exposure remains uncovered.


class CostModel:
    def __init__(self, fit):
        self.reference = fit['backend_cpu_at_125000']
        self.samples = 0
        self.fixed = fit['per_call_overhead_cpu']
        self.certificate = fit['independent_certificate_cpu_prior']

    def estimate(self, height, probability=0):
        return self.fixed + self.reference*(height/125000)**2 + probability*self.certificate

    def observe(self, height, seconds, complete):
        if complete:
            normalized = max(.001, seconds-self.fixed)/(height/125000)**2
            # Shrink measured complete-call cost toward four development calls.
            self.reference = (self.reference*(4+self.samples)+normalized)/(5+self.samples)
            self.samples += 1


def choose(anchors, cost, fit, next_index, count):
    options = [o for anchor in anchors.values() for o in anchor.options(cost)]
    best = max(options, key=lambda o:(o[0], -o[1], -o[2])) if options else None
    probability = prior_for(fit, next_index)*sum(h <= HEIGHTS[0] for h in fit['radial_heights'])/len(fit['radial_heights'])
    # A new anchor includes a fixed development preparation allowance. Actual
    # preparation is always charged; this prior never supplies a CPU receipt.
    fresh_score = probability/(cost.estimate(HEIGHTS[0], probability)+.008)
    if next_index < count and (best is None or fresh_score > best[0]):
        return {'kind':'prepare', 'index':next_index, 'score':fresh_score}
    if best is None:
        return {'kind':'exhausted'}
    score, index, job, model, height, probability, price = best
    return {'kind':'search', 'index':index, 'job':job, 'model':model, 'height':height,
            'score':score, 'probability':probability, 'predicted_cpu':price}
