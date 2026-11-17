"""Small exact permutation model and canonical tuple augmentation."""

from itertools import combinations, permutations


IDENTITY = tuple(range(8))


def compose(a, b):
    assert len(a) == len(b)
    return tuple(a[b[i]] for i in range(len(a)))


def inverse(a):
    result = [0] * len(a)
    for i, j in enumerate(a):
        result[j] = i
    return tuple(result)


def power(a, n):
    value = tuple(range(len(a)))
    for _ in range(n):
        value = compose(value, a)
    return value


def model(degree=8):
    elements = []
    for support in combinations(range(degree), 4):
        a, b, c, d = support
        for edges in [((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c))]:
            p = list(range(degree))
            for u, v in edges:
                p[u], p[v] = v, u
            elements.append(tuple(p))
    assert len(elements) == len(set(elements)) == 3 * len(list(combinations(range(degree), 4)))
    index = {p: i for i, p in enumerate(elements)}
    table = [[index[compose(compose(a, b), a)] for b in elements] for a in elements]
    tau = elements[0]
    assert tau == (1, 0, 3, 2) + tuple(range(4, degree))
    centralizer = [p for p in permutations(range(degree)) if compose(p, tau) == compose(tau, p)]
    assert len(centralizer) == (192 if degree == 8 else 960)
    actions = [[index[compose(compose(p, a), inverse(p))] for a in elements] for p in centralizer]
    return elements, table, centralizer, actions


def tuple_representatives(actions):
    # At each level choose the least member of an orbit under the stabilizer
    # of the entire retained prefix. Orbit-stabilizer supplies raw coverage.
    result = []
    size = len(actions[0])
    level_counts = [0, 0, 0, 0]
    raw_weight = 0

    def descend(prefix, stabilizer):
        nonlocal raw_weight
        depth = len(prefix) - 1
        level_counts[depth] += 1
        if len(prefix) == 4:
            weight = len(actions) // len(stabilizer)
            raw_weight += weight
            result.append((tuple(prefix), len(stabilizer), weight))
            return
        seen = set()
        for candidate in range(size):
            if candidate in seen:
                continue
            orbit = {actions[g][candidate] for g in stabilizer}
            assert candidate == min(orbit)
            seen.update(orbit)
            smaller = [g for g in stabilizer if actions[g][candidate] == candidate]
            descend(prefix + [candidate], smaller)
        assert len(seen) == size

    descend([0], list(range(len(actions))))
    assert raw_weight == size ** 3
    return result, level_counts


def support_connected(generators):
    reached = {0}
    while True:
        new = reached | {p[i] for p in generators for i in reached}
        if new == reached:
            return len(reached) == len(generators[0])
        reached = new


def generated_group(generators):
    group = [IDENTITY]
    seen = {IDENTITY}
    for p in group:
        for q in generators:
            value = compose(p, q)
            if value not in seen:
                assert len(group) < 20160
                seen.add(value)
                group.append(value)
    return group
