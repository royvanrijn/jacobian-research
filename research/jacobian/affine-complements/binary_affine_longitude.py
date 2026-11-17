"""Small exact controls for PF125; all-degree proof is in the companion note."""
from collections import Counter


def pairing(mask, x):
    return (mask & x).bit_count() % 2


def normal_form(a, v, u, mask):
    """Coordinates l(x)=x_0, q(x)=x+x_0*v."""
    if a < 2 or not 0 < v < 2**a or v & 1 or u & 1:
        raise ValueError("invalid hyperplane data")
    if not 0 <= u < 2**a or not 0 <= mask < 2**a or pairing(mask, v):
        raise ValueError("invalid affine parameters")
    return [x ^ (u if not x & 1 else 0) ^ (v if pairing(mask, x) else 0)
            for x in range(2**a)]


def compose(p, q):
    return [p[q[x]] for x in range(len(p))]


def cycle_profile(p, domain):
    unseen = set(domain)
    counts = Counter()
    while unseen:
        x = min(unseen)
        y = p[x]
        orbit = {x}
        while y != x:
            if y in orbit or y not in unseen:
                raise ValueError("not an invariant permutation domain")
            orbit.add(y)
            y = p[y]
        unseen.difference_update(orbit)
        counts[len(orbit)] += 1
    return dict(sorted(counts.items()))


def peripheral_data(a, v, u, mask):
    p = normal_form(a, v, u, mask)
    representatives = [x for x in range(2**a) if x & 1 and x < (x ^ v)]
    critical = [(2, 1, pairing(mask, x)) for x in representatives]
    residual = cycle_profile(p, range(0, 2**a, 2))
    return critical + [(1, k, 0) for k, n in residual.items() for _ in range(n)]
