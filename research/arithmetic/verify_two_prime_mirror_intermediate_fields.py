"""Exact affine group, subgroup mass and conductor fixtures for AFM1–AFM5."""

from collections import Counter
from itertools import product

from sage.all import CyclotomicField, GF, ZZ, prod, valuation


def element_order_exponent(element, ell, exponents):
    orders = (exponent - min(valuation(ZZ(coordinate), ell), exponent)
              if coordinate else 0 for coordinate, exponent in zip(element, exponents))
    return max(orders, default=0)


def group_elements(ell, exponents):
    return list(product(*(range(ell**exponent) for exponent in exponents)))


def phi_power(ell, exponent):
    return (ell-1)*ell**(exponent-1) if exponent else 1


def automorphism_count(ell, a, b):
    if a == b:
        return ell**(4*(a-1))*(ell**2-1)*(ell**2-ell)
    return phi_power(ell, a)*phi_power(ell, b)*ell**(2*b)


def check_profiles_and_conductors(ell, level):
    for power_sum in range(1, 2*level+1):
        exponents_by_type = []
        for a in range(1, level+1):
            b = power_sum-a
            if not 0 <= b <= a:
                continue
            elements = group_elements(ell, (a, b))
            orders = {x: element_order_exponent(x, ell, (a, b)) for x in elements}
            profiles = Counter()
            for x, y in product(elements, repeat=2):
                if b:
                    generates = (x[0]*y[1]-x[1]*y[0]) % ell != 0
                else:
                    generates = x[0] % ell != 0 or y[0] % ell != 0
                if generates:
                    profiles[orders[x], orders[y]] += 1
            expected = {}
            for i, j in product(range(a+1), repeat=2):
                c = min(i, j)
                if max(i, j) == a and c >= b:
                    expected[i, j] = 1 if c == b else phi_power(ell, c-b)
            automorphisms = automorphism_count(ell, a, b)
            assert profiles == {key: count*automorphisms for key, count in expected.items()}

            counts = Counter(orders.values())
            root_orbit_sum = 0
            for j in range(1, a+1):
                orbit_count, remainder = divmod(counts[j], phi_power(ell, j))
                assert remainder == 0
                root_disc = CyclotomicField(ell**j).discriminant()
                root_orbit_sum += orbit_count*valuation(ZZ(root_disc), ell)
            assert root_orbit_sum == sum(orders.values()) - (len(elements)-1)//(ell-1)
            exponents_by_type.append((a, root_orbit_sum))
        assert exponents_by_type == sorted(set(exponents_by_type), key=lambda pair: pair[0])
        exponents = [exponent for _, exponent in exponents_by_type]
        assert exponents == sorted(set(exponents))


def check_fixed_resolvent_mass():
    ell, invariants = 3, (2, 2, 1)
    elements = group_elements(ell, invariants)
    orders = {x: element_order_exponent(x, ell, invariants) for x in elements}
    sizes = [prod(ell**min(j, power) for power in invariants) for j in range(3)]

    def multiply(x, integer):
        return tuple(integer*coordinate % ell**power
                     for coordinate, power in zip(x, invariants))

    for a, b in ((1, 1), (2, 1), (2, 2)):
        injections = 0
        for x in elements:
            if orders[x] != a:
                continue
            line = {multiply(x, scalar*ell**(a-1)) for scalar in range(ell)}
            for y in elements:
                if orders[y] == b and multiply(y, ell**(b-1)) not in line:
                    injections += 1
        assert injections == (sizes[a]-sizes[a-1])*(sizes[b]-ell*sizes[b-1])
        assert injections % automorphism_count(ell, a, b) == 0


def check_prime_level_local_profiles():
    for ell in (3, 5, 7, 11, 13):
        # At an ell-ramified quadratic prime the inertia projection sends
        # a to its Legendre symbol times a.
        scalars = GF(ell)
        projected = {a if a.is_square() else -a for a in scalars if a}
        inertia_order = ell-1 if ell % 4 == 1 else (ell-1)//2
        assert len(projected) == inertia_order
        assert (ell**2-1) - (ell**2-1)//inertia_order == (
            (ell+1)*(ell-2 if ell % 4 == 1 else ell-3)
        )


if __name__ == "__main__":
    check_profiles_and_conductors(3, 2)
    check_profiles_and_conductors(5, 1)
    check_fixed_resolvent_mass()
    check_prime_level_local_profiles()
    print("AFM1–AFM5 affine order profiles, masses and conductor fixtures verified.")
