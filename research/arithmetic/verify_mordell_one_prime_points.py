#!/usr/bin/env python3
"""Exact controls for MCT17--MCT19; standard library, no rank oracle.

Checks points and their dual-isogeny images, a governing norm, all four
symbol embeddings, and the complete declared finite parameter scan.
Infinitude and the Cassels--Tate group structure remain written inputs.
"""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path

import verify_mordell_ray_pairing as ray


ROOT = Path(__file__).resolve().parent


def verify_norm(data):
    a, b = map(ray.quadratic, (data["a"], data["b"]))
    A, B, C = map(ray.quadratic, data["coefficients"])
    value = ray.add(ray.add(ray.power(A, 3), ray.mul(a, ray.power(B, 3))),
                    ray.mul(ray.power(a, 2), ray.power(C, 3)))
    value = ray.add(value, ray.scale(ray.mul(ray.mul(a, A), ray.mul(B, C)), -3))
    assert value == b, "Wrong exact governing norm"
    denominators = [c.denominator for v in (A, B, C) for c in v]
    assert all(d in (1, 3) for d in denominators)
    assert ray.norm(a) == -5**3 and ray.norm(b) == -71**3
    return a, (A, B, C)


def verify(data):
    assert data["family"] == dict(R0=355, B=5, C=71,
                                  a="71*t+1", q="71*t^2+2*t-19")
    assert data["bounded_scan"] == dict(t_min=2, t_max=200, t_step=2,
                                        t_mod_9=[0, 2])
    a5, coefficients = verify_norm(data["pair_norm"])
    base = next(c for c in json.loads(
        (ROOT / "mordell_ray_full_pairing_certificates.json").read_text())
                ["cases"] if int(c["R"]) == 355)
    assert base["basis"] == [5, 71]
    assert base["representatives"] == [["-35", "15"], ["5183", "2130"]]
    assert base["matrix"] == [[0, 2], [1, 0]]
    assert data["reference"]["entry"] == 2

    expected = []
    for t in range(2, 201, 2):
        if t % 9 in (0, 2) and ray.allowed(71*t*t + 2*t - 19):
            expected.append(t)
    assert [c["t"] for c in data["cases"]] == expected
    nonzero = []
    for case in data["cases"]:
        t, q, R, a = (case[k] for k in ("t", "q", "R", "a"))
        assert a == 71*t + 1 and q == 71*t*t + 2*t - 19
        assert R == 355*q and ray.allowed(q) and q % 72 == 53
        assert q not in (5, 71)
        X, Y = map(F, case["Q"])
        assert (X, Y) == (71*q, a*71*q)
        assert Y*Y == X**3 + 54*R*R, "Wrong point on E prime"
        x, y = map(F, case["point_on_E"])
        assert x == (X**3 + 216*R*R)/(9*X*X)
        assert y == Y*(X**3 - 432*R*R)/(27*X**3)
        assert y*y == x**3 - 2*R*R, "Wrong point on E"
        assert x.denominator == 9 and y.denominator == 27
        values = {ray.pair_symbol(q, a5, coefficients, -5,
                                 negative_s=ns, inverse_zeta=iz)
                  for ns in (False, True) for iz in (False, True)}
        assert values == {case["kappa"]}
        assert case["entry"] == (2 + case["kappa"]) % 3
        if case["entry"]:
            nonzero.append(R)
    return nonzero


def main():
    data = json.loads((ROOT / "mordell_one_prime_point_certificates.json").read_text())
    nonzero = verify(data)
    # These finite residue checks have exactly the universal local scopes
    # stated in the written proof, without an inference to prime infinitude.
    units = {b*(3*a*a - 2*b*b) % 9 for a in range(9) for b in range(9)
             if b % 3}
    assert units == {1, 2, 7, 8}
    for p in (11, 17):
        for R in range(1, p):
            for sign in (-2, 54):
                count = 1 + sum((y*y - x*x*x - sign*R*R) % p == 0
                                for x in range(p) for y in range(p))
                assert count == p + 1
    bad_norm = deepcopy(data)
    bad_norm["pair_norm"]["coefficients"][0][0] = str(
        F(bad_norm["pair_norm"]["coefficients"][0][0]) + 1)
    bad_point = deepcopy(data)
    bad_point["cases"][1]["point_on_E"][1] = str(
        F(bad_point["cases"][1]["point_on_E"][1]) + 1)
    bad_symbol = deepcopy(data)
    bad_symbol["cases"][1]["entry"] = 0
    for bad in (bad_norm, bad_point, bad_symbol):
        try:
            verify(bad)
        except AssertionError:
            pass
        else:
            raise AssertionError("A corrupted certificate was accepted")
    print("PASS: exact norm; four embeddings per symbol; five complete bounded "
          "prime controls; rational points and dual-isogeny images; local "
          "residue counts; three corruptions rejected.")
    print("Nonzero fixed blocks certify rank one and Sha[3^infinity]=(Z/3)^2 at R:",
          nonzero)
    print("No prime-value infinitude, CM converse, or independent review is certified.")


if __name__ == "__main__":
    main()
