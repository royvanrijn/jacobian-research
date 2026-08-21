from sage.all import *
from itertools import product
import argparse


parser = argparse.ArgumentParser(description="Search Möbius sections on the exact MW3 quartic modulo p.")
parser.add_argument("--p", type=int, default=23)
args = parser.parse_args()

k = GF(args.p)
K = FunctionField(k, "U")
U = K.gen()
RT = PolynomialRing(K, "t")
t = RT.gen()


def chord_discriminant(z):
    return (
        t**5 + (k(21)/50*z - k(347)/200)*t**4
        + (-k(651)/625*z + k(53)/625)*t**3
        + (-k(3969)/62500*z**2 + k(13251)/31250*z + k(30471)/62500)*t**2
        + (k(3528)/390625*z**2 + k(7056)/390625*z + k(3528)/390625)*t
        + k(194481)/78125000*z**4 - k(583443)/39062500*z**2
        - k(194481)/9765625*z - k(583443)/78125000
    )


z = t * (t - 1) * U - 1 - k(25)/21 * t
completed = chord_discriminant(z)
quartic, remainder = completed.quo_rem(t**2 * (t - 1)**2)
assert remainder == 0 and quartic.degree() == 4 and quartic[0] == 0


def normalized_projective_tuples():
    for pivot in range(4):
        prefix = [k.zero()] * pivot + [k.one()]
        for suffix in product(k, repeat=3 - pivot):
            yield tuple(prefix + list(suffix))


hits = []
tested = 0
for a, b, c, d in normalized_projective_tuples():
    if a * d - b * c == 0:
        continue
    tested += 1
    denominator = c * U + d
    if denominator == 0:
        continue
    t_value = (a * U + b) / denominator
    value = K(quartic(t_value))
    if value.is_square():
        y_value = value.sqrt()
        hits.append((a, b, c, d, t_value, y_value))
        print(
            f"MW3MOBIUS|matrix={int(a)},{int(b)},{int(c)},{int(d)}"
            f"|t={t_value}|y={y_value}",
            flush=True,
        )

print(f"MW3MOBIUS|p={args.p}|tested={tested}|hits={len(hits)}|status=PASS", flush=True)
