#!/usr/bin/env python3
"""Exact terminal-fourteen rays at small outer mean."""

from fractions import Fraction as Q
from sys import argv

from verify_multivariate_terminal_fourteen_retained_profiles import configure_shared_model

SCALE = 400
CAP = Q(1, 8)
RAYS = {
    0: (972935, 876387, 765171, 591197, 394244, 177302, -60268, -277210, -500974, -674948, -833644, -930192),
    1: (972657, 875256, 763044, 588047, 390016, 172668, -65119, -281566, -504528, -677490, -835060, -930821),
    2: (972351, 874039, 760773, 584750, 385664, 167977, -69941, -285838, -507963, -679922, -836406, -931413),
    3: (972019, 872751, 758385, 581345, 381232, 163262, -74720, -290031, -511301, -682269, -837698, -931978),
    4: (971660, 871394, 755888, 577847, 376738, 158538, -79454, -294152, -514554, -684544, -838944, -932522),
    5: (971271, 869972, 753292, 574269, 372196, 153812, -84141, -298206, -517731, -686757, -840153, -933048),
    6: (970852, 868487, 750603, 570619, 367616, 149089, -88782, -302197, -520840, -688914, -841327, -933557),
    7: (970400, 866940, 747826, 566904, 363005, 144375, -93377, -306127, -523885, -691021, -842470, -934052),
    8: (969914, 865332, 744966, 563131, 358370, 139671, -97927, -310001, -526872, -693082, -843585, -934535),
}
RECTANGLES = {
    0: (123, 182),
    1: (207, 434),
    2: (340, 724),
    3: (370, 1294),
    4: (628, 1352),
    5: (704, 2568),
    6: (708, 4097),
    7: (796, 16396),
}


def main():
    model = configure_shared_model()
    model.verify_setup(RAYS, {}, CAP, True, knot_scale=SCALE)
    chosen = tuple(map(int, argv[1:])) if len(argv) > 1 else tuple(RECTANGLES)
    for left in chosen:
        accepted = tuple(model.verify_small_band(
            order, CAP, RAYS, left, left+1, knot_scale=SCALE
        ) for order in model.ORDERS)
        expected = RECTANGLES[left]
        model.require(accepted == expected, ("small outer mean replay", left, accepted, expected))
        print("Verified", Q(left, SCALE), Q(left+1, SCALE), accepted, flush=True)


if __name__ == "__main__":
    main()
