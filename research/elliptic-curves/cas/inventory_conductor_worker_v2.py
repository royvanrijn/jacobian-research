#!/usr/bin/env python3
"""Verify a denominator-clearing Q-isomorphism, then use the frozen V1 worker."""
import argparse
from pathlib import Path
from sage.all import QQ, ZZ, EllipticCurve
from v3_warm_support import read, require
from inventory_conductor_worker import run as integral_run


def verify_transport(source):
    original = list(map(QQ, source['original_curve']))
    integral = list(map(QQ, source['curve']))
    scale = ZZ(source['integral_model_scale'])
    require(scale > 0 and len(original) == len(integral) == 5, 'invalid model transport')
    require(integral == [a*scale**w for a,w in zip(original, (1,2,3,4,6))],
            'model transport coefficients differ')
    require(all(a.denominator() == 1 for a in integral), 'transport is not integral')
    old, new = EllipticCurve(QQ, original), EllipticCurve(QQ, integral)
    require(new.discriminant() == old.discriminant()*scale**12, 'transport discriminant differs')
    require(old.is_isomorphic(new), 'rational model isomorphism failed')
    # Explicit coordinate map: (x,y) -> (scale^2*x, scale^3*y).
    # This changes no rational group or conductor; no minimality is claimed.


def run(packet, output, check=False):
    verify_transport(read(packet))
    integral_run(packet, output, check)


if __name__ == '__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--packet',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--check',action='store_true')
    a=p.parse_args();run(a.packet,a.output,a.check)
