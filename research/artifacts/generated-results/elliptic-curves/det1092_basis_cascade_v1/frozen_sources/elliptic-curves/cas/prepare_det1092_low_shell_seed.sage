#!/usr/bin/env sage-python
"""Specialize the fixed MW17 determinant-1092 parent for one low-shell row."""
import argparse
import sys
from fractions import Fraction as F
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ

CAS=Path(__file__).resolve().parent; sys.path.insert(0,str(CAS))
import det1092_low_shell_cascade as control
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache
from research_runtime.search_state import raw_state


def main(index):
    row=control.campaign()["rows"][index]; parent=control.cert.read(control.selection.SOURCE)
    ring=PolynomialRing(QQ,"t"); t=QQ(row["parameter"]); denominator=t.denominator()
    value=lambda record:ring(record["numerator"])(t)/ring(record["denominator"])(t)
    a=[value(record) for record in parent["a_invariants"]]
    model=[0,0,0,a[3]*denominator**8,a[4]*denominator**12]
    if list(map(str,model)) != row["model"]: raise ArithmeticError("frozen low-shell specialization model differs")
    curve=EllipticCurve(QQ,model)
    points=[curve([value(point[0])*denominator**4,value(point[1])*denominator**6]) for point in parent["basis_weierstrass_coordinates"]]
    if len(points)!=17 or any(point.is_zero() for point in points): raise ArithmeticError("generic MW17 specialization malformed")
    model_q=tuple(F(str(value)) for value in model); points_q=tuple(tuple(F(str(value)) for value in point.xy()) for point in points)
    state=raw_state(model_q,points_q,cache=QuotientOnlyReductionCache(MemoryFactStore()),prime_bound=1000)
    if tuple(tuple(map(F,point)) for point in state.basis) != points_q: raise ArithmeticError("seed changed before exact certificate")
    proof=control.checked_rank(model_q,points_q,state.reductions.primes,state.no_two_torsion_prime)
    out=control.BATCH/row["id"]/"initial-seed.json"
    if out.exists(): raise FileExistsError("preserve frozen low-shell initial seed")
    control.checkpoint(out,{"family":"det1092-reduced","parameter":row["parameter"],"curve":row["model"],"points":[list(map(str,point)) for point in points_q],"generic_points":[list(map(str,point)) for point in points_q],"rank_certificate":proof,"source_parent_sha256":control.cert.hashed(control.selection.SOURCE),"source_selected_model":row["model"],"homogeneous_scale":str(denominator**2)})
    print("FROZEN LOW-SHELL MW17 SEED",row["id"],flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--index",type=int,required=True);args=parser.parse_args();main(args.index)
