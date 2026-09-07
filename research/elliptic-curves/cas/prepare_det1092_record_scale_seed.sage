#!/usr/bin/env sage-python
"""Exact homogeneous specialization of the seventeen generic sections."""
import sys,argparse
from pathlib import Path
from fractions import Fraction as F
from sage.all import QQ,PolynomialRing,EllipticCurve
CAS=Path(__file__).resolve().parent;sys.path.insert(0,str(CAS))
import det1092_record_scale_points as trial
from research_runtime.search_state import raw_state
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache
from research_runtime.memory_store import MemoryFactStore
def main(index):
    row=trial.campaign()['rows'][index];parent=trial.cert.read(trial.selection.SOURCE)
    R=PolynomialRing(QQ,'t');t=QQ(row['parameter']);q=t.denominator()
    value=lambda v:R(v['numerator'])(t)/R(v['denominator'])(t)
    a=[value(v) for v in parent['a_invariants']];assert a[:3]==[0,0,0]
    model=[0,0,0,a[3]*q**8,a[4]*q**12];assert list(map(str,model))==row['model']
    E=EllipticCurve(QQ,model)
    points=[E([value(P[0])*q**4,value(P[1])*q**6]) for P in parent['basis_weierstrass_coordinates']]
    assert len(points)==17 and all(not P.is_zero() for P in points)
    m=tuple(F(str(v)) for v in model);pts=tuple(tuple(F(str(v)) for v in P.xy()) for P in points)
    state=raw_state(m,pts,cache=QuotientOnlyReductionCache(MemoryFactStore()),prime_bound=1000)
    assert tuple(tuple(map(F,P)) for P in state.basis)==pts
    proof=trial.checked_rank(m,pts,state.reductions.primes,state.no_two_torsion_prime)
    xy=[list(map(str,P)) for P in pts]
    out=trial.BATCH/row['id']/'initial-seed.json';assert not out.exists()
    trial.checkpoint(out,dict(family='det1092-reduced',parameter=row['parameter'],curve=row['model'],points=xy,generic_points=xy,rank_certificate=proof,
        source_parent_sha256=trial.cert.hashed(trial.selection.SOURCE),source_selected_model=row['model'],homogeneous_scale=str(q**2)))
    print('EXACT17 GENERIC SECTIONS',row['id'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--index',type=int,required=True);a=ap.parse_args();main(a.index)
