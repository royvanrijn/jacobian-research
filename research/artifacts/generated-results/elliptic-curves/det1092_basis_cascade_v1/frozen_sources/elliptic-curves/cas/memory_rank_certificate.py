"""Exact stand-alone mod2 certificates with an explicit ephemeral finite cache.

The returned points/signatures/torsion proof are replay inputs. No persistent
number-field fact or worker state is replaced. This avoids per-point disk
fact publication when independently reconstructing a finite certificate.
"""
from dataclasses import asdict
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from mod2_reduction_independence import combined_mod2_rank
from elliptic_candidate_record import is_on_weierstrass_curve
from certify_compact_r17_candidates import short_curve_has_no_rational_2_torsion_modular_certificate

def checked_rank(model,points,primes,torsion_prime):
    if not points or any(not is_on_weierstrass_curve(model,p) for p in points):raise ArithmeticError('point membership failed')
    cache=ReductionCache(MemoryFactStore());sigs=tuple(cache.signature(model,points,p) for p in primes)
    if combined_mod2_rank(sigs,len(points))!=len(points):raise ArithmeticError('independent finite columns not established')
    if not short_curve_has_no_rational_2_torsion_modular_certificate(model,torsion_prime):raise ArithmeticError('2-torsion witness failed')
    return {'rank_lower_bound':len(points),'no_rational_2_torsion_prime':torsion_prime,'signatures':[asdict(s) for s in sigs],'argument':'Every integral relation is divisible by 2; E(Q)[2]=0 permits infinite descent. Hence all listed points are independent.'}
