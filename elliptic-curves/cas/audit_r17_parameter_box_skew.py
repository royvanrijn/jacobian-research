#!/usr/bin/env python3
"""Exact finite optimization of triangle coefficient bounds for equal-area parameter boxes."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';ATLAS=ART/'compact_six_r17_atlas_v1.json';OUT=ART/'r17_parameter_box_skew_v1.json'

def expected():
    rows=[]
    for f in cert.read(ATLAS)['families']:
        A=list(map(F,f['A_coefficients_low_to_high']));B=list(map(F,f['B_coefficients_low_to_high']))
        if len(A)!=9 or len(B)!=13:raise ArithmeticError('fixed degrees8 and12 required')
        candidates=[]
        for k in range(-4,5):
            ua=sum(abs(a)*F(2)**(k*(2*i-8)) for i,a in enumerate(A));ub=sum(abs(b)*F(2)**(k*(2*i-12)) for i,b in enumerate(B));w=max(4*ua**3,27*ub**2)
            candidates.append({'k':k,'numerator_bound':int(32768*F(2)**k),'denominator_bound':int(32768/F(2)**k),'A_triangle_constant':str(ua),'B_triangle_constant':str(ub),'weighted_triangle_bound':str(w)})
        chosen=min(candidates,key=lambda r:(F(r['weighted_triangle_bound']),abs(r['k']),r['k']));old=next(r for r in candidates if r['k']==0)
        if any(r['numerator_bound']*r['denominator_bound']!=32768**2 for r in candidates):raise ArithmeticError('equal rectangle areas required')
        rows.append({'family':f['family'],'candidates':candidates,'selected_k':chosen['k'],'selected_numerator_bound':chosen['numerator_bound'],'selected_denominator_bound':chosen['denominator_bound'],'weighted_bound_improvement':str(F(old['weighted_triangle_bound'])/F(chosen['weighted_triangle_bound']))})
    return {'schema':'elliptic-curves.r17-parameter-box-skew.v1','status':'PASS_EXACT_FINITE_BOUND_OPTIMIZATION','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),ATLAS)},'height':32768,'k_values':list(range(-4,5)),'rows':rows,'identity':'For A_h(n,d)=sum A_i*n^i*d^(8-i) and B_h(n,d)=sum B_i*n^i*d^(12-i), a rectangle |n|<=H*2^k,|d|<=H*2^-k has |A_h|<=H^8*U_A(k) and |B_h|<=H^12*U_B(k), where U_A=sum |A_i|*2^(k*(2i-8)) and U_B=sum |B_i|*2^(k*(2i-12)). Hence max(4|A_h|^3,27|B_h|^2)<=H^24*max(4U_A^3,27U_B^2). All values and the minimum over the nine declared k are exact rational arithmetic.','claim_boundary':'A conservative coefficient-bound comparison of equal-area rectangles, not minimal discriminant, conductor, height-density, candidate-rank or point-visibility calibration. Equal geometric area does not imply equal primitive counts in fixed denominator slices. A new rectangle is a new parameter population and requires its own frozen scanner and proof gates. No parameter scan or point search is performed here.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=d:raise ArithmeticError('exact parameter skew bound differs')
    else:
        if OUT.exists():raise FileExistsError('preserve parameter-bound audit')
        checkpoint(OUT,d)
    for r in d['rows']:print('EXACT SKEW BOUND',r['family'],r['selected_k'],r['selected_numerator_bound'],r['selected_denominator_bound'],'improvement bits approx',F(r['weighted_bound_improvement']).numerator.bit_length()-F(r['weighted_bound_improvement']).denominator.bit_length(),flush=True)
