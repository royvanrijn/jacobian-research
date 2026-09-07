#!/usr/bin/env python3
"""Sage-free replay of twelve exact simultaneous lifts and height obstructions."""
import argparse,json
from fractions import Fraction as F
from math import gcd,isqrt
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
INPUT=ART/'native_rank3_carrier_images_v1.json';OUT=ART/'native_rank3_carrier_images_replay_v1.json'
def polynomial(coefficients,t):return sum(F(a)*t**i for i,a in enumerate(coefficients))
def expected():
    d=cert.read(INPUT);g=cert.read(ART/'rank_jump_minimal_native_block_carrier_v1.json')['geometry'];ci=cert.read(ART/'rank_jump_minimal_native_block_carrier_inputs_v1.json');inp=cert.read(ART/'rank_jump_native_triple_intersection_inputs_v1.json');family=next(f for f in cert.read(spec.ATLAS)['families'] if f['family']=='08234');covers={c['label']:c for c in inp['covers']};seen=set();rows=[]
    for name,h in d['sources'].items():
        if cert.hashed(ROOT/name)!=h:raise ArithmeticError('source input differs')
    for r in d['rows']:
        if r['status']!='EXACT_TWO_POINT_IMAGE':raise ArithmeticError('complete twelve finite images required')
        z,w=map(F,r['quartic_point']);N=polynomial(list(reversed(g['conic_parameter_numerator'])),z);D=polynomial(list(reversed(g['conic_parameter_denominator'])),z);U=polynomial(list(reversed(g['conic_root_numerator'])),z)
        if not D or w*w!=polynomial(list(reversed(g['quartic_coefficients'])),z):raise ArithmeticError('quartic point differs')
        t=N/D;s=-26*t-50;roots=[U/D,w/(F(g['quartic_square_scale'])*D)];height=max(abs(s.numerator),s.denominator)
        if str(t)!=r['published_parameter'] or str(s)!=r['compact_parameter'] or list(map(str,roots))!=r['square_roots'] or r['parameter_height']!=str(height) or r['parameter_bits']!=height.bit_length() or r['beyond_compact_region']!=(height>4096) or r['is_known_origin']!=(t==F(ci['anchor']['t'])):raise ArithmeticError('parameter projection differs')
        model=[F(0)]*3+[sum(F(a)*s.numerator**k*s.denominator**(weight-k) for k,a in enumerate(family[name])) for name,weight in [('A_coefficients_low_to_high',8),('B_coefficients_low_to_high',12)]]
        if list(map(str,model))!=r['curve'] or any(a.denominator!=1 for a in model):raise ArithmeticError('independent homogeneous compact model differs')
        pub=[F(0)]*3+[polynomial(inp['A'],t),polynomial(inp['B'],t)];native=[]
        for row,form,root in zip(ci['covers'],g['forms'],roots):
            if root*root!=polynomial(form,t):raise ArithmeticError('simultaneous square condition differs')
            c=covers[row['label']];v=root*F(row['removed_rational_square_root']);lift=c['lifted_section'];P=tuple(polynomial(lift[k+'0_coefficients'],t)+v*polynomial(lift[k+'1_coefficients'],t) for k in ('x','y'))
            if not cert.is_on_weierstrass_curve(pub,P):raise ArithmeticError('native point membership failed')
            native.append(P)
        a,b=cert.weierstrass_invariants(pub),cert.weierstrass_invariants(model);u=cert.square_root(a['c6']*b['c4']/(a['c4']*b['c6']))
        if u is None or not u or a['c4']!=u**4*b['c4'] or a['c6']!=u**6*b['c6']:raise ArithmeticError('exact rational model scale differs')
        extra=[tuple(map(F,P)) for P in r['supplied_points']]
        if not any(extra==[(x/u**2,y/(sign*u)**3) for x,y in native] for sign in (1,-1)) or any(not cert.is_on_weierstrass_curve(model,P) for P in extra):raise ArithmeticError('exact native point transport differs')
        c4,c6=int(b['c4']),int(b['c6']);G=gcd(abs(c4),abs(c6));num=c6*c6;den=G**3;M=2**400-1;excluded=num>den*(1224*M+521)**2
        q=(num+den-1)//den;lower=isqrt(q);lower+=int(lower*lower<q);bits=max(0,(lower-521+1223)//1224).bit_length()
        if not excluded or r['all_normalized_400_bit_integral_models_excluded']!=excluded or r['normalized_coefficient_bits_lower_bound']!=bits or r['c4']!=str(c4) or r['c6']!=str(c6) or r['invariant_gcd']!=str(G):raise ArithmeticError('exact normalized height inequality differs')
        j=b['c4']**3/b['discriminant']
        if j in seen:raise ArithmeticError('same j needs a separate equation comparison')
        seen.add(j);rows.append({'word':r['word'],'compact_parameter':str(s),'normalized_coefficient_bits_lower_bound':bits,'beyond_compact_region':height>4096,'known_origin':r['is_known_origin']})
    if len(rows)!=12:raise ArithmeticError('fixed twelve words required')
    paths=[Path(__file__).resolve(),INPUT,Path(cert.__file__),spec.ATLAS,ART/'rank_jump_minimal_native_block_carrier_v1.json',ART/'rank_jump_minimal_native_block_carrier_inputs_v1.json',ART/'rank_jump_native_triple_intersection_inputs_v1.json']
    return {'schema':'elliptic-curves.native-rank3-carrier-images-replay.v1','status':'PASS','sources':{str(q.relative_to(ROOT)):cert.hashed(q) for q in paths},'rows':rows,'distinct_j_invariants':len(seen),'minimum_normalized_coefficient_bits_lower_bound':min(r['normalized_coefficient_bits_lower_bound'] for r in rows),'boundaries':'Independent rational arithmetic verifies quartic images, both square roots, primitive compact parameters, homogeneous models, both native point identities and exact model transports. All twelve j-invariants differ and every image exceeds normalized400-bit integral size. No independence/rank claim for these specializations, all-carrier height bound, auxiliary rank replay or universal novelty.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('independent image replay differs')
    else:
        if OUT.exists():raise FileExistsError('preserve independent image replay')
        checkpoint(OUT,d)
    print('INDEPENDENT NATIVE CARRIER IMAGE REPLAY',d['distinct_j_invariants'],d['minimum_normalized_coefficient_bits_lower_bound'],'PASS',flush=True)
