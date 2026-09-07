#!/usr/bin/env python3
"""Exact height inequalities for all finite images, including earlier size censoring."""
import argparse,json
from fractions import Fraction as F
from math import gcd,isqrt
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';INPUT=ART/'native_rank3_fullspan_images_v1.json';OUT=ART/'native_rank3_fullspan_heights_v1.json'
def poly(coeff,t):return sum(F(a)*t**i for i,a in enumerate(coeff))
def expected():
    d=cert.read(INPUT);source=cert.read(ART/'rank_jump_minimal_native_block_carrier_v1.json');g=source['geometry'];ci=cert.read(ART/'rank_jump_minimal_native_block_carrier_inputs_v1.json');inp=cert.read(ART/'rank_jump_native_triple_intersection_inputs_v1.json');family=next(f for f in cert.read(spec.ATLAS)['families'] if f['family']=='08234');covers={c['label']:c for c in inp['covers']};old={r['compact_parameter'] for r in cert.read(ART/'native_rank3_carrier_images_v1.json')['rows']};rows=[];seen={};exceptions=[]
    if d['status']!='PASS_FIXED_AUDIT' or len(d['rows'])!=125 or any(cert.hashed(ROOT/n)!=h for n,h in d['sources'].items()):raise ArithmeticError('complete fixed125 image audit required')
    for r in d['rows']:
        if r['status']=='BIRATIONAL_EXCEPTION':exceptions.append(r['word']);continue
        z,w=map(F,r['quartic_point']);D=poly(list(reversed(g['conic_parameter_denominator'])),z);N=poly(list(reversed(g['conic_parameter_numerator'])),z);U=poly(list(reversed(g['conic_root_numerator'])),z)
        if not D or w*w!=poly(list(reversed(g['quartic_coefficients'])),z):raise ArithmeticError('quartic image identity differs')
        t=N/D;s=-26*t-50;roots=[U/D,w/(F(g['quartic_square_scale'])*D)];height=max(abs(s.numerator),s.denominator)
        if str(t)!=r['published_parameter'] or str(s)!=r['compact_parameter'] or list(map(str,roots))!=r['square_roots'] or r['parameter_bits']!=height.bit_length():raise ArithmeticError('finite parameter projection differs')
        if r['status']=='PARAMETER_HEIGHT_CENSORED' and height.bit_length()<=512:raise ArithmeticError('censoring threshold differs')
        model=[F(0)]*3+[sum(F(a)*s.numerator**k*s.denominator**(weight-k) for k,a in enumerate(family[name])) for name,weight in [('A_coefficients_low_to_high',8),('B_coefficients_low_to_high',12)]]
        if any(a.denominator!=1 for a in model):raise ArithmeticError('integral compact model required')
        if 'curve' in r and tuple(map(F,r['curve']))!=tuple(model):raise ArithmeticError('homogeneous specialization differs')
        pub=[F(0)]*3+[poly(inp['A'],t),poly(inp['B'],t)];native=[]
        for form,root,c in zip(g['forms'],roots,ci['covers']):
            if root*root!=poly(form,t):raise ArithmeticError('simultaneous rational square differs')
            lift=covers[c['label']]['lifted_section'];v=root*F(c['removed_rational_square_root']);P=tuple(poly(lift[k+'0_coefficients'],t)+v*poly(lift[k+'1_coefficients'],t) for k in ('x','y'))
            if not cert.is_on_weierstrass_curve(pub,P):raise ArithmeticError('native point identity failed')
            native.append(P)
        a,b=cert.weierstrass_invariants(pub),cert.weierstrass_invariants(model);u=cert.square_root(a['c6']*b['c4']/(a['c4']*b['c6']))
        if u is None or not u or a['c4']!=u**4*b['c4'] or a['c6']!=u**6*b['c6']:raise ArithmeticError('exact model scale differs')
        moved=[(x/u**2,y/u**3) for x,y in native]
        if any(not cert.is_on_weierstrass_curve(model,P) for P in moved):raise ArithmeticError('compact native point identity failed')
        if 'supplied_points' in r and not any([tuple(map(F,P)) for P in r['supplied_points']]==[(x,sign*y) for x,y in moved] for sign in (1,-1)):raise ArithmeticError('recorded native transport differs')
        c4,c6=int(b['c4']),int(b['c6']);G=gcd(abs(c4),abs(c6));num=c6*c6;den=G**3;M=2**400-1;excluded=num>den*(1224*M+521)**2
        ceil=(num+den-1)//den;lower=isqrt(ceil);lower+=int(lower*lower<ceil);bits=max(0,(lower-521+1223)//1224).bit_length();origin=t==F(ci['anchor']['t']);j=b['c4']**3/b['discriminant']
        if not origin and not excluded:raise ArithmeticError('a new affordable image requires explicit review')
        if 'all_normalized_400_bit_integral_models_excluded' in r and r['all_normalized_400_bit_integral_models_excluded']!=excluded:raise ArithmeticError('earlier height gate differs')
        seen.setdefault(j,set()).add(str(s));rows.append({'word':r['word'],'compact_parameter':str(s),'known_origin':origin,'earlier_twelve_image':str(s) in old,'previously_size_censored':r['status']=='PARAMETER_HEIGHT_CENSORED','parameter_bits':height.bit_length(),'normalized_coefficient_bits_lower_bound':bits,'normalized_400_bit_models_excluded':excluded,'c4_hex':hex(c4),'c6_hex':hex(c6),'invariant_gcd_hex':hex(G),'both_native_points_verified':True})
    parameters={r['compact_parameter'] for r in rows};nonanchor={r['compact_parameter'] for r in rows if not r['known_origin']}
    if len(rows)!=124 or len(exceptions)!=1:raise ArithmeticError('fixed finite/exception census differs')
    paths=[Path(__file__).resolve(),INPUT,Path(cert.__file__),spec.ATLAS,ART/'rank_jump_minimal_native_block_carrier_v1.json',ART/'rank_jump_minimal_native_block_carrier_inputs_v1.json',ART/'rank_jump_native_triple_intersection_inputs_v1.json',ART/'native_rank3_carrier_images_v1.json']
    return {'schema':'elliptic-curves.native-rank3-fullspan-heights.v1','status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'rows':rows,'source_birational_exceptions':exceptions,'finite_images':len(rows),'distinct_parameters':len(parameters),'distinct_j_invariants':len(seen),'nonanchor_distinct_parameters_excluded':len(nonanchor),'additional_parameters_beyond_first_twelve_and_origin':len(nonanchor-old),'minimum_nonanchor_normalized_bits_lower_bound':min(r['normalized_coefficient_bits_lower_bound'] for r in rows if not r['known_origin']),'claim_boundary':'Independent exact arithmetic verifies all124 finite parameter images, including45 earlier size-censored cases, both native point identities and rational model transports. Hexadecimal invariants preserve large exact integers without decimal conversion limits. All nonanchor finite images exceed the400-bit normalized integral model gate by c6^2>gcd(c4,c6)^3*(1224*(2^400-1)+521)^2. The single birational exception is retained from the source. This fixed cube is not a full-carrier height theorem, saturation, global point enumeration, specialized rank proof or universal novelty claim.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('fullspan exact height replay differs')
    else:
        if OUT.exists():raise FileExistsError('preserve fullspan height proof')
        checkpoint(OUT,d)
    print('FULLSPAN EXACT HEIGHTS',d['finite_images'],d['distinct_parameters'],d['nonanchor_distinct_parameters_excluded'],d['additional_parameters_beyond_first_twelve_and_origin'],'PASS',flush=True)
