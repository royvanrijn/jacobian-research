#!/usr/bin/env sage-python
"""A fixed125-word grid using three certified independent auxiliary points."""
import argparse,json,sys
from itertools import product
from math import gcd,isqrt
from pathlib import Path
from importlib.machinery import SourceFileLoader
from sage.all import QQ,PolynomialRing,EllipticCurve
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
from research_runtime.store import checkpoint
from research_runtime.search_state import raw_state
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from memory_rank_certificate import checked_rank
ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/native-rank3-fullspan-images-v1';OUT=ART/'native_rank3_fullspan_images_v1.json'
HELPER=ROOT/'elkies-k3/scripts/prepare_r17_norm12_11952_alternate_v4_rank_one_bases.sage'
NAMES=['rank_jump_minimal_native_block_carrier_v1.json','rank_jump_minimal_native_block_carrier_inputs_v1.json','rank_jump_native_triple_intersection_inputs_v1.json','rank_jump_native_triple_carrier_verification_v1.json','native_rank3_carrier_marked_point_replay_v1.json']
WORDS=[list(w) for w in product(range(-2,3),repeat=3)]
def sources():
    paths=[Path(__file__).resolve(),HELPER,spec.ATLAS,*(ART/n for n in NAMES),Path(spec.__file__),Path(cert.__file__),CAS/'memory_rank_certificate.py',CAS/'research_runtime/search_state.py',CAS/'research_runtime/finite_reduction.py',CAS/'research_runtime/store.py']
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve native carrier image protocol')
    v=cert.read(ART/NAMES[3]);explicit=cert.read(ART/NAMES[4])
    if explicit['status']!='PASS' or explicit['rank_lower_bound']!=3 or len(explicit['three_independent_Jacobian_points'])!=3:raise ArithmeticError('three explicit independent auxiliary points required')
    if v['status']!='PASS' or v['auxiliary_Jacobian_exact_rank']!=3 or v['specialized_pair_quotient_rank']!=2:raise ArithmeticError('marked carrier proof gate required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.native-rank3-fullspan-images.v1','sources':sources(),'words':WORDS,'maximum_parameter_bits':512,'model_bits_gate':400,'compact_inner_height':4096,'seconds_per_stage':300,'rss_bytes':2147483648,'gate':'The earlier twelve-word construction sampled only two returned auxiliary points despite the proven Jacobian rank3. The opposite marked quartic point now supplies a third explicit direction, with independent exact mod3 and mod5 rank3 proofs. Explore one fixed coefficient cube in this full-rank subgroup to test the previously omitted auxiliary direction.','scope':'Exactly125 triples in[-2,2]^3, including the identity, in the three certified independent auxiliary points. All words freeze before any parameter image. Exact pointed-quartic inversion, simultaneous roots and native points, compact08234 transports and normalized400-bit height inequalities. Parameters over512 bits are explicitly censored before specialization. Known origins and duplicates remain in the roster without refill. Finite19-point certificates only for new outer-compact parameters with displayed models at most400 bits. No new auxiliary or original point search, saturation, change of basis, full-carrier height completeness, rank-density claim or automatic expansion.'})
def expected():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['words']!=WORDS:raise ArithmeticError('native carrier source binding changed')
    carrier,ci,inp,verification,explicit=[cert.read(ART/n) for n in NAMES]
    for data in (carrier,ci,inp,verification):
        for name,h in data['bindings'].items():
            if cert.hashed(ROOT/name)!=h:raise ArithmeticError('upstream exact carrier input changed')
    for name,h in explicit['sources'].items():
        if cert.hashed(ROOT/name)!=h:raise ArithmeticError('explicit auxiliary proof source changed')
    g=carrier['geometry'];des=carrier['descent']
    if carrier['execution']['status']!='COMPLETE' or des['rank_lower_bound']!=3 or des['rank_upper_bound']!=3 or len(des['points'])!=2:raise ArithmeticError('rank3 auxiliary input required')
    helper=SourceFileLoader('native_carrier_pointed_helpers',str(HELPER)).load_module();R=PolynomialRing(QQ,'x');x=R.gen()
    quartic=R(list(reversed(g['quartic_coefficients'])));q=quartic[0].sqrt()
    if q not in QQ or not q or not quartic.is_squarefree():raise ArithmeticError('pointed rational infinity required')
    _,rebuilt,base,opp,constants=helper.pointed_curve(list(quartic),QQ(0),q,'x')
    J=EllipticCurve(QQ,g['minimal_Jacobian_model']);iso=J.isomorphism_to(base);generators=[J([QQ(v) for v in P]) for P in explicit['three_independent_Jacobian_points']]
    N=R(list(reversed(g['conic_parameter_numerator'])));den=R(list(reversed(g['conic_parameter_denominator'])));U=R(list(reversed(g['conic_root_numerator'])));scale=QQ(g['quartic_square_scale']);forms=[R(v) for v in g['forms']]
    if rebuilt!=quartic or U*U!=sum(forms[0][i]*N**i*den**(2-i) for i in range(3)) or quartic!=scale**2*sum(forms[1][i]*N**i*den**(2-i) for i in range(3)):raise ArithmeticError('carrier identity differs')
    A=R(inp['A']);B=R(inp['B']);covers={c['label']:c for c in inp['covers']};lifts=[]
    for c,form in zip(ci['covers'],forms):
        raw=covers[c['label']];factor=QQ(c['removed_rational_square_root']);pol=form*factor**2
        if pol!=R(raw['residual_chord']['q_coefficients']):raise ArithmeticError('primitive square scaling differs')
        x0,x1,y0,y1=[R(raw['lifted_section'][k+'_coefficients']) for k in ('x0','x1','y0','y1')]
        if y0*y0+pol*y1*y1-x0**3-3*pol*x0*x1*x1-A*x0-B or 2*y0*y1-3*x0*x0*x1-pol*x1**3-A*x1:raise ArithmeticError('generic native lift identity differs')
        lifts.append((factor,x0,x1,y0,y1))
    family=next(f for f in cert.read(spec.ATLAS)['families'] if f['family']=='08234');rows=[]
    for word in WORDS:
        point=sum((a*Q for a,Q in zip(word,generators)),J(0));image=helper.inverse_pointed(iso(point),constants);row={'word':word,'Jacobian_point':[str(v) for v in point.xy()] if point else None}
        if image is None:row['status']='BIRATIONAL_EXCEPTION';rows.append(row);continue
        z,w=image
        if w*w!=quartic(z):raise ArithmeticError('pointed inverse differs')
        if not den(z):row['status']='PARAMETER_POLE';rows.append(row);continue
        t=N(z)/den(z);roots=[U(z)/den(z),w/(scale*den(z))];s=-26*t-50
        if any(u*u!=f(t) for u,f in zip(roots,forms)):raise ArithmeticError('simultaneous square roots differ')
        height=int(max(abs(s.numerator()),s.denominator()));row.update(quartic_point=[str(z),str(w)],published_parameter=str(t),compact_parameter=str(s),square_roots=list(map(str,roots)),parameter_height=str(height),parameter_bits=height.bit_length(),beyond_compact_region=height>p['compact_inner_height'],is_known_origin=t==QQ(ci['anchor']['t']))
        if height.bit_length()>p['maximum_parameter_bits']:row['status']='PARAMETER_HEIGHT_CENSORED';rows.append(row);continue
        model,generic=spec.specialize(family,str(s));E=EllipticCurve(QQ,[QQ(str(v)) for v in model]);pub=EllipticCurve(QQ,[A(t),B(t)]);transport=pub.isomorphism_to(E);extra=[]
        for (factor,x0,x1,y0,y1),u in zip(lifts,roots):
            v=u*factor;Q=transport(pub(x0(t)+v*x1(t),y0(t)+v*y1(t)));extra.append(tuple(cert.F(str(a)) for a in Q.xy()))
        if any(a.denominator!=1 for a in model):raise ArithmeticError('integral compact model required')
        inv=cert.weierstrass_invariants(model);c4,c6=int(inv['c4']),int(inv['c6']);G=gcd(abs(c4),abs(c6));num=c6*c6;denominator=G**3;M=2**p['model_bits_gate']-1
        excluded=num>denominator*(1224*M+521)**2
        ceiling=(num+denominator-1)//denominator;lower=isqrt(ceiling);lower+=int(lower*lower<ceiling);coefficient=max(0,(lower-521+1223)//1224)
        bits=max(abs(int(a)).bit_length() for a in model)
        row.update(status='EXACT_TWO_POINT_IMAGE',curve=list(map(str,model)),supplied_points=[[str(v) for v in Q] for Q in extra],model_coefficient_bits=bits,normalized_coefficient_bits_lower_bound=coefficient.bit_length(),all_normalized_400_bit_integral_models_excluded=excluded,invariant_gcd=str(G),c4=str(c4),c6=str(c6))
        if height>p['compact_inner_height'] and bits<=p['model_bits_gate']:
            state=raw_state(model,tuple(generic)+tuple(extra),cache=ReductionCache(MemoryFactStore()),prime_bound=997);proof=checked_rank(model,state.basis,state.reductions.primes,state.no_two_torsion_prime)
            row.update(independent_points=[[str(v) for v in Q] for Q in state.basis],rank_certificate=proof,rank_lower_bound=state.rank)
        rows.append(row)
    return {'schema':'elliptic-curves.native-rank3-fullspan-images-result.v1','status':'PASS_FIXED_AUDIT','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'rows':rows,'labels':[c['label'] for c in ci['covers']],'pointed_quartic':list(map(str,quartic)),'actual_carrier_model':list(map(str,base.a_invariants())),'height_lemma':'For integral source invariants let G=gcd(abs(c4),abs(c6)). Any isomorphic integral model has abs(c6_target)^2 >= c6^2/G^3. A normalized integral equation with a1,a3 in{0,1}, a2 in{-1,0,1} and abs(a4),abs(a6)<=M has abs(c6_target)<=1224M+521. The recorded strict integer inequality excludes every such400-bit model, including normalized global minimal and integral short models.','claim_boundary':p['scope']}
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','build','check']);a=parser.parse_args()
    if a.stage=='prepare':prepare()
    else:
        if a.stage=='build' and OUT.exists():raise FileExistsError('preserve fixed image audit')
        d=expected()
        if a.stage=='check':
            if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('native image replay differs')
        else:checkpoint(OUT,d)
        print('NATIVE FULLSPAN IMAGES',[(r['word'],r['status'],r.get('parameter_bits'),r.get('normalized_coefficient_bits_lower_bound'),r.get('rank_lower_bound')) for r in d['rows']],flush=True)
