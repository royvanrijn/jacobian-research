#!/usr/bin/env sage-python
"""Twelve fixed actual-carrier images, with exact lifts and a model-height gate."""
import argparse,json,sys
from pathlib import Path
from importlib.machinery import SourceFileLoader
from sage.all import QQ,EllipticCurve,PolynomialRing
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
from research_runtime.store import checkpoint
from research_runtime.search_state import raw_state
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from memory_rank_certificate import checked_rank
ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/soluble-pair-carrier-height-v1';OUT=ART/'soluble_pair_carrier_height_v1.json'
HELPER=ROOT/'elkies-k3/scripts/prepare_r17_norm12_11952_alternate_v4_rank_one_bases.sage'
WORDS=[[-1,0],[1,0],[0,-1],[0,1],[-1,-1],[-1,1],[1,-1],[1,1],[-2,0],[2,0],[0,-2],[0,2]]

def sources():
    paths=[Path(__file__).resolve(),HELPER,spec.ATLAS,ART/'rank_jump_global_pair_solubility_v1.json',ART/'rank_jump_global_carrier_verification_v1.json',ART/'rank_jump_solubility_first_inputs_v1.json',ART/'rank_jump_local_solubility_block_inputs_v1.json',CAS/'compact_atlas_specialization.py',CAS/'memory_rank_certificate.py',CAS/'research_runtime/search_state.py']
    return {str(q.relative_to(ROOT)):cert.hashed(q) for q in paths}

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve actual-carrier height protocol')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.soluble-pair-carrier-height.v1','sources':sources(),'case':'observed_positive','labels':['orbit-1795d','orbit-0911e'],'words':WORDS,'wall_seconds':300,'rss_bytes':2147483648,'maximum_compact_parameter_bits':512,'search_model_bits_gate':360,'finite_certificate_prime_bound':997,'scope':'Use the actual globally soluble positive pair carrier from the newly completed auxiliary-Jacobian proof. Reverse its quartic at the retained rational infinity point, verify the pointed Jacobian and exact birational map, and map exactly twelve fixed combinations of the two already recorded Jacobian points. Verify both generic lift identities before specialization, transport to compact08234, and certify the rank of the17 generic sections plus the two supplied points by finite reductions through997. This is a bounded construction and height audit, not a score or point-search campaign. The pair and its origin are retrospective; its source known fibre is not a new discovery. No public target ranks enter the twelve words. No automatic parameter expansion, point search, descent, saturation, density law, record or universal novelty.'})

def expected():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['words']!=WORDS:raise ArithmeticError('frozen carrier inputs changed')
    result=cert.read(ART/'rank_jump_global_pair_solubility_v1.json');row=next(r for r in result['rows'] if r['id']==p['case']);g=row['geometry'];desc=row['descent']
    if row['execution']['status']!='COMPLETE' or not desc['global_carrier_solubility_proved'] or (desc['rank_lower_bound'],desc['rank_upper_bound'],desc['Sha_2_dimension'])!=(2,2,0) or len(desc['points'])!=2:raise ArithmeticError('completed soluble rank2 actual-carrier gate required')
    helper=SourceFileLoader('soluble_pair_pointed_helpers',str(HELPER)).load_module();R=PolynomialRing(QQ,'x');x=R.gen()
    quartic=R(list(reversed(g['quartic_coefficients'])));q=quartic[0].sqrt()
    if q not in QQ or not q or not quartic.is_squarefree():raise ArithmeticError('retained positive infinity point required')
    _,rebuilt,base,opp,constants=helper.pointed_curve(list(quartic),QQ(0),q,'x')
    J=EllipticCurve(QQ,g['minimal_Jacobian_model']);iso=J.isomorphism_to(base)
    generators=[J([QQ(v) for v in P]) for P in desc['points']]
    maps=cert.read(ART/'rank_jump_solubility_first_inputs_v1.json');covers={c['label']:c for c in cert.read(ART/'rank_jump_local_solubility_block_inputs_v1.json')['covers']};A=R(maps['A']);B=R(maps['B']);lifts=[]
    for label,form in zip(p['labels'],g['forms']):
        c=covers[label];lift=maps['split_lift_maps'][label]
        if c['form']!=form:raise ArithmeticError('carrier cover label differs')
        factor=QQ(c['removed_rational_square_root']);pol=R(form)*factor**2
        x0,x1,y0,y1=[R(lift[k+'_coefficients']) for k in ('x0','x1','y0','y1')]
        if y0*y0+pol*y1*y1-x0**3-3*pol*x0*x1*x1-A*x0-B or 2*y0*y1-3*x0*x0*x1-pol*x1**3-A*x1:raise ArithmeticError('generic lift polynomial identity failed')
        lifts.append((factor,x0,x1,y0,y1))
    N=R(list(reversed(g['parameter_numerator'])));den=R(list(reversed(g['parameter_denominator'])));U=R(list(reversed(g['conic_root_numerator'])));scale=QQ(g['quartic_square_scaling']);f1,f2=[R(v) for v in g['forms']]
    if rebuilt!=quartic or U*U!=sum(f1[i]*N**i*den**(2-i) for i in range(3)) or quartic!=scale**2*sum(f2[i]*N**i*den**(2-i) for i in range(3)):raise ArithmeticError('actual-carrier polynomial identities differ')
    family=next(f for f in cert.read(spec.ATLAS)['families'] if f['family']=='08234');rows=[]
    for word in p['words']:
        P=sum((a*Q for a,Q in zip(word,generators)),J(0));image=helper.inverse_pointed(iso(P),constants);item={'word':word,'Jacobian_point':[str(v) for v in P.xy()] if P else None}
        if image is None:item['status']='BIRATIONAL_EXCEPTION';rows.append(item);continue
        z,w=image
        if w*w!=quartic(z):raise ArithmeticError('inverse pointed map failed')
        if not den(z):item['status']='PARAMETER_POLE';rows.append(item);continue
        t=N(z)/den(z);roots=[U(z)/den(z),w/(scale*den(z))];s=-26*t-50
        if any(u*u!=f(t) for u,f in zip(roots,[f1,f2])):raise ArithmeticError('simultaneous rational split failed')
        item.update(quartic_point=[str(z),str(w)],published_parameter=str(t),compact_parameter=str(s),square_roots=list(map(str,roots)),compact_parameter_bits=int(max(abs(s.numerator()),s.denominator()).nbits()))
        if item['compact_parameter_bits']>p['maximum_compact_parameter_bits']:item['status']='DECLARED_PARAMETER_HEIGHT_CENSORED';rows.append(item);continue
        model,generic=spec.specialize(family,str(s));E=EllipticCurve(QQ,[QQ(str(v)) for v in model]);pub=EllipticCurve(QQ,[A(t),B(t)]);transport=pub.isomorphism_to(E);points=[]
        for (factor,x0,x1,y0,y1),u in zip(lifts,roots):
            v=u*factor;Q=pub(x0(t)+v*x1(t),y0(t)+v*y1(t));Q=transport(Q);points.append(tuple(cert.F(str(v)) for v in Q.xy()))
        state=raw_state(model,tuple(generic)+tuple(points),cache=ReductionCache(MemoryFactStore()),prime_bound=p['finite_certificate_prime_bound']);proof=checked_rank(model,state.basis,state.reductions.primes,state.no_two_torsion_prime)
        bits=max(max(abs(v.numerator),v.denominator).bit_length() for v in model)
        item.update(status='EXACT_TWO_LIFTS_AND_FINITE_POINT_PROOF',curve=list(map(str,model)),generic_points=[[str(v) for v in Q] for Q in generic],supplied_points=[[str(v) for v in Q] for Q in points],independent_points=[[str(v) for v in Q] for Q in state.basis],rank_certificate=proof,rank_lower_bound=state.rank,model_coefficient_bits=bits,passes_declared_height_gate=bits<=p['search_model_bits_gate'])
        rows.append(item)
    return {'schema':'elliptic-curves.soluble-pair-carrier-height-result.v1','status':'PASS_FIXED_AUDIT','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'labels':p['labels'],'actual_carrier_curve':list(map(str,base.ainvs())),'pointed_quartic':list(map(str,quartic)),'rows':rows,'claim_boundary':p['scope']}

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['prepare','build','check']);v=a.parse_args()
    if v.stage=='prepare':prepare()
    else:
        if v.stage=='build' and OUT.exists():raise FileExistsError('preserve carrier height audit')
        d=expected()
        if v.stage=='check':
            if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('carrier height replay differs')
        else:checkpoint(OUT,d)
        print('ACTUAL CARRIER HEIGHT AUDIT',[(r['word'],r['status'],r.get('rank_lower_bound'),r.get('model_coefficient_bits')) for r in d['rows']],flush=True)
