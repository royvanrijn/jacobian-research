#!/usr/bin/env python3
"""X1092 Elkies-pair scheduler v2: choose usable carriers, never starve V3.

This is a thin policy layer over run_x1092_elkies_pairs.py.  The certified
conic construction, carrier geometry, M19 admission and V3 implementation are
unchanged.  New campaigns select among a bounded set of successful carriers
using only induced parameter heights, and reduce the requested control count
when necessary so at least one ranked fibre survives pre-score attrition.
"""
from __future__ import annotations

from fractions import Fraction as F
import json
from pathlib import Path

import run_x1092_elkies_pairs as base
from elkies_pair_policy import (
    carrier_quality_key,
    carrier_quality_summary,
    effective_control_count,
    final_selection,
    mestre_score,
    pair_order,
    rational_height,
    require,
    short_isomorphic,
    shortlist,
    split_controls,
)

SELF=Path(__file__).resolve()
QUALITY_ADDRESSES=32
MAX_SAME_CARRIER_SUCCESSES=1
MAX_ODD_CARRIER_SUCCESSES=8
QUALITY_SECONDS=30


def carrier_quality(data, plan):
    """Measure only the parameter-height yield of a certified generic carrier."""
    heights=[];errors=0
    for unused,value,error in base.carrier_addresses(
            data,min(QUALITY_ADDRESSES,plan['addresses']),plan['parameter_bits']):
        if error or value is None:
            errors+=1;continue
        t,unused_w1,unused_w2=value
        heights.append(rational_height(t))
    summary=carrier_quality_summary(heights,plan['parameter_bits'])
    summary.update(requested=min(QUALITY_ADDRESSES,plan['addresses']),errors=errors)
    return summary


def carrier_phase(folder):
    """Prefer low-height yield; same-cover wins only an outcome-free tie."""
    from sage.all import QQ,matrix,ZZ,vector
    import elkies_pair_geometry as geo

    plan,rt,out,parent=base.context(folder);E,points,G=geo.load_parent(parent)
    cs=base.conics(out,geo,E,points,G);by={c['mask']:c for c in cs}
    order=pair_order(cs,parent['generic_height_gram'],plan['pair_limit'])
    base.write(out/'pair-order.json',{'rows':order})
    attempts=[];successes=[];same_successes=0;odd_successes=0

    for i,row in enumerate(order):
        if row['kind']=='same-cover' and same_successes>=MAX_SAME_CARRIER_SUCCESSES:
            continue
        if row['kind']=='odd-intersection' and odd_successes>=MAX_ODD_CARRIER_SUCCESSES:
            break
        a,b=(by[m] for m in row['masks']);dest=out/'pairs'/f'{i:04d}'
        try:
            if row['kind']=='same-cover':
                field=E.base_ring();T=geo.function(field,a['T']);W=geo.function(field,a['W'])
                scale=QQ(row['scale']);data=None
                for j,u in enumerate(plan['same_cover_addresses']):
                    try:t,w=T(QQ(u)),W(QQ(u))
                    except (ValueError,ZeroDivisionError):continue
                    packet=base.bounded_math(
                        90,lambda:base.certify_seed(parent,t,w,scale*w,a,b,dest/f'probe-{j}'))
                    if packet:
                        data={'kind':'rational','masks':row['masks'],'T':a['T'],'W1':a['W'],
                              'scale':str(scale),'function_field_rank_lower_bound':19,
                              'independence':'Exact generic section identities and one specialization with certified independent 19 images.',
                              'witness_seed':str((dest/f'probe-{j}/seed.json').relative_to(rt)),
                              'witness_seed_sha256':base.sha(dest/f'probe-{j}/seed.json')}
                        break
                require(data is not None,'no independent same-cover specialization in the fixed probe window')
            else:
                data=base.bounded_math(90,lambda:geo.carrier(a,b))
                w,v=vector(ZZ,a['word']),vector(ZZ,b['word'])
                gram=matrix(QQ,19,19);gram[:17,:17]=4*G
                for j,z in enumerate((w,v)):
                    for k,value in enumerate(2*G*z):gram[k,17+j]=gram[17+j,k]=value
                    gram[17+j,17+j]=16
                gram[17,18]=gram[18,17]=w*G*v
                require(gram.is_positive_definite() and gram.det()==4**17*1092*36,
                        'incorrect biquadratic height Gram')
                data['generic_height_gram']=[list(map(str,r)) for r in gram.rows()]
                data['generic_height_determinant']=str(gram.det())

            data.update(
                conic_sha256={str(c['mask']):base.sha(out/'conics'/f"{c['mask']}.json") for c in (a,b)},
                pair_rule=row,uses_known_fibre=False)
            try:
                quality=base.bounded_math(QUALITY_SECONDS,lambda:carrier_quality(data,plan))
            except (ValueError,ArithmeticError,ZeroDivisionError,TimeoutError):
                quality=carrier_quality_summary([],plan['parameter_bits'])
                quality.update(requested=min(QUALITY_ADDRESSES,plan['addresses']),errors='QUALITY_UNKNOWN')
            data['carrier_quality']=quality
            base.write(dest/'carrier.json',data)
            result={'index':i,**row,'status':'PASS_GENERIC_M19_CARRIER','quality':quality}
            attempts.append(result);successes.append((i,row,data))
            if row['kind']=='same-cover':same_successes+=1
            else:odd_successes+=1
            base.write(out/'carrier-attempts.json',{'results':attempts},False)
            print('PAIR_CARRIER_CANDIDATE',data['kind'],row['masks'],json.dumps(quality,sort_keys=True),flush=True)
        except (ValueError,ArithmeticError,ZeroDivisionError,TimeoutError) as e:
            result={'index':i,**row,'status':'UNKNOWN_BOUNDED_PAIR','reason':str(e)}
            base.write(dest/'unknown.json',result);attempts.append(result)
            base.write(out/'carrier-attempts.json',{'results':attempts},False)
            print('PAIR_UNKNOWN',json.dumps(result),flush=True)

    if not successes:
        base.write(out/'carrier-status.json',{
            'status':'NO_CERTIFIED_CARRIER_IN_FINITE_WINDOW','attempts':len(attempts),
            'rank_or_existence_upper_bound':None})
        return

    def key(item):
        i,row,data=item
        return (*carrier_quality_key(data['carrier_quality'],row['kind']=='same-cover'),-i)
    i,row,data=max(successes,key=key)
    base.write(out/'carrier.json',data)
    base.write(out/'carrier-selection.json',{
        'status':'PASS_OUTCOME_FREE_LOW_HEIGHT_CARRIER_SELECTION',
        'selected_index':i,'selected_masks':row['masks'],'selected_kind':row['kind'],
        'quality_addresses':min(QUALITY_ADDRESSES,plan['addresses']),
        'primary_cap':plan['parameter_bits'],'extended_cap':4*plan['parameter_bits'],
        'same_cover_tie_break_only':True,
        'candidates':[{'index':j,'kind':r['kind'],'masks':r['masks'],
                       'quality':d['carrier_quality']} for j,r,d in successes],
        'boundary':'Carrier scheduling uses only map-defined parameter heights; no Mestre score, specialized rank outcome or record fibre enters selection.'})
    print('PAIR_CARRIER_SELECTED',data['kind'],row['masks'],json.dumps(data['carrier_quality'],sort_keys=True),flush=True)


def score_phase(folder):
    """v1 scoring with deterministic graceful degradation of control count."""
    from sage.all import prime_range

    plan,rt,out,parent=base.context(folder);data=base.read(out/'carrier.json');rows=[];attrition=[];seen={}
    exclusions=[]
    for aa in base.read(rt/'pair-inputs/exclusions.json')['models']:
        a1,a2,a3,a4,a6=map(F,aa);b2=a1*a1+4*a2;b4=a1*a3+2*a4;b6=a3*a3+4*a6
        exclusions.append((-(b2*b2-24*b4)/48,-(-b2**3+36*b2*b4-216*b6)/864))

    for i,value,error in base.carrier_addresses(data,plan['addresses'],plan['parameter_bits']):
        if error:attrition.append({'index':i,'reason':error});continue
        t,w1,w2=value
        if rational_height(t)>plan['parameter_bits']:
            attrition.append({'index':i,'reason':'PARAMETER_HEIGHT_CAP','bits':rational_height(t)});continue
        try:a,b=base.equation(parent,t)
        except (ValueError,ZeroDivisionError) as e:attrition.append({'index':i,'reason':str(e)});continue
        j=F(1728*4*a**3,4*a**3+27*b*b)
        if any(short_isomorphic(a,b,c,d) for c,d in seen.get(str(j),[])):
            attrition.append({'index':i,'reason':'EXACT_Q_ISOMORPHISM_DUPLICATE'});continue
        if any(4*c**3+27*d*d!=0 and F(1728*4*c**3,4*c**3+27*d*d)==j and
               short_isomorphic(a,b,c,d) for c,d in exclusions):
            attrition.append({'index':i,'reason':'PINNED_EQUATION_EXCLUSION'});continue
        seen.setdefault(str(j),[]).append((a,b))
        rows.append({'id':str(i),'parameter':str(t),'w1':str(w1),'w2':str(w2),
                     'model':[str(a),str(b)],'j':str(j),'j_bits':rational_height(j)})

    base.write(out/'address-roster.json',{
        'rows':rows,'attrition':attrition,'requested_addresses':plan['addresses'],
        'boundary':'Fixed generator multiples with exact-isomorphism/height attrition, not uniformly random fibres.'})
    if not rows:
        base.write(out/'selection.json',{'rows':[],'status':'NO_FINITE_ADDRESSES_WITHIN_FIXED_HEIGHT_CAP'})
        return

    requested=plan['controls'];effective=effective_control_count(len(rows),requested)
    controls,_=split_controls(rows,effective)
    base.write(out/'controls.json',{
        'ids':[r['id'] for r in controls],'requested':requested,'effective':effective,
        'degraded':effective!=requested,
        'rule':'Controls are reduced only by pre-score eligible population size so at least one ranked fibre survives.'})

    cache={};small=list(prime_range(5,plan['shallow_bound']+1));large=list(prime_range(plan['shallow_bound']+1,plan['deep_bound']+1))
    for r in rows:
        r['shallow_traces']=base.traces(tuple(map(int,r['model'])),small,cache)
        r['shallow']=mestre_score(r['shallow_traces'])
    base.write(out/'shallow.json',{'rows':rows,'bound':plan['shallow_bound']})
    deep=shortlist(rows,controls,plan['deep_keep'])
    for r in deep+controls:
        r['extension_traces']=base.traces(tuple(map(int,r['model'])),large,cache)
        r['deep']=mestre_score(r['shallow_traces']+r['extension_traces'])
    base.write(out/'deep.json',{'rows':deep+controls,'bound':plan['deep_bound'],'validation_primes_read':[]})
    base.write(out/'selection.json',{
        'rows':final_selection(deep,controls,plan['ranked']),
        'status':'PASS_FIXED_TWO_STAGE_SELECTION',
        'requested_controls':requested,'effective_controls':effective,
        'control_degraded':effective!=requested})


# Make every supervised _phase child execute this entrypoint and these two policies.
base.SELF=SELF
base.carrier_phase=carrier_phase
base.score_phase=score_phase


if __name__=='__main__':
    base.main()
