#!/usr/bin/env python3
"""Bounded lazy half-lattice search consuming and returning immutable MWState.

The request fixes the initial point basis, a positive definite scoring metric,
chart policy and finite budgets. No minimal model, factorization, BNF or full
Selmer computation is a prerequisite. Replay checks retained hits and points;
regeneration is a separate discovery invocation.
"""
import argparse
from dataclasses import asdict
from fractions import Fraction
import json
from pathlib import Path
import sys

from research_runtime.chart_policy import ChartPolicy
from research_runtime.cvp import VoronoiIterator
from research_runtime.mw_state import MWState
from research_runtime.pruning import PruningRegistry, SearchRequest
from research_runtime.search_state import raw_state, reduction_cache
from research_runtime.search_state import input_generators
from research_runtime.arithmetic import CurveModel, rationals
from research_runtime.store import checkpoint, digest, default_store
from research_runtime.supervisor import Limits, capture, preserve_previous
from pointed_quartic_search import PointedQuarticSearch, ROOT, sources


def search_request(request, output, *, retained=None):
    if request.get('schema') != 'elliptic-curves.lazy-mw-search.v1':
        raise ValueError('unknown lazy MW search schema')
    cache=reduction_cache()
    if retained is not None:
        if retained['request_hash']!=digest(request):raise ArithmeticError('search request changed')
        cache.store.import_snapshot(retained['arithmetic_facts'])
        state=MWState.from_record(retained['initial_state'],cache=cache)
    else:
        state=raw_state(request['curve'],request['points'],cache=cache,prime_bound=request.get('prime_bound',1000))
    if state.rank!=len(request['points']):
        raise ArithmeticError('initial list is not certified independent at the declared primes')
    model=tuple(request['curve'])
    if len(model)==2:model=(0,0,0,*model)
    points=tuple(rationals((p['x'],p['y']) if isinstance(p,dict) else p) for p in request['points'])
    if state.model!=CurveModel(model) or input_generators(state)!=points:
        raise ArithmeticError('initial MWState differs from the requested curve or generators')
    policy=ChartPolicy(**request.get('policy',{}))
    coordinate={'kind':policy.chart_metric_kind,'weight':policy.chart_metric_weight,'matrix':[1,0,0,1]}
    if policy.model_normalization!='raw' or policy.chart_parameterization!='pointed-quartic':
        raise ValueError('this production adapter uses raw pointed charts; normalization comparators are separate policies')
    if policy.enumeration_backend != 'gmp-pointed-sieve':
        raise ValueError('unsupported pointed enumeration backend')
    policy_record=asdict(policy)
    count=request['next_holes'];height=request['height'];seconds=float(request['seconds_per_chart'])
    if any(type(n) is not int or n<1 for n in (count,height,request['cvp_node_budget'])) or seconds<=0:
        raise ValueError('positive declared chart, height, time and CVP budgets are required')
    gram=request['metric_gram']
    if len(gram)!=state.rank:raise ValueError('metric dimensions differ from the initial basis')
    registry=PruningRegistry(ROOT)
    scope=(('curve',state.model.key),)
    initial_gate=registry.decision(SearchRequest('rank_target',scope,target_rank=state.rank+1))
    if retained is None:
        binding={'state':state.key,'policy':policy.key}
        iterator=(VoronoiIterator.resume(request['cvp_checkpoint'],binding=binding)
                  if 'cvp_checkpoint' in request else
                  VoronoiIterator(gram,seen=request.get('seen_masks',()),binding=binding))
        if iterator.gram!=tuple(tuple(Fraction(str(v)) for v in row) for row in gram):
            raise ArithmeticError('resumed CVP metric differs from the request')
        cvp_status='COMPLETE_REQUEST'
        try:
            holes=(iterator.next_holes(count,diversity_window=policy.diversity_window,node_budget=request['cvp_node_budget'])
                   if initial_gate['search_allowed'] else [])
        except TimeoutError:
            holes=[];cvp_status='CVP_BUDGET_REACHED'
        frontier=iterator.checkpoint()
    else:
        # The metric is a scheduling policy, not an exclusion theorem. Replay
        # need not revisit the CVP frontier to verify every searched cover.
        from research_runtime.cvp import Hole
        holes=[Hole.from_record(h) for h in retained['holes']]
        frontier=retained['cvp_checkpoint']
        VoronoiIterator.resume(frontier,binding={'state':state.key,'policy':policy.key})
        cvp_status=retained['cvp_status']
    initial=state
    result={'schema':'elliptic-curves.lazy-mw-search-result.v1','request':request,'request_hash':digest(request),
        'initial_state':initial.record(),'policy':policy_record,'metric_is_height_certificate':False,
        'holes':[h.record() for h in holes],'cvp_checkpoint':frontier,'cvp_status':cvp_status,'charts':[],'status':'RUNNING',
        'source_hashes':sources(),'full_bnf_requested':False}
    if retained is None and not initial_gate['search_allowed']:
        result['status']='EXCLUDED_BY_THEOREM';result['pruning']=initial_gate
    elif retained is not None and retained['status']=='EXCLUDED_BY_THEOREM' and not holes:
        result['status']=retained['status'];result['pruning']=retained['pruning']
    for i,hole in enumerate(holes):
        if retained is not None and i>=len(retained['charts']):
            if retained['status']!='EXCLUDED_BY_THEOREM':raise ArithmeticError('missing retained chart')
            result['status']=retained['status'];result['pruning']=retained['pruning'];break
        gate=registry.decision(SearchRequest('rank_target',scope,target_rank=state.rank+1))
        if retained is None and not gate['search_allowed']:
            result['status']='EXCLUDED_BY_THEOREM';result['pruning']=gate;break
        # Point admission only appends to this basis. Old planned centres retain
        # their exact meaning while discoveries are accumulated incrementally.
        representative=(*hole.doubled_coordinates,*(0 for _ in range(state.rank-initial.rank)))
        search=PointedQuarticSearch(state=state,centre={'coefficients':representative},coordinate_policy=coordinate)
        if retained is None:
            updated,outcome=search.search_state(height,seconds,
                extra_primes=request.get('extra_primes',()),checkpoint_dir=Path(output).parent/'charts')
        else:
            old=retained['charts'][i]
            if old['hole']!=hole.record():raise ArithmeticError('changed chart centre')
            outcome=search.verify_record(old['search'])
            updated=state
            for point in outcome.curve_points:
                updated=updated.adjoin(point,cache=cache,extra_primes=request.get('extra_primes',()))
        row={'hole':hole.record(),'search':outcome.record,'state_before':state.key,'state_after':updated.key,
             'certified_rank_gain':updated.rank-state.rank}
        if retained is not None and row!=retained['charts'][i]:raise ArithmeticError('point-admission replay differs')
        result['charts'].append(row);state=updated
        result.update({'final_state':state.record(),'arithmetic_facts':cache.store.snapshot()})
        if retained is None:checkpoint(output,result)
    if result['status']=='RUNNING':
        result['status']=cvp_status if cvp_status!='COMPLETE_REQUEST' else (
            'COMPLETE_DECLARED_CHARTS' if all(r['search']['status']=='bounded_search_complete' for r in result['charts']) else 'INCOMPLETE_CHARTS')
    result.update({'final_state':state.record(),'arithmetic_facts':cache.store.snapshot(),
                   'certified_rank_gain':state.rank-initial.rank,
                   'requires_new_metric_after_basis_growth':state.rank!=initial.rank,
                   'claim_boundary':'Exact independent points and bounded coverage only; a miss is not a rank upper bound.'})
    if retained is not None and result['final_state']!=retained['final_state']:
        raise ArithmeticError('final MWState did not replay')
    checkpoint(output,result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--request',type=Path)
    group.add_argument('--verify',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--wall-seconds',type=int,default=120)
    parser.add_argument('--rss-bytes',type=int,default=1_073_741_824)
    parser.add_argument('--worker',action='store_true',help=argparse.SUPPRESS)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('use a new output; retained discovery witnesses are immutable')
    if not args.worker:
        result=capture([sys.executable,str(Path(__file__).resolve()),*sys.argv[1:],'--worker'],
            limits=Limits(args.wall_seconds,args.rss_bytes),log_path=args.output.with_suffix('.log'))
        print(result.stdout,end='');return
    retained=json.loads(args.verify.read_text()) if args.verify else None
    request=retained['request'] if retained else json.loads(args.request.read_text())
    result=search_request(request,args.output,retained=retained)
    print(f"LAZY_MW_SEARCH|status={result['status']}|charts={len(result['charts'])}|rank_gain={result['certified_rank_gain']}")


if __name__=='__main__':main()
