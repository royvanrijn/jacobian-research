#!/usr/bin/env python3
"""Production local/subspace-first arithmetic on the MW18 priority cohort.

The old unconditional ellrank campaign is retained with its original sources
and certificate. This runner defaults to cheap local scheduling features.
Full Selmer is explicit and reserved for an upper-bound requirement; supplied
subspace requests use the shared vectorized local/cover/CT pipeline.
"""
import argparse
from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path
import sys

from research_runtime.arithmetic import CurveModel
from research_runtime.store import checkpoint, digest as content_digest
from research_runtime.supervisor import Limits, supervise_source, preserve_previous

ROOT=Path(__file__).resolve().parents[2]
CAS=ROOT/'elliptic-curves/cas'
SPECIALIZATIONS=ROOT/'artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-specializations-h1000-v1.json'
OUTPUT=ROOT/'artifacts/local/elkies-k3/r17-extreme-anchored-mw18-production-arithmetic-v2.json'
CHECKPOINT_DIRECTORY=ROOT/'artifacts/local/elkies-k3/r17-extreme-anchored-mw18-production-arithmetic'
DEFAULT_SAGE=Path('/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python')
PROTOCOL='R17ANCHORMW18ARITHMETIC'

WORKER=r'''import json,sys
from pathlib import Path
sys.path.insert(0,CAS_PATH)
from run_arithmetic_pipeline import worker
from research_runtime.store import default_store,checkpoint
from research_runtime.supervisor import Limits
payload=json.loads(Path(INPUT_PATH).read_text())
result=worker(payload['request'],default_store(),Limits(**payload['limits']),retained=payload.get('retained'))
checkpoint(OUTPUT_PATH,result)
print(result['status'],flush=True)
'''


def digest(path):return sha256(Path(path).read_bytes()).hexdigest()
def relative(path):
    try:return str(Path(path).resolve().relative_to(ROOT))
    except ValueError:return str(Path(path).resolve())


def raw_model_key(candidate) -> str:
    return sha256(
        json.dumps(candidate["raw_short_model"], separators=(",", ":")).encode()
    ).hexdigest()


def priority_cohort(candidates):
    certified = [
        candidate
        for candidate in candidates
        if candidate["independence"]["status"]
        == "CERTIFIED_INTEGRALLY_INDEPENDENT_RANK_AT_LEAST_18"
    ]
    by_cover = {}
    for candidate in certified:
        label = candidate["cover_label"]
        if label not in by_cover or candidate["nagao"]["total_score"] > by_cover[label]["nagao"]["total_score"]:
            by_cover[label] = candidate
    selected = [(candidate, "top_nagao_score_on_cover") for candidate in by_cover.values()]
    selected.extend(
        (candidate, "known_extreme_anchor_surviving_nagao")
        for candidate in certified
        if candidate["nagao"]["is_certified_anchor"]
    )
    groups = {}
    for candidate, reason in selected:
        key = raw_model_key(candidate)
        group = groups.setdefault(
            key,
            {
                "representative": candidate,
                "candidate_ids": [],
                "selection_reasons": [],
            },
        )
        if candidate["candidate_id"] not in group["candidate_ids"]:
            group["candidate_ids"].append(candidate["candidate_id"])
        if reason not in group["selection_reasons"]:
            group["selection_reasons"].append(reason)
    result = list(groups.values())
    result.sort(
        key=lambda group: (
            -float(group["representative"]["nagao"]["total_score"]),
            group["representative"]["candidate_id"],
        )
    )
    return result



def arithmetic_request(candidate,mode,subspaces):
    model=candidate['raw_short_model']
    if len(model)==2:model=[0,0,0,*model]
    if mode=='subspace':
        request=subspaces[candidate['candidate_id']]
        target=request.get('target',request['source'])
        if request.get('mode')!='subspace' or CurveModel(target['model'])!=CurveModel(model):
            raise ValueError('subspace request targets a different MW18 curve')
        return request
    result={'mode':mode,'curve':{'model':model},'feature_primes':[2,3,5,7,11,13]}
    if mode=='complete-selmer':result['requirement']='unconditional-upper-bound'
    return result


def gate(worker):
    # These curves have no rational 2-torsion, certified in the source audit.
    total=(worker['selmer']['full_selmer_dimension'] if worker and worker['mode']=='complete-selmer' else None)
    if total is not None and total<18:raise ArithmeticError('completed upper bound contradicts the certified basis')
    excluded=total is not None and total<32
    return {'known_rank_lower_bound':18,'target_rank':32,'rank_upper_bound':total,
            'residual_dimension_upper_bound':None if total is None else total-18,
            'mathematically_excluded':excluded,
            'bounded_point_search_possible':not excluded,
            'full_selmer_required_for_bounded_search':False,
            'decision':'EXCLUDED_BY_COMPLETE_SELMER' if excluded else 'UNKNOWN',
            'unbounded_search_authorized':False}


def supervise(candidate,request,args,retained=None):
    directory=args.checkpoint_directory/(candidate['candidate_id']+'-'+content_digest(request)[:16])
    result_path=directory/('replayed.json' if retained else 'result.json')
    limits=Limits(args.timeout_seconds,args.rss_limit_bytes,pari_stack_bytes=args.pari_stack_bytes)
    receipt=supervise_source(str(args.sage_python),f'CAS_PATH={str(CAS)!r}\n'+WORKER,
        {'request':request,'limits':asdict(limits),'pari_stack_bytes':args.pari_stack_bytes,'retained':retained},
        result_path,directory/'worker.log',timeout=args.timeout_seconds,rss_limit_bytes=args.rss_limit_bytes)
    result=json.loads(result_path.read_text()) if receipt['outcome']=='completed' else None
    return {'candidate_id':candidate['candidate_id'],'request':request,'arithmetic':result,
            'supervisor':receipt,'gate':gate(result)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--specializations',type=Path,default=SPECIALIZATIONS)
    parser.add_argument('--output',type=Path,default=OUTPUT)
    parser.add_argument('--checkpoint-directory',type=Path,default=CHECKPOINT_DIRECTORY)
    parser.add_argument('--sage-python',type=Path,default=DEFAULT_SAGE)
    parser.add_argument('--timeout-seconds',type=float,default=60)
    parser.add_argument('--rss-limit-bytes',type=int,default=1073741824)
    parser.add_argument('--pari-stack-bytes',type=int,default=256000000)
    parser.add_argument('--mode',choices=['features','subspace','complete-selmer'],default='features')
    parser.add_argument('--subspace-requests',type=Path,help='candidate-id mapping to known-squareclass pipeline requests')
    parser.add_argument('--candidate',action='append',help='restrict to declared priority candidate ids')
    parser.add_argument('--no-resume',action='store_true')
    parser.add_argument('--verify',type=Path,help='replay retained local/subspace witnesses; full-Selmer checks retained integrity')
    args=parser.parse_args()
    if args.mode=='subspace' and not args.subspace_requests and not args.verify:parser.error('subspace mode requires supplied global squareclasses')
    source=json.loads(args.specializations.read_text())
    if source['status']!='COMPLETE_EXACT_MW18_FINALIST_SPECIALIZATION_AUDIT' or source['certified_rank_at_least_18_count']!=178:
        raise ArithmeticError('MW18 specialization source is incomplete')
    groups=priority_cohort(source['candidates'])
    if args.candidate:
        selected=set(args.candidate)
        groups=[g for g in groups if g['representative']['candidate_id'] in selected]
        if {g['representative']['candidate_id'] for g in groups}!=selected:raise ValueError('candidate is not in the priority cohort')
    subspaces=json.loads(args.subspace_requests.read_text()) if args.subspace_requests else {}
    previous=json.loads(args.verify.read_text()) if args.verify else None
    inputs={relative(args.specializations):digest(args.specializations)}
    if previous and previous['inputs']!=inputs:raise ArithmeticError('retained cohort source changed')
    old={r['candidate_id']:r for r in previous['records']} if previous else {}
    records=[]
    for group in groups:
        candidate=group['representative'];cid=candidate['candidate_id']
        if previous and cid not in old:continue
        request=old[cid]['request'] if previous else arithmetic_request(candidate,args.mode,subspaces)
        path=args.checkpoint_directory/(cid+'-'+content_digest(request)[:16])/'row.json'
        if not previous and not args.no_resume and path.exists():
            record=json.loads(path.read_text())
            if record['request']!=request:raise ArithmeticError('cached request changed')
        else:
            retained=old[cid]['arithmetic'] if previous else None
            if previous and retained is None:
                records.append(old[cid]);continue
            record=supervise(candidate,request,args,retained)
            if not previous:preserve_previous(path);checkpoint(path,record)
        records.append(record)
        result={'schema':'elliptic-curves.mw18-production-arithmetic.v2','inputs':inputs,'records':records,
                'status':'COMPLETE_REQUESTED_COHORT' if len(records)==len(groups) else 'PARTIAL',
                'claim_boundary':'Local/subspace dimensions are not full Selmer bounds. Missing arithmetic is UNKNOWN and does not require a full BNF before a separately bounded point search.'}
        preserve_previous(args.output);checkpoint(args.output,result)
        print(f"{PROTOCOL}|candidate={cid}|mode={request['mode']}|outcome={record['supervisor']['outcome']}|upper={record['gate']['rank_upper_bound']}",flush=True)

if __name__=='__main__':main()
