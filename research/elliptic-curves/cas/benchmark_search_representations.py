#!/usr/bin/env python3
"""Bounded independent normalization/chart/enumerator comparisons.

A request supplies one exact known subgroup, centre and finite budgets. Each
combination runs in its own supervised worker. Censored normalization or search
is a measurement, never a reason to call the curve insoluble or rank bounded.
"""
import argparse
from dataclasses import asdict
from fractions import Fraction
from math import isqrt,gcd
import json
from pathlib import Path
import sys

from research_runtime.arithmetic import ArithmeticContext, CurveModel
from research_runtime.chart_policy import ChartPolicy, RepresentationPipeline
from research_runtime.mw_state import MWState
from research_runtime.search_state import raw_state,reduction_cache
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import Limits,capture,run
from pointed_quartic_search import PointedQuarticSearch


def cell(request,policy,retained=None):
    cache=reduction_cache()
    if retained:
        if retained['request_hash']!=digest(request) or retained['policy']!=asdict(policy):raise ArithmeticError('benchmark replay request differs')
        cache.store.import_snapshot(retained['arithmetic_facts'])
        state=MWState.from_record(retained['initial_state'],cache=cache)
        model=tuple(request['curve']);model=(0,0,0,*model) if len(model)==2 else model
        if state.model!=CurveModel(model) or state.basis!=tuple(tuple(str(Fraction(str(c))) for c in p) for p in request['points']):raise ArithmeticError('benchmark state differs from requested basis')
    else:state=raw_state(request['curve'],request['points'])
    if state.rank!=len(request['points']):raise ArithmeticError('benchmark basis is not certified independent')
    evidence={}
    def normalize(state,limits):
        if policy.model_normalization=='raw':
            context=state.arithmetic;model=state.model;transport=(Fraction(1),Fraction(0),Fraction(0),Fraction(0))
        else:
            if retained:context=ArithmeticContext.from_record(retained['normalization_context'])
            else:
                from research_runtime.sage_arithmetic import SageArithmetic
                context=SageArithmetic(cache.store).prepare_context(state.arithmetic,factor_primes=request.get('factor_primes',()),discover=True)
            if context.model!=state.model:raise ArithmeticError('normalization changed input model')
            context.require_prepared();model=context.minimal_model;transport=tuple(map(Fraction,context.minimal_to_input))
        u,r,s,t=transport
        points=[((Fraction(x)-r)/u**2,(Fraction(y)-s*(Fraction(x)-r)-t)/u**3) for x,y in state.basis]
        evidence['normalization_context']=context.record()
        return model,points,transport
    def parameterize(state,normalized,centre,policy,limits):
        model,points,transport=normalized
        search=PointedQuarticSearch(model.coefficients,points,centre,
            {'kind':policy.chart_metric_kind,'weight':policy.chart_metric_weight,'matrix':[1,0,0,1]})
        return search,transport
    def enumerate_points(state,chart,limits):
        search,transport=chart;u,r,s,t=transport
        if policy.enumeration_backend=='gmp-pointed-sieve':
            result=(search.verify_record(retained['witnesses']['search']) if retained else
                    search.search(request['height'],float(request['seconds'])))
            points=result.curve_points;witness={'search':result.record}
        else:
            if retained:
                witness=retained['witnesses']
                if witness['chart']!=search.chart_record():raise ArithmeticError('retained PARI chart differs')
                hits=witness['primitive_square_hits']
            else:
                polynomial='+'.join(f'({c})*x^{i}' for i,c in enumerate(search.coefficients))
                program=f'P=hyperellratpoints({polynomial},{request["height"]});for(i=1,#P,print("HIT|",P[i][1],"|",P[i][2]));quit\n'
                result=capture(['gp','-fq'],input_text=program,limits=Limits(float(request['seconds']),request['rss_bytes']))
                hits=[]
                for line in result.stdout.splitlines():
                    if not line.startswith('HIT|'):continue
                    _,xx,yy=line.split('|');xx=Fraction(xx);yy=Fraction(yy)
                    rr=yy*xx.denominator**2
                    if rr.denominator!=1:raise ArithmeticError('nonintegral homogenized root')
                    hits.append([str(xx.numerator),str(xx.denominator),str(abs(rr.numerator))])
                lead=search.coefficients[4]
                if lead>=0 and isqrt(lead)**2==lead:hits.append(['1','0',str(isqrt(lead))])
                witness={'chart':search.chart_record(),'primitive_square_hits':hits,
                         'supervisor':result.supervision,'bounded_census_replayed':False}
            points=set()
            for n,d,root in (tuple(map(int,h)) for h in hits):
                if (d==0 and n!=1) or (d!=0 and (not 1<=d<=request['height'] or abs(n)>request['height'] or gcd(n,d)!=1)):
                    raise ArithmeticError('retained PARI hit lies outside its primitive search box')
                if root<0 or root*root!=sum(c*n**i*d**(4-i) for i,c in enumerate(search.coefficients)):
                    raise ArithmeticError('retained PARI square identity failed')
                for sign in {root,-root}:
                    point=search.map_hit(n,d,sign)
                    if point is not None:points.add(point)
        updated=state
        for x,y in sorted(points):
            x,y=Fraction(x),Fraction(y)
            updated=updated.adjoin((u*u*x+r,u**3*y+s*u*u*x+t),cache=cache,extra_primes=request.get('extra_primes',()))
        return updated,witness
    pipeline=RepresentationPipeline(normalizers={'raw':normalize,'minimal':normalize},
        parameterizations={'pointed-quartic':parameterize},
        enumerators={'gmp-pointed-sieve':enumerate_points,'pari-hyperellratpoints':enumerate_points})
    final,row=pipeline.run(state,policy,centre=request['centre'],limits=request)
    result={**row,**evidence,'request_hash':digest(request),'request':request,'initial_state':state.record(),
            'final_state':final.record(),'arithmetic_facts':cache.store.snapshot(),'certified_rank_gain':final.rank-state.rank}
    if retained and result['final_state']!=retained['final_state']:raise ArithmeticError('benchmark point admission differs')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--request',type=Path);p.add_argument('--verify',type=Path)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--cell',help=argparse.SUPPRESS)
    p.add_argument('--sage-python',default='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python')
    a=p.parse_args()
    if a.cell:
        retained=json.loads(a.verify.read_text()) if a.verify else None
        request=retained['request'] if retained else json.loads(a.request.read_text())
        result=cell(request,ChartPolicy(**json.loads(a.cell)),retained)
        checkpoint(a.output,result);print('REPRESENTATION_CELL|PASS',flush=True);return
    if bool(a.request)==bool(a.verify):p.error('choose --request or --verify')
    if a.verify:
        old=json.loads(a.verify.read_text());root=a.verify.parent;records=old['records'];request=old['request']
        policies=[r['policy'] for r in records]
    else:
        request=json.loads(a.request.read_text());root=a.output.parent
        policies=[asdict(ChartPolicy(model_normalization=n,enumeration_backend=e,chart_metric_kind=k,chart_metric_weight=w))
                  for n in ('raw','minimal') for e in ('gmp-pointed-sieve','pari-hyperellratpoints')
                  for k,w in [('raw','1'),('metric','16')]]
    out=[]
    for index,policy in enumerate(policies):
        result_path=a.output.parent/('replay-cells' if a.verify else 'cells')/(digest(policy)+'.json')
        args=[a.sage_python,str(Path(__file__).resolve()),'--cell',json.dumps(policy),'--output',str(result_path)]
        if a.verify:
            row=records[index]
            if row['outcome']!='completed':out.append(row);continue
            args+=['--verify',str(root/row['result'])]
        else:args+=['--request',str(a.request)]
        receipt=run(args,limits=Limits(request['cell_wall_seconds'],request['rss_bytes']),
                    log_path=result_path.with_suffix('.log'),result_path=result_path,
                    checkpoint_path=result_path.with_suffix('.supervisor.json'))
        row={'policy':policy,'outcome':receipt['outcome'],'result':str(result_path.relative_to(a.output.parent)),
             'supervisor':receipt}
        if receipt['outcome']=='completed':
            value=json.loads(result_path.read_text());row.update({'measurements':value['measurements'],'certified_rank_gain':value['certified_rank_gain']})
        out.append(row)
        checkpoint(a.output,{'schema':'elliptic-curves.representation-benchmark.v1','request':request,'records':out,
                             'claim_boundary':'Each stage is timed separately. Censoring or bounded misses do not prove rank bounds or optimal policies.'})
        print(f"REPRESENTATION|cell={index+1}/{len(policies)}|outcome={receipt['outcome']}",flush=True)

if __name__=='__main__':main()
