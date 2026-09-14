#!/usr/bin/env sage -python
"""Independent Sage quotient arithmetic, transcript replay and timing accounting."""
from collections import defaultdict
from fractions import Fraction as F
import gzip
import json
from pathlib import Path
import time

from sage.all import EllipticCurve, GF, QQ
from finite_cancellation_corpus import OUT, canonical, digest, write, short, pointkey
from finite_cancellation_features import prepare
from half_lattice_pointed_sieve import linear_combination
from importlib.machinery import SourceFileLoader
from v3_warm_engine import certified_state
from pointed_quartic_search import PointedQuarticSearch
from pari_pointed_backend import replay
from memory_rank_certificate import checked_rank

CPU=OUT/'cpu'

def sage_rank(packet):
    """Construct complete cosets E(Fp)/2E(Fp) with Sage point arithmetic."""
    curve=list(map(QQ,packet['curve']));points=[tuple(map(QQ,p)) for p in packet['points']]
    E=EllipticCurve(QQ,curve);assert all(E(p) for p in points)
    pivots={};used=[]
    for sig in packet['proof']['signatures']:
        p=sig['prime'];Fp=GF(p);Ep=EllipticCurve(Fp,curve)
        allpoints=Ep.points();doubled={2*P for P in allpoints};classes={P:0 for P in doubled};reps=[Ep(0)];dim=0
        for P in allpoints:
            if P in classes:continue
            current=list(reps)
            for i,R in enumerate(current):
                representative=P+R;code=i+(1<<dim)
                for T in doubled:classes[representative+T]=code
                reps.append(representative)
            dim+=1
        assert len(classes)==len(allpoints) and dim<=2
        columns=[]
        for x,y in points:
            P=Ep(0) if x.denominator()%p==0 or y.denominator()%p==0 else Ep(Fp(x),Fp(y))
            columns.append(classes[P])
        for j in range(dim):
            row=sum(((c>>j)&1)<<i for i,c in enumerate(columns))
            while row:
                bit=row.bit_length()-1
                if bit not in pivots:pivots[bit]=row;break
                row^=pivots[bit]
        used.append(p)
        if len(pivots)==len(points):break
    tp=packet['proof']['no_rational_2_torsion_prime'];Et=EllipticCurve(GF(tp),curve)
    assert all(2*P!=Et(0) for P in Et.points() if P!=Et(0))
    assert len(pivots)==len(points)
    return {'rank':len(pivots),'primes_used':used,'no_two_torsion_prime':tp}

def main():
    start=time.process_time();plan=json.loads((CPU/'protocol.json').read_text());ph=digest((CPU/'protocol.json').read_bytes())
    inputs=json.loads((CPU/'inputs.json').read_text());supervision=json.loads((CPU/'supervision.json').read_text())
    assert supervision['status']=='COMPLETE' and len(supervision['records'])==24
    corpus={r['id']:r for r in json.loads(gzip.decompress((OUT/'corpus.json.gz').read_bytes()))}
    mapper=SourceFileLoader('cpu_replay_mapper',str(Path(__file__).with_name('lean_factor_free_pari_mapping.sage'))).load_module()
    mapper.pari.allocatemem(256000000,silent=True);rows=[];totals=defaultdict(lambda:defaultdict(float));pair={}
    for case in inputs:
        seed=case['seed'];curve=tuple(map(F,seed['curve']));basis=[tuple(map(F,p)) for p in seed['points']]
        state=certified_state(curve,basis,seed['proof']);_,normalized=short(seed['curve'],seed['points'])
        kept={pointkey(p) for p in normalized};withheld={pointkey(p) for p in corpus[case['id']]['targets']} - kept
        assert withheld
        for arm in plan['arms']:
            folder=CPU/'arms'/case['id']/arm;supervisor=json.loads((folder/'supervisor.json').read_text());result=json.loads((folder/'result.json').read_text())
            assert supervisor['status']=='COMPLETE' and supervisor['result_sha256']==digest((folder/'result.json').read_bytes())
            assert result['protocol_sha256']==supervisor['protocol_sha256']==ph
            oracle_hits=set();chart_count=0;new_points=[];selected=[]
            for record in result['charts']:
                p=folder/f'chart-{record["index"]:03d}.json';row=json.loads(p.read_text());word=case['centres'][record['index']]
                assert row['centre']==word
                mapping=row['mapping']
                if arm=='finite_selector':
                    Q=linear_combination(curve,basis,word);prepared=prepare(seed['curve'],[list(map(str,Q))],0,mapper)
                    assert row['selection']['prepared_sha256']==digest(canonical(prepared))
                    fit=plan['fit'];base=prepared['models'][0]['features'];scores=[]
                    for M in prepared['models']:
                        z=[(M['features'][k]-base[k]-mu)/sd for k,mu,sd in zip(fit['features'],fit['means'],fit['scales'])]
                        scores.append(fit['coefficients'][0]+sum(a*b for a,b in zip(fit['coefficients'][1:],z)))
                    mi=min(range(len(scores)),key=lambda i:scores[i]);assert row['selection']['scores']==scores
                    assert mapping==prepared['models'][mi]['mapping']
                else:assert mapping==mapper.mapping(curve,basis,{'representative':word})
                search=PointedQuarticSearch(state=state,centre={'coefficients':word},coordinate_policy=mapping['coordinate_policy'])
                points=replay(search,mapping,row['search']);_,normalized=short(seed['curve'],points)
                oracle_hits|={pointkey(p) for p in normalized}&withheld
                if row.get('new_point'):new_points.append(row['new_point'])
                chart_count+=1;selected.append(row['selection']['selected'])
            assert chart_count==result['point_search_calls']
            rank_replay=None
            if result['success']:
                packet=json.loads((folder/'rank-input.json').read_text());assert packet['points'][:len(basis)]==seed['points']
                assert packet['points'][len(basis):]==new_points and len(new_points)==1
                proof=packet['proof'];actual=checked_rank(curve,[tuple(map(F,p)) for p in packet['points']],
                    [s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
                assert json.loads(json.dumps(actual))==proof
                rank_replay=sage_rank(packet)
            else:assert not new_points
            r={'case':case['id'],'family':case['family'],'arm':arm,'initial_rank':len(basis),'rank_lower_bound':result['rank_lower_bound'],
                'success':result['success'],'literal_withheld_points_recovered':len(oracle_hits),'point_calls':chart_count,
                'charged_cpu_seconds':supervisor['charged_cpu_seconds'],'rank_replay':rank_replay,
                'result_sha256':supervisor['result_sha256'],'models_selected':dict(__import__('collections').Counter(selected))}
            rows.append(r);pair[case['id'],arm]=r
            totals[arm]['recoveries']+=r['success'];totals[arm]['cpu_seconds']+=r['charged_cpu_seconds'];totals[arm]['point_calls']+=chart_count
            totals[arm]['cases_with_literal_withheld_recovery']+=bool(oracle_hits)
    ratio=totals['finite_selector']['cpu_seconds']/totals['factor_free']['cpu_seconds']
    gate=totals['finite_selector']['recoveries']>=totals['factor_free']['recoveries'] and ratio<=.9
    pairs=[{'case':c['id'],'baseline_success':pair[c['id'],'factor_free']['success'],'selector_success':pair[c['id'],'finite_selector']['success'],
        'cpu_ratio':pair[c['id'],'finite_selector']['charged_cpu_seconds']/pair[c['id'],'factor_free']['charged_cpu_seconds']} for c in inputs]
    write(CPU/'verified.json',{'status':'PASS_TRANSCRIPTS_SELECTOR_AND_TWO_FINITE_IMPLEMENTATIONS','strict_cpu_gate_passed':gate,
        'cpu_ratio':ratio,'cpu_saving_percent':100*(1-ratio),'totals':dict(totals),'rows':rows,'pairs':pairs,
        'protocol_sha256':ph,'verifier_sha256':digest(Path(__file__).read_bytes()),'verification_cpu_seconds':time.process_time()-start,
        'boundary':'Newly executed retrospective controls on12 curves, cached V3 banks, M18-M20 inputs. Every claimed gain independently replayed using complete Sage finite-group cosets and the portable implementation. Literal withheld-point hits are counted separately. Runtime measurements and conditional examples do not prove a population speedup or a new record.'})
    print(json.dumps({'status':'PASS_TRANSCRIPTS_SELECTOR_AND_TWO_FINITE_IMPLEMENTATIONS','strict_cpu_gate_passed':gate,'ratio':ratio,'totals':dict(totals),'cpu_seconds':time.process_time()-start}),flush=True)

if __name__=='__main__':main()
