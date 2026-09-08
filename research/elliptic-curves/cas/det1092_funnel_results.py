#!/usr/bin/env python3
"""Post-search rank harvest and bounded conductor checks; no selection feedback.

Harvest accepts completed, source-bound V3 epochs, even while the rest of a
cascade is live. Every exported lower bound is recomputed from rational points.
Full search-policy replay and mathematical point independence are separate fields.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys

import det1092_funnel as f
from v3_warm_support import atomic, bindings, curve_tuple, point_tuple, read, require, sha


def canonical(value):
    return json.loads(f.packed(value))


def compare_bound(lower, upper, benchmark):
    require(0 < lower <= upper, 'invalid conductor interval')
    if benchmark is None:
        return 'NO_RECORDED_BENCHMARK'
    require(benchmark > 0, 'invalid benchmark')
    if lower > benchmark:
        return 'PROVED_ABOVE_RECORDED_MINIMUM'
    if upper < benchmark:
        return 'PROVED_BELOW_RECORDED_MINIMUM'
    return 'UNRESOLVED_COMPARISON'


def benchmarks(catalogue, rank):
    rows = [r for r in catalogue if r['rank_lower_bound'] >= rank]
    listed = [r for r in rows if r.get('conductor') is not None]
    require(all(str(r['conductor']).isdigit() and int(r['conductor']) > 0 for r in listed),
            'invalid recorded conductor')
    best = min(listed,key=lambda r:int(r['conductor'])) if listed else None
    return dict(rank_threshold=rank, id=best['id'] if best else None,
                reported_minimum=best['conductor'] if best else None,
                missing_conductor_ids=[r['id'] for r in rows if r.get('conductor') is None],
                scope='Frozen catalogue metadata; not a proof of a universal record.')


def harvest(run, out):
    import certify_compact_r17_candidates as cert
    from det1092_v3_contract import specialize
    from refresh_icarm_local_database import load_catalogue, load_inventory
    p = read(run/'protocol.json'); bindings(f.ROOT,p['sources']); bindings(f.ROOT,p['inputs'])
    out.mkdir(parents=True,exist_ok=True)
    catalogue_path = out/'comparison-snapshot.json'
    if not catalogue_path.exists():
        catalogue = load_catalogue()
        snapshot = dict(public=[{k:r[k] for k in ('id','ainvs','rank_lower_bound','conductor')} for r in catalogue['curves']],
                        inventory=[{'id':r['id'],'curve':r['curve'],'rank_lower_bound':r['rank_lower_bound']} for r in load_inventory()],
                        source_manifest_sha256=sha(f.ROOT/'elliptic-curves/data/icarm_current.json'))
        atomic(catalogue_path,snapshot,immutable=True)
    snapshot = read(catalogue_path); public = snapshot['public']
    by_j = {}
    for kind,values,model_key in [('public',public,'ainvs'),('inventory',snapshot['inventory'],'curve')]:
        for r in values:
            inv=cert.weierstrass_invariants(tuple(map(F,r[model_key])))
            by_j.setdefault(inv['c4']**3/inv['discriminant'],[]).append((kind,r['id'],r[model_key]))
    reports=[]
    for folder in sorted((run/'amplifiers').glob('*')):
        if not folder.is_dir():
            continue
        row=read(run/'seeds'/folder.name/'intake.json')
        model,base=specialize(read(f.PARENT),row)
        verified_path=folder/'trial-verified.json'
        policy_replayed=False
        if verified_path.exists():
            verified=read(verified_path);bindings(f.ROOT,verified['bindings'])
            require(verified['status']=='PASS_INDEPENDENT_DET1092_REPLAY','failed policy replay')
            policy_replayed=True
        for stage_path in sorted((folder/'replay-M17').glob('epoch-*/stage.json')):
            stage=read(stage_path)
            if stage['after'] < 23:
                continue
            audit_path=stage_path.parent/stage['audit']; audit=read(audit_path)
            require(sha(audit_path)==stage['audit_sha256'],'epoch/audit binding differs')
            require(audit['status']=='COMPLETE_DECLARED_FINITE_AUDIT','incomplete point certificate')
            cloud=f.ROOT/audit['input_path']
            require(cloud.resolve().is_relative_to(folder.resolve()) and sha(cloud)==audit['input_sha256'],
                    'unbound point cloud')
            points=point_tuple(audit['independent_points'])
            require(points[:17]==base and curve_tuple(audit['curve'])==model,'specialization prefix or curve differs')
            require(len(points)==stage['after'],'epoch rank differs from point list')
            proof=audit['rank_certificate']
            exact=canonical(cert.checked_rank(model,points,[r['prime'] for r in proof['signatures']],proof['no_rational_2_torsion_prime']))
            rank=exact['rank_lower_bound']
            inv=cert.weierstrass_invariants(model)
            matches={'public':[],'inventory':[]}
            require(inv['c4'] and inv['c6'],'special-j novelty requires a separate exact isomorphism test')
            for kind,identifier,other in by_j.get(inv['c4']**3/inv['discriminant'],[]):
                if cert.isomorphic(model,other):
                    matches[kind].append(identifier)
            result=dict(status='REPLAY_CERTIFIED_SUBGROUP',case=folder.name,epoch=stage['epoch'],
                        rank_lower_bound=rank,parameter=row['parameter'],curve=list(map(str,model)),
                        points=[list(map(str,q)) for q in points],rank_certificate=exact,
                        full_search_policy_replayed=policy_replayed,
                        catalogue_matches=matches,comparison_snapshot_sha256=sha(catalogue_path),
                        benchmarks={str(r):benchmarks(public,r) for r in range(23,rank+1)},
                        bindings={str(q.relative_to(f.ROOT)):sha(q) for q in [stage_path,audit_path,cloud,run/'protocol.json']},
                        checker_sha256=sha(Path(__file__)),
                        boundary='Exact subgroup lower bound on the retained determinant1092 fibre. Catalogue nonmatch is snapshot-relative, not literature-wide novelty. No conductor or exact rank follows from point count.')
            target=out/folder.name/f'rank-{rank:02d}-epoch-{stage["epoch"]:02d}.json'
            # Search-policy completion can change later; preserve each certificate version.
            target=target.with_name(target.stem+('-policy-replayed' if policy_replayed else '-subgroup')+'.json')
            atomic(target,result,immutable=True)
            reports.append(dict(case=folder.name,rank_lower_bound=rank,certificate=str(target.relative_to(f.ROOT)),sha256=sha(target)))
    atomic(out/'summary.json',dict(status='HARVEST_SNAPSHOT',certificates=reports,
                                  highest_certified_rank=max((r['rank_lower_bound'] for r in reports),default=None),
                                  comparison_snapshot_sha256=sha(catalogue_path)))
    return reports


def local_bounds(model, prime_bound, progress=None):
    from sage.all import QQ, ZZ, EllipticCurve, pari, prime_range
    require(type(prime_bound) is int and 3 <= prime_bound <= 10000,'local prime bound outside protocol')
    model=curve_tuple(model)
    require(all(q.denominator==1 for q in model),'integral short equation required')
    E=EllipticCurve(QQ,[QQ(str(q)) for q in model]); ep=pari.ellinit(E.a_invariants())
    delta=ZZ(E.discriminant());require(delta and ep.disc()==delta,'independent discriminants disagree')
    remainder=abs(delta);divisor=ZZ(1);rows=[]
    for prime in prime_range(prime_bound+1):
        exponent=0
        while remainder % prime==0:
            remainder//=prime;exponent+=1
        if not exponent and prime not in (2,3):
            continue
        # Over Q every prime ideal is principal. This requests its generator
        # and avoids Sage's number-field-only negative-uniformizer branch on
        # a nonminimal rational equation; it performs no global factorization.
        local=E.local_data(prime,algorithm='generic',proof=True,globally=True)
        fp=int(local.conductor_valuation());independent=ep.elllocalred(prime)
        require(fp==int(independent[0]),'Sage/PARI local conductor disagreement')
        divisor*=prime**fp
        rows.append(dict(prime=int(prime),input_discriminant_exponent=exponent,
                         minimal_discriminant_exponent=int(local.discriminant_valuation()),conductor_exponent=fp))
        if progress:
            progress(dict(local_data=rows,unresolved_cofactor=str(remainder)))
    require(remainder*ZZ.prod(ZZ(r['prime'])**r['input_discriminant_exponent'] for r in rows)==abs(delta),'discriminant reconstruction differs')
    return dict(status='EXACT' if remainder==1 else 'UNKNOWN',curve=list(map(str,model)),discriminant=str(delta),
                local_data=rows,unresolved_cofactor=str(remainder),conductor_divisor=str(divisor),
                conductor_upper_bound=str(divisor*remainder),exact_conductor=str(divisor) if remainder==1 else None,
                argument='The verified local product divides N. The remaining cofactor is prime to2,3 and bounds remaining conductor exponents by the integral discriminant. No squarefreeness or primality assumption is made.')


def conductor(packet_path,output,check=False):
    import certify_compact_r17_candidates as cert
    from sage.all import pari
    import sage.version
    packet=read(packet_path);bindings(f.ROOT,packet['bindings'])
    require(packet['status']=='REPLAY_CERTIFIED_SUBGROUP' and packet['rank_lower_bound']>=23,'conductor gate requires certified rank above22')
    model=curve_tuple(packet['curve']);points=point_tuple(packet['points']);proof=packet['rank_certificate']
    exact=canonical(cert.checked_rank(model,points,[r['prime'] for r in proof['signatures']],proof['no_rational_2_torsion_prime']))
    require(exact==proof and len(points)==packet['rank_lower_bound'],'rank certificate differs')
    result=local_bounds(model,10000,lambda row:atomic(output.with_suffix('.progress.json'),row))
    result.update(packet_sha256=sha(packet_path),rank_lower_bound=len(points),checker_sha256=sha(Path(__file__)),
                  software=dict(sage=sage.version.version,pari=str(pari.version())),comparisons={})
    for rank,record in packet['benchmarks'].items():
        benchmark=int(record['reported_minimum']) if record['reported_minimum'] else None
        result['comparisons'][rank]=dict(benchmark=record,status=compare_bound(int(result['conductor_divisor']),int(result['conductor_upper_bound']),benchmark))
    if check:
        require(read(output)==result,'conductor replay changed')
    else:
        atomic(output,result,immutable=True)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['harvest','conductor','check-conductor','local-control'])
    p.add_argument('--run',type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--packet',type=Path)
    a=p.parse_args()
    if a.action=='harvest':
        require(a.run,'--run required');reports=harvest(a.run.resolve(),a.output.resolve());print('HARVESTED',len(reports))
    elif a.action=='local-control':
        exact=local_bounds([0,0,0,-16,16],100)
        partial=local_bounds([0,0,0,-16,16],3)
        require(exact['exact_conductor']=='37' and partial['status']=='UNKNOWN' and partial['conductor_upper_bound']=='37','conductor37 control failed')
        atomic(a.output,dict(status='PASS_NONMINIMAL_LOCAL_CONDUCTOR_CONTROL',exact=exact,partial=partial),immutable=True)
    else:
        require(a.packet,'--packet required');result=conductor(a.packet.resolve(),a.output.resolve(),a.action=='check-conductor');print(result['status'])


if __name__=='__main__':
    main()
