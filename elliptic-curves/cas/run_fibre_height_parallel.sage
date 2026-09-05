#!/usr/bin/env sage-python
"""Execute only the frozen selected chart jobs, then certify with the runner.

Parallelism changes throughput, not the selected jobs or per-chart budgets.
The shared backend validates and reuses completed content-addressed charts.
This execution manifest supplements the frozen preparation provenance.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from importlib.machinery import SourceFileLoader
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[2]
runner=SourceFileLoader('fibre_population_parallel',str(
    ROOT/'elliptic-curves/cas/fibre_height_population.sage')).load_module()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers',type=int,default=8)
    args=parser.parse_args()
    if not 1<=args.workers<=8:
        raise ValueError('bounded executor permits one through eight workers')
    population=runner.read(runner.POPULATION)
    protocol=population['protocol']
    if population['status']!='FROZEN_BEFORE_PROSPECTIVE_SEARCH' or population['inputs']!=runner.provenance(protocol):
        raise ValueError('unfrozen population or source drift')
    jobs=[]
    for name,family in population['families'].items():
        if family['arms']!=runner.select_arms(family['rows'],protocol):
            raise ValueError('frozen arms changed')
        selected=set(family['arms']['nagao']+family['arms']['height_cost'])
        for row in family['rows']:
            if row['id'] in selected:
                for index,centre in enumerate(protocol['centres'][name]):
                    jobs.append((row,index,centre))
    jobs.sort(key=lambda j:runner.content_hash(['parallel-execution-v1',j[0]['id'],j[1]]))
    path=runner.ART/'fibre_height_execution_v1.json'
    identity=dict(schema='elliptic-curves.fibre-height-execution.v1',
        population_sha256=runner.digest(runner.POPULATION),
        executor_sha256=runner.digest(Path(__file__)),workers=args.workers,
        jobs=[dict(candidate_id=r['id'],centre_index=i) for r,i,_ in jobs],
        chart_height=protocol['search_height'],seconds_per_chart=protocol['search_seconds'])
    payload=runner.read(path) if path.exists() else {**identity,'status':'RUNNING','completed_jobs':[]}
    for key,value in identity.items():
        if payload[key]!=value:
            raise ValueError('executor checkpoint drift: '+key)
    runner.write(path,payload)
    done={tuple(r) for r in payload['completed_jobs']}

    def work(job):
        row,index,centre=job
        chart=runner.detector.PointedQuarticSearch(curve=row['search_model'],subgroup=row['subgroup'],
            centre={'coefficients':centre},coordinate_policy='metric:16')
        transcript=chart.search(protocol['search_height'],protocol['search_seconds'],
            checkpoint_dir=runner.LOCAL/'search-charts').record
        if transcript['status']!='bounded_search_complete':
            raise RuntimeError(f'incomplete chart: {row["id"]}/{index}; no automatic budget extension')
        return [row['id'],index]

    started=time.perf_counter()
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        pending=[j for j in jobs if (j[0]['id'],j[1]) not in done]
        for result in pool.map(work,pending):
            payload['completed_jobs'].append(result)
            runner.write(path,payload)
            print(f'FROZEN_CHARTS|completed={len(payload["completed_jobs"])}/{len(jobs)}',flush=True)
    payload['status']='CHARTS_COMPLETE'
    payload['batch_wall_seconds']=time.perf_counter()-started
    runner.write(path,payload)
    summary=runner.search(runner.POPULATION,runner.RESULTS)
    summary.update(population_sha256=runner.digest(runner.POPULATION),results_sha256=runner.digest(runner.RESULTS),
        execution_manifest=str(path.relative_to(ROOT)),execution_manifest_sha256=runner.digest(path))
    runner.write(runner.SUMMARY,summary)


if __name__=='__main__':
    main()
