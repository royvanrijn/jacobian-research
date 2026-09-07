#!/usr/bin/env python3
"""One immutable completed chart; exact replay parity and CPU profile, no point search."""
import cProfile,pstats,sys
from pathlib import Path
import det1092_reduced_chart_point_pilot_v2 as trial
import retained_native19_trial_v3 as engine
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];D=ROOT/'artifacts/local/elliptic-curves/recorded-admission-profile-v1'
def main():
 p=cert.read(D/'protocol.json');assert cert.hashed(Path(__file__))==p['script_sha256'];trial.configure(0);folder=trial.D
 assert all(cert.hashed(ROOT/n)==h for n,h in p['inputs'].items())
 transcript=cert.read(folder/'result.json');maps=cert.read(folder/'maps.json');row=transcript['charts'][0];cache=engine.ReductionCache(engine.MemoryFactStore());seed,state=engine.initial(cache);state,archive=engine.rotate(state)
 assert archive==cert.read(ROOT/row['archive_path'])
 points=[(cert.F(P['x']),cert.F(P['y'])) for P in row['search']['finite_curve_points']];model=tuple(map(cert.F,seed['curve']));compression=engine.compress(model,state.basis,maps['rows'][0]['centre']['representative'],points);assert compression==row['admission_compression']
 profiler=cProfile.Profile();profiler.enable()
 for i in compression['kept_indices']:state=state.adjoin(points[i],cache=cache)
 profiler.disable();assert state.key==row['state_key'] and state.rank==row['rank_lower_bound'] and state.record()['state']['observations']==row['admission_observations']
 profiler.dump_stats(str(D/'admission.pstats'));stats=pstats.Stats(profiler);rows=[]
 for (file,line,name),(primitive,total,own,cumulative,callers) in stats.stats.items():rows.append(dict(file=file,line=line,function=name,primitive_calls=primitive,total_calls=total,own_seconds=own,cumulative_seconds=cumulative))
 checkpoint(D/'result.json',dict(status='PASS',raw_points=len(points),admission_points=len(compression['kept_indices']),rank_lower_bound=state.rank,total_profile_seconds=stats.total_tt,top_cumulative=sorted(rows,key=lambda r:-r['cumulative_seconds'])[:25],top_own=sorted(rows,key=lambda r:-r['own_seconds'])[:25],boundary='One completed first-chart admission replay with profiling overhead; exact state key, observations and rank match. Not an unprofiled benchmark, a new search, or an optimization claim.'))
 print('PASS exact recorded admissions',len(compression['kept_indices']),'profiled seconds',stats.total_tt,flush=True)
if __name__=='__main__':main()
