#!/usr/bin/env python3
"""A bounded, explicitly noncertifying triangulated-cover cue."""
import hashlib
import json
from pathlib import Path
import signal
import time
import snappy

HERE=Path(__file__).resolve().parent


def timeout(*_):raise TimeoutError('frozen geometry cue time limit')


def main():
    plan=json.loads((HERE/'compact_cover_geometry_plan.json').read_text())
    compact=json.loads((HERE/'compact_cover_results.json').read_text())
    files=['compact_cover_geometry_plan.json','probe_compact_cover_geometry.py','compact_cover_results.json','compact_cover_verification.json']
    result={'status':'RUNNING','software':{'snappy':snappy.__version__},'results':[],
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files}}
    path=HERE/'compact_cover_geometry_cues.json';signal.signal(signal.SIGALRM,timeout);started=time.monotonic()
    cases=[('4_1',snappy.Link('4_1').PD_code())]+[(r['name'],r['pd']) for r in compact['targets']]
    for name,pd in cases:
        assert time.monotonic()-started<plan['limits']['seconds_total']
        signal.alarm(plan['limits']['seconds_each'])
        try:
            exterior=snappy.Link(pd).exterior();covers=exterior.covers(2,cover_type='cyclic');assert len(covers)==1
            M=covers[0];cover_info=M.cover_info();M.dehn_fill((1,0));G=M.fundamental_group()
            row={'name':name,'status':'CUE_ONLY','pd':pd,'cover_info':cover_info,'candidate_filling':[1,0],
                'homology':str(M.homology()),'solution_type':M.solution_type(),'numerical_volume':str(M.volume()),
                'generators':G.generators(),'relators':G.relators(),'original_generators':G.original_generators(),
                'generators_in_originals':G.generators_in_originals(),'word_moves':G._word_moves(),
                'triangulation':M._to_string(),
                'assurance':'Exploratory cue: peripheral filling provenance and exact transport to core group not certified; no geometry/rank theorem.'}
        except TimeoutError as error:row={'name':name,'status':'CAPPED','reason':str(error)}
        finally:signal.alarm(0)
        result['results'].append(row);path.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps({k:row[k] for k in ['name','status','homology','solution_type','numerical_volume','generators','relators'] if k in row}),flush=True)
    result['status']='COMPLETE_CUES';path.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
