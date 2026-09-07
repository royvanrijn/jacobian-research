#!/usr/bin/env sage-python
"""Budget-only continuation of the sealed cascade; identical selector."""
import argparse,json,time
from fractions import Fraction as F
from pathlib import Path
from importlib.machinery import SourceFileLoader
CAS=Path(__file__).resolve().parent
b=SourceFileLoader('cascade_base',str(CAS/'adaptive_visibility_cascade.sage')).load_module()
D=b.D/'continuation-v1'
from research_runtime.store import checkpoint

def freeze():
    if D.exists():raise FileExistsError('preserve continuation')
    p=b.protocol();D.mkdir()
    checkpoint(D/'protocol.json',{'schema':'cascade-budget-continuation.v1','base_protocol_sha256':b.sha(b.D/'protocol.json'),
       'sources':{**b.sources(),b.rel(Path(__file__)):b.sha(Path(__file__))},'maximum_total_waves':12,
       'amendment':'Only the four-wave budget is extended to twelve. Same selector, chart and height budgets, no-gain stop, and score-stratified roster. The completed initial prefix is retained verbatim.',
       'calibration_observation':'Initial calibration run reached rank26 after four gaining waves. No residual identities or future point coordinates are inputs.'})

def run(index):
    b.guard();p=b.protocol();cp=b.read(D/'protocol.json')
    if cp['sources']!={**b.sources(),b.rel(Path(__file__)):b.sha(Path(__file__))} or cp['base_protocol_sha256']!=b.sha(b.D/'protocol.json'):
        raise ArithmeticError('continuation freeze changed')
    row=b.read(b.D/'roster.json')[index];prior=b.D/row['id'];out=D/row['id'];out.mkdir(exist_ok=True)
    if (out/'terminal.json').exists():return
    initial=b.read(prior/'terminal.json');stages=initial['stages'];prefix=list(stages)
    if stages[-1]['after']==stages[-1]['before']:
        checkpoint(out/'terminal.json',{'status':'PRIOR_NO_GAIN_STOP_RETAINED','id':row['id'],'prefix_sha256':b.sha(prior/'terminal.json'),'stages':stages,'final_rank_lower_bound':initial['final_rank_lower_bound']});return
    from research_runtime.memory_store import MemoryFactStore
    from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as Cache
    from research_runtime.search_state import raw_state
    from research_runtime.pointed_orbit_compression import compress
    from pointed_quartic_search import PointedQuarticSearch
    import pari_pointed_backend as backend
    import audit_recorded_point_mod2_rank_v3 as mod2
    import audit_retained_cloud_modl as modl
    seen=set()
    for s in prefix:
        for kind in ('mod2','modl'):
            if b.sha(b.ROOT/s[kind])!=s[kind+'_sha256']:raise ArithmeticError('prefix certificate changed')
        mod2.check(b.ROOT/s['mod2']);modl.check(b.ROOT/s['modl'])
        for c in b.read(prior/f"wave-{s['wave']:02d}"/'selection.json')['centres']:
            seen.add((c['lane'],c['orbit'] if c['lane']=='canonical' else c['parity']))
    cloud=b.read(b.ROOT/stages[-1]['mod2']);model=tuple(map(F,cloud['curve']));basis=tuple(tuple(map(F,x)) for x in cloud['independent_points'])
    cache=Cache(MemoryFactStore());state=raw_state(model,basis,cache=cache,prime_bound=1000)
    mapper=b.load('factor_free_pari_mapping.sage');mapper.pari.allocatemem(256000000,silent=True)
    for wave in range(len(stages),cp['maximum_total_waves']):
        wd=out/f'wave-{wave:02d}';wd.mkdir(exist_ok=True)
        if (wd/'result.json').exists():raise FileExistsError('preserve continuation wave')
        before=state.rank;basis=tuple(tuple(map(F,x)) for x in state.basis);centres,metric=b.candidates(model,basis,seen)
        checkpoint(wd/'selection.json',{'protocol_sha256':b.sha(D/'protocol.json'),'rank':before,'basis':[list(map(str,x)) for x in basis],'centres':centres,'metric':metric})
        charts=[];started=time.monotonic()
        for j,c in enumerate(centres):
            mapping=mapper.mapping(model,basis,c)
            search=PointedQuarticSearch(state=state,centre={'coefficients':c['representative']+[0]*(state.rank-before)},coordinate_policy=mapping['coordinate_policy'])
            transcript,points=backend.execute(search,mapping,p['height'],p['seconds_per_chart'],p['gp_sha256'])
            if backend.replay(search,mapping,transcript)!=points:raise ArithmeticError('exact map replay differs')
            compression=compress(model,state.basis,c['representative']+[0]*(state.rank-before),points)
            charts.append({'index':j,'centre':c,'mapping':mapping,'search':transcript,'admission_compression':compression})
            for k in compression['kept_indices']:state=state.adjoin(points[k],cache=cache)
            seen.add((c['lane'],c['orbit'] if c['lane']=='canonical' else c['parity']))
            checkpoint(wd/f'chart-{j:03d}.json',charts[-1])
            if (j+1)%20==0:print('CONTINUATION',row['id'],wave,j+1,state.rank,flush=True)
        checkpoint(wd/'result.json',{'status':'COMPLETE_DECLARED_WAVE','family':'det1092-cascade','parameter':row['parameter'],'curve':list(map(str,model)),'charts':charts,'final_state':state.record(),'rank_lower_bound':state.rank})
        mod2.build(wd/'result.json',wd/'mod2.json',1000,b.sha(wd/'result.json'));mod2.check(wd/'mod2.json')
        modl.build(wd/'mod2.json',wd/'modl.json');modl.check(wd/'modl.json')
        cloud=b.read(wd/'mod2.json');basis=tuple(tuple(map(F,x)) for x in cloud['independent_points']);state=raw_state(model,basis,cache=cache,prime_bound=1000)
        if state.rank!=cloud['rank_lower_bound']:raise ArithmeticError('continuation seed rank mismatch')
        stage={'wave':wave,'before':before,'after':state.rank,'charts':len(charts),'completed':sum(c['search']['status']=='bounded_search_complete' for c in charts),'seconds':time.monotonic()-started,
          'mod2':b.rel(wd/'mod2.json'),'mod2_sha256':b.sha(wd/'mod2.json'),'modl':b.rel(wd/'modl.json'),'modl_sha256':b.sha(wd/'modl.json')}
        stages.append(stage);checkpoint(out/'stages.json',stages);print('CONTINUATION STAGE',row['id'],before,state.rank,flush=True)
        if state.rank==before or state.rank>=32:break
    checkpoint(out/'terminal.json',{'status':'BOUNDED_CONTINUATION_TERMINAL','id':row['id'],'prefix_sha256':b.sha(prior/'terminal.json'),'protocol_sha256':b.sha(D/'protocol.json'),'stages':stages,'final_rank_lower_bound':state.rank,'read_paths':sorted(b.READS)})

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['freeze','run']);ap.add_argument('--index',type=int,default=0);a=ap.parse_args()
    freeze() if a.stage=='freeze' else run(a.index)
