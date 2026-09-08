"""Sage worker for the determinant-1092 V4 wide-bootstrap experiment.

The bootstrap stage is new: 512 fresh M17 parity classes per fibre from the
complete nonzero quotient.  All prior V3 M17 classes are excluded.  Once an
independent point is certified, the code calls the unchanged V3 landscape for
all subsequent epochs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
from pathlib import Path
import shutil
import sys
import tempfile
import time
import uuid

import det1092_v4_bootstrap_contract as c
from v3_warm_support import (atomic, assert_basis, bindings, check_chart, curve_tuple,
    indexed_paths, point_tuple, read, require, sha)


def copy_immutable(source, target):
    source, target = Path(source), Path(target); target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        require(sha(source) == sha(target), 'immutable copy differs: '+str(target)); return
    fd, name = tempfile.mkstemp(prefix='.copy-', dir=target.parent)
    try:
        with os.fdopen(fd,'wb') as out, source.open('rb') as inp:
            shutil.copyfileobj(inp,out); out.flush(); os.fsync(out.fileno())
        try: os.link(name,target)
        except FileExistsError: require(sha(source)==sha(target),'concurrent immutable copy differs')
    finally:
        try: os.unlink(name)
        except FileNotFoundError: pass


def prepare():
    from sage.all import pari
    import sage.version
    from v3_warm_engine import certified_state, load

    v3roster, rows = c.validate_v3_panel()
    require(not (c.D/'roster.json').exists(), 'V4 roster already prepared')
    c.D.mkdir(parents=True, exist_ok=True)
    impl = c.own_sources(); jobs = {}; cases = []
    top_inputs = {str((c.V3/'roster.json').relative_to(c.ROOT)):sha(c.V3/'roster.json')}
    for row in rows:
        case = row['id']; src = c.V3/case; folder = c.D/case; folder.mkdir(exist_ok=True)
        for name in ('seed-input.json','seed-proof.json','orbits.tsv'):
            copy_immutable(src/name, folder/name)
        copy_immutable(src/'replay-M17/epoch-00/selection.json', folder/'prior-v3-selection.json')
        copy_immutable(src/'trial-verified.json', folder/'prior-v3-verified.json')
        seed, proof = read(folder/'seed-input.json'), read(folder/'seed-proof.json')
        model, points = curve_tuple(seed['curve']), point_tuple(seed['points'])
        state = certified_state(model, points, proof); assert_basis(state,model,points,17)
        prior_policy = read(src/'protocol.json')
        numerical = load('v4_prepare_numeric_'+case.replace('-','_'), c.CAS/'adaptive_visibility_cascade_v3.sage')
        require(numerical.sources() == prior_policy['sources'], case+': frozen V3 numerical sources changed')
        policy = c.policy_from_v3(case)
        policy.update(
            sources=numerical.sources(), implementation_sources=impl,
            software={'sage':sage.version.version,'pari':str(pari.version()),'python':sys.version},
            inputs={str(p.relative_to(c.ROOT)):sha(p) for p in
                    (folder/'seed-input.json',folder/'seed-proof.json',folder/'orbits.tsv',
                     folder/'prior-v3-selection.json',folder/'prior-v3-verified.json')},
            resource_limits=c.RESOURCE,
            prior_v3_protocol_sha256=sha(src/'protocol.json'),
            prior_v3_terminal_sha256=row['prior_terminal_sha256'])
        atomic(folder/'protocol.json',policy,immutable=True)
        jobs[case]=sha(folder/'protocol.json'); cases.append(row)
        for p in (src/'trial-verified.json',src/'replay-M17/terminal.json',src/'replay-M17/epoch-00/selection.json'):
            top_inputs[str(p.relative_to(c.ROOT))]=sha(p)
        print('V4_PREPARED',case,'fresh',c.FRESH_CHARTS,flush=True)
    atomic(c.D/'roster.json',{'status':'READY_V4_EIGHT','cases':cases,'jobs':jobs,
           'inputs':top_inputs,'implementation_sources':impl,'claim_boundary':c.CLAIM},immutable=True)
    c.validate_roster(); print('V4_PREPARED_ALL_EIGHT',flush=True)


def setup(case):
    from sage.all import pari
    import sage.version
    from v3_warm_engine import Context, certified_state, load
    roster=c.validate_roster(); row=next((r for r in roster['cases'] if r['id']==case),None)
    require(row is not None,'case outside frozen V4 roster')
    folder=c.D/case; p=read(folder/'protocol.json')
    for key in ('inputs','sources','implementation_sources'): bindings(c.ROOT,p[key])
    require(p['implementation_sources']==c.own_sources(),'V4 implementation changed after preparation')
    require(p['software']['sage']==sage.version.version and p['software']['pari']==str(pari.version()),
            'Sage/PARI changed since V4 preparation')
    require(sha(Path('/usr/bin/gp'))==p['gp_sha256'],'GP binary changed')
    seed,proof=read(folder/'seed-input.json'),read(folder/'seed-proof.json')
    model,points=curve_tuple(seed['curve']),point_tuple(seed['points'])
    state=certified_state(model,points,proof); assert_basis(state,model,points,17)
    engine=load('det1092_v4_numeric_'+case.replace('-','_'),c.CAS/'adaptive_visibility_cascade_v3.sage')
    engine.D,engine.ORBITS=folder,folder/'orbits.tsv'; engine.v1.D,engine.v1.ORBITS=folder,engine.ORBITS
    require(engine.sources()==p['sources'],'frozen V3 engine changed')
    sources={**p['sources'],**p['implementation_sources']}
    engine.v1.guard()
    return Context(case,folder,c.ROOT,p,engine,model,state,sources)


def _prior(ctx):
    selection=read(ctx.folder/'prior-v3-selection.json')
    masks={c.parity_mask(row['representative']) for row in selection['centres']}
    points={tuple(row['point']) for row in selection['centres']}
    require(len(masks)==82 and len(points)==82,'prior V3 exclusion set changed')
    return masks,points


def _even_positions(order,count):
    if count<=0:return []
    if len(order)<=count:return list(order)
    return [order[min(len(order)-1,(2*j+1)*len(order)//(2*count))] for j in range(count)]


def bootstrap_selection(ctx, publish=False):
    """Recompute the complete nonzero M17 atlas and choose exactly 512 fresh classes."""
    import numpy as np
    from sage.all import ZZ,matrix
    from visibility_selection_v3 import chart_profile

    require(ctx.state.rank==17,'wide bootstrap is M17-only')
    prior_masks,prior_points=_prior(ctx)
    geo=ctx.engine.load('prospective_half_lattice_v3.sage')
    hg,asym=geo.canonical_height_gram(ctx.model,point_tuple(ctx.state.basis))
    g=matrix(ZZ,geo.rounded_gram(hg,1000000)); require(g.is_positive_definite(),'nonpositive V4 metric')
    rows=[]
    with ctx.engine.ORBITS.open() as stream:
        for raw in csv.DictReader(stream,delimiter='\t'):
            shell=int(raw['minimum_norm'])
            if shell not in c.SHELLS: continue
            word=list(map(int,raw['parent_MW17_w'].split()))
            require(len(word)==17,'orbit representative dimension changed')
            mask=c.parity_mask(word)
            if mask==0 or mask in prior_masks: continue
            rows.append({'orbit':int(raw['orbit_mask']),'shell':shell,'representative':word,'parity':mask})
    require(rows,'empty V4 atlas')
    require(len({r['parity'] for r in rows})==len(rows),'duplicate parity in complete quotient')
    words=np.asarray([r['representative'] for r in rows],dtype=np.int64)
    gg=np.asarray(g.rows(),dtype=np.int64)
    require(17**2*int(abs(words).max())**2*int(abs(gg).max())<2**62,'V4 norm overflow')
    norms=np.einsum('ij,jk,ik->i',words,gg,words,optimize=True)
    for r,n in zip(rows,norms): r['metric_norm']=int(n); r['lanes']=[]
    chosen=set()
    def add(indices,lane):
        for i in indices:
            chosen.add(i)
            if lane not in rows[i]['lanes']: rows[i]['lanes'].append(lane)
    deep=sorted(range(len(rows)),key=lambda i:(-rows[i]['metric_norm'],rows[i]['parity']))[:64]
    shallow=sorted(range(len(rows)),key=lambda i:(rows[i]['metric_norm'],rows[i]['parity']))[:64]
    add(deep,'deep64'); add(shallow,'shallow64')
    global_order=sorted(range(len(rows)),key=lambda i:(rows[i]['metric_norm'],rows[i]['parity']))
    add(_even_positions(global_order,128),'global-quantile128')
    for shell in c.SHELLS:
        order=sorted((i for i,r in enumerate(rows) if r['shell']==shell),
                     key=lambda i:(rows[i]['metric_norm'],rows[i]['parity']))
        require(order,'missing declared shell '+str(shell))
        add(_even_positions(order,min(32,len(order))),'shell-'+str(shell)+'-quantile')
    hashed=sorted(range(len(rows)),key=lambda i:(c.hash_key(ctx.case,rows[i]['parity']),rows[i]['parity']))
    add(hashed[:96],'sha96')
    for i in hashed:
        if len(chosen)>=c.FRESH_CHARTS: break
        add([i],'sha-fill')
    require(len(chosen)==c.FRESH_CHARTS,'V4 selector did not produce exactly 512 classes')
    mapper=ctx.engine.load('factor_free_pari_mapping.sage'); mapper.pari.allocatemem(256000000,silent=True)
    centres=[]
    for i in chosen:
        r=rows[i]; key,word=ctx.engine.point_key(ctx.model,point_tuple(ctx.state.basis),r['representative'])
        require(c.parity_mask(word)==r['parity'],'sign normalization changed parity')
        centre={**r,'representative':word,'point':list(map(str,key)),'lane':'v4-wide-atlas',
                'lanes':sorted(r['lanes'])}
        require(tuple(centre['point']) not in prior_points,'V4 repeated a prior V3 centre')
        profile=chart_profile(mapper.mapping(ctx.model,point_tuple(ctx.state.basis),centre)); centre.update(profile)
        centres.append(centre)
    require(len({r['parity'] for r in centres})==len({tuple(r['point']) for r in centres})==c.FRESH_CHARTS,
            'V4 selected duplicate centres')
    centres.sort(key=lambda r:(r['quartic_bits'],r['quartic_max_bits'],c.hash_key(ctx.case,r['parity']),r['parity']))
    result={'schema':'det1092-v4-wide-bootstrap-selection.v1','rank':17,
            'basis':[list(map(str,p)) for p in ctx.state.basis],
            'rounded_gram':[list(map(int,row)) for row in g.rows()], 'height_asymmetry':str(asym),
            'shells':list(c.SHELLS),'prior_v3_parities':sorted(prior_masks),
            'eligible_parity_classes':len(rows),'selected_count':len(centres),
            'selector':{'deep':64,'shallow':64,'global_quantiles':128,'shell_quantiles_each':32,'sha':96,
                        'fill':'SHA until exactly 512; union deduplicated before fill',
                        'order':'quartic_bits, quartic_max_bits, SHA(domain:case:parity), parity'},
            'centres':centres,'claim_boundary':c.CLAIM}
    if publish:
        atomic(ctx.folder/'bootstrap/selection.json',result,immutable=True)
    return result


def preflight(ctx):
    from pointed_quartic_search import PointedQuarticSearch
    import pari_pointed_backend as backend
    selection=bootstrap_selection(ctx,publish=False); centre=selection['centres'][0]
    mapper=ctx.engine.load('factor_free_pari_mapping.sage'); mapper.pari.allocatemem(256000000,silent=True)
    mapping=mapper.mapping(ctx.model,point_tuple(ctx.state.basis),centre)
    search=PointedQuarticSearch(state=ctx.state,centre={'coefficients':centre['representative']},
                               coordinate_policy=mapping['coordinate_policy'])
    backend.validate_map(search,mapping)
    print('V4_PREFLIGHT_PASS',ctx.case,'eligible',selection['eligible_parity_classes'],'selected',512,flush=True)


def _snapshot(ctx,path,charts,state,family):
    atomic(path,{'status':'INCREMENTAL_RETAINED_CLOUD','family':family,'parameter':ctx.case,
           'curve':list(map(str,ctx.model)),'charts':charts,'final_state':state.record(),
           'rank_lower_bound':state.rank},immutable=True)


def _audit_snapshot(path,audit):
    from v3_warm_engine import _audit
    return _audit(path,audit)


def _odd(audit,path):
    from v3_warm_engine import _odd
    return _odd(audit,path)


def _returned_new(chart,seen):
    from fractions import Fraction as F
    new=False
    for raw in chart['search']['finite_curve_points']:
        p=F(raw['x']),F(raw['y']); key=p[0],abs(p[1])
        if key not in seen: seen.add(key); new=True
    return new


def run_bootstrap(ctx):
    from pointed_quartic_search import PointedQuarticSearch
    from v3_warm_engine import certified_state, restore_state
    import pari_pointed_backend as backend
    out=ctx.folder/'bootstrap'; out.mkdir(exist_ok=True)
    if (out/'terminal.json').exists(): return read(out/'terminal.json')
    selection=bootstrap_selection(ctx,publish=not (out/'selection.json').exists())
    if (out/'selection.json').exists(): require(read(out/'selection.json')==selection,'saved V4 selection differs')
    state=ctx.state
    if (out/'seed-state.json').exists(): state=restore_state(read(out/'seed-state.json'),ctx.model,point_tuple(ctx.state.basis))
    else: atomic(out/'seed-state.json',state.record(),immutable=True)
    mapper=ctx.engine.load('factor_free_pari_mapping.sage'); mapper.pari.allocatemem(256000000,silent=True)
    existing=indexed_paths(out); require(len(existing)<=c.FRESH_CHARTS,'too many V4 bootstrap charts')
    charts=[]; seen={(x,abs(y)) for x,y in point_tuple(state.basis)}; censored=0; gain_audit=None
    for j,centre in enumerate(selection['centres']):
        path=existing[j] if j<len(existing) else out/f'chart-{j:03d}.json'
        if path.exists():
            chart=read(path); check_chart(chart,selection,j)
            search=PointedQuarticSearch(state=state,centre={'coefficients':centre['representative']},
                                       coordinate_policy=chart['mapping']['coordinate_policy'])
            backend.replay(search,chart['mapping'],chart['search'])
        else:
            mapping=mapper.mapping(ctx.model,point_tuple(state.basis),centre)
            search=PointedQuarticSearch(state=state,centre={'coefficients':centre['representative']},
                                       coordinate_policy=mapping['coordinate_policy'])
            transcript,points=backend.execute(search,mapping,ctx.policy['height'],ctx.policy['seconds_per_chart'],ctx.policy['gp_sha256'])
            require(backend.replay(search,mapping,transcript)==points,'V4 bootstrap witness replay changed')
            chart={'index':j,'centre':centre,'mapping':mapping,'search':transcript}; atomic(path,chart,immutable=True)
        charts.append(chart); censored += chart['search']['status']!='bounded_search_complete'
        novel=_returned_new(chart,seen)
        if novel:
            snap=out/f'candidate-cloud-{j:04d}.json'; audit=out/f'candidate-mod2-{j:04d}.json'
            if not snap.exists(): _snapshot(ctx,snap,charts,state,'det1092-v4-bootstrap')
            cloud=_audit_snapshot(snap,audit)
            if cloud['rank_lower_bound']>17:
                enlarged=point_tuple(cloud['independent_points']); require(enlarged[:17]==point_tuple(ctx.state.basis),'generic prefix lost')
                state=certified_state(ctx.model,enlarged,cloud['rank_certificate']); gain_audit=audit
                print('V4_BOOTSTRAP_GAIN',ctx.case,'chart',j+1,'rank',state.rank,flush=True); break
        if (j+1)%16==0: print(ctx.case,'V4 chart',j+1,'/',c.FRESH_CHARTS,'certified',state.rank,flush=True)
    final_snap=out/'final-cloud.json'; final_audit=out/'final-mod2.json'
    if not final_snap.exists(): _snapshot(ctx,final_snap,charts,ctx.state,'det1092-v4-bootstrap')
    cloud=_audit_snapshot(final_snap,final_audit)
    require(cloud['rank_lower_bound']==state.rank,'V4 final audit differs from incremental gain')
    if state.rank>17:
        state=certified_state(ctx.model,point_tuple(cloud['independent_points']),cloud['rank_certificate'])
    odd=_odd(final_audit,out/'modl.json'); require(odd=={'3':state.rank,'5':state.rank},'V4 bootstrap odd ranks disagree')
    if state.rank>17: reason='BOOTSTRAP_GAIN'
    elif censored: reason='CENSORED_BOOTSTRAP'
    else:
        require(len(charts)==c.FRESH_CHARTS,'uncensored V4 no-gain did not finish 512 charts')
        reason='COMPLETE_FRESH_BOOTSTRAP_NO_GAIN'
    terminal={'status':'V4_BOOTSTRAP_TERMINAL','case':ctx.case,'initial_rank':17,'rank_lower_bound':state.rank,
              'charts':len(charts),'censored_charts':censored,'stop_reason':reason,
              'selection_sha256':sha(out/'selection.json'),'final_audit':final_audit.name,
              'final_audit_sha256':sha(final_audit),'scope':c.CLAIM}
    atomic(out/'terminal.json',terminal,immutable=True)
    if state.rank>17:
        atomic(out/'bootstrap-seed.json',{'curve':list(map(str,ctx.model)),
               'points':cloud['independent_points'],'rank_certificate':cloud['rank_certificate']},immutable=True)
    return terminal


def _all_tested(ctx):
    _,prior_points=_prior(ctx); tested=set(prior_points)
    sel=read(ctx.folder/'bootstrap/selection.json')
    for row in sel['centres']: tested.add(tuple(row['point']))
    return tested


def run_cascade(ctx):
    """After V4 finds the first point, use unchanged V3 landscape and budgets."""
    from pointed_quartic_search import PointedQuarticSearch
    from v3_warm_engine import certified_state
    from v3_warm_replay import verify_landscape
    import pari_pointed_backend as backend
    boot=read(ctx.folder/'bootstrap/bootstrap-seed.json')
    state=certified_state(ctx.model,point_tuple(boot['points']),boot['rank_certificate'])
    require(state.rank>17,'cascade requires a certified V4 bootstrap gain')
    out=ctx.folder/'cascade'; out.mkdir(exist_ok=True)
    if (out/'terminal.json').exists(): return read(out/'terminal.json')
    policy=ctx.policy['v3_policy']; mapper=ctx.engine.load('factor_free_pari_mapping.sage'); mapper.pari.allocatemem(256000000,silent=True)
    tested=_all_tested(ctx); stages=[]; total=0; reason='EPOCH_BUDGET_EXHAUSTED'
    for epoch in range(policy['max_epochs']):
        if total>=policy['max_charts']: reason='CHART_BUDGET_EXHAUSTED'; break
        wd=out/f'epoch-{epoch:02d}'; wd.mkdir(exist_ok=True); basis=point_tuple(state.basis); before=state.rank
        if (wd/'selection.json').exists():
            verify_landscape(ctx.engine,ctx.model,basis,tested,wd,policy); selection=read(wd/'selection.json')
        else:
            ctx.engine.landscape(ctx.model,basis,tested,wd,policy); selection=read(wd/'selection.json')
        centres=selection['centres']; require(centres,'V3 cascade produced no finalists')
        existing=indexed_paths(wd); require(len(existing)<=len(centres),'extra cascade charts')
        charts=[]; seen={(x,abs(y)) for x,y in basis}; censored=0; gained=False
        for j,centre in enumerate(centres):
            if total>=policy['max_charts']: break
            path=existing[j] if j<len(existing) else wd/f'chart-{j:03d}.json'
            if path.exists():
                chart=read(path); check_chart(chart,selection,j)
                search=PointedQuarticSearch(state=state,centre={'coefficients':centre['representative']},coordinate_policy=chart['mapping']['coordinate_policy'])
                backend.replay(search,chart['mapping'],chart['search'])
            else:
                mapping=mapper.mapping(ctx.model,basis,centre)
                search=PointedQuarticSearch(state=state,centre={'coefficients':centre['representative']},coordinate_policy=mapping['coordinate_policy'])
                transcript,points=backend.execute(search,mapping,policy['height'],policy['seconds_per_chart'],policy['gp_sha256'])
                require(backend.replay(search,mapping,transcript)==points,'post-gain V3 witness replay changed')
                chart={'index':j,'centre':centre,'mapping':mapping,'search':transcript}; atomic(path,chart,immutable=True)
            charts.append(chart); tested.add(tuple(centre['point'])); total+=1; censored += chart['search']['status']!='bounded_search_complete'
            if _returned_new(chart,seen):
                snap=wd/f'candidate-cloud-{j:04d}.json'; audit=wd/f'candidate-mod2-{j:04d}.json'
                if not snap.exists(): _snapshot(ctx,snap,charts,state,'det1092-v4-postgain-v3')
                cloud=_audit_snapshot(snap,audit)
                if cloud['rank_lower_bound']>before:
                    enlarged=point_tuple(cloud['independent_points']); require(enlarged[:before]==basis,'cascade prefix lost')
                    state=certified_state(ctx.model,enlarged,cloud['rank_certificate']); gained=True; break
        snap=wd/'final-cloud.json'; audit=wd/'final-mod2.json'
        if not snap.exists(): _snapshot(ctx,snap,charts,certified_state(ctx.model,basis,read(snap)['rank_certificate']) if False else ctx.state,'unused')
        # Rewrite the final snapshot correctly if this is its first publication: the seed state is reconstructed from basis.
        if read(snap).get('family')=='unused':
            raise RuntimeError('internal V4 snapshot sentinel should be unreachable')
        cloud=_audit_snapshot(snap,audit)
        # Candidate snapshot may have found a gain before final snapshot existed. Build a canonical final snapshot using the epoch seed record.
        # Existing code path above is intentionally replaced below on first execution by final-cloud-v2.
        final2=wd/'final-cloud-v2.json'; audit2=wd/'final-mod2-v2.json'
        if not final2.exists():
            from research_runtime.memory_store import MemoryFactStore
            from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as Cache
            from research_runtime.search_state import raw_state
            epoch_state=raw_state(ctx.model,basis,cache=Cache(MemoryFactStore()),prime_bound=1000)
            _snapshot(ctx,final2,charts,epoch_state,'det1092-v4-postgain-v3')
        cloud=_audit_snapshot(final2,audit2)
        if cloud['rank_lower_bound']>before:
            state=certified_state(ctx.model,point_tuple(cloud['independent_points']),cloud['rank_certificate']); gained=True
        else:
            require(state.rank==before,'cascade incremental/final audit disagreement')
        odd=_odd(audit2,wd/'modl.json'); require(odd=={'3':state.rank,'5':state.rank},'cascade odd ranks disagree')
        if state.rank>=policy['target_rank']: reason='TARGET_LOWER_BOUND_REACHED'
        elif total>=policy['max_charts']: reason='CHART_BUDGET_EXHAUSTED'
        elif censored: reason='CENSORED_SEARCH'
        elif gained: reason='REBUILD_AFTER_CERTIFIED_GAIN'
        else: reason='COMPLETE_FINITE_NO_GAIN'
        stage={'epoch':epoch,'before':before,'after':state.rank,'charts':len(charts),'censored_charts':censored,
               'stale_charts_cancelled':len(centres)-len(charts) if gained else 0,
               'audit':audit2.name,'audit_sha256':sha(audit2),'full_cosets_scored':selection['full_cosets_scored'],
               'stop_reason':reason}
        atomic(wd/'stage.json',stage,immutable=True); stages.append(stage); atomic(out/'stages.json',stages)
        print('V4_POSTGAIN_V3',ctx.case,before,'->',state.rank,'charts',len(charts),reason,flush=True)
        if reason!='REBUILD_AFTER_CERTIFIED_GAIN': break
    terminal={'status':'V4_POSTGAIN_V3_TERMINAL','case':ctx.case,'initial_rank':read(ctx.folder/'bootstrap/terminal.json')['rank_lower_bound'],
              'rank_lower_bound':state.rank,'charts':total,'stages':stages,'stop_reason':reason,'scope':c.CLAIM}
    atomic(out/'terminal.json',terminal,immutable=True); return terminal


def run_search(ctx):
    boot=run_bootstrap(ctx)
    if boot['rank_lower_bound']>17: run_cascade(ctx)
    print('V4_SEARCH_COMPLETE',ctx.case,'bootstrap',boot['rank_lower_bound'],flush=True)


def main():
    p=argparse.ArgumentParser(); p.add_argument('action',choices=['prepare','preflight','search','replay']); p.add_argument('--case'); a=p.parse_args()
    require(os.environ.get('DET1092_V4_SUPERVISED')=='1','use the supervised V4 controller')
    if a.action=='prepare': prepare(); return
    require(a.case is not None,'case required'); ctx=setup(a.case)
    if a.action=='preflight': preflight(ctx)
    elif a.action=='search': run_search(ctx)
    else:
        from det1092_v4_bootstrap_replay import replay_case
        replay_case(ctx)

if __name__=='__main__': main()
