"""Sage worker for the determinant-1092 V4 wide-bootstrap experiment.

V4 searches exactly 512 fresh M17 parity classes per fibre. It excludes all 82
classes searched by the completed V3 pilot. It stops on a certified 18th point;
that point is exported for a later unchanged-V3 cascade rather than embedding a
second cascade implementation here.
"""
from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import shutil
import sys
import tempfile

import det1092_v4_bootstrap_contract as c
from v3_warm_support import (atomic, assert_basis, bindings, check_chart, curve_tuple,
    indexed_paths, point_tuple, read, require, sha)


def copy_immutable(source, target):
    source,target=Path(source),Path(target); target.parent.mkdir(parents=True,exist_ok=True)
    if target.exists(): require(sha(source)==sha(target),'immutable copy differs: '+str(target)); return
    fd,name=tempfile.mkstemp(prefix='.copy-',dir=target.parent)
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

    _,rows=c.validate_v3_panel()
    require(not (c.D/'roster.json').exists(),'V4 roster already prepared')
    c.D.mkdir(parents=True,exist_ok=True); impl=c.own_sources(); jobs={}; cases=[]
    top_inputs={str((c.V3/'roster.json').relative_to(c.ROOT)):sha(c.V3/'roster.json')}
    for row in rows:
        case=row['id']; src=c.V3/case; folder=c.D/case; folder.mkdir(exist_ok=True)
        for name in ('seed-input.json','seed-proof.json','orbits.tsv'):
            copy_immutable(src/name,folder/name)
        copy_immutable(src/'replay-M17/epoch-00/selection.json',folder/'prior-v3-selection.json')
        copy_immutable(src/'trial-verified.json',folder/'prior-v3-verified.json')
        seed,proof=read(folder/'seed-input.json'),read(folder/'seed-proof.json')
        model,points=curve_tuple(seed['curve']),point_tuple(seed['points'])
        state=certified_state(model,points,proof); assert_basis(state,model,points,17)
        prior=read(src/'protocol.json')
        numerical=load('v4_prepare_'+case.replace('-','_'),c.CAS/'adaptive_visibility_cascade_v3.sage')
        require(numerical.sources()==prior['sources'],case+': V3 numerical sources changed')
        policy=c.policy_from_v3(case)
        policy.update(sources=numerical.sources(),implementation_sources=impl,
            software={'sage':sage.version.version,'pari':str(pari.version()),'python':sys.version},
            inputs={str(p.relative_to(c.ROOT)):sha(p) for p in
                    (folder/'seed-input.json',folder/'seed-proof.json',folder/'orbits.tsv',
                     folder/'prior-v3-selection.json',folder/'prior-v3-verified.json')},
            resource_limits=c.RESOURCE,prior_v3_protocol_sha256=sha(src/'protocol.json'),
            prior_v3_terminal_sha256=row['prior_terminal_sha256'])
        atomic(folder/'protocol.json',policy,immutable=True); jobs[case]=sha(folder/'protocol.json'); cases.append(row)
        for p in (src/'trial-verified.json',src/'replay-M17/terminal.json',src/'replay-M17/epoch-00/selection.json'):
            top_inputs[str(p.relative_to(c.ROOT))]=sha(p)
        print('V4_PREPARED',case,flush=True)
    atomic(c.D/'roster.json',{'status':'READY_V4_EIGHT','cases':cases,'jobs':jobs,'inputs':top_inputs,
           'implementation_sources':impl,'claim_boundary':c.CLAIM},immutable=True)
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
    require(p['software']['sage']==sage.version.version and p['software']['pari']==str(pari.version()),'Sage/PARI changed')
    require(sha(Path('/usr/bin/gp'))==p['gp_sha256'],'GP binary changed')
    seed,proof=read(folder/'seed-input.json'),read(folder/'seed-proof.json')
    model,points=curve_tuple(seed['curve']),point_tuple(seed['points'])
    state=certified_state(model,points,proof); assert_basis(state,model,points,17)
    engine=load('det1092_v4_'+case.replace('-','_'),c.CAS/'adaptive_visibility_cascade_v3.sage')
    engine.D,engine.ORBITS=folder,folder/'orbits.tsv'; engine.v1.D,engine.v1.ORBITS=folder,engine.ORBITS
    require(engine.sources()==p['sources'],'frozen V3 numerical engine changed')
    sources={**p['sources'],**p['implementation_sources']}; engine.v1.guard()
    return Context(case,folder,c.ROOT,p,engine,model,state,sources)


def prior_sets(ctx):
    selection=read(ctx.folder/'prior-v3-selection.json')
    masks={c.parity_mask(row['representative']) for row in selection['centres']}
    points={tuple(row['point']) for row in selection['centres']}
    require(len(masks)==len(points)==82,'prior V3 exclusion set changed')
    return masks,points


def even_positions(order,count):
    if len(order)<=count:return list(order)
    return [order[min(len(order)-1,(2*j+1)*len(order)//(2*count))] for j in range(count)]


def bootstrap_selection(ctx,publish=False):
    """Score the complete nonzero M17 quotient, then select exactly 512 fresh classes."""
    import numpy as np
    from sage.all import ZZ,matrix
    from visibility_selection_v3 import chart_profile
    require(ctx.state.rank==17,'V4 selector is M17-only')
    prior_masks,prior_points=prior_sets(ctx)
    geo=ctx.engine.load('prospective_half_lattice_v3.sage')
    hg,asym=geo.canonical_height_gram(ctx.model,point_tuple(ctx.state.basis))
    g=matrix(ZZ,geo.rounded_gram(hg,1000000)); require(g.is_positive_definite(),'nonpositive V4 metric')
    rows=[]
    with ctx.engine.ORBITS.open() as stream:
        for raw in csv.DictReader(stream,delimiter='\t'):
            shell=int(raw['minimum_norm'])
            if shell not in c.SHELLS: continue
            word=list(map(int,raw['parent_MW17_w'].split())); require(len(word)==17,'orbit dimension changed')
            mask=c.parity_mask(word)
            if mask and mask not in prior_masks:
                rows.append({'orbit':int(raw['orbit_mask']),'shell':shell,'representative':word,'parity':mask,'lanes':[]})
    require(rows and len({r['parity'] for r in rows})==len(rows),'V4 atlas empty/duplicated')
    words=np.asarray([r['representative'] for r in rows],dtype=np.int64); gg=np.asarray(g.rows(),dtype=np.int64)
    require(17**2*int(abs(words).max())**2*int(abs(gg).max())<2**62,'V4 norm overflow')
    norms=np.einsum('ij,jk,ik->i',words,gg,words,optimize=True)
    for r,n in zip(rows,norms): r['metric_norm']=int(n)
    chosen=set()
    def add(indices,lane):
        for i in indices:
            chosen.add(i)
            if lane not in rows[i]['lanes']: rows[i]['lanes'].append(lane)
    add(sorted(range(len(rows)),key=lambda i:(-rows[i]['metric_norm'],rows[i]['parity']))[:64],'deep64')
    add(sorted(range(len(rows)),key=lambda i:(rows[i]['metric_norm'],rows[i]['parity']))[:64],'shallow64')
    order=sorted(range(len(rows)),key=lambda i:(rows[i]['metric_norm'],rows[i]['parity']))
    add(even_positions(order,128),'global-quantile128')
    for shell in c.SHELLS:
        order=sorted((i for i,r in enumerate(rows) if r['shell']==shell),key=lambda i:(rows[i]['metric_norm'],rows[i]['parity']))
        require(order,'missing shell '+str(shell)); add(even_positions(order,min(32,len(order))),'shell-'+str(shell)+'-quantile')
    hashed=sorted(range(len(rows)),key=lambda i:(c.hash_key(ctx.case,rows[i]['parity']),rows[i]['parity']))
    add(hashed[:96],'sha96')
    for i in hashed:
        if len(chosen)>=c.FRESH_CHARTS: break
        add([i],'sha-fill')
    require(len(chosen)==c.FRESH_CHARTS,'V4 selector did not end at 512')
    mapper=ctx.engine.load('factor_free_pari_mapping.sage'); mapper.pari.allocatemem(256000000,silent=True)
    centres=[]
    for i in chosen:
        r=rows[i]; key,word=ctx.engine.point_key(ctx.model,point_tuple(ctx.state.basis),r['representative'])
        require(c.parity_mask(word)==r['parity'],'point sign changed parity')
        centre={**r,'representative':word,'point':list(map(str,key)),'lane':'v4-wide-atlas','lanes':sorted(r['lanes'])}
        require(tuple(centre['point']) not in prior_points,'V4 repeated a V3 centre')
        centre.update(chart_profile(mapper.mapping(ctx.model,point_tuple(ctx.state.basis),centre))); centres.append(centre)
    require(len({r['parity'] for r in centres})==len({tuple(r['point']) for r in centres})==c.FRESH_CHARTS,'duplicate V4 centres')
    centres.sort(key=lambda r:(r['quartic_bits'],r['quartic_max_bits'],c.hash_key(ctx.case,r['parity']),r['parity']))
    result={'schema':'det1092-v4-wide-bootstrap-selection.v1','rank':17,
            'basis':[list(map(str,p)) for p in ctx.state.basis],
            'rounded_gram':[list(map(int,row)) for row in g.rows()],'height_asymmetry':str(asym),
            'shells':list(c.SHELLS),'prior_v3_parities':sorted(prior_masks),'eligible_parity_classes':len(rows),
            'selected_count':len(centres),'selector':{'deep':64,'shallow':64,'global_quantiles':128,
            'shell_quantiles_each':32,'sha':96,'fill':'SHA to exactly 512 after union/dedup',
            'order':'quartic_bits, quartic_max_bits, SHA(domain:case:parity), parity'},
            'centres':centres,'claim_boundary':c.CLAIM}
    if publish: atomic(ctx.folder/'bootstrap/selection.json',result,immutable=True)
    return result


def preflight(ctx):
    from pointed_quartic_search import PointedQuarticSearch
    import pari_pointed_backend as backend
    selection=bootstrap_selection(ctx,False); centre=selection['centres'][0]
    mapper=ctx.engine.load('factor_free_pari_mapping.sage'); mapper.pari.allocatemem(256000000,silent=True)
    mapping=mapper.mapping(ctx.model,point_tuple(ctx.state.basis),centre)
    search=PointedQuarticSearch(state=ctx.state,centre={'coefficients':centre['representative']},coordinate_policy=mapping['coordinate_policy'])
    backend.validate_map(search,mapping)
    print('V4_PREFLIGHT_PASS',ctx.case,'eligible',selection['eligible_parity_classes'],'fresh',512,flush=True)


def snapshot(ctx,path,charts):
    atomic(path,{'status':'INCREMENTAL_RETAINED_CLOUD','family':'det1092-v4-bootstrap','parameter':ctx.case,
           'curve':list(map(str,ctx.model)),'charts':charts,'final_state':ctx.state.record(),'rank_lower_bound':17},immutable=True)


def returned_new(chart,seen):
    from fractions import Fraction as F
    novel=False
    for raw in chart['search']['finite_curve_points']:
        p=F(raw['x']),F(raw['y']); key=p[0],abs(p[1])
        if key not in seen: seen.add(key); novel=True
    return novel


def run_search(ctx):
    from pointed_quartic_search import PointedQuarticSearch
    from v3_warm_engine import certified_state,_audit,_odd
    import pari_pointed_backend as backend
    out=ctx.folder/'bootstrap'; out.mkdir(exist_ok=True)
    if (out/'terminal.json').exists(): print('V4_SEALED_SEARCH_REUSED',ctx.case,flush=True); return
    selection=bootstrap_selection(ctx,publish=not (out/'selection.json').exists())
    require(read(out/'selection.json')==selection,'saved V4 selection differs')
    existing=indexed_paths(out); require(len(existing)<=c.FRESH_CHARTS,'too many bootstrap charts')
    mapper=ctx.engine.load('factor_free_pari_mapping.sage'); mapper.pari.allocatemem(256000000,silent=True)
    charts=[]; seen={(x,abs(y)) for x,y in point_tuple(ctx.state.basis)}; censored=0; rank=17
    for j,centre in enumerate(selection['centres']):
        path=existing[j] if j<len(existing) else out/f'chart-{j:03d}.json'
        if path.exists():
            chart=read(path); check_chart(chart,selection,j)
            search=PointedQuarticSearch(state=ctx.state,centre={'coefficients':centre['representative']},coordinate_policy=chart['mapping']['coordinate_policy'])
            backend.replay(search,chart['mapping'],chart['search'])
        else:
            mapping=mapper.mapping(ctx.model,point_tuple(ctx.state.basis),centre)
            search=PointedQuarticSearch(state=ctx.state,centre={'coefficients':centre['representative']},coordinate_policy=mapping['coordinate_policy'])
            transcript,points=backend.execute(search,mapping,ctx.policy['height'],ctx.policy['seconds_per_chart'],ctx.policy['gp_sha256'])
            require(backend.replay(search,mapping,transcript)==points,'V4 witness replay changed')
            chart={'index':j,'centre':centre,'mapping':mapping,'search':transcript}; atomic(path,chart,immutable=True)
        charts.append(chart); censored += chart['search']['status']!='bounded_search_complete'
        if returned_new(chart,seen):
            sp=out/f'candidate-cloud-{j:04d}.json'; ap=out/f'candidate-mod2-{j:04d}.json'
            if not sp.exists(): snapshot(ctx,sp,charts)
            cloud=_audit(sp,ap)
            if cloud['rank_lower_bound']>17:
                rank=cloud['rank_lower_bound']; print('V4_BOOTSTRAP_GAIN',ctx.case,'chart',j+1,'rank',rank,flush=True); break
        if (j+1)%16==0: print(ctx.case,'V4 chart',j+1,'/',c.FRESH_CHARTS,'certified',rank,flush=True)
    final_cloud=out/'final-cloud.json'; final_mod2=out/'final-mod2.json'
    if not final_cloud.exists(): snapshot(ctx,final_cloud,charts)
    cloud=_audit(final_cloud,final_mod2); require(cloud['rank_lower_bound']==rank,'incremental/final V4 rank differs')
    if rank>17:
        state=certified_state(ctx.model,point_tuple(cloud['independent_points']),cloud['rank_certificate']); rank=state.rank
    odd=_odd(final_mod2,out/'modl.json'); require(odd=={'3':rank,'5':rank},'V4 odd-prime ranks disagree')
    if rank>17: reason='CERTIFIED_BOOTSTRAP_GAIN'
    elif censored: reason='CENSORED_BOOTSTRAP'
    else:
        require(len(charts)==c.FRESH_CHARTS,'clean no-gain did not finish 512 fresh charts')
        reason='COMPLETE_FRESH_BOOTSTRAP_NO_GAIN'
    terminal={'status':'V4_BOOTSTRAP_TERMINAL','case':ctx.case,'initial_rank':17,'rank_lower_bound':rank,
              'charts':len(charts),'censored_charts':censored,'stop_reason':reason,
              'selection_sha256':sha(out/'selection.json'),'final_audit':final_mod2.name,
              'final_audit_sha256':sha(final_mod2),'scope':c.CLAIM}
    atomic(out/'terminal.json',terminal,immutable=True)
    if rank>17:
        atomic(out/'bootstrap-seed.json',{'curve':list(map(str,ctx.model)),'points':cloud['independent_points'],
               'rank_certificate':cloud['rank_certificate'],'rank_lower_bound':rank,
               'status':'VERIFIED_SEARCH_SEED_PENDING_INDEPENDENT_V4_REPLAY'},immutable=True)
    print('V4_SEARCH_TERMINAL',ctx.case,rank,reason,flush=True)


def main():
    p=argparse.ArgumentParser(); p.add_argument('action',choices=['prepare','preflight','search','replay']); p.add_argument('--case'); a=p.parse_args()
    require(os.environ.get('DET1092_V4_SUPERVISED')=='1','use supervised V4 controller')
    if a.action=='prepare': prepare(); return
    require(a.case is not None,'case required'); ctx=setup(a.case)
    if a.action=='preflight': preflight(ctx)
    elif a.action=='search': run_search(ctx)
    else:
        from det1092_v4_bootstrap_replay import replay_case
        replay_case(ctx)

if __name__=='__main__': main()
