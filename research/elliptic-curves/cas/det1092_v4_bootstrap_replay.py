"""Independent replay for the determinant-1092 V4 wide-bootstrap experiment."""
from __future__ import annotations

from fractions import Fraction as F

import det1092_v4_bootstrap_contract as c
from v3_warm_support import (atomic, bindings, check_chart, curve_tuple, indexed_paths,
                             point_tuple, read, require, sha, within)


def replay_case(ctx):
    from det1092_v4_bootstrap_worker import bootstrap_selection
    from v3_warm_engine import certified_state
    from pointed_quartic_search import PointedQuarticSearch
    import pari_pointed_backend as backend
    import audit_recorded_point_mod2_rank_v3 as mod2
    import audit_retained_cloud_modl as modl

    folder,root=ctx.folder,ctx.root; out=folder/'bootstrap'
    terminal=read(out/'terminal.json')
    require(terminal.get('status')=='V4_BOOTSTRAP_TERMINAL' and terminal.get('case')==ctx.case,
            'missing/unrecognized V4 terminal')
    require(terminal['selection_sha256']==sha(out/'selection.json'),'V4 terminal selection binding changed')
    expected_selection=bootstrap_selection(ctx,False)
    require(read(out/'selection.json')==expected_selection,'independent V4 atlas selection differs')
    selection=expected_selection
    charts=indexed_paths(out,expected=terminal['charts'])
    if terminal['stop_reason']=='COMPLETE_FRESH_BOOTSTRAP_NO_GAIN':
        require(len(charts)==c.FRESH_CHARTS and terminal['censored_charts']==0,
                'clean V4 no-gain is not a complete 512-chart exposure')
    require(0<len(charts)<=c.FRESH_CHARTS,'invalid V4 chart count')
    mapper=ctx.engine.load('factor_free_pari_mapping.sage'); mapper.pari.allocatemem(256000000,silent=True)
    expected_points=list(point_tuple(ctx.state.basis)); seen={(x,abs(y)) for x,y in expected_points}
    censored=0
    for j,path in enumerate(charts):
        chart=read(path); check_chart(chart,selection,j)
        require(chart['search']['height_bound']==ctx.policy['height'] and
                chart['search']['timeout_seconds']==ctx.policy['seconds_per_chart'] and
                chart['search']['gp_binary_sha256']==ctx.policy['gp_sha256'],
                'V4 chart budget/backend differs')
        require(mapper.mapping(ctx.model,point_tuple(ctx.state.basis),chart['centre'])==chart['mapping'],
                'V4 quartic mapping differs')
        search=PointedQuarticSearch(state=ctx.state,centre={'coefficients':chart['centre']['representative']},
                                   coordinate_policy=chart['mapping']['coordinate_policy'])
        points=backend.replay(search,chart['mapping'],chart['search'])
        censored += chart['search']['status']!='bounded_search_complete'
        for p in points:
            key=p[0],abs(p[1])
            if key not in seen: seen.add(key); expected_points.append(p)
        if j%32==0 or j+1==len(charts): print('V4_REPLAY_CHART',ctx.case,j+1,'/',len(charts),flush=True)
    require(censored==terminal['censored_charts'],'V4 censored count changed')
    cloud_path,audit_path=out/'final-cloud.json',out/'final-mod2.json'
    snapshot,audit=read(cloud_path),read(audit_path)
    require(snapshot['charts']==[read(p) for p in charts] and curve_tuple(snapshot['curve'])==ctx.model and
            snapshot['final_state']==ctx.state.record(),'V4 final snapshot differs from executed seed/charts')
    require(audit.get('status')=='COMPLETE_DECLARED_FINITE_AUDIT' and
            audit['input_sha256']==sha(cloud_path) and within(root,audit['input_path'])==cloud_path,
            'V4 mod2 audit is not bound to final snapshot')
    require(curve_tuple(audit['curve'])==ctx.model and point_tuple(audit['points'])==tuple(expected_points),
            'V4 certified point cloud contains missing/unreturned points')
    mod2.check(audit_path)
    enlarged=point_tuple(audit['independent_points'])
    require(enlarged[:17]==point_tuple(ctx.state.basis),'V4 generic prefix changed')
    state=certified_state(ctx.model,enlarged,audit['rank_certificate'])
    require(state.rank==terminal['rank_lower_bound'],'V4 replay rank differs from terminal')
    odd_path=out/'modl.json'; odd=read(odd_path)
    require(odd['input_sha256']==sha(audit_path) and odd['points']==audit['points'],
            'V4 odd-prime audit bound to different point cloud')
    modl.check(odd_path)
    ranks={str(a['modulus']):a['finite_column_rank'] for a in odd['audits']}
    require(ranks=={'3':state.rank,'5':state.rank},'V4 mod3/mod5 replay disagrees')
    replay={'schema':'det1092-v4-wide-bootstrap-replay.v1','case':ctx.case,'sources':ctx.sources,
            'protocol_sha256':sha(folder/'protocol.json'),'terminal_sha256':sha(out/'terminal.json'),
            'selection_sha256':sha(out/'selection.json'),'charts_replayed':len(charts),
            'rank_lower_bound':state.rank,'odd_ranks':ranks,'claim_boundary':c.CLAIM}
    atomic(folder/'v4-replay.json',replay,immutable=True)
    final_paths=[folder/'v4-replay.json',folder/'protocol.json',folder/'seed-input.json',folder/'seed-proof.json',
                 out/'terminal.json',out/'selection.json',audit_path,odd_path]
    verified={'schema':'det1092-v4-wide-bootstrap-result.v1','status':'PASS_INDEPENDENT_DET1092_V4_REPLAY',
              'case':ctx.case,'parameter':read(folder/'seed-input.json')['parameter'],'initial_rank':17,
              'rank_lower_bound':state.rank,'gain':state.rank-17,'bootstrap_charts':len(charts),
              'censored_charts':censored,'stop_reason':terminal['stop_reason'],'sources':ctx.sources,
              'bindings':{str(p.relative_to(root)):sha(p) for p in final_paths},'claim_boundary':c.CLAIM}
    atomic(folder/'v4-verified.json',verified,immutable=True)
    if state.rank>17:
        discovery={'status':'VERIFIED_V4_BOOTSTRAP_GAIN_ELIGIBLE_FOR_FROZEN_V3_CASCADE',
                   'case':ctx.case,'parameter':verified['parameter'],'curve':list(map(str,ctx.model)),
                   'points':[list(map(str,p)) for p in state.basis],'rank_lower_bound':state.rank,
                   'rank_certificate':audit['rank_certificate'],'verified_sha256':sha(folder/'v4-verified.json'),
                   'claim_boundary':'Certified new lower bound for an existing prospective fibre; no catalogue novelty or rank-record claim.'}
        atomic(folder/'v4-discovery.json',discovery,immutable=True)
    print('INDEPENDENT_V4_REPLAY_COMPLETE',ctx.case,state.rank,flush=True)
    return verified
