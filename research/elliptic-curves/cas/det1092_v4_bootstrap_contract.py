"""Contracts for the wide-atlas determinant-1092 V4 bootstrap experiment.

V4 changes only the M17 bootstrap exposure. It searches 512 fresh parity
classes per fibre, excluding every M17 parity already searched by the completed
V3 eight-fibre pilot. A certified gain is exported as a new independently
verified seed for a later unchanged-V3 cascade; V4 itself stops there. Null
results are bounded visibility experiments, never saturation or rank upper bounds.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

from v3_warm_support import bindings, read, require, sha
import det1092_v3_contract as v3c

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
LOCAL = ROOT/'artifacts/local/elliptic-curves'
D = LOCAL/'det1092-v4-wide-bootstrap'
AUTO = LOCAL/'det1092-v4-wide-controller'
V3 = v3c.D
DOMAIN = 'det1092-v4-wide-bootstrap-v1'
FRESH_CHARTS = 512
SHELLS = (4, 6, 8, 10, 12)
RESOURCE = {'prepare':3600, 'preflight':180, 'search':21600, 'replay':14400,
            'rss_bytes':3221225472}
CLAIM = ('Wide M17 visibility follow-up on the same eight frozen determinant-1092 fibres. '
         'Exactly 512 previously unsearched parity classes are selected per fibre from the '
         'complete nonzero degree-two quotient using equation/lattice metrics only. A finite '
         'no-gain result is not a rank upper bound, saturation proof, or proof that no jump exists. '
         'Any certified gain is independently replayed and exported only as an eligible seed for '
         'the already-frozen V3 cascade; V4 does not retune or execute that cascade.')


def own_sources():
    names = ('det1092_v4_bootstrap_contract.py','det1092_v4_bootstrap_worker.py',
             'det1092_v4_bootstrap_replay.py','run_det1092_v4_bootstrap.py',
             'v3_warm_support.py','v3_warm_engine.py','det1092_v3_contract.py',
             'visibility_selection_v3.py','research_runtime/supervisor.py')
    return {str((CAS/n).relative_to(ROOT)):sha(CAS/n) for n in names}


def parity_mask(word):
    mask = 0
    for i, value in enumerate(word):
        if int(value) & 1:
            mask |= 1 << i
    return mask


def hash_key(case, mask):
    return hashlib.sha256((DOMAIN+':'+case+':'+str(int(mask))).encode()).hexdigest()


def prior_case(case):
    """Validate the completed 82-chart V3 null and return its frozen centres."""
    folder = V3/case
    verified = read(folder/'trial-verified.json')
    require(verified.get('status') == 'PASS_INDEPENDENT_DET1092_V3_REPLAY',
            case+': missing independent V3 verification')
    require(verified.get('initial_rank') == verified.get('rank_lower_bound') == 17 and
            verified.get('gain') == 0 and verified.get('charts') == 82 and
            verified.get('stop_reason') == 'COMPLETE_FINITE_NO_GAIN',
            case+': V4 requires the clean 82-chart 17->17 V3 terminal')
    bindings(ROOT, verified['bindings']); bindings(ROOT, verified['sources'])
    selection = read(folder/'replay-M17/epoch-00/selection.json')
    require(selection['rank'] == 17 and len(selection['centres']) == 82,
            case+': unexpected V3 M17 schedule')
    masks = [parity_mask(row['representative']) for row in selection['centres']]
    require(len(set(masks)) == 82, case+': prior V3 schedule repeats a parity class')
    return verified, selection, set(masks)


def validate_v3_panel():
    roster = v3c.validate_roster()
    require(len(roster['pilot']) == 8, 'V3 pilot roster changed')
    rows = []
    for row in roster['pilot']:
        _, _, masks = prior_case(row['id'])
        rows.append(dict(row, prior_verified_sha256=sha(V3/row['id']/'trial-verified.json'),
                         prior_terminal_sha256=sha(V3/row['id']/'replay-M17/terminal.json'),
                         prior_selection_sha256=sha(V3/row['id']/'replay-M17/epoch-00/selection.json'),
                         prior_parity_count=len(masks)))
    return roster, rows


def policy_from_v3(case):
    p = dict(read(V3/case/'protocol.json'))
    require(p['initial_rank'] == 17 and p['target_rank'] == 32, 'V3 rank policy changed')
    for key in ('height','seconds_per_chart','gp_sha256'):
        require(key in p, 'missing frozen V3 policy key '+key)
    return {
        'schema':'det1092-v4-wide-bootstrap-job.v1',
        'initial_rank':17, 'generic_rank':17, 'bootstrap_target_rank':18,
        'bootstrap_fresh_charts':FRESH_CHARTS, 'bootstrap_shells':list(SHELLS),
        'bootstrap_domain':DOMAIN,
        'height':p['height'], 'seconds_per_chart':p['seconds_per_chart'],
        'gp_sha256':p['gp_sha256'], 'frozen_v3_policy':p,
        'scope':CLAIM,
    }


def validate_roster():
    r = read(D/'roster.json')
    require(r.get('status') == 'READY_V4_EIGHT' and len(r['cases']) == 8,
            'V4 preparation incomplete')
    bindings(ROOT, r['inputs']); bindings(ROOT, r['implementation_sources'])
    _, expected = validate_v3_panel()
    require([x['id'] for x in r['cases']] == [x['id'] for x in expected], 'V4 case order changed')
    for row in r['cases']:
        folder = D/row['id']; p = read(folder/'protocol.json')
        require(sha(folder/'protocol.json') == r['jobs'][row['id']], row['id']+': V4 protocol changed')
        bindings(ROOT, p['inputs']); bindings(ROOT, p['sources']); bindings(ROOT, p['implementation_sources'])
        require(p['bootstrap_fresh_charts'] == FRESH_CHARTS and tuple(p['bootstrap_shells']) == SHELLS,
                row['id']+': V4 bootstrap policy changed')
    return r


def result_summary(results):
    return {
        'completed':len(results),
        'fibres_with_certified_gain':sum(r.get('gain',0)>0 for r in results),
        'maximum_rank_lower_bound':max([17]+[r.get('rank_lower_bound',17) for r in results]),
        'fresh_bootstrap_charts':sum(r.get('bootstrap_charts',0) for r in results),
        'automatic_v3_cascade':False,
        'automatic_new_parameter_expansion':False,
        'claim_boundary':CLAIM,
    }
