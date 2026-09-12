"""Validate the original arithmetic evidence despite later live inventory edits.

Only the explicitly pinned inventory is allowed to resolve from Git, and its
bytes must match the original plan SHA256. Never repin old evidence or replace
the user's live inventory. All other frozen sources must still match on disk.
"""
import json
from pathlib import Path
import subprocess

import historical_external_arithmetic as old
import wide_arithmetic_profile_core as core

INVENTORY_SOURCE = 'elliptic-curves/data/research_curves/database.json'
INVENTORY_REVISION = '1959f550ca43e5ff4492e8b5eb383f3104562a1a'


def verify():
    path = old.DEFAULT_OUTPUT
    plan = old.read(path/'plan.json')
    old.require(old.tree_hashes(Path(plan['census'])) == plan['prospective_files'], 'frozen census bytes changed')
    resolution = {}
    inventory = None
    for source, digest in plan['sources'].items():
        p = old.ROOT/source
        if old.sha(p) == digest:
            resolution[source] = {'mode':'live_matching_bytes','sha256':digest}
            if source == INVENTORY_SOURCE:
                inventory = old.read(p)
        else:
            old.require(source == INVENTORY_SOURCE, f'frozen non-inventory source changed: {source}')
            git_path = 'research/'+source
            raw = subprocess.check_output(['git','show',INVENTORY_REVISION+':'+git_path], cwd=old.ROOT.parent)
            old.require(core.sha256_bytes(raw) == digest, 'archived inventory does not match original frozen hash')
            inventory = json.loads(raw)
            resolution[source] = {'mode':'exact_archived_git_bytes','revision':INVENTORY_REVISION,
                                  'sha256':digest,'live_sha256_at_check':old.sha(p)}
    for name, digest in plan['evidence_hashes'].items():
        old.require(old.sha(path/name)==digest, 'historical evidence changed: '+name)
    roster = [{k:r[k] for k in ('id','family','parameter','rank_lower_bound','rank_provenance')} for r in inventory['curves']]
    old.require(roster == old.read(path/'inventory_roster.json')['rows'], 'original inventory roster mismatch')
    old.require({r['id'] for r in old.select(roster)}=={r['id'] for r in plan['rows']}, 'historical selection mismatch')
    for row in plan['rows']:
        inp = path/'inputs'/f"{row['curve_key']}.json"
        old.require(inp.read_text()==core.stable_json(row), 'historical input mismatch')
        for mode in ('base','local'):
            result = old.read(path/mode/inp.name)
            old.require(result['curve_key']==row['curve_key'] and result['input_sha256']==old.sha(inp), 'historical checkpoint binding mismatch')
            for k in ('timeout_seconds','memory_gb'):
                old.require(result['runtime'][k]==plan['policies'][mode][k], 'historical policy mismatch')
    replay = old.read(path/'ARITHMETIC_REPLAY.json')
    for name, digest in replay['checkpoint_sha256'].items():
        old.require(old.sha(path/name)==digest, 'arithmetic replay changed')
    report = old.read(path/'REPORT.json')
    old.require(report['plan_sha256']==old.sha(path/'plan.json'), 'historical plan changed')
    old.require(report['arithmetic_replay_sha256']==old.sha(path/'ARITHMETIC_REPLAY.json'), 'replay receipt changed')
    for name, digest in report['output_sha256'].items():
        old.require(old.sha(path/name)==digest, 'historical analysis changed: '+name)
    return plan, resolution
