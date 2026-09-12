#!/usr/bin/env python3
"""Audit the retained rescue definition and source amendment without a search.

This checks protocol bytes, deterministic assignments, historical chunk hashes
and reported counts. It does not replay point arithmetic or certify inputs for
the amended runtime. Only the standard library and read-only git show are used.
"""

import ast
import hashlib
import json
from pathlib import Path
import runpy
import subprocess


ROOT = Path(__file__).resolve().parents[2]
BUILDER = ROOT / 'elkies-k3/scripts/build_mw17_jump_v2_zero_gain_rescue.py'
MIGRATION = ROOT / 'elliptic-curves/cas/pointed_quartic_migration.py'
GENERATED = ROOT / 'artifacts/generated-results'
PROTOCOL = GENERATED / 'elkies-k3-mw17-jump-v2-zero-gain-rescue-arm-v1.json'
LEDGER = GENERATED / 'elkies-k3-mw17-jump-v2-ledger-v1.json'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def canonical_hash(value):
    return digest(json.dumps(value, sort_keys=True, separators=(',', ':')).encode())


def compare_definition(retained, current, migrated, historical_source):
    """Allow only the already declared implementation amendment, not new data."""
    for doc in (retained, current):
        definition = {k: v for k, v in doc.items() if k != 'protocol_definition_sha256'}
        if canonical_hash(definition) != doc['protocol_definition_sha256']:
            raise ValueError('protocol definition digest differs')
    old = retained['implementation_hashes']
    new = current['implementation_hashes']
    if old.keys() != new.keys():
        raise ValueError('implementation roster differs')
    amended = []
    for name, expected in old.items():
        if new[name] == expected:
            continue
        if name not in migrated or digest(historical_source(name)) != expected:
            raise ValueError('source change outside the pinned amendment: ' + name)
        amended.append(name)
    ignored = {'implementation_hashes', 'protocol_definition_sha256'}
    if {k: v for k, v in retained.items() if k not in ignored} != {
        k: v for k, v in current.items() if k not in ignored
    }:
        raise ValueError('frozen population, assignment, context or budget differs')
    return amended


def audit():
    # Read literal amendment declarations without importing its search backend.
    constants = {}
    for node in ast.parse(MIGRATION.read_text()).body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id in {'MIGRATED', 'REGRESSION_REVISION'}:
                    constants[target.id] = ast.literal_eval(node.value)
    retained = json.loads(PROTOCOL.read_text())
    # The inspected builder uses only deterministic metadata and file hashes.
    current = runpy.run_path(str(BUILDER))['build']()

    def historical_source(name):
        return subprocess.check_output(
            ['git', 'show', constants['REGRESSION_REVISION'] + ':' + name], cwd=ROOT,
        )

    amended = compare_definition(retained, current, constants['MIGRATED'], historical_source)
    raw = LEDGER.read_bytes()
    ledger = json.loads(raw)
    context = retained['known_context_at_freeze']
    if digest(raw) != context['ledger_sha256_at_protocol_freeze']:
        raise ValueError('freeze-time base ledger changed')
    records = []
    inputs = {str(PROTOCOL.relative_to(ROOT)): digest(PROTOCOL.read_bytes()),
              str(LEDGER.relative_to(ROOT)): digest(raw)}
    for chunk in ledger['chunk_provenance']:
        p = ROOT / chunk['path']
        raw = p.read_bytes()
        if digest(raw) != chunk['sha256']:
            raise ValueError('historical chunk changed: ' + chunk['path'])
        rows = json.loads(raw)['records']
        if len(rows) != chunk['record_count']:
            raise ValueError('historical chunk count differs')
        records.extend(rows)
        inputs[chunk['path']] = digest(raw)
    if len({row['campaign_index'] for row in records}) != len(records):
        raise ValueError('duplicate historical campaign index')
    if len(records) != ledger['completed_worker_count']:
        raise ValueError('base ledger count differs')
    assigned = {row['sample_id'] for row in retained['assignments'] if row['assigned_to_rescue_arm']}
    if len(assigned) != retained['assignment']['assigned_candidate_count']:
        raise ValueError('assigned count differs')
    clean = []
    for row in records:
        covers = row['initial']['cover_records']
        if not (row['status'] == 'PASS_EXACT_CERTIFIED_QUOTIENT_GAIN'
                and row['actual_certified_quotient_rank_gain'] == 0
                and row['attempted_chart_count'] == len(covers) == 43
                and row['bounded_cover_timeout_count'] == row['cover_backend_failure_count'] == 0
                and all(c['search']['status'] == 'bounded_search_complete' for c in covers)):
            raise ValueError('historical clean-zero scope differs')
        clean.append(row)
    return dict(
        status='PASS_RETAINED_PROTOCOL_AND_HISTORICAL_CHECKPOINTS',
        candidates=len(retained['assignments']), assigned=len(assigned),
        historical_clean_zeros=len(clean),
        assigned_historical_clean_zeros=sum(row['sample_id'] in assigned for row in clean),
        historical_completed_charts=sum(row['attempted_chart_count'] for row in clean),
        amended_implementation_paths=amended,
        historical_source_revision=constants['REGRESSION_REVISION'], inputs=inputs,
        current_runtime_base_chunks_present=len(list((ROOT / 'artifacts/local/elkies-k3/mw17-jump-v2-pointed-v1').glob('chunk-*-of-*.json'))),
        execution_readiness='NOT_CERTIFIED',
        boundary='Metadata and retained-byte audit only. Old raw-backend results cannot be resumed under changed coordinate boxes. No point, CVP, parameter or descent calculation ran.',
    )


if __name__ == '__main__':
    print(json.dumps(audit(), indent=2, sort_keys=True))
