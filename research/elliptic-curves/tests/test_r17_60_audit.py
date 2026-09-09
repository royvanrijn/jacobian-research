import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'cas'))
import audit_r17_60_panel as audit


@pytest.fixture
def bounded_miss(tmp_path, monkeypatch):
    monkeypatch.setattr(audit.panel, 'D', tmp_path)
    monkeypatch.setattr(audit.panel, 'ROOT', tmp_path)
    def write(path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value))
    row = {'id': 'case', 'family': 'f', 'parameter': '1', 'stratum': 'low'}
    folder = tmp_path / 'cases/case'
    seed = folder / 'seed-search'
    write(tmp_path / 'protocol.json', {})
    write(seed / 'protocol.json', {'max_point_invocations': 86})
    write(seed / 'terminal.json', {'charts': 1, 'rank_lower_bound': 17,
                                  'protocol_sha256': audit.sha(seed / 'protocol.json')})
    write(seed / 'verified.json', {'status': 'PASS_INDEPENDENT_FIRST_SEED_REPLAY',
                                  'terminal_sha256': audit.sha(seed / 'terminal.json')})
    write(seed / 'chart-0000.json', {})
    packet = {'bindings': {}, 'rank_lower_bound': 17, 'point_searches': 0, 'reconciliation_gains': []}
    write(folder / 'seed-reconciled.json', packet)
    result = dict(row, status='PASS_INDEPENDENT_R17_60_CASE', bindings={},
                  protocol_sha256=audit.sha(tmp_path / 'protocol.json'), first_M18_found=False,
                  packet=packet, gain_timeline=[], rank_lower_bound=17, seed_calls=1,
                  complement_calls=0, total_calls=1, added_directions=0, calls_per_added_direction=None,
                  gaining_complement_calls=0, last_complement_gain_call=None, complement_calls_since_last_gain=0)
    write(folder / 'result.json', result)
    return row, folder, result, write


def test_valid_bounded_miss(bounded_miss):
    row, _, _, _ = bounded_miss
    assert audit.case_audit(row)['rank_lower_bound'] == 17


def test_rejects_fabricated_gain_timing(bounded_miss):
    row, folder, result, write = bounded_miss
    result['last_complement_gain_call'] = 99
    write(folder / 'result.json', result)
    with pytest.raises(Exception, match='reported gain timing differs'):
        audit.case_audit(row)


def test_rejects_extra_unaccounted_chart(bounded_miss):
    row, folder, _, write = bounded_miss
    write(folder / 'seed-search/chart-0001.json', {})
    with pytest.raises(Exception, match='seed chart count differs'):
        audit.case_audit(row)


def test_rejects_budget_expansion(bounded_miss):
    row, folder, _, write = bounded_miss
    write(folder / 'complement-00/protocol.json', {'max_charts': 101})
    with pytest.raises(Exception, match='remaining budget differs'):
        audit.case_audit(row)


def test_rejects_stale_replay_receipt(bounded_miss):
    row, folder, _, write = bounded_miss
    write(folder / 'seed-search/terminal.json', {'changed': True})
    with pytest.raises(Exception, match='seed replay differs'):
        audit.case_audit(row)
