import importlib.util
import json
from pathlib import Path
import sys

import pytest

CAS = Path(__file__).resolve().parents[1] / 'cas'
sys.path.insert(0, str(CAS))
spec = importlib.util.spec_from_file_location('seed_panel_parallel', CAS / 'run_curve302_seed_panel_parallel.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_frozen_order_and_occupied_slots():
    order = ['a', 'b', 'c', 'd', 'e']
    assert m.next_cases(order, {'a'}, {'b', 'd'}, 2) == ['c', 'e']
    assert m.next_cases(order, set(), set(), 0) == []


def test_exclusive_lock(tmp_path):
    path = tmp_path / 'lock'
    fd = m.acquire(path)
    try:
        with pytest.raises(RuntimeError, match='already owned'):
            m.acquire(path)
    finally:
        m.os.close(fd)
    m.os.close(m.acquire(path))


def test_reused_pid_is_not_live(monkeypatch):
    monkeypatch.setattr(m, 'identity', lambda pid: {'pid': pid, 'start_token': 'new'})
    assert not m.alive({'pid': 42, 'start_token': 'old'})


def test_receipt_rejects_changed_terminal(tmp_path, monkeypatch):
    monkeypatch.setattr(m.panel, 'D', tmp_path)
    folder = tmp_path / 'seed'
    (folder / 'replay-M17').mkdir(parents=True)
    (folder / 'protocol.json').write_text('{}')
    terminal = folder / 'replay-M17/terminal.json'
    terminal.write_text('{}')
    r = {'status': 'PASS_INDEPENDENT_SEEDED_V3_REPLAY', 'seed_direction': 'seed',
         'initial_rank': 18, 'terminal_sha256': m.sha(terminal),
         'protocol_sha256': m.sha(folder / 'protocol.json')}
    (folder / 'seeded-verified.json').write_text(json.dumps(r))
    assert m.verified('seed') == r
    terminal.write_text('{"changed":true}')
    with pytest.raises(ValueError, match='stale replay'):
        m.verified('seed')


def test_case_uses_original_context_search_and_full_replay(tmp_path, monkeypatch):
    import types
    monkeypatch.setattr(m, 'D', tmp_path)
    monkeypatch.setattr(m, 'MANIFEST', tmp_path / 'manifest.json')
    m.MANIFEST.write_text(json.dumps({'order': ['seed']}))
    monkeypatch.setattr(m, 'check_bindings', lambda _: None)
    completed = set()
    monkeypatch.setattr(m, 'verified', lambda seed: {'rank': 31} if seed in completed else None)
    ctx = types.SimpleNamespace(folder=tmp_path / 'seed')
    calls = []
    def search(c):
        assert c is ctx
        calls.append('search')
    def replay(c):
        assert c is ctx
        calls.append('replay')
        completed.add('seed')
    base = types.SimpleNamespace(context=lambda seed: ctx, verify_case=replay)
    monkeypatch.setattr(m.panel, 'load_base', lambda: base)
    monkeypatch.setitem(sys.modules, 'det1092_v3_worker', types.SimpleNamespace(run_search=search))
    m.run_case('seed')
    assert calls == ['search', 'replay']
    m.run_case('seed')
    assert calls == ['search', 'replay']
    completed.clear()
    calls.clear()
    (ctx.folder / 'replay-M17').mkdir(parents=True)
    (ctx.folder / 'replay-M17/terminal.json').write_text('{}')
    m.run_case('seed')
    assert calls == ['replay']


def test_controller_bounded_concurrency_out_of_order_completion(tmp_path, monkeypatch):
    monkeypatch.setattr(m, 'D', tmp_path)
    monkeypatch.setattr(m, 'MANIFEST', tmp_path / 'manifest.json')
    order = list('abcdef')
    m.MANIFEST.write_text(json.dumps({'order': order, 'workers': 4}))
    monkeypatch.setattr(m, 'check_bindings', lambda _: None)
    monkeypatch.setattr(m.panel, 'sage_launcher', lambda: '/fake/sage')
    monkeypatch.setattr(m.signal, 'signal', lambda *args: None)
    monkeypatch.setattr(m.time, 'sleep', lambda _: None)
    completed = {'a'}
    monkeypatch.setattr(m, 'verified', lambda s: {'rank_lower_bound': 31, 'charts': 5} if s in completed else None)
    launches, active = [], set()
    class Process:
        def __init__(self, argv, **kwargs):
            self.seed = argv[-1]
            self.pid = 100 + len(launches)
            self.remaining = 3 if self.seed == 'b' else 1
            launches.append(self.seed)
            active.add(self.seed)
            assert len(active) <= 4
        def poll(self):
            self.remaining -= 1
            if self.remaining > 0:
                return None
            active.discard(self.seed)
            completed.add(self.seed)
            return 0
    monkeypatch.setattr(m.subprocess, 'Popen', Process)
    published, summaries = [], []
    monkeypatch.setattr(m, 'publish', lambda s: published.append(json.loads(json.dumps(s))))
    monkeypatch.setattr(m, 'summarize', lambda o: summaries.append(list(o)))
    m.controller(m.acquire(tmp_path / 'lock'))
    assert launches == list('bcdef')
    assert summaries == [order]
    assert published[-1]['status'] == 'COMPLETE_SEED_UNIVERSALITY_PANEL'
    assert set(published[-1]['completed']) == set(order)


def test_takeover_refuses_live_search_before_signal(tmp_path, monkeypatch):
    monkeypatch.setattr(m, 'identity', lambda pid: {'pid': pid, 'start_token': '1'})
    with pytest.raises(ValueError, match='sealed search'):
        m.retire_serial({'pid': 123, 'status': 'RUNNING_SEARCH'})
