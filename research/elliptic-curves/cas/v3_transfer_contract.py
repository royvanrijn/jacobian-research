"""Fail-closed input/replay contracts for the first out-of-302 V3 transfer.

No Sage imports, network access, point search, or mathematical-status writes.
All rank assertions are inherited from separately replayed certificates.
"""
from __future__ import annotations

import ast
import gzip
import hashlib
import json
import os
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
from typing import Any

BASE_COMMIT = '305d4e788b370bd9fdce9321dde7441fc68ea7c3'
CASE_SPECS = (
    {'id': 'control-native11952', 'kind': 'positive-control', 'initial_rank': 17},
    {'id': 'warm-11952-41', 'kind': 'warm-inventory', 'initial_rank': 27,
     'parameter': '-2448/11', 'inventory_id': 'new-20260906-41',
     'export': 'new_paired_rank27_curve_11952.sage'},
    {'id': 'warm-11952-72', 'kind': 'warm-inventory', 'initial_rank': 27,
     'parameter': '2012/211', 'inventory_id': 'new-20260906-72',
     'export': 'new_retention_rank27_curve_11952.sage'},
    {'id': 'warm-11952-186', 'kind': 'warm-inventory', 'initial_rank': 27,
     'parameter': '4286/1881', 'inventory_id': 'new-20260906-186',
     'export': 'new_full11952_high_rank_curves.sage'},
)
# Transfer budget only. Mathematical selector and per-chart bounds come from V3.
LIMITS = {'max_charts': 4096, 'max_epochs': 16, 'target_rank': 32,
          'orbit_node_limit': 20000000, 'prepare_wall_seconds': 3600,
          'case_wall_seconds': 7200, 'replay_wall_seconds': 7200,
          'rss_bytes': 3221225472, 'maximum_workers': 1}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def resolve_json(path: Path) -> Path:
    path = Path(path)
    if path.is_file():
        return path
    zipped = Path(str(path) + '.gz')
    if zipped.is_file():
        return zipped
    raise FileNotFoundError(f'Required retained input missing: {path} (or .gz). Restore it; do not regenerate a different experiment.')


def read_json(path: Path) -> Any:
    path = resolve_json(path)
    raw = gzip.decompress(path.read_bytes()) if path.suffix == '.gz' else path.read_bytes()
    return json.loads(raw)


def write_new(path: Path, value: Any) -> None:
    """Publish a complete immutable file; never overwrite an existing endpoint."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = (json.dumps(value, sort_keys=True, indent=2) + '\n').encode()
    temp = path.with_name(path.name + f'.partial-{os.getpid()}')
    try:
        with temp.open('xb') as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temp, path)  # atomic no-clobber publication, including concurrent calls
    finally:
        if temp.exists():
            temp.unlink()


def within(root: Path, name: str) -> Path:
    root = root.resolve()
    path = (root / name).resolve()
    require(path.is_relative_to(root), f'Path escapes research root: {name}')
    return path


def check_bindings(root: Path, bindings: dict[str, str]) -> None:
    require(bool(bindings), 'Empty binding set')
    for name, digest in bindings.items():
        require(sha(within(root, name)) == digest, f'Changed bound input: {name}')


def on_curve(model, point) -> bool:
    a1, a2, a3, a4, a6 = map(F, model)
    x, y = map(F, point)
    return y*y + a1*x*y + a3*y == x*x*x + a2*x*x + a4*x + a6


def extract_sage_literals(text: str):
    """Read the existing single-curve Sage exports WITHOUT executing their code."""
    tree = ast.parse(text)
    models, point_lists = [], []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == 'EllipticCurve':
            require(len(node.args) == 2, 'Unexpected curve export call')
            models.append(ast.literal_eval(node.args[1]))
        if isinstance(node, ast.Assign):
            names = {t.id for t in node.targets if isinstance(t, ast.Name)}
            if 'point_coordinates' in names:
                point_lists.append(ast.literal_eval(node.value))
            elif 'points' in names and isinstance(node.value, ast.ListComp):
                require(len(node.value.generators) == 1, 'Unexpected point comprehension')
                point_lists.append(ast.literal_eval(node.value.generators[0].iter))
    require(len(models) == len(point_lists) == 1, 'Require a single literal curve and point list')
    model = tuple(F(str(x)) for x in models[0])
    points = tuple(tuple(F(str(x)) for x in p) for p in point_lists[0])
    require(len(model) == 5 and all(len(p) == 2 for p in points), 'Malformed model/points')
    require(all(on_curve(model, p) for p in points), 'Export contains an off-curve point')
    return model, points


def sqrt_rational(q: F) -> F:
    q = F(q)
    require(q >= 0, 'Nonrational Weierstrass scale')
    a, b = isqrt(q.numerator), isqrt(q.denominator)
    require(a*a == q.numerator and b*b == q.denominator, 'Nonrational Weierstrass scale')
    return F(a, b)


def transport_to_short(source, target, points):
    """Exact Q-isomorphism to a specified short model; no j-only matching."""
    a1, a2, a3, a4, a6 = map(F, source)
    target = tuple(map(F, target))
    require(len(target) == 5 and target[:3] == (0, 0, 0), 'Expected short target model')
    b2, b4, b6 = a1*a1+4*a2, 2*a4+a1*a3, a3*a3+4*a6
    c4, c6 = b2*b2-24*b4, -b2**3+36*b2*b4-216*b6
    A, B = -27*c4, -54*c6
    require(A*B*target[3]*target[4] != 0, 'Exceptional j needs a separate adapter')
    u = sqrt_rational(B*target[3]/(A*target[4]))
    require(u != 0 and A == u**4*target[3] and B == u**6*target[4], 'Equations are not Q-isomorphic')
    result = tuple(((36*F(x)+3*b2)/u**2, 108*(2*F(y)+a1*F(x)+a3)/u**3) for x, y in points)
    require(all(on_curve(target, p) for p in result), 'Point transport failed')
    return result, {'scale': str(u), 'b2': str(b2), 'source_ainvs': list(map(str, source))}


def check_orbit_rows(gram, rows):
    """Witness validation, not a proof of completeness or coset minimality."""
    n = len(gram)
    require(n == 17 and all(len(r) == n for r in gram), 'This first adapter is rootless MW17 only')
    masks = set()
    for row in rows:
        w = row['word']; norm = row['norm']; mask = row['mask']
        require(len(w) == n and norm in (8, 10), 'Unsupported shell witness')
        require(mask == sum((int(x) % 2) << i for i, x in enumerate(w)), 'Wrong basis/parity transport')
        require(mask != 0 and mask not in masks, 'Duplicate or zero orbit')
        require(sum(w[i]*gram[i][j]*w[j] for i in range(n) for j in range(n)) == norm, 'Wrong generic shell norm')
        masks.add(mask)
    require(bool(rows), 'Empty parent bank')


def validate_gate(root: Path, gate_dir: Path, replay_path: Path, metric_path: Path) -> dict:
    """Bind the COMPLETED M17 replay, not a live progress file or final point list.

    This validates replay provenance/completeness. It does not rerun Sage here.
    Numerical/finite-group calculations are performed by the pinned checkers.
    """
    root, gate_dir = root.resolve(), gate_dir.resolve()
    p = read_json(gate_dir/'protocol.json')
    t = read_json(gate_dir/'replay-M17/terminal.json')
    r, m = read_json(replay_path), read_json(metric_path)
    ph = sha(gate_dir/'protocol.json')
    require(t.get('status') == 'COMPLETE_BOUNDED_CALIBRATION', 'V3 has not terminated normally')
    require(t.get('final_rank_lower_bound', 0) >= 31, 'V3 rank-31 execution is not closed')
    require(r.get('schema') == 'visibility-cascade-v3-replay', 'Need final full replay, not a progress checkpoint')
    require(r.get('rank_lower_bound') == t['final_rank_lower_bound'], 'Final replay is incomplete')
    require(r.get('guard_rejected_unredacted_parent') is True, 'Missing artifact-boundary replay')
    require(all(x.get('protocol_sha256') == ph for x in (r, m, t)), 'Replay protocol mismatch')
    check_bindings(root, p['sources']); check_bindings(root, p['inputs'])
    checker = root/'elliptic-curves/cas/check_visibility_cascade_v3.sage'
    metric_checker = root/'elliptic-curves/cas/check_visibility_metric_v3.sage'
    require(r.get('checker_sha256') == sha(checker), 'V3 checker changed')
    require(m.get('checker_sha256') == sha(metric_checker), 'V3 metric checker changed')
    require(len(r['stages']) == len(t['stages']) == len(m['stages']) > 0, 'Incomplete replay/metric coverage')
    rank, charts = 17, 0
    for i, (s, q, h) in enumerate(zip(t['stages'], r['stages'], m['stages'])):
        require(s['epoch'] == q['epoch'] == i, 'Epoch gap/reordering')
        require(s['before'] == q['before'] == h['rank'] == rank, 'Replay chain is not rooted at M17')
        require(s['after'] == q['after'] and s['after'] > rank, 'Unexpected terminal/no-gain stage')
        require(s['charts'] == q['charts_replayed'], 'Incomplete chart replay')
        require(s['full_cosets_scored'] == q['full_cosets_replayed'], 'Incomplete landscape replay')
        require(all(q['independent_modl_ranks'].get(str(ell)) == s['after'] for ell in (3, 5)), 'Missing odd-prime rank certificate')
        wd = gate_dir/f'replay-M17/epoch-{i:02d}'
        require(sha(wd/'selection.json') == h['sha256'], 'Metric selection mismatch')
        require(within(root, h['selection']) == wd/'selection.json', 'Wrong metric epoch')
        require(sha(wd/s['audit']) == s['audit_sha256'], 'Final epoch cloud changed')
        check_bindings(root, q['checkpoint_hashes'])
        required = {str(x.relative_to(root)) for x in wd.iterdir() if x.suffix in ('.json', '.npz')}
        require(required == set(q['checkpoint_hashes']), 'Replay did not bind the full epoch')
        rank, charts = s['after'], charts+s['charts']
    require(rank == r['rank_lower_bound'] and charts == t['charts'], 'Terminal totals disagree')
    paths = (gate_dir/'protocol.json', gate_dir/'replay-M17/terminal.json', replay_path, metric_path)
    return {'status': 'PASS_COMPLETED_V3_GATE', 'rank_lower_bound': rank, 'charts': charts,
            'bindings': {str(x.resolve().relative_to(root)): sha(x) for x in paths},
            'scope': 'Completed independently replayed 302 calibration. Not a new curve, general sensitivity theorem, or basis-invariance proof.'}


def classify_stop(*, rank: int, target: int, charts: int, max_charts: int,
                  executed: int, scheduled: int, censored: int, gain: bool) -> str:
    if rank >= target:
        return 'TARGET_LOWER_BOUND_REACHED'
    if charts >= max_charts:
        return 'CHART_BUDGET_EXHAUSTED'
    if censored:
        return 'CENSORED_SEARCH'
    if not gain:
        return 'COMPLETE_FINITE_NO_GAIN' if executed == scheduled else 'INCOMPLETE_EPOCH'
    return 'REBUILD_AFTER_CERTIFIED_GAIN'
