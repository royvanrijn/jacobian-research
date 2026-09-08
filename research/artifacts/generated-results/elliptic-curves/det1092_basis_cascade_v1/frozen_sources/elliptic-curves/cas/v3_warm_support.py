"""Sage-free checkpoint contracts for the maintained warm-transfer runner.

File order is numeric, never lexical. These checks authorize continuation or
independent replay; they do not themselves establish a rank lower bound.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import tempfile
from fractions import Fraction
from pathlib import Path

WARM = ('warm-11952-41', 'warm-11952-72', 'warm-11952-186')
TERMINAL_REASONS = {
    'COMPLETE_FINITE_NO_GAIN', 'TARGET_LOWER_BOUND_REACHED',
    'CHART_BUDGET_EXHAUSTED', 'CENSORED_SEARCH', 'INCOMPLETE_EPOCH',
    'EPOCH_BUDGET_EXHAUSTED', 'NO_NEW_FINALISTS',
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    with Path(path).open() as stream:
        return json.load(stream)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def atomic(path, value, *, immutable=False):
    """Publish atomically. Immutable endpoints may only be reused identically."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if immutable and path.exists():
        require(read(path) == value, f'immutable endpoint differs: {path}')
        return
    fd, temp = tempfile.mkstemp(prefix='.' + path.name + '.', dir=path.parent)
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(value, stream, sort_keys=True, indent=2)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
        if immutable:
            try:
                os.link(temp, path)
            except FileExistsError:
                require(read(path) == value, f'concurrent immutable endpoint differs: {path}')
        else:
            os.replace(temp, path)
        directory_fd = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        if os.path.exists(temp):
            os.unlink(temp)


def within(root, name):
    root = Path(root).resolve()
    path = (root / name).resolve()
    require(path.is_relative_to(root), f'path escapes research directory: {name}')
    return path


def bindings(root, values):
    require(isinstance(values, dict) and values, 'empty or invalid input binding set')
    for name, digest in values.items():
        require(sha(within(root, name)) == digest, f'changed bound file: {name}')


def point_tuple(points):
    """Normalize exact VALUES, preserving order and sign; never use floats."""
    result = tuple(tuple(Fraction(x) for x in point) for point in points)
    require(all(len(p) == 2 for p in result), 'expected affine coordinate pairs')
    return result


def curve_tuple(curve):
    result = tuple(Fraction(x) for x in curve)
    if len(result) == 2:
        result = (Fraction(0),) * 3 + result
    require(len(result) == 5 and result[:3] == (0, 0, 0), 'short Weierstrass model required')
    return result


def assert_basis(state, model, points, expected_rank):
    actual = point_tuple(state.basis)
    expected = point_tuple(points)
    require(curve_tuple(state.model.coefficients) == curve_tuple(model), 'seed model differs')
    require(state.rank == expected_rank, f'seed rank {state.rank} != {expected_rank}')
    if actual != expected:
        diff = next((i for i, pair in enumerate(zip(actual, expected)) if pair[0] != pair[1]),
                    min(len(actual), len(expected)))
        raise ValueError(f'ordered seed coordinates differ at index {diff}; lengths {len(actual)}/{len(expected)}')


def indexed_paths(folder, prefix='chart', *, expected=None):
    """Return a contiguous 0-based sequence; reject holes and aliased filenames.

    :03d is a MINIMUM width: chart-1000 sorts before chart-101 as a string.
    Do not rename historical evidence or drop the centre/index assertions.
    """
    pattern = re.compile(re.escape(prefix) + r'-(\d+)\.json\Z')
    found = {}
    for path in Path(folder).glob(prefix + '-*.json'):
        match = pattern.fullmatch(path.name)
        require(match is not None, f'malformed indexed checkpoint: {path}')
        index = int(match.group(1))
        require(index not in found, f'duplicate {prefix} index {index}')
        found[index] = path
    indices = sorted(found)
    require(indices == list(range(len(indices))), f'noncontiguous {prefix} checkpoints in {folder}')
    if expected is not None:
        require(len(indices) == expected, f'{prefix} count {len(indices)} != {expected}')
    return [found[i] for i in indices]


def check_chart(chart, selection, index):
    require(type(chart.get('index')) is int and chart['index'] == index,
            f'chart filename/payload index mismatch at {index}')
    centres = selection['centres']
    require(index < len(centres), f'chart {index} lies outside frozen schedule')
    require(chart['centre'] == centres[index], f'chart centre differs from frozen schedule at index {index}')


def epoch_paths(folder):
    entries = []
    for path in Path(folder).glob('epoch-*'):
        match = re.fullmatch(r'epoch-(\d+)', path.name)
        require(match is not None and path.is_dir(), f'malformed epoch: {path}')
        entries.append((int(match.group(1)), path))
    entries.sort()
    require([i for i, _ in entries] == list(range(len(entries))), 'epoch gap/alias')
    return [p for _, p in entries]


def terminal_structure(folder):
    """Validate a sealed terminal and schedule, WITHOUT promoting its rank."""
    folder = Path(folder)
    endpoint = folder / 'replay-M17/terminal.json'
    if not endpoint.exists():
        return None
    terminal = read(endpoint)
    policy = read(folder / 'protocol.json')
    require(terminal.get('status') == 'TERMINAL_BOUNDED_TRANSFER', 'unexpected terminal schema')
    require(terminal['protocol_sha256'] == sha(folder / 'protocol.json'), 'terminal protocol mismatch')
    require(terminal['initial_rank'] == policy['initial_rank'], 'terminal seed rank mismatch')
    stages = terminal['stages']
    require(bool(stages), 'terminal has no completed stage')
    rank, total = policy['initial_rank'], 0
    paths = epoch_paths(folder / 'replay-M17')
    require(len(paths) == len(stages), 'terminal/epoch count mismatch')
    for i, (stage, wd) in enumerate(zip(stages, paths)):
        require(stage['epoch'] == i and stage['before'] == rank, 'broken terminal rank chain')
        require(read(wd / 'stage.json') == stage, 'terminal differs from sealed stage')
        selection = read(wd / 'selection.json')
        require(selection['rank'] == rank and len(selection['basis']) == rank, 'selection rank mismatch')
        charts = indexed_paths(wd, expected=stage['charts'])
        censored = 0
        for j, path in enumerate(charts):
            chart = read(path)
            check_chart(chart, selection, j)
            censored += chart['search']['status'] != 'bounded_search_complete'
        require(censored == stage.get('censored_charts', censored), 'censored count mismatch')
        require(stage['after'] >= rank, 'certified subgroup cannot shrink')
        require(stage['full_cosets_scored'] == selection['full_cosets_scored'], 'coset total mismatch')
        audit = within(wd, stage['audit'])
        require(sha(audit) == stage['audit_sha256'], 'sealed audit changed')
        if charts:
            require(audit.name == f'mod2-{len(charts)-1:03d}.json', 'stage does not bind last chart audit')
        cancelled = len(selection['centres']) - len(charts) if stage['after'] > rank else 0
        require(stage['stale_charts_cancelled'] == cancelled, 'cancelled schedule mismatch')
        if stage['stop_reason'] == 'COMPLETE_FINITE_NO_GAIN':
            require(stage['after'] == rank and not censored and len(charts) == len(selection['centres']),
                    'no-gain terminal does not cover its complete uncensored schedule')
        if i < len(stages)-1:
            require(stage['stop_reason'] == 'REBUILD_AFTER_CERTIFIED_GAIN' and stage['after'] > rank,
                    'continued after terminal/no-gain epoch')
        rank, total = stage['after'], total + len(charts)
    require(total == terminal['charts'] and rank == terminal['final_rank_lower_bound'], 'terminal totals mismatch')
    require(terminal['stop_reason'] in TERMINAL_REASONS, 'unknown terminal stop reason')
    require(terminal['stop_reason'] == stages[-1]['stop_reason'] or
            terminal['stop_reason'] == 'EPOCH_BUDGET_EXHAUSTED', 'stage/terminal stop mismatch')
    return terminal


def process_info(pid, proc_root=Path('/proc')):
    try:
        pid = int(pid)
        if pid <= 0:
            return None
        text = (Path(proc_root) / str(pid) / 'stat').read_text()
        fields = text.rsplit(')', 1)[1].split()
        return {'pid': pid, 'state': fields[0], 'start_token': fields[19], 'ppid': int(fields[1])}
    except (OSError, ValueError, IndexError, TypeError):
        return None


def same_process(pid, token):
    row = process_info(pid)
    return bool(row and token and row['start_token'] == str(token) and row['state'] not in ('Z', 'X'))


def tail(path, byte_limit=65536):
    path = Path(path)
    if not path.exists():
        return None
    with path.open('rb') as stream:
        stream.seek(0, os.SEEK_END)
        stream.seek(max(0, stream.tell() - byte_limit))
        return stream.read(byte_limit).decode('utf-8', 'replace')
