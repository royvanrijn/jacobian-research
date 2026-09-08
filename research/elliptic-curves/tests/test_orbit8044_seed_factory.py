import json
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

CAS = Path(__file__).resolve().parents[1]/'cas'
sys.path.insert(0, str(CAS))
import det1092_funnel as f
import orbit8044_seed_factory as factory
import run_orbit8044_seed_factory as controller
from v3_warm_support import atomic, read, sha


@pytest.mark.parametrize('height', [1, 2, 3, 7])
def test_complete_height_order(height):
    expected = sorted([(a, b) for b in range(1, height+1) for a in range(-height, height+1)
                       if gcd(a, b) == 1], key=lambda q:(max(abs(q[0]), q[1]), q[1], q[0]))
    actual = list(factory.parameters(height, 100000))
    assert [(r['a'], r['b']) for r in actual] == expected
    assert len({r['u'] for r in actual}) == len(actual)
    assert sum(r['u'] == '0' for r in actual) == 1
    assert list(factory.parameters(height, 2)) == actual[:2]


@pytest.mark.parametrize('args', [(0, 1), (1, 0), (True, 4), (1, 100001)])
def test_enumeration_caps(args):
    with pytest.raises(ValueError):
        list(factory.parameters(*args))


def test_homogeneous_map_and_pole():
    chart = dict(s_numerator=['2', '3', '5'], s_denominator=['1', '0', '-1'])
    assert factory.image(chart, 2, 3) == F(2*9+3*6+5*4, 9-4)
    assert factory.image(chart, 1, 1) is None
    with pytest.raises(ValueError, match='base point'):
        factory.image(dict(s_numerator=['0', '0', '0'], s_denominator=['0', '0', '0']), 1, 1)


def test_retained_chart_polynomial_identity_without_pari():
    chart = read(f.ART/'orbit8044_seed_factory_v1/parameter-chart.json')
    X, Y, Z = [list(map(int, chart[k])) for k in ['s_numerator', 's_denominator', 'z_numerator']]
    q0, q1, q2 = map(int, chart['primitive_conic'])
    def product(a, b):
        return [sum(a[i]*b[k-i] for i in range(3) if 0 <= k-i < 3) for k in range(5)]
    xx, xy, yy, zz = product(X, X), product(X, Y), product(Y, Y), product(Z, Z)
    assert all(z == q0*y+q1*m+q2*x for z, y, m, x in zip(zz, yy, xy, xx))
    assert factory.image(chart, 0, 1) == F(5193, 35630)


@pytest.mark.parametrize('A,B,d', [(-1, 1, 2), (0, 1, 2), (1, 0, 3)])
def test_isomorphism_does_not_merge_twists(A, B, d):
    model = [0, 0, 0, A, B]
    scaled = [0, 0, 0, 16*A, 64*B]
    twist = [0, 0, 0, d*d*A, d**3*B]
    assert factory.rational_isomorphic(model, scaled)
    assert factory.rational_j(model) == factory.rational_j(twist)
    assert not factory.rational_isomorphic(model, twist)


def test_real_dependent_conic_control_and_trace():
    source = f.ART/'det1092_funnel_first_seeds_v1/dependent-conic-witnesses'
    m17 = read(source/'m17.json'); model = tuple(map(F, m17['curve']))
    base = tuple(tuple(map(F, P)) for P in m17['points'])
    points = tuple(tuple(map(F, P)) for P in read(source/'conic-points.json'))
    trace = read(f.COVER)['lift']['trace_word']
    result = factory.obvious_dependencies(model, base, points, trace)
    assert result['trace_identity_checked']
    assert all(r['status'] == 'EXACT_INHERITED_WORD' for r in result['points'])
    for r, point in zip(result['points'], points):
        assert factory.word_point(model, base, r['word']) == point
    with pytest.raises(ValueError, match='trace'):
        factory.obvious_dependencies(model, base, points, [0]*17)


def test_real_independent_conic_is_not_filtered():
    source = f.ART/'det1092_funnel_first_seeds_v1/seeds/funnel-002537010'
    m17 = read(source/'m17.json'); model = tuple(map(F, m17['curve']))
    base = tuple(tuple(map(F, P)) for P in m17['points'])
    points = tuple(tuple(map(F, P)) for P in read(source/'conic-points.json'))
    result = factory.obvious_dependencies(model, base, points, read(f.COVER)['lift']['trace_word'])
    assert [r['status'] for r in result['points']] == ['UNKNOWN', 'UNKNOWN']


def synthetic_verified_seed(tmp_path, monkeypatch, case='seed0', model=None):
    """Synthetic certificate metadata tests queue contracts, not mathematics."""
    monkeypatch.setattr(f, 'ROOT', tmp_path)
    folder = tmp_path/'run'; seed = folder/'seeds'/case
    atomic(folder/'protocol.json', {'test_only':True})
    atomic(seed/'intake.json', {'parameter':'0'})
    packet = dict(status='CERTIFIED_M18', parameter='0', curve=model or [0, 0, 0, -1, 1],
                  intake_sha256=sha(seed/'intake.json'))
    atomic(seed/'m18.json', packet)
    atomic(seed/'standalone-replay.json', dict(status='PASS_STANDALONE_FUNNEL_M18', rank_lower_bound=18,
        matrix_rank=18, inputs={str((seed/'m18.json').relative_to(tmp_path)):sha(seed/'m18.json')}))
    return folder


def test_queue_survives_interrupt_before_admission(tmp_path, monkeypatch):
    folder = synthetic_verified_seed(tmp_path, monkeypatch)
    normal = factory.atomic
    def interrupted(path, value, **kwargs):
        if path.parent.name == 'admissions':
            raise RuntimeError('interruption after queue publication')
        return normal(path, value, **kwargs)
    monkeypatch.setattr(factory, 'atomic', interrupted)
    with pytest.raises(RuntimeError):
        factory.queue_admission(folder, 'seed0', [])
    assert (folder/'v3-queue/seed0.json').exists()
    monkeypatch.setattr(factory, 'atomic', normal)
    result = factory.queue_admission(folder, 'seed0', [])
    assert result == factory.queue_admission(folder, 'seed0', [])
    assert len(list((folder/'v3-queue').glob('*.json'))) == 1


def test_queue_duplicate_and_twist(tmp_path, monkeypatch):
    folder = synthetic_verified_seed(tmp_path, monkeypatch)
    model = [0, 0, 0, -16, 64]
    registry = [dict(case='prior', curve=model, j=str(factory.rational_j(model)))]
    result = factory.queue_admission(folder, 'seed0', registry)
    assert result['status'] == 'DUPLICATE_RATIONAL_ISOMORPHISM'
    assert not (folder/'v3-queue/seed0.json').exists()
    synthetic_verified_seed(tmp_path, monkeypatch, 'twist', [0, 0, 0, -4, 8])
    assert factory.queue_admission(folder, 'twist', registry)['status'] == 'QUEUED_M18'


def test_queue_rejects_changed_seed(tmp_path, monkeypatch):
    folder = synthetic_verified_seed(tmp_path, monkeypatch)
    path = folder/'seeds/seed0/m18.json'; packet = read(path)
    packet['curve'] = [0, 0, 0, -4, 8]; atomic(path, packet)
    with pytest.raises(ValueError, match='changed after proof'):
        factory.queue_admission(folder, 'seed0', [])


def test_receipt_recovery_does_not_rerun_or_reset_budget(tmp_path, monkeypatch):
    terminal = tmp_path/'result.json'; count = []
    def interrupted_run(command, *, limits, checkpoint_path, **kwargs):
        count.append(command)
        atomic(terminal, {'status':'done'})
        atomic(checkpoint_path, dict(command=command, outcome='completed', returncode=0,
            wall_seconds=5, limits={'wall_seconds':limits.wall_seconds}))
        raise RuntimeError('controller dies after supervised completion')
    monkeypatch.setattr(controller, 'run', interrupted_run)
    jobs = controller.Stages(tmp_path, 'factory', 5, 1024)
    with pytest.raises(RuntimeError):
        jobs.stage('first', ['test'], 10, terminal)
    resumed = controller.Stages(tmp_path, 'factory', 5, 1024)
    assert resumed.stage('first', ['test'], 10, terminal)
    assert len(count) == 1
    assert not resumed.stage('second', ['test2'], 10, terminal)
    assert resumed.exhausted
    atomic(terminal, {'status':'changed'})
    with pytest.raises(ValueError, match='changed'):
        resumed.stage('first', ['test'], 10, terminal)


def test_unreceipted_interruption_consumes_reservation(tmp_path):
    atomic(tmp_path/'factory-controller.json', dict(stages={'first':dict(status='RUNNING', reserved_seconds=10)}, cursor=0))
    jobs = controller.Stages(tmp_path, 'factory', 10, 1024)
    assert not jobs.stage('first', ['test'], 10, tmp_path/'missing.json')
    assert jobs.data['stages']['first']['status'] == 'CENSORED_INTERRUPTION'
    assert not jobs.stage('next', ['test2'], 10, tmp_path/'missing.json')
    assert jobs.exhausted


def test_direct_quadratic_branches_match_independent_retained_controls():
    sage = shutil.which('sage')
    if not sage:
        pytest.skip('Sage integration requires installed runtime')
    code = '''
import sys
sys.path.insert(0,sys.argv[1])
import det1092_funnel as f
from v3_warm_support import read,point_tuple
from orbit8044_seed_factory_worker import both_points
paths=[f.ART/'det1092_funnel_first_seeds_v1/seeds'/c for c in ['funnel-019791672','funnel-002537010','funnel-020887854']]
paths += [f.ART/'det1092_funnel_first_seeds_v1/dependent-conic-witnesses',f.ART/'det1092_small_conic_seed_v1/seeds/conic-small-01']
for path in paths:
    actual,error=both_points(read(path/'intake.json'))
    if path.name=='funnel-019791672':
        assert error=='UNRESOLVED_TWO_DISTINCT_RATIONAL_BRANCHES'
    else:
        assert error is None and tuple(actual)==point_tuple(read(path/'conic-points.json'))
print('PASS_DIRECT_RESIDUAL_BRANCH_CONTROLS')
'''
    result = subprocess.run([sage, '-python', '-c', code, str(CAS)], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout+result.stderr
