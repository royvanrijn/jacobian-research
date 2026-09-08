import importlib.util
import json
from pathlib import Path
import sqlite3
import sys

import pytest

CAS = Path(__file__).resolve().parents[1]/'cas'
sys.path.insert(0, str(CAS))
import det1092_funnel as f


def policy():
    return dict(population=97, maximum_draws=2000, block_size=13,
                shell_bits=[4, 6, 8], domain='test-funnel-resume', primes=[5,7,13],
                height_bands=list(range(32)), retain_strong_per_band=5,
                retain_control_per_band=3, seed_strong_per_band=2,
                seed_control_per_band=1, arithmetic_cap=20, seed_cap=7,
                resources={'storage_bytes':10**8})


@pytest.mark.parametrize('when', ['before_commit', 'after_commit'])
def test_crash_resume_matches_uninterrupted(tmp_path, when):
    p = policy(); arithmetic = f.Arithmetic(p)
    uninterrupted = f.scan(tmp_path/'normal', p, arithmetic)
    def crash(place):
        if place == when:
            raise RuntimeError('injected process interruption')
    with pytest.raises(RuntimeError, match='injected'):
        f.scan(tmp_path/'interrupted', p, arithmetic, failpoint=crash)
    resumed = f.scan(tmp_path/'interrupted', p, arithmetic)
    assert resumed == uninterrupted
    assert f.selection(resumed,p,arithmetic) == f.selection(uninterrupted,p,arithmetic)
    with f.database(tmp_path/'normal') as left, f.database(tmp_path/'interrupted') as right:
        assert list(left.execute('SELECT * FROM blocks')) == list(right.execute('SELECT * FROM blocks'))
    f.audit(tmp_path/'interrupted',p,arithmetic)


def test_changed_protocol_and_tampered_block_rejected(tmp_path):
    p = policy(); arithmetic = f.Arithmetic(p)
    f.scan(tmp_path,p,arithmetic,max_blocks=1)
    with pytest.raises(ValueError,match='protocol mismatch'):
        f.scan(tmp_path,dict(p,domain='changed'),arithmetic)
    with f.database(tmp_path) as db:
        db.execute("UPDATE blocks SET payload=x'0001' WHERE id=0"); db.commit()
    with pytest.raises(ValueError,match='block replay differs'):
        f.audit(tmp_path,p,arithmetic)


def test_finite_draw_cap_preserves_partial_intake(tmp_path):
    p = dict(policy(),maximum_draws=2)
    arithmetic = f.Arithmetic(p)
    state = f.scan(tmp_path,p,arithmetic)
    assert state['cursor'] == 2 and state['accepted'] < p['population']
    with pytest.raises(ValueError,match='incomplete'):
        f.selection(state,p,arithmetic)


def test_primitive_duplicate_and_singular_handling(tmp_path, monkeypatch):
    p = dict(policy(),population=2,maximum_draws=5)
    arithmetic = f.Arithmetic(p)
    monkeypatch.setattr(f,'draw',lambda p,i: (1,1) if i < 3 else (2,1))
    arithmetic.A = [0]*9; arithmetic.B = [0]*13
    state = f.scan(tmp_path,p,arithmetic)
    assert state['accepted'] == state['singular'] == 2
    assert state['rejections']['duplicate'] == 2
    assert f.selection(state,p,arithmetic)['seed_inputs'] == []


def test_nonminimal_residues_replayed_on_exact_equation():
    p = policy(); arithmetic = f.Arithmetic(p)
    scaled = []
    for m,n in [(1,1),(-3,1),(5,7),(1009,101),(7,13),(32,1)]:
        row = arithmetic.record(0,m,n); A,B = arithmetic.model(m,n)
        for prime,record in zip(p['primes'],row[8]):
            a,b,k = A,B,0
            while a % prime**4 == b % prime**6 == 0:
                a //= prime**4; b //= prime**6; k += 1
            expected = arithmetic.local(a % prime,b % prime,prime)
            assert record == [expected[0],expected[1],k]
            if k:
                scaled.append(prime)
    # Force both historical lost-good-prime branches on scaled short models.
    for prime in [5,13]:
        arithmetic = f.Arithmetic(dict(p,primes=[prime]))
        arithmetic.A = [prime**4]+[0]*8
        arithmetic.B = [prime**6]+[0]*12
        arithmetic.tables = [[arithmetic.local(0,0,prime)]*prime+
                             [arithmetic.local(0,0,prime)]]
        row = arithmetic.record(0,1,1)
        expected = arithmetic.local(1,1,prime)
        assert row[8] == [[expected[0],expected[1],1]]
        assert expected[0] is not None


def test_local_trace_and_cubic_splitting_by_counting():
    for prime in [5,7,13,23]:
        for A,B in [(1,1),(2,3),(0,1),(-1,0)]:
            trace,roots,_ = f.Arithmetic.local(A,B,prime)
            if trace is None:
                assert (4*A**3+27*B*B) % prime == 0
                continue
            affine = sum((y*y-x**3-A*x-B) % prime == 0
                         for x in range(prime) for y in range(prime))
            assert trace == prime-affine
            assert roots == sum((x**3+A*x+B) % prime == 0 for x in range(prime))


def test_height_formula_and_conic_transport():
    from fractions import Fraction as F
    from det1092_v3_contract import evaluate
    arithmetic = f.Arithmetic(policy()); parent=f.read(f.PARENT)
    m,n = -17,9; row=arithmetic.record(0,m,n)
    A,B=arithmetic.model(m,n)
    assert F(A,n**8) == evaluate(parent['a_invariants'][3],F(m,n))
    assert F(B,n**12) == evaluate(parent['a_invariants'][4],F(m,n))
    j=F(6912*A**3,4*A**3+27*B*B)
    assert row[4:7] == [abs(j.numerator).bit_length(),j.denominator.bit_length(),
                         max(abs(j.numerator).bit_length(),j.denominator.bit_length())//64]
    aa,bb,cc,dd=arithmetic.matrix; tau=F(aa*m+bb*n,cc*m+dd*n)
    q=list(map(int,f.read(f.COVER)['curve_over_Q']['q_coefficients']))
    assert F(f.homogeneous(arithmetic.conic,m,n),(cc*m+dd*n)**2) == sum(v*tau**i for i,v in enumerate(q))


def test_selection_caps_unknowns_and_controls(tmp_path):
    p=policy(); arithmetic=f.Arithmetic(p)
    state=f.scan(tmp_path,p,arithmetic); result=f.selection(state,p,arithmetic)
    assert len(result['arithmetic_candidates']) <= p['arithmetic_cap']
    assert len(result['seed_inputs']) <= p['seed_cap']
    ids=[r['id'] for r in result['seed_inputs']]
    assert len(ids) == len(set(ids))
    for row in result['seed_inputs']:
        assert all(row[k] == 'UNKNOWN' for k in ['strict','selmer','kummer_extra_direction','xi'])
        assert row['id'] in {r['id'] for r in result['arithmetic_candidates']}


@pytest.mark.parametrize('rank,expected',[(18,'RETAIN'),(20,'INTERESTING'),(24,'INTERESTING'),
    (27,'INTERESTING'),(28,'DEEP'),(29,'DEEP'),(30,'AGGRESSIVE'),(31,'AGGRESSIVE'),(32,'TARGET_LOWER_BOUND')])
def test_rank_queue(rank,expected):
    assert f.queue(rank) == expected


def test_unproved_rank_cannot_enter_queue():
    for rank in [17,None,32.0,'32',True]:
        with pytest.raises(ValueError):
            f.queue(rank)


def test_real_seed_stops_before_consuming_later_candidates():
    from fractions import Fraction as F
    from det1092_funnel_worker import first_seed
    folder = f.ART/'det1092_funnel_validation_v1/seeds/conic-control'
    old, packet = f.read(folder/'m17.json'), f.read(folder/'m18.json')
    model = tuple(map(F,packet['curve']))
    points = tuple(tuple(map(F,p)) for p in packet['points'])
    consumed = []
    def candidates():
        consumed.append('inherited'); yield points[0]
        consumed.append('extra'); yield points[17]
        raise AssertionError('consumed another candidate after certified M18')
    found = first_seed(model,points[:17],old['proof'],candidates())
    assert found is not None and found[0] == points[17]
    assert found[1]['rank_lower_bound'] == 18
    assert consumed == ['inherited','extra']


def test_real_seed_gate_rejects_off_curve_point():
    from fractions import Fraction as F
    from det1092_funnel_worker import first_seed
    folder = f.ART/'det1092_funnel_validation_v1/seeds/conic-control'
    old = f.read(folder/'m17.json')
    model = tuple(map(F,old['curve']))
    points = tuple(tuple(map(F,p)) for p in old['points'])
    with pytest.raises(ValueError,match='off the exact curve'):
        first_seed(model,points,old['proof'],[(points[0][0],points[0][1]+1)])
