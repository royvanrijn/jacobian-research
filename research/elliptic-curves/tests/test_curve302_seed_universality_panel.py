import ast
import importlib.util
import inspect
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
SCRIPT=ROOT/'elliptic-curves/cas/run_curve302_seed_universality_panel.py'
spec=importlib.util.spec_from_file_location('curve302_seed_universality_panel',SCRIPT)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def test_panel_is_exact_remaining_twelve():
    assert len(m.EXPECTED_DIRECTIONS)==14
    assert m.KNOWN_WINNERS==('recovered-strict-02','recovered-strict-03')
    assert set(m.KNOWN_WINNERS)<set(m.EXPECTED_DIRECTIONS)
    order=m.panel_order()
    assert len(order)==12
    assert set(order)==set(m.EXPECTED_DIRECTIONS)-set(m.KNOWN_WINNERS)
    assert not set(order)&set(m.KNOWN_WINNERS)


def test_order_comes_from_frozen_census_not_outcomes():
    source=inspect.getsource(m.panel_order)
    assert "single_seed_ranking" in source
    assert "rank_lower_bound" not in source
    assert "seeded-verified" not in source
    assert m.CLOSURE.name=='curve302_unlock_seed_closure_v2.json'


def test_separate_evidence_namespace():
    assert m.D.name=='curve302-seed-universality-panel-v1'
    assert 'seeded-v3-amplifier-v1' not in str(m.D)


def test_foreground_zero_chart_preflight_is_required():
    source=inspect.getsource(m.preflight_all)
    assert "ctx.state.rank!=18" in source
    assert "'charts':0" in source
    launch=inspect.getsource(m.launch)
    assert "preflight.json" in launch
    assert "run foreground preflight" in launch


def test_worker_reuses_unchanged_base_search_and_replay():
    source=inspect.getsource(m.worker)
    assert 'searcher.run_search(ctx)' in source
    assert 'base.verify_case(ctx)' in source
    assert "rank_lower_bound" in source
    # No outcome-conditioned replacement/refill roster.
    assert 'append(' not in source
    assert 'panel_order()' in source


def test_controller_import_has_no_sage_dependency():
    assert callable(m.preflight_all)
    assert callable(m.launch)
    assert callable(m.status)
    assert callable(m.worker)


def test_claim_boundary_is_retrospective():
    source=SCRIPT.read_text()
    assert 'Retrospective known-seed amplification panel' in source
    assert 'not a prospective selector or new rank claim' in source
