import ast
import importlib.util
import inspect
from fractions import Fraction
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
SCRIPT=ROOT/'elliptic-curves/cas/run_curve302_seeded_v3_amplifier.py'
spec=importlib.util.spec_from_file_location('curve302_seeded_v3_amplifier',SCRIPT)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def test_two_frozen_single_seeds_only():
    assert m.SEEDS==('recovered-strict-02','recovered-strict-03')
    assert len(set(m.SEEDS))==2


def test_rank_boundary_is_m18_to_m31():
    assert m.INITIAL_RANK==18
    assert m.TARGET_RANK==31


def test_evidence_namespace_is_separate():
    assert m.D.name=='curve302-seeded-v3-amplifier-v1'
    assert 'det1092-v4' not in str(m.D)
    assert m.D != m.V3


def test_controller_import_has_no_sage_dependency():
    assert callable(m.launch)
    assert callable(m.status)
    assert callable(m.worker)
    assert callable(m.ordered_seed_state)


def test_oracle_inputs_are_explicit():
    assert m.VISIBILITY.name=='curve302_residual_visibility_geometry_v1.json'
    assert m.M24.name=='curve302_recovered_followup_wave_03_mod2_v1.json'
    assert m.ORBITS.name=='curve302_parent_degree2_multisection_orbits_v1.tsv'


def test_ordered_seed_state_uses_existing_normalized_trust_boundary():
    tree=ast.parse(inspect.getsource(m.ordered_seed_state))
    calls=[]
    for node in ast.walk(tree):
        if not isinstance(node,ast.Call):
            continue
        if isinstance(node.func,ast.Name):
            calls.append(node.func.id)
        elif isinstance(node.func,ast.Attribute):
            calls.append(node.func.attr)
    assert 'certified_state' in calls
    assert 'assert_basis' in calls
    assert 'raw_state' not in calls
    source=inspect.getsource(m.ordered_seed_state)
    assert 'tuple(state.basis)' not in source
    assert 'reductions.points' not in source


def test_basis_comparison_normalizes_storage_representation():
    # MWState/research_runtime stores exact rationals canonically as strings;
    # callers frequently hold the same values as Fraction/Sage rationals.
    import sys
    sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
    from v3_warm_support import point_tuple
    stored=(('1/2','-7/3'),('123','5/11'))
    caller=((Fraction(1,2),Fraction(-7,3)),(Fraction(123),Fraction(5,11)))
    assert point_tuple(stored)==point_tuple(caller)
