import ast
import importlib.util
import inspect
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
SCRIPT=ROOT/'elliptic-curves/cas/analyze_curve302_seed_pair_arithmetic.py'
spec=importlib.util.spec_from_file_location('curve302_seed_pair_arithmetic',SCRIPT)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def test_frozen_pair_and_bounds():
    assert m.SEEDS==('recovered-strict-02','recovered-strict-03')
    assert m.PRIME_BOUND==1000
    assert m.OUT.name=='curve302_seed_pair_arithmetic_v1.json'


def test_no_point_search_or_cascade_execution():
    source=SCRIPT.read_text()
    forbidden=('PointedQuarticSearch','run_search(','backend.execute(','adaptive_visibility_cascade_v3')
    assert not any(token in source for token in forbidden)


def test_joint_rank19_is_required_before_distinct_claim():
    source=inspect.getsource(m.build)
    assert "require(joint_rank == 19" in source
    assert "'distinct_over_generic_rational_span': True" in source
    assert "'distinct_mod2_cosets_over_generic_subgroup': True" in source
    assert source.index("require(joint_rank == 19") < source.index("'distinct_over_generic_rational_span': True")


def test_exact_mod2_profile_is_frozen():
    source=inspect.getsource(m.build)
    for expected in (
        "'generic_M17': 17",
        "'M17_plus_recovered_strict_02': 18",
        "'M17_plus_recovered_strict_03': 18",
        "'M17_plus_both': 19",
    ):
        assert expected in source


def test_crop_rank_is_exact_gf2_not_numeric():
    tree=ast.parse(inspect.getsource(m.crop_rank))
    calls=[]
    for node in ast.walk(tree):
        if isinstance(node,ast.Call):
            if isinstance(node.func,ast.Name): calls.append(node.func.id)
            elif isinstance(node.func,ast.Attribute): calls.append(node.func.attr)
    assert 'gf2_rank' in calls


def test_completed_amplifier_evidence_is_required():
    source=inspect.getsource(m.build)
    assert "COMPLETE_TWO_SEED_AMPLIFIER" in source
    assert "PASS_INDEPENDENT_SEEDED_V3_REPLAY" in source
    assert "row['rank_lower_bound'] == 31" in source


def test_json_native_roundtrip_removes_tuple_list_false_mismatches():
    value={'rows':[(1,0),(0,1)],'nested':{'x':('1','2')}}
    native=m.json_native(value)
    assert native==json.loads(json.dumps(native))
    assert native=={'rows':[[1,0],[0,1]],'nested':{'x':['1','2']}}


def test_build_returns_json_native_before_immutable_replay():
    source=inspect.getsource(m.build)
    assert 'return json_native(result)' in source
