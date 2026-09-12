#!/usr/bin/env python3
"""Package this completed experiment; never restart arithmetic or overwrite it."""
import hashlib
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'artifacts/generated-results/elliptic-curves/constructed_class_blind_recovery_v1'
LOCAL = ROOT/'artifacts/local/elliptic-curves'
V3 = LOCAL/'constructed-class-blind-v3-v2'
COVER = LOCAL/'constructed-class-blind-cover-v1'


def read(p):
    return json.loads(p.read_text())


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def write(p, d):
    with p.open('x') as f:
        json.dump(d, f, sort_keys=True, indent=2)
        f.write('\n')


def package():
    worker = read(OUT/'cover-worker-report.json')
    terminal = read(V3/'terminal.json')
    replay = read(V3/'verified-v2.json')
    lifts = read(OUT/'v3-fixed-cover-lifts.json')
    assert worker['status'] == 'BOTH_COVERS_RECOVERED_EXACT'
    assert replay['terminal_sha256'] == sha(V3/'terminal.json')
    assert read(OUT/'independent-cover-replay.json')['status'] == 'PASS'
    assert read(OUT/'independent-quartic-replay.json')['status'] == 'PASS'
    assert read(OUT/'blind-lift-independence.json')['combined_rank_lower_bound'] == 18
    assert read(OUT/'v3-fixed-cover-replay.json')['status'] == 'PASS'
    first_failure = read(LOCAL/'constructed-class-blind-v3-v1/prepare-supervisor.json')['wall_seconds']
    prep = read(V3/'prepare-supervisor.json')['wall_seconds']
    search = read(V3/'search-supervisor.json')['wall_seconds']
    checks = read(V3/'verify-supervisor.json')['wall_seconds']+read(V3/'verify-v2-supervisor.json')['wall_seconds']
    evaluation = read(V3/'fixed-lifts-supervisor.json')['wall_seconds']
    result = {'schema': 'constructed-class-blind-comparison.v1', 'status': 'PASS',
        'family': 'MW16-05', 'parameter': '3/17', 'fixed_compaction_columns': [6, 7],
        'cover_arm': {'rational_covers_recovered': 2, 'independent_rank_over_generic': 2,
            'initial_generic_rank': 16, 'combined_rank_lower_bound': 18,
            'arithmetic_wall_seconds_charged': worker['arithmetic_wall_seconds_charged'],
            'quartic_search_wall_seconds': worker['local_quartic_search_wall_seconds_sum'],
            'raw_coefficient_bits': [2464, 1766], 'reduced_degree4_coefficient_bits': [39, 44],
            'quartic_parameters': ['1/2', '3/10'], 'quartic_search_bound': 1000,
            'direct_degree4_search_bounds_without_hit': [1000, 100000]},
        'v3_arm': {'initial_generic_rank': 16, 'final_rank_lower_bound': 22,
            'charts': terminal['charts'], 'completed_charts': replay['completed_charts'],
            'chart_bound': 125000, 'stop_reason': terminal['stop_reason'],
            'preparation_wall_seconds': prep, 'retained_failed_preparation_wall_seconds': first_failure,
            'full_search_wall_seconds': search, 'backend_search_wall_seconds': replay['backend_search_wall_seconds'],
            'landscape_seconds': sum(s['landscape_seconds'] for s in terminal['stages']),
            'map_seconds': sum(a['wall_seconds'] for s in terminal['stages'] for a in s['attempts']),
            'replay_wall_seconds_including_serialization_failure': checks,
            'fixed_class_evaluation_wall_seconds': evaluation,
            'fixed_class_independent_replay_seconds': read(OUT/'v3-fixed-cover-replay.json')['wall_seconds'],
            'gain_times': [{k:v for k,v in g.items() if k != 'point'} for g in terminal['gain_times']],
            'first_two_directions_seconds_including_preparation': first_failure+prep+terminal['gain_times'][1]['search_arm_seconds'],
            'rank22_prefix_seconds_including_preparation': first_failure+prep+terminal['gain_times'][-1]['search_arm_seconds'],
            'fixed_covers_recovered_from_new_subgroup': 2,
            'fixed_cover_earliest_prefix_ranks': [r['earliest_certified_prefix_rank'] for r in lifts['cases']]},
        'comparison_boundary': 'The cover arm starts from pre-existing fixed classes; their upstream construction cost is excluded. V3 starts only from equation and the same generic subgroup. Remote Magma versus local Sage/PARI, independent preparation choices and implementation time prevent an end-to-end or general speed-superiority claim. This is a retrospective fixed-fibre recovery experiment, not a parameter-selection mechanism, new rank record, exact-rank proof or global class/Selmer computation.',
        'bindings': {str(p.relative_to(ROOT)): sha(p) for p in [OUT/'protocol.json', OUT/'covers.json', OUT/'generic.json',
            OUT/'cover-worker-report.json', OUT/'cover-worker-witness.json', OUT/'independent-cover-replay.json',
            OUT/'independent-quartic-replay.json', OUT/'blind-lift-independence.json', OUT/'v3-fixed-cover-lifts.json',
            OUT/'v3-fixed-cover-replay.json', V3/'protocol.json', V3/'terminal.json', V3/'verified-v2.json']}}
    write(OUT/'results.json', result)
    paths = set()
    for base in [LOCAL/'constructed-class-blind-v3-v1', V3, COVER]:
        paths.update(p for p in base.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix != '.pdf')
    paths.update(p for p in OUT.glob('*.json'))
    own = ['prepare_constructed_class_blind.py', 'run_constructed_class_v3.py', 'replay_constructed_class_v3.py',
           'evaluate_constructed_class_v3.py', 'verify_constructed_class_blind.py', 'replay_constructed_class_quartic.py',
           'certify_constructed_class_blind.py', 'test_constructed_class_blind.py', 'package_constructed_class_blind.py']
    paths.update(Path(__file__).with_name(n) for n in own)
    paths.update(Path(__file__).parent.glob('blind_constructed_cover_*.py'))
    for name in read(V3/'protocol.json')['sources']:
        paths.add(ROOT/name)
    for name in read(V3/'protocol.json')['inputs']:
        paths.add(ROOT/name)
    for p in V3.glob('epoch-*/chart-*.json'):
        paths.update(ROOT/n for n in read(p)['search']['source_hashes'])
    for name in ['rank_jump_reference_strict_class_construction_inputs_v1.json', 'rank_jump_constructed_class_compaction_v1.json']:
        paths.add(OUT.parent/name)
    manifest = {str(p.relative_to(ROOT)): {'sha256': sha(p), 'bytes': p.stat().st_size} for p in sorted(paths)}
    target = OUT/'evidence.zip'
    with zipfile.ZipFile(target, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for p in sorted(paths):
            info = zipfile.ZipInfo(str(p.relative_to(ROOT)), date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info, p.read_bytes())
    with zipfile.ZipFile(target) as archive:
        assert set(archive.namelist()) == set(manifest)
        for name, record in manifest.items():
            assert hashlib.sha256(archive.read(name)).hexdigest() == record['sha256']
    write(OUT/'evidence-manifest.json', {'schema': 'constructed-class-blind-evidence.v1',
        'archive_sha256': sha(target), 'archive_bytes': target.stat().st_size,
        'files': manifest, 'file_count': len(manifest),
        'restore': 'Extract only into an empty replay checkout, with paths rooted at research/. Never overwrite an active research tree.',
        'boundary': 'Retained arithmetic inputs, programs, raw responses, failures, maps, searches and proof records. Hashes certify byte integrity; mathematical replay is separate. The referenced public paper PDF is omitted.'})
    print('PACKAGED', len(manifest), target.stat().st_size, sha(target))


if __name__ == '__main__':
    package()
