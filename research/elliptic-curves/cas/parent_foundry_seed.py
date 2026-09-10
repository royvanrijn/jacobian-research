"""Prepare generic-only or seeded parent searches without changing old checkers.

The search and replay engine is the unchanged V3 implementation. This preparer
permits a certified generic-rank input and pins its own source in each protocol.
The launched v2 snapshot has the equivalent opt-in gate inside its frozen V3
copy; this extraction preserves historical checker-source pins in the checkout.
"""
from pathlib import Path
import shutil
from run_complement_seed_v3 import (CAS, ROOT, load, read, require, sha,
                                    checkpoint, MAPPERS)


def freeze(folder, prep):
    PREP = BANK = prep
    import sage.version
    from sage.all import pari
    import numpy
    old = load('productive_sources', CAS/'adaptive_visibility_cascade_v3.sage')
    bank = read(BANK/'anchor-bank.json')
    require(bank['status'] in ('COMPLETE_FROZEN_PRODUCTIVE_SUBSET', 'COMPLETE_FROZEN_PAIRWISE_PARENT_SUBSET', 'COMPLETE_FROZEN_COMPLEMENT_PARENT_SUBSET'), 'exact parent subset required')
    require(max(sum(r['norm'] == q for r in bank['rows']) for q in bank['shells']) <= 16, 'all-parent subset exceeds existing V3 shell quota')
    origin = 'gains' if bank['status'] == 'COMPLETE_FROZEN_PRODUCTIVE_SUBSET' else 'derivation'
    require(bank['protocol_sha256'] == sha(BANK/'protocol.json') and
            bank[origin+'_sha256'] == sha(BANK/(origin+'.json')), 'productive bank seal differs')
    seed_names = sorted(PREP.glob('seed-M*.json'))
    require(len(seed_names) == 1, 'exactly one prepared seed required')
    seed_path = seed_names[0]
    seed = read(seed_path)
    initial_rank = seed['rank_lower_bound']
    minimum_rank = seed['generic_rank']
    require(seed['generic_rank'] in (16, 17) and minimum_rank <= initial_rank < 32, 'native MW16/R17 seed below target required')
    prepared = read(PREP/'prepared.json')
    require(prepared['files'][seed_path.name] == sha(seed_path), 'seed changed')
    require(all(sha(PREP/n) == h for n,h in prepared['files'].items()), 'prepared native inputs changed')
    names = ('parent_foundry_seed.py', 'run_complement_seed_v3.py', 'lean_preconditioned_map_receipts.py', 'lean_preconditioned_map_worker.py',
             'lean_preconditioned_full_pari_mapping.sage',
             'lean_factor_free_pari_mapping.sage', 'visibility_complement_subset.py',
             'visibility_lattice_fast.py', 'pointed_box_equivalence.py', 'future_point_admission.py',
             'memory_rank_certificate.py', 'v3_warm_engine.py', 'v3_warm_support.py',
             'prepare_extended20_mw16_pari_batch.sage')
    inputs = [seed_path, PREP/'prepared.json', BANK/'anchor-bank.json',
              BANK/'protocol.json', BANK/(origin+'.json'), PREP/'generic-cvp-proofs.json']
    protocol = {'schema': 'prepared-parent-complement-boxes.v1', 'initial_rank': initial_rank,
        'anchor_count': len(bank['rows']),
        'family': read(PREP/'protocol.json')['family'],
        'parameter': read(PREP/'protocol.json')['parameter'],
        'target_rank': 32, 'generic_rank': seed['generic_rank'], 'scaled_shells': bank['shells'],
        'anchors_per_shell': 16, 'canonical_per_shell': 25,
        'exact_cvp_node_limit': 2000000, 'prime_bound': 1000,
        'height': 125000, 'seconds_per_chart': 10, 'max_charts': 100,
        'seconds_per_map': 5, 'map_rss_bytes': 1024**3,
        'map_python': shutil.which('sage'), 'map_python_sha256': sha(Path(shutil.which('sage'))),
        'max_epochs': 33 - initial_rank, 'mapping_order': list(MAPPERS),
        'search_wall_limit_seconds': 1800, 'replay_wall_limit_seconds': 1800,
        'rss_limit_bytes': 3*1024**3, 'maximum_workers': 1,
        'software': {'sage': sage.version.version, 'pari': str(pari.version()), 'numpy': numpy.__version__},
        'gp_sha256': sha(Path('/usr/bin/gp')),
        'inputs': {str(p.relative_to(ROOT)): sha(p) for p in inputs},
        'sources': {**old.sources(), **{str((CAS/n).relative_to(ROOT)): sha(CAS/n) for n in names}},
        'selection': 'V3 diversified extension shortlist over every frozen exact parent class; no per-shell anchor truncation. '
            'Full extension scoring, exact CVP. Factor-free chart profiles '
            'determine centre order. For each centre try factor-free-preconditioned full PARI minimization, then factor-free '
            'coordinates, skipping identical height boxes under signed coordinate permutations including infinity. Each map has its own5-second '
            'and1-GiB worker bound; resource-limited maps have explicit incomplete receipts. '
            'Rebuild immediately after the first standalone certified gain.',
        'scope': 'New directions beyond the supplied certified native MW16/R17 seed. Exact parent-class '
            'selection follows the bound preparation derivation; productivity is not assumed. No withheld/higher-rank point or visibility '
            'labels enter execution. The curve90 masked pointwise nulls are not passed controls '
            'or rank upper bounds. Map replay verifies exact rational identities without rerunning minimization; '
            'no canonical minimality claim. Resource limits imply no point absence. '
            'This is a bounded discovery experiment, not a success prediction.'}
    folder.mkdir(exist_ok=False)
    checkpoint(folder/'seed.json', seed)
    checkpoint(folder/'bank.json', bank)
    protocol['inputs'].update({str((folder/n).relative_to(ROOT)): sha(folder/n)
                               for n in ('seed.json', 'bank.json')})
    checkpoint(folder/'protocol.json', protocol)
    return protocol


