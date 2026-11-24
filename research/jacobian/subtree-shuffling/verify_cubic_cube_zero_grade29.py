#!/usr/bin/env python3
"""Replay the pinned prefix and a selected grade-29 root-relation certificate."""
if not __debug__:
    raise SystemExit('Exact checks require ordinary Python execution.')

from argparse import ArgumentParser
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import gzip
import json
import sys
import time


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load(name, path):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    start = time.monotonic()
    folder = args.repo_root / 'research/jacobian/subtree-shuffling'
    ledger_path = folder / 'cubic_cube_zero_shuffle_source_ledger.json'
    ledger = json.loads(ledger_path.read_text())
    for entry in ledger['inputs']:
        base = args.repo_root if entry['scope'] == 'repository' else folder
        path = base / entry['path']
        require(sha256(path.read_bytes()).hexdigest() == entry['sha256'],
                'Changed pinned prefix source: ' + entry['path'])
    prefix_path = folder / 'cubic_cube_zero_shuffle_certificates.json.gz'
    prefix = json.loads(gzip.decompress(prefix_path.read_bytes()))
    cert = json.loads(gzip.decompress(args.input.read_bytes()))
    require(cert['schema_version'] == 1 and cert['status'] == 'FULL', 'Incomplete certificate')
    require(cert['grade'] == 29 and cert['prime'] == 1000000009, 'Wrong problem or field')
    require(cert['columns'] == cert['rank'] == len(cert['relations']) == 46928,
            'Wrong certificate size')
    require(cert['prefix_certificate_sha256'] == sha256(prefix_path.read_bytes()).hexdigest(),
            'Wrong prefix certificate')
    require([c['grade'] for c in prefix['exact_grades']] == list(range(1, 19)),
            'Incomplete exact prefix')
    require([c['grade'] for c in prefix['modular_zero_grades']] == list(range(19, 29)),
            'Incomplete zero prefix')

    s = load('grade29_replay_kernel', folder / 'verify_subtree_shuffle.py')
    v = load('grade29_replay_prefix', folder / 'verify_cubic_cube_zero_shuffle.py')
    q = s.ShuffleQuotient(3, 3)
    exact, modular_pivots = [], [{}]
    for item in prefix['exact_grades']:
        row, pivots = v.verify_exact(q, item, s)
        exact.append(row)
        modular_pivots.append(pivots)
        print('EXACT', json.dumps(row), flush=True)

    class ModularQuotient(s.ShuffleQuotient):
        @lru_cache(maxsize=250000)
        def lower_shuffle(self, side, bottom):
            if not side:
                return {bottom: 1}
            grade = self.grade[bottom] + sum(self.grade[i] for i in side) + len(side) // 2
            require(grade < len(self.basis), 'Using an unproved lower grade')
            if grade >= 19:
                require(not self.basis[grade], 'Zero shortcut on a nonzero piece')
                return {}
            return v.nf_mod(self.root_shuffles(side, bottom), self.pivots[grade])

    mq = ModularQuotient(3, 3)
    mq.grade = list(q.grade)
    mq.tree_ids = dict(q.tree_ids)
    mq.basis = list(q.basis)
    mq.pivots = modular_pivots
    zeros = []
    for item in [*prefix['modular_zero_grades'], cert]:
        k = item['grade']
        require(item['prime'] == v.P and k == len(mq.basis), 'Nonconsecutive zero piece')
        children = v.make_children(mq, k)
        require(len(item['relations']) == len(children), 'Wrong number of relations')
        pivots = {}
        last_log = time.monotonic()
        for side, bottom in item['relations']:
            side = tuple(side)
            v.validate_substitution(mq, k, side, bottom)
            row = mq.root_shuffles(side, bottom)
            require(all(i in children for i in row), 'Relation outside current prequotient')
            require(v.insert_mod(row, pivots), 'Dependent selected relation')
            if k == 29 and time.monotonic() - last_log >= 15:
                print('GRADE29', json.dumps({'rank': len(pivots),
                      'columns': len(children), 'seconds': time.monotonic() - start}), flush=True)
                last_log = time.monotonic()
        require(len(pivots) == len(children), 'Non-full zero-piece rank')
        mq.pivots.append({})
        mq.basis.append(())
        row = {'grade': k, 'columns': len(children),
               'independent_valid_relations': len(pivots), 'dimension': 0}
        zeros.append(row)
        print('ZERO', json.dumps(row), flush=True)

    inherited_controls = v.controls()
    result = {'status': 'PASS', 'arithmetic': 'QQ and F_1000000009; standard library only',
              'exact_prefix': exact, 'modular_zero_grades': zeros,
              'conclusion': 'A_(3,3),29 = 0 over every characteristic-zero field',
              'whole_algebra_dimension_claimed': False,
              'input_sha256': sha256(args.input.read_bytes()).hexdigest(),
              'prefix_certificate_sha256': sha256(prefix_path.read_bytes()).hexdigest(),
              'prefix_ledger_sha256': sha256(ledger_path.read_bytes()).hexdigest(),
              'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
              'inherited_controls': inherited_controls, 'seconds': time.monotonic() - start}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print('PASS', json.dumps({'grade': 29, 'columns': 46928, 'seconds': result['seconds']}), flush=True)


if __name__ == '__main__':
    main()
