#!/usr/bin/env python3
"""Publish compact views of the frozen experiment; never rerun a search."""
import csv
import hashlib
import io
import json
import subprocess
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'artifacts/generated-results/elliptic-curves/rank_triangle_v1'


def read(p):
    return json.loads(p.read_text())


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def emit(path, content):
    data = content if isinstance(content, bytes) else content.encode()
    if path.exists():
        assert path.read_bytes() == data, path
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)


def save(name, data):
    emit(OUT / name, json.dumps(data, sort_keys=True, indent=2) + '\n')


def main():
    verified = read(OUT / 'verified.json')
    assert verified['status'] == 'PASS_BOUNDED_TRIANGLE_REPLAY'
    for path, digest in verified['bindings'].items():
        assert sha(Path(path)) == digest, path
    for name, key in [('geometry-protocol-v2.json', 'sha256'),
                      ('native-carrier-protocol.json', 'sources'), ('strict.json', 'bindings')]:
        for path, digest in read(OUT / name)[key].items():
            assert sha(Path(path)) == digest, path
    geometry = read(OUT / 'geometry.json')
    native = read(OUT / 'native-carriers.json')
    strict = read(OUT / 'strict.json')
    classes = {r['target']: r for r in strict['individual_exceptional_classes']}
    assert len(geometry['cover_equivalence_classes']) == 784
    target_rows = []
    for row in geometry['rows']:
        covers = row['covers'] + ([c for c in native['rows'] if c['target'] == row['target']] if row['fibre'] == '11952' else [])
        best = [c for c in covers if c['normalization_genus'] == 1]
        assert len(best) == 1
        target_rows.append({'fibre': row['fibre'], 'target': row['target'],
                            'source_point': row['source_point'], 'parent_point': row['parent_point'],
                            'source_to_parent_isomorphism': row['fibre_isomorphism'],
                            'genus_histogram': dict(Counter(c['normalization_genus'] for c in covers)),
                            'best_tested_carrier': best[0],
                            'global_minimum_genus': 'UNKNOWN',
                            'degree_minimum': 2,
                            'degree_lower_bound_reason': 'A degree-one curve is a generic section; the certified target is outside the rational span of the full-rank specialized M17.',
                            'strict_dictionary': classes[row['target']] if row['fibre'] == '302' else None})
    save('targets.json', target_rows)
    panel = read(OUT / 'panel.json')['rows']
    arithmetic = {r['id']: r for r in read(OUT / 'arithmetic.json')['results']}
    arithmetic['302'] = read(OUT / 'arithmetic/302-cached-support/result.json')
    buf = io.StringIO()
    names = ['id', 'family', 'parameter', 'rank_lower_bound', 'jump_lower_bound', 'calls',
             'bk_offset', 'field_discriminant', 'signature', 'root_number', 'conductor',
             'g_upper', 'rank_upper', 'status']
    writer = csv.DictWriter(buf, fieldnames=names); writer.writeheader()
    for row in panel:
        rec = arithmetic[row['id']]
        merged = dict(row, **{k: rec.get(k, 'UNKNOWN') for k in names[6:]})
        writer.writerow({k: 'UNKNOWN' if merged.get(k) is None else merged.get(k, '') for k in names})
    emit(OUT / 'panel.csv', buf.getvalue())
    # Preserve the small parent input otherwise located only in the live runtime.
    parent = ROOT / 'artifacts/local/elliptic-curves/broad-rank-v1/runtime/research/broad-inputs/parents/11952.json'
    emit(OUT / 'inputs/11952-parent.json', parent.read_bytes())
    save('inputs/README.json', {'source': str(parent), 'sha256': sha(parent),
                               'usage': 'Byte-for-byte frozen input. Original producer paths and hashes are retained in geometry-protocol-v2.json.'})
    sources = list(Path(__file__).parent.glob('rank_triangle*')) + [Path(__file__).with_name('verify_rank_triangle.sage'), Path(__file__)]
    inputs = [p for p in OUT.rglob('*') if p.is_file() and p.name != 'manifest.json']
    save('manifest.json', {'status': 'FROZEN_BOUNDED_TRIANGLE',
                           'software': subprocess.check_output(['sage', '--version'], text=True).strip(),
                           'files': {str(p.relative_to(ROOT)): sha(p) for p in sorted(set(sources + inputs))},
                           'boundary': 'Finite fixed-target atlas and equation panel, not a complete multisection census, exact rank result, or independently implemented full descent.'})
    print('Published 22 target records, 18 panel rows, and hashed manifest.')


if __name__ == '__main__':
    main()
