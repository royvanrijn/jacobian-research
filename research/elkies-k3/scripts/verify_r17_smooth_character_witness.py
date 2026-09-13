#!/usr/bin/env python3
"""Independent, standard-library replay of two fixed smooth-character atlases.

Distinct irreducible quadratic polynomial factors give independent valuations
in Q(u)^*/Q(u)^*2, regardless of their rational scalar factors. This checks
injectivity and absence of internal character triples, plus attachment to the
retained priority tables. It does not construct sections, prove completeness of
the lattice census, transport old-base characters, or treat singular pencils.

The default reads only the compact packet. --check-source-projection also reads
the four original files (about 692 MB); --export-inputs writes a NEW packet from
those exact files. Neither mode reconstructs missing inputs or launches a CAS.
"""

import argparse
import csv
from fractions import Fraction
import gzip
import hashlib
import json
from math import gcd, isqrt, lcm
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[2]
RESULTS = 'artifacts/generated-results/'
PACKET = ROOT / (RESULTS + 'elkies-k3-r17-smooth-character-replay-inputs-v1.json.gz')
PACKET_SHA256 = 'fa46f7e1e344e4b9003a4b686783f88c023509ab98f5abef44cb6143f36a0450'
SCHEMA = 'elkies-k3.r17-smooth-character-replay-inputs.v1'
SPECS = [
    dict(chart='11952', count=39147,
         atlas=RESULTS+'elkies-k3-r17-norm12-11952-alternate-bisections-full-v1.json',
         atlas_sha256='fbf979bfe7d92528405c62330a80dbfd7742dc27c41a0426dcd4014f6865c8ce',
         atlas_schema='elkies-k3.r17-norm12-11952-alternate-bisections-full.v1',
         vector_field='alternate_rank17_w',
         priority=RESULTS+'elkies-k3-r17-norm12-11952-alternate-bisection-priority-v1.tsv',
         priority_sha256='d052174f241d9ea786570138bda6b32498c3526eaf0f3e8ee99f65f7e8925256',
         priority_vector_field='historical_alternate_w'),
    dict(chart='103b2', count=39120,
         atlas=RESULTS+'elkies-k3-r17-norm12-103b2-bisections-full-v1.json',
         atlas_sha256='faa98745f9bc8fb304493533ca93ecbb4b36606981ed2899282de381305fdeed',
         atlas_schema='elkies-k3.r17-norm12-103b2-hidden-bisections-full.v1',
         vector_field='direct_hidden_w',
         priority=RESULTS+'elkies-k3-r17-norm12-103b2-bisection-priority-v1.tsv',
         priority_sha256='c43750377c5772b3fabd4ed0bcd5e9950f6e3aa1c3fba3eee1ef41a608922898',
         priority_vector_field='direct_hidden_w'),
]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            value.update(block)
    return value.hexdigest()


def quadratic_atom(numerator, denominator):
    """Primitive polynomial, up to a nonzero rational scalar; fail if reducible."""
    require(isinstance(numerator, list) and len(numerator) == 3,
            'expected exactly three quadratic coefficients')
    require(isinstance(denominator, list) and len(denominator) == 1,
            'nonconstant branch denominator requires a different proof')
    require(all(type(x) is str for x in numerator + denominator),
            'coefficients must be exact rational strings')
    divisor = Fraction(denominator[0])
    require(divisor != 0, 'zero branch denominator')
    q = [Fraction(x)/divisor for x in numerator]
    require(q[2] != 0, 'branch polynomial is not quadratic')
    scale = lcm(*(x.denominator for x in q))
    integers = [int(x*scale) for x in q]
    content = gcd(*integers)
    if integers[2] < 0:
        content = -content
    a, b, c = (x//content for x in integers)  # low to high
    discriminant = b*b-4*a*c
    require(discriminant < 0 or isqrt(discriminant)**2 != discriminant,
            'quadratic is reducible or has a repeated root')
    return a, b, c


def validate_identity(mask, word):
    require(type(mask) is int and 0 <= mask < 2**17, 'invalid orbit mask')
    require(isinstance(word, list) and len(word) == 17
            and all(type(x) is int for x in word), 'invalid fixed-frame word')


def verify_atlas(atlas, expected_count):
    rows, priority = atlas['rows'], atlas['priority']
    require(len(rows) == len(priority) == expected_count, 'incomplete atlas or priority table')
    expected = {}
    for mask, word in priority:
        validate_identity(mask, word)
        require(mask not in expected, 'duplicate priority mask')
        expected[mask] = word
    masks, labels, atoms = set(), set(), set()
    for mask, label, word, numerator, denominator in rows:
        validate_identity(mask, word)
        require(isinstance(label, str) and label, 'missing branch label')
        require(mask not in masks and label not in labels, 'duplicate branch identity')
        require(mask in expected and word == expected[mask], 'priority/frame attachment mismatch')
        atom = quadratic_atom(numerator, denominator)
        require(atom not in atoms,
                'repeated polynomial support: full rational squareclasses must be compared')
        masks.add(mask)
        labels.add(label)
        atoms.add(atom)
    require(masks == expected.keys(), 'priority coverage mismatch')
    return dict(count=len(rows), distinct_irreducible_quadratic_atoms=len(atoms),
                priority_attachments=len(masks), internal_three_character_relations=0)


def verify(packet):
    require(packet['schema'] == SCHEMA, 'wrong packet schema')
    require(len(packet['atlases']) == len(SPECS), 'wrong atlas collection')
    results = {}
    for spec, atlas in zip(SPECS, packet['atlases']):
        require(atlas['source'] == spec, 'source identity or schema changed')
        results[spec['chart']] = verify_atlas(atlas, spec['count'])
    return results


def project_source(spec, root=ROOT):
    """Copy raw coefficients and both frame words; no derived math or status flags."""
    for kind in ('atlas', 'priority'):
        path = root / spec[kind]
        require(path.is_file(), f'missing original input: {path}')
        require(digest(path) == spec[kind+'_sha256'], f'changed original input: {path}')
    with (root / spec['atlas']).open() as stream:
        original = json.load(stream)
    require(original['artifact_schema'] == spec['atlas_schema'], 'wrong atlas schema')
    require(original['base_parameter'] == 'u', 'wrong base parameter')
    rows = [[r['lattice_orbit_mask'], r['label'], r[spec['vector_field']],
             r['branch']['numerator_coefficients'], r['branch']['denominator_coefficients']]
            for r in original['bisections']]
    del original
    with (root / spec['priority']).open(newline='') as stream:
        priority = [[int(r['orbit_mask']), [int(x) for x in r[spec['priority_vector_field']].split()]]
                    for r in csv.DictReader(stream, delimiter='\t')]
    return dict(source=spec, rows=rows, priority=priority)


def check_source_projection(packet, root=ROOT):
    for spec, atlas in zip(SPECS, packet['atlases']):
        require(project_source(spec, root) == atlas, f"source projection mismatch: {spec['chart']}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--check-source-projection', action='store_true')
    mode.add_argument('--export-inputs', type=Path,
                      help='write a new deterministic gzip packet from the pinned originals')
    args = parser.parse_args()
    started = time.monotonic()
    if args.export_inputs:
        require(not args.export_inputs.exists(), 'refusing to overwrite a retained packet')
        packet = dict(schema=SCHEMA, atlases=[project_source(s) for s in SPECS])
    else:
        require(digest(PACKET) == PACKET_SHA256, 'retained packet changed')
        with gzip.open(PACKET, 'rt') as stream:
            packet = json.load(stream)
    results = verify(packet)
    if args.check_source_projection:
        check_source_projection(packet)
    if args.export_inputs:
        raw = (json.dumps(packet, separators=(',', ':'), ensure_ascii=True)+'\n').encode()
        # No filename or timestamp in gzip headers; the JSON content is deterministic.
        with args.export_inputs.open('xb') as stream:
            with gzip.GzipFile(fileobj=stream, mode='wb', filename='', mtime=0,
                               compresslevel=6) as zipped:
                zipped.write(raw)
    path = args.export_inputs or PACKET
    print(json.dumps(dict(status='PASS_RETAINED_SMOOTH_CHARACTER_WITNESS', atlases=results,
                          packet_sha256=digest(path), packet_bytes=path.stat().st_size,
                          source_projection_checked=bool(args.export_inputs or args.check_source_projection),
                          elapsed_seconds=round(time.monotonic()-started, 3)), sort_keys=True))


if __name__ == '__main__':
    main()
