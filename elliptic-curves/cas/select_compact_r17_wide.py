#!/usr/bin/env python3
"""Balanced H4096 selector on all six compact R17 families, using retained traces."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
from research_runtime.store import checkpoint, digest
from research_runtime.supervisor import capture, Limits

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT/'elliptic-curves/cas'
PARENT = ROOT/'artifacts/local/elliptic-curves/compact-six-r17-h1024-v2'
DIRECTORY = ROOT/'artifacts/local/elliptic-curves/compact-six-r17-h4096-v1'
BINARY = ROOT/'artifacts/local/elliptic-curves/compact-r17-wide-v1/scanner'
YIELD = ROOT/'artifacts/generated-results/elliptic-curves/compact_atlas_new_curves_v1.json'

def sources():
    return {str(p.relative_to(ROOT)): cert.hashed(p) for p in (
        Path(__file__).resolve(), spec.ATLAS, Path(spec.__file__).resolve(), Path(cert.__file__).resolve(),
        YIELD, PARENT/'trace-table-check.json', CAS/'newfamily/scan_rational_nagao_tables.cpp',
        CAS/'research_runtime/store.py', CAS/'research_runtime/supervisor.py')}

def prepare(directory):
    if (directory/'protocol.json').exists():
        raise FileExistsError('preserve wide R17 protocol')
    if max(r['rank_certificate']['rank_lower_bound'] for r in cert.read(YIELD)['curves']) < 25:
        raise ArithmeticError('prospective positive-yield gate failed')
    families = [f['family'] for f in cert.read(spec.ATLAS)['families']]
    tables = {f: {str(sign): {'path': str((PARENT/f/f'tables-{sign}.txt').relative_to(ROOT)),
                             'sha256': cert.hashed(PARENT/f/f'tables-{sign}.txt')}
                  for sign in (-1, 1)} for f in families}
    checkpoint(directory/'protocol.json', {
        'schema': 'elliptic-curves.compact-six-r17-wide-selection.v1', 'sources': sources(),
        'scanner_binary_sha256': cert.hashed(BINARY), 'families': families,
        'height': 4096, 'prime_bound': 4093, 'retained_per_family': 128, 'finalists_per_family': 4,
        'trace_tables': tables,
        'score': 'All 562 primes 5..4093: sum round(1e12*(2-a_p)*log(p)/(p+1-a_p)) at good residues.',
        'order': 'score descending, good-prime count descending, denominator ascending, signed numerator ascending',
        'population': 'All 20,400,078 signed primitive nonzero n/d per family with abs(n),d<=4096; six families,122,400,468 addresses. Zero and infinity excluded. This contains the earlier H1024 box.',
        'mathematical_gate': 'The completed compact six-family pilot produced six new curves including a certified rank25 example, and the full-rational incidence audit excludes transfer of further generic directions from the other existing presentations. A balanced wider rational box tests candidate incidence on these same six exact compact families; generic rank alone is not the gate.',
        'scanner_wall_seconds': 180, 'selection_worker_wall_seconds': 400,
        'selection_worker_rss_bytes': 1073741824, 'maximum_selection_workers': 2,
        'point_followup_scope': 'A separate fixed24-address batch after all six fresh generic parity censuses pass exact replay. Every address receives a new bounded attempt, including repeated addresses; no replacement or catalogue prefilter.43 generic charts, height100000,4seconds per chart,300seconds/1.5GiB per worker, at most4 concurrently. Numerical heights explicitly use384bits; exact membership and quotient-only caches preserve proof gates.',
        'target_free_boundary': 'No known-record parameters, ranks, target j-invariants or public points are read by selection or prospective point workers. Only prior prospective yield supplies this scheduling gate. Post-batch comparison determines which curves are new.',
        'claim_boundary': 'This is a finite incidence experiment, not a rank classifier or upper bound. Reused trace bytes are hash-bound; the full candidate ranking uses every declared prime.'})

def parse(text, sign):
    rows, summary = [], None
    for line in text.splitlines():
        v = line.split()
        if not v:
            continue
        if v[0] == 'C':
            _, n, d, a, b, g, h = v
            n, d = sign*int(n), int(d)
            rows.append({'numerator': n, 'denominator': d, 'parameter': str(F(n, d)),
                         'score_units': int(F(a)*10**12), 'good_primes': int(g)})
        elif v[0] == 'S':
            summary = list(map(int, v[1:]))
    if summary is None:
        raise ArithmeticError('scanner did not finish')
    return rows, summary

def select(directory, family):
    protocol = cert.read(directory/'protocol.json')
    if protocol['sources'] != sources() or protocol['scanner_binary_sha256'] != cert.hashed(BINARY):
        raise ArithmeticError('frozen selection source differs')
    folder = directory/family
    if (folder/'population.json').exists():
        raise FileExistsError('preserve wider population')
    shards = []
    for sign in (-1, 1):
        source = protocol['trace_tables'][family][str(sign)]
        table = ROOT/source['path']
        if cert.hashed(table) != source['sha256']:
            raise ArithmeticError('retained trace table changed')
        result = capture([str(BINARY), str(table), '4096', '4096', '128', '0', '1'],
                         limits=Limits(180, 536870912), log_path=folder/f'scan-{sign}.log')
        rows, summary = parse(result.stdout, sign)
        shard = {'rows': rows, 'summary': summary, 'supervision': result.supervision,
                 'table_source': source, 'protocol_hash': digest(protocol)}
        checkpoint(folder/f'scan-{sign}.json', shard)
        shards.append(shard)
    rows = [r for s in shards for r in s['rows']]
    rows.sort(key=lambda r: (-r['score_units'], -r['good_primes'], r['denominator'], r['numerator']))
    if not all(rows[3]['score_units'] > s['rows'][-1]['score_units'] for s in shards):
        raise ArithmeticError('finalist boundary requires tie-complete enumeration')
    count = sum(s['summary'][3] for s in shards)
    if count != 20400078 or len({r['parameter'] for r in rows}) != len(rows):
        raise ArithmeticError('primitive population count or deduplication differs')
    checkpoint(folder/'population.json', {'family': family, 'protocol_hash': digest(protocol),
        'candidate_count': count, 'retained_candidates': rows[:128], 'finalists': rows[:4],
        'unused_H_band': 'duplicate p5; unused, not validation', 'target_free': True})
    print('R17 H4096 SELECTED', family, [(r['parameter'], r['score_units']/10**12) for r in rows[:4]], flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage', choices=['prepare', 'select'])
    p.add_argument('--directory', type=Path, default=DIRECTORY)
    p.add_argument('--family')
    a = p.parse_args()
    prepare(a.directory) if a.stage == 'prepare' else select(a.directory, a.family)
