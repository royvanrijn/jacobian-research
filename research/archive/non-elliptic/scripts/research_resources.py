"""Read bounded projections of maintained JSON ledgers, never generated workers."""

from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

SOURCES = {
    'k3': 'elkies-k3/data/process/elkies_k3_process_ledger.json',
    'support': 'verified/SUPPORT_SATURATION_PATHS.json',
    'curve': 'elliptic-curves/data/research_curves/database.json',
}
KINDS = ('mechanism', 'stage', 'transition', 'event', 'literature', 'support', 'curve')
BOUNDARY = ('Structured reference record, not a new mathematical status. Read related claims in '
            'MATH_STATUS.json; a historical process stage does not supersede its current scope.')


def records(root: Path, entries: list[dict]) -> list[dict]:
    """Only three maintained inputs; no recursive scan of computation output."""
    claims_by_path = defaultdict(set)
    known_ids = {e['id'] for e in entries}
    for entry in entries:
        for path in [entry['canonical_source'], *entry.get('software_lock', [])]:
            claims_by_path[path].add(entry['id'])
    result = []

    def add(identifier, kind, heading, source, payload, group, *, summary=None):
        body = json.dumps(payload if summary is None else summary, ensure_ascii=False, sort_keys=True)
        related = set()
        for reference in payload.get('evidence', []):
            if reference in known_ids:
                related.add(reference)
            related.update(claims_by_path[reference])
        for key in ('theorem_status_id', 'formal_core_status_id', 'status_id'):
            if payload.get(key) in known_ids:
                related.add(payload[key])
        # Model-frontier references are precise claim IDs in the support ledger.
        for frontier in payload.get('model_frontiers', []):
            if frontier.get('status_id') in known_ids:
                related.add(frontier['status_id'])
        if payload.get('canonical_source'):
            related.update(claims_by_path[payload['canonical_source']])
        for key in ('rank_source_certificate', 'conductor_certificate'):
            if payload.get(key):
                related.update(claims_by_path[payload[key]])
        result.append(dict(id=identifier, kind=kind, title=heading, source=source,
                           area=group, record=payload, summary=body,
                           related_claims=sorted(related), boundary=BOUNDARY,
                           temporal_role='historical' if kind in {'event', 'stage', 'transition'} else 'reference'))

    process_path = root / SOURCES['k3']
    if process_path.is_file():
        process = json.loads(process_path.read_text())
        for plural, kind in [('mechanisms', 'mechanism'), ('stages', 'stage'),
                             ('transitions', 'transition'), ('events', 'event'), ('literature', 'literature')]:
            for record in process[plural]:
                heading = next((record[key] for key in ('name', 'label', 'behavior', 'title', 'operation')
                                if record.get(key)), record['id'])
                # Stage IDs overlap event/transition syntax; include the kind.
                identifier = 'K3-' + kind.upper() + '-' + record['id']
                add(identifier, kind, heading, SOURCES['k3'], record, 'elkies-k3')
    support_path = root / SOURCES['support']
    if support_path.is_file():
        support = json.loads(support_path.read_text())
        for identifier, record in support['programmes'].items():
            payload = {**record, 'common_gate': support['common_gate'],
                       'theorem_status_id': support['theorem_status_id'],
                       'formal_core_status_id': support['formal_core_status_id']}
            add('SUPPORT-' + identifier, 'support', identifier + ': ' + record['target'],
                SOURCES['support'], payload, {'plane-jc': 'plane-jc', 'cubic-keller': 'core'}.get(identifier, 'geometry'))
    curve_path = root / SOURCES['curve']
    if curve_path.is_file():
        inventory = json.loads(curve_path.read_text())
        for record in inventory['curves']:
            # The export's certified bound can be stronger than the earlier
            # local-search bound. Keep both with their provenance; never use max
            # or a public rank field to invent an exact rank.
            summary = {k: v for k, v in record.items()
                       if k not in {'points', 'rank_certificate', 'bad_primes', 'known_bad_primes', 'discriminant'}}
            summary['point_count'] = len(record['points'])
            heading = (f"{record['id']}: recorded rank >= {record['rank_lower_bound']}; "
                       f"conductor {record['conductor_status']}; family {record.get('family', 'unspecified')}; "
                       f"parameter {record.get('parameter', 'unspecified')}; ICARM {record.get('icarm_ids', [])}")
            add('CURVE-' + record['id'], 'curve', heading, SOURCES['curve'], record,
                'elliptic-curves', summary=summary)
    ids = [record['id'] for record in result]
    assert len(ids) == len(set(ids)), 'duplicate structured resource ID'
    return result


def filtered(rows: list[dict], kind: str | None = None, unknown_only: bool = False) -> list[dict]:
    if kind:
        rows = [row for row in rows if row['kind'] == kind]
    if unknown_only:
        # Explicit machine state only. Word matching in prose would confuse an
        # excluded "unknown" possibility with an actual missing certificate.
        rows = [row for row in rows if row['kind'] == 'curve' and row['record']['conductor_status'] == 'UNKNOWN']
    return rows


def render(rows: list[dict], root: Path) -> dict[Path, str]:
    lines = ['# Structured research records', '',
             '<!-- Generated by scripts/research.py render; do not edit. -->', '',
             'These maintained JSON ledgers contain process lessons, support gates and curve '
             'evidence that a Markdown-only search misses. Retrieval reads them directly. '
             'Their records are references; [MATH_STATUS.json](../MATH_STATUS.json) retains '
             'mathematical authority.', '',
             '| Source | Records | SHA-256 |', '|---|---:|---|']
    counts = Counter(row['source'] for row in rows)
    for path in SOURCES.values():
        source = root / path
        assert source.is_file(), f'missing maintained structured source: {path}'
        lines.append(f'| [{Path(path).name}](../{path}) | {counts[path]} | `{hashlib.sha256(source.read_bytes()).hexdigest()}` |')
    kinds = Counter(row['kind'] for row in rows)
    unknown = len(filtered(rows, 'curve', True))
    lines += ['', 'Available records: ' + ', '.join(f'{value} {key}' for key, value in sorted(kinds.items())) + '.', '',
              f'The curve export records **{unknown} UNKNOWN conductors**. This is a certificate '
              'availability frontier, not a command to factor them. Its rank column is a lower '
              'bound with its own provenance; earlier local-search bounds are kept separately.', '',
              '```sh', 'python3 research/scripts/research.py search "full denominator" --kind resource',
              'python3 research/scripts/research.py resources --kind mechanism',
              'python3 research/scripts/research.py show K3-MECHANISM-MECH-05',
              'python3 research/scripts/research.py show SUPPORT-cubic-keller',
              'python3 research/scripts/research.py show CURVE-new-20260906-188',
              'python3 research/scripts/research.py resources --kind curve --unknown-only', '```', '',
              '`show` returns the complete selected record, including evidence. Search excludes '
              'large point-coordinate and finite-signature arrays from snippets; they remain '
              'available in that exact record. Process events retain their historical timestamps '
              'and reported stages without promoting them to current status.', '',
              '[Work ledger](../knowledge/WORK_LEDGER.md) · [All research](README.md)', '']
    return {root / 'index/resources.md': '\n'.join(lines)}
