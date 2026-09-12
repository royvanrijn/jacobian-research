"""Source reviews of partial results and entries without executable checkers."""

from __future__ import annotations

from collections import defaultdict
import hashlib
import json
from pathlib import Path

from research import cell
from research_programmes import is_active


def entry_fingerprint(entry: dict) -> str:
    """Bind dependency, replay and assurance metadata as well as the prose scope."""
    encoded = json.dumps(entry, sort_keys=True, separators=(',', ':'), ensure_ascii=False)
    return hashlib.sha256(encoded.encode()).hexdigest()


def load(root: Path) -> dict:
    return json.loads((root / 'knowledge/partial_reviews.json').read_text())


def stale(review: dict, entries: list[dict], root: Path) -> list[str]:
    known = {e['id']: e for e in entries}
    changed = [identifier for identifier, digest in review['reviewed_claims'].items()
               if identifier not in known or entry_fingerprint(known[identifier]) != digest]
    changed += [path for path, digest in review['reviewed_files'].items()
                if not (root / path).is_file() or hashlib.sha256((root / path).read_bytes()).hexdigest() != digest]
    # A receipt records bytes inspected in the original checkout. Absence of an
    # ignored cache is reported separately; present bytes must still match.
    changed += [path for path, digest in review.get('reviewed_local_files', {}).items()
                if (root / path).exists() and (not (root / path).is_file()
                or hashlib.sha256((root / path).read_bytes()).hexdigest() != digest)]
    # A newly recovered file still needs identity/schema review before the old
    # missing-input assessment can be counted as current.
    changed += [path for path in review.get('missing_inputs', []) if (root / path).exists()]
    return changed


def unavailable_local_evidence(review: dict, root: Path) -> list[str]:
    return [path for path in review.get('reviewed_local_files', {}) if not (root / path).is_file()]


def validate(data: dict, entries: list[dict], root: Path, *, check_reviews: bool = False,
             require_complete: bool = False, require_checker_review: bool = False) -> None:
    assert data['schema'] == 1
    known = {e['id']: e for e in entries}
    absence = data.get('checker_absence_reviews', [])
    absence_ids = {review['id'] for review in absence}
    ids = [review['id'] for review in data['reviews']]
    all_ids = ids + [review['id'] for review in absence]
    assert len(all_ids) == len(set(all_ids)), 'duplicate source review'
    for review in [*data['reviews'], *absence]:
        assert review['id'] in known
        assert review['id'] in review['reviewed_claims']
        assert set(review['reviewed_claims']) <= set(known)
        assert review['source'] in review['reviewed_files']
        assert review['source'] == known[review['id']]['canonical_source'], 'partial review source changed'
        for field in ('completed', 'remaining', 'source_review', 'input_review', 'dependency_review'):
            assert isinstance(review[field], str) and review[field].strip(), f"{review['id']}: missing {field}"
        if review['id'] in absence_ids:
            entry = known[review['id']]
            assert entry['checker'] is None, f"{entry['id']}: checker added; reconcile absence review"
            if review['disposition'] == 'written-proof':
                assert entry['state'] == 'proved' and entry['kind'] != 'open_problem'
            elif review['disposition'] == 'research-objective':
                assert entry['kind'] == 'open_problem' and entry['state'] in {'open', 'parked'}
            else:
                assert review['disposition'] == 'missing-witness', 'invalid checker-absence disposition'
        else:
            assert review['disposition'] in {'open-obligation', 'bounded-only', 'superseded-test-family', 'historical-route'}
        assert review['evidence'], 'a semantic review needs source passages'
        for evidence in review['evidence']:
            assert evidence['path'] in review['reviewed_files']
            assert evidence['quote'] in (root / evidence['path']).read_text(), f"{review['id']}: missing review evidence passage"
        for path in review['reviewed_files']:
            p = Path(path)
            assert not p.is_absolute() and '..' not in p.parts
            assert (root / p).is_file()
        for path, digest in review.get('reviewed_local_files', {}).items():
            p = Path(path)
            assert not p.is_absolute() and '..' not in p.parts
            assert p.parts[:2] in {('artifacts', 'local'), ('artifacts', 'generated-results')}, (
                'local receipts are only for computation artifacts, never canonical source code or prose'
            )
            assert path not in review['reviewed_files']
            assert len(digest) == 64 and all(c in '0123456789abcdef' for c in digest)
        missing = review.get('missing_inputs', [])
        assert isinstance(missing, list) and len(missing) == len(set(missing))
        for path in missing:
            p = Path(path)
            assert not p.is_absolute() and '..' not in p.parts
            assert path not in review['reviewed_files'] and path not in review.get('reviewed_local_files', {}), (
                'input cannot be both present and missing at review time'
            )
        if check_reviews or require_complete or require_checker_review:
            changed = stale(review, entries, root)
            assert not changed, f"{review['id']}: source-level review needs reconciliation: {', '.join(changed)}"
    if require_complete:
        remaining = {e['id'] for e in entries if is_active(e) and e['state'] == 'partial'} - set(ids)
        assert not remaining, 'partial results not yet reviewed: ' + ', '.join(sorted(remaining))
    if require_checker_review:
        remaining = {e['id'] for e in entries if is_active(e) and e['checker'] is None} - absence_ids
        assert not remaining, 'checker absence not yet reviewed: ' + ', '.join(sorted(remaining))


def detail(identifier: str, data: dict, entries: list[dict], root: Path, *,
           collection: str = 'reviews') -> dict | None:
    review = next((r for r in data.get(collection, []) if r['id'] == identifier), None)
    return None if review is None else {**review, 'review_needed': stale(review, entries, root),
                                       'unavailable_local_evidence': unavailable_local_evidence(review, root)}


def render(data: dict, entries: list[dict], root: Path) -> dict[Path, str]:
    partials = [e for e in entries if is_active(e) and e['state'] == 'partial']
    by_id = {e['id']: e for e in entries}
    reviews = {r['id']: r for r in data['reviews']}
    pending = [e for e in partials if e['id'] not in reviews]
    outdated = [e for e in partials if e['id'] in reviews and stale(reviews[e['id']], entries, root)]
    reviewed = len(partials) - len(pending) - len(outdated)
    lines = ['# Source-level review of partial results', '',
             '<!-- Generated by scripts/research.py render; do not edit. -->', '',
             f'**{reviewed}/{len(partials)} current partial results reviewed; {len(pending)} remain unreviewed; {len(outdated)} reviews need reconciliation.** '
             'This covers elliptic curves and supporting K3 work. Archived programme reviews remain in the JSON records. '
             'This is progress on the source review, not a completion certificate. '
             'Mathematical states remain in [MATH_STATUS.json](../MATH_STATUS.json).', '',
             'Each record states the exact surviving obligation, the later theorem or experiment '
             'already available, the inspected implementation/input boundary and source passages. '
             'Maintained source files and claim scopes are hash-bound. `research.py show ID` includes the full '
             'review; rendering never approves changed evidence.', '',
             '`reviewed_local_files` retains hashes of ignored calculation artifacts inspected in the '
             'review checkout. When present, their hashes must still match. When absent, `show ID` '
             'reports unavailable local evidence. A source-review PASS does not certify that every '
             'replay input is distributed or independently replayed. Missing caches never trigger reconstruction.', '',
             'Where present, `missing_inputs` names inspected replay prerequisites, not expected '
             'outputs or an automatic rebuild queue. Recovering one requires a fresh identity/schema '
             'review. The list is scoped to the replay discussed in each record.', '',
             '[Work ledger](WORK_LEDGER.md) · [Review records](partial_reviews.json)', '']
    grouped = defaultdict(list)
    for review in data['reviews']:
        if not is_active(by_id[review['id']]):
            continue
        grouped[review['source']].append(review)
    for source, group in sorted(grouped.items()):
        lines += [f'## {Path(source).stem}', '', f'[Canonical source](../{source}).', '',
                  '| Claim | Disposition | Precise remaining obligation |', '|---|---|---|']
        for review in group:
            changed = stale(review, entries, root)
            state = by_id[review['id']]['state']
            disposition = 'REVIEW NEEDED' if changed else review['disposition']
            lines.append(f"| `{review['id']}` ({state}) | {disposition} | {cell(review['remaining'])} |")
        lines.append('')
    if pending:
        lines += ['## Unreviewed source groups', '',
                  'These rows expose the remaining work. They have not received a source-level '
                  'review merely because their claim scopes were indexed.', '',
                  '| Source | Pending claim IDs |', '|---|---|']
        groups = defaultdict(list)
        for entry in pending:
            groups[entry['canonical_source']].append(entry['id'])
        for source, identifiers in sorted(groups.items()):
            lines.append(f"| [{Path(source).name}](../{source}) | " + ', '.join(f'`{i}`' for i in identifiers) + ' |')
    output = {root / 'knowledge/PARTIAL_REVIEW.md': '\n'.join(lines) + '\n'}
    if 'checker_absence_reviews' in data:
        absence = {r['id']: r for r in data['checker_absence_reviews']}
        targets = [e for e in entries if is_active(e) and e['checker'] is None]
        current = sum(e['id'] in absence and not stale(absence[e['id']], entries, root) for e in targets)
        lines = ['# Entries without executable checkers', '',
                 '<!-- Generated by scripts/research.py render; do not edit. -->', '',
                 f'**{current}/{len(targets)} current EC/K3 entries reviewed.** '
                 'This classifies the absence of a checker; it grants no independent replay, '
                 'formal verification or external review. Mathematical state remains in '
                 '[MATH_STATUS.json](../MATH_STATUS.json).', '',
                 'Written deductions can use computational premises certified by other claims. '
                 'Open objectives have no completed theorem to replay. A missing finite witness, '
                 'if identified, remains an obligation. Read `research.py show ID` for the exact '
                 'proof sections, input boundary and hash-bound source review. These records share '
                 '[partial_reviews.json](partial_reviews.json); they are excluded from partial-result counts.', '',
                 '[Partial-result review](PARTIAL_REVIEW.md) · [Work ledger](WORK_LEDGER.md)', '',
                 '| Entry | Review outcome | Checker purpose / absence |', '|---|---|---|']
        for entry in sorted(targets, key=lambda e: e['id']):
            review = absence.get(entry['id'])
            outcome = 'UNREVIEWED' if review is None else 'REVIEW NEEDED' if stale(review, entries, root) else review['disposition']
            reason = 'Inspect the canonical source.' if review is None else review['input_review']
            lines.append(f"| [{entry['id']}](../{entry['canonical_source']}) | {outcome} | {cell(reason)} |")
        output[root / 'knowledge/CHECKER_REVIEW.md'] = '\n'.join(lines) + '\n'
    return output
