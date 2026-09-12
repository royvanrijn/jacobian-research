"""Proposed work and inherited obligations; mathematical state is read live."""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
import re

from research import AREAS, area, cell, claim_fingerprint, stale_claims
from research_programmes import is_active, validate_programme


def load_work(root: Path) -> tuple[dict, dict]:
    return (json.loads((root / 'knowledge/work_items.json').read_text()),
            json.loads((root / 'knowledge/legacy_work_review.json').read_text()))


def original_checkboxes(text: str) -> list[dict]:
    """Retain every unchecked item, including its indented completed subsets."""
    found, current = [], None
    for number, line in enumerate(text.splitlines(), 1):
        match = re.match(r'^- \[ \] (.*)', line)
        if match:
            current = {'line': number, 'text': match[1]}
            found.append(current)
        elif current and line.startswith('  '):
            current['text'] += ' ' + line.strip()
        else:
            current = None
    return found


def validate(data: dict, legacy: dict, entries: list[dict], root: Path,
             *, check_reviews: bool = False, require_coverage: bool = True) -> None:
    assert data['schema'] == legacy['schema'] == 1
    assert data['authority'] == 'MATH_STATUS.json'
    known = {e['id']: e for e in entries}
    actions = data['actions']
    ids = {a['id'] for a in actions}
    assert len(ids) == len(actions), 'duplicate work ID'
    targets = [a['target'] for a in actions if a['target']]
    assert len(targets) == len(set(targets)), 'duplicate primary work target'
    active = {e['id'] for e in entries if e['kind'] == 'open_problem' and e['state'] == 'open'}
    if require_coverage:
        assert active <= set(targets), 'open problems missing work proposals: ' + ', '.join(sorted(active - set(targets)))
    required = {'id', 'kind', 'area', 'title', 'target', 'claims', 'sources',
                'next_step', 'done_when', 'prerequisites', 'compute', 'reviewed_claims'}
    for action in actions:
        assert required <= set(action) <= required | {'programme_status'}, f"{action.get('id')}: invalid work schema"
        validate_programme(action)
        assert action['id'].startswith('WORK-') and action['area'] in AREAS
        assert action['kind'] in {'research', 'review', 'replay', 'provenance', 'maintenance'}
        assert action['compute'] in {'none', 'bounded-replay', 'new-compute'}
        for field in ('title', 'next_step', 'done_when', 'prerequisites'):
            assert isinstance(action[field], str) and action[field].strip()
        assert isinstance(action['claims'], list) and len(action['claims']) == len(set(action['claims']))
        assert set(action['claims']) <= set(known), f"{action['id']}: unknown claim"
        assert set(action['reviewed_claims']) == set(action['claims'])
        assert all(re.fullmatch(r'[0-9a-f]{16}', value) for value in action['reviewed_claims'].values())
        if action['target']:
            assert action['target'] in action['claims']
            assert known[action['target']]['kind'] == 'open_problem'
        if check_reviews:
            stale = stale_claims(action, entries)
            assert not stale, f"{action['id']}: review changed scopes before updating fingerprints: {', '.join(stale)}"
        assert action['sources'] and len(action['sources']) == len(set(action['sources']))
        for source in action['sources']:
            path = Path(source.split('#', 1)[0])
            assert not path.is_absolute() and '..' not in path.parts
            assert (root / path).is_file(), f"{action['id']}: missing source {source}"
    path = Path(legacy['source'])
    assert not path.is_absolute() and '..' not in path.parts
    source = root / path
    assert hashlib.sha256(source.read_bytes()).hexdigest() == legacy['source_sha256'], 'historical checklist changed'
    assert original_checkboxes(source.read_text()) == [
        {'line': item['line'], 'text': item['text']} for item in legacy['items']
    ], 'every original unchecked item must be accounted for exactly once, without losing completed subsets'
    legacy_ids = [item['id'] for item in legacy['items']]
    assert len(legacy_ids) == len(set(legacy_ids))
    for item in legacy['items']:
        assert item['id'] == 'LEGACY-20260904-' + str(item['line'])
        assert item['actions'] and set(item['actions']) <= ids
        assert set(item.get('claims', [])) <= set(known), 'unknown legacy reconciliation evidence'
        assert item['rationale'].strip()
        assert item['disposition'] in {
            'maintained-check', 'superseded-navigation', 'retained-review', 'unscheduled-replay',
            'needs-reconciliation', 'retained-unknown', 'partly-superseded', 'maintained-open'}
        if 'resolution' in item:
            resolution = item['resolution']
            assert resolution['outcome'] == 'completed-at-snapshot'
            assert resolution['summary'].strip() and resolution['checks']
            assert re.fullmatch(r'\d{4}-\d{2}-\d{2}', resolution['date'])
            assert resolution['evidence'], 'a legacy completion needs inspected evidence'
            for source, digest in resolution['evidence'].items():
                path = Path(source)
                assert not path.is_absolute() and '..' not in path.parts
                assert (root / path).is_file()
                assert re.fullmatch(r'[0-9a-f]{64}', digest)


def action_record(action: dict, entries: list[dict], legacy: dict) -> dict:
    by_id = {e['id']: e for e in entries}
    stale = stale_claims(action, entries)
    target = by_id.get(action['target'])
    disposition = 'review-needed' if stale else 'proposed'
    if target and target['state'] != 'open':
        disposition = 'target-' + target['state']
    if not is_active(action):
        disposition = 'archived-programme'
    return {**action, 'disposition': disposition, 'review_needed': stale,
            'claim_states': {i: by_id[i]['state'] for i in action['claims']},
            'replacement_edges': {i: by_id[i]['replaced_by'] for i in action['claims'] if by_id[i]['replaced_by']},
            'current_scope_fingerprints': {i: claim_fingerprint(by_id[i]) for i in action['claims']},
            'inherited_items': [item['id'] for item in legacy['items'] if action['id'] in item['actions']],
            'execution': 'Unscheduled. This work item does not authorize a calculation or change mathematical status.'}


def records(data: dict, legacy: dict, entries: list[dict], kind: str = 'actions',
            selected_area: str | None = None, *, history: bool = False) -> list[dict]:
    if kind == 'actions':
        rows = [action_record(action, entries, legacy) for action in data['actions']]
    else:
        rows = [{**e, 'area': area(e['canonical_source']), 'disposition': e['state']}
                for e in entries if e['state'] == kind]
    return [row for row in rows if (history or is_active(row))
            and (not selected_area or row['area'] == selected_area)]


def show(identifier: str, data: dict, legacy: dict, entries: list[dict]) -> dict | None:
    for action in data['actions']:
        if action['id'] == identifier:
            result = action_record(action, entries, legacy)
            result['supporting_claims'] = [e for e in entries if e['id'] in action['claims']]
            result['inherited_items'] = [item for item in legacy['items'] if action['id'] in item['actions']]
            return result
    for item in legacy['items']:
        if item['id'] == identifier:
            return {**item, 'source': legacy['source'], 'source_sha256': legacy['source_sha256'],
                    'boundary': legacy['boundary'],
                    'current_claims': [e for e in entries if e['id'] in item.get('claims', [])],
                    'current_actions': [action_record(a, entries, legacy) for a in data['actions'] if a['id'] in item['actions']]}
    return None


def render(data: dict, legacy: dict, entries: list[dict], root: Path) -> dict[Path, str]:
    actions = data['actions']
    active_actions = [a for a in actions if is_active(a)]
    by_id = {e['id']: e for e in entries}
    counts = Counter(e['state'] for e in entries if is_active(e))
    header = '<!-- Generated by scripts/research.py render; do not edit. -->'
    intro = ['# Unknowns and proposed work', '', header, '',
             'Mathematical states come from [MATH_STATUS.json](../MATH_STATUS.json). '
             'This ledger makes remaining obligations and suggested next steps discoverable. '
             '**All actions are unscheduled**; a work item does not authorize a search or replay.', '',
             f"There are **{counts['open']} open problems**, **{counts['partial']} partial results** and "
             f"**{counts['parked']} parked problems** in elliptic curves and supporting K3 work. "
             f"The {len(active_actions)} active work proposals include shared repository maintenance. "
             f"The {len(actions)-len(active_actions)} proposals for other programmes are archived. "
             f"All {len(legacy['items'])} inherited checklist records remain preserved; "
             'archiving is not completion.', '',
             'Start with the exact source and [method lessons](ALGORITHMS.md). The proposals summarize '
             'a useful next gate; they do not replace the full scope. Partial-result rows guarantee '
             'coverage, while the [source-level partial review](PARTIAL_REVIEW.md) records inspected '
             'obligations and explicitly lists the unreviewed remainder. '
             'A parked replacement can retire a route without proving its original statement.', '',
             '| Area | Open | Partial | Parked | Proposals |', '|---|---:|---:|---:|---:|']
    outputs = {}
    for key, label in AREAS.items():
        group = [e for e in entries if area(e['canonical_source']) == key]
        group_actions = [a for a in actions if a['area'] == key]
        count = Counter(e['state'] for e in group if is_active(e))
        active_group = [a for a in group_actions if is_active(a)]
        if active_group or any(is_active(e) for e in group):
            intro.append(f"| [{label}](work/{key}.md) | {count['open']} | {count['partial']} | {count['parked']} | {len(active_group)} |")
        lines = [f'# Work: {label}', '', header, '',
                 '[All work](../WORK_LEDGER.md) · [All claims](../../index/' + key + '.md)', '',
                 'Proposals are unscheduled. Read the full claim scope before acting; '
                 'supporting claims retain their literal state and replacement edges.', '']
        for action in sorted(group_actions, key=lambda a: not is_active(a)):
            record = action_record(action, entries, legacy)
            lines += [f"## {action['id']}", '', f"**{action['title']}**", '',
                      f"**Disposition:** {record['disposition']}. **Kind:** {action['kind']}. "
                      f"**Compute class:** {action['compute']}.", '',
                      f"**Next step:** {action['next_step']}", '',
                      f"**Completion evidence:** {action['done_when']}", '',
                      f"**Prerequisites and boundary:** {action['prerequisites']}", '']
            if record['review_needed']:
                lines += ['**Review needed:** supporting scopes changed: ' + ', '.join(record['review_needed']) + '.', '']
            if action['claims']:
                lines += ['Authority: ' + '; '.join(f"[{i}](../../{by_id[i]['canonical_source']}) ({by_id[i]['state']})" for i in action['claims']) + '.', '']
            if record['replacement_edges']:
                lines += ['Recorded replacements: ' + '; '.join(f"`{k}` → " + ', '.join(f'`{i}`' for i in v) for k, v in record['replacement_edges'].items()) + '.', '']
            lines += ['Sources: ' + '; '.join(f'[{Path(p).name}](../../{p})' for p in action['sources']) + '.', '']
            if record['inherited_items']:
                lines += ['Inherited checklist items: ' + ', '.join(f'`{i}`' for i in record['inherited_items']) + '. Full wording: `research.py show LEGACY-ID`.', '']
        for state, label_state in [('partial', 'Partial-result register'), ('parked', 'Parked and replaced problems')]:
            subset = sorted((e for e in group if e['state'] == state), key=lambda e: (e['canonical_source'], e['id']))
            if not subset:
                continue
            lines += [f'## {label_state}', '',
                      'These are literal authority records, not newly scheduled tasks. '
                      '`research.py show ID` returns the full current scope, evidence and dependencies.', '',
                      '| ID | Result / canonical source | Recorded replacements |', '|---|---|---|']
            for e in subset:
                replacements = ', '.join(f'`{i}`' for i in e['replaced_by']) or '—'
                lines.append(f"| `{e['id']}` | [{cell(e['title'])}](../../{e['canonical_source']}) | {replacements} |")
            lines.append('')
        outputs[root / f'knowledge/work/{key}.md'] = '\n'.join(lines) + '\n'
    intro += ['', '## Suggested starting points', '',
              '| Action | Remaining gate |', '|---|---|']
    for action in actions:
        if is_active(action) and action['target']:
            intro.append(f"| [{action['id']}](work/{action['area']}.md#{action['id'].lower()}) | {cell(action['title'])} |")
    intro += ['', '## Retrieve the details', '', '```sh',
              'python3 research/scripts/research.py work --area elkies-k3',
              'python3 research/scripts/research.py work --kind partial --area elliptic-curves',
              'python3 research/scripts/research.py show WORK-EC-NEXT',
              'python3 research/scripts/research.py show LEGACY-20260904-556', '```', '',
              '[Historical checklist reconciliation](LEGACY_WORK_REVIEW.md) retains every original unchecked '
              'item and its completed subsets. [work_items.json](work_items.json) stores only proposed '
              'actions and reviewed scope fingerprints. A changed supporting scope makes an item require '
              'review; rendering cannot approve it.', '',
              'Unknowns inside proved bounded claims remain in their full scopes. Curve conductor '
              'unknowns are exposed by `research.py resources --kind curve --unknown-only`; '
              'archived support-gate records are available through `--kind support --history`. '
              'Neither resource view strengthens the authority.', '']
    outputs[root / 'knowledge/WORK_LEDGER.md'] = '\n'.join(intro)
    review_lines = ['# Historical checklist reconciliation', '', header, '',
                    f"All **{len(legacy['items'])} unchecked items** in the [preserved retrospective](../{legacy['source']}) "
                    'are retained in [legacy_work_review.json](legacy_work_review.json), including the long '
                    'completed-subset annotations. These dispositions describe cleanup decisions; they '
                    'do not certify completion of a mathematical audit or schedule a calculation.', '',
                    f"**{sum('resolution' in item for item in legacy['items'])}/{len(legacy['items'])} items have explicit completion records** "
                    'with the inspected snapshot hashes and checks performed. Completion applies to '
                    'that recorded snapshot; the maintained checks must still pass on current work. '
                    'An item without a completion record remains unfinished, regardless of its destination.', '',
                    f"{sum(not is_active(item) for item in legacy['items'])} inherited items now belong only to archived programmes. "
                    'Their completion status remains as recorded; they are outside the active EC/K3 cleanup scope.', '',
                    'Use `research.py show LEGACY-20260904-LINE` for the full text, rationale and current '
                    'destination. [Current work ledger](WORK_LEDGER.md).', '',
                    '| Original line / ID | Original item (excerpt) | Disposition | Completion | Current destination |', '|---|---|---|---|---|']
    action_map = {a['id']: a for a in actions}
    for item in legacy['items']:
        preview = item['text'][:135] + ('…' if len(item['text']) > 135 else '')
        links = ', '.join(f"[{i}](work/{action_map[i]['area']}.md#{i.lower()})" for i in item['actions'])
        resolution = item.get('resolution')
        completion = 'completed ' + resolution['date'] if resolution else 'unfinished'
        if not is_active(item):
            completion += '; archived programme'
        review_lines.append(f"| `{item['id']}` | {cell(preview)} | {item['disposition']} | {completion} | {links} |")
    outputs[root / 'knowledge/LEGACY_WORK_REVIEW.md'] = '\n'.join(review_lines) + '\n'
    return outputs
