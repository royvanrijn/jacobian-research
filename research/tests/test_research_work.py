"""Unknowns must stay visible without reviving completed work or running CAS."""

import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import research
import research_resources
import research_work


def claim(identifier='OP-TEST', state='open', **updates):
    entry = dict(id=identifier, state=state, kind='open_problem', title='Remaining marked gate',
                 canonical_source='canonical.md', scope='The full marked model remains UNKNOWN.',
                 dependencies=[], replaced_by=[], software_lock=[], checker='never_execute.py')
    entry.update(updates)
    return entry


def fixture(root, entry):
    text = '- [ ] Reconcile parent recovery.\n  Completed subset: the alternative family.\n'
    (root / 'old.txt').write_text(text)
    (root / 'canonical.md').write_text('# Marked gate\n')
    action = dict(id='WORK-TEST', kind='research', area='core', title='Review the marked gate',
                  target=entry['id'], claims=[entry['id']], sources=['canonical.md'],
                  next_step='Identify the exact missing marking.', done_when='A complete marking is certified.',
                  prerequisites='The alternative family is already complete.', compute='none',
                  reviewed_claims={entry['id']: research.claim_fingerprint(entry)})
    legacy = dict(schema=1, source='old.txt', source_sha256=hashlib.sha256(text.encode()).hexdigest(),
                  boundary='Historical review only.', items=[dict(
                      id='LEGACY-20260904-1', line=1,
                      text='Reconcile parent recovery. Completed subset: the alternative family.',
                      actions=['WORK-TEST'], disposition='partly-superseded',
                      rationale='Keep the unresolved marking separate from the completed family.')])
    return dict(schema=1, authority='MATH_STATUS.json', purpose='Unscheduled proposals.', actions=[action]), legacy


class WorkTests(unittest.TestCase):
    def test_mapping_an_inherited_item_is_not_completing_it(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entry = claim()
            data, legacy = fixture(root, entry)
            text = research_work.render(data, legacy, [entry], root)[root / 'knowledge/LEGACY_WORK_REVIEW.md']
            self.assertIn('0/1', text)
            self.assertIn('unfinished', text)
            legacy['items'][0]['resolution'] = dict(
                outcome='completed-at-snapshot', date='2026-09-12',
                summary='Claimed completed.', checks=['metadata check'], evidence={})
            with self.assertRaisesRegex(AssertionError, 'completion needs inspected evidence'):
                research_work.validate(data, legacy, [entry], root)

    def test_new_open_problem_requires_a_proposal(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entries = [claim()]
            data, legacy = fixture(root, entries[0])
            research_work.validate(data, legacy, entries, root, check_reviews=True)
            with self.assertRaisesRegex(AssertionError, 'missing work proposals'):
                research_work.validate(data, legacy, entries + [claim('OP-NEW')], root)
            # Retrieval must remain available while a newly added frontier is
            # waiting for its proposal; generation and CI still reject the gap.
            research_work.validate(data, legacy, entries + [claim('OP-NEW')], root, require_coverage=False)

    def test_render_does_not_approve_a_changed_or_closed_target(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            original = claim()
            data, legacy = fixture(root, original)
            closed = claim(state='parked', scope='Alternative parent complete; original provenance UNKNOWN.',
                           replaced_by=['RECOVERED'])
            replacement = claim('RECOVERED', state='proved', kind='theorem')
            entries = [closed, replacement]
            views = research_work.render(data, legacy, entries, root)
            self.assertIn('target-parked', views[root / 'knowledge/work/core.md'])
            self.assertIn('Review needed', views[root / 'knowledge/work/core.md'])
            row = research_work.records(data, legacy, entries)[0]
            self.assertEqual(row['replacement_edges'], {'OP-TEST': ['RECOVERED']})
            with self.assertRaisesRegex(AssertionError, 'review changed scopes'):
                research_work.validate(data, legacy, entries, root, check_reviews=True)

    def test_every_legacy_item_and_completed_subset_survives(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entry = claim()
            data, legacy = fixture(root, entry)
            for mutation in ('drop', 'duplicate', 'strip-subset'):
                damaged = copy.deepcopy(legacy)
                if mutation == 'drop':
                    damaged['items'] = []
                elif mutation == 'duplicate':
                    damaged['items'] *= 2
                else:
                    damaged['items'][0]['text'] = 'Reconcile parent recovery.'
                with self.assertRaisesRegex(AssertionError, 'every original unchecked item'):
                    research_work.validate(data, damaged, [entry], root)
            detail = research_work.show('LEGACY-20260904-1', data, legacy, [entry])
            self.assertIn('Completed subset', detail['text'])
            self.assertEqual(detail['current_actions'][0]['disposition'], 'proposed')

    def test_missing_legacy_destination_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entry = claim()
            data, legacy = fixture(root, entry)
            legacy['items'][0]['actions'] = ['WORK-NOT-FOUND']
            with self.assertRaises(AssertionError):
                research_work.validate(data, legacy, [entry], root)

    def test_every_partial_and_parked_record_has_one_destination(self):
        entries = [claim('P1', 'partial'), claim('P2', 'partial', canonical_source='HC4_NOTE.md'),
                   claim('OLD', 'parked', canonical_source='elkies-k3/old.md', replaced_by=['P1'])]
        views = research_work.render({'actions': []}, {'items': [], 'source': 'old.txt'}, entries, Path('/virtual'))
        pages = [value for path, value in views.items() if path.parent.name == 'work']
        for entry in entries:
            self.assertEqual(sum(line.startswith(f"| `{entry['id']}` |")
                                 for page in pages for line in page.splitlines()), 1)
        self.assertEqual(research_work.records({'actions': []}, {}, entries, 'partial', 'hessian')[0]['id'], 'P2')

    def test_work_lookup_never_runs_its_checker(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entry = claim(checker=str(root / 'never_execute.py'))
            (root / 'never_execute.py').write_text("raise RuntimeError('Must only read metadata')\n")
            data = fixture(root, entry)
            detail = research.show('WORK-TEST', [entry], [], work_data=data)
            self.assertEqual(detail['supporting_claims'][0]['checker'], entry['checker'])
            self.assertIn('Unscheduled', detail['execution'])


class StructuredRetrievalTests(unittest.TestCase):
    def make_curve_database(self, root):
        curve = dict(id='test-curve', rank_lower_bound=28, local_search_rank_lower_bound=27,
                     rank_provenance='PUBLIC_POINT_REPRODUCTION', conductor_status='UNKNOWN',
                     conductor=None, family='11952', parameter='110314/102227', icarm_ids=[619],
                     points=[['1', '2']], rank_certificate={'rank_lower_bound': 28, 'signatures': ['large payload']})
        path = root / research_resources.SOURCES['curve']
        path.parent.mkdir(parents=True)
        path.write_text(json.dumps({'curves': [curve]}))
        return path

    def test_curve_bounds_and_unknowns_retain_their_provenance(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = self.make_curve_database(root)
            rows = research_resources.records(root, [])
            row = rows[0]
            self.assertIn('rank >= 28', row['title'])
            self.assertEqual(row['record']['local_search_rank_lower_bound'], 27)
            self.assertNotIn('large payload', row['summary'])
            self.assertEqual(len(research_resources.filtered(rows, 'curve', True)), 1)
            detail = research.show(row['id'], [], [], resource_rows=rows)
            self.assertEqual(detail['record']['points'], [['1', '2']])
            # A live update, not a generated navigation cache, controls filtering.
            database = json.loads(path.read_text())
            database['curves'][0]['conductor_status'] = 'EXACT'
            path.write_text(json.dumps(database))
            self.assertEqual(research_resources.filtered(research_resources.records(root, []), 'curve', True), [])

    def test_process_mechanism_is_searchable_but_not_a_current_theorem(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / research_resources.SOURCES['k3']
            path.parent.mkdir(parents=True)
            path.write_text(json.dumps(dict(mechanisms=[dict(id='MECH-05', name='Full denominator',
                                                            resolution='Cancel the entire pole divisor.', evidence=['OP-TEST'])],
                                            stages=[], transitions=[], events=[], literature=[])))
            entries = [claim()]
            rows = research_resources.records(root, entries)
            with patch.object(research, 'ROOT', root), patch.object(research, 'source_documents', return_value=[]):
                hits = research.search('full denominator', entries, [], resource_rows=rows)
            self.assertEqual(hits[0]['id'], 'K3-MECHANISM-MECH-05')
            self.assertEqual(hits[0]['kind'], 'resource')
            self.assertNotIn('state', hits[0])
            self.assertEqual(hits[0]['related_claims'], ['OP-TEST'])

    def test_canonical_tex_body_is_retrievable(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'proof.tex').write_text('\\section{Conductor gluing}\nA tagged obstruction survives at infinity.')
            entries = [claim(canonical_source='proof.tex')]
            with patch.object(research, 'ROOT', root), patch.object(research, 'source_documents', return_value=[]):
                hits = research.search('tagged obstruction', entries, [])
            self.assertEqual(hits[0]['source'], 'proof.tex')
            self.assertEqual(hits[0]['kind'], 'source')


if __name__ == '__main__':
    unittest.main()
