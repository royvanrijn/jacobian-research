"""A rendered catalogue must never certify an unfinished or stale source review."""

import hashlib
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import research_partial_reviews as reviews


def fixture(root):
    (root / 'proof.md').write_text('The tested chart is closed; the full problem remains UNKNOWN.\n')
    (root / 'checker.py').write_text("raise RuntimeError('retrieval must not execute research')\n")
    entry = dict(id='P', state='partial', title='A bounded chart', canonical_source='proof.md',
                 scope='The full problem remains UNKNOWN.', dependencies=[], checker='checker.py',
                 software_lock=[], independent_replay=False)
    review = dict(id='P', source='proof.md', completed='One chart is closed.',
                  remaining='Classify the other charts.', disposition='bounded-only',
                  source_review='Read the proof.', input_review='Inspected the source inputs.',
                  dependency_review='No external theorem is assumed.',
                  evidence=[dict(path='proof.md', quote='the full problem remains UNKNOWN')],
                  reviewed_claims={'P': reviews.entry_fingerprint(entry)},
                  reviewed_files={p: hashlib.sha256((root / p).read_bytes()).hexdigest()
                                  for p in ('proof.md', 'checker.py')})
    return entry, dict(schema=1, reviews=[review])


class PartialReviewTests(unittest.TestCase):
    def test_complete_gate_rejects_an_unreviewed_partial(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entry, data = fixture(root)
            entries = [entry, {**entry, 'id': 'NEW'}]
            reviews.validate(data, entries, root, check_reviews=True)
            text = reviews.render(data, entries, root)[root / 'knowledge/PARTIAL_REVIEW.md']
            self.assertIn('1/2', text)
            self.assertIn('`NEW`', text)
            with self.assertRaisesRegex(AssertionError, 'not yet reviewed: NEW'):
                reviews.validate(data, entries, root, require_complete=True)

    def test_dependency_input_and_assurance_changes_invalidate_the_review(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entry, data = fixture(root)
            for field, value in [('dependencies', ['external: new lemma']),
                                 ('software_lock', ['new-input.json']), ('independent_replay', True)]:
                with self.subTest(field=field):
                    changed = {**entry, field: value}
                    text = reviews.render(data, [changed], root)[root / 'knowledge/PARTIAL_REVIEW.md']
                    self.assertIn('0/1', text)
                    self.assertIn('REVIEW NEEDED', text)
                    with self.assertRaisesRegex(AssertionError, 'needs reconciliation'):
                        reviews.validate(data, [changed], root, require_complete=True)

    def test_quote_match_cannot_hide_a_changed_source(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entry, data = fixture(root)
            with (root / 'proof.md').open('a') as stream:
                stream.write('A later counterexample invalidates the earlier exclusion.\n')
            with self.assertRaisesRegex(AssertionError, 'proof.md'):
                reviews.validate(data, [entry], root, check_reviews=True)

    def test_missing_passage_is_rejected_and_checker_is_never_imported(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entry, data = fixture(root)
            reviews.validate(data, [entry], root, require_complete=True)
            self.assertEqual(reviews.detail('P', data, [entry], root)['review_needed'], [])
            data['reviews'][0]['evidence'][0]['quote'] = 'The full problem is solved.'
            with self.assertRaisesRegex(AssertionError, 'missing review evidence passage'):
                reviews.validate(data, [entry], root)

    def test_recovered_input_requires_identity_review(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entry, data = fixture(root)
            data['reviews'][0]['missing_inputs'] = ['old-marking.json']
            reviews.validate(data, [entry], root, require_complete=True)
            (root / 'old-marking.json').write_text('{"schema": "different-marking"}')
            self.assertIn('old-marking.json', reviews.detail('P', data, [entry], root)['review_needed'])
            with self.assertRaisesRegex(AssertionError, 'old-marking.json'):
                reviews.validate(data, [entry], root, require_complete=True)

    def test_missing_input_cannot_also_be_hash_reviewed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entry, data = fixture(root)
            data['reviews'][0]['missing_inputs'] = ['checker.py']
            with self.assertRaisesRegex(AssertionError, 'both present and missing'):
                reviews.validate(data, [entry], root)

    def test_local_receipt_reports_absence_and_still_rejects_changed_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entry, data = fixture(root)
            name = 'artifacts/local/exact-input.json'
            original = b'{"bounded_result": "UNKNOWN"}'
            data['reviews'][0]['reviewed_local_files'] = {name: hashlib.sha256(original).hexdigest()}
            reviews.validate(data, [entry], root, require_complete=True)
            self.assertEqual(reviews.detail('P', data, [entry], root)['unavailable_local_evidence'], [name])
            path = root / name
            path.parent.mkdir(parents=True)
            path.write_bytes(original)
            reviews.validate(data, [entry], root, require_complete=True)
            path.write_text('{"unjustified_result": "PROVED"}')
            with self.assertRaisesRegex(AssertionError, 'needs reconciliation'):
                reviews.validate(data, [entry], root, require_complete=True)

    def test_canonical_sources_cannot_be_downgraded_to_optional_receipts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entry, data = fixture(root)
            review = data['reviews'][0]
            review['reviewed_local_files'] = {'checker.py': review['reviewed_files'].pop('checker.py')}
            with self.assertRaisesRegex(AssertionError, 'local receipts are only'):
                reviews.validate(data, [entry], root)


if __name__ == '__main__':
    unittest.main()
