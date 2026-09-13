"""Small adversarial controls; default maintenance never reads a full atlas."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    'smooth_character_witness', ROOT/'elkies-k3/scripts/verify_r17_smooth_character_witness.py')
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)


def fixture():
    w = [0]*17
    return dict(rows=[[1, 'first', w[:], ['-2', '0', '1'], ['1']],
                      [2, 'second', w[:], ['1', '0', '1'], ['1']]],
                priority=[[1, w[:]], [2, w[:]]])


class CharacterWitnessTests(unittest.TestCase):
    def test_distinct_atoms_are_independent(self):
        self.assertEqual(replay.verify_atlas(fixture(), 2), dict(
            count=2, distinct_irreducible_quadratic_atoms=2,
            priority_attachments=2, internal_three_character_relations=0))

    def test_scalar_and_denominator_do_not_change_polynomial_atom(self):
        self.assertEqual(replay.quadratic_atom(['2/3', '0', '-1/3'], ['-7/11']), (-2, 0, 1))

    def test_reducible_repeated_and_linear_polynomials_fail(self):
        for q in (['-1', '0', '1'], ['1', '-2', '1'], ['1', '2', '0']):
            with self.subTest(q=q), self.assertRaises(ValueError):
                replay.quadratic_atom(q, ['1'])

    def test_nonconstant_and_zero_denominators_fail(self):
        for d in (['0'], ['1', '1'], []):
            with self.subTest(d=d), self.assertRaises(ValueError):
                replay.quadratic_atom(['-2', '0', '1'], d)

    def test_floats_are_not_exact_input(self):
        with self.assertRaisesRegex(ValueError, 'exact rational strings'):
            replay.quadratic_atom([-2.0, '0', '1'], ['1'])

    def test_equal_support_requires_full_scalar_comparison(self):
        data = fixture()
        # 2*(u^2-2) is a DISTINCT character, but unsupported by this shortcut.
        data['rows'][1][3] = ['-4', '0', '2']
        with self.assertRaisesRegex(ValueError, 'full rational squareclasses'):
            replay.verify_atlas(data, 2)

    def test_wrong_frame_word_fails(self):
        data = fixture()
        data['rows'][1][2][0] = 1
        with self.assertRaisesRegex(ValueError, 'attachment mismatch'):
            replay.verify_atlas(data, 2)

    def test_unknown_mask_fails(self):
        data = fixture()
        data['rows'][1][0] = 3
        with self.assertRaisesRegex(ValueError, 'attachment mismatch'):
            replay.verify_atlas(data, 2)

    def test_duplicate_label_and_priority_fail(self):
        data = fixture()
        data['rows'][1][1] = 'first'
        with self.assertRaisesRegex(ValueError, 'duplicate branch identity'):
            replay.verify_atlas(data, 2)
        data = fixture()
        data['priority'][1][0] = 1
        with self.assertRaisesRegex(ValueError, 'duplicate priority mask'):
            replay.verify_atlas(data, 2)

    def test_missing_row_fails(self):
        data = fixture()
        data['rows'].pop()
        with self.assertRaisesRegex(ValueError, 'incomplete'):
            replay.verify_atlas(data, 2)

    def test_wrong_source_identity_fails(self):
        packet = dict(schema=replay.SCHEMA, atlases=[
            dict(fixture(), source=copy.deepcopy(s)) for s in replay.SPECS])
        packet['atlases'][0]['source']['vector_field'] = 'direct_hidden_w'
        with self.assertRaisesRegex(ValueError, 'source identity'):
            replay.verify(packet)

    def test_source_projection_uses_exact_fields_and_rejects_changed_input(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = fixture()
            original = dict(artifact_schema='fixture.v1', base_parameter='u', bisections=[
                dict(lattice_orbit_mask=m, label=l, fixed_w=w,
                     branch=dict(numerator_coefficients=n, denominator_coefficients=d),
                     irrelevant_status='PASS') for m, l, w, n, d in data['rows']])
            (root/'atlas.json').write_text(json.dumps(original))
            (root/'priority.tsv').write_text('orbit_mask\tpriority_w\n'+''.join(
                str(m)+'\t'+' '.join(map(str, w))+'\n' for m, w in data['priority']))
            spec = dict(atlas='atlas.json', priority='priority.tsv', atlas_schema='fixture.v1',
                        vector_field='fixed_w', priority_vector_field='priority_w')
            for kind in ('atlas', 'priority'):
                spec[kind+'_sha256'] = hashlib.sha256((root/spec[kind]).read_bytes()).hexdigest()
            self.assertEqual(replay.project_source(spec, root), dict(data, source=spec))
            original['bisections'][0]['branch']['denominator_coefficients'] = ['1', '1']
            (root/'atlas.json').write_text(json.dumps(original))
            with self.assertRaisesRegex(ValueError, 'changed original input'):
                replay.project_source(spec, root)
            (root/'atlas.json').unlink()
            with self.assertRaisesRegex(ValueError, 'missing original input'):
                replay.project_source(spec, root)


if __name__ == '__main__':
    unittest.main()
