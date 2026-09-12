"""A failed CAS record must not retain parseable apparent successes."""

from pathlib import Path
import re
import subprocess
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from singular_result_runner import run_singular_result


class SingularResultTests(unittest.TestCase):
    source = 'print("RESULT_BEGIN"); print("unit=1"); print("RESULT_END"); quit;\n'

    def execute(self, output, *, finish=True):
        def process(command, **options):
            marker = re.search(r'print\("(SINGULAR_RUN_COMPLETE_[a-f0-9]+)"\)', options['input'])[1]
            return subprocess.CompletedProcess(command, 0, output + (marker + '\n' if finish else ''), '')
        with patch('singular_runner.subprocess.run', side_effect=process):
            return run_singular_result(self.source)

    def test_clean_result_preserves_the_original_output(self):
        output = 'RESULT_BEGIN\nunit=1\nRESULT_END\n'
        result = self.execute(output)
        self.assertEqual((result.returncode, result.stdout, result.stderr), (0, output, ''))

    def test_error_before_apparent_unit_returns_no_parseable_success(self):
        result = self.execute(' ? failed coefficient\nRESULT_BEGIN\nunit=1\nRESULT_END\n')
        self.assertEqual(result.returncode, 1)
        self.assertEqual(result.stdout, '')
        self.assertIn('failed coefficient', result.stderr)

    def test_incomplete_or_repeated_result_fails(self):
        for output, finish in [('RESULT_BEGIN\nunit=1\nRESULT_END\n', False), ('RESULT_BEGIN\nRESULT_END\nRESULT_END\n', True)]:
            with self.subTest(output=output, finish=finish):
                self.assertEqual(self.execute(output, finish=finish).returncode, 1)

    def test_timeout_propagates_and_unexpected_program_shape_is_rejected(self):
        with patch('singular_runner.subprocess.run', side_effect=subprocess.TimeoutExpired('Singular', 1)):
            with self.assertRaises(subprocess.TimeoutExpired):
                run_singular_result(self.source, timeout=1)
        with self.assertRaisesRegex(ValueError, 'explicit quit'):
            run_singular_result('print("RESULT_END");')


if __name__ == '__main__':
    unittest.main()
