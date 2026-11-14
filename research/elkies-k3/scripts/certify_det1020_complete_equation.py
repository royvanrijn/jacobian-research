"""Replay the final section and source-inverse certificates.

The prior rational-pencil and cubic checkers are named proof dependencies.
No discovery lifting, enumeration, specialization sweep or rank search runs.
"""
from pathlib import Path
import runpy
import sys

P = Path(__file__).resolve().parents[2] / 'artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1'
sys.path.insert(0, str(P))
for name in [
    'compile_trace_components.py',
    'build_trace_basis_transport.py',
    'verify_marked_sections.py',
    'verify_marked_section_reduction.py',
    'prepare_inverse_divisors.py',
    'check_source_inverse.py',
]:
    sys.argv = [str(P / name), '--check']
    runpy.run_path(str(P / name), run_name='__main__')
print('PASS_DET1020_SEVENTEEN_MARKED_SECTIONS_AND_DETERMINANTAL_SOURCE_INVERSE')
