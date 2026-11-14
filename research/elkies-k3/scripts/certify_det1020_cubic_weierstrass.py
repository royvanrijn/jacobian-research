"""Replay the fixed1020 cubic, Weierstrass maps and four known sections."""
from pathlib import Path
import runpy
import sys

P = Path(__file__).resolve().parents[2] / 'artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1'
sys.path.insert(0, str(P))
for name, args in [
    ('verify_trisection_numerators.py', ['1024', '--check']),
    ('build_trisection_quartic.py', ['--check']),
    ('extract_trisection_invariants.py', ['--check']),
    ('convert_pointed_quartic_exact.py', ['--check']),
    ('compile_component_sections.py', ['--check']),
]:
    sys.argv = [str(P / name)] + args
    runpy.run_path(str(P / name), run_name='__main__')
print('PASS_DET1020_EXACT_CUBIC_WEIERSTRASS_AND_FOUR_SECTIONS')
