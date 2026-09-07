#!/usr/bin/env sage-python
"""Replay the complete frozen triangle result, including its QQ inverse.

Rebuild the dictionary, ranks and rational pencils, then verify the stored
birational conversion and independently reconstruct the sixteen modular
j-maps by generic plane-cubic algebra. No new parameter or point search.
"""
from pathlib import Path
import json
import runpy
import signal

ROOT = Path(__file__).resolve().parents[2]
signal.alarm(180)
runpy.run_path(str(ROOT/'elkies-k3/scripts/search_curve302_shortword_triangles.sage'))['replay']()
result = runpy.run_path(str(ROOT/'elkies-k3/scripts/certify_curve302_triangle_generic_conversion.sage'))['build'](True)
branches = runpy.run_path(str(ROOT/'elkies-k3/scripts/certify_curve302_ten_triangle_branches.sage'))['build']()
assert branches == json.loads((ROOT/'artifacts/generated-results/elkies-k3-curve302-ten-triangle-branches-v1.json').read_text())
print('PASS_COMPLETE_TEN_TRIANGLE_REPLAY',result['QQ_inverse']['prime'],flush=True)
