"""Replay the exact QQ pencil, plane presentation and convenient origin."""
from pathlib import Path
import runpy,sys
P=Path(__file__).resolve().parents[2]/'artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1'
sys.path.insert(0,str(P))
for name,args in [('verify_rr_rational_numerators.py',['1024','--check']),('build_rr_plane_and_origin.py',['--check'])]:
 sys.argv=[str(P/name)]+args;runpy.run_path(str(P/name),run_name='__main__')
print('PASS_DET1020_EXACT_QQ_PENCIL_PLANE_AND_ORIGIN')
