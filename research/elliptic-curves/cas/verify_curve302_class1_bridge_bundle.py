#!/usr/bin/env python3
"""Read-only portable replay of the class1/core/carrier proof bundle."""
import subprocess
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
CHECKS=[
    [sys.executable,'research/elliptic-curves/cas/package_curve302_class1_bridge.py','--check'],
    ['timeout','180','sage','-python','research/elliptic-curves/cas/verify_curve302_class1_bridge.sage','--check'],
    ['timeout','120','sage','-python','research/elliptic-curves/cas/curve302_class1_prescribed_core_glue.sage','--check'],
    ['timeout','120','sage','-python','research/elliptic-curves/cas/curve302_class1_twist_sections.sage','--check'],
    ['timeout','60','sage','-python','research/elliptic-curves/cas/curve302_class1_branch_veronese.sage','--check'],
]
if __name__=='__main__':
    for command in CHECKS:
        subprocess.run(command,cwd=ROOT,check=True)
    print('PASS complete class1 bridge/core/carrier bundle',flush=True)
