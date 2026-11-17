"""Replay the optional large literal certificate in six checkpointed workers."""
from pathlib import Path
import subprocess
import sys

root=Path(__file__).resolve().parent
for column in range(1,7):
    subprocess.run([sys.executable,str(root/'check_lift_column.py'),str(column)],check=True,timeout=30)
print('PASS ALL SIX EXACT REVERSE-IDEAL CERTIFICATES')
