"""Run byte-preserved ancillary source using the installed full Sage distribution.

The source expects the smaller passagemath-singular module. Only that import
is aliased; no coefficient, equation, branch or algebra algorithm is changed.
"""
import json
from pathlib import Path
import runpy
import sys

sys.path.insert(0, '/tmp/jacobian-hc4-sage-deps')
import sage.all
import flint
import sympy
import sage.version

sys.modules['sage.all__sagemath_singular'] = sage.all
print('[AUDIT_RUNTIME] '+json.dumps(dict(python=sys.version.split()[0],
      sage=sage.version.version, flint=flint.__version__, sympy=sympy.__version__,
      compatibility_shim='sage.all__sagemath_singular -> sage.all')), flush=True)
source=Path(sys.argv[1]).resolve()
runpy.run_path(str(source), run_name='__main__')
