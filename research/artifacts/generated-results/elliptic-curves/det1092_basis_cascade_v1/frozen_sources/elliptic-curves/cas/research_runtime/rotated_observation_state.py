"""Archive observation metadata outside the next live chart state.

The caller must persist the returned full archive before using the new state.
MWState.adjoin uses observations only to append history and to derive a parent
content key; point membership, basis, reductions, torsion and height decisions
use the other fields, all preserved here. The original source hash is enforced.
"""
from dataclasses import replace
from hashlib import sha256
from pathlib import Path
from . import mw_state
PIN='e53c24a786afe1096a3bc90eff94396b1b3b6b1fd3ec52c24e51c69622b4596d'
if sha256(Path(mw_state.__file__).read_bytes()).hexdigest()!=PIN:raise RuntimeError('observation rotation needs reviewed MWState source')

def rotate(state):
    archive=state.record()
    return replace(state,observations=(),parent_state=archive['key']),archive
