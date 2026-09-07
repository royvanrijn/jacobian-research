#!/usr/bin/env sage-python
"""V3 wrapper: preserves failed transfer v1/v2 and fixes preparation supervision path.

Builds on v2's cache-identity guard. New artifacts live under
v3-transfer-11952-v3. The only additional change is that preparation logs and
supervisor records are derived from the active campaign directory instead of
the preserved v1 literal path.
"""
from importlib.machinery import SourceFileLoader
from pathlib import Path
import os
import sys

SELF = Path(__file__).resolve()
CAS = SELF.parent
v2 = SourceFileLoader('v3_transfer_campaign_v2_preserved', str(CAS/'v3_transfer_campaign_v2.sage')).load_module()
base = v2.base
base.D = base.ROOT/'artifacts/local/elliptic-curves/v3-transfer-11952-v3'
base.__file__ = str(SELF)

# Bind v1, v2 and this wrapper into the new protocol.
_prev_sources = base.driver_sources
def driver_sources():
    rows = dict(_prev_sources())
    for p in (CAS/'v3_transfer_campaign_v2.sage', SELF):
        rows[str(p.relative_to(base.ROOT))] = base.sha(p)
    return rows
base.driver_sources = driver_sources

# v1 hard-coded `v3-transfer-11952-v1-preparation` for prepare-worker. That
# makes preserved attempts collide across wrappers. Keep case supervision
# unchanged, but derive preparation storage from the active D.
_original_supervise = base.supervise
def supervise(action, case=None, v3_replay=None):
    if action != 'prepare-worker':
        return _original_supervise(action, case=case, v3_replay=v3_replay)
    from research_runtime.supervisor import Limits, run
    logdir = base.D.parent/(base.D.name+'-preparation')
    logdir.mkdir(parents=True, exist_ok=True)
    stem = 'prepare'
    base.require(not (logdir/(stem+'.supervisor.json')).exists(),
                 'Preserve previous supervised preparation attempt. Review before retry.')
    command = [sys.executable, str(SELF), action]
    if v3_replay:
        command += ['--v3-replay', str(v3_replay.resolve())]
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
               MKL_NUM_THREADS='1', V3_TRANSFER_SUPERVISED='1')
    report = run(command,
        limits=Limits(base.LIMITS['prepare_wall_seconds'], base.LIMITS['rss_bytes']),
        cwd=base.ROOT, env=env, log_path=logdir/(stem+'.log'),
        checkpoint_path=logdir/(stem+'.supervisor.json'))
    base.require(report['outcome'] == 'completed',
                 f"prepare-worker stopped: {report['outcome']}; inspect {logdir/'prepare.log'}")
base.supervise = supervise

if __name__ == '__main__':
    base.main()
