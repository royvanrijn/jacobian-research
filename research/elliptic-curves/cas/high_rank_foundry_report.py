"""Human-readable live view of certified endpoints; no scheduling decisions."""
from datetime import datetime,timezone
import os
from pathlib import Path
import tempfile


def write_report(path,status,config):
    timestamp=datetime.fromtimestamp(status['updated_at'],timezone.utc).isoformat()
    scope='this run and its preserved predecessor' if status.get('inherited_from') else 'this run'
    lines=['# High-rank search foundry', '', f"**{status['status']}** · updated {timestamp}", '',
           f"{status['new_curves']} new certified curves across {scope}; {status['point_calls']} completed point calls. "
           f"Jobs: {status['jobs']}. Up to {status['workers']} workers.", '',
           'Ranks are certified lower bounds. Imported historical packets are marked below. '
           'Novelty is relative to the frozen catalogue; exact rank and worldwide record status remain unknown.', '',
           '| Curve | Family | Parameter | Certified lower bound | Calls since gain | State | Certificate |',
           '| --- | --- | --- | ---: | ---: | --- | --- |']
    frozen=Path(config['frozen_root'])
    for c in status['best']:
        label=c['id']+(' (imported)' if c.get('historical_baseline') else '')
        packet=c.get('certificate') or frozen/c['packet']
        lines.append(f"| {label} | {c['family']} | {c['parameter']} | {c['rank']} | {c['stale_calls']} | {c['state']} | [points and proof]({packet}) |")
    lines+=['','Active jobs:']
    for j in status['active_jobs']:
        lines.append(f"- {j['id']}: {j['cid']} ({j['kind']}); [job directory]({frozen/j['path']})")
    lines+=['','[Full machine-readable status](LIVE_STATUS.json). '
            'The accompanying immutable job exports contain portable certificates and complete trajectories. '
            'Family/selection-arm aggregates are descriptive and adaptively censored.','']
    path=Path(path);fd,temp=tempfile.mkstemp(prefix='.report-',dir=path.parent)
    try:
        with os.fdopen(fd,'w') as stream:
            stream.write('\n'.join(lines));stream.flush();os.fsync(stream.fileno())
        os.replace(temp,path)
    finally:
        if os.path.exists(temp):os.unlink(temp)
