#!/usr/bin/env sage -python
"""Bounded negative replays on copies of one completed basis-refresh arm."""
import json
from pathlib import Path
import resource
import shutil
import signal
import tempfile
import time

from cancellation_basis_amplification import OUT, RAW, guard
from finite_cancellation_corpus import canonical, digest, write
from finite_cancellation_features import alarm
from verify_cancellation_basis_amplification import verify_arm


def read(path): return json.loads(path.read_text())
def sha(path): return digest(path.read_bytes())


def main():
    plan, cases = guard()
    if read(OUT/'supervision.json')['status'] != 'COMPLETE':
        raise ArithmeticError('finish timed arms before supplementary replays')
    if (OUT/'transition-regression.json').exists():
        raise FileExistsError('preserve transition replay')
    resource.setrlimit(resource.RLIMIT_CPU, (30, 30))
    signal.signal(signal.SIGALRM, alarm); signal.alarm(60)
    recipe = read(RAW/'transition-regression-recipe.json')
    case = next(c for c in cases if c['id'] == recipe['case'])
    source = RAW/'arms'/case['id']/recipe['arm']; result = read(source/'result.json')
    tick = time.process_time()
    positive = verify_arm(case, recipe['arm'], source, plan, result)
    checks = []
    for test, expected in [('missing_rebuild', 'search continued on obsolete basis'),
                           ('discarded_cloud', 'bank does not use complete cloud')]:
        with tempfile.TemporaryDirectory(prefix='basis-transition-replay-') as tmp:
            dest = Path(tmp)/'arm'; shutil.copytree(source, dest)
            events = read(dest/'events.json')
            index = next(i for i, e in enumerate(events) if e['kind'] == 'bank_start' and e['epoch'] > 0)
            if test == 'missing_rebuild':
                del events[index]
            else:
                events[index]['rank'] = len(case['seed']['points'])
                events[index]['seed_sha256'] = digest(canonical(case['seed']))
            write(dest/'events.json', events)
            # Rebind only the local test result so the semantic invariant,
            # rather than a trivial changed-file hash, must reject the copy.
            altered = {**result, 'events_sha256': sha(dest/'events.json')}
            try:
                verify_arm(case, recipe['arm'], dest, plan, altered)
            except ArithmeticError as error:
                if expected not in str(error):
                    raise
                checks.append({'test': test, 'status': 'REJECTED_EXPECTED_INVARIANT', 'reason': str(error)})
            else:
                raise ArithmeticError('invalid transition accepted: '+test)
    signal.alarm(0)
    report = {'status': 'PASS_POSITIVE_AND_TWO_NEGATIVE_TRANSITION_REPLAYS', 'case': case['id'],
              'positive_rank': positive['rank'], 'checks': checks, 'cpu_seconds': time.process_time()-tick,
              'source_sha256': sha(Path(__file__)), 'original_result_sha256': sha(source/'result.json'),
              'protocol_sha256': sha(OUT/'protocol.json'), 'point_search_calls': 0,
              'boundary': 'Supplementary verifier regressions after the fixed trial, on temporary copies. No experimental input, source seal, or outcome was altered.'}
    write(OUT/'transition-regression.json', report)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
