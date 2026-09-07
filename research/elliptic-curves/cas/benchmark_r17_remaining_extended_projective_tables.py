#!/usr/bin/env python3
"""Balanced fifteen-table cost gate for longer scoring before R17 retention."""
import argparse
import time
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
import benchmark_extended_projective_trace_cache_v2 as engine
from research_runtime.store import checkpoint, digest
from research_runtime.supervisor import capture, Limits

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / 'elliptic-curves/cas'
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
D = ROOT / 'artifacts/local/elliptic-curves/r17_remaining-extended-projective-benchmark-v1'
OUT = ART / 'r17_remaining_extended_projective_benchmark_v1.json'

def sources():
    paths = [Path(__file__).resolve(), spec.ATLAS, Path(spec.__file__),
             Path(engine.__file__), Path(engine.scalar.__file__), Path(cert.__file__),
             CAS/'mod2_reduction_independence.py', CAS/'research_runtime/store.py',
             CAS/'research_runtime/supervisor.py',
             ART/'r17_extended_score_prime_minimality_v1.json',
             ART/'r17_extended_score_prime_minimality_sage_replay_v1.json',
             ART/'broad60_mw16_experiment_v1.json']
    return {str(p.relative_to(ROOT)): cert.hashed(p) for p in paths}

def prepare():
    if (D/'protocol.json').exists():
        raise FileExistsError('preserve fifteen-table protocol')
    for name in ['r17_extended_score_prime_minimality_v1.json', 'r17_extended_score_prime_minimality_sage_replay_v1.json', 'broad60_mw16_experiment_v1.json']:
        if cert.read(ART/name)['status'] != 'PASS':
            raise ArithmeticError('complete arithmetic and population gates required')
    primes = [p for p in engine._primes_up_to(32749) if p >= 4099]
    cases = [primes[0], primes[(len(primes)-1)//2], primes[-1]]
    families = []
    for row in cert.read(spec.ATLAS)['families']:
        if row['family'] == '11952':
            continue
        model = {k: row[k]+['0']*(n-len(row[k])) for k,n in
                 [('A_coefficients_low_to_high',9),('B_coefficients_low_to_high',13)]}
        if any(cert.F(c).denominator != 1 for a in model.values() for c in a):
            raise ArithmeticError('integer atlas coefficients required')
        families.append({'family':row['family'], 'model':model, 'model_hash':digest(model)})
    if len(families) != 5 or len({r['family'] for r in families}) != 5:
        raise ArithmeticError('all five independent atlas presentations required')
    checkpoint(D/'protocol.json', {
        'schema':'elliptic-curves.r17_remaining-extended-projective-benchmark.v1',
        'sources':sources(), 'gp_sha256':cert.hashed(engine.scalar.GP),
        'families':families, 'primes':cases, 'full_prime_roster':primes,
        'gp_seconds_per_case':20, 'rss_bytes':536870912, 'outer_seconds':300,
        'projected_serial_seconds_per_family_gate':1800,
        'gate':'The user requests broader earlier selection on new territory. The completed broader MW16 cohort changes initial coverage while preserving downstream budgets and produces three new high-rank inventory curves. Five of the six existing compact R17 families still have only short-prime projective tables and scalar extension after short retention;11952 already has its complete extended cache. Test the remaining five equally, without a new catalogue calibration or target input. The exact six-family arithmetic audit excludes removable scaling at every prime4099through32749. Complete height-independent extended tables would enable a later new-territory scan to apply the existing longer score before discarding candidates. Generic sections and fibration equations are already certified; this benchmark changes neither.',
        'scope':'Exactly three predetermined prime tables per family, all five remaining families equally represented;11952 is excluded because its complete cache already exists; compute every projective residue, including infinity. One worker,20 seconds per case and300 seconds total; preserve failure without retry. Verify complete framing, singularity markers and Hasse bounds, plus five independent character sums per table. No catalogue, known record equation/parameter/point/rank/j-invariant or jump label is read. A full cache or fresh parameter/point campaign requires a separate finite protocol after this cost gate; this benchmark launches none. Runtime extrapolations are descriptive, not mathematical bounds.'})

def protocol():
    p = cert.read(D/'protocol.json')
    if p['sources'] != sources() or p['gp_sha256'] != cert.hashed(engine.scalar.GP):
        raise ArithmeticError('frozen benchmark inputs changed')
    return p

def table(f, prime, p, create):
    folder = D/f['family']/str(prime)
    rawpath = folder/'raw.json'
    code = engine.program(f['model'], prime)
    if create:
        if rawpath.exists():
            raise FileExistsError('preserve attempted table')
        result = capture([str(engine.scalar.GP),'-q','-s','256000000'],
                         input_text=code, limits=Limits(p['gp_seconds_per_case'],p['rss_bytes']),
                         log_path=folder/'gp.log', separate_stderr=True, check=False)
        checkpoint(rawpath, {'program':code,'stdout':result.stdout,'stderr':result.stderr,
                             'supervision':result.supervision})
    raw = cert.read(rawpath)
    if raw['program'] != code or raw['stderr'] or raw['supervision']['outcome'] != 'completed' or raw['supervision']['returncode'] != 0:
        raise ArithmeticError('table failed or censored')
    lines = raw['stdout'].splitlines()
    if len(lines) != prime+3 or lines[-1] != 'DONE' or not lines[-2].startswith('MS|'):
        raise ArithmeticError('incomplete projective table')
    a,b = ([int(c)%prime for c in f['model'][k]] for k in
           ('A_coefficients_low_to_high','B_coefficients_low_to_high'))
    values = []
    for i,line in enumerate(lines[:-2]):
        aa = a[-1] if i == prime else engine.evaluate(a,i,prime)
        bb = b[-1] if i == prime else engine.evaluate(b,i,prime)
        bad = (4*aa**3+27*bb**2)%prime == 0
        fields = line.split('|')
        if len(fields) != (2 if bad else 3) or fields[0] != ('B' if bad else 'T') or int(fields[1]) != i:
            raise ArithmeticError('residue or discriminant mismatch')
        value = None if bad else int(fields[2])
        if value is not None and value*value > 4*prime:
            raise ArithmeticError('Hasse bound violated')
        values.append(value)
    direct = []
    for i in [0,prime//3,2*prime//3,prime-1,prime]:
        aa = a[-1] if i == prime else engine.evaluate(a,i,prime)
        bb = b[-1] if i == prime else engine.evaluate(b,i,prime)
        value = engine.scalar.direct([0,0,0,aa,bb],prime)
        if value != values[i]:
            raise ArithmeticError('independent character sum differs')
        direct.append([i,value])
    saved = {'input':{'family':f['family'],'model_hash':f['model_hash'],'prime':prime},
             'traces':[v or 0 for v in values], 'good':[v is not None for v in values]}
    path = folder/'table.json'
    if create:
        if path.exists(): raise FileExistsError('preserve prime table')
        checkpoint(path,saved)
    elif cert.read(path) != saved:
        raise ArithmeticError('table replay differs')
    return {'family':f['family'],'prime':prime,'projective_rows':prime+1,
            'wall_seconds':raw['supervision']['wall_seconds'],'direct_checks':direct,
            'raw_sha256':cert.hashed(rawpath),'table_sha256':cert.hashed(path)}

def run(check=False):
    p = protocol()
    if not check and (D/'result.json').exists(): raise FileExistsError('preserve benchmark attempt')
    data = {'status':'RUNNING','protocol_sha256':cert.hashed(D/'protocol.json'),'rows':[]}
    start = time.monotonic()
    try:
        for f in p['families']:
            for prime in p['primes']:
                if time.monotonic()-start > p['outer_seconds']: raise TimeoutError('fixed outer deadline')
                row = table(f,prime,p,not check);data['rows'].append(row)
                if not check: checkpoint(D/'result.json',data)
                print('R17 PROJECTIVE TABLE',f['family'],prime,row['wall_seconds'],'PASS',flush=True)
        full_rows = sum(q+1 for q in p['full_prime_roster'])
        projections = {f['family']:max(r['wall_seconds']/r['projective_rows'] for r in data['rows'] if r['family']==f['family'])*full_rows for f in p['families']}
        data.update(status='PASS',projected_serial_seconds=projections,
                    cost_gate_passed=all(t<=p['projected_serial_seconds_per_family_gate'] for t in projections.values()),
                    direct_character_sums=75)
        artifact = {'schema':'elliptic-curves.r17_remaining-extended-projective-benchmark-result.v1',**data,'sources':sources(),'claim_boundary':p['scope']}
        if check:
            if cert.read(D/'result.json') != data or cert.read(OUT) != artifact: raise ArithmeticError('benchmark replay differs')
        else:
            if OUT.exists(): raise FileExistsError('preserve benchmark certificate')
            checkpoint(D/'result.json',data);checkpoint(OUT,artifact)
    except Exception as exc:
        if not check:
            data.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(D/'result.json',data)
        raise
    print('FIVE-FAMILY CACHE COST GATE',data['cost_gate_passed'],projections,flush=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage',choices=['prepare','run','check'])
    a = parser.parse_args()
    prepare() if a.stage == 'prepare' else run(a.stage == 'check')
