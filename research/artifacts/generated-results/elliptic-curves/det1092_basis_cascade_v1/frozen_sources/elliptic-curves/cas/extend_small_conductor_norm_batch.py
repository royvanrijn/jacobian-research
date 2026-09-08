#!/usr/bin/env python3
"""Fixed 512-box continuation, with product-tree smoothness and scalar replay."""
import argparse
import json
from math import gcd, prod
from pathlib import Path
import sys
import time
import pilot_small_conductor_norm_smoothness as prior
from research_runtime.store import checkpoint
from research_runtime.supervisor import Limits, run

forms, cert = prior.forms, prior.forms.target.original.cert
ROOT, ART = forms.ROOT, forms.ART
D = ROOT / 'artifacts/local/elliptic-curves/small-conductor-norm-batch-v1'
OUT = ART / 'small_conductor_norm_batch_v1.json'


def protocol():
    p = cert.read(D / 'protocol.json')
    for name, digest in p['sources'].items():
        if cert.hashed(ROOT / name) != digest:
            raise ArithmeticError('frozen input changed: ' + name)
    return p


def prepare():
    if (D / 'protocol.json').exists():
        raise FileExistsError('preserve protocol')
    proof = ART / 'small_conductor_norm_relations_v1.json'
    if cert.read(proof)['additional_relation_rank'] != 18:
        raise ArithmeticError('previous productive pilot required')
    checkpoint(D / 'protocol.json', {
        'schema': 'elliptic-curves.small-conductor-norm-batch-protocol.v1',
        'sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in
                    [Path(__file__).resolve(), forms.OUT, proof, Path(prior.__file__).resolve()]},
        'box': 512, 'smooth_bound': 400000, 'm_per_chunk': 16,
        'workers': 1, 'rss_bytes': 1610612736,
        'stage_seconds': {'worker': 180, 'check': 300, 'audit': 180, 'audit-check': 180},
        'gate': 'Prior fixed pilot produced18 exact noncanonical relation rows. Test sustained yield on only this field and measure the full factor-base matrix gap.',
        'population': 'Every coprime(m,n) with -512<=m<=512 and1<=n<=512, in lexicographic order. The old64-box is included to measure incremental rank exactly.',
        'method': 'Product/remainder tree computes primorial modulo each norm; exact gcd removes all supported prime powers. Replay uses scalar gcd with the primorial on every value.',
        'factor_base': 'All prime ideals above every rational prime<=400000, with canonical rational-prime relations; generation is justified only under the stated Bach GRH hypothesis.',
        'claim_boundary': 'Finite relation experiment. No unconditional generation, new rational point, full descent, or exact rank claim. The measured quotient gives a class-2-rank upper bound only under GRH.'})


def residues(values, modulus):
    """Exact simultaneous modulus % values via a product/remainder tree."""
    levels = [values]
    while len(levels[-1]) > 1:
        old = levels[-1]
        levels.append([prod(old[i:i+2]) for i in range(0, len(old), 2)])
    current = [modulus % levels[-1][0]]
    for level in reversed(levels[:-1]):
        current = [current[i//2] % value for i, value in enumerate(level)]
    return current


def strip(value, support):
    while True:
        common = gcd(value, support)
        if common == 1:
            return value
        value //= common


def factor_smooth(value, primes):
    factors = []
    for p in primes:
        exponent = 0
        while value % p == 0:
            value //= p
            exponent += 1
        if exponent:
            factors.append([p, exponent])
        if value == 1:
            return factors
    raise ArithmeticError('smooth factorization incomplete')


def calculate(check=False):
    p = protocol()
    primes = prior.primes_to(p['smooth_bound'])
    primorial = prod(primes)
    c = list(map(int, cert.read(forms.OUT)['reduced_binary_cubic_descending']))
    chunks, smooth = [], []
    started = time.monotonic()
    count = 0
    for lo in range(-p['box'], p['box']+1, p['m_per_chunk']):
        hi = min(p['box']+1, lo+p['m_per_chunk'])
        pairs = [(m,n) for m in range(lo,hi) for n in range(1,p['box']+1) if gcd(m,n)==1]
        values = [sum(a*m**(3-i)*n**i for i,a in enumerate(c)) for m,n in pairs]
        if any(v == 0 for v in values):
            raise ArithmeticError('irreducible form vanishes')
        if check:
            remainders = [strip(abs(v), primorial) for v in values]
        else:
            supports = [gcd(abs(v),r) for v,r in zip(values,residues(list(map(abs,values)),primorial))]
            remainders = [strip(abs(v),s) for v,s in zip(values,supports)]
        records = []
        for (m,n), value, remainder in zip(pairs,values,remainders):
            record = {'m':m,'n':n,'value':str(value),'remainder':str(remainder)}
            if remainder == 1:
                record['factorization'] = factor_smooth(abs(value),primes)
                smooth.append(record)
            records.append(record)
        path = D / ('chunk_%04d.json' % len(chunks))
        data = {'m_start':lo,'m_stop_exclusive':hi,'records':records}
        if check:
            if cert.read(path) != data:
                raise ArithmeticError('scalar replay differs: '+str(path))
        else:
            if path.exists():
                raise FileExistsError('preserve completed chunks')
            checkpoint(path,data)
        chunks.append({'path':str(path.relative_to(ROOT)),'sha256':cert.hashed(path),'pairs':len(pairs)})
        count += len(pairs)
        if not check:
            checkpoint(D/'progress.json',{'completed_m':hi-1,'pairs':count,'smooth_values':len(smooth),'wall_seconds':time.monotonic()-started})
        if len(chunks)%8==0:
            print('CHUNKS',len(chunks),'PAIRS',count,'SMOOTH',len(smooth),flush=True)
    seconds = time.monotonic()-started
    record = {'schema':'elliptic-curves.small-conductor-norm-batch.v1','status':'PASS',
              'protocol':p,'chunks':chunks,'primitive_pairs':count,'smooth_records':smooth,
              'smooth_values':len(smooth),'old_box_smooth_values':sum(abs(r['m'])<=64 and r['n']<=64 for r in smooth),
              'wall_seconds':seconds,'claim_boundary':p['claim_boundary']}
    if check:
        old = cert.read(OUT)
        record['wall_seconds'] = old['wall_seconds']
        if record != old:
            raise ArithmeticError('batch summary differs')
        checkpoint(D/'scalar_replay.json',{'status':'PASS','source_sha256':cert.hashed(OUT),'wall_seconds':seconds})
    else:
        if OUT.exists():
            raise FileExistsError('preserve summary')
        checkpoint(OUT,record)
    print('NORM BATCH', 'CHECK' if check else 'BUILD', 'PASS',count,len(smooth),'SECONDS',seconds,flush=True)


def launch(stage):
    p = protocol()
    path = D/(stage+'.supervisor.json')
    if path.exists():
        raise FileExistsError('preserve supervisor evidence')
    command = [sys.executable,str(Path(__file__).resolve()),stage]
    if stage.startswith('audit'):
        command = ['/home/royvanrijn/.local/bin/sage','-python',str(ROOT/'elliptic-curves/cas/audit_small_conductor_norm_batch.sage')]
        if stage == 'audit-check':
            command.append('--check')
    result = run(command,limits=Limits(p['stage_seconds'][stage],p['rss_bytes']),cwd=ROOT,
                 log_path=D/(stage+'.log'),checkpoint_path=path)
    print(stage,result['outcome'],result['returncode'],flush=True)
    if result['outcome'] != 'completed' or result['returncode'] != 0:
        raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['prepare','worker','check','launch-worker','launch-check','launch-audit','launch-audit-check'])
    args=parser.parse_args()
    if args.stage=='prepare': prepare()
    elif args.stage.startswith('launch-'): launch(args.stage[7:])
    else: calculate(args.stage=='check')
