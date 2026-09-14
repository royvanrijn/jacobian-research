#!/usr/bin/env sage -python
"""Bounded enlarged-basis integration gate; no rational-point search.

One previously discovered fresh M18 is a labelled development input. Rebuild
its landscape from the retained generic parent subset, then independently
replay rational CVP and all selected anchor points. Charge the entire process.
"""
import argparse
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time

from finite_cancellation_corpus import ROOT, LOCAL, digest, write

CAS = ROOT/'elliptic-curves/cas'
OUT = ROOT/'artifacts/generated-results/elliptic-curves/cancellation_basis_refresh_v1'
RAW = LOCAL/'cancellation-basis-refresh-v1'
FRESH = ROOT/'artifacts/generated-results/elliptic-curves/cancellation_scheduler_fresh_v1'


def read(path): return json.loads(path.read_text())
def sha(path): return digest(path.read_bytes())
def need(value, message):
    if not value: raise ArithmeticError(message)
def cpu():
    r = resource.getrusage(resource.RUSAGE_CHILDREN)
    return time.process_time()+r.ru_utime+r.ru_stime


def freeze():
    from cancellation_scheduler_prepare import new_write
    from cancellation_scheduler_fresh import source_closure
    rows = [r for r in read(FRESH/'summary.json')['rows'] if r['arm'] == 'adaptive_local']
    # Highest already-certified fresh cloud, then case hash. This is explicitly
    # an integration/cost control, not a new predictive holdout.
    row = min(rows, key=lambda r: (-r['cloud_rank_lower_bound'], r['case']))
    source = LOCAL/'cancellation-scheduler-fresh-v1/arms'/row['case']/'adaptive_local'
    seed = read(source/'cloud-certificate.json'); generic = read(source/'cold/bank-00/seed-M16.json')
    bank = read(source/'cold/bank-00/anchor-bank.json')
    need(len(seed['points']) == 18 and seed['points'][:16] == generic['points'], 'development basis differs')
    need(read(source/'cloud-verification.json')['certificate_sha256'] == sha(source/'cloud-certificate.json'), 'cloud certification differs')
    inputs = dict(case=row['case'], seed=seed, generic_points=generic['points'], bank=bank,
        provenance_sha256={str(p.relative_to(ROOT)):sha(p) for p in (
            source/'cloud-certificate.json', source/'cloud-verification.json', source/'cold-receipt.json',
            source/'cold/bank-00/anchor-bank.json', source/'cold/bank-00/generic-cvp-proofs.json',
            source/'cold/bank-00/seed-M16.json', FRESH/'summary.json', FRESH/'replay-manifest.json')})
    new_write(OUT/'input.json', inputs)
    policy = dict(generic_rank=16, scaled_shells=bank['shells'], anchor_count=len(bank['rows']),
        anchors_per_shell=16, canonical_per_shell=25, exact_cvp_node_limit=2000000)
    prior = read(ROOT/'artifacts/generated-results/elliptic-curves/cancellation_scheduler_v4/protocol.json')
    sources = source_closure([ROOT/p for p in prior['source_sha256']] + [Path(__file__),
        CAS/'visibility_complement_subset.py', CAS/'run_complement_seed_v3.py', CAS/'verify_finite_cancellation_cpu.sage'])
    new_write(OUT/'protocol.json', dict(status='FROZEN_BASIS_REFRESH_INTEGRATION',
        input_sha256=sha(OUT/'input.json'), source_sha256=sources, policy=policy,
        hard_process_cpu_seconds=60, wall_seconds=90, point_search_calls=0,
        prospective_search_cpu_seconds=40, maximum_refresh_fraction=.25,
        selection='Highest previously certified adaptive cloud rank in the four-fibre pilot, then case hash; bank00 only. Labelled development/cost input, not an independent efficacy test.',
        gate='Exact producer/reference landscape equality and independently certified starting subgroup and anchor identities; complete isolated CPU must be at most10 seconds, a quarter of the existing40-second search allowance. Passing permits designing a separately frozen basis-aware control, not a fresh search or campaign restart.',
        mathematical_boundary='The input subgroup is certified. CVP is exact for the recorded rounded positive numerical height metric; this is not a certified canonical-height matrix. Only the retained16-row generic subset is extended, not the full generic shell.',
        continuation='Any later epoch starts from the entire reconciled certified cloud, discards the previous box exposures, and requires its own compatible verified landscape. This gate makes no point calls.'))
    print(json.dumps(dict(status='FROZEN_BASIS_REFRESH_INTEGRATION', case=row['case'], rank=18,
        generic_anchors=len(bank['rows']), extension_classes=len(bank['rows'])*4)), flush=True)


def guard():
    plan = read(OUT/'protocol.json'); inputs = read(OUT/'input.json')
    need(sha(OUT/'input.json') == plan['input_sha256'], 'input changed')
    for name,h in {**plan['source_sha256'], **inputs['provenance_sha256']}.items():
        need(sha(ROOT/name) == h, 'bound bytes changed: '+name)
    return plan, inputs


def worker():
    plan, inputs = guard()
    resource.setrlimit(resource.RLIMIT_CPU, (plan['hard_process_cpu_seconds'],)*2)
    need(not (RAW/'verification.json').exists(), 'preserve verification')
    from sage.all import EllipticCurve, QQ, pari
    from sage.env import SAGE_VERSION
    from visibility_complement_subset import landscape
    from run_complement_seed_v3 import normalized_selection
    from v3_warm_engine import certified_state
    seed = inputs['seed']; model = tuple(map(F, seed['curve'])); basis = tuple(tuple(map(F,p)) for p in seed['points'])
    tick = cpu(); certified_state(model, basis, seed['proof'])
    rank_checker = SourceFileLoader('refresh_independent_rank', str(CAS/'verify_finite_cancellation_cpu.sage')).load_module().sage_rank
    independent_rank = rank_checker(seed)
    need(independent_rank['rank'] == len(basis), 'input rank not independently certified')
    write(RAW/'rank-gate.json', dict(status='PASS', rank=independent_rank, cpu_seconds=cpu()-tick))
    tick = cpu(); producer = landscape(model,basis,set(),RAW/'producer',plan['policy'],inputs['bank'])
    production_cpu = cpu()-tick
    tick = cpu(); reference = landscape(model,basis,set(),RAW/'reference',plan['policy'],inputs['bank'],reference=True)
    reference_cpu = cpu()-tick
    need(normalized_selection(producer) == normalized_selection(reference), 'rational CVP replay differs')
    tick = cpu(); E = EllipticCurve(QQ,list(model)); native = [E(list(p)) for p in basis]
    for centre in producer['centres']:
        P = sum((int(c)*point for c,point in zip(centre['representative'],native)),E(0))
        need(not P.is_zero() and list(map(str,P.xy())) == centre['point'], 'independent anchor differs')
    native_cpu = cpu()-tick
    write(OUT/'prepared-bank.json', dict(status='VERIFIED_ENLARGED_SUBGROUP_BANK', seed=seed,
        centres=[r['representative'] for r in producer['centres'][:48]],
        input_sha256=sha(OUT/'input.json'), selection_sha256=sha(RAW/'producer/selection.json'),
        boundary='Certified M18 basis with an independently replayed landscape for16 retained generic anchors and all4 extension classes each. A later basis change invalidates compatibility. No point search has run on this bank.'))
    write(RAW/'verification.json', dict(status='PASS_NEW_BASIS_LANDSCAPE', rank=independent_rank,
        generic_anchors=len(inputs['bank']['rows']), full_cosets_scored=producer['full_cosets_scored'],
        selected_centres=len(producer['centres']), exported_centres=min(48,len(producer['centres'])),
        production_cpu_seconds=production_cpu, reference_cpu_seconds=reference_cpu,
        native_anchor_cpu_seconds=native_cpu, total_internal_cpu_seconds=cpu(),
        runtime=dict(sage=SAGE_VERSION, pari=str(pari('version()')), python=sys.version),
        protocol_sha256=sha(OUT/'protocol.json'), point_search_calls=0,
        selection_sha256=sha(RAW/'producer/selection.json'), reference_sha256=sha(RAW/'reference/selection.json'),
        prepared_bank_sha256=sha(OUT/'prepared-bank.json')))


def run():
    plan,_ = guard(); need(not (OUT/'supervision.json').exists() and not (RAW/'start.json').exists(), 'preserve started attempt')
    RAW.mkdir(parents=True,exist_ok=True); write(RAW/'start.json', dict(protocol_sha256=sha(OUT/'protocol.json'), status='STARTED'))
    start=time.monotonic(); before=resource.getrusage(resource.RUSAGE_CHILDREN)
    with (RAW/'worker.log').open('w') as log:
        proc=subprocess.Popen([sys.executable,str(Path(__file__).resolve()),'worker'],stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        timed_out=False
        try: code=proc.wait(timeout=plan['wall_seconds'])
        except subprocess.TimeoutExpired:
            timed_out=True; os.killpg(proc.pid,signal.SIGKILL); code=proc.wait()
    after=resource.getrusage(resource.RUSAGE_CHILDREN); charged=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime
    complete=code==0 and (RAW/'verification.json').exists()
    record=dict(status='COMPLETE' if complete else 'UNKNOWN_INCOMPLETE_RETAINED', returncode=code,
        wall_seconds=time.monotonic()-start, charged_cpu_seconds=charged, timed_out=timed_out,
        protocol_sha256=sha(OUT/'protocol.json'), verification_sha256=sha(RAW/'verification.json') if complete else None,
        cost_gate_passed=complete and charged<=plan['prospective_search_cpu_seconds']*plan['maximum_refresh_fraction'],
        point_search_calls=0, boundary='One development integration/cost input. No efficacy or fresh-point claim. Failed limits do not exclude directions or trigger a rerun.')
    write(OUT/'supervision.json',record); print(json.dumps(record,indent=2),flush=True)
    need(complete, 'integration gate incomplete; checkpoint retained')


def check():
    plan,inputs=guard(); record=read(OUT/'supervision.json')
    need(record['status']=='COMPLETE' and record['verification_sha256']==sha(RAW/'verification.json'), 'worker receipt differs')
    v=read(RAW/'verification.json')
    need(v['status']=='PASS_NEW_BASIS_LANDSCAPE' and v['point_search_calls']==0, 'integration scope differs')
    for folder,key in (('producer','selection_sha256'),('reference','reference_sha256')):
        need(sha(RAW/folder/'selection.json')==v[key], 'selection bytes differ')
    need(sha(OUT/'prepared-bank.json')==v['prepared_bank_sha256'], 'exported bank differs')
    need(record['cost_gate_passed']==(record['charged_cpu_seconds']<=10), 'cost gate differs')
    print(json.dumps(dict(status='PASS_RETAINED_BASIS_REFRESH_RECEIPTS', **{k:v[k] for k in ('generic_anchors','full_cosets_scored','selected_centres','point_search_calls')},
        charged_cpu_seconds=record['charged_cpu_seconds'],cost_gate_passed=record['cost_gate_passed'])),flush=True)


def pack():
    check()
    import pack_cancellation_scheduler as bundle
    bundle.OUT=OUT; bundle.RAW=RAW; bundle.pack()


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('command',choices=['freeze','worker','run','check','pack'])
    globals()[p.parse_args().command]()
