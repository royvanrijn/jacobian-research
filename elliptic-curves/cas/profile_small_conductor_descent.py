#!/usr/bin/env python3
"""Bounded, staged descent diagnostics; no incomplete BNF is an upper bound."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import certify_small_conductor_curve as proof_checker
from research_runtime.store import checkpoint
from research_runtime.supervisor import Limits, run

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
D = ROOT / 'artifacts/local/elliptic-curves/small-conductor-descent-profile-v1'
PROOF = ART / 'small_conductor_rank22_proof_v1.json'
GP = Path('/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/gp')
HISTORY = [ROOT / 'elliptic-curves/notes/R17_MW29_RELATIVE_SCLASS_RETRY_2026-09-04.md',
           ROOT / 'elliptic-curves/notes/BNF_FREE_RESIDUAL_2SELMER.md']


def prepare():
    if (D / 'protocol.json').exists():
        raise FileExistsError('preserve profile protocol')
    proof = proof_checker.cert.read(PROOF)
    proof_checker.verify(proof)
    a, b = map(int, proof['integral_model'][3:])
    f = f'x^3+x^2+({16*a})*x+({64*b})'
    primes = '[' + ','.join(p for p, _ in proof['discriminant_factorization']) + ']'
    curve = '[' + ','.join(proof['integral_model']) + ']'
    common = f'default(realprecision,80);setrand(1);f={f};E=ellinit({curve});P={primes};\n'
    programs = {}
    programs['field'] = common + f'''
addprimes(P);gettime();print("STAGE|nfinit_start");
nf=nfinit([f,P]);print("MS|nfinit|",gettime());
print("FIELD_DISCRIMINANT|",nf.disc);print("FIELD_SIGNATURE|",nf.sign);
print("FIELD_INDEX|",nf.index);print("FIELD_BASIS|",nf.zk);
print("CERTIFY_ORDER|",nfcertify(nf));print("MS|order_certify|",gettime());
writebin("{D / 'nf.bin'}",nf);
R=polredabs(nf,1);print("REDUCED|",R);print("MS|polredabs|",gettime());
write("{D / 'reduced.gp'}",R);
for(i=1,#P,Q=idealprimedec(nf,P[i]);print("SPLIT|",P[i],"|",vector(#Q,j,[Q[j].e,Q[j].f])));
print("DONE|field");quit
'''
    for name, hinted in [('cold_rankinit', False), ('hinted_rankinit', True)]:
        programs[name] = common + ('addprimes(P);\n' if hinted else '') + f'''
setdebug("bnf",3);gettime();print("STAGE|ellrankinit_start");
r=ellrankinit(E);print("MS|ellrankinit|",gettime());
writebin("{D / (name + '.bin')}",r);print("DONE|{name}");quit
'''
    programs['prepared_bnf'] = common + f'''
addprimes(P);nf=read("{D / 'nf.bin'}");setdebug("bnf",3);
gettime();print("STAGE|bnfinit_start");b=bnfinit(nf,0);
print("MS|bnfinit|",gettime());print("COMPUTED_CYCLIC_FACTORS|",b.cyc);
writebin("{D / 'uncertified_bnf.bin'}",b);
print("STAGE|bnfcertify_quotient_start");C=bnfcertify(b,1);
print("CLASS_QUOTIENT_CERTIFIED|",C);print("MS|bnfcertify_quotient|",gettime());
print("DONE|prepared_bnf");quit
'''
    for name, program in programs.items():
        (D / (name + '.gp')).parent.mkdir(parents=True, exist_ok=True)
        (D / (name + '.gp')).write_text(program)
    inputs = [PROOF, *HISTORY, Path(__file__).resolve(), ROOT / 'elliptic-curves/cas/research_runtime/supervisor.py', GP]
    checkpoint(D / 'protocol.json', {
        'schema': 'elliptic-curves.small-conductor-descent-profile.v1',
        'inputs': {str(p): proof_checker.cert.hashed(p) for p in inputs},
        'gp_version': subprocess.run([str(GP), '--version-short'], capture_output=True, text=True, check=True).stdout.strip(),
        'programs': {name: proof_checker.cert.hashed(D / (name + '.gp')) for name in programs},
        'limits_seconds': {'field': 30, 'cold_rankinit': 30, 'hinted_rankinit': 30, 'prepared_bnf': 60},
        'rss_bytes': 1610612736, 'pari_stack_bytes': 512000000, 'maximum_workers': 1,
        'selection': 'Only the user-selected MW16 family05 fibre3/17. No previous-curve sweep.',
        'gate': 'Audit prior generic and BNF-free failures first. Profile maximal-order setup with certified factor hints; compare short cold/hinted initialization solely to locate cost. At most one prepared BNF flag0 canary tests whether skipping factorization and later unit certification changes the front-end bottleneck.',
        'claim_boundary': 'Only a completed exact order calculation or separately certified class quotient can support its corresponding claim. No timeout, incomplete relation count, or computed but uncertified class group is a Selmer/rank upper bound. No point search.'})
    print('PREPARED fixed profile:30+30+30+60 seconds, one worker')


def execute(stage):
    protocol = proof_checker.cert.read(D / 'protocol.json')
    if any(proof_checker.cert.hashed(Path(p)) != h for p, h in protocol['inputs'].items()):
        raise ArithmeticError('profile input changed')
    path = D / (stage + '.gp')
    if proof_checker.cert.hashed(path) != protocol['programs'][stage]:
        raise ArithmeticError('GP input changed')
    out = D / (stage + '.json')
    if out.exists():
        raise FileExistsError('preserve profile result')
    if stage == 'prepared_bnf':
        parent = proof_checker.cert.read(D / 'field.json')
        if parent['supervision']['outcome'] != 'completed' or 'CERTIFY_ORDER|[]' not in (D / 'field.log').read_text():
            raise ArithmeticError('certified field preflight required')
    result = run([str(GP), '-q', '-f', '-s', str(protocol['pari_stack_bytes']), str(path)],
                 limits=Limits(protocol['limits_seconds'][stage], protocol['rss_bytes']),
                 log_path=D / (stage + '.log'), checkpoint_path=D / (stage + '.supervisor.json'), cwd=ROOT)
    text = (D / (stage + '.log')).read_text()
    checkpoint(out, {'schema': 'elliptic-curves.small-conductor-profile-stage.v1',
                     'stage': stage, 'supervision': result,
                     'program_sha256': protocol['programs'][stage],
                     'log_sha256': proof_checker.cert.hashed(D / (stage + '.log')),
                     'gp_done_marker': 'DONE|' + stage in text,
                     'gp_error': '***' in text, 'unconditional_rank_upper_bound': None})
    print(stage, result['outcome'], result['wall_seconds'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['prepare', 'field', 'cold_rankinit', 'hinted_rankinit', 'prepared_bnf'])
    args = parser.parse_args()
    prepare() if args.stage == 'prepare' else execute(args.stage)
