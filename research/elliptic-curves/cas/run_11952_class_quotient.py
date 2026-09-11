#!/usr/bin/env python3
"""Separate large-stack, fail-closed quotient-only rank certificate lane.

Stages are one-shot, preserve all failure logs, and never use known points in
BNF arithmetic. Provisional GRH and unconditional bounds have distinct fields.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import socket
import subprocess

from research_runtime.store import checkpoint
from research_runtime.supervisor import Limits, run

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
PRIOR = ART / '11952_rank_bound_audit_v1'
OUT = ART / '11952_class_quotient_v1'
LOCAL = ROOT / 'artifacts/local/elliptic-curves/11952-class-quotient-v1'
NF = ROOT / 'artifacts/local/elliptic-curves/11952-rank-bound-audit-v1/nf.bin'
CASE = ROOT / 'artifacts/local/elliptic-curves/broad-rank-v1/runtime/research/broad-cases/b-7bb187bc9254e81c6283/batch-002/search-00'
GP = Path('/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/gp')
STACK = 8 * 2**30
RSS = 12 * 2**30


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def read(p):
    return json.loads(p.read_text())


def memory_available():
    return int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines()
                    if line.startswith('MemAvailable:'))) * 1024


def fields(log):
    result = {}
    for line in log.splitlines():
        if '|' in line:
            key, value = line.split('|', 1)
            if key in result and key != 'STAGE':
                raise ValueError('duplicate output marker ' + key)
            result[key] = value
    return result


def classify(stage, log, outcome, binary_present, provisional=None):
    """Pure status gate, including GP's zero-exit-after-error pitfall."""
    result = {'status': 'INCOMPLETE', 'provisional_class_2rank': None,
              'grh_conditional_rank_upper': None, 'grh_conditional_exact_rank': None,
              'unconditional_class_2rank_upper': None, 'unconditional_rank_upper': None,
              'unconditional_exact_rank': None}
    if stage == 'quotient' and provisional:
        for key in ['provisional_class_2rank', 'grh_conditional_rank_upper', 'grh_conditional_exact_rank']:
            result[key] = provisional[key]
    f = fields(log)
    if outcome != 'completed' or 'FAIL' in f or any(line.lstrip().startswith('***') for line in log.splitlines()):
        result['status'] = 'RESOURCE_LIMIT' if outcome in ['strict_wall_timeout', 'rss_limit'] else 'BACKEND_FAILURE'
        return result
    if f.get('DONE') != stage:
        return result
    if stage == 'provisional':
        if not binary_present or f.get('CERTIFY_ORDER') != '[]':
            return result
        cyc = json.loads(f['PROVISIONAL_CYCLIC_FACTORS'])
        if not isinstance(cyc, list) or any(type(n) is not int or n < 2 for n in cyc):
            raise ArithmeticError('invalid provisional invariants')
        if any(cyc[i] % cyc[i+1] for i in range(len(cyc)-1)):
            raise ArithmeticError('invariant factor divisibility')
        g = sum(n % 2 == 0 for n in cyc)
        if g != int(f['PROVISIONAL_2RANK']):
            raise ArithmeticError('2-rank disagreement')
        result.update(provisional_cyclic_factors=list(map(str, cyc)), provisional_class_2rank=g)
        if g + 9 < 25:
            result['status'] = 'INCONSISTENT_WITH_CERTIFIED_LOWER_BOUND'
            return result
        result.update(grh_conditional_rank_upper=g+9,
                      grh_conditional_exact_rank=25 if g+9 == 25 else None,
                      status='PROVISIONAL_GRH_SUFFICIENT' if g <= 16 else 'NOT_SUFFICIENT')
    else:
        if not provisional or provisional['status'] != 'PROVISIONAL_GRH_SUFFICIENT':
            raise ArithmeticError('quotient gate requires sufficient provisional bound')
        g = provisional['provisional_class_2rank']
        if json.loads(f['PROVISIONAL_CYCLIC_FACTORS']) != list(map(int, provisional['provisional_cyclic_factors'])):
            raise ArithmeticError('BNF checkpoint invariants changed')
        if f.get('CLASSGROUP_QUOTIENT_CERTIFIED') != '1':
            result['status'] = 'QUOTIENT_NOT_CERTIFIED'
            return result
        if int(f['G_UPPER']) != g or g > 16:
            raise ArithmeticError('quotient dimension gate')
        result.update(status='CLASSGROUP_QUOTIENT_CERTIFIED', unconditional_class_2rank_upper=g,
                      unconditional_rank_upper=g+9, unconditional_exact_rank=25)
    return result


def prepare():
    if OUT.exists() or LOCAL.exists():
        raise FileExistsError('preserve all prior attempts')
    prior = read(PRIOR / 'field.json')
    assert prior['maximal_order_certified'] and prior['bound_offset'] == 9
    lower = read(CASE / 'verified.json')
    assert lower['status'] == 'PASS_INDEPENDENT_PRODUCTIVE_V3_REPLAY'
    assert lower['rank_lower_bound'] == 25 and digest(CASE / 'terminal.json') == lower['terminal_sha256']
    assert digest(CASE / 'terminal.json') == read(PRIOR / 'equation.json')['source_sha256']
    assert memory_available() >= RSS + 4*2**30, 'leave headroom for other work'
    OUT.mkdir(parents=True)
    LOCAL.mkdir(parents=True)
    common = f'default(realprecision,80);default(parisizemax,{STACK});default(nbthreads,1);setrand(1);\n'
    polynomial = 'x^3+x^2-4060408041437129659386958883089829522172534180*x+69302996254333102097575929282499740768182420645699427883330979661600'
    provisional = common + f'''
main()={{
  my(nf,b,g,C);
  nf=read("{NF}");
  if(nf.pol!={polynomial},error("cubic mismatch"));
  if(nf.disc!={prior['field_discriminant']},error("field discriminant mismatch"));
  C=nfcertify(nf);if(C!=[],error("maximal order not certified"));print("CERTIFY_ORDER|",C);
  print("STAGE|bnfinit");gettime();setdebug("bnf",1);
  b=bnfinit(nf,0);setdebug("bnf",0);print("BNFINIT_MS|",gettime());
  print("PROVISIONAL_CYCLIC_FACTORS|",b.cyc);
  g=sum(i=1,#b.cyc,b.cyc[i]%2==0);print("PROVISIONAL_2RANK|",g);
  writebin("{LOCAL / 'provisional_bnf.bin'}",b);
  print("DONE|provisional");
}};
iferr(main(),e,print("FAIL|",e);quit(1));quit(0)
'''
    quotient = common + f'''
main()={{
  my(b,g,C);
  b=read("{LOCAL / 'provisional_bnf.bin'}");
  if(b.nf.pol!={polynomial},error("cubic mismatch"));
  if(b.nf.disc!={prior['field_discriminant']},error("field discriminant mismatch"));
  print("PROVISIONAL_CYCLIC_FACTORS|",b.cyc);
  g=sum(i=1,#b.cyc,b.cyc[i]%2==0);if(g>16,error("NOT_SUFFICIENT"));
  print("STAGE|bnfcertify_flag1");gettime();C=bnfcertify(b,1);
  print("CLASSGROUP_QUOTIENT_CERTIFIED|",C);if(C!=1,error("quotient certification failed"));
  print("CERTIFY_MS|",gettime());print("G_UPPER|",g);print("DONE|quotient");
}};
iferr(main(),e,print("FAIL|",e);quit(1));quit(0)
'''
    for s, program in [('provisional', provisional), ('quotient', quotient)]:
        (OUT / (s+'.gp')).write_text(program)
    sources = [Path(__file__).resolve(), GP, NF, PRIOR / 'field.json', PRIOR / 'field.log', PRIOR / 'equation.json',
               CASE / 'terminal.json', CASE / 'verified.json', ROOT / 'elliptic-curves/cas/research_runtime/supervisor.py']
    checkpoint(OUT / 'protocol.json', {
        'schema': '11952.classgroup-quotient-lane.v1', 'host': socket.gethostname(),
        'available_memory_at_prepare': memory_available(),
        'input_sha256': {str(p): digest(p) for p in sources},
        'gp_version': subprocess.check_output([str(GP), '--version-short'], text=True).strip(),
        'program_sha256': {s: digest(OUT / (s+'.gp')) for s in ['provisional', 'quotient']},
        'wall_seconds': {'provisional': 600, 'quotient': 300},
        'pari_stack_bytes': STACK, 'rss_bytes': RSS, 'maximum_workers': 1,
        'known_points_passed_to_arithmetic': False, 'point_search': False,
        'bnf_flag': 0, 'bnfcertify_flag': 1, 'class_2rank_target': 16,
        'automatic_retry': False, 'automatic_budget_extension': False,
        'claim_gate': 'Completed sufficient provisional BNF gives only a GRH-conditional rank25 statement. Quotient certification ==1 plus g_comp<=16 plus the exact g+9 local bound and hash-matched independent rank>=25 replay gives unconditional rank25. Resource failures give no new bound. g_comp>16 stops before certification.'})
    print('PREPARED: provisional600s, quotient300s; PARI8GiB; RSS12GiB; one worker.')


def execute(stage):
    p = read(OUT / 'protocol.json')
    for name, expected in p['input_sha256'].items():
        assert digest(Path(name)) == expected, name
    assert digest(OUT / (stage+'.gp')) == p['program_sha256'][stage]
    if (OUT / (stage+'.json')).exists() or (OUT / (stage+'.supervisor.json')).exists():
        raise FileExistsError('one-shot stage; prior attempt preserved')
    provisional = None
    if stage == 'quotient':
        provisional = read(OUT / 'provisional.json')
        assert provisional['status'] == 'PROVISIONAL_GRH_SUFFICIENT'
        assert digest(LOCAL / 'provisional_bnf.bin') == provisional['bnf_sha256']
    assert memory_available() >= RSS + 2*2**30, 'insufficient current memory headroom'
    result = run([str(GP), '-q', '-f', '-s', str(STACK), str(OUT / (stage+'.gp'))],
                 limits=Limits(p['wall_seconds'][stage], RSS, pari_stack_bytes=STACK),
                 log_path=OUT / (stage+'.log'), checkpoint_path=OUT / (stage+'.supervisor.json'), cwd=ROOT)
    binary = LOCAL / 'provisional_bnf.bin'
    record = classify(stage, (OUT / (stage+'.log')).read_text(), result['outcome'], binary.exists(), provisional)
    record.update(stage=stage, supervision=result, protocol_sha256=digest(OUT / 'protocol.json'),
                  log_sha256=digest(OUT / (stage+'.log')), bnf_sha256=digest(binary) if binary.exists() else None)
    checkpoint(OUT / (stage+'.json'), record)
    print(json.dumps(record, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['prepare', 'provisional', 'quotient'])
    args = parser.parse_args()
    prepare() if args.stage == 'prepare' else execute(args.stage)
