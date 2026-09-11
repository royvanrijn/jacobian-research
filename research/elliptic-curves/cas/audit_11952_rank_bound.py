#!/usr/bin/env python3
"""One equation-only, bounded cubic-class-group probe for 11952 at 921/653.

No point coordinates enter either arithmetic worker. Preserve unsuccessful
attempts; a timeout or an uncertified BNF never becomes an upper bound.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess

from research_runtime.store import checkpoint
from research_runtime.supervisor import Limits, run

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'artifacts/generated-results/elliptic-curves/11952_rank_bound_audit_v1'
LOCAL = ROOT / 'artifacts/local/elliptic-curves/11952-rank-bound-audit-v1'
SOURCE = ROOT / 'artifacts/local/elliptic-curves/broad-rank-v1/runtime/research/broad-cases/b-7bb187bc9254e81c6283/batch-002/search-00/terminal.json'
GP = Path('/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/gp')
PRIMES = [2, 3, 5, 7, 13, 17, 19, 47, 83, 7568348148323,
          1023034923031265640083839109496335411267837093185764838343605884505996669629389738604483673]


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def read(p):
    return json.loads(p.read_text())


def prepare():
    if OUT.exists() or LOCAL.exists():
        raise FileExistsError('preserve the frozen attempt')
    curve = read(SOURCE)['curve']
    assert curve == ['0', '0', '0', '-84196621147240320617047979399750704971769668763392',
                     '206937637967498573533732205060490773473330035432339348502618816917963546624']
    OUT.mkdir(parents=True)
    LOCAL.mkdir(parents=True)
    checkpoint(OUT / 'equation.json', {'fibre': '11952', 't': '921/653', 'a_invariants': curve,
                                     'source_sha256': digest(SOURCE), 'point_data_imported': False})
    common = ('default(realprecision,80);setrand(1);\n'
              f'P={PRIMES};for(i=1,#P,if(!isprime(P[i]),error("nonprime hint")));addprimes(P);\n')
    field = common + f'''
E=ellinit([{','.join(curve)}]);v=0;M=ellminimalmodel(E,&v);
print("MINIMAL_MODEL|",vector(5,i,M[i]));print("CHANGE_TO_MINIMAL|",v);
print("MINIMAL_DISCRIMINANT|",M.disc);print("C4|",M.c4);
D=factor(M.disc);if(prod(i=1,matsize(D)[1],D[i,1]^D[i,2])!=M.disc,error("factorization"));
for(i=1,matsize(D)[1],if(!isprime(D[i,1]),error("nonprime factor"));print("DISC_FACTOR|",[D[i,1],D[i,2]]));
if(M.a1!=0 || M.a3!=0,error("unsupported model"));
f=x^3+M.a2*x^2+M.a4*x+M.a6;
print("CUBIC_COEFFICIENTS_ASCENDING|",vector(4,i,polcoef(f,i-1)));
print("IRREDUCIBLE_MOD_23|",polisirreducible(Mod(1,23)*f));
print("STAGE|nfinit");gettime();nf=nfinit([f,P]);
print("MS_NFINIT|",gettime());print("CERTIFY_ORDER|",nfcertify(nf));
print("FIELD_DISCRIMINANT|",nf.disc);print("FIELD_SIGNATURE|",nf.sign);
print("FIELD_INDEX|",nf.index);print("FIELD_BASIS|",nf.zk);
n=0;
for(i=1,matsize(D)[1],p=D[i,1];L=elllocalred(M,p);Q=idealprimedec(nf,p);mult=(L[1]==1);add=(L[1]>1);contribution=if(mult,D[i,2]%2==0,if(add,#Q-1,0));n+=contribution;print("LOCAL|",[p,D[i,2],L[1],L[2],#Q,contribution]);print("SPLIT|",p,"|",vector(#Q,j,[Q[j].e,Q[j].f])));
u=if(M.disc>0,2,1);print("BOUND_OFFSET_U_PLUS_N|",u+n);
print("ROOT_NUMBER|",ellrootno(M));
writebin("{LOCAL / 'nf.bin'}",nf);print("DONE|field");quit
'''
    bnf = common + f'''
nf=read("{LOCAL / 'nf.bin'}");setdebug("bnf",3);gettime();
print("STAGE|bnfinit");b=bnfinit(nf,0);print("MS_BNFINIT|",gettime());
print("COMPUTED_CYCLIC_FACTORS|",b.cyc);writebin("{LOCAL / 'uncertified_bnf.bin'}",b);
print("STAGE|bnfcertify_quotient");C=bnfcertify(b,1);
print("CLASS_QUOTIENT_CERTIFIED|",C);print("MS_CERTIFY|",gettime());
print("CLASS_2RANK_UPPER|",sum(i=1,#b.cyc,b.cyc[i]%2==0));
print("DONE|class_group");quit
'''
    for name, program in [('field', field), ('class_group', bnf)]:
        (OUT / (name + '.gp')).write_text(program)
    paths = [Path(__file__).resolve(), GP, OUT / 'equation.json',
             ROOT / 'elliptic-curves/cas/research_runtime/supervisor.py']
    checkpoint(OUT / 'protocol.json', {
        'schema': '11952.equation-only-rank-bound.v1',
        'inputs': {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p): digest(p) for p in paths},
        'gp_version': subprocess.check_output([str(GP), '--version-short'], text=True).strip(),
        'programs': {s: digest(OUT / (s + '.gp')) for s in ['field', 'class_group']},
        'wall_seconds': {'field': 30, 'class_group': 60}, 'maximum_workers': 1,
        'rss_bytes': 1610612736, 'pari_stack_bytes': 512000000,
        'point_search': False, 'automatic_retry': False,
        'arithmetic_gate': 'Klagsbrun-Sherman-Weigandt Proposition 3.1: rank <= dim Sel2 <= g+u+n, for no rational 2-torsion. Certify maximal order and local offset, then one class-quotient canary.',
        'claim_boundary': 'No known point/Kummer generator is supplied to the measured class-group object. The existing 25-point certificate is used only after an independent upper bound. BNF invariants alone are GRH-conditional; flag-1 bnfcertify success certifies a quotient and hence an unconditional 2-rank upper bound. No algebraic parity assumption.'})
    print('Prepared one 30-second field diagnostic and one 60-second class-group probe.')


def execute(stage):
    p = read(OUT / 'protocol.json')
    for name, value in p['inputs'].items():
        assert digest(ROOT / name) == value, name
    assert digest(OUT / (stage + '.gp')) == p['programs'][stage]
    if (OUT / (stage + '.json')).exists():
        raise FileExistsError('no automatic repeat of a bounded attempt')
    if stage == 'class_group':
        prior = read(OUT / 'field.json')
        assert prior['valid_completion'] and prior['maximal_order_certified']
    result = run([str(GP), '-q', '-f', '-s', str(p['pari_stack_bytes']), str(OUT / (stage + '.gp'))],
                 limits=Limits(p['wall_seconds'][stage], p['rss_bytes']),
                 log_path=OUT / (stage + '.log'), checkpoint_path=OUT / (stage + '.supervisor.json'), cwd=ROOT)
    log = (OUT / (stage + '.log')).read_text()
    valid = result['outcome'] == 'completed' and f'DONE|{stage}' in log and '***' not in log
    record = {'stage': stage, 'supervision': result, 'valid_completion': valid,
              'program_sha256': p['programs'][stage], 'log_sha256': digest(OUT / (stage + '.log')),
              'unconditional_rank_upper_bound': None}
    if stage == 'field':
        record['maximal_order_certified'] = valid and 'CERTIFY_ORDER|[]' in log
        if record['maximal_order_certified']:
            fields = dict(line.split('|', 1) for line in log.splitlines() if '|' in line)
            assert fields['IRREDUCIBLE_MOD_23'] == '1'
            offset = int(fields['BOUND_OFFSET_U_PLUS_N'])
            record.update(bound_offset=offset, sufficient_class_2rank_upper_for_rank25=25-offset,
                          field_discriminant=fields['FIELD_DISCRIMINANT'], rational_2torsion=False)
    else:
        certified = valid and 'CLASS_QUOTIENT_CERTIFIED|1' in log
        record['class_quotient_certified'] = certified
        if certified:
            fields = dict(line.split('|', 1) for line in log.splitlines() if '|' in line)
            record['class_2rank_upper'] = int(fields['CLASS_2RANK_UPPER'])
            record['unconditional_rank_upper_bound'] = record['class_2rank_upper'] + read(OUT / 'field.json')['bound_offset']
    checkpoint(OUT / (stage + '.json'), record)
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['prepare', 'field', 'class_group'])
    args = parser.parse_args()
    prepare() if args.stage == 'prepare' else execute(args.stage)
