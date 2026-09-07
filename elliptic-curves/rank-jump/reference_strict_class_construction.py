#!/usr/bin/env python3
"""Checkpoint the equation-only reference class/unit construction attempt."""
import argparse
from pathlib import Path
import subprocess
import sys
import time
import retrospective as r

PROTOCOL=Path(__file__).with_name('REFERENCE_STRICT_CLASS_CONSTRUCTION_PROTOCOL.json')
SOURCE=r.OUT/'rank_jump_prime_square_class_constructor_inputs_v1.json'
INPUT=r.OUT/'rank_jump_reference_strict_class_construction_inputs_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-reference-strict-class-construction-v1'


def export():
    source=r.read(SOURCE)
    data={k:source[k] for k in ['cubic_ascending','field_discriminant','S_finite','generic_classes']}
    assert len(data['generic_classes'])==16
    data['bindings']={str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in [SOURCE,PROTOCOL,Path(__file__)]}
    r.write_new(INPUT,data)


def worker():
    from sage.all import QQ,PolynomialRing,pari
    p=r.read(PROTOCOL);c=r.read(INPUT)
    pari.allocatemem(p['bounds']['initial_stack_bytes'],p['bounds']['pari_stack_bytes'],silent=True)
    pari.setrand(p['bounds']['random_seed']);pari.default('debug',1)
    R=PolynomialRing(QQ,'z');f=R(c['cubic_ascending'])
    nf=pari.nfinit([pari(f),c['S_finite']])
    assert str(nf.disc())==c['field_discriminant']
    r.write_new(WORK/'setup.json',{'status':'PASS','field_discriminant':str(nf.disc()),
        'basis':[str(x) for x in nf.nf_get_zk()],'signature':list(map(int,nf.nf_get_sign())),
        'pari_version':str(pari.version())})
    print('SETUP_COMPLETE',flush=True)
    try:
        bnf=pari.bnfinit(nf,1)
        r.write_new(WORK/'bnf.json',{'status':'RETURNED_UNCERTIFIED','compact_bnf':str(bnf),
            'cyclic_factors':list(map(str,bnf.bnf_get_cyc()))})
        result={'status':'BNF_RETURNED_UNCERTIFIED','explicit_new_strict_class':'NOT_YET_CONSTRUCTED'}
    except Exception as exc:
        result={'status':'UNKNOWN','exception_type':type(exc).__name__,'reason':str(exc)}
    r.write_new(WORK/'terminal.json',result);print(result['status'],flush=True)


def capture():
    WORK.mkdir(parents=True,exist_ok=True)
    assert not (WORK/'terminal.json').exists() and not (WORK/'worker.log').exists()
    start=time.monotonic()
    with (WORK/'worker.log').open('x') as stream:
        try:
            proc=subprocess.run([sys.executable,__file__,'worker'],stdout=stream,stderr=stream,
                                timeout=r.read(PROTOCOL)['bounds']['reference_bnf_seconds'])
            reason=None if proc.returncode==0 else 'worker failure'
        except subprocess.TimeoutExpired:reason='300-second timeout'
    if reason and not (WORK/'terminal.json').exists():r.write_new(WORK/'terminal.json',{'status':'UNKNOWN','reason':reason})
    r.write_new(WORK/'execution.json',{'elapsed_seconds':time.monotonic()-start,'process_error':reason})
    print(r.read(WORK/'terminal.json'),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','worker','capture']);args=p.parse_args()
    globals()[args.mode]()
