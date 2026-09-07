#!/usr/bin/env python3
"""Replay fixed ideal-reduction identities and actual repeated-ideal units."""
import argparse
from pathlib import Path
import retrospective as r
import curve302_arithmetic_unit_class as source
OUTPUT=r.OUT/'rank_jump_curve302_arithmetic_unit_class_verification_v1.json'

def compute():
    from sage.all import QQ,PolynomialRing,pari
    d=r.read(source.OUTPUT);ar=r.read(source.ARITH)
    for path,sha in d['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    R=PolynomialRing(QQ,'z');nf=pari.nfinit([pari(R(ar['cubic_ascending'])),ar['S_finite']])
    checked=0;repeats=0
    for item in d['trials']:
        path=r.ROOT/item['path'];assert r.digest(path.read_bytes())==item['sha256'];row=r.read(path)
        I=pari(row['starting_ideal_hnf']);seen={}
        assert len(row['trials'])==81
        for trial in row['trials']:
            J=pari(trial['reduced_ideal_hnf']);a=pari(trial['multiplier_GP'])
            assert pari.idealmul(nf,J,a)==I;checked+=1;key=str(J)
            if key in seen:
                u=pari.nfbasistoalg(nf,a)/pari.nfbasistoalg(nf,seen[key])
                assert u in [1,-1];repeats+=1
            else:seen[key]=a
    assert checked==567 and repeats==81 and not d['retained_units'] and not d['strict_candidates']
    return {'schema':'rank-jump.curve302-arithmetic-unit-class-verification.v1','status':'PASS',
        'exact_ideal_identities':checked,'repeated_ideal_pairs':repeats,'nontrivial_norm_one_units':0,
        'bindings':source.r.read(source.OUTPUT)['bindings']|{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in [Path(__file__),source.OUTPUT]},
        'boundary':'All repeated ideals in this fixed run give only units+/-1. This is not an exclusion of other unit algorithms or other302 classes.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();result=compute()
    if a.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS567 exact ideal identities; repeated-ideal units trivial')
