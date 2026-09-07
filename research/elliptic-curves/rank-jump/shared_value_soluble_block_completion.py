#!/usr/bin/env python3
"""Preserve the collapsed anchor; one declared replacement closes the S3 rank proof."""
import argparse
from pathlib import Path
from types import SimpleNamespace
import retrospective as r
import shared_value_soluble_block as first
import shared_value_halving as half

PROTOCOL=Path(__file__).with_name('SHARED_VALUE_COMPLETION_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_shared_value_soluble_block_completion_v1.json'


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
        (Path(__file__),PROTOCOL,Path(first.__file__),first.PROTOCOL,Path(half.__file__),half.OUTPUT)}


def compute():
    from sage.all import QQ,PolynomialRing,EllipticCurve
    from sage.version import version
    spec=r.read(PROTOCOL);effective=r.read(first.PROTOCOL)
    effective['specializations']['independence_anchor_n']=spec['replacement_anchor']['n']
    saved=first.r
    proxy=SimpleNamespace(**vars(r))
    proxy.read=lambda p:effective if Path(p)==first.PROTOCOL else r.read(p)
    try:
        first.r=proxy
        replacement=first.compute(spec['replacement_anchor']['c'])
    finally:first.r=saved
    replacement['frozen_producer_bindings']=replacement.pop('bindings')
    replacement['effective_protocol']=effective
    original=r.read(first.WORK/'c-1.json');assert original['bindings']==first.bindings()
    # This relation is an identity along a whole parameter locus, not a guessed
    # explanation of a failed fingerprint at one specialization.
    R=PolynomialRing(QQ,'c');c=R.gen();n=c+2;m=(n*n-c*c-2)/2
    E=EllipticCurve(R.fraction_field(),[0,m,0,-m-3,c*c])
    P0=E(0,c);Pm=E(-1,n);Pp=E(2,n)
    assert P0+Pm+2*Pp==E(0)
    halves=r.read(half.OUTPUT)
    assert halves['subset_divisions'][2]['mask']==3
    assert halves['subset_divisions'][2]['halves']==[['2','-5','1']]
    assert halves['original_fingerprint_rank']==2
    return {'schema':'rank-jump.shared-value-soluble-block-completion.v1','bindings':bindings(),
        'software':{'sage':version},'rows':[original,replacement],
        'collapsed_S3_anchor':{'c':3,'n':5,'m':7,'retained_subgroup_rank':2,
            'relation':'R+P_-1+2P_2=0','lower_bound':'two independent Kummer fingerprints in the retained exact halving certificate'},
        'symbolic_collapse_locus':{'equation':'n=c+2','parameter_m':'2*c+1','relation':'R+P_-1+2P_2=0'},
        'initial_failure_preserved':{'path':str((first.WORK/'c-3.log').relative_to(r.ROOT)),
            'sha256':r.digest((first.WORK/'c-3.log').read_bytes()),
            'reason':'The original three-point fingerprint rank is two because of an actual dependence.'}}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args()
    data=compute()
    if args.mode=='check':assert r.read(OUTPUT)==data;print('PASS completed simultaneous-solubility mechanism')
    else:r.write_new(OUTPUT,data);print('PASS two generic rank-three families and both collapse controls')
