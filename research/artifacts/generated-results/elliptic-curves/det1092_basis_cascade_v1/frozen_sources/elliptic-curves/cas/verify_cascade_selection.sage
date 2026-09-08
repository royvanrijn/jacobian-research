#!/usr/bin/env sage-python
"""Recompute every frozen selection using only generic inputs and own seeds."""
from pathlib import Path
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
CAS=Path(__file__).resolve().parent
b=SourceFileLoader('selector_replay_base',str(CAS/'adaptive_visibility_cascade.sage')).load_module()
from research_runtime.store import checkpoint

def main():
    b.guard();b.protocol()
    try:
        (b.ART/'curve302_residual_visibility_geometry_v1.json').read_text()
    except PermissionError:pass
    else:raise ArithmeticError('target artifact guard failed')
    checked=0
    for row in b.read(b.D/'roster.json'):
        model,initial=b.seed(row)
        if [list(map(str,p)) for p in initial]!=b.read(b.D/row['id']/'seed.json')['points']:raise ArithmeticError('generic seed changed')
        terminal=b.read(b.D/'continuation-v1'/row['id']/'terminal.json');seen=set()
        for s in terminal['stages']:
            wd=(b.ROOT/s['mod2']).parent;selection=b.read(wd/'selection.json')
            basis=tuple(tuple(map(F,p)) for p in selection['basis'])
            centres,metric=b.candidates(model,basis,seen)
            if centres!=selection['centres'] or metric!=selection['metric']:raise ArithmeticError('target-free selector did not reproduce')
            for c in centres:seen.add((c['lane'],c['orbit'] if c['lane']=='canonical' else c['parity']))
            checked+=len(centres)
        print('SELECTION VERIFIED',row['id'],len(terminal['stages']),'waves',flush=True)
    checkpoint(b.D/'selection-replay.json',{'status':'PASS_GENERIC_ONLY_SELECTION_REPLAY','charts':checked,'forbidden_target_artifact_read_rejected':True,'source_sha256':b.sha(Path(__file__)),'read_paths':sorted(b.READS)})

if __name__=='__main__':main()
