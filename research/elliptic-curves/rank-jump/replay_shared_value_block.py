#!/usr/bin/env python3
"""Portable replay of the frozen mechanism and its initial failure evidence."""
import argparse
from pathlib import Path
from unittest.mock import patch
import retrospective as r
import shared_value_soluble_block as first
import shared_value_soluble_block_completion as completed

INPUT=r.OUT/'rank_jump_shared_value_initial_inputs_v1.json'


def capture():
    paths=[first.WORK/'c-1.json',first.WORK/'c-3.log']
    r.write_new(INPUT,{'schema':'rank-jump.shared-value-initial-inputs.v1',
        'bindings':first.bindings(),'files':[{'path':str(p.relative_to(r.ROOT)),
            'sha256':r.digest(p.read_bytes()),'text':p.read_text()} for p in paths]})


def replay():
    data=r.read(INPUT);assert data['bindings']==first.bindings()
    files={r.ROOT/row['path']:row['text'].encode() for row in data['files']}
    for row in data['files']:assert r.digest(row['text'].encode())==row['sha256']
    original=Path.read_bytes
    def read(path):return files[path] if path in files else original(path)
    with patch.object(Path,'read_bytes',read):
        assert first.compute(1)==r.read(first.WORK/'c-1.json')
        assert completed.compute()==r.read(completed.OUTPUT)
    print('PASS portable symbolic, arithmetic and collapse-control replay')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','check']);args=p.parse_args()
    capture() if args.mode=='capture' else replay()
