#!/usr/bin/env python3
"""One documented memory correction; preserve both original seeded attempts."""
from pathlib import Path
import retrospective as r
import seeded_reference_relations as run
import seeded_reference_class as previous

if __name__=='__main__':
    run.WORK=r.ROOT/'artifacts/local/rank-jump-seeded-reference-relations-memory-v1'
    run.PROTOCOL=Path(__file__).with_name('SEEDED_REFERENCE_RELATIONS_MEMORY_PROTOCOL.json')
    run.OUTPUT=r.OUT/'rank_jump_seeded_reference_relations_memory_v1.json'
    assert previous.MAIN.count('268435456')==1
    previous.MAIN=previous.MAIN.replace('268435456','1073741824')
    run.WORK.mkdir(parents=True,exist_ok=True)
    r.write_new(run.WORK/'launch.json',{'bindings':previous.bindings([Path(__file__),run.PROTOCOL,
        Path(run.__file__),Path(previous.__file__),previous.INPUT]),
        'memory_only_substitution':['268435456','1073741824']})
    run.capture()
