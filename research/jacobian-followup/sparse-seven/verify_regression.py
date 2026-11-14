#!/usr/bin/env python3
"""Regression of seeded closure against the original unseeded enumeration."""
import itertools
import json
import tempfile
from pathlib import Path
import search

def main():
    degree,bound=3,7
    compiler=search.SupportEnumerator(degree,bound)
    original,stats=compiler.enumerate()
    expected={s for row in original.values() for s in row
              if s.collision_axis_unimodularity_possible()}
    actual=set()
    with tempfile.TemporaryDirectory() as directory:
        for p in range(2,degree+1):
            shard=search.Shard(degree,bound,p,2,3,compiler)
            state=shard.run(Path(directory)/f"{p}.json",seconds=60,node_limit=500000)
            assert state["status"]=="complete"
            actual.update(search.read_support(rows) for rows in state["supports"])
    # Orientation picks F2 to carry >=2 axis terms. Include its residual swap
    # before comparing to the unseeded labelled enumeration.
    actual|={s.swapped_23() for s in list(actual)}
    assert actual==expected,(len(actual),len(expected),len(actual^expected))
    # Validate resumption across the node limit instead of relying on a
    # checkpoint's existence as evidence of completion.
    with tempfile.TemporaryDirectory() as directory:
        path=Path(directory)/"resume.json"
        shard=search.Shard(3,7,2,2,3,compiler)
        state=shard.run(path,60,1)
        assert state["status"]=="incomplete" and state["stack"]
        state=shard.run(path,60,500000)
        assert state["status"]=="complete" and not state["stack"]
        resumed={search.read_support(rows) for rows in state["supports"]}
        direct={s for s in expected if s.rows[0] and
                min(e[0] for e in s.rows[0] if e[1:]==(0,0))==2 and
                s.collision_axis_counts()[1]>=2}
        assert resumed==direct
    receipt=dict(status="pass",degree=degree,bound=bound,
                 original_balanced=sum(map(len,original.values())),
                 axis_admissible_labelled=len(expected),labelled_matched=len(actual),
                 original_statistics=stats,resume_verified=True)
    search.write_json(search.HERE/"regression.json",receipt)
    print(json.dumps(receipt))
if __name__=="__main__":
    main()

