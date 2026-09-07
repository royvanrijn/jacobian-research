#!/usr/bin/env python3
"""Replay the sealed producer with canonical JSON comparison of basis keys."""
import json
import retrospective as r
import small_jacobian_selmer as ex


def replay():
    retained=r.read(ex.OUTPUT);write=r.write_new;checked=[]
    def compare(path,value):
        assert path==ex.OUTPUT
        # JSON object keys are strings, including integer pivot indices.
        assert json.loads(json.dumps(value))==retained
        checked.append(True)
    try:
        r.write_new=compare
        ex.build(check=False)
    finally:
        r.write_new=write
    assert checked==[True]
    print("PASS sealed producer replay with canonical JSON basis indices")


if __name__=="__main__":
    replay()
