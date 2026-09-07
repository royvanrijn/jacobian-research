#!/usr/bin/env python3
"""Read-only replay of the frozen capacity object through its JSON boundary.
The producer's direct check compares Python tuples with loaded JSON lists.
Normalize only JSON representation; never change values or frozen sources.
"""
import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'rank-jump'))
import native_common_class_capacity as capacity
actual=json.loads(json.dumps(capacity.compute()))
assert actual==capacity.r.read(capacity.OUTPUT)
print('PASS native common-class capacity after canonical JSON tuple/list normalization')
