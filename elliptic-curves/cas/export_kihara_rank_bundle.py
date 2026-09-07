#!/usr/bin/env python3
"""Portable inputs for independent first-parent Picard and saturation replay."""
import json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
files={'height':ART/'kihara_section_involution_v1.json','rank':ART/'kihara_first_parent_rank_v1.json','counts':ROOT/'artifacts/local/elliptic-curves/kihara-picard-count-v1/counts.json','odd':ART/'kihara_fresh_point_pilot_fibre0_modl_v1.json','parents':ART/'kihara_five_parent_distinctness_v1.json','seeds':ART/'kihara_parent_replay_bundle_v1.json'}
d={k:json.loads(p.read_text()) for k,p in files.items()};d['parent']=d.pop('parents')['rows'][1];d['seed']=d.pop('seeds')['rows'][1];d['sources']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files.values()}
out=ART/'kihara_first_parent_rank_bundle_v1.json'
with out.open('x') as f:json.dump(d,f,indent=2);f.write('\n')
print('EXPORTED first-parent rank witnesses')
