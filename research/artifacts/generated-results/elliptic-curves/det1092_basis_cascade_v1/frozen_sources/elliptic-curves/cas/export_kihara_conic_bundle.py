#!/usr/bin/env python3
"""Portable conic witnesses, with prior full-parent theorem as a dependency."""
import json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
files={'construction':'kihara_split_infinity_v1.json','control':'kihara_conic_control_v1.json','parent_rank':'kihara_first_parent_rank_v1.json','parent_replay':'kihara_first_parent_rank_replay_v1.json','parents':'kihara_five_parent_distinctness_v1.json'}
d={k:json.loads((ART/n).read_text()) for k,n in files.items()};d['parent']=d.pop('parents')['rows'][1];d['sources']={str((ART/n).relative_to(ROOT)):hashlib.sha256((ART/n).read_bytes()).hexdigest() for n in files.values()}
with (ART/'kihara_conic_replay_bundle_v1.json').open('x') as f:json.dump(d,f,indent=2);f.write('\n')
print('EXPORTED conic proof and control')
