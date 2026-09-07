#!/usr/bin/env python3
"""Export retained Kihara witnesses for a repository-independent Sage replay."""
import json,sys,hashlib
from pathlib import Path
from fractions import Fraction
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';sys.path.insert(0,str(ROOT/'elliptic-curves'))
from ecsearch.kihara import kihara_rank14_replay
files=['kihara_global_parent_geometry_v1.json','kihara_fixed_parent_geometry_v1.json','kihara_five_parent_distinctness_v1.json','curve302_inverse_kihara_and_rank14_16_intake_v1.json','mestre_parent_portfolio_intake_v1.json','kihara_rank14_t2_v1.json']
paths=[ART/n for n in files];data={n:json.loads((ART/n).read_text()) for n in files};rows=[]
for i,t in enumerate(['2','3/2','5/2','11/3','1009/101']):
 c=kihara_rank14_replay(Fraction(t));row={'parameter':t,'old_model':list(map(str,c.weierstrass_coefficients)),'old_points':[[str(x) for x in p] for p in c.weierstrass_points]}
 if i:
  for key,path in [('seed',LOCAL/'kihara-fresh-fibres-v2'/f'fibre{i-1}'/'seed.json'),('old_intake',LOCAL/'kihara-fresh-fibres-v1'/f'fibre{i-1}'/'seed.json'),('finite',ART/f'kihara_fresh_point_pilot_fibre{i-1}_mod2_v1.json')]:
   row[key]=json.loads(path.read_text());paths.append(path)
 rows.append(row)
out=ART/'kihara_parent_replay_bundle_v1.json';assert not out.exists()
out.write_text(json.dumps({'schema':'elliptic-curves.kihara-parent-replay-bundle.v1','certificates':data,'rows':rows,'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}},indent=2)+'\n')
print('EXPORTED five parents and four corrected fibre seeds')
