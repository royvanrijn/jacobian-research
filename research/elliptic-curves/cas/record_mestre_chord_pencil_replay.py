#!/usr/bin/env python3
"""Bind the completed standalone nine-pencil verification to portable files."""
import argparse,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/mestre-chord-pencils-standalone-v1';OUT=ART/'mestre_chord_pencils_portable_replay_v1.json'
def compute():
 r=json.loads((D/'supervisor.json').read_text());log=(D/'replay.log').read_bytes()
 assert r['outcome']=='completed' and r['returncode']==0 and r['failure_reason'] is None
 assert hashlib.sha256(log).hexdigest()==r['log_sha256']
 assert log.decode().splitlines()[0]=='PASS all781 words,151 eligible divisors and fixed eight-pencil admission'
 assert log.decode().splitlines()[-1]=='PASS9 INDEPENDENT BIRATIONAL MAPS AND EXACT GENERIC RANKS' and len(log.decode().splitlines())==11
 files={}
 for p in [ROOT/'elliptic-curves/cas/verify_mestre_chord_pencils.sage',ART/'mestre_chord_pencils_v1.json']:
  data=p.read_bytes();assert data==(D/p.name).read_bytes();files[str(p.relative_to(ROOT))]={'sha256':hashlib.sha256(data).hexdigest(),'bytes':len(data)}
 return {'schema':'elliptic-curves.mestre-chord-pencils-portable-replay.v1','status':'PASS','files':files,'supervisor':r,'transcript':log.decode(),'command':'sage -python verify_mestre_chord_pencils.sage --input mestre_chord_pencils_v1.json','scope':'Completed exact standalone replay, not a new point search, specialized-rank upper bound, or proof that these presentations have inferior discovery yield.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();r=compute()
 if a.check:assert r==json.loads(OUT.read_text())
 else:
  assert not OUT.exists();OUT.write_text(json.dumps(r,indent=2)+'\n')
 print('PASS nine-pencil portable bindings')
