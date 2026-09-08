#!/usr/bin/env python3
"""Lossless evidence archive + compact algebraic certificates after all checks.

Expanded residual/discriminant polynomials are deterministic from RR lines and
retained as exact raw bytes in the ignored local evidence archive. All candidate
coordinates, all RR lines and all failures/censors remain in compact certificates.
"""
import copy,hashlib,json,os,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4];PKG=ROOT/'research/artifacts/generated-results/elliptic-curves/det1092_blind_mw16_v1'
LOCAL=ROOT/'research/artifacts/local/elliptic-curves/det1092_blind_mw16_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def put(p,d):
 with p.open('x') as f:json.dump(d,f,indent=2,sort_keys=True);f.write('\n')
 p.chmod(0o444)
replay=json.loads((PKG/'replay.json').read_text());assert replay['status']=='PASS_INDEPENDENT_GENERIC_REPLAY'
assert (PKG/'diagnostic.json').exists();commit=replay['starting_commit']
barrier=json.loads((PKG/'primary_complete.json').read_text())
for run in barrier['runs']:
 for p,h in run['files'].items():assert sha(PKG/p)==h
(PKG/'arms').mkdir();records=[]
for run in barrier['runs']:
 aid=run['arm_id'];d=json.loads((PKG/'runs'/aid/'result.json').read_text());d['original_result_sha256']=sha(PKG/'runs'/aid/'result.json')
 for stage in d.get('stages',[]):
  for member in stage.get('members',[]):
   removed={k:member.pop(k) for k in ['residual_coefficients','discriminant','square_root'] if k in member}
   if removed:
    member['derived_expansions_sha256']=hashlib.sha256(json.dumps(removed,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    member['derived_fields']=list(removed)
 d['serialization_note']='Only expanded residual coefficients, discriminant and square root replaced by a canonical JSON digest. Recompute exactly from the retained RR line, centre word and fixture. Original bytes including every field remain in the bound local raw archive.'
 target=PKG/'arms'/aid;target.mkdir();put(target/'certificate.json',d)
 for name in ['selection.json','accounting.json','stdout.txt']:
  src=PKG/'runs'/aid/name
  if src.exists():dest=target/name;dest.write_bytes(src.read_bytes());dest.chmod(0o444)
 events=[]
 for line in (PKG/'runs'/aid/'events.jsonl').read_text().splitlines():
  event=json.loads(line)
  if 'stage' in event:
   st=event.pop('stage');event['stage_status']=st['status'];event['centre']=st['centre'];event['n']=st['n'];event['recorded_members']=len(st.get('members',[]))
  if 'candidate' in event:event.pop('candidate')
  events.append(event)
 put(target/'stage_ledger.json',{'starting_commit':commit,'arm_id':aid,'raw_events_sha256':sha(PKG/'runs'/aid/'events.jsonl'),'events':events})
 records.append({'arm_id':aid,'compact_certificate':str((target/'certificate.json').relative_to(PKG)),'compact_sha256':sha(target/'certificate.json'),'raw_files':run['files']})
LOCAL.mkdir(parents=True,exist_ok=True);assert not (LOCAL/'runs').exists()
(PKG/'runs').rename(LOCAL/'runs')
put(PKG/'raw_archive.json',{'starting_commit':commit,'primary_barrier_sha256':sha(PKG/'primary_complete.json'),'raw_root_relative_to_repository':str(LOCAL.relative_to(ROOT)),'raw_bytes_preserved':True,'raw_files_moved_without_editing':True,'compact_reconstruction':'RR elimination from lines regenerates omitted expansions. Verify compact certificate digests with replay.py --compact. Full replay.py can still read the local originals.','arms':records,'packager_sha256':sha(Path(__file__))})
print('PACKAGED; raw bytes preserved under',LOCAL,flush=True)
