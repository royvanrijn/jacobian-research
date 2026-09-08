#!/usr/bin/env python3
"""Archive the immutable selector and verified nine-curve intake without live outcomes."""
import argparse,gzip
from pathlib import Path
import mixed_reduced_parent_exposure as batch
from research_runtime.store import checkpoint
cert=batch.cert;D=batch.BATCH;OUT=batch.ART/'mixed_reduced_parent_intake_v1'
def expected():
 p=batch.protocol();s=batch.L/'det1092-reduced-score-strata-v1'
 files={n:D/n for n in ['protocol.json','prepare-protocol.json','intake-input.json','intake-result.json','run-protocol.json']}
 files.update({'seed-replay-'+p.name:p for p in (D/'seed-replay').iterdir() if p.is_file()})
 files.update({'score-protocol.json':s/'protocol.json','score-result.json':s/'result.json'})
 for row in p['rows']:files[row['id']+'-seed.json']=D/row['id']/'seed.json'
 return files
def main(check):
 files=expected()
 if not check:
  assert not OUT.exists();OUT.mkdir()
  for n,p in files.items():
   data=p.read_bytes();dest=OUT/(n+'.gz' if n=='score-result.json' else n);dest.write_bytes(gzip.compress(data,mtime=0) if n=='score-result.json' else data)
  checkpoint(OUT/'manifest.json',dict(status='PASS_IMMUTABLE_INTAKE',files={n:cert.hashed(p) for n,p in files.items()},rows=len(batch.protocol()['rows']),maximum_point_boxes=441,point_outcomes='Not part of this immutable intake package; consult the live ledger or completed exposure certificate.',boundary='Frozen source-only score strata and unchanged retained26 source gaps. Exact specialized seed transports replay. This package proves no completed point-search outcome or new rank gain.'))
 manifest=cert.read(OUT/'manifest.json');assert manifest['files']=={n:cert.hashed(p) for n,p in files.items()}
 for n,p in files.items():
  raw=(OUT/(n+'.gz' if n=='score-result.json' else n)).read_bytes();assert (gzip.decompress(raw) if n=='score-result.json' else raw)==p.read_bytes()
 print('PASS immutable mixed-parent intake and selector package',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');main(p.parse_args().check)
