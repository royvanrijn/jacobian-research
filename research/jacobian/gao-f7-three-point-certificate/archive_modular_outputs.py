"""Lossless compression of owned modular logs, retaining hashes and mapping."""
from pathlib import Path
import gzip,hashlib,json,os
HERE=Path(__file__).resolve().parent
manifest_path=HERE/'modular-output-archive.json'
rows=json.loads(manifest_path.read_text())['files'] if manifest_path.exists() else []
for source in sorted(HERE.glob('cubic-six-dag-*.txt')):
 assert source.resolve().parent==HERE.resolve()
 data=source.read_bytes();digest=hashlib.sha256(data).hexdigest();target=source.with_suffix('.txt.gz');temporary=target.with_suffix('.gz.research/.scratch')
 with temporary.open('wb') as raw:
  with gzip.GzipFile(filename='',mode='wb',fileobj=raw,compresslevel=1,mtime=0) as zipped:zipped.write(data)
 assert hashlib.sha256(gzip.open(temporary,'rb').read()).hexdigest()==digest
 os.replace(temporary,target)
 rows.append({'original':source.name,'archive':target.name,'uncompressed_sha256':digest,'archive_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'uncompressed_bytes':len(data),'archive_bytes':target.stat().st_size})
 # The byte-for-byte checked gzip retains the entire original evidence.
 source.unlink()
 manifest_path.write_text(json.dumps({'format':'gzip; original bytes checked before removing the uncompressed duplicate','files':rows},indent=2)+'\n')
print(json.dumps({'files':len(rows),'uncompressed_MB':sum(x['uncompressed_bytes'] for x in rows)/1e6,'archived_MB':sum(x['archive_bytes'] for x in rows)/1e6}))
