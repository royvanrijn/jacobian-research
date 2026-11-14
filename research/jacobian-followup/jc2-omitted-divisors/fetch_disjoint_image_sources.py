"""Explicit network refresh of the versioned source audit; not a replay step."""
from pathlib import Path
import hashlib
import json
import urllib.request

HERE=Path(__file__).resolve().parent/'sources'
rows=[]
for name,url,scope in [
    ('jelonek-2011.03472v3.html','https://arxiv.org/html/2011.03472v3',
     'Theorem1.2 retains smoothness; proof of Theorem1.1 contains the elementary disjoint-curve constant restriction argument. No connected singular-set Euler conclusion imported.'),
    ('jelonek-2011.03472-version-audit.html','https://arxiv.org/abs/2011.03472',
     'Version history: author rejects v4/v5; latest withdrawal is not used as a theorem. v3 is the explicitly selected text.')]:
    request=urllib.request.Request(url,headers={'User-Agent':'jacobian-research source custody'})
    data=urllib.request.urlopen(request,timeout=30).read()
    (HERE/name).write_bytes(data)
    rows.append(dict(path=name,url=url,bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),scope=scope))
manifest=dict(as_of='2026-11-14',sources=rows,
    boundary='Version/scope audit only. The disjoint-curve lemma, pole-module argument and critical-value obstruction are proved in DISJOINT_IMAGE_GATE.md. No nonsmooth connected Euler theorem is assumed.')
(HERE/'disjoint_image_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(rows))
