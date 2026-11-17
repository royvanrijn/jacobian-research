#!/usr/bin/env python3
"""Seal the read source scope without rewriting the older deposited manifest."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
COMMIT = '06cbbb7e2737e1a11085b64e0ead021d63c550c3'


def main():
    files = []
    for name in ['glsearch.py', 'w3.py']:
        path = HERE/'sources/current13/search'/name
        raw = path.read_bytes()
        files.append({'path':str(path.relative_to(HERE)),
                      'url':f'https://github.com/anaxia-research/conjecture-to-choose/blob/{COMMIT}/search/{name}',
                      'sha256':hashlib.sha256(raw).hexdigest(),
                      'computed_git_blob_sha1':hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest(),
                      'provenance':'Pinned GitHub connector raw content; blob hash computed locally, not supplied by connector'})
    code = (HERE/'sources/current13/search/glsearch.py').read_text()
    assert 'seed vectors e_1..e_4' in code and 'keep=det4(A,q)!=0' in code
    random_code = (HERE/'sources/current13/search/w3.py').read_text()
    assert '.backtrack(random.randint(2,40))' in random_code
    report = {'status':'PASS source seals and explicit search guards','commit':COMMIT,'files':files,
              'coverage_finding':'Four seed image directions are assumed independent, and singular seed pairings discarded. Homogeneous affine3 extensions have rank at most3 for every four-direction pairing and lie outside this search even when the common fixed space is zero.',
              'certificate_validity':'No challenge to any243 replayed deposited certificates or the rank-one fixed-space criterion.',
              'three_seed_boundary':'Bounded random diagram moves, not an exhaustive three-bridge exclusion.',
              'ATLAS_import':'Primary L32 page lists natural3,dual3,Steinberg8 in characteristic2. Full matrix-algebra spans and four2-regular classes are independently checked; Brauer class-count theorem imported.',
              'assurance':'Source content inspected and mathematical scope argued in the canonical proof; text guards are not a substitute for that argument.',
              'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in ['audit_gl3_sources.py','sources/atlas-l32.html','sources/current13/manifest.json']}}
    (HERE/'gl3_source_scope_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'source_files':len(files)}))


if __name__=='__main__':main()
