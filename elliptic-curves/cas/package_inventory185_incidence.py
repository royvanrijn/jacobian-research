#!/usr/bin/env python3
"""Standalone all185 rational incidence and generic duplicate-subgroup proof."""
from hashlib import sha256
from pathlib import Path
import zipfile
import certify_compact_r17_candidates as cert
import certify_inventory185_incidence as cohort
from package_recorded_mod2_audit import dependencies
ROOT=cohort.ROOT;ART=cohort.ART;CAS=ROOT/'elliptic-curves/cas'
OUT=ART/'inventory185_incidence_evidence_v1.json'

def main():
    archive=OUT.with_suffix('.zip')
    if OUT.exists() or archive.exists():raise FileExistsError('preserve incidence evidence')
    proof=ART/'inventory185_incidence_v1.json'
    if cert.read(proof)!=cohort.result():raise ArithmeticError('exact185 aggregate differs')
    paths={Path(__file__).resolve(),CAS/'verify_inventory185_incidence.py',CAS/'replay_inventory185_all_incidence.py',CAS/'certify_inventory185_incidence.py',proof,cohort.INDEX,cohort.TRANSPORT,ROOT/'elliptic-curves/notes/INVENTORY185_INCIDENCE_2026-09-06.md',*cohort.PROOFS}
    for p in [proof,cohort.TRANSPORT,*cohort.PROOFS]:
        for name,h in cert.read(p)['sources'].items():
            q=ROOT/name
            if cert.hashed(q)!=h:raise ArithmeticError('frozen incidence input changed')
            paths.add(q)
    paths=dependencies(paths);members=[]
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for p in sorted(paths):
            raw=p.read_bytes();name=str(p.relative_to(ROOT));info=zipfile.ZipInfo(name,date_time=(2026,9,6,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o644<<16
            z.writestr(info,raw);members.append({'path':name,'sha256':sha256(raw).hexdigest(),'bytes':len(raw)})
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None or any(sha256(z.read(r['path'])).hexdigest()!=r['sha256'] for r in members):raise ArithmeticError('archive integrity differs')
    cert.write(OUT,{'schema':'elliptic-curves.inventory185-incidence-evidence.v1','archive':str(archive.relative_to(ROOT)),'archive_sha256':cert.hashed(archive),'archive_bytes':archive.stat().st_size,'files':members,'required_base_archives':[],
        'claim_boundary':'Standalone independent replay of all2220 rational j-incidence comparisons for185 equations across twelve family presentations, the nine-cohort aggregate and exact Q(t) model/unimodular17-section transport for all25 duplicated presentations. No point discovery or rank-proof replay; historical rank metadata only must not exceed current lower bounds. No unlisted family, nongeneric point, higher rank or universal novelty is excluded.'})
    print('ALL185 INCIDENCE BUNDLE',len(members),'files',archive.stat().st_size,'bytes',flush=True)
if __name__=='__main__':main()
