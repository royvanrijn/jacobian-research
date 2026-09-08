#!/usr/bin/env python3
"""Preserve the five pinned population inputs omitted from the initial archive."""
from pathlib import Path
import zipfile
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';OUT=ART/'million_height_mw16_population_supplement_v1.json'
def main():
    archive=OUT.with_suffix('.zip')
    if OUT.exists() or archive.exists():raise FileExistsError('preserve population supplement')
    p=cert.read(LOCAL/'mw16-extended-prime-benchmark-v1/protocol.json');paths=[Path(__file__).resolve(),ROOT/'elliptic-curves/cas/verify_million_height_mw16_bundle_v2.py']
    for row in p['rows']:
        path=LOCAL/'prospective-mw16-h4096-v1'/row['family']/'population.json'
        if cert.hashed(path)!=row['population_sha256'] or cert.read(path)['retained_candidates'][4]['parameter']!=row['parameter']:raise ArithmeticError('frozen benchmark population differs')
        paths.append(path)
    members=[]
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for path in paths:
            name=str(path.relative_to(ROOT));raw=path.read_bytes();info=zipfile.ZipInfo(name,date_time=(2026,9,6,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o644<<16;z.writestr(info,raw);members.append({'path':name,'sha256':cert.hashed(path),'bytes':len(raw)})
    checkpoint(OUT,{'schema':'elliptic-curves.million-mw16-population-supplement.v1','archive':str(archive.relative_to(ROOT)),'archive_sha256':cert.hashed(archive),'archive_bytes':archive.stat().st_size,'base_manifest_sha256':cert.hashed(ART/'million_height_mw16_evidence_v1.json'),'files':members,'reason':'The first isolated trace benchmark replay exposed five omitted population.json inputs. Their hashes had been pinned and checked locally, but the initial packager did not discover the implicit paths. Include the exact five original inputs and a version2 extractor. Original archive, failed replay and frozen arithmetic sources remain unchanged. No new search, trace, candidate selection or curve proof.'})
    print('PACKAGED FIVE PINNED MW16 POPULATIONS',archive.stat().st_size,flush=True)
if __name__=='__main__':main()
