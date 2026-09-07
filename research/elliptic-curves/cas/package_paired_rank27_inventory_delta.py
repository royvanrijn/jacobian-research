#!/usr/bin/env python3
"""Immutable delta retaining the timed-out inventory attempt and its replacement."""
from pathlib import Path
import zipfile
from hashlib import sha256
import certify_compact_r17_candidates as cert
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';OUT=ART/'paired_rank27_inventory_replay_evidence_v2.json'

def main():
    archive=OUT.with_suffix('.zip')
    if OUT.exists() or archive.exists():raise FileExistsError('preserve inventory replay delta')
    r=cert.read(LOCAL/'inventory-v4-memory-replay-v1/replay.supervisor.json')
    if r['outcome']!='completed' or r['returncode']!=0:raise ArithmeticError('replacement inventory replay incomplete')
    prior=cert.read(ART/'paired_rank27_portable_replay_v1.json')
    if len(prior['ledger']['rows'])!=35 or sum(r['status']=='PASS' for r in prior['ledger']['rows'])!=34:raise ArithmeticError('original isolated record differs')
    main_manifest=ART/'paired_rank27_discovery_evidence_v1.json';m=cert.read(main_manifest);bases=list(m['required_base_archives'])
    for path in (main_manifest,ART/'latest8_cross_family_evidence_v1.json'):
        b=cert.read(path);bases.append({'manifest':str(path.relative_to(ROOT)),'manifest_sha256':cert.hashed(path),'archive':b['archive'],'archive_sha256':b['archive_sha256']})
    for b in bases:
        if cert.hashed(ROOT/b['archive'])!=b['archive_sha256']:raise ArithmeticError('base archive changed')
    paths={Path(__file__).resolve(),ROOT/'elliptic-curves/cas/replay_inventory_v4_memory.py',ROOT/'elliptic-curves/cas/complete_paired_rank27_inventory_replay.py',ROOT/'elliptic-curves/notes/EXTENDED_PRIME_RANK27_DISCOVERIES_2026-09-06.md',ART/'paired_rank27_portable_replay_v1.json',ART/'new_high_rank_curve_index_v4_memory_replay_v1.json',ART/'latest8_incidence_portable_replay_v1.json'}
    for folder in ('inventory-v4-memory-replay-v1','paired-rank27-portable-v1/verification'):paths.update(p for p in (LOCAL/folder).rglob('*') if p.is_file())
    members=[]
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for p in sorted(paths):
            raw=p.read_bytes();name=str(p.relative_to(ROOT));info=zipfile.ZipInfo(name,date_time=(2026,9,6,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o644<<16;z.writestr(info,raw);members.append({'path':name,'bytes':len(raw),'sha256':sha256(raw).hexdigest()})
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None or any(sha256(z.read(r['path'])).hexdigest()!=r['sha256'] for r in members):raise ArithmeticError('delta integrity failure')
    cert.write(OUT,{'schema':'elliptic-curves.paired-rank27-inventory-replay-evidence.v2','archive':str(archive.relative_to(ROOT)),'archive_sha256':cert.hashed(archive),'archive_bytes':archive.stat().st_size,'required_base_archives':bases,'files':members,'claim_boundary':'Preserves34 passed original isolated checks and the inventory180-second timeout. Adds a separately versioned memory-cache inventory/CSV replayer and exact latest8 incidence evidence. Completion requires the new isolated inventory replay; no original outcome is rewritten.'});print('PACKAGED INVENTORY REPLAY DELTA',len(members),archive.stat().st_size,flush=True)
if __name__=='__main__':main()
