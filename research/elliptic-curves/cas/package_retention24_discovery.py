#!/usr/bin/env python3
"""Portable wider-retention discovery and exact control-diagnostic evidence."""
from pathlib import Path
from hashlib import sha256
import zipfile
import certify_compact_r17_candidates as cert
from package_recorded_mod2_audit import dependencies
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';OUT=ART/'retention24_discovery_evidence_v1.json'
FOLDERS=('r17-retention512-benchmark-v1','compact-six-r17-retention512-v1','r17-retention512-extended-v1','retention24-r17-pari-v1','native11952-rank28-coset-visibility-v1','retention-cross-family-v1','retention24-current-catalogue-v1')
def main():
    archive=OUT.with_suffix('.zip')
    if OUT.exists() or archive.exists():raise FileExistsError('preserve height/discarded evidence')
    if cert.read(LOCAL/'retention24-r17-pari-v1/verification-ledger.json')['status']!='PASS' or cert.read(ART/'native11952_rank28_coset_visibility_replay_v1.json')['status']!='PASS' or cert.read(ART/'new_high_rank_curve_index_v7_memory_replay_v1.json')['status']!='PASS':raise ArithmeticError('required finite checks incomplete')
    prior=ART/'height_and_discarded_evidence_v1.json';pm=cert.read(prior);base_names=[r['manifest'] for r in pm['required_base_archives']]+[str(prior.relative_to(ROOT))];bases=[];base_members={}
    for name in base_names:
        p=ROOT/name;m=cert.read(p)
        if cert.hashed(ROOT/m['archive'])!=m['archive_sha256']:raise ArithmeticError('base archive changed')
        base_members.update({r['path']:r['sha256'] for r in m['files']});bases.append({'manifest':name,'manifest_sha256':cert.hashed(p),'archive':m['archive'],'archive_sha256':m['archive_sha256']})
    paths={Path(__file__).resolve(),ROOT/'elliptic-curves/notes/WIDER_RETENTION_DISCOVERIES_2026-09-06.md'}
    for name in FOLDERS:paths.update(p for p in (LOCAL/name).rglob('*') if p.is_file())
    names=('native11952_rank28_coset_visibility_v1.json','native11952_rank28_coset_visibility_replay_v1.json','retention24_r17_results_v1.json','retention_high_rank_minimal_proofs_v1.json','new_high_rank_curve_index_v7.json','new_high_rank_curve_index_v7.csv','new_high_rank_curve_index_v7_memory_replay_v1.json','height_and_discarded_portable_replay_v1.json','retention_refreshed_catalogue_comparison_v1.json','retention_cross_family_j_incidence_v1.json','retention_cross_family_j_incidence_replay_v1.json','inventory89_cross_family_incidence_v1.json')
    paths.update(ART/n for n in names);paths.update(ART.glob('retention24_r17_*_mod2_v1.json'));paths.update(ART.glob('new_retention_rank*_curve*.sage'))
    names=('audit_native11952_rank28_coset_visibility.sage','replay_native11952_rank28_coset_visibility.py','benchmark_r17_retention512.py','retain512_compact_r17.py','extend_retention512_r17_scores.py','retention24_r17_pari_batch.py','prepare_retention24_r17_pari_batch.sage','verify_retention24_r17_pari_batch.py','certify_retention24_r17_results.py','certify_retention_high_rank_minimal.py','export_new_high_rank_curve_index_v7.py','replay_inventory_v7_memory.py','replay_retention24_geometry.py','verify_retention24_discovery_bundle.py','refresh_retention_catalogue.py','audit_retention_cross_family_incidence.sage','replay_retention_cross_family_incidence.py','certify_inventory89_incidence.py','replay_compact_cross_family_incidence.py','replay_latest7_cross_family_incidence.py','replay_latest8_cross_family_incidence.py','replay_latest23_cross_family_incidence.py')
    paths.update(ROOT/'elliptic-curves/cas'/n for n in names);seen=set()
    while True:
        paths=dependencies(paths);todo=[p for p in paths-seen if p.suffix=='.json']
        if not todo:break
        for p in todo:
            seen.add(p);d=cert.read(p)
            if not isinstance(d,dict):continue
            for key in ('sources','source_hashes','inputs','bindings','shard_hashes','canonical_table_hashes','old_scored_hashes'):
                value=d.get(key,{})
                if isinstance(value,dict):paths.update(ROOT/name for name in value if (ROOT/name).is_file())
            for name in d.get('source_certificate_hashes',{}):
                if (ART/name).is_file():paths.add(ART/name)
            for key in ('input_path','seed_path','maps_path'):
                if isinstance(d.get(key),str) and (ROOT/d[key]).is_file():paths.add(ROOT/d[key])
    members=[];inherited=[]
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for p in sorted(paths):
            raw=p.read_bytes();name=str(p.relative_to(ROOT));h=sha256(raw).hexdigest();row={'path':name,'sha256':h,'bytes':len(raw)}
            if base_members.get(name)==h:inherited.append(row);continue
            info=zipfile.ZipInfo(name,date_time=(2026,9,6,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o644<<16;z.writestr(info,raw);members.append(row)
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None or any(sha256(z.read(r['path'])).hexdigest()!=r['sha256'] for r in members):raise ArithmeticError('archive integrity differs')
    cert.write(OUT,{'schema':'elliptic-curves.retention24-discovery-evidence.v1','archive':str(archive.relative_to(ROOT)),'archive_sha256':cert.hashed(archive),'archive_bytes':archive.stat().st_size,'required_base_archives':bases,'files':members,'inherited_exact_members':inherited,'claim_boundary':'Extract ten pinned bases in order then this supplement. Retains6144 wider short-score finalists from the same H4096 population,4608 new extended-score trace rosters, all1080 charts and histories on the fixed24 selected candidates, the current inventory and high-rank minimal proofs, plus98 separate control-only oracle translations. All rank claims are independent-point lower bounds; no exact rank, universal novelty or general selector guarantee.'})
    print('PACKAGED RETENTION24',len(members),'new',len(inherited),'inherited',archive.stat().st_size,'bytes',flush=True)
if __name__=='__main__':main()
