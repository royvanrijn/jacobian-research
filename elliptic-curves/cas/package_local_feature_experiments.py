#!/usr/bin/env python3
"""Portable local-feature cohorts and the new101st curve with25 independent points."""
from pathlib import Path
from hashlib import sha256
import zipfile
import certify_compact_r17_candidates as cert
from package_recorded_mod2_audit import dependencies
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';OUT=ART/'local_feature_evidence_v1.json'
FOLDERS=('higher13-scaled-selection-v1','higher-split-node-score-v1','scaled13_24-r17-pari-v1','splitnode24-r17-pari-v1','scaled13-rank25-11952-300-adaptive-v1','inventory101-promotion-v1')

def main():
    archive=OUT.with_suffix('.zip')
    if OUT.exists() or archive.exists():raise FileExistsError('preserve height/discarded evidence')
    if cert.read(ART/'local_feature_experiments_v1.json')['status']!='PASS':raise ArithmeticError('exact new26 and inventory gates required')
    prior=ART/'product_first_evidence_v1.json';pm=cert.read(prior);base_names=[r['manifest'] for r in pm['required_base_archives']]+[str(prior.relative_to(ROOT))];bases=[];base_members={}
    for name in base_names:
        p=ROOT/name;m=cert.read(p)
        if cert.hashed(ROOT/m['archive'])!=m['archive_sha256']:raise ArithmeticError('base archive changed')
        base_members.update({r['path']:r['sha256'] for r in m['files']});bases.append({'manifest':name,'manifest_sha256':cert.hashed(p),'archive':m['archive'],'archive_sha256':m['archive_sha256']})
    paths={Path(__file__).resolve(),ROOT/'elliptic-curves/notes/LOCAL_FEATURE_RANK25_2026-09-06.md'}
    for name in FOLDERS:paths.update(p for p in (LOCAL/name).rglob('*') if p.is_file())
    for pattern in ('scaled13*_*v1.json','splitnode24_*_v1.json','local_feature_experiments_v1.json','inventory101_*_v1.json','new_high_rank_curve_index_v12.*','new_high_rank_curve_index_v12_memory_replay_v1.json','new_scaled13_rank25_curve_11952.sage'):paths.update(ART.glob(pattern))
    names=('select_higher13_scaled_candidates.py','score_higher_split_nodes.py','scaled13_24_r17_pari_batch.py','prepare_scaled13_24_r17_pari_batch.sage','verify_scaled13_24_r17_pari_batch.py','certify_scaled13_24_r17_results.py','replay_scaled13_24_geometry.py','audit_scaled13_24_clouds_modl.py','splitnode24_r17_pari_batch.py','prepare_splitnode24_r17_pari_batch.sage','verify_splitnode24_r17_pari_batch.py','certify_splitnode24_r17_results.py','replay_splitnode24_geometry.py','audit_splitnode24_clouds_modl.py','followup_scaled13_rank25.py','prepare_scaled13_rank25_adaptive.sage','audit_scaled13_rank25_followup.py','replay_scaled13_rank25_followup_geometry.py','certify_scaled13_25_followup_results.py','certify_scaled13_24_rank25_minimal.py','export_scaled13_rank25_sage.py','export_new_high_rank_curve_index_v12.py','replay_inventory_v12_memory.py','audit_inventory101_added_incidence.sage','replay_inventory101_added_incidence.py','certify_inventory101_incidence.py','report_local_feature_experiments.py','verify_local_feature_bundle.py')
    paths.update(ROOT/'elliptic-curves/cas'/n for n in names);seen=set()
    while True:
        paths=dependencies(paths);todo=[p for p in paths-seen if p.suffix=='.json']
        if not todo:break
        for p in todo:
            seen.add(p);d=cert.read(p)
            if not isinstance(d,dict):continue
            if p.name=='elkies-k3-r17-norm12-icarm-calibration-dataset-v1.json':continue # Retrospective metadata only; this audit does not reprove its historical ranks or geometry.
            for key in ('sources','source_hashes','inputs','bindings','shard_hashes','canonical_table_hashes','old_scored_hashes','short_table_hashes'):
                value=d.get(key,{})
                if isinstance(value,dict):paths.update(ROOT/name for name in value if (ROOT/name).is_file())
            if isinstance(d.get('table'),dict) and (ROOT/d['table'].get('path','')).is_file():paths.add(ROOT/d['table']['path'])
            for name in d.get('source_certificate_hashes',{}):
                if (ART/name).is_file():paths.add(ART/name)
            for key in ('input_path','seed_path','maps_path','generic_census_path','supersedes_failed_protocol'):
                if isinstance(d.get(key),str) and (ROOT/d[key]).is_file():paths.add(ROOT/d[key])
    members=[];inherited=[]
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for p in sorted(paths):
            raw=p.read_bytes();name=str(p.relative_to(ROOT));h=sha256(raw).hexdigest();row={'path':name,'sha256':h,'bytes':len(raw)}
            if base_members.get(name)==h:inherited.append(row);continue
            info=zipfile.ZipInfo(name,date_time=(2026,9,6,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o644<<16;z.writestr(info,raw);members.append(row)
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None or any(sha256(z.read(r['path'])).hexdigest()!=r['sha256'] for r in members):raise ArithmeticError('archive integrity differs')
    cert.write(OUT,{'schema':'elliptic-curves.local-feature-evidence.v1','archive':str(archive.relative_to(ROOT)),'archive_sha256':cert.hashed(archive),'archive_bytes':archive.stat().st_size,'required_base_archives':bases,'files':members,'inherited_exact_members':inherited,'claim_boundary':'Extract nineteen pinned base archives and this supplement. Retains two frozen local-feature selectors,48 point curves with2160 completed initial boxes,301 completed follow-up boxes and its3260-point union, all full-cloud proofs modulo2,3,5, the new25 minimal equation and101-curve inventory with1212 incidence pairs. No exact rank, new near-record, upper rank, saturation, selector superiority or universal novelty.'})
    print('PACKAGED LOCAL FEATURES',len(members),'new',len(inherited),'inherited',archive.stat().st_size,'bytes',flush=True)
if __name__=='__main__':main()
