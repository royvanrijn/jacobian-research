#!/usr/bin/env python3
"""Portable wider-retention discovery and exact control-diagnostic evidence."""
from pathlib import Path
from hashlib import sha256
import zipfile
import certify_compact_r17_candidates as cert
from package_recorded_mod2_audit import dependencies
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';OUT=ART/'higher_r17_rank26_evidence_v1.json'
FOLDERS=('higher32768-r17-stratified-v1','higher32768-r17-extended-v1','higher24-r17-pari-v1','higher24-incidence-v1','higher-rank25-11952-069-adaptive-v1')

def main():
    archive=OUT.with_suffix('.zip')
    if OUT.exists() or archive.exists():raise FileExistsError('preserve height/discarded evidence')
    if cert.read(ART/'higher26_minimal_proof_v1.json')['status']!='PASS' or cert.read(ART/'new_high_rank_curve_index_v11_memory_replay_v1.json')['status']!='PASS':raise ArithmeticError('exact new26 and inventory gates required')
    prior=ART/'new_mw16_rank27_evidence_v1.json';pm=cert.read(prior);base_names=[r['manifest'] for r in pm['required_base_archives']]+[str(prior.relative_to(ROOT))];bases=[];base_members={}
    for name in base_names:
        p=ROOT/name;m=cert.read(p)
        if cert.hashed(ROOT/m['archive'])!=m['archive_sha256']:raise ArithmeticError('base archive changed')
        base_members.update({r['path']:r['sha256'] for r in m['files']});bases.append({'manifest':name,'manifest_sha256':cert.hashed(p),'archive':m['archive'],'archive_sha256':m['archive_sha256']})
    paths={Path(__file__).resolve(),ROOT/'elliptic-curves/notes/HIGHER_PARAMETER_RANK26_2026-09-06.md'}
    for name in FOLDERS:paths.update(p for p in (LOCAL/name).rglob('*') if p.is_file())
    for pattern in ('higher24_r17*_v1.json','higher24_rank25_minimal_proof_v1.json','higher25_followup_results_v1.json','higher26_minimal_proof_v1.json','higher_rank25_*_v1.json','higher24_visibility_cost_v1.json','higher24_cross_family*_v1.json','inventory100_cross_family_incidence_v*.json','new_high_rank_curve_index_v10*','new_high_rank_curve_index_v11*','new_higher_rank26_curve_11952.sage','new_mw16_rank27_portable_replay_v1.json'):
        paths.update(ART.glob(pattern))
    names=('scan_higher_r17_stratified.py','extend_higher_r17_stratified.py','replay_higher_r17_stratified_portable.py','higher24_r17_pari_batch.py','prepare_higher24_r17_pari_batch.sage','verify_higher24_r17_pari_batch.py','certify_higher24_r17_results.py','replay_higher24_geometry.py','audit_higher24_visibility_cost.py','certify_higher24_rank25_minimal.py','followup_higher_rank25.py','prepare_higher_rank25_adaptive.sage','audit_higher_rank25_followup.py','replay_higher_rank25_followup_geometry.py','certify_higher25_followup_results.py','certify_higher26_minimal.py','export_new_high_rank_curve_index_v10.py','replay_inventory_v10_memory.py','export_new_high_rank_curve_index_v11.py','replay_inventory_v11_memory.py','audit_higher24_incidence.sage','replay_higher24_incidence.py','certify_inventory100_incidence.py','certify_inventory100_incidence_v2.py','export_higher26_sage.py','verify_higher_r17_rank26_bundle.py','report_new_mw16_rank27_portable.py')
    paths.update(ROOT/'elliptic-curves/cas'/n for n in names);seen=set()
    while True:
        paths=dependencies(paths);todo=[p for p in paths-seen if p.suffix=='.json']
        if not todo:break
        for p in todo:
            seen.add(p);d=cert.read(p)
            if not isinstance(d,dict):continue
            if p.name=='elkies-k3-r17-norm12-icarm-calibration-dataset-v1.json':continue # Retrospective metadata only; this audit does not reprove its historical ranks or geometry.
            for key in ('sources','source_hashes','inputs','bindings','shard_hashes','canonical_table_hashes','old_scored_hashes'):
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
    cert.write(OUT,{'schema':'elliptic-curves.higher-r17-rank26-evidence.v1','archive':str(archive.relative_to(ROOT)),'archive_sha256':cert.hashed(archive),'archive_bytes':archive.stat().st_size,'required_base_archives':bases,'files':members,'inherited_exact_members':inherited,'claim_boundary':'Extract fifteen pinned base archives and this supplement. Retains the122368792-address twelve-slice population,6144 saved short and extended-prime score rosters,24 generic17-only point attempts and full clouds, a301-centre adaptive gain25to26,100-curve inventory and1200 incidence pairs, exact minimal26 model and standalone point file, and descriptive coefficient/visibility comparisons. Higher26 new-direction experiments have separate evidence. No full32768 coverage, exact rank, new28/32, conductor record, calibrated score predictor or universal novelty.'})
    print('PACKAGED HIGHER R17 RANK26',len(members),'new',len(inherited),'inherited',archive.stat().st_size,'bytes',flush=True)
if __name__=='__main__':main()
