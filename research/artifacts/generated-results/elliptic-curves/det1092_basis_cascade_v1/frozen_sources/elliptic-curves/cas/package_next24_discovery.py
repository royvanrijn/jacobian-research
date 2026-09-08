#!/usr/bin/env python3
"""Portable supplement for the third new rank27 curve and their finite discovery experiment."""
from pathlib import Path
from hashlib import sha256
import zipfile
import certify_compact_r17_candidates as cert
from package_recorded_mod2_audit import dependencies
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';OUTPUT=ART/'next24_discovery_evidence_v1.json'
FOLDERS=('next24-r17-extended-pari-v1','memory-cloud-audit-regression-v3','native11952-published-visibility-v1','native11952-published-visibility-v2','next24-high-rank-minimal-v2','next24-high-rank-minimal-v3','next24-high-rank-minimal-v4','next24-rank27-adaptive-v1','next24-certification-v1')


def main():
    archive=OUTPUT.with_suffix('.zip')
    if OUTPUT.exists() or archive.exists():raise FileExistsError('preserve fast pipeline bundle')
    if cert.read(LOCAL/'next24-r17-extended-pari-v1/verification-ledger.json')['status']!='PASS' or cert.read(LOCAL/'next24-rank27-adaptive-v1/cloud-verification-ledger.json')['status']!='PASS':raise ArithmeticError('histories/clouds incomplete')
    base_members={};bases=[]
    for name in ('compact_r17_wide_evidence_v1.json','new_rank26_followup_evidence_v1.json','exact_parity_coordinate_evidence_v1.json','small_conductor_followup_evidence_v1.json','fast_point_pipeline_evidence_v1.json','paired_rank27_discovery_evidence_v1.json','latest8_cross_family_evidence_v1.json','paired_rank27_inventory_replay_evidence_v2.json'):
        path=ART/name;m=cert.read(path)
        if cert.hashed(ROOT/m['archive'])!=m['archive_sha256']:raise ArithmeticError('base archive changed')
        base_members.update({r['path']:r['sha256'] for r in m['files']});bases.append({'manifest':str(path.relative_to(ROOT)),'manifest_sha256':cert.hashed(path),'archive':m['archive'],'archive_sha256':m['archive_sha256']})
    paths={Path(__file__).resolve(),ROOT/'elliptic-curves/notes/NEXT24_RANK27_DISCOVERY_2026-09-06.md'}
    for name in FOLDERS:paths.update(p for p in (LOCAL/name).rglob('*') if p.is_file())
    names=('next24_r17_results_v1.json','next24_high_rank_minimal_proofs_v3.json','next24_high_rank_minimal_proofs_v4.json','new_high_rank_curve_index_v5.json','new_high_rank_curve_index_v5.csv','new_high_rank_curve_index_v5_memory_replay_v1.json','native11952_published_visibility_v1.json','native11952_published_visibility_v2.json','next24_rank27_adaptive_coverage_v1.json','next24_rank27_all_retained_mod2_v1.json','next24_rank27_all_retained_modl_v1.json','paired_rank27_portable_completion_v2.json')
    paths.update(ART/name for name in names);paths.update(ART.glob('next24_r17_extended_*_mod2_v1.json'));paths.update(ART.glob('new_next24_*_curve.sage'))
    scripts=('next24_r17_pari_batch.py','prepare_next24_r17_pari_batch.sage','verify_next24_r17_pari_batch.py','audit_recorded_point_mod2_rank_v3.py','certify_next24_r17_results.py','export_new_high_rank_curve_index_v5.py','replay_inventory_v5_memory.py','certify_next24_high_rank_minimal.py','certify_next24_high_rank_minimal_v2.py','certify_next24_high_rank_minimal_v3.py','certify_next24_high_rank_minimal_v4.py','audit_native11952_published_visibility.py','audit_native11952_published_visibility_v2.py','followup_next24_rank27.py','prepare_next24_rank27_adaptive.sage','audit_next24_rank27_followup.py','replay_next24_point_geometry.py','verify_next24_discovery_bundle.py')
    paths.update(ROOT/'elliptic-curves/cas'/name for name in scripts)
    paths.add(ROOT/'elliptic-curves/tests/test_next24_local_minimality.py')
    seen=set()
    while True:
        paths=dependencies(paths);todo=[p for p in paths-seen if p.suffix=='.json']
        if not todo:break
        for path in todo:
            seen.add(path);data=cert.read(path)
            if not isinstance(data,dict):continue
            for key in ('sources','source_hashes','inputs','bindings'):
                value=data.get(key,{})
                if isinstance(value,dict):
                    paths.update(ROOT/p for p in value if (ROOT/p).is_file())
            for name in data.get('source_certificate_hashes',{}):
                if (ART/name).is_file():paths.add(ART/name)
            for key in ('input_path','seed_path','maps_path'):
                if isinstance(data.get(key),str) and (ROOT/data[key]).is_file():paths.add(ROOT/data[key])
            for row in data.get('source_inputs',[]):
                if (ROOT/row['path']).is_file():paths.add(ROOT/row['path'])
    members=[];inherited=[]
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for path in sorted(paths):
            raw=path.read_bytes();name=str(path.relative_to(ROOT));h=sha256(raw).hexdigest();row={'path':name,'bytes':len(raw),'sha256':h}
            if base_members.get(name)==h:inherited.append(row);continue
            info=zipfile.ZipInfo(name,date_time=(2026,9,6,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o644<<16;z.writestr(info,raw);members.append(row)
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None or any(sha256(z.read(r['path'])).hexdigest()!=r['sha256'] for r in members):raise ArithmeticError('archive integrity failure')
    cert.write(OUTPUT,{'schema':'elliptic-curves.next24-discovery-evidence.v1','archive':str(archive.relative_to(ROOT)),'archive_sha256':cert.hashed(archive),'archive_bytes':archive.stat().st_size,'required_base_archives':bases,'files':members,'inherited_exact_members':inherited,'claim_boundary':'Extract the eight pinned bases in order, then this supplement. Retains all24 fixed candidate attempts and1080 charts, one301-chart adaptive follow-up, all state archives and local replay outcomes,62-curve inventory, three newly certified globally minimal rank26/27 models, and a separate published29-point visibility audit. Failed minimality builders and the JSON tuple/list checker failure remain retained. No exact ranks, universal novelty or new rank28/32 claim.'})
    print('PACKAGED NEXT24 DISCOVERY',len(members),'new members',len(inherited),'inherited',archive.stat().st_size,'bytes',flush=True)
if __name__=='__main__':main()
