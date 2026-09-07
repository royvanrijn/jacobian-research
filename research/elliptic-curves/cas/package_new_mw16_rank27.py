#!/usr/bin/env python3
"""Portable wider-retention discovery and exact control-diagnostic evidence."""
from pathlib import Path
from hashlib import sha256
import zipfile
import certify_compact_r17_candidates as cert
from package_recorded_mod2_audit import dependencies
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';OUT=ART/'new_mw16_rank27_evidence_v1.json'
FOLDERS=('mw16-new26-followups-v1','mw16-new26-a1-fibration-01-015-adaptive-v1','mw16-new26-a1-fibration-01-052-adaptive-v1','mw16-new26-a1-fibration-02-014-adaptive-v1','mw16-new27-a1-fibration-01-052-newdirection-v1','periodic-nagao-scanner-benchmark-v1','periodic-nagao-strict-v1')

def main():
    archive=OUT.with_suffix('.zip')
    if OUT.exists() or archive.exists():raise FileExistsError('preserve height/discarded evidence')
    if cert.read(ART/'new_mw16_rank27_minimal_proof_v1.json')['status']!='PASS' or cert.read(ART/'new_high_rank_curve_index_v9_memory_replay_v1.json')['status']!='PASS':raise ArithmeticError('exact new27 and inventory gates required')
    prior=ART/'million_height_mw16_evidence_v1.json';pm=cert.read(prior);base_names=[r['manifest'] for r in pm['required_base_archives']]+[str(prior.relative_to(ROOT)),str((ART/'million_height_mw16_population_supplement_v1.json').relative_to(ROOT))];bases=[];base_members={}
    for name in base_names:
        p=ROOT/name;m=cert.read(p)
        if cert.hashed(ROOT/m['archive'])!=m['archive_sha256']:raise ArithmeticError('base archive changed')
        base_members.update({r['path']:r['sha256'] for r in m['files']});bases.append({'manifest':name,'manifest_sha256':cert.hashed(p),'archive':m['archive'],'archive_sha256':m['archive_sha256']})
    paths={Path(__file__).resolve(),ROOT/'elliptic-curves/notes/NEW_MW16_RANK27_2026-09-06.md'}
    for name in FOLDERS:paths.update(p for p in (LOCAL/name).rglob('*') if p.is_file())
    for pattern in ('new_mw16_followup_results_v1.json','new_mw16_rank27_minimal_proof_v1.json','mw16_new26_*_v1.json','mw16_new27_direction_*_v1.json','new_high_rank_curve_index_v9*','new_mw16_rank27_curve_a1_01.sage','inventory98_cross_family_incidence_v2.json','public_compact_parameter_heights_v1.json','periodic_nagao_scanner*_v1.json'):
        paths.update(ART.glob(pattern))
    names=('followup_new_mw16_rank26.py','prepare_new_mw16_rank26_adaptive.sage','run_new_mw16_rank26_followups.py','audit_new_mw16_rank26_followup.py','replay_new_mw16_rank26_followup_geometry.py','certify_new_mw16_followup_results.py','certify_new_mw16_rank27_minimal.py','export_new_high_rank_curve_index_v9.py','replay_inventory_v9_memory.py','certify_inventory98_incidence_v2.py','export_new_mw16_rank27_sage.py','followup_new_mw16_rank27_direction.py','prepare_new_mw16_rank27_direction.sage','audit_new_mw16_rank27_direction.py','replay_new_mw16_rank27_direction_geometry.py','audit_public_compact_parameter_heights.py','benchmark_periodic_nagao_scanner.py','verify_periodic_nagao_scanner.py','replay_periodic_nagao_scanner_portable.py','newfamily/scan_rational_nagao_tables_v2.cpp','verify_new_mw16_rank27_bundle.py')
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
    cert.write(OUT,{'schema':'elliptic-curves.new-mw16-rank27-evidence.v1','archive':str(archive.relative_to(ROOT)),'archive_sha256':cert.hashed(archive),'archive_bytes':archive.stat().st_size,'required_base_archives':bases,'files':members,'inherited_exact_members':inherited,'claim_boundary':'Extract fourteen pinned base archives in order and this supplement. Retains903 three-curve adaptive boxes,301 new27-direction boxes, exact histories and complete cloud proofs, one minimal27-point curve, unchanged98 IDs with six27s, and incidence rank-metadata binding. Also contains the retrospective67-of69 compact coordinate audit and strict periodic Nagao benchmark. The historical calibration dataset is opaque pinned metadata: only this audit coordinate inversion and Q-isomorphism are rechecked, not historical ranks or parent geometry. The higher-height prospective population and point searches have separate evidence and are excluded. No exact rank, new28/32, conductor record or universal novelty.'})
    print('PACKAGED NEW MW16 RANK27',len(members),'new',len(inherited),'inherited',archive.stat().st_size,'bytes',flush=True)
if __name__=='__main__':main()
