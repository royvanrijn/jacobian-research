#!/usr/bin/env python3
"""Portable supplement for calibrated PARI searches and exact admission improvements."""
from pathlib import Path
from hashlib import sha256
import zipfile
import certify_compact_r17_candidates as cert
from package_recorded_mod2_audit import dependencies
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';OUTPUT=ART/'fast_point_pipeline_evidence_v1.json'
FOLDERS=('new-rank26-engine-comparison-v1','new-rank26-engine-comparison-v2','new-rank26-adaptive-pari-maps-v1','new-rank26-fast301-v1','rank26-primebank-continuation-v1','rank26-complete-cloud-modl-v1','small-conductor-adaptive-pari-maps-v1','small-conductor-fast301-v1','small-conductor-primebank-continuation-v1','preloaded-prime-state-benchmark-v1','small-conductor-complete-cloud-modl-v1')


def main():
    archive=OUTPUT.with_suffix('.zip')
    if OUTPUT.exists() or archive.exists():raise FileExistsError('preserve fast pipeline bundle')
    for case in ('rank26','small-conductor'):
        d=LOCAL/(case+'-primebank-continuation-v1')
        for name in ('charts-result.json','completion-result.json'):
            r=cert.read(d/'verification'/name)
            if r['outcome']!='completed' or r['returncode']!=0:raise ArithmeticError('continuation replay incomplete')
    for name in ('rank26-complete-cloud-modl-v1','small-conductor-complete-cloud-modl-v1'):
        r=cert.read(LOCAL/name/'check-result.json')
        if r['outcome']!='completed' or r['returncode']!=0:raise ArithmeticError('odd-modulus replay incomplete')
    base_members={};bases=[]
    for name in ('compact_r17_wide_evidence_v1.json','new_rank26_followup_evidence_v1.json','exact_parity_coordinate_evidence_v1.json','small_conductor_followup_evidence_v1.json'):
        path=ART/name;m=cert.read(path)
        if cert.hashed(ROOT/m['archive'])!=m['archive_sha256']:raise ArithmeticError('base archive changed')
        base_members.update({r['path']:r['sha256'] for r in m['files']});bases.append({'manifest':str(path.relative_to(ROOT)),'manifest_sha256':cert.hashed(path),'archive':m['archive'],'archive_sha256':m['archive_sha256']})
    paths={Path(__file__).resolve(),ROOT/'elliptic-curves/notes/FAST_POINT_PIPELINE_AUDIT_2026-09-06.md'}
    for name in FOLDERS:paths.update(p for p in (LOCAL/name).rglob('*') if p.is_file())
    names=('new_rank26_engine_comparison_v1.json','new_rank26_engine_comparison_v2.json','pari_height_boundary_regression_v1.json','new_rank26_fast301_recorded_mod2_v1.json','small_conductor_fast301_recorded_mod2_v1.json','rank26_fast301_completion_v1.json','small_conductor_fast301_completion_v1.json','rank26_all_current_retained_mod2_v1.json','small_conductor_all_current_retained_mod2_v1.json','rank26_all_retained_modl_v1.json','small_conductor_all_retained_modl_v1.json','strict_backend_calibration_replay_v1.json')
    paths.update(ART/name for name in names)
    scripts=('compare_new_rank26_search_engines.py','compare_new_rank26_search_engines_v2.py','pari_pointed_backend.py','prepare_new_rank26_adaptive_pari_maps.sage','prepare_small_conductor_pari_maps.sage','search_new_rank26_fast301.py','replay_new_rank26_fast301.py','verify_new_rank26_fast301.py','search_small_conductor_fast301.py','replay_small_conductor_fast301.py','verify_small_conductor_fast301.py','continue_fixed_pari_search.py','audit_fixed_pari_completion.py','benchmark_preloaded_prime_state.py','audit_retained_cloud_modl.py','check_pari_height_boundaries.py','replay_backend_calibrations_strict.py','replay_all_fast_backend_geometry.py','verify_fast_point_pipeline_bundle.py')
    paths.update(ROOT/'elliptic-curves/cas'/name for name in scripts)
    paths.update(ROOT/'elliptic-curves/tests'/name for name in ('test_pari_pointed_backend.py','test_preloaded_prime_state.py','test_rotated_observation_state.py','test_retained_cloud_modl.py','test_strict_calibration_roster.py'))
    seen=set()
    while True:
        paths=dependencies(paths);todo=[p for p in paths-seen if p.suffix=='.json']
        if not todo:break
        for path in todo:
            seen.add(path);data=cert.read(path)
            if not isinstance(data,dict):continue
            for key in ('sources','source_hashes','inputs'):
                value=data.get(key,{})
                if isinstance(value,dict):
                    paths.update(ROOT/p for p in value if (ROOT/p).is_file())
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
    cert.write(OUTPUT,{'schema':'elliptic-curves.fast-point-pipeline-evidence.v1','archive':str(archive.relative_to(ROOT)),'archive_sha256':cert.hashed(archive),'archive_bytes':archive.stat().st_size,'required_base_archives':bases,'files':members,'inherited_exact_members':inherited,'claim_boundary':'Extract the four pinned bases in order, then this supplement. Retains failed comparator, timed-out and intentionally stopped experiments, exact continuations, archived states, full point-cloud proofs, backend boundary checks and admission tests. No new rank28/32 or rank23 small-conductor result is implied.'})
    print('PACKAGED FAST POINT PIPELINE',len(members),'new members',len(inherited),'inherited',archive.stat().st_size,'bytes',flush=True)
if __name__=='__main__':main()
