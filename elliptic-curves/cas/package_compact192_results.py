#!/usr/bin/env python3
"""Portable fixed192 compact point cohort, inventory and universal13 proof."""
from pathlib import Path
from hashlib import sha256
import zipfile
import certify_compact_r17_candidates as cert
from package_recorded_mod2_audit import dependencies
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';OUT=ART/'compact192_evidence_v1.json'
FOLDERS=('compact192-unsearched-selection-v1','compact192-r17-pari-v1','r17-13-scaling-classification-v1')

def main():
    archive=OUT.with_suffix('.zip')
    if OUT.exists() or archive.exists():raise FileExistsError('preserve height/discarded evidence')
    import export_new_high_rank_curve_index_v13 as inventory
    inventory.promotion_gate()
    replay=cert.read(ART/'new_high_rank_curve_index_v13_memory_replay_v1.json')
    if replay['status']!='PASS' or any(cert.hashed(ROOT/n)!=h for n,h in replay['sources'].items()):raise ArithmeticError('completed boundV13 inventory replay required')
    prior=ART/'endpoint_point_trial_evidence_v1.json';pm=cert.read(prior);base_names=[r['manifest'] for r in pm['required_base_archives']]+[str(prior.relative_to(ROOT))];bases=[];base_members={}
    for name in base_names:
        p=ROOT/name;m=cert.read(p)
        if cert.hashed(ROOT/m['archive'])!=m['archive_sha256']:raise ArithmeticError('base archive changed')
        base_members.update({r['path']:r['sha256'] for r in m['files']});bases.append({'manifest':name,'manifest_sha256':cert.hashed(p),'archive':m['archive'],'archive_sha256':m['archive_sha256']})
    paths={Path(__file__).resolve(),ROOT/'elliptic-curves/cas/finish_compact192_evidence.py',LOCAL/'compact192-evidence-finalization-v1/protocol.json',ROOT/'elliptic-curves/notes/COMPACT192_UNSEARCHED_TRIAL_2026-09-06.md',ROOT/'elliptic-curves/notes/R17_INTEGRAL_13_PARAMETER_CHARTS_2026-09-06.md'}
    for name in FOLDERS:paths.update(p for p in (LOCAL/name).rglob('*') if p.is_file())
    for pattern in ('compact192_*_v1.json','compact192_first26_candidate.sage','new_high_rank_curve_index_v13.*','new_high_rank_curve_index_v13_memory_replay_v1.json','r17_integral13_parameter_charts_v1.json','r17_13_scaling_*_v*.json','endpoint_point_trial_portable_replay_v1.json'):paths.update(ART.glob(pattern))
    names=('select_compact192_unsearched.py', 'compact192_r17_pari_batch.py', 'prepare_compact192_r17_pari_batch.sage', 'verify_compact192_r17_pari_batch.py', 'stream_compact192_verification.py', 'certify_compact192_r17_results.py', 'replay_compact192_geometry.py', 'audit_compact192_clouds_modl.py', 'report_compact192_experiment.py', 'finalize_compact192_points.py', 'export_new_high_rank_curve_index_v13.py', 'replay_inventory_v13_memory.py', 'export_compact192_first26_candidate.py', 'audit_r17_integral13_charts.py', 'verify_r17_integral13_charts.sage', 'classify_r17_13_scaling.py', 'verify_r17_13_scaling_classification.py', 'report_r17_13_scaling_geometry.py', 'report_r17_13_scaling_geometry_v2.py', 'verify_compact192_bundle.py', 'report_compact192_portable.py')
    paths.update(ROOT/'elliptic-curves/cas'/n for n in names);seen=set()
    while True:
        paths=dependencies(paths);todo=[p for p in paths-seen if p.suffix=='.json']
        if not todo:break
        for p in todo:
            seen.add(p);d=cert.read(p)
            if not isinstance(d,dict):continue
            if p.name=='elkies-k3-r17-norm12-icarm-calibration-dataset-v1.json':continue # Retrospective metadata only; this audit does not reprove its historical ranks or geometry.
            for key in ('sources','source_hashes','inputs','bindings','shard_hashes','canonical_table_hashes','old_scored_hashes','short_table_hashes','cache_result_hashes','seed_hashes'):
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
    cert.write(OUT,{'schema':'elliptic-curves.compact192-evidence.v1','archive':str(archive.relative_to(ROOT)),'archive_sha256':cert.hashed(archive),'archive_bytes':archive.stat().st_size,'required_base_archives':bases,'files':members,'inherited_exact_members':inherited,'claim_boundary':'Extract twenty-four pinned base archives and this supplement. Retain the fixed192 distinct compactH4096 candidate selection, all pre-search generic maps and bounded point attempts, exact point histories and full clouds modulo2,3,5, post-terminal pinned catalogue comparisons and V13 inventory proofs. Also retain the exact84-cell parameter-chart proof and complete universal13-scaling and bad-reduction classification, with both aggregate versions. No new point search occurs in replay. No whole-curve rank upper bound, full saturation, complete parameter coverage, score optimality, record claim or universal novelty.'})
    print('PACKAGED COMPACT192 RESULTS',len(members),'new',len(inherited),'inherited',archive.stat().st_size,'bytes',flush=True)
if __name__=='__main__':main()
