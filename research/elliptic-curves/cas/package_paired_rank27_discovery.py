#!/usr/bin/env python3
"""Portable supplement for two new rank27 curves and their finite discovery experiment."""
from pathlib import Path
from hashlib import sha256
import zipfile
import certify_compact_r17_candidates as cert
from package_recorded_mod2_audit import dependencies
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';OUTPUT=ART/'paired_rank27_discovery_evidence_v1.json'
FOLDERS=('r17-extended-prime-benchmark-v1','r17-retained-extended-primes-v1','fresh-r17-paired-pari-v1','fresh-r17-paired-certification-v1','new-paired-rank27-certification-v1','paired-high-rank-minimal-v1','paired-high-rank-minimal-v2','paired-rank27-conductor-bound-v1','paired-rank27-adaptive-v1','paired-second27-adaptive-v1','paired-point-geometry-replay-v1')


def main():
    archive=OUTPUT.with_suffix('.zip')
    if OUTPUT.exists() or archive.exists():raise FileExistsError('preserve fast pipeline bundle')
    if cert.read(LOCAL/'fresh-r17-paired-pari-v1/verification-ledger.json')['status']!='PASS':raise ArithmeticError('paired histories/clouds not verified')
    for folder,name in [('r17-retained-extended-primes-v1','replay.supervisor.json'),('paired-high-rank-minimal-v2','check.supervisor.json'),('fresh-r17-paired-certification-v1','check.supervisor.json'),('paired-point-geometry-replay-v1','replay.supervisor.json')]:
        r=cert.read(LOCAL/folder/name)
        if r['outcome']!='completed' or r['returncode']!=0:raise ArithmeticError('required exact replay incomplete')
    for folder in ('paired-rank27-adaptive-v1','paired-second27-adaptive-v1'):
        if cert.read(LOCAL/folder/'cloud-verification-ledger.json')['status']!='PASS':raise ArithmeticError('adaptive clouds incomplete')
    base_members={};bases=[]
    for name in ('compact_r17_wide_evidence_v1.json','new_rank26_followup_evidence_v1.json','exact_parity_coordinate_evidence_v1.json','small_conductor_followup_evidence_v1.json','fast_point_pipeline_evidence_v1.json'):
        path=ART/name;m=cert.read(path)
        if cert.hashed(ROOT/m['archive'])!=m['archive_sha256']:raise ArithmeticError('base archive changed')
        base_members.update({r['path']:r['sha256'] for r in m['files']});bases.append({'manifest':str(path.relative_to(ROOT)),'manifest_sha256':cert.hashed(path),'archive':m['archive'],'archive_sha256':m['archive_sha256']})
    paths={Path(__file__).resolve(),ROOT/'elliptic-curves/notes/EXTENDED_PRIME_RANK27_DISCOVERIES_2026-09-06.md'}
    for name in FOLDERS:paths.update(p for p in (LOCAL/name).rglob('*') if p.is_file())
    names=('r17_retained_prime_extension_diagnostics_v1.json','fresh_r17_paired_results_v1.json','new_paired_rank27_proof_v1.json','paired_high_rank_minimal_proofs_v2.json','new_high_rank_curve_index_v4.json','new_high_rank_curve_index_v4.csv','paired_rank27_adaptive_coverage_v1.json','paired_second27_adaptive_coverage_v1.json','paired_rank27_all_retained_mod2_v1.json','paired_rank27_all_retained_modl_v1.json','paired_second27_all_retained_mod2_v1.json','paired_second27_all_retained_modl_v1.json','new_paired_rank27_curve.sage','new_paired_rank27_curve_11952.sage','new_paired_rank26_curve.sage','fast_point_pipeline_portable_replay_v1.json')
    paths.update(ART/name for name in names);paths.update(ART.glob('fresh_r17_paired_*_mod2_v1.json'))
    scripts=('benchmark_r17_extended_prime_traces.py','extend_retained_r17_prime_scores.py','analyze_r17_prime_extension.py','fresh_r17_pari_batch.py','prepare_fresh_r17_pari_batch.sage','verify_fresh_r17_pari_batch.py','certify_fresh_r17_paired_results.py','certify_paired_rank27.py','certify_paired_high_rank_minimal.py','certify_paired_high_rank_minimal_v2.py','memory_rank_certificate.py','export_new_high_rank_curve_index_v4.py','audit_paired_rank27_conductor_bound.sage','followup_paired_rank27.py','prepare_paired_rank27_adaptive.sage','audit_paired_rank27_followup.py','followup_paired_second27.py','prepare_paired_second27_adaptive.sage','audit_paired_second27_followup.py','replay_paired_point_geometry.py','verify_paired_rank27_discovery_bundle.py')
    paths.update(ROOT/'elliptic-curves/cas'/name for name in scripts)
    paths.update(ROOT/'elliptic-curves/tests'/name for name in ('test_r17_prime_extension.py','test_memory_rank_certificate.py'))
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
    cert.write(OUTPUT,{'schema':'elliptic-curves.paired-rank27-discovery-evidence.v1','archive':str(archive.relative_to(ROOT)),'archive_sha256':cert.hashed(archive),'archive_bytes':archive.stat().st_size,'required_base_archives':bases,'files':members,'inherited_exact_members':inherited,'claim_boundary':'Extract the five pinned bases in order, then this supplement. Retains768 extended trace rosters,1037 initial charts,602 adaptive charts, all state archives and local replay records,47-curve inventory, two new minimal rank27 curves and a new minimal rank26 curve. The rejected four-new-curves assumption and corrected known390 label remain. Bounded misses do not establish rank upper bounds; universal novelty and new rank28/32 remain open.'})
    print('PACKAGED PAIRED RANK27 DISCOVERY',len(members),'new members',len(inherited),'inherited',archive.stat().st_size,'bytes',flush=True)
if __name__=='__main__':main()
