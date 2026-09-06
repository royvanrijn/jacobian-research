#!/usr/bin/env python3
"""Portable wider-retention discovery and exact control-diagnostic evidence."""
from pathlib import Path
from hashlib import sha256
import zipfile
import certify_compact_r17_candidates as cert
from package_recorded_mod2_audit import dependencies
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';OUT=ART/'million_height_mw16_evidence_v1.json'
FOLDERS=('native29-million-chart-benchmark-v1','native29-cvp-neighbours-v1','pointed-sieve-v2-benchmark-v1','mw16-top25-pari-followup-v1','mw16-extended-prime-benchmark-v1','mw16-retained-extended-primes-v1','new27-million-height-pilot-v1','extended20-mw16-pari-v1','extended20-mw16-point-unions-v1','mw16-exact-maximum-parities-v1','mw16-exact-maximum-parities-v2','extended20-mw16-incidence-v1','million-mw16-diagnostic-replay-v1')

def main():
    archive=OUT.with_suffix('.zip')
    if OUT.exists() or archive.exists():raise FileExistsError('preserve height/discarded evidence')
    if cert.read(LOCAL/'extended20-mw16-point-unions-v1/ledger.json')['status']!='PASS' or cert.read(ART/'new_high_rank_curve_index_v8_memory_replay_v1.json')['status']!='PASS' or cert.read(LOCAL/'mw16-exact-maximum-parities-v2/ledger.json')['status']!='PASS':raise ArithmeticError('completed exact gates required')
    prior=ART/'retention_rank27_followup_evidence_v1.json';pm=cert.read(prior);base_names=[r['manifest'] for r in pm['required_base_archives']]+[str(prior.relative_to(ROOT))];bases=[];base_members={}
    for name in base_names:
        p=ROOT/name;m=cert.read(p)
        if cert.hashed(ROOT/m['archive'])!=m['archive_sha256']:raise ArithmeticError('base archive changed')
        base_members.update({r['path']:r['sha256'] for r in m['files']});bases.append({'manifest':name,'manifest_sha256':cert.hashed(p),'archive':m['archive'],'archive_sha256':m['archive_sha256']})
    paths={Path(__file__).resolve(),ROOT/'elliptic-curves/notes/MILLION_HEIGHT_AND_MW16_EXTENSION_2026-09-06.md'}
    for name in FOLDERS:paths.update(p for p in (LOCAL/name).rglob('*') if p.is_file())
    for pattern in ('native29_million_chart_benchmark_v1.json','native29_cvp_neighbours*_v1.json','pointed_sieve_v2_benchmark_v1.json','mw16_top25_*_v1.json','new27_million_*_v1.json','extended20_mw16*_v*.json','mw16_exact_maximum_parities_*_v2.json','new_high_rank_curve_index_v8*','new_mw16_rank26_*.sage','inventory98_cross_family_incidence_v1.json'):
        paths.update(ART.glob(pattern))
    names=('benchmark_native29_million_chart.py','audit_native29_cvp_neighbours.py','benchmark_pointed_sieve_v2.py','pointed_quartic_sieve_v2.cpp','mw16_top25_pari_followup.py','prepare_mw16_top25_pari.sage','benchmark_mw16_extended_prime_traces.py','extend_mw16_retained_prime_scores.py','new27_million_height_pilot.py','audit_mw16_and_million_followups.py','replay_million_and_mw16_diagnostics.py','extended20_mw16_pari_batch.py','prepare_extended20_mw16_pari_batch.sage','verify_extended20_mw16_pari_batch.py','certify_extended20_mw16_results.py','replay_extended20_mw16_geometry.py','combine_extended20_mw16_retained_points.py','certify_extended20_mw16_minimal.py','export_new_high_rank_curve_index_v8.py','replay_inventory_v8_memory.py','export_extended20_mw16_sage.py','audit_mw16_exact_maximum_parities.sage','audit_mw16_exact_maximum_parities_v2.sage','replay_mw16_exact_maximum_parities.py','audit_extended20_mw16_incidence.sage','replay_extended20_mw16_incidence.py','replay_extended20_mw16_incidence_v2.py','certify_inventory98_incidence.py','verify_million_height_mw16_bundle.py')
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
    cert.write(OUT,{'schema':'elliptic-curves.million-height-mw16-evidence.v1','archive':str(archive.relative_to(ROOT)),'archive_sha256':cert.hashed(archive),'archive_bytes':archive.stat().st_size,'required_base_archives':bases,'files':members,'inherited_exact_members':inherited,'claim_boundary':'Extract twelve pinned base archives in order and then this supplement. Retains the native29 retrospective control,5488 exact translates, unpromoted GMP prototype, MW16 selector and20 point attempts, complete point unions, three new minimal26-point curves,98-curve inventory, parity and incidence certificates, and the completed new27 million-height pilot. Failures and corrections remain. No exact rank, new28/32, conductor record or universal novelty.'})
    print('PACKAGED MILLION HEIGHT AND MW16',len(members),'new',len(inherited),'inherited',archive.stat().st_size,'bytes',flush=True)
if __name__=='__main__':main()
