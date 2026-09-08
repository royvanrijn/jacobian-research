#!/usr/bin/env python3
"""Portable height calibration and saved-shard discovery evidence supplement."""
from pathlib import Path
from hashlib import sha256
import zipfile
import certify_compact_r17_candidates as cert
from package_recorded_mod2_audit import dependencies
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';OUT=ART/'height_and_discarded_evidence_v1.json'
FOLDERS=('native11952-translated-visibility-v1','native11952-height-pair-v1','native11952-height125-control-v1','new27-height125-followup-v1','r17-discarded-shards-extended-v1','discarded12-r17-pari-v1')
def main():
    archive=OUT.with_suffix('.zip')
    if OUT.exists() or archive.exists():raise FileExistsError('preserve height/discarded evidence')
    if cert.read(LOCAL/'discarded12-r17-pari-v1/verification-ledger.json')['status']!='PASS' or cert.read(ART/'new27_height125_followup_coverage_v1.json')['status']!='PASS' or cert.read(ART/'native29_height_control_v2.json')['status']!='PASS':raise ArithmeticError('required finite checks incomplete')
    prior=ART/'next24_discovery_evidence_v1.json';pm=cert.read(prior);base_names=[r['manifest'] for r in pm['required_base_archives']]+[str(prior.relative_to(ROOT))];bases=[];base_members={}
    for name in base_names:
        p=ROOT/name;m=cert.read(p)
        if cert.hashed(ROOT/m['archive'])!=m['archive_sha256']:raise ArithmeticError('base archive changed')
        base_members.update({r['path']:r['sha256'] for r in m['files']});bases.append({'manifest':name,'manifest_sha256':cert.hashed(p),'archive':m['archive'],'archive_sha256':m['archive_sha256']})
    paths={Path(__file__).resolve(),ROOT/'elliptic-curves/notes/TRANSLATED_HEIGHT_AND_DISCARDED_SHARDS_2026-09-06.md'}
    for name in FOLDERS:paths.update(p for p in (LOCAL/name).rglob('*') if p.is_file())
    names=('native11952_translated_visibility_v1.json','native11952_translated_visibility_replay_v1.json','native29_height_control_v1.json','native29_height_control_v2.json','native11952_height100000_mod2_v1.json','native11952_height125000_mod2_v1.json','native11952_height1000000_mod2_v1.json','new27_height125_followup_coverage_v1.json','discarded12_r17_results_v1.json','discarded_rank26_minimal_proof_v1.json','new_high_rank_curve_index_v6.json','new_high_rank_curve_index_v6.csv','new_high_rank_curve_index_v6_memory_replay_v1.json','new_discarded_rank26_curve.sage','next24_portable_replay_v1.json')
    paths.update(ART/n for n in names);paths.update(ART.glob('discarded12_r17_*_mod2_v1.json'));paths.update(ART.glob('new27_height125_*_mod*_v1.json'))
    names=('audit_native11952_translated_visibility.sage','replay_native11952_translated_visibility.py','replay_native11952_translated_visibility_v2.py','native11952_height_pair.py','native11952_height125_control.py','verify_native11952_height_controls.py','certify_native29_height_control.py','certify_native29_height_control_v2.py','new27_height125_followup.py','audit_new27_height125_followup.py','extend_r17_discarded_shard_scores.py','discarded12_r17_pari_batch.py','prepare_discarded12_r17_pari_batch.sage','verify_discarded12_r17_pari_batch.py','certify_discarded12_r17_results.py','certify_discarded_rank26_minimal.py','export_new_high_rank_curve_index_v6.py','replay_inventory_v6_memory.py','replay_height_and_discarded_geometry.py','verify_height_and_discarded_bundle.py')
    paths.update(ROOT/'elliptic-curves/cas'/n for n in names);seen=set()
    while True:
        paths=dependencies(paths);todo=[p for p in paths-seen if p.suffix=='.json']
        if not todo:break
        for p in todo:
            seen.add(p);d=cert.read(p)
            if not isinstance(d,dict):continue
            for key in ('sources','source_hashes','inputs','bindings','shard_hashes'):
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
    cert.write(OUT,{'schema':'elliptic-curves.height-and-discarded-evidence.v1','archive':str(archive.relative_to(ROOT)),'archive_sha256':cert.hashed(archive),'archive_bytes':archive.stat().st_size,'required_base_archives':bases,'files':members,'inherited_exact_members':inherited,'claim_boundary':'Extract nine pinned bases in order then this supplement. Retains196 exact retrospective translations,147 matched/capped control charts,141 new27 follow-up charts,768 saved-address trace rosters,540 fresh generic boxes and all archived histories,70-curve inventory and a new globally minimal rank26 curve. The49 million-height timeouts remain censored. No exact rank, universal novelty or new rank28/32 curve.'})
    print('PACKAGED HEIGHT/DISCARDED',len(members),'new',len(inherited),'inherited',archive.stat().st_size,'bytes',flush=True)
if __name__=='__main__':main()
