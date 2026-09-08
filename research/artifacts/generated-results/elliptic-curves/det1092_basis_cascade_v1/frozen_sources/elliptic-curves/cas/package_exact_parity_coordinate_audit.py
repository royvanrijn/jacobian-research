#!/usr/bin/env python3
"""Package the exact geometry and completed search-gap supplement over two pinned bases."""
from hashlib import sha256
from pathlib import Path
import zipfile
import certify_compact_r17_candidates as cert
from package_recorded_mod2_audit import dependencies
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves'
OUTPUT=ART/'exact_parity_coordinate_evidence_v1.json'
FOLDERS=('r17-norm12-exact-minima-v1','r17-norm12-exact-minima-remaining4-v1','r17-exact-parity-dataset-v1','r17-omitted-generic-classes-v1','new-rank26-fixed-tails-v1','native11952-metric49-control-v1','native11952-metric49-tails-v1','native11952-control-union-v1','native11952-pari-maps-v1','native11952-pari49-control-v1','new-rank26-pari-maps-v1','new-rank26-pari43-v1','native11952-engine-comparison-v1')

def main():
    archive=OUTPUT.with_suffix('.zip')
    if OUTPUT.exists() or archive.exists():raise FileExistsError('preserve supplement')
    for name in ('r17-omitted-generic-classes-v1',):
        d=cert.read(LOCAL/name/'verification/ledger.json')
        if len(d['rows'])!=8 or any(r['status']!='PASS' for r in d['rows']):raise ArithmeticError('omitted-class replay missing')
    for name in ('new-rank26-fixed-tails-v1','native11952-metric49-control-v1','native11952-metric49-tails-v1','new-rank26-pari43-v1'):
        d=cert.read(LOCAL/name/'verification/result.json')
        if len(d['rows'])!=3 or any(r['status']!='PASS' for r in d['rows']):raise ArithmeticError('chart/cloud replay missing')
    v=cert.read(LOCAL/'native11952-pari49-control-v1/verification-v2/result.json')
    if v['outcome']!='completed' or v['returncode']!=0:raise ArithmeticError('v2 map replay missing')
    v=cert.read(LOCAL/'native11952-pari49-control-v1/verification/result.json')
    if any(r['status']!='PASS' for r in v['rows'] if r['stage']!='charts'):raise ArithmeticError('PARI map cloud replay missing')
    base_members={};bases=[]
    for name in ('compact_r17_wide_evidence_v1.json','new_rank26_followup_evidence_v1.json'):
        path=ART/name;m=cert.read(path);zpath=ROOT/m['archive']
        if cert.hashed(zpath)!=m['archive_sha256']:raise ArithmeticError('base archive changed')
        base_members.update({r['path']:r['sha256'] for r in m['files']});bases.append({'manifest':str(path.relative_to(ROOT)),'manifest_sha256':cert.hashed(path),'archive':m['archive'],'archive_sha256':m['archive_sha256']})
    paths={Path(__file__).resolve(),ROOT/'elliptic-curves/cas/verify_exact_parity_coordinate_bundle.py',ROOT/'elliptic-curves/notes/EXACT_PARITY_AND_COORDINATE_AUDIT_2026-09-06.md',ROOT/'elliptic-curves/tests/test_exact_parity_ellipsoid.py'}
    for name in FOLDERS:paths.update(p for p in (LOCAL/name).rglob('*') if p.is_file())
    for pattern in ('r17_exact_parity_radius_*_v1.json','r17_omitted_classes_recorded_mod2_*_v1.json','native11952_*_v1.json'):
        paths.update(ART.glob(pattern))
    paths.update(ART/name for name in ('r17_exact_maximum_parity_classes_v1.json','new_rank26_tail_coverage_v1.json','new_rank26_tail_coverage_v2.json','new_rank26_all_retained_mod2_v1.json','new_rank26_pari43_recorded_mod2_v1.json'))
    paths.update(ROOT/'elliptic-curves/cas'/name for name in ('export_r17_exact_maximum_parity_classes.py','audit_r17_norm12_minima.sage','audit_r17_norm12_minima_remaining4.sage','replay_r17_norm12_minima.py','search_r17_omitted_classes.sage','run_r17_omitted_classes.py','replay_r17_omitted_classes.py','verify_r17_omitted_classes.py','search_new_rank26_tails.py','replay_new_rank26_tails.py','replay_new_rank26_tails_v2.py','verify_new_rank26_tails.py','calibrate_native11952_metric49.sage','replay_native11952_metric49.py','verify_native11952_metric49.py','search_native11952_tails.py','replay_native11952_tails.py','verify_native11952_tails.py','audit_native11952_visibility.py','audit_native11952_control_union.py','prepare_native11952_pari_maps.sage','search_native11952_pari49.py','replay_native11952_pari49.py','verify_native11952_pari49.py','replay_native11952_pari49_v2.py','prepare_new_rank26_pari_maps.sage','search_new_rank26_pari43.py','replay_new_rank26_pari43.py','verify_new_rank26_pari43.py','compare_native11952_search_engines.py'))
    seen=set()
    while True:
        paths=dependencies(paths);todo=[p for p in paths-seen if p.suffix=='.json']
        if not todo:break
        for p in todo:
            seen.add(p);data=cert.read(p)
            if not isinstance(data,dict):continue
            for key in ('sources','source_hashes','inputs'):
                mapping=data.get(key,{})
                if isinstance(mapping,dict):
                    for name in mapping:
                        q=ROOT/name
                        if q.is_file():paths.add(q)
            for key in ('input_path','parent_input_path','generic_census_path'):
                if isinstance(data.get(key),str) and (ROOT/data[key]).is_file():paths.add(ROOT/data[key])
    members=[];inherited=[]
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for p in sorted(paths):
            raw=p.read_bytes();name=str(p.relative_to(ROOT));h=sha256(raw).hexdigest();row={'path':name,'bytes':len(raw),'sha256':h}
            if base_members.get(name)==h:inherited.append(row);continue
            info=zipfile.ZipInfo(name,date_time=(2026,9,6,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o644<<16;z.writestr(info,raw);members.append(row)
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None or any(sha256(z.read(r['path'])).hexdigest()!=r['sha256'] for r in members):raise ArithmeticError('supplement integrity failed')
    cert.write(OUTPUT,{'schema':'elliptic-curves.exact-parity-coordinate-evidence.v1','archive':str(archive.relative_to(ROOT)),'archive_sha256':cert.hashed(archive),'archive_bytes':archive.stat().st_size,'required_base_archives':bases,'files':members,'inherited_exact_members':inherited,'claim_boundary':'Extract the two pinned base archives in listed order, then this supplement. Inherited member bytes are hash-identical; changed sources and notes are included explicitly. Contains exact lattice proofs, completed finite search-gap records, versioned coverage correction and known-curve coordinate controls. Does not certify a new rank28/32 curve.'})
    print('PACKAGED PARITY/COORDINATE SUPPLEMENT',len(members),'new members',len(inherited),'inherited',archive.stat().st_size,'bytes',flush=True)
if __name__=='__main__':main()
