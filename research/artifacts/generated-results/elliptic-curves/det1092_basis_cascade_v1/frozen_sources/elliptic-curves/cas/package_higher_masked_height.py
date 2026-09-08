#!/usr/bin/env python3
"""Portable wider-retention discovery and exact control-diagnostic evidence."""
from pathlib import Path
from hashlib import sha256
import zipfile
import certify_compact_r17_candidates as cert
from package_recorded_mod2_audit import dependencies
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';OUT=ART/'higher_masked_height_evidence_v1.json'
FOLDERS=('higher-r17-masked-controls-v1','higher26-million-height-pilot-v1')

def main():
    archive=OUT.with_suffix('.zip')
    if OUT.exists() or archive.exists():raise FileExistsError('preserve height/discarded evidence')
    if cert.read(ART/'higher26_million_higher_11952_coverage_v1.json')['status']!='PASS' or any(r['recovery_status']!='RECOVERED_KNOWN_DIRECTION' for r in cert.read(ART/'higher_r17_masked_relations_v1.json')['rows']):raise ArithmeticError('exact new26 and inventory gates required')
    prior=ART/'product22_comparison_evidence_v1.json';pm=cert.read(prior);base_names=[r['manifest'] for r in pm['required_base_archives']]+[str(prior.relative_to(ROOT))];bases=[];base_members={}
    for name in base_names:
        p=ROOT/name;m=cert.read(p)
        if cert.hashed(ROOT/m['archive'])!=m['archive_sha256']:raise ArithmeticError('base archive changed')
        base_members.update({r['path']:r['sha256'] for r in m['files']});bases.append({'manifest':name,'manifest_sha256':cert.hashed(p),'archive':m['archive'],'archive_sha256':m['archive_sha256']})
    paths={Path(__file__).resolve(),ROOT/'elliptic-curves/notes/HIGHER_R17_MASKED_CONTROLS_2026-09-06.md',ROOT/'elliptic-curves/notes/HIGHER26_MILLION_HEIGHT_2026-09-06.md',LOCAL/'higher-rank26-11952-069-newdirection-v2/all-retained-point-cloud-only.json'}
    for name in FOLDERS:paths.update(p for p in (LOCAL/name).rglob('*') if p.is_file())
    for pattern in ('higher26_million_*_v1.json','higher_r17_masked_relations_v1.json'):paths.update(ART.glob(pattern))
    names=('higher_r17_masked_controls.py','prepare_higher_r17_masked_controls.sage','audit_higher_r17_masked_relations.sage','higher26_million_height_pilot.py','audit_higher26_million_height.py','verify_higher_masked_height_bundle.py')
    paths.update(ROOT/'elliptic-curves/cas'/n for n in names);seen=set()
    while True:
        paths=dependencies(paths);todo=[p for p in paths-seen if p.suffix=='.json']
        if not todo:break
        for p in todo:
            seen.add(p);d=cert.read(p)
            if not isinstance(d,dict):continue
            if p.name=='elkies-k3-r17-norm12-icarm-calibration-dataset-v1.json':continue # Retrospective metadata only; this audit does not reprove its historical ranks or geometry.
            for key in ('sources','source_hashes','inputs','bindings','shard_hashes','canonical_table_hashes','old_scored_hashes','height_control_bindings','selection_proof_hashes'):
                value=d.get(key,{})
                if isinstance(value,dict):paths.update(ROOT/name for name in value if (ROOT/name).is_file())
            if isinstance(d.get('table'),dict) and (ROOT/d['table'].get('path','')).is_file():paths.add(ROOT/d['table']['path'])
            for name in d.get('source_certificate_hashes',{}):
                if (ART/name).is_file():paths.add(ART/name)
            for key in ('input_path','seed_path','maps_path','generic_census_path','supersedes_failed_protocol'):
                if isinstance(d.get(key),str) and (ROOT/d[key]).is_file():paths.add(ROOT/d[key])
    for row in cert.read(LOCAL/'higher26-million-height-pilot-v1/protocol.json')['rows']:
        paths.update(ROOT/row[k] for k in ('seed_path','maps_path','proof_path'))
    members=[];inherited=[]
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for p in sorted(paths):
            raw=p.read_bytes();name=str(p.relative_to(ROOT));h=sha256(raw).hexdigest();row={'path':name,'sha256':h,'bytes':len(raw)}
            if base_members.get(name)==h:inherited.append(row);continue
            info=zipfile.ZipInfo(name,date_time=(2026,9,6,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o644<<16;z.writestr(info,raw);members.append(row)
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None or any(sha256(z.read(r['path'])).hexdigest()!=r['sha256'] for r in members):raise ArithmeticError('archive integrity differs')
    cert.write(OUT,{'schema':'elliptic-curves.higher-masked-height-evidence.v1','archive':str(archive.relative_to(ROOT)),'archive_sha256':cert.hashed(archive),'archive_bytes':archive.stat().st_size,'required_base_archives':bases,'files':members,'inherited_exact_members':inherited,'claim_boundary':'Extract seventeen pinned base archives and this supplement. Retains all72 masked-control boxes, exact recovered-known-direction relations, and the fixed49-map million-height pilot on the new26 curve with its complete prior and current point-cloud union. Exact finite quotient certificates are lower bounds. These controls do not calibrate exceptional-direction sensitivity or establish upper ranks, saturation or universal novelty.'})
    print('PACKAGED HIGHER MASKED/HEIGHT',len(members),'new',len(inherited),'inherited',archive.stat().st_size,'bytes',flush=True)
if __name__=='__main__':main()
