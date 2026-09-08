#!/usr/bin/env python3
"""Bind the degree-three construction attempt and its independent proof."""
import argparse,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves'
OUT=ART/'mestre_degree_three_gate_v1.json'
def read(p):return json.loads(p.read_text())
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def compute():
    stages=[];paths=[Path(__file__).resolve()]
    jobs=[('bounded-triangles','mestre-degree-three-pencils-v1','roster.supervisor.json'),
          ('height-domain','mestre-all-section-triangles-v1','enumeration.supervisor.json'),
          ('alternate-section-admission','mestre-rational-bisections-v1','admission.supervisor.json'),
          ('rational-curves','mestre-rational-bisections-v1','curves.supervisor.json'),
          ('basis-translations','mestre-translated-bisections-v1','supervisor.json'),
          ('class-selection','mestre-bisection-pencil-admission-v1','supervisor.json'),
          ('complete-frames','mestre-degree-three-frames-v1','supervisor.json'),
          ('independent-proof','mestre-degree-three-standalone-v1','supervisor.json')]
    for name,folder,supervisor in jobs:
        d=LOCAL/folder;s=read(d/supervisor);assert s['outcome']=='completed' and s['returncode']==0 and s['failure_reason'] is None
        assert digest(Path(s['log']))==s['log_sha256'];p=read(d/'protocol.json')
        for filename,h in p.get('sources',{}).items():assert digest(ROOT/filename)==h;paths.append(ROOT/filename)
        stages.append({'name':name,'wall_seconds':s['wall_seconds'],'supervisor':str((d/supervisor).relative_to(ROOT))})
        paths.extend([d/'protocol.json',d/supervisor,Path(s['log'])])
    bundle=ART/'mestre_degree_three_proof_bundle_v1.json';b=read(bundle)
    for filename,h in b['sources'].items():assert digest(ROOT/filename)==h;paths.append(ROOT/filename)
    verifier=ROOT/'elliptic-curves/cas/verify_mestre_degree_three_gate.sage';fresh=LOCAL/'mestre-degree-three-standalone-v1'
    for p in (bundle,verifier):assert p.read_bytes()==(fresh/p.name).read_bytes();paths.append(p)
    transcript=(fresh/'replay.log').read_text();assert transcript.strip()=='PASS standalone exact bisections, primitive degree3 pencils and exhaustive roots [[8, 9, 20, 776], [10, 11, 14, 486], [10, 11, 14, 578]]'
    bounded=read(LOCAL/'mestre-degree-three-pencils-v1/roster.json');admission=read(LOCAL/'mestre-bisection-pencil-admission-v1/result.json')
    assert bounded['eligible_pencils']==1726 and len(admission['rows'])==554
    curves=read(LOCAL/'mestre-rational-bisections-v1/result.json');translations=read(LOCAL/'mestre-translated-bisections-v1/result.json')
    assert len(curves['rows'])==50 and len(translations['rows'])==46
    assert [r['index'] for r in curves['rows'] if r.get('old_base_degree')==2]==[27,35]
    rows=[{k:r[k] for k in ['source_index','subtract_section_word','fibre_class','visible_MW_upper_bound','generic_Q_MW_rank','generic_geometric_MW_rank','geometric_root_count','geometric_root_rank','rational_root_rank']} for r in b['frames']['rows']]
    return {'schema':'mestre-degree-three-gate.v1','status':'PASS','rows':rows,
        'bounded_three_section_pencils':1726,'alternate_section_words':50,'rational_bisections':2,'basis_translations':46,'class_selected_pencils':554,
        'rational_MW_lower_upper_equal':[8,10,10],'geometric_MW_lower_upper_equal':[9,11,11],
        'point_search_boxes':0,'new_parents':0,'inventory_additions':0,'stages':stages,'total_supervised_seconds':sum(s['wall_seconds'] for s in stages),
        'standalone_transcript':transcript,'sources':{str(p.relative_to(ROOT)):digest(p) for p in paths},
        'scope':'Three exact primitive old-degree3 rational Jacobian pencils on the explicit u11 Mestre K3, constructed as O plus translates of two rational bisections. Independent rational LDL enumeration closes all geometric frame roots and gives exact generic Q ranks8,10,10. No pairwise inequivalence, new Weierstrass equations or point bases, high-rank fibre, search-yield estimate or parent expansion is claimed. Visible rank14 upper bounds were loose by four or six directions. The retained height domain is selection input; this report does not certify completeness of that separate height-six enumeration or all degree3 fibrations. Two initial launcher API errors preceded the first CAS worker; no arithmetic rerun occurred. Supervised stage accounting excludes launcher, source-review and report CPU.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();r=compute()
    if a.check:assert r==read(OUT)
    else:
        with OUT.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print('PASS degree3 exact generic Q ranks8,10,10;',r['total_supervised_seconds'],'supervised seconds')
