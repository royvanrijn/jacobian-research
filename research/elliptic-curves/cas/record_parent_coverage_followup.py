#!/usr/bin/env python3
"""Bind completed retained-frame and Kihara section proofs, including failures."""
import argparse,ast,hashlib,json
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def collect(jobs,artifacts):
    paths={Path(__file__).resolve()};stages=[]
    def bind(obj):
        for name,h in obj.get('sources',{}).items():
            p=ROOT/name;assert sha(p)==h,(name,'source changed');paths.add(p)
    for folder,expected in jobs:
        d=LOCAL/folder;s=read(d/'supervisor.json');p=read(d/'protocol.json')
        assert s['outcome']==expected and (s['returncode']==0)==(expected=='completed')
        log=Path(s['log']);assert sha(log)==s['log_sha256'];bind(p)
        for name,h in p.get('files',{}).items():
            f=d/name;assert sha(f)==h;paths.add(f)
            canonical=(ART if name.endswith('.json') else ROOT/'elliptic-curves/cas')/name
            assert canonical.read_bytes()==f.read_bytes();paths.add(canonical)
        paths.update(d.glob('*.json'));paths.add(log)
        stages.append({'stage':folder,'outcome':expected,'wall_seconds':s['wall_seconds'],'supervisor':str((d/'supervisor.json').relative_to(ROOT))})
    for name in artifacts:
        p=ART/name;paths.add(p);bind(read(p))
    return stages,paths
def main(which):
    if which=='mestre':
        artifacts=['mestre_degree_three_gate_v1.json','mestre_degree_three_proof_bundle_v1.json','mestre_retained_pencil_proof_bundle_v1.json','mestre_retained_visible_proof_bundle_v1.json']
        jobs=[('mestre-retained-pencil-frames-v1','completed'),('mestre-retained-pencil-standalone-v1','backend_failure'),('mestre-retained-pencil-standalone-v2','completed'),('mestre-retained-visible-standalone-v1','completed')]
        stages,paths=collect(jobs,artifacts)
        b=read(ART/artifacts[2]);v=read(ART/artifacts[3]);prior=read(ART/artifacts[1]);rows=b['frames']['rows']
        assert rows==read(LOCAL/jobs[0][0]/'result.json')['rows']==v['new_frames']
        assert v['prior_frames']==prior['frames']['rows'] and len(rows)==200 and len(v['admission_rows'])==554
        prefix='PASS standalone exact bisections, primitive degree3 pencils and exhaustive roots '
        transcript=(LOCAL/jobs[2][0]/'replay.log').read_text().strip();assert transcript.startswith(prefix)
        replay=ast.literal_eval(transcript[len(prefix):]);assert len(replay)==200
        assert [r[0] for r in replay]==[r['generic_Q_MW_rank'] for r in rows]
        assert [r[1] for r in replay]==[r['generic_geometric_MW_rank'] for r in rows]
        visible=(LOCAL/jobs[3][0]/'replay.log').read_text().strip()
        assert visible=='PASS554 distinct fibre classes;203 exact-frame inputs;351 visible ceilings; all MW_Q<=11'
        histogram=dict(sorted(Counter(r['generic_Q_MW_rank'] for r in rows).items()));assert histogram=={6:5,7:28,8:69,9:71,10:27}
        result={'new_exact_frame_count':200,'new_exact_Q_rank_histogram':{str(k):v for k,v in histogram.items()},'previous_exact_frame_count':3,'retained_distinct_fibre_classes':554,'visible_only_count':351,'all_retained_generic_Q_MW_upper_bound':11,
            'scope':'All554 retained Mestre u11 O-plus-translated-bisection pencils have generic rational MW rank at most11. Exactly200 additional full geometric frames were independently checked; three earlier proofs are reused and351 lower visible ceilings independently recomputed. Distinct divisor classes need not be inequivalent fibrations. No all-degree3 classification, height-domain completeness, specialized-rank or yield bound. The first standalone attempt failed an overrestrictive translated-class uniqueness assertion; V2 uses exact component rotation without rerunning frame production. Failed cost included.'}
    else:
        artifacts=['kihara_quadratic_section_bundle_v1.json','kihara_quadratic_section_full_NS_v1.json','kihara_universal_anti_bundle_v1.json','kihara_universal_anti_replay_v1.json']
        jobs=[('kihara-quotient-characters-v1','completed'),('kihara-quotient-minimal-sections-v1','strict_wall_timeout'),('kihara-quotient-section-extraction-v1','completed'),('kihara-quadratic-section-standalone-v1','completed'),('kihara-universal-anti-v1','backend_failure'),('kihara-universal-anti-extraction-v1','completed'),('kihara-universal-anti-standalone-v1','completed')]
        stages,paths=collect(jobs,artifacts)
        full=read(ART/artifacts[1]);universal=read(ART/artifacts[3]);assert full['status']==universal['status']=='PASS'
        assert full['geometric_NS_determinant']==-756 and full['constant_field']=='Q(sqrt(-3))'
        for folder,prefix in [(jobs[3][0],'PASS first Kihara full geometric MW13'),(jobs[6][0],'PASS unrestricted Kihara quadratic section')]:assert (LOCAL/folder/'replay.log').read_text().startswith(prefix)
        result={'first_parent_path_parameter':'3/2','full_geometric_NS_rank':18,'full_geometric_NS_determinant':-756,'full_geometric_NS_discriminant_group':[3,6,42],'exact_first_parent_constant_field':'Q(sqrt(-3))','first_parent_rational_MW_rank':12,'first_parent_geometric_MW_rank':13,'unrestricted_ratio_anti_section':'Q(v)(sqrt(-3))(T)',
            'path_coordinate_coverage':universal['path_coverage'],
            'scope':'First Kihara parent full geometric basis and field proved using direct section, Galois trace and height identities plus prior full rational basis and exact ranks. Separate unrestricted-ratio anti-section identity does not determine every specialized full lattice, field or rank. The path covers only -1/2<v<0; omitted coordinate regions are not automatically distinct moduli. Five retained quotient character audits are diagnostic, not field proofs. Preserved FGLM timeout and overly linear universal-extraction assumption are charged; saved bases are reused. No point search or new parent values.'}
    result.update({'schema':which+'-parent-coverage-followup.v1','status':'PASS','stages':stages,'total_supervised_seconds':sum(s['wall_seconds'] for s in stages),'new_parents':0,'point_search_boxes':0,'inventory_additions':0,'sources':{str(p.relative_to(ROOT)):sha(p) for p in sorted(paths)}})
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('which',choices=['mestre','kihara']);p.add_argument('--check',action='store_true');a=p.parse_args();result=main(a.which);out=ART/(a.which+'_parent_coverage_followup_v1.json')
    if a.check:assert result==read(out)
    else:
        with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print('PASS',a.which,'parent coverage follow-up;',result['total_supervised_seconds'],'supervised seconds including failed stages')
