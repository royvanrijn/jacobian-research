#!/usr/bin/env python3
"""Seal the completed bounded result and allowlisted provenance; no discovery."""
import hashlib
import json
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
import zipfile

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/'artifacts/local/elliptic-curves/marked-two-class-v1'
GEN=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=GEN/'marked_two_class_propagation_v1'


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def write(p,x):
    with p.open('x') as f:json.dump(x,f,indent=2,sort_keys=True);f.write('\n')
def ev(cs,t):return sum((Q(c)*t**j for j,c in enumerate(cs)),Q(0))
def sqrt(q):
    q=Q(q);a=isqrt(q.numerator);b=isqrt(q.denominator)
    assert a*a==q.numerator and b*b==q.denominator
    return Q(a,b)


def main():
    OUT.mkdir(exist_ok=False)
    geometry=read(LOCAL/'geometry.json');family=geometry['family'];t0=Q(3,17)
    compact=read(GEN/'rank_jump_constructed_class_compaction_v1.json')
    atoms=read(GEN/'rank_jump_reference_additional_strict_verification_v1.json')
    ideals=read(GEN/'rank_jump_constructed_class_half_ideal_v1.json')
    generic=read(GEN/'constructed_class_blind_recovery_v1/generic.json')
    covers=read(GEN/'constructed_class_blind_recovery_v1/covers.json')
    original=read(GEN/'rank_jump_reference_strict_class_construction_inputs_v1.json')
    assert original['cubic_ascending']==covers['cubic_ascending']
    aa=ev(family['A_coefficients_low_to_high'],t0);bb=ev(family['B_coefficients_low_to_high'],t0)
    targetA,targetB=map(Q,generic['curve'][3:]);u2=(targetB/bb)/(targetA/aa);u=sqrt(u2)
    assert u**4*aa==targetA and u**6*bb==targetB
    def val(r):return ev(r['numerator_coefficients_low_to_high'],t0)/ev(r['denominator_coefficients_low_to_high'],t0)
    signs=[]
    for section,point,gamma in zip(family['sections'],generic['original_points'],original['generic_classes']):
        x,y=map(Q,point);xx=u2*val(section['X']);yy=u**3*val(section['Y'])
        assert xx==36*x+3 and yy**2==(216*y+108*x)**2
        signs.append(1 if yy==216*y+108*x else -1)
        beta=list(map(Q,gamma['beta_ascending']));scale=-beta[1]
        sqrt(scale);assert beta==[4*x*scale,-scale,Q(0)]
    source_names=['rank_jump_constructed_class_compaction_v1.json',
        'rank_jump_reference_additional_strict_verification_v1.json',
        'rank_jump_constructed_class_half_ideal_v1.json',
        'rank_jump_reference_strict_class_construction_inputs_v1.json',
        'constructed_class_blind_recovery_v1/covers.json',
        'constructed_class_blind_recovery_v1/generic.json']
    sources={}
    for name in source_names:
        p=GEN/name;dst=LOCAL/'block-provenance'/name;dst.parent.mkdir(parents=True,exist_ok=True)
        with dst.open('xb') as stream:stream.write(p.read_bytes())
        sources[str(p.relative_to(ROOT))]=digest(p)
    circuits=[]
    archive=GEN/'rank_jump_class_block_lift_evidence_v1.zip'
    with zipfile.ZipFile(archive) as z:
        for case in compact['cases']:
            index=case['column'];assert index in (6,7)
            source=case['square_equivalence_factorization'];data=z.read(source)
            assert hashlib.sha256(data).hexdigest()==case['square_equivalence_sha256']
            circuit=json.loads(data);assert circuit['column']==index
            assert ideals['columns'][index]['factor_labels']==atoms['class_representatives'][index-6]['factor_labels']
            count=sum(f['kind']=='projected_atom' for f in ideals['columns'][index]['factor_labels'])
            assert count==[1676,1572][index-6]
            target=LOCAL/'block-provenance'/f'column_{index:02d}_square_equivalence.json'
            with target.open('xb') as stream:stream.write(data)
            circuits.append({'column':index,'principal_atom_count':count,'source_member':source,
                'member_sha256':digest(target),'retained_path':str(target.relative_to(ROOT))})
            blind=next(r for r in covers['cases'] if r['column']==index)
            assert case['beta_ascending']==blind['beta_ascending']
            assert Q(blind['positive_norm_square_root'])**2==Q(blind['norm'])==Q(case['norm'])
    write(OUT/'frozen-block.json',{'schema':'marked-two-class.frozen-block.v1',
        'columns':[6,7],'control_parameter':'3/17','circuits':circuits,'source_bindings':sources,
        'circuit_archive_sha256':digest(archive),
        'generic_control_transport':{'short_scale_u':str(u),'formula':'Xshort=u^2*Xcompact=36*x_original+3; Yshort=u^3*Ycompact=216*y_original+108*x_original',
            'generic_ordinate_signs_relative_to_positive_control':signs,'generic_classes_checked':16},
        'boundary':'Retains full original atom dependencies and norm/square-equivalence circuits. Their earlier certificates are inherited, not recomputed. Geometry was already frozen independently; no carrier splits at the control, so field-square label matching is not reached. No historical exceptional points or oracle lifts are included.'})
    for name in ['protocol.json','geometry.json','frozen-candidates.json','independent-replay.json']:
        (OUT/name).write_bytes((LOCAL/name).read_bytes())
    results=read(LOCAL/'independent-replay.json');compiled=read(LOCAL/'compilation-v2.json')
    enum=read(LOCAL/'enumeration.json');proto=read(LOCAL/'protocol.json')
    write(OUT/'carriers.json',{'schema':'marked-two-class.carriers.v1','parameter':'compact MW16-05 t',
        'geometry_sha256':digest(OUT/'geometry.json'),
        'records':[read(LOCAL/'compiled-export-v2'/f'{i:03d}.json') for i in results['exact_smooth_geometrically_rational_bisections']]})
    measured={'marking':proto['preparation_cpu_seconds'],'enumeration':enum['cpu_seconds'],
        'compilation_children_including_prior_and_conservative_interrupt_charge':compiled['total_charged_child_cpu_seconds'],
        'export_repair_parent':compiled['parent_cpu_seconds'],'independent_replay':results['cpu_seconds']}
    write(OUT/'result.json',{'schema':'marked-two-class.result.v1','status':'PASS_BOUNDED_VISIBILITY_MISS',
        'theorem':'Valid sufficient simultaneous-splitting implication; explicit MW16-05 pair NOT constructed.',
        'enumeration_nodes':enum['nodes'],'shells_complete':False,'retained_shell_counts':enum['counts'],
        'candidate_count':256,'exact_bisections':182,'nonsplit_at_control':182,'reducible_residuals':73,
        'unknown_candidates':[171],'exact_nonzero_labels':0,'common_bases':0,'fresh_fibres':0,
        'strict_transfer':'UNKNOWN_NO_FRESH_FIBRES','ideal_transfer':'UNKNOWN_NO_FRESH_FIBRES',
        'nonsplitting_progression':results['nonsplitting_progression'],
        'recorded_or_charged_cpu_seconds':measured,'recorded_or_charged_sum':sum(measured.values()),
        'timing_boundary':'Includes the whole30-second allowance for an interrupted child. Process import/startup, the interrupted first parent overhead, authoring, source retrieval, packaging and navigation are not component-timed. No end-to-end speed comparison is claimed.',
        'independent_replay_sha256':digest(OUT/'independent-replay.json'),
        'tests':'18 small exact regressions passed; source included for replay.',
        'next_gate':'Two actual class-labelled carriers and an infinite rational common base. Any new search, shell continuation or interrupted-candidate work requires a new explicitly scoped protocol; no automatic enlargement.'})
    own_sources=[ROOT/'elliptic-curves/rank-jump'/n for n in [
        'prepare_marked_two_class.sage','enumerate_marked_two_class.py','compile_marked_two_class.sage',
        'finish_marked_two_class_export.sage','verify_marked_two_class.sage','test_marked_two_class.sage',
        'package_marked_two_class.py','test_shared_value_block.py','test_oblique_split_cubic.py',
        'retrospective.py','shared_value_soluble_block.py','shared_value_soluble_block_completion.py',
        'shared_value_halving.py','oblique_split_cubic.py']]
    controls=[GEN/n for n in ['rank_jump_shared_value_soluble_block_completion_v1.json','rank_jump_oblique_split_cubic_v1.json']]
    paths=sorted(set([p for p in LOCAL.rglob('*') if p.is_file() and '__pycache__' not in p.parts]+own_sources+controls+
        [p for p in OUT.iterdir() if p.is_file()]),key=lambda p:str(p.relative_to(ROOT)))
    records=[]
    with zipfile.ZipFile(OUT/'evidence.zip','x',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
        for p in paths:
            data=p.read_bytes();rel=str(p.relative_to(ROOT));info=zipfile.ZipInfo(rel,date_time=(1980,1,1,0,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED;z.writestr(info,data)
            records.append({'path':rel,'sha256':hashlib.sha256(data).hexdigest(),'bytes':len(data)})
    with zipfile.ZipFile(OUT/'evidence.zip') as z:
        assert len(z.namelist())==len(records)
        for r in records:assert hashlib.sha256(z.read(r['path'])).hexdigest()==r['sha256']
    write(OUT/'evidence-manifest.json',{'archive_sha256':digest(OUT/'evidence.zip'),'files':records,
        'restore':'Restore archive paths relative to research/ in an EMPTY replay checkout; never overwrite the live shared tree. Preserve the restored independent-replay.json before explicitly rerunning its exclusive-output checker.',
        'boundary':'Only allowlisted arithmetic/geometry, generic sections, frozen class provenance, failed receipts and this experiment are packaged. No old exceptional-point or V3-cloud payload is included.'})
    print(json.dumps({'archive_bytes':(OUT/'evidence.zip').stat().st_size,'members':len(records),
        'sha256':digest(OUT/'evidence.zip'),'recorded_or_charged_cpu_seconds':sum(measured.values())}),flush=True)


if __name__=='__main__':main()
