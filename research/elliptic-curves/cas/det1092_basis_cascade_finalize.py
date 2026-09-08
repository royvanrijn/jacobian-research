#!/usr/bin/env python3
"""Seal the completed experiment and render its terminal result table."""
import hashlib,json,datetime,subprocess,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];PKG=ROOT/'artifacts/generated-results/elliptic-curves/det1092_basis_cascade_v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
s=read(PKG/'summary.json');replay=read(PKG/'compact_replay.json');groups=read(PKG/'group_diagnostic.json');roster=read(PKG/'roster.json');protocol=read(PKG/'protocol.json');ledger=read(PKG/'presentation_ledger.json');commit=s['starting_commit']
assert replay['status']=='PASS_PORTABLE_EXACT_REPLAY' and groups['status']=='COMPLETE_POST_RUN_EXACT_DIAGNOSTIC'
assert replay['summary_sha256']==groups['summary_sha256']==sha(PKG/'summary.json')
assert protocol['roster_sha256']==sha(PKG/'roster.json') and len(s['arms'])==len(roster['arms'])==4
for row in s['arms']:
 assert row['status']=='TERMINAL' and row['replay_status']=='PASS_EXACT_REPLAY'
 assert sha(PKG/f"evidence-{row['arm_id']}.json.gz")==row['evidence_sha256']
 terminal=read(ROOT/'artifacts/local/elliptic-curves/det1092-basis-cascade-v1'/row['arm_id']/'run/terminal.json')
 allowed='artifacts/local/elliptic-curves/det1092-basis-cascade-v1/'+row['arm_id']+'/'
 assert all(p.startswith(allowed) or p.startswith('artifacts/local/elliptic-curves/pointed-sieve-build/') for p in terminal.get('read_paths',[]))
 assert next(x for x in replay['arms'] if x['arm_id']==row['arm_id'])['rank_lower_bound']==row['final_rank_lower_bound']
for path,h in protocol['all_source_bindings'].items():assert sha(PKG/'frozen_sources'/path)==h
ranks=[r['final_rank_lower_bound'] for r in s['arms']]
classification='FOUR_BASES_REACH31_WITH_DIFFERENT_PATHS' if all(r>=31 for r in ranks) else 'BASIS_DEPENDENT_BOUNDED_CHAIN_RECOVERY' if any(r>=31 for r in ranks) else 'BOUNDED_CHAIN_MISS'
post_sources={}
for p in (ROOT/'elliptic-curves/cas').glob('det1092_basis_cascade*'):
 if p.is_file():
  shutil.copyfile(p,PKG/'frozen_orchestration'/p.name);post_sources[str(p.relative_to(ROOT))]=sha(p)
files=[p for p in PKG.iterdir() if p.is_file() and p.suffix in ('.json','.gz','.svg','.png')]
result={'starting_commit':commit,'completion_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'COMPLETE_INDEPENDENTLY_REPLAYED_EXPERIMENT','classification':classification,'artifact_bindings':{p.name:sha(p) for p in files},'postprocessing_and_orchestration_sources':post_sources,'primary_arms':4,'primary_launches':4,'primary_reruns':0,'recorded_worker_artifact_access_stayed_within_own_arm_or_build_cache':True,'independently_replayed_chart_records':sum(r['charts'] for r in s['arms']),'new_elliptic_curve_isomorphism_classes':ledger['new_elliptic_curve_isomorphism_classes'],'maximum_post_run_union_rank_lower_bound':groups['maximum_union_rank_lower_bound'],'above31_review_required':groups['above31_independent_review_required'],'finalizer_sha256':sha(Path(__file__)),'claim_boundary':s['claim_boundary']}
with (PKG/'completion.json').open('x') as f:json.dump(result,f,sort_keys=True,indent=2);f.write('\n')
note=ROOT/'elliptic-curves/notes/DET1092_BASIS_CASCADE_2026-09-08.md';text=note.read_text();text=text.replace('**Experiment in progress; no terminal conclusion yet.**',f"**Complete: `{classification}`.**")
section='\n## Terminal results\n\n![Certified rank versus completed charts](../../artifacts/generated-results/elliptic-curves/det1092_basis_cascade_v1/rank_path_comparison.svg)\n\nAll ranks below are independently certified lower bounds.\n\n| Basis | Final rank | Charts | Worker wall (min) | Total CPU (min) | Final subgroup versus original |\n|---|---:|---:|---:|---:|---|\n'
for row in s['arms']:
 r=row['resources'];group=next(x for x in groups['final_subgroup_comparisons'] if x['arm_id']==row['arm_id'])
 section+=f"| {row['arm_id']} | ≥{row['final_rank_lower_bound']} | {row['charts']} | {r['wall_seconds']/60:.2f} | {(r['own_cpu_seconds']+r['child_cpu_seconds'])/60:.2f} | {group['status']} |\n"
same=next(r for r in s['arms'] if r['arm_id']=='original')['chart_centre_path_sha256']==next(r for r in s['arms'] if r['arm_id']=='recovered-P14')['chart_centre_path_sha256']
section+='\nThe original and recovered-P14 arms '+('tested identical ordered centre paths.' if same else 'tested different ordered centre paths.')+' No general basis invariance theorem follows.\nThe timing is one concurrent run per arm, without randomized timing trials;\nchart totals and exact paths are the stronger comparison. Once different\npoints are admitted, the arms have different adaptive reference lists, so\nlater costs include their effects.\n\n'
section+='Charts used for each successive gain:\n\n| Basis | Gain-stage chart counts, starting at17 |\n|---|---|\n'
for row in s['arms']:section+=f"| {row['arm_id']} | "+', '.join(map(str,row['charts_per_stage']))+' |\n'
section+=f"\nAcross the runs, {len(ledger['equations'])} distinct exact reduced quartic\nequations were tested; {ledger['new_equation_presentations_in_compared_catalogue']} are absent from the compared historical\nV1/V2/V3 chart files as exact coefficient pairs. **None is a new elliptic\ncurve:** every one is exactly birational over `Q` to302. This follows from\nverified PGL2/square-scaling identities and the universal inverse chord map,\nnot merely equality of j-invariants. All four conductors remain the same\n163-digit integer. The equation and point ledgers preserve all occurrences.\n\n"
section+='The [post-run group diagnostic](../../artifacts/generated-results/elliptic-curves/det1092_basis_cascade_v1/group_diagnostic.json)\ncertifies unions of the separately recovered intermediate groups and exact\nfinal group relations when found. These unions are post-run diagnostics;\ntheir points never feed a primary arm and their ranks do not measure a\nprospectively scheduled hybrid search. Height arithmetic only proposes\nrational words; exact rational group identities decide their validity.\n\n'
section+='The [completion seal](../../artifacts/generated-results/elliptic-curves/det1092_basis_cascade_v1/completion.json)\nbinds the [summary](../../artifacts/generated-results/elliptic-curves/det1092_basis_cascade_v1/summary.json),\nall four portable evidence archives, both replay layers, source/roster records,\nresource accounts, raw archive hashes and novelty ledgers. One evaluator-only\nuniversal-map preflight used the wrong centre sign; its failed assertion and\ncorrected exact identity check are retained in `verification_attempts.json`.\nNo primary arm was rerun or given a changed policy.\n\n'
section+='This result concerns representations and paths on one known curve. The next\ntransfer experiment must freeze genuinely different parent/fibre inputs,\nexclude exact isomorphic duplicates, record their conductor status, and\nkeep known exceptional points out of selection and search. Success here\ndoes not establish discovery power on that separate population.\n'
text+=section;note.write_text(text)
print(json.dumps(result,indent=2))
