#!/usr/bin/env python3
"""Portable source-only reduced chart and standalone generic identity replay."""
import argparse,json,hashlib,shutil,sys
from pathlib import Path
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];L=ROOT/'artifacts/local/elliptic-curves';ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ART/'det1092_reduced_parameter_chart_v1'
FILES={
 'source-parent.json':L/'det1092-point-pilot-v2/parent-sections.json',
 'chart-search.json':L/'det1092-parameter-lattice-reduction-v1/result.json',
 'reduced-parent.json':L/'det1092-reduced-chart-proof-v2/parent-sections.json',
 'proof-protocol.json':L/'det1092-reduced-chart-proof-v2/protocol.json',
 'generic-proof.json':L/'det1092-reduced-chart-proof-v2/result.json',
 'verify.sage':ROOT/'elliptic-curves/cas/verify_det1092_reduced_parameter_chart_v2.sage',
 'unit-relations.json':L/'det1092-unit-section-relations-v2/result.json',
 'integral-source.json':L/'det1092-integral-parameter-model-v1/result.json',
 'hessian-support.json':L/'det1092-parameter-hessian-support-v1/result.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_bytes())
def expected():
 for n,p in FILES.items():assert sha(D/n)==sha(p)
 source=read(D/'source-parent.json');original=ART/'curve302_recovered_mw17_parent_v1.json';a=read(original)
 assert source['source_sha256']==sha(original)
 assert set(source)=={'source_sha256','a_invariants','basis_weierstrass_coordinates'}
 assert all(source[k]==a[k] for k in ['a_invariants','basis_weierstrass_coordinates'])
 trial=read(D/'chart-search.json');proof=read(D/'generic-proof.json');assert trial['status']==proof['status']=='PASS'
 stages=[]
 for name in ['det1092-integral-parameter-model-v1','det1092-parameter-hessian-support-v1','det1092-parameter-lattice-reduction-v1','det1092-reduced-chart-proof-v1','det1092-reduced-chart-proof-v2','det1092-reduced-seed-diagnosis-v1','det1092-reduced-seed-diagnosis-v2','det1092-unit-section-relations-v1','det1092-unit-section-relations-v2']:
  p=L/name/'supervisor.json';r=read(p);stages.append(dict(name=name,outcome=r['outcome'],returncode=r['returncode'],wall_seconds=r['wall_seconds'],supervisor_sha256=sha(p)))
 return dict(schema='elliptic-curves.det1092-reduced-parameter-chart.v1',status='PASS',files={n:sha(D/n) for n in FILES},source_parent_sha256=sha(original),strict_contractions=len(trial['steps']),weighted_lattice_charts=len(trial['charts']),selected_coefficient_bits=proof['coefficient_bits'],selected_parameter_matrix=proof['parameter_matrix'],section_count=17,unit_section_span=15,supervised_stages=stages,total_derivation_diagnosis_seconds=sum(s['wall_seconds'] for s in stages),boundary='Exact Q(s) base and Weierstrass isomorphism and seventeen section transports; same parent and fibration. Nineteen strict local contractions followed by seventeen weighted two-dimensional LLL bases selected from source equation coefficients only. No global minimality, optimality, or completeness of partial-collapse primes claimed. Unit span15 follows from exact relations and independently certified fifteen-point subset. The first section proof wrongly required polynomial section coordinates; the first diagnosis constructed invalid Fraction(Sage QQ) values and incorrectly reported an empty mod2 audit. Those attempts are superseded, preserved and counted; their erroneous outputs are not evidence. Signed-kernel relation proposals were insufficient; bounded numerical proposals are accepted only after exact group equality. No new rank direction is claimed by the coordinate reduction.')
def main(check):
 if check:assert read(D/'manifest.json')==expected();print('PASS reduced chart portable manifest');return
 assert not D.exists();D.mkdir()
 for n,p in FILES.items():shutil.copy2(p,D/n)
 checkpoint(D/'manifest.json',expected())
 fresh=L/'det1092-reduced-chart-portable-replay-v1';fresh.mkdir(exist_ok=False)
 for n in ['source-parent.json','chart-search.json','proof-protocol.json','verify.sage']:shutil.copy2(D/n,fresh/n)
 cmd=['/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python',str(fresh/'verify.sage'),'--protocol',str(fresh/'proof-protocol.json'),'--parent',str(fresh/'source-parent.json'),'--chart',str(fresh/'chart-search.json'),'--output',str(fresh/'reduced-parent.json'),'--certificate',str(fresh/'generic-proof.json')]
 checkpoint(fresh/'replay-protocol.json',dict(files={p.name:sha(p) for p in fresh.iterdir()},seconds=120,rss_bytes=2147483648,workers=1))
 r=run(cmd,limits=Limits(120,2147483648),log_path=fresh/'replay.log',checkpoint_path=fresh/'supervisor.json',cwd=fresh)
 assert r['outcome']=='completed' and r['returncode']==0
 assert all(sha(fresh/n)==sha(D/n) for n in ['reduced-parent.json','generic-proof.json'])
 checkpoint(D/'portable-replay.json',dict(status='PASS',supervision=r,manifest_sha256=sha(D/'manifest.json'),matched_files={n:sha(fresh/n) for n in ['reduced-parent.json','generic-proof.json']}));print('PASS isolated generic chart replay',r['wall_seconds'])
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');main(p.parse_args().check)
