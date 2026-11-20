#!/usr/bin/env python3
"""Complete source support for the uniform quadratic ordinary-SFC tail.

The written analytic/root-count proof supplies the theorem. Replays validate
its source freshness and literal algebra/rational premises, not a formal proof.
"""
from __future__ import annotations
if not __debug__:raise SystemExit('Refusing -O/-OO: all exact checks must execute.')
from pathlib import Path
from hashlib import sha256
import argparse,copy,json,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent
PREFIX='research/factorial/image-rigidity/'
LOCAL={'UNIFORM_ACTUAL_PHASE_AND_CONTOUR_ESTIMATES.md','verify_uniform_actual_phase_and_contour_estimates.py','uniform_actual_phase_and_contour_verification.json','UNIFORM_QUADRATIC_STRONG_FACTORIAL_TAIL.md','verify_uniform_quadratic_strong_factorial_tail.py','uniform_quadratic_strong_factorial_tail_verification.json','verify_uniform_quadratic_strong_factorial_completion.py','uniform_quadratic_strong_factorial_completion_inputs.json'}
REPO={PREFIX+p for p in{'GAP_TWO_EFFECTIVE_LARGE_LOWER_EXPONENT_THEOREM.md','verify_gap_two_effective_large_lower_completion.py','gap_two_effective_large_lower_completion_source_ledger.json','gap_two_effective_large_lower_completion_verification.json','GAMMA_PRODUCT_ROOT_MODULUS.tex'}}
REQUIRED={('local',p)for p in LOCAL}|{('repo',p)for p in REPO}
FORBIDDEN_400={'ITERATED_GRACE_POLYNOMIAL_ROOT_STRIP.md','verify_iterated_grace_polynomial_root_strip.py','iterated_grace_polynomial_root_strip_source_ledger.json','iterated_grace_polynomial_root_strip_inputs.json','iterated_grace_polynomial_root_strip_verification.json','UNIFORM_QUARTIC_STRONG_FACTORIAL_TAIL.md'}
EXPECTED={'schema':1,'domain':{'integer_d_min':2,'integer_N_min':2,'ordinary_integer_a_lower':'10^24*d^2','homotopy_real_a':'d^2/delta^2','homotopy_delta_max':'1/1000000000000','homotopy_m_min':'10000000000000000','homotopy_requires_hard_order':False,'uniform_phase_actual_delta_max':'1/10000000000','uniform_phase_outer_delta_max':'1/100000'},'actual_moment':'A_(m;a,d)(Z)=sum binom(m,k)(a*m+1)_(d*k) Z^(m-k)','rectangle':{'real_abs_max':'1/10000000','imaginary_abs_max':'21/10'},'quad_expected':{'rational_guards':18,'universal_identities':3,'exact_laurent_bell_cases':27,'negative_controls':4},'phase_expected':{'rational_guards':70,'universal_identities':3,'falling_coefficient_controls':2079,'negative_controls':2},'coverage':[{'declaration':'Theorem UniformQuadraticStrongFactorialTail','path':'UNIFORM_QUADRATIC_STRONG_FACTORIAL_TAIL.md'},{'declaration':'Lemma UniformActualMomentRootRectangle','path':'UNIFORM_QUADRATIC_STRONG_FACTORIAL_TAIL.md'},{'declaration':'Corollary QuadraticBadLowerExponentEquivalence','path':'UNIFORM_QUADRATIC_STRONG_FACTORIAL_TAIL.md'},{'declaration':'Lemma UniformActualPhaseMajorants','path':'UNIFORM_ACTUAL_PHASE_AND_CONTOUR_ESTIMATES.md'},{'declaration':'Lemma UniformAllHeightOnePacket','path':'UNIFORM_ACTUAL_PHASE_AND_CONTOUR_ESTIMATES.md'}],'scope':'Every characteristic-zero field, nonzero ordinary binomials, all integer gaps d>=2 and every terminal N>=2 at a>=10^24*d^2. Bad_d is bounded by the same quadratic cutoff.','excludes':['full ordinary binomial SFC','Bad_d empty','a>=100*d^2 universal claim','any bound on bad N below the cutoff','formal verification'],'main_route_excludes_400':True}
def sha(path):return sha256(path.read_bytes()).hexdigest()
def validate(ledger,repo):
 assert isinstance(ledger,dict)and ledger.get('schema')==1
 assert isinstance(ledger.get('files'),list)and len(ledger['files'])==13
 assert all(isinstance(v,dict)and set(v)=={'root','path','sha256'}for v in ledger['files'])
 assert{(v['root'],v['path'])for v in ledger['files']}==REQUIRED
 for v in ledger['files']:
  assert v['root']in{'local','repo'}and isinstance(v['sha256'],str)and len(v['sha256'])==64 and all(c in'0123456789abcdef'for c in v['sha256'])
  p=(HERE if v['root']=='local'else repo)/v['path'];assert p.is_file()and sha(p)==v['sha256'],v['path']
def reject(v,repo):
 try:validate(v,repo)
 except(AssertionError,KeyError,TypeError):return
 raise AssertionError('malformed ledger accepted')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,default=HERE/'uniform_quadratic_strong_factorial_completion_verification.json');args=ap.parse_args();repo=args.repo_root.resolve()
 lp=HERE/'uniform_quadratic_strong_factorial_completion_source_ledger.json';ledger=json.loads(lp.read_text());validate(ledger,repo)
 inp=json.loads((HERE/'uniform_quadratic_strong_factorial_completion_inputs.json').read_text());assert inp==EXPECTED
 assert len(LOCAL)==8 and len(REPO)==5 and not any(Path(p).name in FORBIDDEN_400 for _,p in REQUIRED)
 # Verify the entire source-ledger import closure also avoids400. Some399
 # ledgers pin older classical391 inputs; those are legitimate retained sources.
 start=repo/PREFIX/'gap_two_effective_large_lower_completion_source_ledger.json'
 seen=set();stack=[start]
 while stack:
  path=stack.pop().resolve()
  if path in seen:continue
  seen.add(path);data=json.loads(path.read_text())
  # Retained398 uses its pinned schema_version=1 spelling; new ledgers
  # and every other imported ledger use schema=1. Its own fresh checker
  # validates that exact legacy schema rather than silently defaulting.
  assert data.get('schema')==1 or(path.name=='explicit_airy_contour_envelope_source_ledger.json'and data.get('schema_version')==1 and 'schema'not in data)
  for row in data.get('files',[]):
   assert Path(row['path']).name not in FORBIDDEN_400
   q=(path.parent if row['root']=='local'else repo)/row['path']
   if q.name.endswith('_source_ledger.json'):stack.append(q)
 replays=[]
 with tempfile.TemporaryDirectory(prefix='uniform_quadratic_completion_')as td:
  td=Path(td)
  for j,(script,receipt)in enumerate([('verify_uniform_actual_phase_and_contour_estimates.py','uniform_actual_phase_and_contour_verification.json'),('verify_uniform_quadratic_strong_factorial_tail.py','uniform_quadratic_strong_factorial_tail_verification.json')]):
   out=td/f'local{j}.json';run=subprocess.run([sys.executable,str(HERE/script),'--output',str(out)],capture_output=True,text=True);assert run.returncode==0,(script,run.stdout,run.stderr)
   assert out.read_bytes()==(HERE/receipt).read_bytes(),receipt
   replays.append(dict(root='local',script=script,receipt=receipt,sha256=sha(out),data=json.loads(out.read_text())))
  out=td/'canonical399.json';script=repo/PREFIX/'verify_gap_two_effective_large_lower_completion.py';receipt=repo/PREFIX/'gap_two_effective_large_lower_completion_verification.json'
  run=subprocess.run([sys.executable,str(script),'--repo-root',str(repo),'--output',str(out)],capture_output=True,text=True);assert run.returncode==0,(run.stdout,run.stderr)
  assert out.read_bytes()==receipt.read_bytes();replays.append(dict(root='repo',script=PREFIX+script.name,receipt=PREFIX+receipt.name,sha256=sha(out),data=json.loads(out.read_text())))
 phase,quad,support=[r['data']for r in replays]
 assert len(phase['rational_guards'])==70 and phase['universal_identity_count']==3 and phase['falling_coefficient_controls']==2079 and phase['negative_controls']==2
 assert len(quad['rational_guards'])==18 and quad['universal_identity_count']==3 and len(quad['exact_laurent_bell_cases'])==27 and quad['negative_controls']==4
 assert support['mandatory_local_pins']==31 and support['mandatory_repo_pins']==22 and support['fresh_local_child_replays']==4 and support['fresh_imported_replays']==5
 failures=0
 for j in range(13):
  bad=copy.deepcopy(ledger);del bad['files'][j];reject(bad,repo);failures+=1
 for key in['schema','files']:
  bad=copy.deepcopy(ledger);del bad[key];reject(bad,repo);failures+=1
 for field,value in[('schema',2),('files',[])]:
  bad=copy.deepcopy(ledger);bad[field]=value;reject(bad,repo);failures+=1
 bad=copy.deepcopy(ledger);bad['files'].append(copy.deepcopy(bad['files'][0]));reject(bad,repo);failures+=1
 bad=copy.deepcopy(ledger);bad['files'][0]['sha256']='0'*64;reject(bad,repo);failures+=1
 bad=copy.deepcopy(ledger);bad['files'][0]['root']='other';reject(bad,repo);failures+=1
 bad=copy.deepcopy(ledger);bad['files'][0]['path']='missing';reject(bad,repo);failures+=1
 for path,value in[('homotopy_requires_hard_order',True),('homotopy_real_a','a must be an integer'),('ordinary_integer_a_lower','100*d^2')]:
  wrong=copy.deepcopy(inp);wrong['domain'][path]=value;assert wrong!=EXPECTED;failures+=1
 for opt in['-O','-OO']:
  run=subprocess.run([sys.executable,opt,str(Path(__file__))],capture_output=True,text=True);assert run.returncode!=0 and 'Refusing -O/-OO'in run.stderr;failures+=1
 result=dict(schema=1,status='PASS_SOURCE_GUARDED_UNIFORM_QUADRATIC_SUPPORT',scope=inp['scope'],coverage=inp['coverage'],excluded_conclusions=inp['excludes'],mandatory_local_pins=8,mandatory_repo_pins=5,main_route_imports_400=False,checked_source_ledger_closure_count=len(seen),fresh_replays=3,source_failure_and_scope_controls=failures,pins=ledger['files'],source_ledger_sha256=sha(lp),replays=replays,formal_verification=False,external_review=False,independent_replay=False)
 args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(dict(status='PASS',mandatory_pins=13,fresh_replays=3,source_scope_controls=failures,source_ledger_closure=len(seen),imports_400=False)))
if __name__=='__main__':main()
