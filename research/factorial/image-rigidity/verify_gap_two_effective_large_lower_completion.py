#!/usr/bin/env python3
"""Source-guarded complete support for the effective gap-two theorem.

The written proof supplies the analytic argument. This program checks all
source freshness, exact certificate replays and rational/algebraic budgets;
it is not a formal verification or an inference from sampled moments.
"""
from __future__ import annotations
if not __debug__:raise SystemExit('Refusing -O/-OO: complete exact checks are required.')
from pathlib import Path
from hashlib import sha256
import argparse,copy,json,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent
PREFIX='research/factorial/image-rigidity/'
LOCAL={
 'GAP_TWO_EFFECTIVE_LARGE_LOWER_EXPONENT_THEOREM.md',
 'verify_gap_two_effective_large_lower_exponent.py',
 'gap_two_effective_large_lower_exponent_verification.json',
 'verify_gap_two_effective_large_lower_completion.py',
 'gap_two_effective_large_lower_completion_inputs.json',
 'GAP_TWO_HORIZONTAL_ENDPOINT_CUBIC_LOCATION.md',
 'verify_gap_two_horizontal_endpoint_cubic_balls.py',
 'gap_two_horizontal_endpoint_cubic_balls_verification.json',
 'verify_gap_two_horizontal_endpoint_location.py',
 'gap_two_horizontal_endpoint_inputs.json',
 'gap_two_horizontal_endpoint_source_ledger.json',
 'gap_two_horizontal_endpoint_location_verification.json',
 'gap_two_horizontal_endpoint_primary_audit.md',
 'NEGATIVE_AIRY_CURVE_OUTER_TAIL.md',
 'verify_negative_airy_curve_outer_tail.py',
 'negative_airy_curve_outer_tail_verification.json',
 'POSITIVE_AIRY_THIN_OUTER_TAIL_AND_CUT_GEOMETRY.md',
 'verify_positive_airy_thin_outer_tail.py',
 'positive_airy_thin_outer_tail_verification.json',
 'COMPACT_AIRY_PHYSICAL_CUT_ADDENDUM.md',
 'verify_compact_airy_physical_cut.py',
 'compact_airy_physical_cut_verification.json',
 'verify_airy_outer_tail_and_cut_package.py',
 'airy_outer_tail_and_cut_inputs.json',
 'airy_outer_tail_and_cut_source_ledger.json',
 'airy_outer_tail_and_cut_verification.json',
 'ACTUAL_SHIFTED_HERMITE_NONLOCAL_BOUND.md',
 'verify_actual_shifted_hermite_nonlocal.py',
 'actual_shifted_hermite_nonlocal_inputs.json',
 'actual_shifted_hermite_nonlocal_source_ledger.json',
 'actual_shifted_hermite_nonlocal_verification.json',
}
REPO_BASES={
 'QUANTITATIVE_GAP_TWO_ENDPOINT_DATA.md',
 'verify_quantitative_gap_two_endpoint_data.py',
 'quantitative_gap_two_endpoint_data_source_ledger.json',
 'quantitative_gap_two_endpoint_data_verification.json',
 'HERMITE_CRITICAL_CUBIC_INVERSE_DISK.md',
 'GAP_TWO_FULL_ENDPOINT_UNIFORM_CUBIC_INVERSE.md',
 'verify_gap_two_full_endpoint_uniform_inverse.py',
 'gap_two_full_endpoint_inverse_source_ledger.json',
 'gap_two_full_endpoint_inverse_verification.json',
 'EXPLICIT_AIRY_CONTOUR_ENVELOPE.md',
 'verify_explicit_airy_contour_envelope.py',
 'explicit_airy_contour_envelope_source_ledger.json',
 'explicit_airy_contour_envelope_verification.json',
 'GAP_TWO_EXTENDED_BULK_ENDPOINT_REDUCTION.md',
 'verify_gap_two_endpoint_reduction_guards.py',
 'gap_two_endpoint_reduction_source_ledger.json',
 'gap_two_endpoint_reduction_guards_verification.json',
 'GAP_TWO_LARGE_LOWER_EXPONENT_ENDPOINT_THEOREM.md',
 'verify_gap_two_endpoint_contour_identities.py',
 'gap_two_large_lower_endpoint_source_ledger.json',
 'gap_two_endpoint_contour_identities_verification.json',
 'GAMMA_PRODUCT_ROOT_MODULUS.tex',
}
REPO={PREFIX+s for s in REPO_BASES}
REQUIRED={('local',p)for p in LOCAL}|{('repo',p)for p in REPO}
assert len(LOCAL)==31 and len(REPO)==22 and len(REQUIRED)==53
EXPECTED_INPUT={
 'schema':1,'theorem':'Theorem EffectiveGapTwoLargeLowerExponent',
 'domain':{'a_min':'40000000000000000','N_min':2,'hard_order':'N>a/4+3/2','local_m_min':'10000000000000000','delta_max':'1/100000000','real_z_abs_max':'1/100000','imaginary_z_bounds':['18999/10000','21001/10000']},
 'remainder':{'compact_xi_abs_max':4,'compact_constant':'100000000000000','large_constant':'1000000000','normalization':'sigma*2*pi*i*m^(-1/3)*exp(m*Phi) after exact nonzero actual moment scalars are removed'},
 'actual_polynomial':'A_m(Z)=sum binom(m,k)(a*m+1)_(2k) Z^(m-k)',
 'coefficient_scope':'Every characteristic-zero field; NONZERO ordinary X^a(A+B*X^2)',
 'expected_arithmetic':{'rational_guards':44,'universal_identities':3,'normalization_controls':4,'negative_controls':3},
 'coverage':[{'declaration':'Theorem EffectiveGapTwoLargeLowerExponent','support':'Written complete actual-period proof, separately proved analytic inputs, exact outward-ball endpoint certificates, and rational/identity budget controls.'}],
 'excludes':['all a>=40000 gap-two claim','uniform complete-period constant C_R=10^10','all-gap Strong Factorial conjecture','formal verification'],
}
def sha(path):return sha256(path.read_bytes()).hexdigest()
def validate_ledger(ledger,base,repo):
 assert isinstance(ledger,dict)and ledger.get('schema')==1
 assert isinstance(ledger.get('files'),list)and len(ledger['files'])==len(REQUIRED)
 assert all(isinstance(row,dict)and set(row)=={'root','path','sha256'}for row in ledger['files'])
 assert {(row['root'],row['path'])for row in ledger['files']}==REQUIRED
 for row in ledger['files']:
  assert row['root']in{'local','repo'}
  assert isinstance(row['sha256'],str)and len(row['sha256'])==64 and all(c in'0123456789abcdef'for c in row['sha256'])
  path=(base if row['root']=='local'else repo)/row['path']
  assert path.is_file()and sha(path)==row['sha256'],row['path']
def must_reject(ledger,base,repo):
 try:validate_ledger(ledger,base,repo)
 except(AssertionError,KeyError,TypeError):return
 raise AssertionError('Malformed source ledger was accepted')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,default=HERE/'gap_two_effective_large_lower_completion_verification.json');args=ap.parse_args();repo=args.repo_root.resolve()
 lp=HERE/'gap_two_effective_large_lower_completion_source_ledger.json';ledger=json.loads(lp.read_text());validate_ledger(ledger,HERE,repo)
 inp=json.loads((HERE/'gap_two_effective_large_lower_completion_inputs.json').read_text());assert inp==EXPECTED_INPUT
 local_children=[
  ('verify_gap_two_horizontal_endpoint_location.py','gap_two_horizontal_endpoint_location_verification.json',True),
  ('verify_airy_outer_tail_and_cut_package.py','airy_outer_tail_and_cut_verification.json',True),
  ('verify_actual_shifted_hermite_nonlocal.py','actual_shifted_hermite_nonlocal_verification.json',True),
  ('verify_gap_two_effective_large_lower_exponent.py','gap_two_effective_large_lower_exponent_verification.json',False),
 ]
 imported_children=[
  ('verify_gap_two_endpoint_reduction_guards.py','gap_two_endpoint_reduction_guards_verification.json'),
  ('verify_gap_two_endpoint_contour_identities.py','gap_two_endpoint_contour_identities_verification.json'),
  ('verify_quantitative_gap_two_endpoint_data.py','quantitative_gap_two_endpoint_data_verification.json'),
  ('verify_gap_two_full_endpoint_uniform_inverse.py','gap_two_full_endpoint_inverse_verification.json'),
  ('verify_explicit_airy_contour_envelope.py','explicit_airy_contour_envelope_verification.json'),
 ]
 results=[]
 with tempfile.TemporaryDirectory(prefix='effective_gap_two_complete_replay_')as temp:
  work=Path(temp)
  for index,(script,receipt,needs_repo)in enumerate(local_children):
   output=work/f'local_{index}.json';cmd=[sys.executable,str(HERE/script)]
   if needs_repo:cmd+=['--repo-root',str(repo)]
   run=subprocess.run(cmd+['--output',str(output)],capture_output=True,text=True)
   assert run.returncode==0,(script,run.stdout,run.stderr)
   assert output.read_bytes()==(HERE/receipt).read_bytes(),('stale local receipt',receipt)
   data=json.loads(output.read_text());results.append(dict(root='local',script=script,receipt=receipt,receipt_sha256=sha(output),data=data))
  for index,(script,receipt)in enumerate(imported_children):
   output=work/f'repo_{index}.json';cmd=[sys.executable,str(repo/PREFIX/script),'--repo-root',str(repo),'--output',str(output)]
   run=subprocess.run(cmd,capture_output=True,text=True);assert run.returncode==0,(script,run.stdout,run.stderr)
   assert output.read_bytes()==(repo/PREFIX/receipt).read_bytes(),('stale imported receipt',receipt)
   data=json.loads(output.read_text());results.append(dict(root='repo',script=PREFIX+script,receipt=PREFIX+receipt,receipt_sha256=sha(output),data=data))
 endpoint,tails,nonlocal_data,arith=[r['data']for r in results[:4]]
 assert endpoint['independent_ball_leaf_checks']==7262 and endpoint['independent_shared_endpoint_discs']==7264
 assert len(endpoint['rational_guards'])==14 and endpoint['negative_controls']==5
 assert len(tails['replays'])==3 and tails['status']=='PASS_CANONICAL_AIRY_TAILS_AND_CUTS'
 assert len(nonlocal_data['rational_guards'])==32 and nonlocal_data['negative_controls']==2
 assert len(arith['rational_guards'])==44 and arith['universal_identities']==3 and len(arith['actual_moment_normalization_controls'])==4 and arith['negative_controls']==3
 # Source failure controls are complete over the REQUIRED set: withdrawing
 # any one proof/source/input/ledger/receipt is rejected, not silently ignored.
 rejected=0
 for j in range(len(ledger['files'])):
  bad=copy.deepcopy(ledger);del bad['files'][j];must_reject(bad,HERE,repo);rejected+=1
 for key in ['schema','files']:
  bad=copy.deepcopy(ledger);del bad[key];must_reject(bad,HERE,repo);rejected+=1
 bad=copy.deepcopy(ledger);bad['schema']=2;must_reject(bad,HERE,repo);rejected+=1
 bad=copy.deepcopy(ledger);bad['files']=[];must_reject(bad,HERE,repo);rejected+=1
 bad=copy.deepcopy(ledger);bad['files'].append(copy.deepcopy(bad['files'][0]));must_reject(bad,HERE,repo);rejected+=1
 bad=copy.deepcopy(ledger);bad['files'][0]['root']='untrusted';must_reject(bad,HERE,repo);rejected+=1
 bad=copy.deepcopy(ledger);bad['files'][0]['sha256']='0'*64;must_reject(bad,HERE,repo);rejected+=1
 bad=copy.deepcopy(ledger);bad['files'][0]['path']='missing_source.md';must_reject(bad,HERE,repo);rejected+=1
 bad=copy.deepcopy(ledger);del bad['files'][0]['sha256'];must_reject(bad,HERE,repo);rejected+=1
 wrong=copy.deepcopy(inp);wrong['domain']['a_min']='40000';assert wrong!=EXPECTED_INPUT;rejected+=1
 wrong=copy.deepcopy(inp);wrong['remainder']['large_constant']='100000000000000';assert wrong!=EXPECTED_INPUT;rejected+=1
 for opt in['-O','-OO']:
  run=subprocess.run([sys.executable,opt,str(Path(__file__))],capture_output=True,text=True)
  assert run.returncode!=0 and 'Refusing -O/-OO'in run.stderr;rejected+=1
 result=dict(schema=1,status='PASS_SOURCE_GUARDED_EFFECTIVE_GAP_TWO_SUPPORT',scope='Supporting certificates for the full written theorem: all integer a>=4e16 and actual ordinary gap-two consecutive coprimality over every characteristic-zero field. The executable is not a formal verification of analytic arguments.',coverage=inp['coverage'],excluded_conclusions=inp['excludes'],mandatory_local_pins=len(LOCAL),mandatory_repo_pins=len(REPO),source_failure_and_scope_controls=rejected,optimized_modes_rejected=['-O','-OO'],fresh_local_child_replays=4,fresh_imported_replays=5,replays=results,pins=ledger['files'],source_ledger_sha256=sha(lp),formal_verification=False,external_review=False,independent_replay=False)
 args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
 print(json.dumps(dict(status='PASS',mandatory_pins=len(REQUIRED),fresh_child_replays=len(results),source_and_scope_controls=rejected,main_arithmetic_guards=44)))
if __name__=='__main__':main()
