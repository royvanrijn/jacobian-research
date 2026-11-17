"""Complete replay for GAOF7FOUR2; mathematical norm lemma in COMPLETE_PROOF.md.
Run without -O. Each component has an explicit bounded subprocess timeout.
"""
from pathlib import Path
import hashlib,json,platform,subprocess,sys,time
HERE=Path(__file__).resolve().parent;START=time.monotonic()
if not __debug__:raise RuntimeError('Run without -O')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
components=[('maximum_multiplicity',HERE.parent/'gao-f7-linear-cube/verify_high_root.py',120,HERE.parent/'gao-f7-linear-cube/high-root-receipt.json'),('fivefold',HERE.parent/'gao-f7-quintuple/verify.py',120,HERE.parent/'gao-f7-quintuple/receipt.json'),('quadruple_source',HERE/'verify.py',120,HERE/'receipt.json'),('integer_norm',HERE/'replay_norm.py',240,HERE/'norm-grid-receipt.json'),('source_discriminant',HERE/'verify_discriminant.py',180,HERE/'discriminant-receipt.json'),('remaining_charts',HERE/'verify_remaining_charts.py',120,HERE/'remaining-charts-receipt.json')]
records=[]
with (HERE/'complete-replay.log').open('w') as fd:
 for name,checker,cap,receiptpath in components:
  print(json.dumps({'stage':name,'seconds':time.monotonic()-START}),flush=True)
  proc=subprocess.run([sys.executable,str(checker)],stdout=fd,stderr=subprocess.STDOUT,timeout=cap)
  assert proc.returncode==0,name
  data=json.loads(receiptpath.read_text());assert data['result']=='PASS' and data['checker_sha256']==sha(checker),name
  records.append({'name':name,'checker':str(checker.relative_to(HERE.parent)),'checker_sha256':sha(checker),'receipt':str(receiptpath.relative_to(HERE.parent)),'receipt_sha256':sha(receiptpath),'checks':data.get('checks'),'seconds':data.get('seconds')})
fac=json.loads((HERE/'crt-norm-factors-candidate.json').read_text());cp=json.loads((HERE/'coordinate-poles-candidate.json').read_text());assert list(map(int,fac['coordinate_pole_factor']['coefficients']))==list(map(int,cp['projection_coefficients']))
five=json.loads((HERE.parent/'gao-f7-quintuple/receipt.json').read_text());assert five['claim_id']=='GAOF7FIVE1' and five['normalized_target_count']==118
chart=json.loads((HERE/'remaining-charts-receipt.json').read_text());assert [r['rank'] for r in chart['charts']]==[174,21,21]
disc=json.loads((HERE/'discriminant-receipt.json').read_text());assert disc['residual_degree']==912
receipt={'result':'PASS','claim_id':'GAOF7FOUR2','scope':'Complete normalized F7 quadruple-or-higher strata over C on c!=0: exactly912 targets of multiplicity type4+2+1^6 (six source points), exactly118 of type5+1^7 (seven points), all others type4+1^8 (eight points). Hence exactly1030 nonordinary normalized quadruple targets. No4+3+2 type occurs. Combined with prior exclusions, only3+2+2+2 remains for a three-point F7 fibre.','normalized_counts':{'4+2+1^6':912,'5+1^7':118,'nonordinary_total':1030},'remaining_three_point_types':['3+2+2+2'],'components':records,'mathematical_argument':'COMPLETE_PROOF.md: discriminant singularity implies norm multiplicity >=2; squarefree degree912 residual and all-pole finite-free coverage give complete strata and counting.','checker_sha256':sha(Path(__file__)),'proof_sha256':sha(HERE/'COMPLETE_PROOF.md'),'source_input_sha256':sha(HERE/'quadruple-input.json'),'python':platform.python_version(),'seconds':time.monotonic()-START,'assurance':'Self-replay with two determinant/ring arithmetic engines and explicit written norm/counting argument. No independent end-to-end, formal or external review, author endorsement or novelty priority.'}
(HERE/'complete-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
