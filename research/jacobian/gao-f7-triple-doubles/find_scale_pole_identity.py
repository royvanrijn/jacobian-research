"""Bounded weighted Macaulay x^3 certificate, adapted from find_identity.py.
Modular row-space calculations choose support only. Acceptance requires a
complete literal characteristic-zero polynomial identity, not modular rank.
"""
from pathlib import Path
import hashlib,json,resource,subprocess,sys,time
HERE=Path(__file__).resolve().parent;START=time.monotonic()
def emit(x):print(json.dumps(x,default=int),flush=True)
def worker():
 from flint import nmod_mat,fmpq_mat,fmpq_mpoly_ctx,fmpq
 data=json.loads((HERE/'two-scale-pole-power-input.json').read_text());variables=data['variables'];weights=[1,1,3,3];generators=[{tuple(ex):int(c) for ex,c in r['terms']} for r in data['equations']];degrees=[max(sum(e*w for e,w in zip(ex,weights)) for ex in f) for f in generators];order=sorted(range(len(generators)),key=lambda i:degrees[i]);target={(0,0,0,3):1};prime=1000000007;rounds=[]
 def monomials(bound,i=0):
  if bound<0:return
  if i==len(weights):yield ();return
  for k in range(bound//weights[i]+1):
   for rest in monomials(bound-k*weights[i],i+1):yield (k,)+rest
 def pivots(M,rank):return [next(j for j in range(M.ncols()) if M[i,j]) for i in range(rank)]
 for bound in [12,14,16,18,20,22,24]:
  columns=[(i,m) for i in order for m in monomials(bound-degrees[i])];sparse=[];rowset=set(target)
  for i,m in columns:
   col={tuple(a+b for a,b in zip(ex,m)):c for ex,c in generators[i].items()};sparse.append(col);rowset.update(col)
  rows=sorted(rowset,key=lambda m:(sum(e*w for e,w in zip(m,weights)),m));where={m:i for i,m in enumerate(rows)};report={'bound':bound,'rows':len(rows),'columns':len(columns),'seconds':time.monotonic()-START}
  if len(rows)*(len(columns)+1)>15000000:report['outcome']='matrix_entry_cap';rounds.append(report);break
  emit(report);M=nmod_mat(len(rows),len(columns)+1,prime)
  for j,col in enumerate(sparse):
   for ex,c in col.items():M[where[ex],j]=c%prime
  M[where[(0,0,0,3)],len(columns)]=1
  reduced,rank=M.rref(inplace=True);pivot=pivots(reduced,rank);base=[j for j in pivot if j<len(columns)];report['modular_rank']=len(base)
  if len(columns) in pivot:
   report['outcome']='no_identity_at_bound_mod_prime';rounds.append(report);emit(report);del M,reduced;continue
  support=[base[i] for i in range(len(base)) if reduced[i,len(columns)]];report['support']=len(support);emit(report)
  if len(support)>2500:report['outcome']='exact_dimension_cap';rounds.append(report);break
  del M,reduced
  transpose=nmod_mat(len(support),len(rows),prime)
  for k,j in enumerate(support):
   for ex,c in sparse[j].items():transpose[k,where[ex]]=c%prime
  rr,rk=transpose.rref(inplace=True);selected=pivots(rr,rk);assert rk==len(support);del transpose,rr
  checkpoint={'status':'EXACT_SOLVE_PENDING','source_sha256':hashlib.sha256((HERE/'two-scale-pole-power-input.json').read_bytes()).hexdigest(),'bound':bound,'support':[{'generator':columns[j][0],'multiplier_exponents':list(columns[j][1])} for j in support],'selected_monomials':[list(rows[i]) for i in selected]}
  (HERE/'scale-pole-macaulay-checkpoint.json').write_text(json.dumps(checkpoint,separators=(',',':'))+'\n');emit({'stage':'exact_solve','dimension':len(support),'seconds':time.monotonic()-START})
  A=fmpq_mat([[sparse[j].get(rows[i],0) for j in support] for i in selected]);B=fmpq_mat([[target.get(rows[i],0)] for i in selected]);solution=A.solve(B,algorithm='dixon');ctx=fmpq_mpoly_ctx.get(variables);multipliers=[{} for _ in generators]
  for k,j in enumerate(support):
   if solution[k,0]:i,m=columns[j];multipliers[i][m]=solution[k,0]
  total=sum((ctx.from_dict(m)*ctx.from_dict({ex:fmpq(c) for ex,c in f.items()}) for m,f in zip(multipliers,generators)),ctx.constant(0));assert total==ctx.from_dict({ex:fmpq(c) for ex,c in target.items()}),'full rational identity failed'
  cert={'status':'EXACT_IDENTITY_CANDIDATE_REQUIRES_INDEPENDENT_REPLAY','variables':variables,'source_sha256':checkpoint['source_sha256'],'target_terms':[[list(ex),str(c)] for ex,c in target.items()],'weights':weights,'bound':bound,'multipliers':[[[list(ex),str(c)] for ex,c in sorted(m.items())] for m in multipliers]}
  (HERE/'scale-pole-identity-candidate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n');report.update(outcome='EXACT_LITERAL_IDENTITY_PASS',multiplier_terms=sum(len(m) for m in multipliers));rounds.append(report);emit(report);break
 (HERE/'scale-pole-macaulay-search.json').write_text(json.dumps({'prime':prime,'weights':weights,'degrees':degrees,'rounds':rounds,'seconds':time.monotonic()-START,'boundary':'Only the exported literal rational identity has mathematical force; independent source/chart and arithmetic replay remains required.'},indent=2)+'\n')
if '--worker' in sys.argv:worker()
else:
 def limits():
  resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));resource.setrlimit(resource.RLIMIT_FSIZE,(50*1024**2,50*1024**2))
 rec={'cap_seconds':120,'address_cap_GiB':4,'per_file_cap_MiB':50,'matrix_entry_cap':15000000,'exact_dimension_cap':2500,'weighted_bounds':[12,14,16,18,20,22,24],'theorem_exclusion':False,'motivation':'Replace the timed-out rational Groebner lift with a finite weighted Macaulay support solve for x^3; dimensions are preflighted and exact acceptance is a full identity.'}
 with (HERE/'scale-pole-macaulay.log').open('wb') as fd:
  try:
   proc=subprocess.run([sys.executable,str(Path(__file__)),'--worker'],stdout=fd,stderr=subprocess.STDOUT,timeout=120,preexec_fn=limits);rec.update(outcome='completed',returncode=proc.returncode)
  except subprocess.TimeoutExpired:rec.update(outcome='timeout',returncode=None)
 rec['seconds']=time.monotonic()-START
 for fn in ['scale-pole-identity-candidate.json','scale-pole-macaulay-checkpoint.json','scale-pole-macaulay-search.json']:
  p=HERE/fn
  if p.exists():rec[fn]={'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
 (HERE/'scale-pole-macaulay-run.json').write_text(json.dumps(rec,indent=2)+'\n');emit(rec);print((HERE/'scale-pole-macaulay.log').read_text()[-4500:])
