"""Bounded global delta-zero normalized carrier-power Macaulay search.
Modular pivots select support; only an exact full rational identity is accepted.
The toric relation is omitted, so any accepted x^k identity is a stronger
obstruction on the full normalized coefficient family.
"""
from pathlib import Path
import hashlib,json,resource,subprocess,sys,time
HERE=Path(__file__).resolve().parent;START=time.monotonic()
def emit(x):print(json.dumps(x,default=int),flush=True)
def worker():
 from flint import nmod_mat,fmpq_mat,fmpq_mpoly_ctx,fmpq
 source=HERE/'two-scale-normalized-input.json';data=json.loads(source.read_text());variables=data['variables'];weights=data['weights'];generators=[{tuple(ex):int(c) for ex,c in r['terms']} for r in data['equations']];degrees=[max(sum(e*w for e,w in zip(ex,weights)) for ex in f) for f in generators];order=sorted(range(len(generators)),key=lambda i:degrees[i]);prime=1000000007;rounds=[]
 def monomials(B,i=0):
  if B<0:return
  if i==len(weights):yield ();return
  for a in range(B//weights[i]+1):
   for tail in monomials(B-a*weights[i],i+1):yield (a,)+tail
 def pivots(M,rank):return [next(j for j in range(M.ncols()) if M[i,j]) for i in range(rank)]
 for bound in [16,18,20,22]:
  powers=[k for k in range(3,8) if 3*k<=bound];targets=[{(0,0,0,0,k):1} for k in powers];columns=[(i,m) for i in order for m in monomials(bound-degrees[i])];sparse=[];rowset=set().union(*targets)
  for i,m in columns:
   col={tuple(a+b for a,b in zip(ex,m)):c for ex,c in generators[i].items()};sparse.append(col);rowset.update(col)
  rows=sorted(rowset,key=lambda m:(sum(e*w for e,w in zip(m,weights)),m));where={m:i for i,m in enumerate(rows)};report={'bound':bound,'rows':len(rows),'columns':len(columns),'target_powers':powers,'seconds':time.monotonic()-START}
  if len(rows)*(len(columns)+len(targets))>50000000:report['outcome']='matrix_entry_cap';rounds.append(report);break
  emit(report);M=nmod_mat(len(rows),len(columns)+len(targets),prime)
  for j,col in enumerate(sparse):
   for ex,c in col.items():M[where[ex],j]=c%prime
  for k,target in enumerate(targets):
   for ex,c in target.items():M[where[ex],len(columns)+k]=c
  R,rank=M.rref(inplace=True);pivot=pivots(R,rank);base=[j for j in pivot if j<len(columns)];base_rank=len(base);consistent=[k for k in range(len(targets)) if all(R[i,len(columns)+k]==0 for i in range(base_rank,rank))];report['modular_rank']=base_rank;report['consistent_powers']=[powers[k] for k in consistent]
  if not consistent:report['outcome']='no_identity_at_bound_mod_prime';rounds.append(report);emit(report);del M,R;continue
  chosen=min(consistent,key=lambda k:sum(bool(R[i,len(columns)+k]) for i in range(base_rank)));support=[base[i] for i in range(base_rank) if R[i,len(columns)+chosen]];report['support']=len(support);report['chosen_power']=powers[chosen];emit(report)
  if len(support)>5000:report['outcome']='exact_dimension_cap';rounds.append(report);break
  del M,R
  transpose=nmod_mat(len(support),len(rows),prime)
  for k,j in enumerate(support):
   for ex,c in sparse[j].items():transpose[k,where[ex]]=c%prime
  rr,rk=transpose.rref(inplace=True);selected=pivots(rr,rk);assert rk==len(support);del transpose,rr
  checkpoint={'status':'EXACT_SOLVE_PENDING','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'bound':bound,'target_power':powers[chosen],'support':[{'generator':columns[j][0],'multiplier_exponents':list(columns[j][1])} for j in support],'selected_monomials':[list(rows[i]) for i in selected]}
  (HERE/'normalized-macaulay-checkpoint.json').write_text(json.dumps(checkpoint,separators=(',',':'))+'\n');emit({'stage':'exact_solve','dimension':len(support),'power':powers[chosen],'seconds':time.monotonic()-START})
  A=fmpq_mat([[sparse[j].get(rows[i],0) for j in support] for i in selected]);B=fmpq_mat([[targets[chosen].get(rows[i],0)] for i in selected]);solution=A.solve(B,algorithm='dixon');ctx=fmpq_mpoly_ctx.get(variables);multipliers=[{} for _ in generators]
  for k,j in enumerate(support):
   if solution[k,0]:i,m=columns[j];multipliers[i][m]=solution[k,0]
  total=sum((ctx.from_dict(m)*ctx.from_dict({ex:fmpq(c) for ex,c in f.items()}) for m,f in zip(multipliers,generators)),ctx.constant(0));assert total==ctx.from_dict({ex:fmpq(c) for ex,c in targets[chosen].items()}),'full rational identity failed'
  certificate={'status':'EXACT_IDENTITY_CANDIDATE_REQUIRES_INDEPENDENT_REPLAY','variables':variables,'source_sha256':checkpoint['source_sha256'],'target_power':powers[chosen],'target_terms':[[list(ex),str(c)] for ex,c in targets[chosen].items()],'weights':weights,'bound':bound,'multipliers':[[[list(ex),str(c)] for ex,c in sorted(m.items())] for m in multipliers]}
  (HERE/'normalized-identity-candidate.json').write_text(json.dumps(certificate,separators=(',',':'))+'\n');report.update(outcome='EXACT_LITERAL_IDENTITY_PASS',multiplier_terms=sum(len(m) for m in multipliers));rounds.append(report);emit(report);break
 (HERE/'normalized-macaulay-search.json').write_text(json.dumps({'prime':prime,'weights':weights,'degrees':degrees,'rounds':rounds,'seconds':time.monotonic()-START,'boundary':'Finite-field selection alone proves nothing about complex solutions; only full literal rational identities are accepted.'},indent=2)+'\n')
if '--worker' in sys.argv:worker()
else:
 def limits():
  resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));resource.setrlimit(resource.RLIMIT_FSIZE,(100*1024**2,100*1024**2))
 rec={'cap_seconds':180,'address_cap_GiB':4,'per_file_cap_MiB':100,'matrix_entry_cap':50000000,'exact_dimension_cap':5000,'weighted_bounds':[16,18,20,22],'target_powers':[3,4,5,6,7],'theorem_exclusion':False,'motivation':'The four-variable pole x^3 target yielded a literal identity at bound22. The corrected global two-scale normalization retains positive weights[1,2,1,3,3] and replaces rational unit lifting by a fixed bounded carrier-power support solve.'}
 with (HERE/'normalized-macaulay.log').open('wb') as fd:
  try:
   proc=subprocess.run([sys.executable,str(Path(__file__)),'--worker'],stdout=fd,stderr=subprocess.STDOUT,timeout=180,preexec_fn=limits);rec.update(outcome='completed',returncode=proc.returncode)
  except subprocess.TimeoutExpired:rec.update(outcome='timeout',returncode=None)
 rec['seconds']=time.monotonic()-START
 for fn in ['normalized-identity-candidate.json','normalized-macaulay-checkpoint.json','normalized-macaulay-search.json']:
  p=HERE/fn
  if p.exists():rec[fn]={'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
 (HERE/'normalized-macaulay-run.json').write_text(json.dumps(rec,indent=2)+'\n');emit(rec);print((HERE/'normalized-macaulay.log').read_text()[-4500:])
