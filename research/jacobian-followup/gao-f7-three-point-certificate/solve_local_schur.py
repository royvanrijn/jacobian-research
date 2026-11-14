"""Exact triangular elimination before small support-guided rational solves.

Upper columns have distinct leading monomials. Their triangular block is
eliminated exactly, reducing the dense Schur system. Every exported node
is still checked as a literal combination of the original parent columns.
"""
from pathlib import Path
import collections,gzip,json,pickle,time,os,resource,math,concurrent.futures,multiprocessing
from flint import fmpq,fmpq_mat,fmpq_mpoly_ctx,nmod_mat
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
HERE=Path(__file__).resolve().parent;START=time.monotonic();V=['u','v','w','z','h'];ctx=fmpq_mpoly_ctx.get(V)
state=pickle.loads((HERE/'reconstruction-state.pkl').read_bytes());modulus=state['modulus'];prime=state['primes'][0]
partial=HERE/'partial-rational-dag.json.gz'
with gzip.open(partial,'rt') as fd:verified={r['id']:r for r in json.load(fd)['nodes']}
for path in sorted(HERE.glob('local-proof-journal-*.jsonl')):
 for line in path.open():
  try:record=json.loads(line)
  except json.JSONDecodeError:continue
  verified[record['id']]=record
polys={i:ctx.from_dict({tuple(t['exponents']):fmpq(t['coefficient']) for t in row['polynomial']}) for i,row in verified.items()}
journal=HERE/('local-proof-journal-'+str(time.time_ns())+'.jsonl');reports=[]
def cfmod(c,p):return c%p if isinstance(c,int) else int(c.p)%p*pow(int(c.q)%p,-1,p)%p
def order(e):return (sum(e),tuple(-x for x in reversed(e)))
def pivots(matrix,rank):return [next(j for j in range(matrix.ncols()) if matrix[i,j]) for i in range(rank)]
def checkpoint():
 (HERE/'schur-solve-progress.json').write_text(json.dumps({'proved_nodes':len(verified),'needed_nodes':len(state['rows']),'seconds':time.monotonic()-START,'journal':journal.name,'reports':reports},indent=2)+'\n')
def normal(poly):
 used={}
 for u in UPPER:
  head=HEADS[u];c=poly[head]
  if c:
   c/=COLS[u][head];poly-=c*COLS[u];used[u]=c
 assert all(poly[HEADS[u]]==0 for u in UPPER)
 return poly,used

def solve_one(i):
 row=state['rows'][i];leading=row['poly_keys'][0]
 b,beta=TARGETS[i]
 rhs=fmpq_mat([[b[SMALL_ROWS[r]]] for r in SELECTED])
 solution=SMALL_MATRIX.solve(rhs,algorithm='dixon')
 lam=[fmpq(0) for c in COLS]
 for u,c in beta.items():lam[u]=c
 for j,colid in enumerate(CHOSEN):
  mu=solution[j,0]*SCALES[colid];original_column=LOWER[colid]
  lam[original_column]+=mu
  for u,c in GAMMA[colid].items():lam[u]-=mu*c
 total=ctx.constant(0);relations=[]
 for j,c in enumerate(lam):
  if c:
   total+=c*COLS[j];parent,e=KEYS[j]
   relations.append({'parent':parent,'exponents':list(e),'coefficient':str(c)})
 actual=total.to_dict()
 if not set(actual).issubset(set(row['poly_keys'])) or actual.get(leading)!=1:return None
 try:consistent=all(cfmod(actual.get(e,fmpq(0)),modulus)==c for e,c in zip(row['poly_keys'],row['residues']))
 except ValueError:consistent=False
 if not consistent:return None
 return {'id':i,'polynomial':[{'exponents':list(map(int,e)),'coefficient':str(c)} for e,c in actual.items() if c],'relations':relations}

def solve_group(ids):
 global KEYS,COLS,HEADS,UPPER,LOWER,GAMMA,TARGETS,SMALL_ROWS,SELECTED,CHOSEN,SCALES,SMALL_MATRIX
 start=time.monotonic();KEYS=state['rows'][ids[0]]['rel_keys'];COLS=[];allmonoms=set()
 for parent,e in KEYS:
  poly=ctx.from_dict({e:1})*polys[parent];COLS.append(poly);allmonoms.update(tuple(map(int,m)) for m in poly.to_dict())
 constraints={i:(allmonoms-set(state['rows'][i]['poly_keys']))|{state['rows'][i]['poly_keys'][0]} for i in ids}
 common=set.intersection(*(constraints[i] for i in ids))
 HEADS=[max((tuple(map(int,e)) for e in poly.to_dict()),key=order) for poly in COLS]
 chosen_heads={}
 for j,head in enumerate(HEADS):
  if head in common and head not in chosen_heads:chosen_heads[head]=j
 UPPER=sorted(chosen_heads.values(),key=lambda j:order(HEADS[j]),reverse=True)
 upper_set=set(UPPER);LOWER=[j for j in range(len(COLS)) if j not in upper_set]
 SMALL_ROWS=sorted(common-set(chosen_heads));index={e:j for j,e in enumerate(SMALL_ROWS)}
 report={'nodes':ids,'original_columns':len(COLS),'triangular_columns':len(UPPER),'small_columns':len(LOWER),'small_rows':len(SMALL_ROWS)}
 print('TRIANGULAR REDUCTION',report,flush=True)
 residuals=[];GAMMA=[]
 for j in LOWER:
  residual,used=normal(COLS[j]);residuals.append(residual);GAMMA.append(used)
 TARGETS={i:normal(ctx.from_dict({state['rows'][i]['poly_keys'][0]:1})) for i in ids}
 integer_columns=[];SCALES=[]
 for residual in residuals:
  terms=residual.to_dict();d=math.lcm(*(int(c.q) for c in terms.values())) if terms else 1
  integer_columns.append({tuple(map(int,e)):int(c.p)*(d//int(c.q)) for e,c in terms.items()});SCALES.append(d)
 k=len(LOWER)
 if not k:report['outcome']='no lower columns';reports.append(report);return 0
 matrix=nmod_mat(len(SMALL_ROWS),k+len(ids),prime)
 for j,col in enumerate(integer_columns):
  for e,c in col.items():
   if e in index:matrix[index[e],j]=c%prime
 for j,i in enumerate(ids):
  for e,c in TARGETS[i][0].to_dict().items():
   if e in index:matrix[index[e],k+j]=cfmod(c,prime)
 reduced,rank=matrix.rref(inplace=True);pv=pivots(reduced,rank);CHOSEN=[j for j in pv if j<k]
 if any(j>=k for j in pv):report['outcome']='inconsistent modular Schur support';reports.append(report);return 0
 transposed=nmod_mat(len(CHOSEN),len(SMALL_ROWS),prime)
 for j,colid in enumerate(CHOSEN):
  for e,c in integer_columns[colid].items():
   if e in index:transposed[j,index[e]]=c%prime
 rr,rk=transposed.rref(inplace=True);SELECTED=pivots(rr,rk);assert rk==len(CHOSEN)
 SMALL_MATRIX=fmpq_mat([[integer_columns[j].get(SMALL_ROWS[r],0) for j in CHOSEN] for r in SELECTED])
 report['exact_dimension']=len(CHOSEN);report['reduction_seconds']=time.monotonic()-start
 print('SMALL EXACT SYSTEM',report,flush=True)
 (HERE/'schur-solve-active.json').write_text(json.dumps(report)+'\n')
 accepted=0
 with concurrent.futures.ProcessPoolExecutor(max_workers=min(8,len(ids)),mp_context=multiprocessing.get_context('fork')) as executor:
  futures={executor.submit(solve_one,i):i for i in ids}
  for future in concurrent.futures.as_completed(futures):
   rec=future.result()
   if rec is None:continue
   i=rec['id'];verified[i]=rec;polys[i]=ctx.from_dict({tuple(t['exponents']):fmpq(t['coefficient']) for t in rec['polynomial']});accepted+=1
   with journal.open('a') as fd:fd.write(json.dumps(rec,separators=(',',':'))+'\n')
   checkpoint();print('VERIFIED',len(verified),'of',len(state['rows']),'node',i,flush=True)
 report['accepted']=accepted;report['seconds']=time.monotonic()-start;report['outcome']='exact local identities replayed';reports.append(report);checkpoint()
 return accepted

while len(verified)<len(state['rows']):
 groups=collections.defaultdict(list)
 for i,row in state['rows'].items():
  if i not in verified and all(j in polys for j,e in row['rel_keys']):groups[tuple(row['rel_keys'])].append(i)
 if not groups:break
 changed=0
 for ids in groups.values():
  added=solve_group(ids);changed+=added
  if added==0 and len(ids)>1:
   for i in ids:changed+=solve_group([i])
 if not changed:break
unit=max(state['rows']);complete=unit in polys and polys[unit]==ctx.constant(1)
if complete:
 certificate={'variables':V,'source':'cubic-six-input.json','unit_node':unit,'nodes':[verified[i] for i in sorted(verified)]}
 temporary=HERE/'complete-certificate.tmp'
 with gzip.open(temporary,'wt',compresslevel=1) as fd:json.dump(certificate,fd,separators=(',',':'))
 os.replace(temporary,HERE/'cubic-six-rational-dag.json.gz')
checkpoint();print(json.dumps({'proved_nodes':len(verified),'unit_proved':complete,'seconds':time.monotonic()-START}),flush=True)
