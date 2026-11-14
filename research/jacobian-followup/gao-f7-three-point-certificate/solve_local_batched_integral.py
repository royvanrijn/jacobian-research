"""Recover local rational identities by shared exact multi-RHS solves.

All accepted nodes are literal QQ identities, rooted in earlier proved nodes.
No modular inference establishes a theorem. The final checker is separate.
"""
from pathlib import Path
import collections,gzip,json,pickle,time,os,resource,math
from flint import fmpq,fmpq_mat,fmpq_mpoly_ctx,nmod_mat
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
HERE=Path(__file__).resolve().parent;START=time.monotonic();V=['u','v','w','z','h'];ctx=fmpq_mpoly_ctx.get(V)
state=pickle.loads((HERE/'reconstruction-state.pkl').read_bytes());modulus=state['modulus'];prime=state['primes'][0]
partial=HERE/'partial-rational-dag.json.gz'
with gzip.open(partial,'rt') as fd:verified={row['id']:row for row in json.load(fd)['nodes']}
for path in sorted(HERE.glob('local-proof-journal-*.jsonl')):
 for line in path.open():
  try:record=json.loads(line)
  except json.JSONDecodeError:continue
  verified[record['id']]=record
polys={i:ctx.from_dict({tuple(t['exponents']):fmpq(t['coefficient']) for t in row['polynomial']}) for i,row in verified.items()}
journal=HERE/('local-proof-journal-'+str(time.time_ns())+'.jsonl');reports=[]
def cfmod(c,p):return c%p if isinstance(c,int) else int(c.p)%p*pow(int(c.q)%p,-1,p)%p
def pivots(matrix,rank):return [next(j for j in range(matrix.ncols()) if matrix[i,j]) for i in range(rank)]
def checkpoint():
 (HERE/'batched-solve-progress.json').write_text(json.dumps({'proved_nodes':len(verified),'needed_nodes':len(state['rows']),'seconds':time.monotonic()-START,'journal':journal.name,'reports':reports},indent=2)+'\n')
def solve_group(ids):
 start=time.monotonic();keys=state['rows'][ids[0]]['rel_keys'];columns=[];allmonoms=set()
 scales=[]
 for parent,mult in keys:
  parent_terms=polys[parent].to_dict()
  denominator=math.lcm(*(int(c.q) for c in parent_terms.values()))
  col={tuple(int(a)+b for a,b in zip(ex,mult)):int(c.p)*(denominator//int(c.q)) for ex,c in parent_terms.items()}
  scales.append(denominator)
  columns.append(col);allmonoms.update(col)
 constraints={i:(allmonoms-set(state['rows'][i]['poly_keys']))|{state['rows'][i]['poly_keys'][0]} for i in ids}
 common=set.intersection(*(constraints[i] for i in ids));rows=sorted(common);index={e:j for j,e in enumerate(rows)};k=len(columns)
 report={'nodes':ids,'rows':len(rows),'columns':k,'rhs':len(ids)}
 if any(state['rows'][i]['poly_keys'][0] not in common for i in ids):report['outcome']='common support omits a required pivot';reports.append(report);return 0
 if len(rows)*(k+len(ids))>12_000_000:report['outcome']='matrix-entry cap';reports.append(report);return 0
 matrix=nmod_mat(len(rows),k+len(ids),prime)
 for j,col in enumerate(columns):
  for e,c in col.items():
   if e in index:matrix[index[e],j]=cfmod(c,prime)
 for j,i in enumerate(ids):matrix[index[state['rows'][i]['poly_keys'][0]],k+j]=1
 reduced,rank=matrix.rref(inplace=True);pv=pivots(reduced,rank);chosen=[j for j in pv if j<k]
 if any(j>=k for j in pv):report['outcome']='inconsistent modular support';reports.append(report);return 0
 report['exact_dimension']=len(chosen)
 if len(chosen)>800:report['outcome']='exact dimension cap';reports.append(report);return 0
 transposed=nmod_mat(len(chosen),len(rows),prime)
 for j,colid in enumerate(chosen):
  for e,c in columns[colid].items():
   if e in index:transposed[j,index[e]]=cfmod(c,prime)
 rr,rk=transposed.rref(inplace=True);selected=pivots(rr,rk);assert rk==len(chosen)
 square=fmpq_mat([[columns[j].get(rows[r],0) for j in chosen] for r in selected])
 rhs=fmpq_mat([[int(rows[r]==state['rows'][i]['poly_keys'][0]) for i in ids] for r in selected])
 del matrix,reduced,transposed,rr
 (HERE/'batched-solve-active.json').write_text(json.dumps(report)+'\n')
 print('EXACT BATCH',report,flush=True)
 solution=square.solve(rhs,algorithm="dixon");accepted=0
 for z,i in enumerate(ids):
  row=state['rows'][i];total=ctx.constant(0);relations=[]
  for j,colid in enumerate(chosen):
   c=solution[j,z]*scales[colid]
   if c:
    parent,e=keys[colid];total+=ctx.from_dict({e:c})*polys[parent]
    relations.append({'parent':parent,'exponents':list(e),'coefficient':str(c)})
  actual=total.to_dict();support=set(row['poly_keys']);leading=row['poly_keys'][0]
  if not set(actual).issubset(support) or actual.get(leading)!=1:continue
  try:consistent=all(cfmod(actual.get(e,fmpq(0)),modulus)==c for e,c in zip(row['poly_keys'],row['residues']))
  except ValueError:consistent=False
  if not consistent:continue
  rec={'id':i,'polynomial':[{'exponents':list(map(int,e)),'coefficient':str(c)} for e,c in actual.items() if c],'relations':relations}
  verified[i]=rec;polys[i]=total;accepted+=1
  with journal.open('a') as fd:fd.write(json.dumps(rec,separators=(',',':'))+'\n')
 report['accepted']=accepted;report['seconds']=time.monotonic()-start;report['outcome']='exact local identities replayed'
 reports.append(report);checkpoint();print('VERIFIED',len(verified),'of',len(state['rows']),'batch_seconds',report['seconds'],flush=True)
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
checkpoint()
print(json.dumps({'proved_nodes':len(verified),'needed_nodes':len(state['rows']),'unit_proved':complete,'seconds':time.monotonic()-START}),flush=True)
