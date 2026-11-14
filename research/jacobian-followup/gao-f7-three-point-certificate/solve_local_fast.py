"""Support-guided exact local solves for an acyclic polynomial certificate.

Modular data selects rows/columns only. A node is accepted solely after
literal QQ identity replay and full-support checks. Earlier proved nodes
are retained; the final verifier independently replays the complete chain.
"""
from pathlib import Path
import collections,gzip,json,pickle,time,os
from flint import fmpq,fmpq_mat,fmpq_mpoly_ctx,nmod_mat
HERE=Path(__file__).resolve().parent
START=time.monotonic();V=['u','v','w','z','h'];ctx=fmpq_mpoly_ctx.get(V)
state=pickle.loads((HERE/'reconstruction-state.pkl').read_bytes());modulus=state['modulus'];prime=state['primes'][0]
partial=HERE/'partial-rational-dag.json.gz'
with gzip.open(partial,'rt') as fd:verified={row['id']:row for row in json.load(fd)['nodes']}
# Each completed line is a previously replayed local identity. An interrupted
# final line is ignored; final verification independently replays everything.
for journal_path in sorted(HERE.glob('local-proof-journal-*.jsonl')):
 with journal_path.open() as fd:
  for line in fd:
   try:record=json.loads(line)
   except json.JSONDecodeError:continue
   verified[record['id']]=record
journal_path=HERE/('local-proof-journal-'+str(time.time_ns())+'.jsonl')

polys={i:ctx.from_dict({tuple(t['exponents']):fmpq(t['coefficient']) for t in row['polynomial']}) for i,row in verified.items()}
reports=[]

def checkpoint():
 (HERE/'local-solve-progress.json').write_text(json.dumps({'proved_nodes':len(verified),'needed_nodes':len(state['rows']),'seconds':time.monotonic()-START,'journal':journal_path.name,'reports':reports},indent=2)+'\n')

def materialize():
 temporary=HERE/'partial-local.tmp'
 with gzip.open(temporary,'wt',compresslevel=1) as fd:json.dump({'variables':V,'nodes':[verified[i] for i in sorted(verified)],'boundary':'Exact local identities only; partial DAG is not a unit proof.'},fd,separators=(',',':'))
 os.replace(temporary,partial)

def pivots(matrix,rank):return [next(j for j in range(matrix.ncols()) if matrix[i,j]) for i in range(rank)]

def cfmod(value,p):return int(value.p)%p*pow(int(value.q)%p,-1,p)%p

for i,row in state['rows'].items():
 if i in verified:continue
 parents={j for j,e in row['rel_keys']}
 if not parents.issubset(polys):continue
 start=time.monotonic();columns=[];monomials=set();support=set(row['poly_keys']);leading=row['poly_keys'][0]
 for parent,mult in row['rel_keys']:
  col={tuple(a+b for a,b in zip(ex,mult)):cf for ex,cf in polys[parent].to_dict().items()}
  columns.append(col);monomials.update(ex for ex in col if ex not in support)
 monomials.add(leading);rows=sorted(monomials);index={ex:j for j,ex in enumerate(rows)};k=len(columns)
 report={'node':i,'rows':len(rows),'columns':k}
 if len(rows)*(k+1)>12_000_000:report['outcome']='matrix-entry cap';reports.append(report);continue
 matrix=nmod_mat(len(rows),k+1,prime)
 for j,col in enumerate(columns):
  for ex,cf in col.items():
   if ex in index:matrix[index[ex],j]=cfmod(cf,prime)
 matrix[index[leading],k]=1
 reduced,rank=matrix.rref(inplace=True);pv=pivots(reduced,rank)
 if k in pv:
  report['outcome']='inconsistent modular support';reports.append(report);continue
 chosen=[pv[j] for j in range(rank) if reduced[j,k]]
 report['exact_dimension']=len(chosen)
 if len(chosen)>800:report['outcome']='exact dimension cap';reports.append(report);continue
 transposed=nmod_mat(len(chosen),len(rows),prime)
 for j,colid in enumerate(chosen):
  for ex,cf in columns[colid].items():
   if ex in index:transposed[j,index[ex]]=cfmod(cf,prime)
 rr,rk=transposed.rref(inplace=True);selected=pivots(rr,rk)
 assert rk==len(chosen)
 square=fmpq_mat([[columns[j].get(rows[r],0) for j in chosen] for r in selected])
 rhs=fmpq_mat([[int(rows[r]==leading)] for r in selected])
 del matrix,reduced,transposed,rr
 (HERE/'local-solve-active.json').write_text(json.dumps(report)+'\n')
 print('EXACT LOCAL SOLVE',report,flush=True)
 solution=square.solve(rhs)
 total=ctx.constant(0);relations=[]
 for j,colid in enumerate(chosen):
  c=solution[j,0]
  if c:
   parent,ex=row['rel_keys'][colid]
   total+=ctx.from_dict({ex:c})*polys[parent]
   relations.append({'parent':parent,'exponents':list(ex),'coefficient':str(c)})
 actual=total.to_dict()
 if not set(actual).issubset(support) or actual.get(leading)!=1:
  report['outcome']='exact constraints failed';reports.append(report);checkpoint();continue
 # This congruence is a consistency check, not the proof of ideal membership.
 # The proof is the preceding literal rational polynomial sum.
 try:consistent=all(cfmod(actual.get(ex,fmpq(0)),modulus)==c for ex,c in zip(row['poly_keys'],row['residues']))
 except ValueError:consistent=False
 if not consistent:
  report['outcome']='exact polynomial differs from retained modular row';reports.append(report);checkpoint();continue
 rec={'id':i,'polynomial':[{'exponents':list(map(int,ex)),'coefficient':str(cf)} for ex,cf in actual.items() if cf],'relations':relations}
 # Independent evaluation from the serialized local multipliers.
 replay=ctx.constant(0)
 groups=collections.defaultdict(dict)
 for relation in rec['relations']:groups[relation['parent']][tuple(relation['exponents'])]=fmpq(relation['coefficient'])
 for parent,terms in groups.items():replay+=ctx.from_dict(terms)*polys[parent]
 assert replay==total
 verified[i]=rec;polys[i]=total
 with journal_path.open('a') as fd:fd.write(json.dumps(rec,separators=(',',':'))+'\n')
 report['outcome']='EXACT LOCAL IDENTITY PASS';report['seconds']=time.monotonic()-start;reports.append(report)
 print('VERIFIED',len(verified),'of',len(state['rows']),'node',i,'seconds',report['seconds'],flush=True)
 checkpoint()
unit=max(state['rows']);complete=unit in polys and polys[unit]==ctx.constant(1)
if complete:
 certificate={'variables':V,'source':'cubic-six-input.json','unit_node':unit,'nodes':[verified[i] for i in sorted(verified)]}
 temporary=HERE/'complete-certificate.tmp'
 with gzip.open(temporary,'wt',compresslevel=1) as fd:json.dump(certificate,fd,separators=(',',':'))
 os.replace(temporary,HERE/'cubic-six-rational-dag.json.gz')
checkpoint()
materialize()
print(json.dumps({'proved_nodes':len(verified),'needed_nodes':len(state['rows']),'unit_proved':complete,'seconds':time.monotonic()-START}),flush=True)
