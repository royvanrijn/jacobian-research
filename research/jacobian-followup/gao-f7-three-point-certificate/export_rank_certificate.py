"""Export the proved predecessor DAG and a full coefficient-space rank witness.

The independent checker must replay the DAG before using the rank argument.
A square matrix covering every actual coefficient row is essential.
"""
from pathlib import Path
import gzip,json,pickle,re,hashlib,time,math,os
from flint import fmpq,fmpz,nmod_mat
HERE=Path(__file__).resolve().parent;START=time.monotonic();prime=1000000007
assert all(prime%d for d in range(2,math.isqrt(prime)+1))
state=pickle.loads((HERE/'reconstruction-state.pkl').read_bytes());unit=max(state['rows']);needed=set(state['rows'])-{unit}
with gzip.open(HERE/'partial-rational-dag.json.gz','rt') as fd:base={r['id']:r for r in json.load(fd)['nodes']}
locations={}
for path in sorted(HERE.glob('local-proof-journal-*.jsonl')):
 with path.open('rb') as fd:
  while True:
   pos=fd.tell();line=fd.readline()
   if not line:break
   m=re.search(rb'"id":(\d+)',line[:80])
   if m:locations[int(m.group(1))]=(path,pos)
def node(i):
 if i in locations:
  path,pos=locations[i]
  with path.open('rb') as fd:fd.seek(pos);record=json.loads(fd.readline())
 else:record=base[i]
 assert record['id']==i
 return record
assert needed.issubset(set(base)|set(locations))
keys=state['rows'][unit]['rel_keys'];parents={};columns=[];rows=set()
for parent,shift in keys:
 if parent not in parents:
  terms={tuple(t['exponents']):fmpq(t['coefficient']) for t in node(parent)['polynomial']}
  denominator=fmpz(1)
  for c in terms.values():denominator=denominator.lcm(c.q)
  assert denominator>0
  integer_terms={ex:c.p*(denominator//c.q) for ex,c in terms.items()}
  assert all(fmpq(c,denominator)==terms[ex] for ex,c in integer_terms.items())
  parents[parent]=integer_terms
 col={tuple(a+b for a,b in zip(ex,shift)):c for ex,c in parents[parent].items()}
 columns.append(col);rows.update(col)
rows=sorted(rows);assert len(rows)==len(columns)==402 and (0,0,0,0,0) in rows
index={e:j for j,e in enumerate(rows)};matrix=nmod_mat(len(rows),len(columns),prime)
for j,col in enumerate(columns):
 for e,c in col.items():matrix[index[e],j]=int(c%prime)
determinant=int(matrix.det());assert determinant!=0
terminal={'kind':'full_integer_coefficient_matrix','prime':prime,'expected_shape':[len(rows),len(columns)],'expected_det_mod':determinant,'columns':[{'parent':j,'exponents':list(e)} for j,e in keys]}
report={'status':'TERMINAL_RANK_CANDIDATE','dag_replay':'PENDING','parent_nodes':len(needed),'shape':[len(rows),len(columns)],'prime':prime,'det_mod':determinant,'seconds':time.monotonic()-START,'lemma':'Every column is a nonzero integer scaling of a monomial times a proved ideal element. The matrix covers ALL actual coefficient monomials including1. A nonzero determinant modulo a proven prime implies nonzero integer determinant, hence full span over Q. A modular unit ideal alone would not imply this.'}
(HERE/'terminal-rank-generation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report),flush=True)
# Stream the existing complete node objects; do not allocate another giant DAG.
header={'variables':['u','v','w','z','h'],'source':'cubic-six-input.json','terminal':terminal}
temporary=HERE/'rank-certificate.tmp';target=HERE/'cubic-six-rank-dag.json.gz'
with gzip.open(temporary,'wt',compresslevel=1) as fd:
 prefix=json.dumps(header,separators=(',',':'));fd.write(prefix[:-1]+',"nodes":[')
 for count,i in enumerate(sorted(needed)):
  if count:fd.write(',')
  json.dump(node(i),fd,separators=(',',':'))
 fd.write(']}')
os.replace(temporary,target)
report['certificate_sha256']=hashlib.sha256(target.read_bytes()).hexdigest();report['certificate_bytes']=target.stat().st_size;report['seconds']=time.monotonic()-START
(HERE/'terminal-rank-generation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
