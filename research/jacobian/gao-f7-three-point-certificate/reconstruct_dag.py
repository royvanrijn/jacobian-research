"""Reconstruct local F4 identities; accept each only after exact QQ replay.

The generation backend is untrusted. No modular unit-ideal inference is used.
"""
from pathlib import Path
import collections,gzip,hashlib,json,math,pickle,time
import sympy as S
from flint import fmpq,fmpq_mpoly_ctx,nmod_mpoly_ctx
HERE=Path(__file__).resolve().parent
V=['u','v','w','z','h'];START=time.monotonic()

def parse(path):
 nodes={};idx=None
 for line in (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()):
  if line.startswith('NODE '):idx=int(line.split()[1]);nodes[idx]=[[],[]]
  elif line.startswith('P '):c,e=line[2:].strip().split('|');nodes[idx][0].append((tuple(map(int,e.split(','))),int(c)))
  elif line.startswith('R '):j,c,e=line[2:].strip().split('|');nodes[idx][1].append(((int(j),tuple(map(int,e.split(',')))),int(c)))
 return nodes

def rr(a,modulus):
 bound=math.isqrt(modulus//2);r0,r1=modulus,a;t0,t1=0,1
 while r1>bound:
  q=r0//r1;r0,r1=r1,r0-q*r1;t0,t1=t1,t0-q*t1
 if t1<0:r1,t1=-r1,-t1
 if t1<=bound and math.gcd(r1,t1)==1 and (a*t1-r1)%modulus==0:return (r1,t1)
 return None

def grouped(rel_keys,values,context):
 groups=collections.defaultdict(dict)
 for (parent,ex),c in zip(rel_keys,values):
  if c:groups[parent][ex]=c
 return {j:context.from_dict(f) for j,f in groups.items()}

state_path=HERE/'reconstruction-state.pkl'
if state_path.exists():
 state=pickle.loads(state_path.read_bytes())
else:
 first=HERE/'cubic-six-dag-1000000007.txt'
 if not first.exists():first=first.with_suffix('.txt.gz')
 nodes=parse(first);keep={max(nodes)};pending=list(keep)
 while pending:
  for (j,ex),cf in nodes[pending.pop()][1]:
   if j not in keep:keep.add(j);pending.append(j)
 data={}
 for i in sorted(keep):
  poly,rels=nodes[i];data[i]={'poly_keys':[key for key,c in poly],'rel_keys':[key for key,c in rels],'residues':[c for key,c in poly]+[c for key,c in rels]}
 state={'rows':data,'modulus':1000000007,'primes':[1000000007],'rejected':[]}
 # First replay every retained local identity over its field.
 context=nmod_mpoly_ctx.get(V,1000000007);polys={};local_count=0
 for i,row in data.items():
  n=len(row['poly_keys']);coeffs=row['residues'];poly=context.from_dict(dict(zip(row['poly_keys'],coeffs[:n])))
  rels=grouped(row['rel_keys'],coeffs[n:],context)
  if rels:
   total=context.constant(0)
   for j,mul in rels.items():total+=mul*polys[j]
   assert total==poly,('modular local identity',i)
   local_count+=1
  polys[i]=poly
 assert polys[max(polys)]==context.constant(1)
 print('MODULAR DAG REPLAY PASS',len(data),local_count,flush=True)
 state_path.write_bytes(pickle.dumps(state,protocol=5))

for path in sorted(HERE.glob('cubic-six-dag-*.txt*'),key=lambda p:int(p.name.split('-')[-1].split('.')[0])):
 prime=int(path.name.split('-')[-1].split('.')[0])
 if prime in state['primes'] or prime in state['rejected']:continue
 nodes=parse(path)
 compatible=all(i in nodes and [k for k,c in nodes[i][0]]==row['poly_keys'] and [k for k,c in nodes[i][1]]==row['rel_keys'] for i,row in state['rows'].items())
 if not compatible:
  print('REJECT INCOMPATIBLE PRIME',prime,flush=True);state['rejected'].append(prime);continue
 modulus=state['modulus'];inv=pow(modulus%prime,-1,prime)
 for i,row in state['rows'].items():
  fresh=[c for key,c in nodes[i][0]]+[c for key,c in nodes[i][1]]
  row['residues']=[old+modulus*(((new-old)%prime)*inv%prime) for old,new in zip(row['residues'],fresh)]
 state['modulus']*=prime;state['primes'].append(prime)
 print('CRT',len(state['primes']),'primes',state['modulus'].bit_length(),'bits',flush=True)
 state_path.write_bytes(pickle.dumps(state,protocol=5))

# Reuse only previously proved local QQ identities, not probabilistic candidates.
partial_path=HERE/'partial-rational-dag.json.gz'
verified={}
if partial_path.exists():
 with gzip.open(partial_path,'rt') as fd:verified={row['id']:row for row in json.load(fd)['nodes']}
ctx=fmpq_mpoly_ctx.get(V);polys={};counts=collections.Counter();max_bits=0
source=json.loads((HERE/'cubic-six-input.json').read_text());symbols=list(map(S.Symbol,V))
originals=[]
for row in source['equations']:
 originals.append(ctx.from_dict({ex:fmpq(int(c)) for ex,c in S.Poly(S.sympify(row['polynomial']),*symbols).terms()}))

def load_polynomial(row):return ctx.from_dict({tuple(term['exponents']):fmpq(term['coefficient']) for term in row['polynomial']})
def dump_checkpoint():
 with gzip.open(partial_path,'wt') as fd:json.dump({'variables':V,'nodes':[verified[i] for i in sorted(verified)],'boundary':'Every saved node has a literal QQ local-identity replay; a partial DAG is not a unit certificate.'},fd,separators=(',',':'))

for i,row in state['rows'].items():
 if i in verified:
  polys[i]=load_polynomial(verified[i]);counts['already_verified']+=1;continue
 parents={j for j,e in row['rel_keys']}
 if not parents.issubset(polys):counts['parent_not_yet_proved']+=1;continue
 values=[rr(c,state['modulus']) for c in row['residues']]
 if any(c is None for c in values):counts['reconstruction_unavailable']+=1;continue
 coeffs=[fmpq(a,b) for a,b in values];n=len(row['poly_keys']);poly=ctx.from_dict(dict(zip(row['poly_keys'],coeffs[:n])))
 rels=grouped(row['rel_keys'],coeffs[n:],ctx)
 rec={'id':i,'polynomial':[{'exponents':list(ex),'coefficient':str(c)} for ex,c in zip(row['poly_keys'],coeffs[:n]) if c]}
 if parents:
  total=ctx.constant(0)
  for j,mul in rels.items():total+=mul*polys[j]
  if total!=poly:counts['candidate_identity_failed']+=1;continue
  rec['relations']=[{'parent':j,'exponents':list(ex),'coefficient':str(c)} for (j,ex),c in zip(row['rel_keys'],coeffs[n:]) if c]
 else:
  choices=[]
  for k,original in enumerate(originals):
   ex=row['poly_keys'][0];scale=poly.to_dict().get(ex,0)/original.to_dict().get(ex,0) if original.to_dict().get(ex,0) else None
   if scale is not None and scale*original==poly:choices.append((k,scale))
  assert len(choices)==1,('initial polynomial',i,choices)
  rec['source_generator'],scale=choices[0];rec['source_scalar']=str(scale)
 polys[i]=poly;verified[i]=rec;counts['newly_verified']+=1
 max_bits=max(max_bits,max(max(abs(a).bit_length(),b.bit_length()) for a,b in values))
 if len(verified)%25==0:
  print('VERIFIED',len(verified),'of',len(state['rows']),'latest',i,'bits',max_bits,flush=True);dump_checkpoint()

dump_checkpoint()
unit=max(state['rows']);complete=unit in polys and polys[unit]==ctx.constant(1)
report={'primes':state['primes'],'rejected_primes':state['rejected'],'modulus_bits':state['modulus'].bit_length(),'needed_nodes':len(state['rows']),'proved_rational_nodes':len(verified),'counts':dict(counts),'max_new_coefficient_bits':max_bits,'unit_proved':complete,'seconds':time.monotonic()-START,'boundary':'Only a complete acyclic chain of exact QQ identities ending in1 excludes the coefficient family.'}
(HERE/'reconstruction-progress.json').write_text(json.dumps(report,indent=2)+'\n')
with (HERE/'reconstruction-history.jsonl').open('a') as fd:fd.write(json.dumps(report)+'\n')
if complete:
 certificate={'variables':V,'source':'cubic-six-input.json','unit_node':unit,'nodes':[verified[i] for i in sorted(verified)]}
 with gzip.open(HERE/'cubic-six-rational-dag.json.gz','wt') as fd:json.dump(certificate,fd,separators=(',',':'))
print(json.dumps(report),flush=True)
