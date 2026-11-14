from pathlib import Path
import gzip,json,re,time
from flint import fmpq,fmpz
HERE=Path(__file__).resolve().parent;started=time.monotonic()
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
  with path.open('rb') as fd:fd.seek(pos);return json.loads(fd.readline())
 return base[i]
denominators={}
def denominator(i):
 if i not in denominators:
  d=fmpz(1)
  for t in node(i)['polynomial']:d=d.lcm(fmpq(t['coefficient']).q)
  denominators[i]=d
 return denominators[i]
i=593;r=node(i);di=denominator(i);alpha=[fmpq(t['coefficient'])*fmpq(di,denominator(t['parent'])) for t in r['relations']]
D=fmpz(1)
for a in alpha:D=D.lcm(a.q)
integers=[]
for a in alpha:
 value=a*D;assert value.q==1;integers.append(value.p)
report={'node':i,'relation_terms':len(alpha),'rational_coefficient_characters':sum(len(t['coefficient']) for t in r['relations']),'integer_coefficient_characters':sum(len(str(n)) for n in integers)+len(str(D)),'divisor_digits':len(str(D)),'max_integer_digits':max(len(str(n)) for n in integers),'seconds':time.monotonic()-started,'identity':'G_i=(sum(N_ij*monomial*G_j))/D_i, where G_j=d_j*g_j and d_j clears all polynomial coefficient denominators.'}
(HERE/'integer-encoding-sample.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
