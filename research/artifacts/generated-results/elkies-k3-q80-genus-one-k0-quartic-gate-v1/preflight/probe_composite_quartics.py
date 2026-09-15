"""Complete composite quartic residue signatures from retained rational/quadratic data."""
import itertools, json, resource, time
from pathlib import Path
from collections import defaultdict
resource.setrlimit(resource.RLIMIT_CPU,(30,35));resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
R=Path(__file__).resolve().parents[4];base=R/'artifacts/generated-results/elkies-k3-q80-degree-two-reciprocity-v1'
rat=json.loads((base/'rational-fibres.json').read_text())['rows'];quad=json.loads((base/'quadratic-fibres.json').read_text())['rows']
start=time.process_time()
rs=[r for r in rat if len(r['rational_roots'])==3]
qs=[r for r in quad if len(r['roots'])==3]
def perm(r,key):
 return set(itertools.permutations(zip(r['roots'],r[key])))
def sig(c):return tuple(sorted(c))
rpairs=[];rindex=defaultdict(list)
for a,b in itertools.combinations(rs,2):
 for aa in perm(a,'rational_codes'):
  for bb in perm(b,'rational_codes'):
   codes=[aa[i][1]^bb[i][1] for i in range(3)]
   row={'bases':[a['t'],b['t']],'roots':[[e for e,c in aa],[e for e,c in bb]],'codes':codes}
   rpairs.append(row);rindex[sig(codes)].append(row)
four=[]
for rows in rindex.values():
 for a,b in itertools.combinations(rows,2):
  if set(a['bases']).isdisjoint(b['bases']):
   four.append([a,b])
mixed=[]
for q in qs:
 key=sig(q['norm_codes'])
 for rr in rindex.get(key,[]):
  mixed.append({'rational':rr,'quadratic':q['t'],'codes':q['norm_codes'],'smooth':q['smooth']})
qindex=defaultdict(list)
for q in qs:qindex[sig(q['norm_codes'])].append(q)
two=[]
for rows in qindex.values():
 for a,b in itertools.combinations(rows,2):
  if a['t']!=b['t']:two.append({'quadratics':[a['t'],b['t']],'codes':[a['norm_codes'],b['norm_codes']]})
print(json.dumps({'rational_pair_choices':len(rpairs),'four_distinct_rational_matches':four,
                  'two_rational_one_quadratic_matches':mixed,'two_distinct_quadratic_matches':two,
                  'rational_signatures':[list(sig(r['rational_codes'])) for r in rs],
                  'cpu_seconds':time.process_time()-start},indent=2))
