from pathlib import Path
from itertools import combinations
from time import monotonic
import json,signal,sys
import sympy as sy
root=Path(__file__).resolve().parents[3]
p=root/'research/jacobian-followup/moment-vanishing-audit';sys.path.insert(0,str(p))
from verify_degree12_sparse_supports import Model

size=int(sys.argv[1]) if len(sys.argv)>1 else 5
start=monotonic();budget=45;cases={}
def alarm(signum,frame):raise TimeoutError()
signal.signal(signal.SIGALRM,alarm)
for name in ['2O','2I']:
 model=Model(name);variables=sy.symbols('t0:'+str(size-2));z=sy.Symbol('z')
 saturation=z*sy.prod(variables)-1
 supports=[q for q in combinations(range(13),size) if q[0]<6<q[-1] and q<=tuple(12-x for x in reversed(q))]
 rows=[]
 for support in supports:
  if monotonic()-start>budget:break
  used=[];status='prefix_unresolved';order=0
  signal.setitimer(signal.ITIMER_REAL,1)
  try:
   for m in range(2,12):
    result=model.moment(support,m)
    if result is None:continue
    val,scalar,poly=result
    used.append(sy.Poly.from_dict({e:sy.Rational(c.numerator,c.denominator) for e,c in poly.items()},variables).as_expr())
    order=m
    if len(poly)==1:status='excluded_monomial';break
    if list(sy.groebner(used+[saturation],z,*variables,order='grevlex'))==[1]:status='excluded_unit';break
  except TimeoutError:status='timeout'
  finally:signal.setitimer(signal.ITIMER_REAL,0)
  rows.append({'support':support,'status':status,'order':order,'polynomials':[str(x) for x in used]})
 cases[name]={'size':size,'total_representatives':len(supports),'rows':rows,
              'unresolved':[r for r in rows if not r['status'].startswith('excluded')]}
 (root/('tmp/degree12_support'+str(size)+'_probe.json')).write_text(json.dumps(cases,indent=2)+'\n')
 print(name,'covered',len(rows),'of',len(supports),'unresolved',len(cases[name]['unresolved']),'seconds',monotonic()-start,flush=True)
 print(cases[name]['unresolved'][:10],flush=True)
