import sys,json,resource,time
from pathlib import Path
from fractions import Fraction as F
resource.setrlimit(resource.RLIMIT_CPU,(20,20));resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
sys.path.insert(0,'research/elkies-k3/scripts')
from verify_q80_branch_trace_specialization import add
from audit_q80_branch_trace_specialization import prime
root=Path('research');src=root/'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json';s=json.load(open(src));b=json.load(open(root/'artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisections-cheapest-1024-v1.json'))['bisections'][0];q=list(map(int,b['branch']['numerator_coefficients']))
def ev(cs,t,p):
 v=0
 for c in reversed(cs):
  z=F(c);v=(v*t+z.numerator*pow(z.denominator,-1,p))%p
 return v
def quotient(a,b,p):
 sq={}
 for y in range(p):sq.setdefault(y*y%p,[]).append(y)
 pts=[None]+[(x,y) for x in range(p) for y in sq.get((x**3+a*x+b)%p,[])]
 triple={add(add(P,P,a,p),P,a,p) for P in pts};labels={P:0 for P in triple};dim=0
 while len(labels)<len(pts):
  P=next(P for P in pts if P not in labels);old=list(labels.items());pp=add(P,P,a,p)
  for Q,v in old:
   for k,R in [(1,P),(2,pp)]:
    z=add(Q,R,a,p);assert z not in labels;labels[z]=v+k*3**dim
  dim+=1
 assert dim<=2 and len(pts)==len(triple)*3**dim
 return labels,dim,len(pts),len(triple)
piv={};records=[];torsion=None
for p in range(5,998):
 if not prime(p) or q[2]%p==0 or (q[1]**2-4*q[0]*q[2])%p==0:continue
 for t in range(p):
  if ev(q,t,p):continue
  try:
   a=ev(s['weierstrass_model']['A_coefficients_low_to_high'],t,p);c=ev(s['weierstrass_model']['B_coefficients_low_to_high'],t,p)
   pts=[]
   for sec in s['sections']['records']:
    pts.append(tuple(ev(sec[k]['numerator_coefficients_low_to_high'],t,p)*pow(ev(sec[k]['denominator_coefficients_low_to_high'],t,p),-1,p)%p for k in ['X','Y']))
  except ValueError:continue
  if (4*a**3+27*c*c)%p==0:continue
  assert all((y*y-x*x*x-a*x-c)%p==0 for x,y in pts)
  labels,dim,n,nt=quotient(a,c,p)
  if n%3 and torsion is None:torsion={'prime':p,'t':t,'order':n}
  rows=[[labels[P]//3**i%3 for P in pts] for i in range(dim)]
  grew=False
  for row in rows:
   z=row[:]
   for i in range(17):
    if not z[i]:continue
    if i in piv:
     h=z[i];z=[(x-h*y)%3 for x,y in zip(z,piv[i])]
    else:
     h=pow(z[i],-1,3);piv[i]=[h*x%3 for x in z];grew=True;break
  if grew:
   records.append({'prime':p,'t':t,'A':a,'B':c,'points':pts,'order':n,'triple_subgroup_order':nt,'rows':rows,'rank':len(piv)})
   print(p,t,len(piv),flush=True)
  if len(piv)==17:break
 if len(piv)==17:break
out={'rank':len(piv),'torsion_witness':torsion,'records':records,'branch':q}
Path('/tmp/q80-branch-mod3.json').write_text(json.dumps(out,indent=2)+'\n');print('RESULT',len(piv),torsion)
