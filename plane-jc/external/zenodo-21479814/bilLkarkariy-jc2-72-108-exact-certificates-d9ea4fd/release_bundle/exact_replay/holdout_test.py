import sympy as sp,pickle,subprocess,sys,re
# find next good prime
S=pickle.load(open('hne0_polred.pkl','rb'));dens={d for e in S for c in e.values() for n,d in c.values()};w=sp.symbols('w');F=w**5-w**4+3*w**3+3*w**2+26
p=500080573
while True:
 p=int(sp.nextprime(p))
 if any(d%p==0 for d in dens):continue
 if sp.Poly(F,w,modulus=p).is_irreducible:break
print('p',p)
subprocess.run([sys.executable,'build_hcert_syz_prime.py',str(p)],check=True)
with open(f'hcert_syz_p{p}.out','w') as f:subprocess.run(['Singular','-q',f'hcert_syz_p{p}.sing'],stdout=f,stderr=subprocess.STDOUT,check=True)
# parser
pat=re.compile(r'^C\|(\d+)\|(\d+)\|(\d+)\|(\d+)\|\((.*)\)$')
def pp(s):
 s=s.replace(' ','');o=[0]*5
 for t in re.findall(r'[+-]?[^+-]+',s):
  if '*w' in t:a,r=t.split('*w',1);e=int(r[1:]) if r.startswith('^') else 1
  elif 'w' in t:
   m=re.fullmatch(r'([+-]?)w(?:\^(\d+))?',t);a='-1' if m.group(1)=='-' else '1';e=int(m.group(2) or 1)
  else:a=t;e=0
  o[e]=int(a)%p
 return o
d=[]
for l in open(f'hcert_syz_p{p}.out'):
 m=pat.match(l)
 if m:d.append((tuple(map(int,m.group(1,2,3,4))),pp(m.group(5))))
base,res=pickle.load(open('rr_cache.pkl','rb'));print('terms',len(d),'same support',[k for k,a in d]==base)
vals=[v for k,a in d for v in a];good=bad=0;examples=[]
for i,z in enumerate(res):
 if z is None:continue
 n,q=z;v=n%p*pow(q%p,-1,p)%p
 if v==vals[i]:good+=1
 else:
  bad+=1
  if len(examples)<5:examples.append((i,z,v,vals[i]))
print('reconstructed check good',good,'bad',bad,'examples',examples)
