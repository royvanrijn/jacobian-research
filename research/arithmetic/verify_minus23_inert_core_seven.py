#!/usr/bin/env python3
"""Exact phase/onset controls and COMPLETE native principal/cube bridges."""
if not __debug__:raise SystemExit('Refusing -O/-OO: exact assertions required')
from pathlib import Path
from fractions import Fraction as F
from math import isqrt,gcd
from hashlib import sha256
import argparse,json,subprocess,tempfile,mmap,struct,sys
HERE=Path(__file__).resolve().parent
PROOF='MINUS23_INERT_CORE_SEVEN_CLOSURE.md'
INPUT='core_seven_completion_inputs.json'
WITNESS='core_seven_cube_atom_witnesses.json'
NATIVE='verify_core_seven_principal_atoms.cpp'
LEDGER='core_seven_completion_source_ledger.json'

def digest(p):
 h=sha256()
 with p.open('rb')as f:
  while(b:=f.read(1<<20)):h.update(b)
 return h.hexdigest()

def Q(x,y):return x*x+x*y+6*y*y

def values(B):
 out={1:(1,0)}
 for y in range(1,isqrt(4*B//23)+1):
  s=isqrt(4*B-23*y*y)
  for x in range(-((s+y)//2),(s-y)//2+1):
   n=Q(x,y)
   if n>0 and gcd(x,y)==1:out.setdefault(n,(x,y))
 return out

def prime_list(B):
 s=bytearray(b'\1')*(B+1);s[0:2]=b'\0\0'
 for p in range(2,isqrt(B)+1):
  if s[p]:s[p*p::p]=b'\0'*((B-p*p)//p+1)
 return s

def root_path(arg):
 if arg:return arg.resolve()
 for p in [HERE,*HERE.parents]:
  if(p/'research/MATH_STATUS.json').is_file():return p
 raise SystemExit('Use --repo-root for scratch or relocated replay')

def controls(inp):
 assert inp['coefficient']==21 and inp['principal_bound']==2**31-1
 assert inp['cube_prime_cutoff']==1024 and inp['cube_count']==57
 pairs=[((0,0),(1,1)),((0,1),(2,0))]
 fails=[];phase_count=0
 for R in [0,1,2]:
  for a in range(3):
   for b in range(3):
    good=[]
    for pair in pairs:
     good.append(all(R+sum((x+y)%3!=0 for x,y in zip((a,b),v))!=1 for v in pair))
    if not any(good):fails.append((R,a,b))
    phase_count+=1
 assert fails==[(0,0,1),(0,1,2),(0,2,0)]
 for n,xy in [(1,(1,0)),(27,(1,2)),(39,(3,2)),(87,(7,2))]:
  assert Q(*xy)==n and gcd(*xy)==1 and n%23
 assert 1+27==28 and 39+87==126
 selected=set()
 for x in range(-10,11):
  for y in [-1,0,1]:
   n=Q(x,y);g=gcd(x,y)
   if 0<n<=21 and g%2 and g%3:selected.add(n)
 assert selected=={1,6,8,12,18}
 assert not any(21-n in selected for n in selected)
 assert pow(7,11,23)==22 and pow(21,11,23)==22 and pow(3,11,23)==1
 M=F(128,77)*(F(1,2)-F(6,25)+F(1,21))
 assert M==F(20672,40425)
 assert F(28,25)**6<2 and 50**6<2**34
 assert F(87,32)**73<1536**10 and F(87,32)**76<2048**10
 P=M*F(896,25)-(F(64)+F(8000,149))/F(73,10)-F(9,8)
 C=M*50-(F(256,3)+F(80000,999))/F(38,5)-F(7,4)
 assert P==F(13571626381,12562935000)>0 and C==F(84031919,40923036)>0
 assert 21*1024**3>2**34
 return {'phase_cases':phase_count,'displayed_pair_failures':fails,'main':str(M),
    'principal_scaled_surplus':str(P),'cube_scaled_surplus':str(C)}

def cube_check():
 d=json.loads((HERE/WITNESS).read_text());flags=values(1023);pr=prime_list(1023)
 expected=[p for p in range(2,1024)if pr[p] and p!=23 and(p==2 or pow(p%23,11,23)==1)and p not in flags]
 assert len(expected)==57 and d['coefficient']==21 and d['prime_cutoff_exclusive']==1024
 assert d['nonprincipal_primes']==expected and [r['prime']for r in d['witnesses']]==expected
 for row in d['witnesses']:
  assert set(row)=={'prime','coordinates'} and len(row['coordinates'])==4
  x,y,u,v=row['coordinates'];assert all(type(z)is int for z in [x,y,u,v])
  a,b=Q(x,y),Q(u,v)
  assert a>0 and b>0 and a+b==21*row['prime']**3 and gcd(x,y)==gcd(u,v)==1
  assert a%23 and b%23
 return len(expected)

def native_check(inp):
 with tempfile.TemporaryDirectory(prefix='minus23_core_seven_')as tmp:
  p=Path(tmp);exe=p/'finite';pre=p/'run'
  subprocess.run(['g++','-O3','-std=c++17',str(HERE/NATIVE),'-o',str(exe)],check=True,capture_output=True,text=True)
  run=subprocess.run([str(exe),str(inp['principal_bound']),str(pre),'5000'],check=True,capture_output=True,text=True)
  counts=json.loads(Path(str(pre)+'_counts.json').read_text());counts.pop('elapsed_seconds')
  assert counts==inp['principal_expected_counts']
  fp=Path(str(pre)+'_primitive_flags_le64.bin');ap=Path(str(pre)+'_assignments.bin')
  assert digest(fp)==inp['principal_flags_sha256'] and digest(ap)==inp['principal_assignments_sha256']
  assert fp.stat().st_size==268435456 and ap.stat().st_size==8*980245
  # Independently authenticate the complete prime domain and each recorded
  # representation via a byte sieve and the exact coordinate-created flags.
  pr=prime_list(counts['prime_cap']);assert pr.count(1)==counts['all_primes_through_cap']
  with fp.open('rb')as f,mmap.mmap(f.fileno(),0,access=mmap.ACCESS_READ)as flags:
   def has(n):return bool(flags[n>>3]&(1<<(n&7)))
   last=0;records=0
   with ap.open('rb')as af:
    while(chunk:=af.read(8*(1<<14))):
     for pp,anchor in struct.iter_unpack('<II',chunk):
      assert last<pp<=counts['prime_cap'] and pr[pp] and pp!=23 and has(pp)
      assert 1<=anchor<=5000 and has(anchor) and anchor%23
      n=21*pp;assert anchor<n<=inp['principal_bound'] and has(n-anchor) and(n-anchor)%23
      last=pp;records+=1
   assert records==980245
   # Any principal split prime in the independently sieved complete domain
   # must occur among the strictly increasing assignment records.
   principal=0;idx=pr.find(1,2)
   while idx!=-1:
    if idx!=23 and has(idx):principal+=1
    idx=pr.find(1,idx+1)
   assert principal==records
   small=values(100000)
   assert all(has(n)==(n in small)for n in range(100001))
  return counts

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path);ap.add_argument('--output',type=Path)
 a=ap.parse_args();root=root_path(a.repo_root);ledger=json.loads((HERE/LEDGER).read_text())
 essential={('packet',n)for n in [PROOF,Path(__file__).name,INPUT,WITNESS,NATIVE]}
 essential|={('repository','research/arithmetic/'+n)for n in ['PRIMITIVE_NORM_CLASS_NUMBER_THREE.md','PRIMITIVE_NORM_MINUS23_EFFECTIVE.md','PRIMITIVE_NORM_MINUS23_SPLIT_TARGETS.md','MINUS23_NEUTRAL_CUBE_ATOM_CLOSURE.md','MINUS23_PRINCIPAL_ATOMS_AND_INERT_CORE_FIVE.md','sources/rosser_schoenfeld_1962.pdf']}
 keys={(v['location'],v['path'])for v in ledger['files']};assert essential<=keys and len(keys)==len(ledger['files'])
 for v in ledger['files']:
  assert v['location']in {'packet','repository'}
  p=(HERE if v['location']=='packet'else root)/v['path'];assert digest(p)==v['sha256'],str(p)
 inp=json.loads((HERE/INPUT).read_text());c=controls(inp);nc=cube_check();n=native_check(inp)
 out={'status':'PASS','scope':'Complete core7 phase/atom proof controls and ACTUAL exhaustive principal/cube finite bridges. No other inert-core or global40 classification. Internal replay only.',
      'controls':c,'cube_witnesses_verified':nc,'principal_counts':n,
      'proof_sha256':digest(HERE/PROOF),'checker_sha256':digest(Path(__file__)),
      'input_sha256':digest(HERE/INPUT),'cube_witness_sha256':digest(HERE/WITNESS),
      'native_sha256':digest(HERE/NATIVE),'ledger_sha256':digest(HERE/LEDGER),
      'versions':{'python':sys.version.split()[0],'g++':subprocess.check_output(['g++','--version'],text=True).splitlines()[0]}}
 dest=a.output or HERE/'core_seven_completion_verification.json'
 dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
 print('PASS all980245 principal targets,57 cubes,27 phase cases; exact onsets2^31/2^34');print(dest)
if __name__=='__main__':main()
