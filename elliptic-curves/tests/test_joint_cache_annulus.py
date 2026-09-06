"""Exhaustive signed/annular/tie-order comparisons against direct modular sums."""
import math,shutil,struct,subprocess,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
@unittest.skipUnless(shutil.which('g++'),'g++ required')
class JointCacheAnnulus(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.temp=tempfile.TemporaryDirectory();cls.d=Path(cls.temp.name);cls.binary=cls.d/'scanner'
  subprocess.run(['g++','-O3','-std=c++17',str(ROOT/'elliptic-curves/cas/newfamily/scan_joint_cache_annulus.cpp'),'-o',str(cls.binary)],check=True,capture_output=True,timeout=20)
 @classmethod
 def tearDownClass(cls):cls.temp.cleanup()
 def caches(self,ties=False):
  tables=[];paths=[]
  for j,primes in enumerate(([5,7],[11,13])):
   data=b'R17XS001'+struct.pack('<I',2)
   for p in primes:
    good=[1 if ties else int(i%3!=0) for i in range(p+1)];units=[0 if ties or not good[i] else (i%5-2)*10**6 for i in range(p+1)];tables.append((p,units,good));data+=struct.pack('<II',p,p+1)+struct.pack('<'+'q'*(p+1),*units)+bytes(good)
   data+=b'ENDXSC01';path=self.d/f'cache{j}.bin';path.write_bytes(data);paths.append(path)
  return paths,tables
 def invoke(self,paths,sign,N,M,K,shard,shards,inner):
  return subprocess.run([str(self.binary),*map(str,paths),*map(str,[sign,N,M,K,shard,shards,inner,4])],text=True,capture_output=True,timeout=5)
 def check_case(self,ties):
  paths,tables=self.caches(ties)
  for sign in (-1,1):
   for N,M,inner,shard,shards in [(31,29,0,0,1),(31,29,7,2,3),(5,29,13,1,7),(31,5,7,4,5)]:
    rows=[]
    for d in range(shard+1,M+1,shards):
     for n in range(1,N+1):
      if max(n,d)<=inner or math.gcd(n,d)!=1:continue
      score=good=0
      for p,u,g in tables:
       t=sign*n*pow(d,-1,p)%p if d%p else p;score+=u[t];good+=g[t]
      rows.append((sign*n,d,score,good))
    rows.sort(key=lambda r:(-r[2],-r[3],r[1],abs(r[0])))
    for K in (1,7,1000):
     r=self.invoke(paths,sign,N,M,K,shard,shards,inner);self.assertEqual(r.returncode,0,r.stderr)
     expected=['JOINT_NAGAO_ANNULUS_V1',f'P {sign} {N} {M} {K} {shard} {shards} {inner} 4']+['C '+' '.join(map(str,row)) for row in rows[:K]]+[f'S {len(rows)} {min(K,len(rows))}']
     self.assertEqual(r.stdout.splitlines(),expected)
 def test_exhaustive_signed_frames(self):self.check_case(False)
 def test_exact_ties(self):self.check_case(True)
 def test_truncated_and_overlapping_cache_fail_before_output(self):
  paths,_=self.caches();paths[1].write_bytes(paths[1].read_bytes()[:-1]);r=self.invoke(paths,1,31,29,7,0,1,7);self.assertNotEqual(r.returncode,0);self.assertEqual(r.stdout,'')
  paths,_=self.caches();r=self.invoke(paths[::-1],1,31,29,7,0,1,7);self.assertNotEqual(r.returncode,0);self.assertEqual(r.stdout,'')
if __name__=='__main__':unittest.main()
