"""Boundary arithmetic and fail-closed framing for the retained-list score engine."""
import shutil
import struct
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]

@unittest.skipUnless(shutil.which('g++'),'g++ required for the compiled score engine')
class RetainedCacheScorer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory();cls.folder=Path(cls.temp.name);cls.binary=cls.folder/'scorer'
        subprocess.run(['g++','-O3','-std=c++17',str(ROOT/'elliptic-curves/cas/newfamily/score_retained_projective_cache.cpp'),'-o',str(cls.binary)],check=True,capture_output=True,timeout=20)
    @classmethod
    def tearDownClass(cls):cls.temp.cleanup()
    def fixture(self):
        primes=[4099,4111];tables=[];data=b'R17XS001'+struct.pack('<I',len(primes))
        for p in primes:
            good=[int(i%11!=0) for i in range(p+1)]
            units=[(i-p//2)*1000000 if good[i] else 0 for i in range(p+1)]
            tables.append((p,units,good));data+=struct.pack('<II',p,p+1)+struct.pack('<'+'q'*(p+1),*units)+bytes(good)
        data+=b'ENDXSC01';cache=self.folder/'cache.bin';cache.write_bytes(data)
        pairs=[(1,1),(-1,1),(4099,1),(1,4099),(-4098,4099),(1,4111),(131071,131072),(-131071,131072)]
        path=self.folder/'candidates.txt';path.write_text('R17-CANDIDATES-V1 '+str(len(pairs))+'\n'+''.join(f'{n} {d}\n' for n,d in pairs))
        return cache,path,tables,pairs
    def test_signed_residues_and_projective_infinity(self):
        cache,path,tables,pairs=self.fixture();r=subprocess.run([str(self.binary),str(cache),str(path)],text=True,capture_output=True,timeout=5)
        self.assertEqual(r.returncode,0,r.stderr);expected=[]
        for i,(n,d) in enumerate(pairs):
            score=good=0
            for p,units,flags in tables:
                t=n*pow(d,-1,p)%p if d%p else p
                score+=units[t];good+=flags[t]
            expected.append(f'R {i} {score} {good}')
        self.assertEqual(r.stdout.splitlines(),expected+[f'S {len(pairs)} {len(tables)}'])
    def test_truncation_yields_no_partial_score_output(self):
        cache,path,_,_=self.fixture();cache.write_bytes(cache.read_bytes()[:-1]);r=subprocess.run([str(self.binary),str(cache),str(path)],text=True,capture_output=True,timeout=5)
        self.assertNotEqual(r.returncode,0);self.assertEqual(r.stdout,'')

if __name__=='__main__':unittest.main()
