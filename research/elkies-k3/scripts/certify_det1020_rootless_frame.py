"""Exact finite witness replay; no genus enumeration or equation assertion."""
from sage.all import ZZ, matrix, lcm, gcd
from pathlib import Path
import argparse, hashlib, json
R = Path(__file__).resolve().parents[2]
P = R/'artifacts/generated-results/elkies-k3-det1020-root-gate-v1'
def rows(M): return [[int(z) for z in r] for r in M.rows()]
def generator(H):
    I = H.inverse(); n = abs(H.det())
    for a in range(H.nrows()):
        for b in range(a, H.nrows()):
            v = I.row(a)+I.row(b)
            if lcm(z.denominator() for z in v) == n:
                return v, v*H*v
    raise AssertionError('No generator in pinned search')
def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--check', action='store_true'); args = parser.parse_args()
    source = P/'4e6-probe.json'; x = json.loads(source.read_text()); h = x['hits'][0]
    G = matrix(ZZ,x['ambient_gram']); A = matrix(ZZ,h['embedding']); B = matrix(ZZ,h['frame_basis']); F = matrix(ZZ,h['frame_gram'])
    assert G.nrows()==24 and G.is_symmetric() and G.is_positive_definite() and G.det()==1
    assert all(z%2==0 for z in G.diagonal())
    assert A.nrows()==7 and A.ncols()==24 and A.rank()==7
    assert A.elementary_divisors()==[1]*7
    K=A*G*A.transpose(); C=K[:6,:6]
    assert C.det()==3 and C.is_positive_definite()
    assert all(C[i,i]==2 for i in range(6)) and all(C[i,j] in [0,-1] for i in range(6) for j in range(i))
    assert sorted(sum(C[i,j]==-1 for j in range(6)) for i in range(6))==[1,1,1,2,2,3]
    assert list(K[6])==[0]*6+[340]
    assert B.nrows()==17 and B.rank()==17 and B.elementary_divisors()==[1]*17
    assert A*G*B.transpose()==0 and B*G*B.transpose()==F
    assert F.det()==1020 and F.is_positive_definite() and all(z%2==0 for z in F.diagonal())
    roots=F.__pari__().qfminim(2,100000,2); assert roots[0]==0
    shell=F.__pari__().qfminim(4,100000,2); assert shell[0]>0
    T=matrix(ZZ,[[-2,1,0],[1,2,0],[0,0,204]])
    assert F.elementary_divisors()==[1]*16+[1020]
    assert T.elementary_divisors()==[1,1,1020]
    v,q=generator(F); w,r=generator(T); u=ZZ(59)
    assert gcd(u,1020)==1 and (q-u*u*r)/2 in ZZ
    result={'status':'PASS_ROOTLESS_FRAME_AND_DISCRIMINANT_MATCH','input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
      'rank':17,'determinant':1020,'minimum':4,'norm4_vectors':int(shell[0]),'norm2_vectors':0,
      'auxiliary_gram':rows(K),'frame_gram':rows(F),'frame_basis':rows(B),
      'cyclic_generator_frame':[str(z) for z in v],'cyclic_generator_T':[str(z) for z in w],
      'q_frame':str(q),'q_T':str(r),'discriminant_isometry_multiplier':59,
      'written_inputs':['Indefinite even lattice uniqueness for rank at least discriminant length plus two identifies U+(-F) with S=(-T)+E8(-1)^2.',
      'Admitted determinant1020 source has geometric NS exactly S and every divisor class descends to Q.',
      'Nef chamber reduction and elliptic K3 theory yield a Q fibration with Q section; Shioda-Tate and primitivity give torsion-free MW rank17 and height lattice F.'],
      'limitations':['No explicit isometry to source NS basis, nef fiber class, K3 equation, period coordinates or section coordinates.','No independent implementation, formal verification or novelty assertion.']}
    target=P/'certificate.json'
    if args.check: assert json.loads(target.read_text())==result
    else: target.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
if __name__=='__main__': main()
