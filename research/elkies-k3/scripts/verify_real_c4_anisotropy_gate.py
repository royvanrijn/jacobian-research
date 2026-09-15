#!/usr/bin/env python3
"""Finite premises for the real-locus obstruction to the order-four mechanism."""
import argparse,hashlib,json,re
from itertools import product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
PACKET=ROOT/'artifacts/generated-results/elkies-k3-real-c4-anisotropy-gate-v1'
CASES=[('948',[[-2,0,1],[0,4,0],[1,0,118]],'elkies-k3/SURFACE_FIBRATION_CHARACTER_GRAPH_2026-09-04.md'),('1092',[[-2,1,0],[1,2,2],[0,2,220]],'elliptic-curves/notes/DET1092_MARKING_AND_CM_GATE_2026-09-07.md'),('1020',[[-2,1,0],[1,2,0],[0,0,204]],'elkies-k3/DET1020_RATIONAL_MARKING_SOURCE_2026-09-15.md')]
def det(m):
 a,b,c=m[0];d,e,f=m[1];g,h,i=m[2]
 return a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)

def run():
 rows=[]
 for name,T,source in CASES:
  text=(ROOT/source).read_text();literal=json.dumps(T,separators=(',',':'))
  assert literal in re.sub(r'\s','',text)
  assert det(T)==-int(name)
  primitive=[v for v in product(range(9),repeat=3) if any(x%3 for x in v)]
  zeros=[v for v in primitive if sum(T[i][j]*v[i]*v[j] for i in range(3) for j in range(3))%9==0]
  assert len(primitive)==702 and not zeros
  rows.append({'determinant':int(name),'T':T,'source':source,'source_sha256':hashlib.sha256((ROOT/source).read_bytes()).hexdigest(),'modulus':9,'primitive_at_3_vectors_checked':702,'isotropic_residues':0})
 # An even rank-two Gram has determinant 0 or3 modulo4, excluding -2.
 assert {(4*a*c-b*b)%4 for a,b,c in product(range(4),repeat=3)}=={0,3}
 assert [n for n in [1,2,4] if (-n)%4 in {0,3}]==[1,4]
 assert 2+(2*1-22)==-18 and (-18//2)%2==1
 # Isotropic boundary: the Inose-type U + <622> has a rational null vector.
 control=[[0,1,0],[1,0,0],[0,0,622]];assert control[0][0]==0 and det(control)==-622
 return {'status':'PASS','checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'cases':rows,'rank2_invariant_lattice_possible_absolute_determinants':[1,4],'anisotropic_full_marking_real_euler_characteristic':-18,'hypothetical_free_orientable_quotient_euler_characteristic':-9,'isotropic_control_T':control,'geometric_argument_formally_verified':False,'positive_correlated_target_complete':False}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();out=run()
 if args.write:
  PACKET.mkdir(exist_ok=True);(PACKET/'result.json').write_text(json.dumps(out,indent=2)+'\n')
 else:assert out==json.loads((PACKET/'result.json').read_text())
 print('PASS: three literal T lattices anisotropic modulo9; determinant and Euler-parity premises checked')
