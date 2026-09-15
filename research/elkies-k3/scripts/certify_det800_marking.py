"""Exact determinant-800 full discriminant-kernel and local-level replay."""
import argparse, hashlib, json
from pathlib import Path
from sage.all import *
from collections import deque
N=50
G=matrix(QQ,[[0,0,2],[0,200,0],[2,0,0]])
J=matrix(QQ,[[0,-1],[1,0]])
D=diagonal_matrix(QQ,[2,200,2])
def action(g, local=False):
 a,b,c,d=(J.inverse()*g*J).list();delta=g.det()
 A=matrix(QQ,[[a*a,-2*a*b,-b*b/N],[-a*c,a*d+b*c,b*d/N],[-N*c*c,2*N*c*d,d*d]])/delta
 assert A.transpose()*G*A==G
 assert all(x in ZZ or (local and (x==0 or x.valuation(2)>=0)) for x in A.list()),(g,A)
 return A
moduli=[2,200,2]
def key(A):
 B=D*A*D.inverse()
 assert all(x in ZZ for x in B.list())
 return tuple(int(B[i,j]%moduli[i]) for i in range(3) for j in range(3))
def mul(a,b):
 A=matrix(ZZ,3,3,a)*matrix(ZZ,3,3,b)
 return tuple(int(A[i,j]%moduli[i]) for i in range(3) for j in range(3))
id=key(identity_matrix(QQ,3))
gens=[key(action(matrix(QQ,g.matrix()))) for g in Gamma0(50).gens()]
image={id};queue=deque([id])
while queue:
 x=queue.popleft()
 for g in gens:
  y=mul(x,g)
  if y not in image:image.add(y);queue.append(y)
assert len(image)==4
coset_records=[]
for label,w in [('1',identity_matrix(QQ,2)),('2',matrix(QQ,[[2,1],[50,26]])),('25',matrix(QQ,[[25,1],[600,25]])),('50',matrix(QQ,[[0,-1],[50,0]]))]:
 A=action(w)
 for sign in [1,-1]:
  stable=id in {mul(key(sign*A),g) for g in image}
  assert stable==(label=='1' and sign==1)
  coset_records.append({'label':label,'sign':sign,'contains_stable':stable,'action':[[int(x) for x in row] for row in A.rows()]})
assert Gamma0(200).genus()==19
Gi=G.inverse();I=identity_matrix(QQ,3)
W=matrix(QQ,[[2,1],[50,26]])
counts={};units=set()
for sign in [1,-1]:
 for coset in [0,1]:
  count=0
  for a in [1,3,5,7]:
   for d in [1,3,5,7]:
    for b in range(8):
     for c in range(8):
      gamma=matrix(QQ,[[a,b],[50*c,d]])
      # Odd determinant means a local unit, not a global integral unit.
      A=action(gamma if not coset else W*gamma, local=True)
      stable=all(x==0 or x.valuation(2)>=0 for x in ((sign*A-I)*Gi).list())
      expected=sign==1 and coset==0 and b%2==0 and c%2==0
      assert stable==expected,(a,b,c,d,sign,coset)
      count+=stable
      if stable:units.add(int(gamma.det()%8))
  counts[str((sign,coset))]=count
assert counts=={'(1, 0)':256,'(1, 1)':0,'(-1, 0)':0,'(-1, 1)':0}
assert units=={1,3,5,7}

C=diagonal_matrix(QQ,[QQ(1)/2,1])
kernel_generators=[]
for generator in Gamma0(200).gens():
    gamma=C.inverse()*matrix(QQ,generator.matrix())*C
    assert all(x in ZZ for x in gamma.list()) and gamma[1,0]%50==0
    assert key(action(gamma))==id
    kernel_generators.append([[int(x) for x in row] for row in gamma.rows()])
assert Gamma0(200).index()==4*Gamma0(50).index()
root=Path(__file__).resolve().parents[2]
input_path=root/'artifacts/generated-results/elkies-k3-rank7-t-arithmetic-v1.json'
row=next(r for r in json.loads(input_path.read_text())['surfaces'] if r['surface_id']=='K3-643eed54b4d256fc')
assert matrix(QQ,row['literal_transcendental_gram'])==G
assert G.det()==-800
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--check',action='store_true')
args=parser.parse_args()
result={'schema':'elkies-k3.det800-full-marking.v1',
        'status':'PASS_EXACT_REPLAY_WITH_WRITTEN_NORMALIZER_AND_PERIOD_INPUTS',
        'surface_id':row['surface_id'],
        'gram':row['literal_transcendental_gram'],
        'input':{'path':str(input_path.relative_to(root)),'sha256':hashlib.sha256(input_path.read_bytes()).hexdigest()},
        'unit_discriminant_image_order':len(image),
        'unit_discriminant_image':[list(g) for g in sorted(image)],
        'signed_normalizer_cosets':coset_records,
        'kernel_generators':kernel_generators,
        'kernel_index_in_Gamma0_50':4,
        'full_marked_curve':'X_0(200)',
        'canonical_model_field':'Q',
        'genus':19,
        'local2_counts':counts,
        'local2_stable_determinant_residues':sorted(units),
        'rational_non_CM_marking_exists':False,
        'theorem_inputs':['Written exhaustive Clifford-order normalizer argument','Written local-level and compatible orientation-character descent','Canonical arithmetic K3 period map','Rational cyclic-isogeny classification excluding degree50'],
        'new_positive_handoff':[]}
output=root/'artifacts/generated-results/elkies-k3-det800-full-marking-v1.json'
if args.check:
    assert json.loads(output.read_text())==result
else:
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print('PASS determinant800: full stable X_0(200), genus19, arithmetic exclusion')
