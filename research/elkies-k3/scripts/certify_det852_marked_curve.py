#!/usr/bin/env sage-python
"""Exact determinant852 full stable marked curve and rational elliptic locus."""
import argparse,hashlib,json,re
from pathlib import Path
from sage.all import *
ROOT=Path(__file__).resolve().parents[2]
PACKET=ROOT/'artifacts/generated-results/elkies-k3-det852-marking-gate-v1'
CAT=ROOT/'artifacts/generated-results/elkies-k3-rank7-auxiliary-catalogue-v1.json'
TINPUT=ROOT/'artifacts/generated-results/elkies-k3-rank7-t-arithmetic-v1.json'
OUTPUT=PACKET/'certificate.json'
sid='K3-4ff75fec54d01662'
row=next(r for r in json.loads(CAT.read_text())['surfaces'] if r['surface_id']==sid)
tr=next(r for r in json.loads(TINPUT.read_text())['surfaces'] if r['surface_id']==sid)
G=matrix(ZZ,tr['literal_transcendental_gram'])
assert G==matrix(ZZ,[[-2,0,1],[0,4,0],[1,0,106]]) and G.det()==-852
assert QuadraticForm(QQ,G).signature_vector()==(2,1,0)
C=CliffordAlgebra(QuadraticForm(QQ,G));v=C.gens()
B=[C.one(),v[0]*v[1],v[0]*v[2],v[1]*v[2]];keys=list(C.basis().keys());inds=[0,4,5,6]
def coords(x):return vector(QQ,[x.monomial_coefficients().get(keys[i],0) for i in inds])
def elem(c):return sum((a*b for a,b in zip(c,B)),C.zero())
def trace(x):return coords(x+x.clifford_conjugate())[0]
def strings(A):return [[str(x) for x in r] for r in A.rows()]
P=matrix(QQ,4,4,lambda i,j:trace(B[i]*B[j]))
assert P.det()==-426**2
assert P==matrix(QQ,tr['clifford']['integral_even_clifford_order']['reduced_trace_pairing'])
assert QuaternionAlgebra(QQ,2,QQ(213)/4).discriminant()==6
for b in B:
 for c in B:assert all(x in ZZ for x in coords(b*c))
vol=v[0]*v[1]*v[2]-G[0,1]/2*v[2]+G[0,2]/2*v[1]-G[1,2]/2*v[0]
images=[vol*x for x in v];M=matrix(QQ,[coords(x) for x in images]).transpose();left=(M.transpose()*M).inverse()*M.transpose()
D,U,V=G.smith_form();assert D.diagonal()==[1,1,852]
dual=G.inverse()*U.inverse()*vector(ZZ,[0,0,1]);q=dual*G*dual/2
ounits=[n for n in range(852) if gcd(n,852)==1 and (n*n-1)*q in ZZ]
assert ounits==[1,143,283,425,427,569,709,851]


records=[]
for label,cc,expected in [(2,[28,-11,-10,-7],427),(6,[20,-11,-10,7],143),(426,[0,1,0,-2],851)]:
 x=elem(cc);norm=coords(x*x.clifford_conjugate())[0];assert norm==label
 xi=x.clifford_conjugate()/norm
 R=matrix(QQ,[coords(x*b*xi) for b in B]).transpose();A=matrix(QQ,[left*coords(x*b*xi) for b in images]).transpose()
 assert all(z in ZZ for z in R.list()+A.list()) and abs(R.det())==1 and A.det()==1
 assert A.transpose()*G*A==G
 act=U*A.inverse().transpose()*U.inverse();mu=int(act[2,2]%852);assert mu==expected
 records.append({'label':label,'coordinates':cc,'norm':int(norm),'order_action':strings(R),'T_action':strings(A),'discriminant_multiplier':mu})
I=identity_matrix(QQ,3);Gi=G.inverse()
assert -A==I+matrix(QQ,[[1],[0],[0]])*matrix(QQ,[G.row(0)])
def diffnum(c):
 x=elem(c);bar=x.clifford_conjugate();norm=coords(x*bar)[0]
 return (matrix(QQ,[left*coords(x*b*bar) for b in images]).transpose()-norm*I)*Gi
basis=list(identity_matrix(QQ,4).rows());polys=[diffnum(c) for c in basis]
polys += [diffnum(basis[i]+basis[j])-diffnum(basis[i])-diffnum(basis[j]) for i in range(4) for j in range(i)]
assert all(z in ZZ for a in polys for z in a.list())
actions={1:1,2:427,6:143,426:851}
actions[3]=actions[2]*actions[6]%852;actions[71]=actions[6]*actions[426]%852;actions[142]=actions[3]*actions[426]%852;actions[213]=actions[2]*actions[426]%852
assert set(actions.values())==set(ounits)
assert sorted(a for a,u in actions.items() if u in [1,851])==[1,426]
import runpy
helper=ROOT/'elkies-k3/scripts/certify_det1236_marked_shimura_curve.sage'
h=runpy.run_path(str(helper));fixed=h['fixed_point_record'](ZZ(6),ZZ(71),ZZ(426));assert fixed['fixed_points']==24
mu=2*(71+1);e2=prod(1-kronecker(-4,p) for p in [2,3])*(1+kronecker(-4,71));e3=prod(1-kronecker(-3,p) for p in [2,3])*(1+kronecker(-3,71));top=1+QQ(mu)/12-QQ(e2)/4-QQ(e3)/3
assert top==13 and (2*top+2-fixed['fixed_points'])/4==1
source=ROOT/'artifacts/generated-results/elkies-k3-shimura-positive-source-preflight-v1'
rankfile=source/'genus_1_AL_quotients_rat_pts_pos_rank.m';modelfile=source/'genus_1_AL_quotient_jacobian_isomorphism_classes.m'
assert re.search(r'\[\*\s*6,\s*71,\s*\{\s*426\s*\}',rankfile.read_text())
assert re.search(r'\[\*\s*6,\s*71,\s*\{\s*426\s*\},\s*"426b1"',modelfile.read_text())
EC=EllipticCurve(QQ,[1,1,0,-286,1780]);point=EC(7,-17)
counts={p:int(EC.change_ring(GF(p)).cardinality()) for p in [5,7,11]};torsion_bound=gcd(list(counts.values()))
assert EC.conductor()==426 and torsion_bound*point!=EC(0)
out={'schema':'elkies-k3.det852-full-marked-curve.v1','status':'PASS_FULL_CURVE_WITH_INFINITE_RATIONAL_LOCUS_EXPLICIT_NONCM_GATE_OPEN','surface_id':sid,'T_gram':strings(G),'Clifford_trace_pairing':strings(P),'quaternion_discriminant':6,'Eichler_level':71,'normalizers':records,'discriminant_orthogonal_units':ounits,'Atkin_Lehner_actions':actions,'local_unit_stability_coefficients':[strings(a) for a in polys],'genus':{'upstairs':int(top),'full_Fr ic ke_fixed_points'.replace(' ',''):fixed,'marked_quotient':1},'marked_curve':'X_0^6(71)/<w_426>','rational_curve_model':{'label':'426b1','ainvs':list(EC.ainvs()),'non_torsion_point':[7,-17,1],'good_reduction_orders':counts,'torsion_bound':int(torsion_bound),'point_proved_non_CM':None,'origin_as_moduli_point':None},'inputs':[{'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in [TINPUT,helper,rankfile,modelfile]],'theorem_inputs':['Squarefree reduced-discriminant Eichler order classification and normalizer theorem','Canonical arithmetic full stable K3 period curve identification','Ogg fixed-point formula','Padurariu-Saia quotient rational-point existence and exact Jacobian isomorphism table'],'boundary':'Full marked curve and infinitely many rational points proved. A particular non-CM period, primitive NS/divisor descent, rootless fibration and MW17 equation are not certified here; no new frame work is authorized by this packet.'}
a=argparse.ArgumentParser();a.add_argument('--check',action='store_true');args=a.parse_args()
# JSON roundtrip normalizes integer dictionary keys for comparison.
out=json.loads(json.dumps(out,default=int))
if args.check:assert json.loads(OUTPUT.read_text())==out
else:OUTPUT.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print('PASS det852 full stable curve =426b1; rational point of infinite order; explicit non-CM marking gate remains open.')
