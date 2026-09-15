#!/usr/bin/env sage-python
"""Exact det1020 full-period point and primitive marking data; descent is a written proof."""
import argparse,hashlib,json,re
from pathlib import Path
from sage.all import *
ROOT=Path(__file__).resolve().parents[2]
PACKET=ROOT/'artifacts/generated-results/elkies-k3-det1020-marking-gate-v1'
CAT=ROOT/'artifacts/generated-results/elkies-k3-rank7-auxiliary-catalogue-v1.json'
TINPUT=ROOT/'artifacts/generated-results/elkies-k3-rank7-t-arithmetic-v1.json'
OUTPUT=PACKET/'certificate.json'
G=matrix(ZZ,[[-2,1,0],[1,2,0],[0,0,204]])
assert G.det()==-1020 and QuadraticForm(QQ,G).signature_vector()==(2,1,0)
C=CliffordAlgebra(QuadraticForm(QQ,G));v=C.gens()
B=[C.one(),v[0]*v[1],v[0]*v[2],v[1]*v[2]];keys=list(C.basis().keys());inds=[0,4,5,6]
def coords(x):return vector(QQ,[x.monomial_coefficients().get(keys[i],0) for i in inds])
def elem(c):return sum((a*b for a,b in zip(c,B)),C.zero())
def trace(x):return coords(x+x.clifford_conjugate())[0]
def strings(A):return [[str(x) for x in r] for r in A.rows()]
P=matrix(QQ,4,4,lambda i,j:trace(B[i]*B[j]))
assert P.det()==-510**2
assert P==matrix(QQ,[[2,1,0,0],[1,3,0,0],[0,0,204,-102],[0,0,-102,-204]])
assert QuaternionAlgebra(QQ,5,102).discriminant()==510
for b in B:
 for c in B:assert all(x in ZZ for x in coords(b*c))
vol=v[0]*v[1]*v[2]-G[0,1]/2*v[2]+G[0,2]/2*v[1]-G[1,2]/2*v[0]
images=[vol*x for x in v];M=matrix(QQ,[coords(x) for x in images]).transpose();left=(M.transpose()*M).inverse()*M.transpose()
D,U,V=G.smith_form();assert D.diagonal()==[1,1,1020]
dual=G.inverse()*U.inverse()*vector(ZZ,[0,0,1]);q=dual*G*dual/2
ounits=[n for n in range(1020) if gcd(n,1020)==1 and (n*n-1)*q in ZZ]
assert len(ounits)==16



records=[];actions={1:1}
for label,cc,expected in [(2,[46,-16,-12,-7],511),(3,[174,-15,-15,4],341),(5,[54,-13,-14,-8],409),(510,[0,0,1,2],1019)]:
 x=elem(cc);norm=coords(x*x.clifford_conjugate())[0];assert norm==label;xi=x.clifford_conjugate()/norm
 R=matrix(QQ,[coords(x*b*xi) for b in B]).transpose();A=matrix(QQ,[left*coords(x*b*xi) for b in images]).transpose()
 assert all(z in ZZ for z in R.list()+A.list()) and abs(R.det())==1 and A.det()==1 and A.transpose()*G*A==G
 act=U*A.inverse().transpose()*U.inverse();mu=int(act[2,2]%1020);assert mu==expected
 records.append({'label':label,'coordinates':cc,'norm':int(norm),'order_action':strings(R),'T_action':strings(A),'discriminant_multiplier':mu})
 for d,a in list(actions.items()):actions[d*label//gcd(d,label)**2]=a*mu%1020
assert len(actions)==16 and set(actions.values())==set(ounits)
assert sorted(d for d,a in actions.items() if a in [1,1019])==[1,510]
I=identity_matrix(QQ,3);Gi=G.inverse();assert -A==I+matrix(QQ,[[1],[0],[0]])*matrix(QQ,[G.row(0)])
def diffnum(c):
 x=elem(c);bar=x.clifford_conjugate();norm=coords(x*bar)[0]
 return (matrix(QQ,[left*coords(x*b*bar) for b in images]).transpose()-norm*I)*Gi
basis=list(identity_matrix(QQ,4).rows());polys=[diffnum(c) for c in basis];polys += [diffnum(basis[i]+basis[j])-diffnum(basis[i])-diffnum(basis[j]) for i in range(4) for j in range(i)]
assert all(z in ZZ for a in polys for z in a.list())
import runpy
helper=ROOT/'elkies-k3/scripts/certify_det1236_marked_shimura_curve.sage';h=runpy.run_path(str(helper));fixed={int(m):h['fixed_point_record'](ZZ(510),ZZ(1),m) for m in divisors(ZZ(510)) if m>1}
assert fixed[510]['fixed_points']==16 and fixed[3]['fixed_points']==8
mu=prod(p-1 for p in [2,3,5,17]);e2=prod(1-kronecker(-4,p) for p in [2,3,5,17]);e3=prod(1-kronecker(-3,p) for p in [2,3,5,17]);top=1+QQ(mu)/12-QQ(e2)/4-QQ(e3)/3;assert top==9 and (2*top+2-16)/4==1
source=ROOT/'artifacts/generated-results/elkies-k3-shimura-positive-source-preflight-v1';rankfile=source/'genus_1_AL_quotients_rat_pts_pos_rank.m';modelfile=source/'genus_1_AL_quotient_jacobian_isomorphism_classes.m'
assert re.search(r'\[\*\s*510,\s*1,\s*\{\s*510\s*\}',rankfile.read_text())
assert re.search(r'\[\*\s*510,\s*1,\s*\{\s*510\s*\},\s*"510d2"',modelfile.read_text())
EC=EllipticCurve(QQ,[1,1,1,-421,-3157]);assert EC.conductor()==510
rank=EC.pari_curve().ellrank();assert int(rank[0])==int(rank[1])==1
Tors=EC.torsion_subgroup();assert Tors.invariants()==(2,2) and gcd([EC.change_ring(GF(p)).cardinality() for p in [7,11,13]])==4

CLASS_NUMBER_ONE_ORDERS = [
    (1, -3), (2, -3), (3, -3), (1, -4), (2, -4),
    (1, -7), (2, -7), (1, -8), (1, -11), (1, -19),
    (1, -43), (1, -67), (1, -163),
]
CLASS_NUMBER_TWO_ORDERS = [
    (4, -3), (5, -3), (7, -3), (3, -4), (4, -4),
    (5, -4), (4, -7), (2, -8), (3, -8), (3, -11),
    (1, -15), (2, -15), (1, -20), (1, -24), (1, -35),
    (1, -40), (1, -51), (1, -52), (1, -88), (1, -91),
    (1, -115), (1, -123), (1, -148), (1, -187),
    (1, -232), (1, -235), (1, -267), (1, -403),
    (1, -427),
]

from sage.quadratic_forms.binary_qf import BinaryQF_reduced_representatives
cm=[]
for hh,orders in [(1,CLASS_NUMBER_ONE_ORDERS),(2,CLASS_NUMBER_TWO_ORDERS)]:
 for f,d in orders:
  disc=f*f*d;assert len(BinaryQF_reduced_representatives(ZZ(disc),primitive_only=True))==hh
  factors=[int(1-(1 if f%p==0 else kronecker(d,p))) for p in [2,3,5,17]];count=hh*prod(factors);ar=prod(p for p in [2,3,5,17] if f%p and kronecker(d,p)==-1);mr=gcd(510,abs(disc));quot=510//mr
  rat=bool(count and ((ar==1 or quot==ar) if hh==1 else (ar==1 and quot==1)))
  if hh==2:assert not(count and ar==1 and quot==1)
  cm.append({'order_discriminant':disc,'class_number':hh,'conductor':f,'local_factors':factors,'D_R':int(ar),'m_R':int(mr),'m_over_m_R':int(quot),'upstairs_count':int(count),'rational':rat,'rational_image_count':int(count//2) if rat else 0})
assert [(r['order_discriminant'],r['rational_image_count']) for r in cm if r['rational']]==[(-3,4),(-163,8)]
assert sum(c['contribution'] for c in fixed[3]['cm_contributions'] if c['order_discriminant']==-3)==8
assert all(c['order_discriminant']!=-163 or c['contribution']==0 for f in fixed.values() for c in f['cm_contributions'])
# Primitive S/T embedding in U^3 + E8(-1)^2, with explicit complementary bases.
H=block_matrix(ZZ,[[zero_matrix(ZZ,3),identity_matrix(ZZ,3)],[identity_matrix(ZZ,3),zero_matrix(ZZ,3)]])
b=matrix(ZZ,3,3,lambda i,j:G[i,i]//2 if i==j else (G[i,j] if i<j else 0));TT=identity_matrix(ZZ,3).augment(b);SS=identity_matrix(ZZ,3).augment(-b.transpose())
assert TT*H*TT.transpose()==G and SS*H*SS.transpose()==-G and SS*H*TT.transpose()==0
C8=matrix(ZZ,CartanMatrix(['E',8]));S=block_diagonal_matrix(-G,-C8,-C8)
assert S.det()==1020 and QuadraticForm(QQ,S).signature_vector()==(1,18,0)
# Identity leading minors prove both rank3 embeddings primitive.
assert TT[:,:3]==SS[:,:3]==identity_matrix(ZZ,3)
# Adjacent E8 roots descend and intersect in a degree-one zero-cycle.
i,j=next((i,j) for i in range(8) for j in range(i) if C8[i,j]==-1)
r=vector(ZZ,[int(k==3+i) for k in range(19)]);s=vector(ZZ,[int(k==3+j) for k in range(19)]);assert r*S*r==s*S*s==-2 and r*S*s==1
# An explicit U exists, for later geometric work after arithmetic admission.
v=vector(ZZ,[1,0,0]+[1]+[0]*15);w=v-vector(ZZ,[0,1,0]+[0]*16);assert v*S*v==w*S*w==0 and v*S*w==1
out={'schema':'elkies-k3.det1020-rational-marking.v1','status':'PASS_FULL_RATIONAL_SATURATED_NS_BY_INTRINSIC_NONCM_AND_DESCENT','T_gram':strings(G),'NS_gram':strings(S),'Clifford_trace_pairing':strings(P),'maximal_order_reduced_discriminant':510,'normalizers':records,'Atkin_Lehner_actions':{int(k):int(v) for k,v in actions.items()},'local_unit_stability_coefficients':[strings(x) for x in polys],'full_marked_curve':'X^510/<w510> overQ, isomorphic to510d2','genus':{'upstairs':int(top),'Fricke_fixed_points':16,'quotient':1},'elliptic_model':{'ainvs':list(EC.ainvs()),'rank_bounds':[int(rank[0]),int(rank[1])],'rational_torsion_invariants':[2,2]},'CM_orders':cm,'fixed_points':fixed,'intrinsic_nonCM_point':'Choose any rational CM(-3) origin O and any rational CM(-163) point R. Then O+2(R-O) is rational non-CM. Its coordinates are uncomputed.','primitive_embedding':{'U3_gram':strings(H),'T_basis':strings(TT),'negative_T_complement_basis':strings(SS)},'divisor_descent_pair':{'r':list(r),'s':list(s),'intersection_gram':[[-2,1],[1,-2]]},'embedded_U':[list(v),list(w)],'inputs':[{'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in [helper,rankfile,modelfile]],'theorem_inputs':['Maximal quaternion-order normalizer and local-unit reduced-norm surjectivity','Canonical arithmetic full stable K3 period identification','Padurariu-Saia rational-point and Jacobian isomorphism tables','Ogg fixed-point formula','Complete small-class-number lists and Gonzalez-Rotger CM residue-field formula','Written residual Atkin-Lehner orbit doubling non-CM argument','Torelli and effective descent from trivial marked inertia','K3 Riemann-Roch and degree-one zero-cycle Picard-Brauer descent'],'boundary':'Outside the retained827-row catalogue. Full rational marking existence is proved by exact inputs and written arguments; period coordinates, K3 equation and MW17 are not computed. No rootless conclusion is asserted.'}
out=json.loads(json.dumps(out,default=int));a=argparse.ArgumentParser();a.add_argument('--check',action='store_true');args=a.parse_args();path=PACKET/'certificate.json'
if args.check:assert json.loads(path.read_text())==out
else:path.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print('PASS1020: full maximal510 marking curve, intrinsic non-CM period, primitive NS and actual divisor descent; rootless/MW17 unresolved.')
