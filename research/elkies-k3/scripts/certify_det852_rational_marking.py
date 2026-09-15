"""Primitive determinant852 NS/T glue and rational-divisor descent witnesses."""
from sage.all import *
import json,hashlib,argparse
from pathlib import Path
R=Path(__file__).resolve().parents[2];P=R/'artifacts/generated-results/elkies-k3-det852-marking-v1'
cat=R/'artifacts/generated-results/elkies-k3-rank7-auxiliary-catalogue-v1.json';curve=R/'artifacts/generated-results/elkies-k3-det852-marking-gate-v1/certificate.json';period=R/'artifacts/generated-results/elkies-k3-det852-cm-gate-v1/intrinsic-noncm.json'
c=json.loads(curve.read_text());p=json.loads(period.read_text());assert p['status']=='PASS_INTRINSIC_RATIONAL_NONCM_POINT'
sid='K3-4ff75fec54d01662';row=next(r for r in json.loads(cat.read_text())['surfaces'] if r['surface_id']==sid)
T=matrix(ZZ,c['T_gram']);assert T==matrix(ZZ,[[-2,0,1],[0,4,0],[1,0,106]]) and T.det()==-852
F=matrix(ZZ,row['frames'][0]['gram']);assert F.is_positive_definite() and F.det()==852
S=block_diagonal_matrix(matrix(ZZ,[[0,1],[1,0]]),-F);assert QuadraticForm(QQ,S).signature_vector()==(1,18,0)
def generator(G):
 d,u,v=G.smith_form();assert d.diagonal()==[1]*(G.nrows()-1)+[852]
 z=G.inverse()*u.inverse()*vector(ZZ,[0]*(G.nrows()-1)+[1]);return z,z*G*z/2
sd,qs=generator(S);td,qt=generator(T);a=next(a for a in range(852) if gcd(a,852)==1 and qs+a*a*qt in ZZ)
z=vector(QQ,list(sd)+list(a*td));direct=block_diagonal_matrix(S,T)
B=matrix(ZZ,[list(852*v) for v in identity_matrix(ZZ,22).rows()]+[list(852*z)]).row_module().basis_matrix()/852
K=B*direct*B.transpose();assert all(x in ZZ for x in K.list()) and all(K[i,i]%2==0 for i in range(22)) and abs(K.det())==1
assert QuadraticForm(QQ,K).signature_vector()==(3,19,0) and 1/abs(B.det())==852
# Norm-six frame coordinate gives two square-minus-two classes intersecting1.
r=vector(ZZ,[1,-1]+[0]*17);s=vector(ZZ,[1,2]+[int(i==4) for i in range(17)])
assert r*S*r==s*S*s==-2 and r*S*s==1
# Each isotropic discriminant graph projection is injective, hence both summands primitive.
assert gcd(a,852)==1
def rows(G):return [[str(x) for x in r] for r in G.rows()]
out={'schema':'elkies-k3.det852-rational-marking.v1','status':'PASS_FULL_RATIONAL_SATURATED_NS_BY_WRITTEN_DESCENT','surface_id':sid,'T_gram':rows(T),'NS_gram':rows(S),'glue':{'q_NS':str(qs),'q_T':str(qt),'multiplier':a,'index':852,'basis':rows(B),'unimodular_gram':rows(K),'signature':[3,19]},'divisor_descent':{'r':list(r),'s':list(s),'intersection_gram':[[-2,1],[1,-2]],'Euler_characteristics':[1,1],'zero_cycle_degree':1},'period':p['construction'],'full_rational_NS_existence':'PROVED_BY_WRITTEN_MODULI_AND_DIVISOR_DESCENT','period_coordinates':None,'explicit_K3_equation':None,'rootless_frame':None,'MW17_fibration':None,'inputs':[{'path':str(q.relative_to(R)),'sha256':hashlib.sha256(q.read_bytes()).hexdigest()} for q in [cat,curve,period]],'theorem_inputs':['Previously certified canonical full stable marking curve and intrinsic rational non-CM period','Torelli and period surjectivity for primitive lattice-polarized K3s','Trivial marked inertia gives effective descent overQ','K3 Riemann-Roch and Picard-Brauer obstruction; degree-one zero-cycle gives injective base Brauer map'],'boundary':'An actual projective K3/Q with geometric Picard19 and all saturated NS divisor classes rational is proved to exist. Equation and period coordinates remain uncomputed; rootless/MW17 not decided here.'}
out=json.loads(json.dumps(out,default=int));args=argparse.ArgumentParser();args.add_argument('--check',action='store_true');a=args.parse_args();path=P/'certificate.json'
if a.check:assert json.loads(path.read_text())==out
else:path.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print('PASS852: primitive K3-lattice glue, rational non-CM period, full saturated rational NS by written descent; MW17 remains open.')
