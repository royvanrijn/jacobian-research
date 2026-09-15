#!/usr/bin/env sage-python
"""Exact det388 full-period point and primitive marking data; descent is a written proof."""
import argparse,hashlib,json,re
from pathlib import Path
from sage.all import *
ROOT=Path(__file__).resolve().parents[2]
PACKET=ROOT/'artifacts/generated-results/elkies-k3-det388-marking-v1'
CAT=ROOT/'artifacts/generated-results/elkies-k3-rank7-auxiliary-catalogue-v1.json'
TINPUT=ROOT/'artifacts/generated-results/elkies-k3-rank7-t-arithmetic-v1.json'
OUTPUT=PACKET/'certificate.json'
sid='K3-f5168aef816b1b25'
row=next(r for r in json.loads(CAT.read_text())['surfaces'] if r['surface_id']==sid)
tr=next(r for r in json.loads(TINPUT.read_text())['surfaces'] if r['surface_id']==sid)
G=matrix(ZZ,tr['literal_transcendental_gram'])
assert G==matrix(ZZ,[[-2,-6,-1],[-6,2,1],[-1,1,10]]) and G.det()==-388
assert QuadraticForm(QQ,G).signature_vector()==(2,1,0)
C=CliffordAlgebra(QuadraticForm(QQ,G));v=C.gens()
B=[C.one(),v[0]*v[1],v[0]*v[2],v[1]*v[2]];keys=list(C.basis().keys());inds=[0,4,5,6]
def coords(x):return vector(QQ,[x.monomial_coefficients().get(keys[i],0) for i in inds])
def elem(c):return sum((a*b for a,b in zip(c,B)),C.zero())
def trace(x):return coords(x+x.clifford_conjugate())[0]
def strings(A):return [[str(x) for x in r] for r in A.rows()]
P=matrix(QQ,4,4,lambda i,j:trace(B[i]*B[j]))
assert P.det()==-194**2
assert P==matrix(QQ,tr['clifford']['integral_even_clifford_order']['reduced_trace_pairing'])
assert QuaternionAlgebra(QQ,10,QQ(97)/20).discriminant()==194
for b in B:
 for c in B:assert all(x in ZZ for x in coords(b*c))
vol=v[0]*v[1]*v[2]-G[0,1]/2*v[2]+G[0,2]/2*v[1]-G[1,2]/2*v[0]
images=[vol*x for x in v];M=matrix(QQ,[coords(x) for x in images]).transpose();left=(M.transpose()*M).inverse()*M.transpose()
D,U,V=G.smith_form();assert D.diagonal()==[1,1,388]
dual=G.inverse()*U.inverse()*vector(ZZ,[0,0,1]);q=dual*G*dual/2
ounits=[n for n in range(388) if gcd(n,388)==1 and (n*n-1)*q in ZZ]
assert ounits==[1,193,195,387]
records=[]
for label,cc in [(2,(-5,3,-2,-4)),(97,(403,56,67,-15)),(194,(-1,1,-6,2))]:
 x=elem(cc);norm=coords(x*x.clifford_conjugate())[0];assert norm==label
 xi=x.clifford_conjugate()/norm
 R=matrix(QQ,[coords(x*b*xi) for b in B]).transpose()
 A=matrix(QQ,[left*coords(x*b*xi) for b in images]).transpose()
 assert all(z in ZZ for z in R.list()+A.list()) and abs(R.det())==1
 assert A.transpose()*G*A==G and A.det()==1
 act=U*A.inverse().transpose()*U.inverse();mu=int(act[2,2]%388)
 assert mu=={2:195,97:193,194:387}[label]
 records.append({'label':label,'coordinates':list(cc),'norm':int(norm),'order_action':strings(R),'T_action':strings(A),'discriminant_multiplier':mu})
# Polynomial numerator coefficients prove stability of ALL local order units.
Gi=G.inverse();I=identity_matrix(QQ,3)
def diffnum(c):
 x=elem(c);bar=x.clifford_conjugate();norm=coords(x*bar)[0]
 return (matrix(QQ,[left*coords(x*b*bar) for b in images]).transpose()-norm*I)*Gi
basis=list(identity_matrix(QQ,4).rows());polys=[diffnum(c) for c in basis]
polys += [diffnum(basis[i]+basis[j])-diffnum(basis[i])-diffnum(basis[j]) for i in range(4) for j in range(i)]
assert all(z in ZZ for a in polys for z in a.list())
assert -matrix(QQ,records[-1]['T_action'])==I+matrix(QQ,[[1],[0],[0]])*matrix(QQ,[G.row(0)])
# Exact published model, with non-CM witness independent of a finite CM table.
raw=(PACKET/'published-model.m').read_text()
coefs=[-19,92,-286,592,-921,1016,-872,-460,1545,-1752,34,1752,1545,460,-872,-1016,-921,-592,-286,-92,-19]
published=re.search(r'Polynomial\(\[RationalField\(\) \| ([^\]]+)\]',raw).group(1)
assert list(map(int,published.split(',')))==coefs
f=PolynomialRing(QQ,'x')(coefs);assert gcd(f,f.derivative())==1 and f.degree()==20
value=f(2);assert value==-315538147 and value==-7411*42577
assert is_prime(7411) and is_prime(42577) and 7411>427 and (7411-1)/3>2
assert not value.is_square()
# Primitive NS/T gluing, using only the retained source frame after period admission.
F=matrix(ZZ,row['frames'][0]['gram']);S=block_diagonal_matrix(matrix(ZZ,[[0,1],[1,0]]),-F)
assert S.det()==388 and QuadraticForm(QQ,S).signature_vector()==(1,18,0)
def dualgen(A):
 dd,uu,vv=A.smith_form();assert abs(dd[-1,-1])==388 and all(abs(dd[i,i])==1 for i in range(A.nrows()-1))
 return A.inverse()*uu.inverse()*vector(ZZ,[0]*(A.nrows()-1)+[1])
sd=dualgen(S);td=dualgen(G);qs=sd*S*sd/2;qt=td*G*td/2
mult=next(u for u in range(388) if gcd(u,388)==1 and qs+u*u*qt in ZZ)
z=vector(QQ,list(sd)+list(mult*td));direct=block_diagonal_matrix(S,G)
change=matrix(ZZ,[list(388*w) for w in identity_matrix(QQ,22).rows()]+[list(388*z)]).row_module().basis_matrix()/388
K=change*direct*change.transpose()
assert all(a in ZZ for a in K.list()) and all(K[i,i]%2==0 for i in range(22))
assert abs(K.det())==1 and QuadraticForm(QQ,K).signature_vector()==(3,19,0)
assert 1/abs(change.det())==388
assert S[2:4,2:4]==matrix(ZZ,[[-2,1],[1,-2]])
result={'schema':'elkies-k3.det388-marking.v1','status':'PASS_FULL_MARKED_NONCM_PERIOD_AND_DESCENT_DATA',
 'surface_id':sid,'T_gram':strings(G),'NS_gram':strings(S),
 'inputs':[{'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in [CAT,TINPUT,PACKET/'published-model.m']],
 'Clifford_trace_pairing':strings(P),'maximal_order_reduced_discriminant':194,
 'normalizers':records,'discriminant_orthogonal_units':ounits,
 'local_unit_stability_polynomial_coefficients':[strings(a) for a in polys],
 'full_marked_curve':'X_0^194(1)/<w_194> = P1 over Q',
 'point':{'x':2,'projective_conic':[4,2,1],'lift_square':int(value),'ramified_primes':[7411,42577],'non_CM':True},
 'glue':{'NS_discriminant_q':str(qs),'T_discriminant_q':str(qt),'multiplier':mult,'basis':strings(change),'unimodular_Gram':strings(K),'index':388,'signature':[3,19]},
 'rational_divisor_descent_pair':[[-2,1],[1,-2]],
 'theorem_inputs':['Maximal quaternion-order normalizer and local reduced-norm surjectivity','Canonical arithmetic K3 period map and Torelli on a fixed ample chamber','Guo-Yang canonical model of X_0^194(1)','Gonzalez-Rotger Theorem5.8: K(P)=ring class field','Complete class-number-one/two imaginary quadratic field classification: absolute discriminant<=427','Ring-class-number and ramification formulas','Trivial marked automorphisms give effective descent; line-bundle obstruction index divides Euler characteristic'],
 'full_rational_NS_existence':'PROVED_BY_WRITTEN_DESCENT_ARGUMENT',
 'explicit_K3_equation':None,'MW17_fibration':None,'rootless_frame':None,
 'boundary':'Exact period and lattice data plus written arithmetic descent. No explicit K3 equation, rootless frame, or17 displayed sections.'}
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');args=parser.parse_args()
if args.check:assert json.loads(OUTPUT.read_text())==result
else:OUTPUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print('PASS det388: full marked rational non-CM period x=2; primitive NS glue and rational-divisor descent pair; MW17 UNKNOWN')
