"""Exact rational Inose model at the published level311 non-CM point; rationalizing twist is existential."""
from sage.all import *
from pathlib import Path
import json,hashlib,argparse
R=Path(__file__).resolve().parents[2];P=R/'artifacts/generated-results/elkies-k3-det622-inose-source-v1';P.mkdir(exist_ok=True)
A=ZZ(31244183594433270730990985793058589729152601677824000000)
B=ZZ(1565810538998051715397339689492195035077551267840000);D=ZZ(39816211853)
assert is_prime(311) and D==11*17*9011*23629 and all(is_prime(p) for p in [11,17,9011,23629])
PR=PolynomialRing(QQ,names=('W','X','Y','Z'));W,X,Y,Z=PR.gens()
f=X**2+W*Y-2*X*Y+2*Y**2+7*X*Z-8*Y*Z+13*Z**2
g=W*X**2-2*W*X*Y+X**2*Y-W*Y**2-X*Y**2-2*Y**3+W**2*Z+6*W*X*Z-X**2*Z-W*Y*Z+5*X*Y*Z+4*Y**2*Z+7*W*Z**2-4*X*Z**2-2*Z**3
point=(6,8,-1,-2);assert f(*point)==g(*point)==0
jac=matrix(QQ,[[h.derivative(v)(*point) for v in PR.gens()] for h in [f,g]]);assert jac.rank()==2
I=QQ(A*A-B*B*D)/1728**2;J=QQ((A-1728)**2-B*B*D)/1728**2
assert I and J and B and I-J+1==QQ(2*A)/1728
S=PolynomialRing(QQ,'t');t=S.gen();a=-3*I*J;c=-2*I*J**2;d=I**2*J**3
a4=a*t**4;a6=t**5*(t**2+c*t+d)
residual=(t**2+c*t+d)**2-4*(I*J)**3*t**2
assert residual.degree()==4 and gcd(residual,residual.derivative())==1 and residual(0)!=0
assert -16*(4*a4**3+27*a6**2)==-432*t**10*residual
assert -a**3/(27*d)==I and c*c/(4*d)==J
assert (I+J-1)**2-4*I*J==QQ(4*B*B*D)/1728**2
E8=matrix(ZZ,CartanMatrix(['E',8]).rows());U=matrix(ZZ,[[0,1],[1,0]])
F=block_diagonal_matrix(E8,E8,matrix(ZZ,[[622]]));NS=block_diagonal_matrix(U,-F);T=block_diagonal_matrix(U,matrix(ZZ,[[622]]))
assert F.det()==622 and NS.det()==622 and T.det()==-622
assert QuadraticForm(QQ,NS).signature_vector()==(1,18,0)
assert NS.smith_form()[0].diagonal()==[1]*18+[622]
def rows(G):return [[int(x) for x in r] for r in G.rows()]
source=R/'artifacts/generated-results/elkies-k3-fricke-source-gates-v1/prime-genus4-6.pdf'
out={'schema':'elkies-k3.det622-inose-source.v1','status':'PASS_RATIONAL_INOSE_MODEL_AND_EXISTENTIAL_FULL_MARKING_TWIST','level':311,'point':list(point),'point_smooth':True,'canonical_curve_equations':[str(f),str(g)],'j_pair':{'trace_half':str(A),'sqrt_coefficient':str(B),'radicand':int(D)},'I':str(I),'J':str(J),'model':{'equation':'y^2 = x^3 + a*t^4*x + t^5*(t^2+c*t+d)','a':str(a),'c':str(c),'d':str(d)},'residual_discriminant_coefficients':[str(v) for v in residual.list()],'fibres':['II* at0','II* atinfinity','fourI1 geometrically'],'NS_gram':rows(NS),'T_gram':rows(T),'geometric_picard_rank':19,'geometric_MW_rank':1,'primitive_generator_height':622,'rationalizing_constant_twist_exists':True,'rationalizing_twist_squareclass':None,'full_rational_NS_exists':True,'rootless_frame':None,'MW17_equation':None,'theorem_inputs':['Adzaga etal Theorem1.1 and Section5.1: point[6:8:-1:-2] on X0+(311) corresponds to the displayed non-CM cyclic311-isogenous j pair','Inose normal form and Shioda isometry MW(F1)=Hom(E1,E2)<2>','Written rank-one Galois-character twisting argument; affineE8 has no graph automorphism'],'published_model_source':{'path':str(source.relative_to(R)),'sha256':hashlib.sha256(source.read_bytes()).hexdigest()},'boundary':'An explicit rational untwisted model and proof that one constant quadratic twist has full saturated rational NS. The twist squareclass and section are not yet computed; rootless existence and MW17 remain UNKNOWN.'}
p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args();output=P/'certificate.json'
if args.check:assert json.loads(output.read_text())==out
else:output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print('PASS det622: exact Inose model, published non-CM311 source, full rational marking after an existential constant twist. Twist value and MW17 UNKNOWN.')
