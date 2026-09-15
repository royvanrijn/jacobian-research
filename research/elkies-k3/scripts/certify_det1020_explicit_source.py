"""Exact equation, full split NS marking and two-prime Picard certificate."""
from sage.all import *
from pathlib import Path
import json, argparse, hashlib, importlib.util
P=Path(__file__).resolve().parents[2]/'artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1'
spec=importlib.util.spec_from_file_location('counts',P/'count_split_sources_python.py');counts=importlib.util.module_from_spec(spec);spec.loader.exec_module(counts)
z=json.loads((P/'split-specializations.json').read_text())['records'][1]
assert z['multiple']==2 and z['parameter']=='203702/35733'
R=PolynomialRing(QQ,'t');t=R.gen()
f,g,A,C=[R(z[key]) for key in ['f','g','A','C']]
l,m=QQ(z['lambda']),QQ(z['mu']);locs=[QQ(0),QQ(1),l,m];orders=[10,3,3,4]
assert len(set(locs))==4 and f.degree()==8 and g.degree()==12
assert A.degree()==6 and A[6]==3 and C.degree()==9 and C[9]==1
assert C*C==(A/3)**3+f*A/3+g
B=A*A/3+f;D=-4*f**3-27*g*g
assert B.valuation(t)>=5 and C.valuation(t)==5
assert C(1)==B(1)==0
fibers=[]
for v,n in zip(locs,orders):
 assert D.valuation(t-v)==n and f(v)!=0
 node=-3*g(v)/(2*f(v));tangent=3*node
 assert 3*node**2+f(v)==0 and node**3+f(v)*node+g(v)==0
 assert tangent!=0 and tangent.is_square()
 assert (A(v)/3==node)==(v in [0,1])
 fibers.append({'location':str(v),'type':'I'+str(n),'node':str(node),'tangent':str(tangent),'tangent_root':str(tangent.sqrt())})
res,rem=D.quo_rem(prod((t-v)**n for v,n in zip(locs,orders)))
assert not rem and res.degree()==4 and res.gcd(res.derivative())==1
assert all(res(v)!=0 for v in locs) and D.degree()==24
frames=json.loads((P/'rank1-frames.json').read_text())
F=matrix(ZZ,next(h['frame_gram'] for h in frames['hits'] if h['root_type']=='A2+A2+A3+A9'))
assert F.det()==1020 and F.elementary_divisors()==[1]*16+[1020]
root=F[:16,:16];pair=F[:16,16]
assert F[16,16]-(pair.transpose()*root.inverse()*pair)[0,0]==QQ(17)/6
assert root==block_diagonal_matrix([CartanMatrix(['A',n]) for n in [2,2,3,9]])
assert list(pair.column(0))==[0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0]
_,_,V=F.smith_form();gen=V.inverse().row(16)*F.inverse();q=gen*F*gen
assert (510**2*q)/2 not in ZZ # sole possible nontrivial even-overlattice index is2
reductions=[]
cpp={tuple(map(int,line.split()[:2])):list(map(int,line.split()[2:])) for line in (P/'count-split-output.txt').read_text().splitlines()}
for p in [41,59]:
 k=GF(p);S=PolynomialRing(k,'t');ff=S(f);gg=S(g);rr=S(res);vv=list(map(k,locs));dd=-4*ff**3-27*gg**2
 assert len(set(vv))==4 and dd.degree()==24 and rr.degree()==4 and rr.gcd(rr.derivative()).degree()==0
 assert all(rr(v)!=0 and ff(v)!=0 for v in vv)
 assert all(dd.valuation(S.gen()-v)==n for v,n in zip(vv,orders))
 assert all(k(QQ(h['tangent']))!=0 and k(QQ(h['tangent'])).is_square() for h in fibers)
 cf=list(map(int,ff.list()));cg=list(map(int,gg.list()))
 N,N2=[ZZ(counts.count(p,j,cf,cg)) for j in [1,2]]
 assert [N,N2]==cpp[(2,p)]
 s=N-1-p*p-19*p;b=N2-1-p**4-18*p*p
 sols=[(e,s-e*p) for e in [-1,1] if (s-e*p)**2==b];assert len(sols)==1
 eps,a=sols[0];assert a not in [-2*p,-p,0,p,2*p]
 disc=a*a-4*p*p
 reductions.append({'p':p,'N_p':int(N),'N_p_squared':int(N2),'T_trace':int(s),'linear_sign':int(eps),'quadratic_trace':int(a),'geometric_reduction_Picard':20,'Artin_Tate_squareclass':int(disc.squarefree_part())})
assert len(set(r['Artin_Tate_squareclass'] for r in reductions))==2
out={'status':'PASS_EXPLICIT_Q_SOURCE_FULL_NS_DET1020_PICARD19','parameter':z['parameter'],'equation':{'form':'y^2=x^3+f(t)*x+g(t)','f':z['f'],'g':z['g']},'section':{'x':[str(c/3) for c in A.list()],'y':z['C'],'P_dot_O':1,'components_at_0_1_lambda_mu':[5,1,0,0],'height':'17/6'},'fibers':fibers,'residual_discriminant':[str(c) for c in res.list()],'Picard_rank':19,'NS_determinant':-1020,'MW_rank_this_fibration':1,'MW_torsion':0,'full_rational_marking':'O, fiber, all16 nonidentity split components and P generate geometric NS. No even index2 overlattice exists.','frame_gram':[[int(c) for c in row] for row in F.rows()],'reductions':reductions,'written_inputs':['Split multiplicative fiber resolution and section component/height formulas','Specialization of NS, Tate theorem and Artin-Tate for K3 over finite fields','Two distinct reduction discriminant squareclasses bound geometric Picard rank by19'],'boundary':'Explicit fully marked source of determinant1020. Rootless-fibration existence is separate; no physical U transport, MW17 equation or17 section coordinates are computed.'}
a=argparse.ArgumentParser();a.add_argument('--check',action='store_true');args=a.parse_args();target=P/'explicit-source-certificate.json'
if args.check:assert json.loads(target.read_text())==out
else:target.write_text(json.dumps(out,indent=2)+'\n')
print(out['status'])
