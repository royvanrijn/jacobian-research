from sage.all import *
from pathlib import Path
import json, argparse, hashlib
P=Path(__file__).resolve().parents[2]/'artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1';x=json.loads((P/'rank1-frames.json').read_text());h=next(h for h in x['hits'] if h['root_type']=='A2+A2+A3+A9');F=matrix(ZZ,h['frame_gram']);R=F[:16,:16];r=F[:16,16];height=F[16,16]-(r.transpose()*R.inverse()*r)[0,0]
assert R==block_diagonal_matrix([CartanMatrix(['A',n]) for n in [2,2,3,9]])
assert F.det()==1020 and F.is_positive_definite() and height==QQ(17)/6 and F[16,16]==6
assert F.elementary_divisors()==[1]*16+[1020]
D,U,V=F.smith_form();g=V.inverse().row(16)*F.inverse();g=vector(QQ,[z-floor(z) for z in g]);q=g*F*g
assert lcm(z.denominator() for z in g)==1020
u=next(u for u in range(1020) if gcd(u,1020)==1 and (q+u*u*QQ(403)/1020)/2 in ZZ)
roots=F.__pari__().qfminim(2,10000,2);assert roots[0]==114
out={'status':'PASS_ALTERNATIVE_FRAME_ARITHMETIC_NOT_EQUATION','frame_gram':[[int(z) for z in a] for a in F.rows()],'root_type':h['root_type'],'component_indices':[0,1,0,5],'height':str(height),'P_dot_O':1,'root_count':114,'q_generator':str(q),'discriminant_generator':[str(z) for z in g],'isometry_multiplier_to_T':u,'construction_tradeoff':'The I10 central component can use local valuations ord(B)>=5 and ord(C)>=5 in y^2=x^3+A*x^2+B*x+C^2 with P=(0,C). This avoids the order9 cancellation needed by the previous component1 presentation, but P.O=1 introduces one pole. The six-parameter chart below incorporates I10 and one I3; the full source equation is not certified.'}

z=json.loads((P/'center-chart.json').read_text());RR=PolynomialRing(QQ,names=z['variables']);T=PolynomialRing(RR,'t');t=T.gen();dd,ee,*rest=RR.gens()
def read(ts):return sum(QQ(c)*prod(v**k for v,k in zip(RR.gens(),m)) for m,c in ts)
AA=T([read(a) for a in z['A']]);CC=T([read(a) for a in z['C']]);ff=T([read(a) for a in z['f']]);gg=T([read(a) for a in z['g']]);BB=AA*AA/3+ff
assert CC*CC==(AA/3)**3+ff*AA/3+gg
assert AA.degree()==6 and AA[6]==3 and CC.degree()==9 and CC[9]==1
assert ff.degree()<=8 and gg.degree()<=12 and AA(0)==dd*dd and AA(1)==ee*ee
assert BB.valuation(t)>=5 and CC.valuation(t)==5 and BB(1)==CC(1)==0
assert BB.derivative()(1)+2*ee*CC.derivative()(1)==0
Delta=-4*ff**3-27*gg**2;assert Delta.valuation(t)==10 and Delta.valuation(t-1)==3 and Delta.degree()==24
sample={RR.gen(i):v for i,v in enumerate([1,2,0,0,0,0])};Q=PolynomialRing(QQ,'t');tt=Q.gen()
fs=Q([a.subs(sample) for a in ff.list()]);gs=Q([a.subs(sample) for a in gg.list()]);ds=-4*fs**3-27*gs**2;res=ds//(tt**10*(tt-1)**3)
assert ds.degree()==24 and res.degree()==11 and res(0)!=0 and res(1)!=0 and res.gcd(res.derivative())==1
out['status']='PASS_CENTER_FRAME_AND_PARTIAL_RATIONAL_CHART'
fe=T([a.subs({ee:0}) for a in ff.list()]);ge=T([a.subs({ee:0}) for a in gg.list()]);ce=T([a.subs({ee:0}) for a in CC.list()]);de=-4*fe**3-27*ge**2
assert fe.valuation(t-1)==2 and ge.valuation(t-1)==2 and de.valuation(t-1)==4
assert ge.derivative(2)(1)/2==ce.derivative()(1)**2
out['nonidentity_IV_branch']={'condition':'e=0','orders_f_g_D_at_one':[2,2,4],'split_witness':'coefficient of (t-1)^2 in g equals C_prime(1)^2','boundary':'Generic branch only; further specializations and remaining source fibers are not certified.'}
out['chart_sha256']=hashlib.sha256((P/'center-chart.json').read_bytes()).hexdigest()
out['partial_chart']={'variables':z['variables'],'equation':'y^2=x^3+f*x+g; P=(A/3,C)','built_fibers':'generic split I10 at0 and split I3 at1','section_components':[5,1],'P_dot_O':1,'sample_parameters':[1,2,0,0,0,0],'sample_f':str(fs),'sample_g':str(gs),'sample_other_fibers':'eleven I1','remaining':'An additional identity I3 or IV and an identity I4 are not constructed. No determinant1020 equation, Picard19 or MW17 equation is asserted.'}
a=argparse.ArgumentParser();a.add_argument('--check',action='store_true');args=a.parse_args();target=P/'center-source-certificate.json'
if args.check:assert json.loads(target.read_text())==out
else:target.write_text(json.dumps(out,indent=2)+'\n')
print(out['status'])
