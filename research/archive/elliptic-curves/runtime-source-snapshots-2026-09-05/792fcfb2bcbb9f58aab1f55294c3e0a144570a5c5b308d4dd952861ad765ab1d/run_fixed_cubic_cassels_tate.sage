#!/usr/bin/env sage
"""Generate exact explicit-cover Cassels--Tate witnesses for the u=-1 subspace.

No class group or ambient Selmer basis is requested. The public certificate
is independently replayed by verify_fixed_cubic_cassels_tate.sage. This
producer uses PARI qfsolve, qfparam, hyperellminimalmodel, hyperellred,
nfroots, and hilbert. Rational local witnesses are accepted only after
exact squareclass checks. Failed or interrupted entries remain UNKNOWN.

The initial run is limited to 18 basis covers, 153 independent pairings,
and the three nonzero radical covers. Every cover and completed pair is
checkpointed. Use a fresh --workdir to reproduce the arithmetic; existing
checkpoints are retained. Regeneration is separate from the cheap replay.

Fisher (2022), Theorem 3.1 and Remark 3.3:
https://antsmath.org/ANTSXV/papers/ANTS-XV_fisher.pdf
"""
from sage.all import *
import json,sys,time,signal
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
from fixed_cubic_field_curve_family import field_multiply,field_product
P=ROOT/'artifacts/local/fixed-cubic-u-minus1-ct'
p=json.load(open(ROOT/'artifacts/generated-results/elliptic-curves/fixed_cubic_field_fermigier_rank20_local_kummer_u2_v1.json'));r=next(row for row in p['runs'] if row['parameter_u']=='-1')
A,B=map(QQ,(p['anchor']['base_polynomial_ascending'][1],p['anchor']['base_polynomial_ascending'][0]))
E=EllipticCurve(list(map(QQ,r['raw_curve_ainvariants'])));I0=E.c4()/16;J0=E.c6()/32
R=PolynomialRing(QQ,'a,b,c');a,b,c=R.gens();T=PolynomialRing(QQ,'x');x=T.gen()
red=pari('(q)->{my(m);my(z=hyperellred(q,&m));[z,m]}')
mini=pari('(q)->{my(m);my(z=hyperellminimalmodel(q,&m));[z,m]}')
def ij(q):
 e,d,c,b,a=q.list();return 12*a*e-3*b*d+c*c,72*a*c*e-27*a*d*d-27*b*b*e+9*b*c*d-2*c**3

def quartic(mask):
 path=P/f'cover-{mask}.json'
 if path.exists():return json.loads(path.read_text())
 start=time.monotonic(); print('COVER_START',mask,flush=True)
 bs=[list(map(QQ,z)) for i,z in enumerate(p['anchor']['known_kummer_basis_beta_power_coordinates']) if mask>>i&1]
 beta=field_product(bs,A,B)
 v=field_multiply(beta,field_multiply((a,b,c),(a,b,c),A,B),A,B)
 Gs=[matrix(QQ,3,3,lambda i,j:v[k].derivative(R.gen(i)).derivative(R.gen(j))/2) for k in (1,2)]
 G=Gs[0]+Gs[1];G=matrix(ZZ,G/QQ(pari.content(G)));fac=G.det().factor()
 fmat=pari.matrix(len(fac),2,[z for pr,e in fac for z in [pr,e]])
 sol=pari.qfsolve([G,fmat]); assert sol.type()=='t_COL'
 par=pari.qfparam(G,sol,1);par=matrix(QQ,par/pari.content(par))
 gm=par*vector([x*x,x,1]);q=-(gm*Gs[0]*gm);assert gm*G*gm==0
 den=q.denominator();q=q*den**2
 assert q.degree()==4
 raw=q;Mtot=identity_matrix(QQ,2)
 for fun in [red,mini,red]:
  z,m=fun(q);q=T(z[0])+T(z[1])**2/4
  Mtot=Mtot*matrix(QQ,m[1])
 ii,jj=ij(q);s=QQ(I0/ii).nth_root(4);assert jj*s**6==J0
 q=q*s*s;assert ij(q)==(I0,J0)
 if q.denominator()!=1: print('RATIONAL',mask,flush=True)
 u,vv=Mtot*vector([x,1]);pull=sum(raw[k]*u**k*vv**(4-k) for k in range(5))
 ratio=QQ(pull.leading_coefficient()/q.leading_coefficient());k=ratio.sqrt();assert k in QQ and pull==k*k*q
 rec={'anchor_mask':mask,'beta':[str(z) for z in beta],'conic_matrix':[[str(z) for z in row] for row in G.rows()], 'conic_point':[str(z) for z in sol], 'parametrization':[[str(z) for z in row] for row in par.rows()], 'initial_y_denominator':str(den),'parameter_transform':[[str(z) for z in row] for row in Mtot.rows()], 'd_over_quartic_y':str(k/den),'quartic':[str(q[k]) for k in range(5)],'seconds':time.monotonic()-start}
 path.write_text(json.dumps(rec,indent=2)+'\n');print('COVER_DONE',mask,rec['seconds'],'digits',max(len(str(z)) for z in q.list()),flush=True);return rec

from collections import deque
Y=PolynomialRing(QQ,'theta'); theta=Y.gen();K=NumberField(theta**3+A*theta+B,'th');th=K.gen();ph=-3*(th-th**2)-E.a2()
V=matrix(QQ,3,3,lambda i,j:(ph**j)[i]);Vinv=V.inverse()
KT=PolynomialRing(K,'x');xx=KT.gen()
pol=pari(f'y^3+({A})*y+({B})');nf=pari.nfinit(pol)

def sqp(a,pr):
 pr=ZZ(pr)
 if not a:return True
 v=a.valuation(pr)
 if v%2:return False
 u=a/pr**v
 return (ZZ(u.numerator()*inverse_mod(u.denominator(),8))%8==1) if pr==2 else kronecker(ZZ(u.numerator()*u.denominator()),pr)==1

def witness(q,pr):
 q=T(q);den=q.denominator();q=ZZ(den)**2*q
 # Every return is a rational x giving an actual Q_p square, verified below.
 if pr=='infinity':
  for j in range(1001):
   for a0 in (j,-j):
    if q(a0)>0:return QQ(a0)
  roots=q.roots(AA,multiplicities=False)
  for lo,hi in zip(roots,roots[1:]):
   a0=QQ(((lo+hi)/2).n(128))
   if q(a0)>0:return a0
  raise RuntimeError('no real witness')
 # Process both projective charts. The reciprocal chart omits t=0;
 # any local neighbourhood of that point contains t!=0.
 queue=deque([(q,ZZ(0),ZZ(1),False),(q.reverse(),ZZ(0),ZZ(1),True)])
 for count in range(10000):
  if not queue:raise RuntimeError('empty local tree')
  g,shift,scale,recip=queue.popleft()
  val=min(z.valuation(pr) for z in g if z)
  gn=g/pr**(2*(val//2))
  picks=range(2) if pr==2 else range(min(pr,64))
  for t0 in picks:
   a0=shift+scale*t0
   if recip and not a0:
    for exponent in range(1,21):
     nearby=ZZ(pr)**exponent
     if gn(nearby) and sqp(QQ(gn(nearby)),pr):
      ans=1/QQ(scale*nearby)
      assert q(ans) and sqp(QQ(q(ans)),pr)
      return ans
    continue
   if gn(t0) and sqp(QQ(gn(t0)),pr):
    ans=1/QQ(a0) if recip else QQ(a0)
    assert q(ans) and sqp(QQ(q(ans)),pr)
    return ans
  if pr==2:
   nexts=[0,1]
  else:
   gg=PolynomialRing(GF(pr),'t')(gn/pr**(min(z.valuation(pr) for z in gn if z)))
   nexts=[ZZ(z) for z in gg.roots(multiplicities=False)]
  for t0 in nexts:
   gg=gn(pr*x+t0)
   # If the constant's squareclass is stable on this disc, it cannot help.
   if gg[0] and all(not z or z.valuation(pr)>gg[0].valuation(pr)+(2 if pr==2 else 0) for z in gg.list()[1:]):
    if not sqp(QQ(gg[0]),pr):continue
   queue.append((gg,shift+scale*t0,scale*pr,recip))
 raise RuntimeError('local witness node cap')

def field_sqrt(z):
 zp=pari.Mod(pari([QQ(z[i]) for i in range(3)]).Polrev('y'),pol)
 rt=pari.nfroots(nf,pari('x')**2-zp)
 if not len(rt):raise ValueError('not a square')
 ans=K(list(QQ(pari.lift(rt[0]).polcoef(i)) for i in range(3)))
 assert ans*ans==z
 return ans

def data(q):
 e,d,c,b,a=q.list();z=(4*a*ph+3*b*b-8*a*c)/3
 h=T([3*d*d-8*c*e,4*(c*d-6*b*e),2*(2*c*c-24*a*e-3*b*d),4*(b*c-6*a*d),3*b*b-8*a*c])
 H=(4*ph*KT(q.derivative(2))+KT(h.derivative(2)))/36+QQ(2)/9*(I0-ph**2)
 assert H[2]==z
 return z,H

def pairing(i,j):
 path=P/f'pair-{i}-{j}.json'
 if path.exists():return json.loads(path.read_text())
 t=time.monotonic();masks=[r['W_u_basis'][k]['mask'] for k in (i,j)];mask=masks[0]^masks[1]
 cs=[quartic(m) for m in masks+[mask]];qs=[T(list(map(QQ,c['quartic']))) for c in cs]
 z1,H=data(qs[0]);z2,_=data(qs[1]);z3,_=data(qs[2]);m=field_sqrt(z1*z2*z3)
 G=(m/z1)*H;gam=T([(Vinv*vector(QQ,[G[k][n] for n in range(3)]))[2] for k in range(3)])
 gam=gam/QQ(pari.content(gam));assert gam
 if gam.leading_coefficient()<0:gam=-gam
 a0=qs[1][4]
 fac=abs(ZZ(a0)).factor();places=sorted(set(r['complete_finite_place_support']+[2,3,5,7]+[int(z) for z,e in fac]))
 # Check the complete discriminant support required by Fisher Remark 3.3.
 DD=qs[0].discriminant();DD=abs(ZZ(DD))
 for pr in places:
  while DD%pr==0:DD//=pr
 assert DD==1
 locals=[];sign=1
 for pr in places+['infinity']:
  w=witness(qs[0],pr)
  if gam(w)==0:raise RuntimeError('zero gamma at witness')
  s=(-1 if a0<0 and gam(w)<0 else 1) if pr=='infinity' else int(pari.hilbert(a0,gam(w),pr))
  sign*=s
  locals.append({'place':pr,'x':str(w),'q_value':str(qs[0](w)),'gamma_value':str(gam(w)),'hilbert_symbol':s})
 rec={'i':i,'j':j,'basis_anchor_masks':masks,'sum_anchor_mask':mask,'square_root_phi_coefficients':[str(z) for z in Vinv*vector(QQ,list(m))], 'gamma':[str(gam[k]) for k in range(3)], 'leading_coefficient_factorization':[[str(z),int(e)] for z,e in fac], 'local_terms':locals,'value':int(sign==-1),'seconds':time.monotonic()-t}
 path.write_text(json.dumps(rec,indent=2)+'\n');print('PAIR_DONE',i,j,rec['value'],'seconds',rec['seconds'],flush=True);return rec


def main():
 global P
 import argparse
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--workdir',type=Path,default=P)
 parser.add_argument('--entry-seconds',type=int,default=60)
 parser.add_argument('--campaign-seconds',type=int,default=900)
 parser.add_argument('--point-height',type=int,default=10000)
 parser.add_argument('--write-certificate',action='store_true',
                     help='write the independently verified compact repository certificate')
 args=parser.parse_args()
 if min(args.entry_seconds,args.campaign_seconds,args.point_height)<=0:
  parser.error('all bounds must be positive')
 P=args.workdir
 P.mkdir(parents=True,exist_ok=True)
 pari.setrand(1)
 def expired(signum,frame):
  raise TimeoutError('per-entry resource limit')
 signal.signal(signal.SIGALRM,expired)
 start=time.monotonic()
 completed=[]
 for i in range(18):
  for j in range(i+1,18):
   if time.monotonic()-start>args.campaign_seconds:
    raise SystemExit('campaign resource limit; uncomputed entries UNKNOWN')
   signal.alarm(args.entry_seconds)
   try:
    completed.append(pairing(i,j))
   except Exception as exc:
    print('PAIR_UNKNOWN',i,j,type(exc).__name__,str(exc),flush=True)
   finally:
    signal.alarm(0)
 if len(completed)!=153:
  raise SystemExit('incomplete matrix; no radical claim')
 M=matrix(GF(2),18)
 for row in completed:
  M[row['i'],row['j']]=M[row['j'],row['i']]=row['value']
 # Check cached entries before using their kernel to schedule point search.
 import importlib.machinery, importlib.util, hashlib
 loader=importlib.machinery.SourceFileLoader('ct_verifier',str(ROOT/'elliptic-curves/cas/verify_fixed_cubic_cassels_tate.sage'))
 spec=importlib.util.spec_from_loader('ct_verifier',loader)
 verifier=importlib.util.module_from_spec(spec)
 loader.exec_module(verifier)
 source=ROOT/'artifacts/generated-results/elliptic-curves/fixed_cubic_field_fermigier_rank20_local_kummer_u2_v1.json'
 evidence={'status':'COMPLETE','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
           'covers':[json.loads(f.read_text()) for f in P.glob('cover-*.json')],
           'pairings':completed,'cross_checks':[]}
 verifier.verify(evidence)
 (P/'matrix.json').write_text(json.dumps([[int(z) for z in row] for row in M.rows()])+'\n')
 print('MATRIX_COMPLETE','rank',M.rank(),'radical',M.right_kernel().dimension(),flush=True)
 rad=list(M.right_kernel().basis())
 if len(rad)>2:
  raise SystemExit('radical exceeds declared three-cover pilot; no larger search launched')
 rows=[]
 for bits in range(1,2**len(rad)):
  v=sum((rad[i] for i in range(len(rad)) if bits>>i&1),M.row_ambient_module().zero())
  mask=0
  for i,bit in enumerate(v):
   if bit:mask^=r['W_u_basis'][i]['mask']
  signal.alarm(args.entry_seconds)
  try:
   cover=quartic(mask)
   q=T(list(map(QQ,cover['quartic'])))
   con=ZZ(pari.content(q))
   scale=prod(pr**(ex//2) for pr,ex in con.factor())
   q=q/scale**2
   points=pari.hyperellratpoints(q,args.point_height)
   row={'W_coordinates':[int(z) for z in v], 'cover':cover,
        'search':{'bound':args.point_height,'quartic':[str(z) for z in q],
                  'quartic_y_rescaling':str(scale),'raw_points':[[str(z) for z in pt] for pt in points],
                  'rational_points_at_infinity':[] if not q[4].is_square() else [str(q[4].sqrt())],
                  'status':'BOUNDED_SEARCH_COMPLETE','rank_claim':None}}
   rows.append(row)
   print('RADICAL_SEARCH',mask,'points',len(points),flush=True)
  finally:
   signal.alarm(0)
 (P/'radical-search.json').write_text(json.dumps(rows,indent=2)+'\n')
 if not args.write_certificate:
  return
 if args.point_height!=10000:
  raise SystemExit('the v1 compact certificate uses the declared height 10000')
 # Arithmetic controls use fresh local evaluations, with a separate cap
 # of eight checks and sixty seconds per check.
 checks=[]
 for i,j in [(6,0),(17,0),(3,2),(10,8)]:
  signal.alarm(args.entry_seconds)
  try:
   row=pairing(i,j)
   row['W_vectors']=[[int(k==i) for k in range(18)],
                     [int(k==j) for k in range(18)]]
   checks.append(row)
  finally:
   signal.alarm(0)
 for k,v in enumerate(rad):
  mask=verifier.coordinates_record(v,r['W_u_basis'][:18])['anchor_mask']
  r['W_u_basis'].append({'mask':mask})
  for j in [0,5]:
   signal.alarm(args.entry_seconds)
   try:
    row=pairing(18+k,j)
    row['W_vectors']=[[int(z) for z in v],[int(i==j) for i in range(18)]]
    checks.append(row)
   finally:
    signal.alarm(0)
 control_quartics=[T([-64,-164,-52,68,-11]),T([-3,-52,-232,-60,-4]),
                    T([-53,102,32,-78,-31])]
 root=['936032/9','-8656/9','20/9']
 gamma=verifier.fisher_gamma(control_quartics,root,44608,18842960)
 terms=[]
 for prime in [2,3,5,7,571,'infinity']:
  xx=witness(control_quartics[0],prime)
  terms.append({'place':prime,'x':str(xx),'q_value':str(control_quartics[0](xx)),
                'gamma_value':str(gamma(xx)),
                'hilbert_symbol':verifier.hilbert_symbol(-4,gamma(xx),prime)})
 control={'name':'Fisher 2022 Example 3.4, 571a1','I':'44608','J':'18842960',
          'quartics':[[str(z) for z in q] for q in control_quartics],
          'square_root_phi_coefficients':root,'gamma':[str(gamma[k]) for k in range(3)],
          'local_terms':terms,'value':1}
 from sage.version import version as sage_version
 evidence.update({'schema':'elliptic-curves.fixed-cubic-cassels-tate-evidence.v1',
                  'covers':[json.loads(f.read_text()) for f in sorted(P.glob('cover-*.json'))],
                  'cross_checks':checks,'published_control':control,'radical_search':rows,
                  'software':{'sage':sage_version,'pari':str(pari.version())},
                  'limits':{'main_pair_count':153,'entry_seconds':args.entry_seconds,
                            'main_campaign_seconds':args.campaign_seconds,
                            'additional_pair_checks':8,'radical_cover_count':3,
                            'point_height':10000}})
 result=verifier.verify(evidence)
 import gzip
 encoded=(json.dumps(evidence,sort_keys=True,separators=(',',':'))+'\n').encode()
 verifier.EVIDENCE.write_bytes(gzip.compress(encoded,mtime=0))
 summary={'schema':'elliptic-curves.fixed-cubic-cassels-tate.v1',
          'status':'PASS_EXACT_RESTRICTED_CASSELS_TATE', 'parameter_u':'-1',
          'source_sha256':verifier.digest(source),
          'evidence':str(verifier.EVIDENCE.relative_to(ROOT)),
          'evidence_sha256':verifier.digest(verifier.EVIDENCE),
          'software':evidence['software'],'limits':evidence['limits'],
          'arithmetic':result,
          'claim_boundary':[
           'This is the pairing on the certified 18-dimensional subspace, not on a complete Selmer basis.',
           'Every class outside the restricted radical maps nontrivially to Sha[2].',
           'The three nonzero radical classes remain UNKNOWN: point classes or further Sha obstructions.',
           'The dimension-two bound concerns rational Kummer classes inside W; it is not a full-curve rank upper bound.',
           'Coefficient height schedules reduced-cover searches and is not a solubility criterion.'],
          'source_hashes':{str(path.relative_to(ROOT)):verifier.digest(path) for path in [
           Path(__file__),ROOT/'elliptic-curves/cas/verify_fixed_cubic_cassels_tate.sage',
           ROOT/'elliptic-curves/tests/test_fixed_cubic_cassels_tate.py']}}
 verifier.SUMMARY.write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
 print('CERTIFICATE_WRITTEN',verifier.SUMMARY,flush=True)


if __name__=='__main__':
 main()
