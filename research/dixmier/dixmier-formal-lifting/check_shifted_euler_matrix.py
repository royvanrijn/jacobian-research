#!/usr/bin/env python3
"""Exact skew-Laurent matrix controls for the shifted Euler companion model."""
from __future__ import annotations
import argparse, hashlib, json, platform, tempfile
from pathlib import Path
import sympy as s
HERE=Path(__file__).resolve().parent
J,r,h,beta=s.symbols('J r h beta',nonzero=True)
eta=4*h/r**2
c=4/r**2
n=5

def simp(x):return s.cancel(x)
def clean(x):return {m:simp(a) for m,a in x.items() if simp(a)!=0}
def add(x,y,fac=1):return clean({m:x.get(m,0)+fac*y.get(m,0) for m in set(x)|set(y)})
def scale(x,a):return clean({m:a*b for m,b in x.items()})
def mul(x,y):
 out={}
 for m,a in x.items():
  for q,b in y.items():
   out[m+q]=out.get(m+q,0)+a.subs(J,J-5*q*eta)*b
 return clean(out)
def invmono(x):
 assert len(x)==1
 m,a=next(iter(x.items()))
 return {-m:simp(1/a.subs(J,J+5*m*eta))}
def scalar(a):return clean({0:s.sympify(a)})
ONE=scalar(1);AA={1:s.Integer(1)};AINV={-1:s.Integer(1)}
def madd(x,y,fac=1):return {ij:v for ij in set(x)|set(y) if (v:=add(x.get(ij,{}),y.get(ij,{}),fac))}
def mscale(x,a):return {ij:v for ij,z in x.items() if (v:=scale(z,a))}
def mmul(x,y):
 out={}
 for (i,j),a in x.items():
  for (k,l),b in y.items():
   if j==k:out[i,l]=add(out.get((i,l),{}),mul(a,b))
 return {ij:v for ij,v in out.items() if v}
def mpow(x,k):
 out=diag([ONE]*n)
 for _ in range(k):out=mmul(out,x)
 return out
def diag(xs):return {(i,i):x for i,x in enumerate(xs) if x}
def mcomm(x,y):return madd(mmul(x,y),mmul(y,x),-1)
def expect(x,y):assert madd(x,y,-1)=={},madd(x,y,-1)
def ds(x):
 # [D_s,X], with [D_s,A]=-4hA/r, [D_s,J]=-2hJ/r, [D_s,r]=h.
 return clean({m:h*(s.diff(a,r)-2*J*s.diff(a,J)/r-4*m*a/r) for m,a in x.items()})
def dsmat(x):return {ij:v for ij,z in x.items() if (v:=ds(z))}
def encode(x):return {str(ij):{str(m):str(a) for m,a in z.items()} for ij,z in sorted(x.items())}

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path);args=ap.parse_args()
 I=diag([ONE]*n)
 H=diag([scalar(J-i*eta) for i in range(n)])
 P={(i+1,i):ONE for i in range(4)};P[0,4]=AA
 Pi={(i,i+1):ONE for i in range(4)};Pi[4,0]=AINV
 expect(mmul(P,Pi),I);expect(mmul(Pi,P),I)
 expect(mpow(P,5),diag([AA]*n))
 expect(mcomm(H,P),mscale(P,-eta))
 Q=mmul(H,Pi)
 expect(mmul(Q,P),H)
 expect(mcomm(P,Q),mscale(I,eta))
 expect(mmul(P,Q),madd(H,mscale(I,eta)))
 # A scalar amplification is different; the proposed step fails already on a unit edge.
 raw=diag([scalar(J)]*n)
 assert madd(mcomm(raw,P),mscale(P,eta))!={}
 badH=diag([scalar(J+i*eta) for i in range(n)])
 assert madd(mcomm(badH,P),mscale(P,eta))!={}

 # The changed diagonal cubic/quintic pair, with A=U*G(J).
 G=J*(J-beta*eta)
 U={1:simp(1/G)}
 Umat=diag([{1:simp(1/G.subs(J,J-i*eta))} for i in range(n)])
 Umi=diag([invmono(Umat[i,i]) for i in range(n)])
 Zmat=mscale(mmul(Umi,H),1/(5*c))
 expect(mcomm(Umat,Zmat),mscale(I,h))
 expect(mscale(mmul(Umat,Zmat),5*c),H)
 assert Umat[0,0]==U
 assert madd(Umat,diag([U]*n),-1)!={}
 expect(mmul(Umat,diag([scalar(G.subs(J,J-i*eta)) for i in range(n)])),mpow(P,5))

 # All remaining rational Long commutators, using a shifted auxiliary derivation.
 N=diag([scalar(i) for i in range(n)])
 kernel=mcomm(N,P)
 assert kernel!={}
 linearized={}
 for i in range(5):linearized=madd(linearized,mmul(mmul(mpow(P,i),kernel),mpow(P,4-i)))
 expect(linearized,{})
 correction=mscale(N,-4*h/(5*r))
 def dhat(x):return madd(dsmat(x),mcomm(correction,x))
 expect(dhat(P),mscale(P,-4*h/(5*r)))
 expect(dhat(H),mscale(H,-2*h/r))
 expect(dhat(Q),mscale(Q,-6*h/(5*r)))
 S=madd(mscale(P,r),mscale(I,4/(3*r)))
 T=madd(madd(mscale(Q,r/4),mscale(S,-1/(6*r))),mscale(I,4/(27*r**2)))
 expect(mcomm(S,T),mscale(I,h))
 expect(dhat(S),mscale(madd(mscale(S,1/(5*r)),mscale(I,-8/(5*r**2))),h))
 expect(dhat(T),mscale(madd(mscale(T,-1/(5*r)),mscale(S,1/(10*r**2))),h))
 F=madd(madd(mscale(mmul(S,T),1/(5*r)),mscale(T,-8/(5*r**2))),mscale(mmul(S,S),-1/(20*r**2)))
 expect(madd(dhat(S),mcomm(F,S)),{})
 expect(madd(dhat(T),mcomm(F,T)),{})
 assert madd(dsmat(P),mscale(P,4*h/(5*r)))!={}
 # [D_s+correction+F,rI]=hI, while all entries in P,Q,S,T,F commute with r.
 assert ds(scalar(r))==scalar(h)

 # An explicit rational-J basis change for beta=1 makes finite kappa differential matrices.
 G1=J*(J-eta)
 Bvals=[1,J,G1,G1,G1]
 B=diag([scalar(x) for x in Bvals]);Bi=diag([scalar(s.Integer(1)/x) for x in Bvals])
 Pg=mmul(mmul(B,P),Bi)
 Qg=mmul(mmul(B,Q),Bi)
 U1={1:simp(1/G1)};U1i=invmono(U1)
 targetP={(1,0):scalar(J),(2,1):scalar(J-eta),(3,2):ONE,(4,3):ONE,(0,4):U1}
 targetQ={(0,1):ONE,(1,2):ONE,(2,3):scalar(J-2*eta),(3,4):scalar(J-3*eta),(4,0):mul(scalar(J-4*eta),U1i)}
 expect(Pg,targetP);expect(Qg,targetQ)
 expect(mcomm(Pg,Qg),mscale(I,eta))
 expect(mmul(mmul(B,H),Bi),H)
 correctiong=madd(mmul(mmul(B,correction),Bi),mmul(dsmat(B),Bi),-1)
 expect(correctiong,diag([scalar(h*z/r) for z in [0,s.Rational(6,5),s.Rational(12,5),s.Rational(8,5),s.Rational(4,5)]]))
 Sg=mmul(mmul(B,S),Bi);Tg=mmul(mmul(B,T),Bi);Fg=mmul(mmul(B,F),Bi)
 for Xg in (Sg,Tg):
  expect(madd(madd(dsmat(Xg),mcomm(correctiong,Xg)),mcomm(Fg,Xg)),{})
 # Direct check in k(kappa)[K; -h*d/dkappa], U rational and J=5c U Z.
 x=s.symbols('kappa',nonzero=True)
 def oclean(z):return {i:simp(a) for i,a in z.items() if simp(a)!=0}
 def om(xo,yo):
  out={}
  for i,a in xo.items():
   for j,b in yo.items():
    for k in range(i+1):out[i+j-k]=out.get(i+j-k,0)+a*s.binomial(i,k)*(-h)**k*s.diff(b,x,k)
  return oclean(out)
 Ux=(x-3)**3/(x-2)**2
 Rx=-(x-3)*(x-2)/x**2
 Jx={1:simp(-5*c*x*Rx)} # H=-5*kappa*K.
 wrapped=om({1:Jx[1],0:-4*eta},{0:1/Ux})
 left=om({0:Ux},wrapped);right=om(wrapped,{0:Ux})
 assert oclean({i:left.get(i,0)-right.get(i,0) for i in set(left)|set(right)})=={0:5*eta}

 # Classical generic branch: a rank-one projector exists at h=0, without a quantum descent claim.
 ell,j0=s.symbols('ell j0',nonzero=True)
 P0=s.zeros(n)
 for i in range(4):P0[i+1,i]=1
 P0[0,4]=ell**5
 X=P0/ell
 E=sum((X**i for i in range(n)),s.zeros(n))/n
 assert s.simplify(E*E-E)==s.zeros(n)
 assert s.trace(E)==1
 assert s.simplify(P0*E-ell*E)==s.zeros(n)
 S0=r*P0+4*s.eye(n)/(3*r)
 assert s.trace(P0)==0 and s.trace(S0)==20/(3*r)
 xx,qq=s.symbols('X Q')
 rglobal=xx*(2-3*xx*qq)
 assert s.diff(rglobal,xx).subs(xx,0)==2
 assert s.diff(rglobal,qq)==-3*xx**2
 # A singular integral gauge can change residue Jordan form but preserves charpoly.
 ts=s.symbols('t')
 singular=s.Matrix([[0,h],[0,0]])
 gauged=s.diag(1/h,1)*singular*s.diag(h,1)
 assert singular.subs(h,0).rank()==0 and gauged.subs(h,0).rank()==1
 assert singular.subs(h,0).charpoly(ts).as_expr()==gauged.subs(h,0).charpoly(ts).as_expr()==ts**2
 # Quantum matrix trace is not cyclic: this integral unit gauge changes it by -h.
 def oa(a,b,fac=1):return oclean({i:a.get(i,0)+fac*b.get(i,0) for i in set(a)|set(b)})
 def omat(a,b):
  out={}
  for (i,j),aa in a.items():
   for (k,l),bb in b.items():
    if j==k:out[i,l]=oa(out.get((i,l),{}),om(aa,bb))
  return out
 unit={0:s.Integer(1)};pos={0:x};mom={1:-s.Integer(1)}
 bunit={(0,0):unit,(1,1):unit,(0,1):pos}
 buniti={(0,0):unit,(1,1):unit,(0,1):{0:-x}}
 witness={(0,1):pos,(1,0):mom}
 witnessg=omat(omat(bunit,witness),buniti)
 assert oa(witnessg.get((0,0),{}),witnessg.get((1,1),{}))=={0:-h}
 files=['SHIFTED_EULER_MATRIX_MODEL.md','check_shifted_euler_matrix.py',
        'QUINTIC_ORDERING_LAURENT_OBSTRUCTION.md','CLASSICAL_QUINTIC_COVER_SQUARE.md',
        'EXACT_POSITION_CENTRALIZER_AND_COMPLETION.md',
        'VALUED_MATRIX_LATTICE_AND_SHIFTED_ADMISSION.md','MATRIX_EULER_SECTORS_AND_SCALAR_RESIDUES.md']
 receipt={'claim':'DIX-SHIFTED-EULER-MATRIX-MODEL','status':'PASS',
  'environment':{'python':platform.python_version(),'sympy':s.__version__},
  'normal_form':'A^m*f(J); f(J)*A^q=A^q*f(J-5*q*eta), eta=4h/r^2',
  'companion_fifth_power_and_inverse':True,'shifted_Euler_relation':True,
  'matrix_Weyl_pair':True,'matrix_fifth_power_linearization_has_nonzero_kernel':True,'scalar_amplification_and_wrong_shift_rejected':True,
  'changed_diagonal_cubic_pair_canonical':True,'original_diagonal_transport_differs':True,
  'shifted_auxiliary_derivation':True,'all_full_rational_Long_commutators':True,
  'beta1_diagonal_gauge':{'B':'diag(1,J,J(J-eta),J(J-eta),J(J-eta))',
    'forward_edges':'J,J-eta,1,1,U','reverse_edges':'1,1,J-2eta,J-3eta,(J-4eta)U^-1'},
  'finite_rational_kappa_chart_control':True,'gauged_full_tuple_connection_without_J_denominators':True,'classical_rank_one_branch_projector':True,'classical_S_trace':'20/(3r)','original_r_divisor_smooth_components_control':True,'singular_gauge_changes_Jordan_not_charpoly':True,'quantum_matrix_trace_cyclicity_counterexample':'trace changes from0 to-h under integral unit gauge',
  'hashes':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files},
  'coverage_boundary':'Exact noncommutative companion, changed transport, all rational matrix Long commutators, one finite kappa-chart gauge and classical branch/trace controls. Irreducibility, universal all-gauge/coherent-lattice admission, and scalar-residue/constant-sector theorems use written proofs; finite controls do not prove their universal scope. No no-go theorem for unrelated matrix representations, globally admitted carrier, scalar realization or bounded Long lift.',
  'independent_replay':False,'formal_verification':False,'external_review':False}
 out=args.output or Path(tempfile.mkdtemp(prefix='shifted-euler-'))/'shifted_euler_matrix_controls.json'
 out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps({'claim':receipt['claim'],'status':'PASS','receipt':str(out)},indent=2))
if __name__=='__main__':main()
