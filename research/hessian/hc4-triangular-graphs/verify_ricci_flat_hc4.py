"""Exact Ricci-flat subclass controls; universal claims use the written proof."""
from pathlib import Path

import tempfile as _check_tempfile

_CHECK_OUTPUT_DIRECTORY = None

def _check_output(original):
    """Write new runs to a temporary directory; preserve deposited receipts."""
    global _CHECK_OUTPUT_DIRECTORY
    if _CHECK_OUTPUT_DIRECTORY is None:
        _CHECK_OUTPUT_DIRECTORY = Path(_check_tempfile.mkdtemp(prefix="math-check-"))
        print("Verification output:", _CHECK_OUTPUT_DIRECTORY)
    return _CHECK_OUTPUT_DIRECTORY / Path(original).name


import hashlib,json,time
import sympy as s
R=Path(__file__).resolve().parents[2];Hdir=Path(__file__).resolve().parent
receipt=Hdir/'receipt-ricci-flat-hc4.json';cp=Hdir/'sources/ricci-flat-hc4-context.json'
c=json.loads(cp.read_text());start=time.perf_counter();checks=[];alg=0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(name,v,a=False):
 global alg
 assert bool(v),name
 checks.append(name);alg+=int(a)
def zero(M):return all(s.cancel(q)==0 for q in M)
for p,d in c['retained_proof_and_source_sha256'].items():check('sealed '+p,sha(R/p)==d)
# Generic trace form and Hessian cotangent block, no positivity assumptions.
h,e,f=s.diag(1,-1),s.Matrix([[0,1],[0,0]]),s.Matrix([[0,0],[1,0]])
basis=[h,e,f];gram=s.Matrix(3,3,lambda i,j:s.trace(basis[i]*basis[j]))
check('sl2 trace form matrix',gram==s.Matrix([[2,0,0],[0,0,1],[0,1,0]]),True)
check('sl2 trace form nondegenerate',gram.det()==-2,True)
a0,b0,c0,d0=s.symbols('a0 b0 c0 d0');A=s.Matrix([[a0,b0],[c0,d0]]);Ai=A.inv()
p,q,r0,pp,qq,rr=s.symbols('p q r0 pp qq rr');L=s.Matrix([[p,q],[r0,-p]]);M=s.Matrix([[pp,qq],[rr,-pp]])
c11,c12,c22,di11,di12,di22,dj11,dj12,dj22=s.symbols('c11 c12 c22 di11 di12 di22 dj11 dj12 dj22')
C=s.Matrix([[c11,c12],[c12,c22]]);Di=s.Matrix([[di11,di12],[di12,di22]]);Dj=s.Matrix([[dj11,dj12],[dj12,dj22]])
Hi=s.zeros(2).row_join(Ai).col_join(Ai.T.row_join(-Ai.T*C*Ai))
Hiz=Di.row_join((A*L).T).col_join((A*L).row_join(s.zeros(2)))
Hjz=Dj.row_join((A*M).T).col_join((A*M).row_join(s.zeros(2)))
Ci=Hi*Hiz;Cj=Hi*Hjz;Cw=Hi*Di.row_join(s.zeros(2)).col_join(s.zeros(2,4))
check('arbitrary cotangent base trace identity',s.cancel(s.trace(Ci*Cj)-2*s.trace(L*M))==0,True)
check('arbitrary cotangent mixed trace zero',s.cancel(s.trace(Ci*Cw))==0,True)
check('arbitrary cotangent vertical trace zero',s.cancel(s.trace(Cw*Cw))==0,True)
check('traceless square is scalar',zero(L*L-s.trace(L*L)*s.eye(2)/2),True)
# General single-shear inverse and determinant; f has arbitrary degree.
x,y,aa,bb=s.symbols('x y aa bb');z=s.Matrix([x,y]);ell=s.Matrix([[aa,bb]]);vv=s.Matrix([bb,-aa]);ff=s.Function('F')
sh=z+vv*ff((ell*z)[0]);J=sh.jacobian(z)
check('shear transverse direction annihilated',ell.dot(vv)==0,True)
check('all-function shear determinant',s.cancel(J.det())==1,True)
check('shear invariant linear form',zero(ell*sh-ell*z),True)
inverse=sh-vv*ff(s.expand((ell*sh)[0]));check('all-function shear inverse',zero(inverse-z),True)
entries=s.Matrix(list(J));dJ=entries.jacobian(z)
check('shear Jacobian image differential rank one',all(s.cancel(dJ.extract([i,j],[0,1]).det())==0 for i in range(4) for j in range(i+1,4)),True)
check('non-Keller rank-one control really nonunit',s.Matrix([x*x,y]).jacobian(z).det()==2*x,True)
# Direct coordinate curvature/Ricci, independently of the trace formula.
rs,ss=s.symbols('rs ss');coords=[x,y,rs,ss]
def tensors(potential,K):
 Hess=s.hessian(potential,coords);Ab=K.jacobian([x,y]);Abinv=Ab.inv();Cb=Hess[:2,:2]
 inv=s.zeros(2).row_join(Abinv).col_join(Abinv.T.row_join(-Abinv.T*Cb*Abinv))
 Cs=[(inv*Hess.diff(v)).applyfunc(s.expand) for v in coords]
 G=[cc/2 for cc in Cs]
 Rc={(i,j):(G[j].diff(coords[i])-G[i].diff(coords[j])+G[i]*G[j]-G[j]*G[i]).applyfunc(s.expand) for i in range(4) for j in range(4)}
 Ric=s.Matrix(4,4,lambda i,j:s.expand(sum(Rc[a,j][a,i] for a in range(4))))
 return Hess,inv,Cs,Rc,Ric
u=x+y*y;K1=s.Matrix([u,y+u*u]);H1,I1,C1,R1,Ric1=tensors(rs*K1[0]+ss*K1[1],K1)
check('two-shear full coordinate Ricci equals trace formula',zero(Ric1-s.Matrix(4,4,lambda i,j:s.trace(C1[i]*C1[j])/4)),True)
check('two-shear Ricci origin nonzero',Ric1.subs({x:0,y:0,rs:0,ss:0})==s.Matrix([[0,2,0,0],[2,0,0,0],[0,0,0,0],[0,0,0,0]]),True)
check('two-shear scalar curvature identically zero',s.cancel(s.trace(I1*Ric1))==0,True)
check('two-shear unit Hessian',s.expand(H1.det())==1,True)
K2=s.Matrix([x,y+x*x]);h2=y**3+rs*x+ss*(y+x*x);H2,I2,C2,R2,Ric2=tensors(h2,K2)
check('curved Ricci-flat control full coordinate Ricci zero',zero(Ric2),True)
check('curved Ricci-flat control trace formula zero',all(s.cancel(s.trace(C2[i]*C2[j]))==0 for i in range(4) for j in range(4)),True)
check('curved Ricci-flat nonzero curvature operator',R2[0,1][2,1].subs({x:0,y:0,rs:0,ss:0})==-3,True)
check('curved Ricci-flat unit Hessian',s.expand(H2.det())==1,True)
# The full gradient inverse for the Ricci-flat hostile control.
X,Y,U,V=s.symbols('X Y U V');iz=s.Matrix([U,V-U*U]);iw=s.Matrix([X,Y])-s.Matrix([0,3*iz[1]**2]);iw=K2.jacobian([x,y]).subs({x:iz[0],y:iz[1]},simultaneous=True).inv().T*iw
grad=s.Matrix([s.diff(h2,v) for v in coords]);ip=list(iz)+list(iw)
check('curved Ricci-flat explicit full gradient inverse',zero(grad.subs(dict(zip(coords,ip)),simultaneous=True)-s.Matrix([X,Y,U,V])),True)
check('no new import or primary cache',c['external_imports']==[] and c['new_primary_caches']==[])
out=dict(outcome='PASS',check_count=len(checks),algebraic_controls=alg,elapsed_seconds=time.perf_counter()-start,checks=checks,proof_sha256=sha(Hdir/'RICCI_FLAT_HC4_AND_ONE_DIMENSIONAL_JACOBIAN_IMAGES.md'),checker_sha256=sha(Path(__file__)),context_sha256=sha(cp),registered_mathematical_checkers_replayed=[],independent_replay=False,formal_verification=False,external_review=False)
_check_output(receipt).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
