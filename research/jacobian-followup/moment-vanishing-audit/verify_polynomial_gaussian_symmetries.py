#!/usr/bin/env python3
"""Exact supplied-inverse certificates; finite controls do not prove the universal theorem."""
from pathlib import Path
import hashlib,json
import sympy as S

class Quotient:
    def __init__(self,epsilon=None,length=1):
        if epsilon is not None and (not isinstance(length,int) or length<2):raise ValueError('nilpotent quotient length must be at least2')
        self.e=epsilon;self.L=length if epsilon is not None else 1
    def red(self,f):
        f=S.expand(f)
        if self.e is None:return f
        return S.expand(sum(c*self.e**k[0] for k,c in S.Poly(f,self.e).terms() if k[0]<self.L))
    def inverse_constant(self,u):
        if self.e is None:
            if u==0:raise ValueError('Jacobian constant is not a unit')
            return 1/u
        u=self.red(u);a=u.subs(self.e,0)
        if a==0:raise ValueError('Jacobian constant is not a unit')
        t=self.red(u/a-1)
        return self.red(sum((-t)**k for k in range(self.L))/a)
    def normalized_log(self,j,x):
        u=self.red(j.subs(dict.fromkeys(x,0)));uinv=self.inverse_constant(u)
        t=self.red(j*uinv-1)
        if self.red(t**self.L)!=0:raise ValueError('Nonconstant Jacobian is not a nilpotent unit')
        ell=self.red(sum(S.Rational((-1)**(k+1),k)*self.red(t**k) for k in range(1,self.L)))
        jinv=self.red(uinv*sum((-t)**k for k in range(self.L)))
        return ell,jinv,u

def compose(f,g,x,ring):return ring.red(f.subs(dict(zip(x,g)),simultaneous=True))
def energy(x,B):return S.expand((S.Matrix(x).T*B.inv()*S.Matrix(x))[0]/2)

def certify_symmetry(phi,psi,x,B,epsilon=None,length=1):
    x=tuple(x);ring=Quotient(epsilon,length)
    if not x or len(set(x))!=len(x) or len(phi)!=len(x) or len(psi)!=len(x):raise ValueError('Matching maps and distinct nonempty coordinates required')
    if B.shape!=(len(x),len(x)) or B!=B.T or B.det()==0:raise ValueError('Fixed symmetric invertible covariance required')
    polyvars=x+((epsilon,) if epsilon is not None else ())
    phi=[S.Poly(f,*polyvars,domain=S.QQ_I).as_expr() for f in phi]
    psi=[S.Poly(f,*polyvars,domain=S.QQ_I).as_expr() for f in psi]
    inv=all(compose(f,psi,x,ring)==q and compose(g,phi,x,ring)==q for f,g,q in zip(phi,psi,x))
    if not inv:raise ValueError('Supplied inverse does not compose to identity in both orders')
    jac=S.Matrix(phi).jacobian(x);j=ring.red(jac.det());ell,jinv,u=ring.normalized_log(j,x)
    en=energy(x,B);defect=ring.red(compose(en,phi,x,ring)-en);res=ring.red(defect-ell)
    preserving=all(ring.red(S.diff(res,q))==0 for q in x)
    # Cofactor/Piola identities and transported drift are independent algebraic controls.
    adj=jac.adjugate();piola=[]
    for i in range(len(x)):
        v=[ring.red(jinv*adj[k,i]) for k in range(len(x))]
        pd=ring.red(sum(S.diff(ring.red(j*v[k]),x[k]) for k in range(len(x))))
        assert pd==0;piola.append(str(pd))
        drift=ring.red(sum(v[k]*S.diff(en,x[k])-S.diff(v[k],x[k]) for k in range(len(x))))
        target=ring.red((B.inv()*S.Matrix(phi))[i])
        assert ring.red(drift-target+sum(v[k]*S.diff(res,x[k]) for k in range(len(x))))==0
    return dict(inverse_verified=True,jacobian=str(j),jacobian_at_origin=str(u),normalized_log_jacobian=str(ell),
        energy_defect=str(defect),residual=str(res),gaussian_preserving=preserving,energy_preserving=defect==0,
        jacobian_constant_in_coordinates=all(ring.red(S.diff(j,q))==0 for q in x),piola_controls=piola)

def wick(f,x,ring):
    # split pair Z,W and optional independent variance-one T
    out=0
    for a,c in S.Poly(ring.red(f),*x).terms():
        if a[0]!=a[1] or any(e%2 for e in a[2:]):continue
        out+=c*S.factorial(a[0])*S.prod(S.factorial2(e-1) if e else 1 for e in a[2:])
    return ring.red(out)

def main():
    Z,W,T,e=S.symbols('Z W T e');x=(Z,W,T);B=S.Matrix([[0,1,0],[1,0,0],[0,0,1]]);q=2*Z*W+T*T;ring=Quotient()
    def invariant_time(f):return [Z,S.expand(W-T*f-Z*f*f/2),S.expand(T+Z*f)]
    nag=invariant_time(q);naginv=invariant_time(-q)
    nc=certify_symmetry(nag,naginv,x,B);assert nc['gaussian_preserving'] and nc['energy_preserving'] and nc['jacobian']=='1'
    # identity and orthogonal reflection are included, with determinant sign unrestricted.
    identity=certify_symmetry(x,x,x,B);reflection=certify_symmetry([-Z,-W,T],[-Z,-W,T],x,B)
    assert identity['gaussian_preserving'] and reflection['gaussian_preserving']
    sign=certify_symmetry([Z,W,-T],[Z,W,-T],x,B);assert sign['gaussian_preserving'] and sign['jacobian']=='-1'
    samples=[S.S.One,Z,W,T,Z*W,T*T,W*W,Z*T,W*T,Z*W*T]
    finite_moments=0
    for f in samples:
        assert wick(compose(f,nag,x,ring),x,ring)==wick(f,x,ring);finite_moments+=1
    # Translation is volume-preserving but not Gaussian-preserving.
    tr=certify_symmetry([Z+1,W,T],[Z-1,W,T],x,B);assert not tr['gaussian_preserving'] and tr['jacobian']=='1'
    # An inverse is required and validated; no constant-Jacobian inversion oracle.
    rejected=False
    try:certify_symmetry(nag,x,x,B)
    except ValueError:rejected=True
    assert rejected
    r2=Quotient(e,2)
    a=[Z-e*Z*Z,W+e*(Z*W-1),T];ainv=[f.subs(e,-e) for f in a]
    ac=certify_symmetry(a,ainv,x,B,e,2);assert ac['gaussian_preserving'] and not ac['energy_preserving']
    assert S.expand(S.sympify(ac['residual']))==0 and S.sympify(ac['jacobian'])==1-e*Z
    for f in samples:
        assert wick(compose(f,a,x,r2),x,r2)==wick(f,x,r2);finite_moments+=1
    bad=[Z,W-e*T*T,T+e*Z*T];badinv=[f.subs(e,-e) for f in bad]
    bc=certify_symmetry(bad,badinv,x,B,e,2);assert not bc['gaussian_preserving'] and bc['energy_preserving']
    badmoment=wick(bad[1],x,r2);assert badmoment==-e
    r3=Quotient(e,3);v=[1-Z*Z-Z*W,-1+Z*W+W*W,0]
    def V(f):return S.expand(sum(a*S.diff(f,b) for a,b in zip(v,x)))
    flow=[r3.red(t+e*V(t)+e*e*V(V(t))/2) for t in x];flowinv=[f.subs(e,-e) for f in flow]
    fc=certify_symmetry(flow,flowinv,x,B,e,3);assert fc['gaussian_preserving'] and S.sympify(fc['residual'])==-e*e
    assert S.sympify(fc['jacobian_at_origin'])==1-e*e
    for f in samples:
        assert wick(compose(f,flow,x,r3),x,r3)==wick(f,x,r3);finite_moments+=1
    # A small exact all-index-formula transcription set, not a degree campaign.
    cutoff_rows=[]
    for j,k,d in [(0,0,1),(0,1,1),(0,1,2),(0,2,1),(1,0,1),(1,0,2),(1,1,1)]:
        f=Z**j*q**k;p=S.expand(W-T*f-Z*f*f/2);m=(2*j+1)*d
        actual=wick(S.expand(W**d*p**m),x,ring)
        expected=(-S.Rational(1,2))**d*S.factorial(m)*2**(2*k*d)*S.rf(S.Rational(2*m+3,2),2*k*d)
        assert S.simplify(actual-expected)==0 and actual!=0
        cutoff_rows.append(dict(j=j,k=k,d=d,m=m,value=str(actual),cutoff=2*j+1))
    # Constant q-degree growth retains cutoff1; sampled pure moments only transcribe the symmetry theorem.
    for m in [1,2,3]:assert wick(S.expand(nag[1]**m),x,ring)==0
    lap=lambda f:S.expand(2*S.diff(f,Z,W)+S.diff(f,T,2))
    lap_output=lap(nag[1]);assert S.expand(lap_output+10*T+14*Z*q)==0 and lap_output!=0
    # Literal marked transport, distinct from same-degree ideal equality.
    p=W;marker=W*T+Z
    marked_transport=[]
    for m in [1,2,3]:
        left=wick(S.expand(marker*nag[1]**m),x,ring)
        right=wick(S.expand(compose(marker,naginv,x,ring)*p**m),x,ring)
        assert left==right;marked_transport.append(str(left))
    receipt=dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),sympy_version=S.__version__,
        nagata_certificate=nc,identity_certificate=identity,reflection_certificate=sign,translation_rejection=tr,
        nilpotent_true_symmetry=ac,quadratic_only_nilpotent_rejection=bc,quadratic_only_first_marked=str(badmoment),
        nonzero_constant_true_symmetry=fc,finite_moment_controls=finite_moments,sharp_cutoff_rows=cutoff_rows,
        nagata_first_Laplacian=str(lap_output),marked_transport_controls=marked_transport,supplied_wrong_inverse_rejected=rejected,
        exact_Q_i_supplied_inverse_classifier_implemented=True,monogenic_nilpotent_quotient_classifier_implemented=True,
        finite_normalized_log_implemented=True,quadratic_invariant_time_constructor_implemented=True,
        general_polynomial_inverse_solver_implemented=False,general_orbit_membership_solver_implemented=False,
        general_quadratic_automorphism_generator_classification_implemented=False,universal_theorem_verified_by_controls=False,
        quadratic_only_criterion_valid_over_all_Q_algebras=False,ordinary_Laplacian_transfer_claimed=False,
        independent_replay=False,external_review=False,formal_verification=False,novelty_claimed=False)
    dest=Path(__file__).with_name('polynomial_gaussian_symmetry_verification.json');dest.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['completed','checker_sha256','finite_moment_controls','quadratic_only_first_marked','nagata_first_Laplacian','marked_transport_controls']}))

if __name__=='__main__':main()
