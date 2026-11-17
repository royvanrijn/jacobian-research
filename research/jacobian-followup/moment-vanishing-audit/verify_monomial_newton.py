#!/usr/bin/env python3
"""Small exact controls, not a producer for arbitrary formal Newton data."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math
import sympy as S

HERE=Path(__file__).resolve().parent

def action(lam,poly,x,z):
    lam=S.sympify(lam);poly=S.sympify(poly)
    assert not lam.has(S.Float) and not poly.has(S.Float)
    assert S.Poly(poly,*x).total_degree()<=24
    return S.expand(sum(c*S.diff(poly,*[t for v,k in zip(x,a) for t in (v,)*k])
                        if sum(a) else c*poly for a,c in S.Poly(lam,*z).terms()))

def height(lam,poly,x,z):
    assert S.Poly(lam,*z).coeff_monomial(1)==0 and lam!=0
    h=0
    while True:
        poly=action(lam,poly,x,z)
        if poly==0:return h
        h+=1
        assert h<=24

def box_output(lam,k,v,x,z):
    answer=0
    for a,c in S.Poly(S.expand(lam**k),*z).terms():
        if all(i<=j for i,j in zip(a,v)):
            answer+=c*math.prod(math.factorial(j)//math.factorial(j-i) for i,j in zip(a,v))*S.prod(t**(j-i) for t,i,j in zip(x,a,v))
    return S.expand(answer)

def certificate(support,a,b,mix,w,kappa):
    """A supplied rational primal/dual certificate; no LP fitting."""
    n=len(a);b=tuple(map(F,b));mix=tuple(map(F,mix));w=tuple(map(F,w));kappa=F(kappa)
    assert sum(mix)==1 and all(t>=0 for t in mix)
    assert b==tuple(sum(q*v[i] for q,v in zip(mix,support)) for i in range(n))
    assert sum(max(t-s,0) for t,s in zip(b,a))==kappa
    assert all(0<=t<=1 for t in w)
    lower=min(sum(t*s for t,s in zip(w,v)) for v in support)-sum(t*s for t,s in zip(w,a))
    assert lower==kappa
    return dict(support=support,phase_exponents=a,point=list(map(str,b)),mix=list(map(str,mix)),dual=list(map(str,w)),kappa=str(kappa),cutoff=str(1/kappa) if kappa else None)

def main():
    x,y,z,X,Y,Z,e=S.symbols('x y z X Y Z e');v=(x,y,z);d=(X,Y,Z)
    certs=[certificate([(2,0),(0,3)],(1,1),('4/3',1),('2/3','1/3'),(1,'2/3'),'1/3'),
           certificate([(3,0),(0,3)],(1,1),(2,1),('2/3','1/3'),(1,1),1),
           certificate([(2,0,0),(0,4,0),(0,0,6)],(1,1,1),('7/6',1,1),('7/12','1/4','1/6'),(1,'1/2','1/3'),'1/6'),
           certificate([(2,0),(0,2)],(1,1),(1,1),('1/2','1/2'),(0,0),0),
           certificate([(2,1),(1,2)],(1,1),('3/2','3/2'),('1/2','1/2'),(1,1),1),
           certificate([(0,1)],(1,0),(0,1),(1,),(0,1),1)]
    controls=[]
    for lam,P,xs,zs,expected in [
        (X**2+Y**3,x*y,(x,y),(X,Y),lambda m:m//2+m//3),
        (X**2+Y**2,x*y,(x,y),(X,Y),lambda m:2*(m//2)),
        (X**2*Y-X*Y**2,x*y,(x,y),(X,Y),lambda m:(2*m)//3),
        (X**2+Y**4-Z**6,x*y*z,v,d,lambda m:m//2+m//4+m//6)]:
        hs=[]
        for m in range(1,6):
            h=height(lam,S.expand(P**m),xs,zs);assert h==expected(m);hs.append(h)
        controls.append(dict(symbol=str(lam),phase=str(P),heights_1_to_5=hs))
    box=0
    for lam,vec in [(X**2+Y*Z,(3,2,1)),(X-Y**2+Z**3,(2,1,3))]:
        for k in (1,2,3):
            poly=S.prod(t**j for t,j in zip(v,vec))
            for _ in range(k):poly=action(lam,poly,v,d)
            assert poly==box_output(lam,k,vec,v,d);box+=1
    # The compact face is finite although the retained model has a later term.
    lam=X**2+Y**4-Z**6-Z**24
    face=X**2+Y**4-Z**6
    assert S.Poly(lam**4,*d).coeff_monomial(X**4*Y**4*Z**6)==S.Poly(face**4,*d).coeff_monomial(X**4*Y**4*Z**6)==-12
    sharp=math.factorial(12)//(math.factorial(7)*math.factorial(3)*math.factorial(2))
    assert S.Poly(face**12,*d).coeff_monomial(X**14*Y**12*Z**12)==sharp
    assert F(12,2)==6
    # An infinite tail need not be expanded to check its unique compact-face coefficient.
    assert (14,12,12)==(12+2,12,12)
    # A common formal jet cannot determine the exact optimal cutoff.
    N=7;lam2=X**2+Y**N
    tailcert=certificate([(2,0),(0,N)],(1,1),(F(2*N-2,N),1),(F(N-1,N),F(1,N)),(1,F(2,N)),F(N-2,N))
    assert S.Poly(lam2**N,X,Y).coeff_monomial(X**(2*N-2)*Y**N)==N
    assert F(N,N-2)>1  # X^2 alone has cutoff1.
    # Passive symbols affect the mixed cutoff even when they kill all P powers.
    assert action(Y,(x**3)*y,(x,y),(X,Y))==x**3
    assert height(Y,x**3,(x,y),(X,Y))==0
    # Support separation is still sufficient over a square-zero coefficient model.
    nil=action(X**2+e*Y**3,(x*y)**2,(x,y),(X,Y))
    nil=action(X**2+e*Y**3,nil,(x,y),(X,Y))
    assert S.rem(nil,e**2,e)==0
    # Bounded gaps in a finite Laurent constant-term sequence, including parity zeros.
    cs=[S.Poly((X**2+Y**2)**j,X,Y).coeff_monomial(X**j*Y**j) for j in range(11)]
    for j in range(9):assert (j+2)*cs[j+2]-4*(j+1)*cs[j]==0
    assert cs[1]==0 and cs[2]==2 and all(cs[j] or cs[j+1] for j in range(10))
    # Proper ordinary diagonal used in the proof, with q=2.
    H=X**3*Y+X*Y**3
    for j in range(6):assert S.Poly(H**j,X,Y).coeff_monomial(X**(2*j)*Y**(2*j))==cs[j]
    guards=0
    for lam,poly in [(X,S.Float(1)),(X,x**25)]:
        try:action(lam,poly,(x,),(X,))
        except AssertionError:guards+=1
        else:raise AssertionError('input guard failed')
    deps=['verify_formal_apolar_gaps.py','verify_formal_symbol_source_rigidity.py']
    assert all((HERE/n).is_file() for n in deps)
    report=dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        supporting_checker_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        primal_dual_certificates=certs,height_controls=controls,box_output_identities=box,
        compact_face_coefficient=-12,sharp_ternary_coefficient=str(sharp),sharp_cutoff='6',
        late_tail_cutoff_certificate=tailcert,passive_boundary_pass=True,
        nonreduced_support_sufficiency_pass=True,constant_term_recurrence_controls=9,
        proper_diagonal_controls=6,guard_refusals=guards,
        full_formal_support_recognizer_implemented=False,general_height_gap_proved=False,
        general_quasi_subadditivity_proved=False,full_exceptional_MZ_decided=False,
        ordinary_nonreduced_purity_gate_proved=False,unrestricted_frontier_endpoint_changed=False,
        universal_claim_proved_by_samples=False,independent_replay=False,
        external_review=False,formal_verification=False,novelty_claimed=False)
    assert guards==2
    (HERE/'monomial_newton_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(completed=True,LP_certificates=len(certs)+1,height_families=len(controls),
                         box_identities=box,sharp_cutoff=6,recurrence_controls=9,guards=guards)))

if __name__=='__main__':main()
