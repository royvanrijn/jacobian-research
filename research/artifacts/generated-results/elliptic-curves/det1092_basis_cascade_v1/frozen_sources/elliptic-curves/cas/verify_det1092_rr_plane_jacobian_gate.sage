#!/usr/bin/env sage-python
"""Independent RR-plane coefficient and fixed-Jacobian proof replay.

Counts Fp2 with integer pairs and norm characters, not Sage field iteration.
Checks extension Frobenius polynomials instead of coefficient exceptions.
"""
import hashlib
import json
from pathlib import Path
from sage.all import QQ, ZZ, GF, PolynomialRing, EllipticCurve, matrix, gcd

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
CERT = ART/'det1092_rr_plane_jacobian_gate_v2.json'
OUT = ART/'det1092_rr_plane_jacobian_gate_replay_v2.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def check():
    d = read(CERT)
    for path, h in d['protocol']['inputs'].items():
        assert sha(ROOT/path) == h
    pp = ART/'curve302_recovered_mw17_parent_v1.json'
    np = ART/'det1092_first_centre_rr_net_v1.json'
    up = ART/'det1092_universal_rr_descent_preflight_v1.json'
    hp = ART/'det1092_rr_net_halving_gate_v1.json'
    bp = ART/'det1092_historical_unlock_obstruction_v1.json'
    parent, net = read(pp), read(np)
    R = PolynomialRing(QQ, 't')
    K = R.fraction_field()
    def dec(v):
        return K(R(v['numerator']))/R(v['denominator'])
    E = EllipticCurve(K, [dec(v) for v in parent['a_invariants']])
    points = [E([dec(v) for v in row]) for row in parent['basis_weierstrass_coordinates']]
    C = -sum((ZZ(n)*P for n,P in zip(net['trace_word'],points)), E(0))
    cx = C[0]+E.b2()/12
    cy = C[1]+(E.a1()*C[0]+E.a3())/2
    aa = -E.c4()/48
    h = R(C[0].denominator().sqrt())
    P = PolynomialRing(QQ, ['X','Y','Z'])
    X,Y,Z = P.gens()
    plane = sum(QQ(row['coefficient'])*X**row['X']*Y**row['Y']*Z**row['Z'] for row in d['plane']['sparse_B'])
    scale = QQ(d['plane']['scale'])
    assert plane.is_homogeneous() and plane.degree()==6
    assert min(i+k for i,j,k in plane.dict())==2
    tangent = sum(c*X**i*Z**k for (i,j,k),c in plane.dict().items() if i+k==2)
    tc = [tangent.monomial_coefficient(v) for v in [X*X,X*Z,Z*Z]]
    assert tc == list(map(QQ,d['plane']['tangent_cone_coefficients_X2_XZ_Z2']))
    assert tc[1]**2-4*tc[0]*tc[2] == QQ(d['plane']['tangent_discriminant']) != 0
    # Compare the entire universal line restriction, not sampled members.
    U = PolynomialRing(QQ, ['T','u','v'])
    T,u,v = U.gens()
    universal = read(up)['universal_genus2']
    previous = sum(QQ(row['coefficient'])*T**row['T']*u**row['u']*v**row['v'] for row in universal['sparse_q'])
    assert scale*plane(T,u+v*T,1) == QQ(universal['scale'])*previous
    # Parent equation checked again with all known denominators cleared.
    def lift(f):
        return U(R(f)(T))
    A,B = [[lift(f) for f in net[key]] for key in ['A','B']]
    f1,f2 = B[1]+u*A[1],B[2]+u*A[2]
    hh,nx,ny = lift(h),lift(cx*h*h),lift(cy*h**3)
    mnum = -f1+lift(E.a1())*f2/2
    raw = hh**4*mnum**4-6*nx*hh**2*mnum*mnum*f2*f2-8*ny*hh*mnum*f2**3-(3*nx*nx+4*lift(aa)*hh**4)*f2**4
    assert raw == hh**10*scale*plane(T,u,1)
    # No vertical content in the affine branch polynomial.
    by_y = [R([plane.monomial_coefficient(X**i*Y**j*Z**(6-i-j)) for i in range(7-j)]) for j in range(5)]
    assert gcd(by_y).degree()==0 and by_y[4].degree()==2
    minimal = R(read(bp)['q_coefficients'])
    ratio = K(scale*by_y[4])/minimal
    assert ratio.numerator().degree()==ratio.denominator().degree()==0
    ratio = QQ(ratio)
    assert ratio.is_square()
    # Existing exact halving branch genus supplies geometric irreducibility
    # and normalization genus9; the node then accounts for all delta=1.
    halves = read(hp)
    assert halves['geometric_halving_cover']['genus']==9
    assert halves['geometric_halving_cover']['degree']==4
    fixed = R([QQ(c) for c in d['fixed_member']['q_coefficients']])
    assert fixed == R(scale*plane(R.gen(),0,1))
    assert fixed.degree()==6 and fixed.gcd(fixed.derivative()).degree()==0
    # This frozen run stopped at its first prime, so there are no earlier
    # exclusions to infer or discard.
    assert len(d['trials'])==1
    trial = d['trials'][0]
    p = int(trial['prime'])
    assert p==d['protocol']['primes'][0]==d['absolute_simplicity_witness_prime']==17
    Fp = PolynomialRing(GF(p),'x')
    fp = Fp(fixed.list())
    assert fp.degree()==6 and fp.discriminant()
    coeff = list(map(int,fp.list()))
    assert coeff==trial['coefficients']
    def chi(a):
        z = pow(int(a)%p,(p-1)//2,p)
        return 0 if z==0 else 1 if z==1 else -1
    def eval1(x):
        y = 0
        for c in reversed(coeff):
            y=(y*x+c)%p
        return y
    n1=p+sum(chi(eval1(x)) for x in range(p))+1+chi(coeff[6])
    nonsquare=next(x for x in range(2,p) if chi(x)==-1)
    def mul(a,b):
        return ((a[0]*b[0]+nonsquare*a[1]*b[1])%p,(a[0]*b[1]+a[1]*b[0])%p)
    def eval2(x):
        y=(0,0)
        for c in reversed(coeff):
            y=mul(y,x)
            y=((y[0]+c)%p,y[1])
        return y
    n2=p*p+2
    for i in range(p):
        for j in range(p):
            a,b=eval2((i,j))
            n2+=chi(a*a-nonsquare*b*b)
    assert [n1,n2]==trial['point_counts_Fp_Fp2']==[30,280]
    s1,s2=p+1-n1,p*p+1-n2
    e2=ZZ((s1*s1-s2)//2)
    frob=R([p*p,-p*s1,e2,-s1,1])
    assert list(map(int,frob.list()))==trial['frobenius_coefficients']
    assert e2%p and frob.is_irreducible()
    comp=matrix(QQ,4,4)
    for i in range(3):
        comp[i+1,i]=1
    for i in range(4):
        comp[i,3]=-frob[i]
    assert R(comp.charpoly())==frob
    extension_polynomials=[]
    for degree in [2,3,4,6]:
        f=R((comp**degree).charpoly())
        assert f.is_irreducible()
        extension_polynomials.append({'extension_degree':degree,'coefficients':list(map(str,f.list()))})
    return {'classification':'verified application and new deduction',
            'status':'PASS_INDEPENDENT_RR_PLANE_AND_ABSOLUTE_JACOBIAN_SIMPLICITY',
            'inputs':{str(path.relative_to(ROOT)):sha(path) for path in [CERT,pp,np,up,hp,bp]},
            'checker_sha256':sha(Path(__file__)),
            'plane_degree':6,'plane_terms':len(plane.dict()),'only_geometric_singularity':'one ordinary node [0:1:0]',
            'minimal_bisection_tangent_square_ratio':str(ratio),
            'curve_mod17':coeff,'counts_F17_F289':[n1,n2],
            'frobenius_polynomial':str(frob),'ordinary':True,
            'extension_polynomials':extension_polynomials,
            'fixed_Jacobian_absolutely_simple':True,'geometric_generic_RR_Jacobian_absolutely_simple':True,
            'point_searches':0,'Selmer_runs':0,'class_group_runs':0,
            'boundary':'No morphism to an elliptic curve from the fixed B member or geometric generic RR member, of any degree. Special RR subfamilies can still have split Jacobians and have not been classified. No genus2 Selmer panel is claimed.'}


if __name__=='__main__':
    d=check()
    text=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if OUT.exists():
        assert OUT.read_text()==text
    else:
        OUT.write_text(text)
    print(d['status'],flush=True)
