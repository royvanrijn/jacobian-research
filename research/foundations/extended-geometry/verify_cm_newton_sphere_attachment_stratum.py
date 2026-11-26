"""Exact and modular checks for the generic-boundary sphere stratum."""
from sage.all import GF, QQ, PolynomialRing, matrix
from math import comb
from collections import defaultdict

def evaluate_pair(field, b_coeffs, c, e, irreducible=None):
    rr=PolynomialRing(field, "z")
    z=rr.gen()
    PP=(z**2+1)**4
    V=z*rr(b_coeffs)
    W=-1+c*z+e*z*z
    a=[18*PP,17*V,16*W]
    g=[z*PP.derivative()-PP,z*V.derivative()-V,z*W.derivative()-W]
    alpha=a[2]*g[1]-g[2]*a[1]
    beta=a[2]*g[0]-g[2]*a[0]
    numerator=beta**2-alpha*(a[1]*g[0]-g[1]*a[0])
    R, remainder=numerator.quo_rem((z*z+1)**3)
    assert not remainder and R.degree()==16 and R[0]==4
    assert R.leading_coefficient()==-1827*e*b_coeffs[5]**2
    assert W.is_squarefree() and (V*V-4*PP*W).is_squarefree()
    assert V.gcd(z*z+1).is_one()
    assert R.is_squarefree()
    assert R.gcd(z*(z*z+1)*W*alpha*beta).is_one()
    R=R.monic()
    if irreducible:
        assert R.is_irreducible()
    A=rr.quotient(R, "Z")
    x=A.gen()
    evaluate=lambda cs, y: A(cs[0])+y*(A(cs[1])+y*A(cs[2]))
    y=-A(beta)/A(alpha)
    assert not evaluate(a,y) and not evaluate(g,y)
    az=evaluate([f.derivative() for f in a],y)
    gz=evaluate([f.derivative() for f in g],y)
    ah=evaluate([a[1],2*a[2],rr.zero()],y)
    gh=evaluate([g[1],2*g[2],rr.zero()],y)
    assert (az*gh-ah*gz).is_unit()
    level=evaluate([PP,V,W],y)/(x*y**18)
    assert level.is_unit()
    vectors=[]
    l=A.one()
    for k in range(16):
        vectors.append([l.lift()[j] for j in range(16)])
        l*=level
    det=matrix(field, vectors).det()
    assert det
    return det

def generic_seed():
    bb=[1,0,0,0,0,1]
    # The two derivative identities are exact over characteristic zero.
    # A finite-field determinant and unit checks certify the same quantities over Q.
    qq=PolynomialRing(QQ,"z")
    z=qq.gen()
    P=(z*z+1)**4
    V=z*(1+z**5)
    W=-1+z+z*z
    a=[18*P,17*V,16*W]
    g=[z*P.derivative()-P,z*V.derivative()-V,z*W.derivative()-W]
    alpha=a[2]*g[1]-g[2]*a[1]
    beta=a[2]*g[0]-g[2]*a[0]
    R,rem=(beta**2-alpha*(a[1]*g[0]-g[1]*a[0])).quo_rem((z*z+1)**3)
    assert not rem
    assert list(R)==[4,64,635,4992,27418,-25152,55283,-87180,25052,-39976,-14940,16420,-2159,29616,18390,-2320,-1827]
    A=qq.quotient(R.monic(),"Z")
    y=-A(beta)/A(alpha)
    assert all(not (A(cs[0])+y*(A(cs[1])+y*A(cs[2]))) for cs in (a,g))
    det=evaluate_pair(GF(103),bb,1,1,irreducible=True)
    assert det == 33
    print("Generic 16-value seed determinant: 33 modulo 103")


def moment_rows(ring, max_k, odd_coefficients, c, e, correction_count, *, include_a0=True):
    phase = [(2*j-1, 0, ring(comb(4,j))) for j in range(5)]
    phase += [(j, 1, ring(value)) for j, value in enumerate(odd_coefficients) if value]
    phase += [(-1, 2, ring(-1)), (0, 2, ring(c)), (1, 2, ring(e))]
    a0 = [(0, 0, 1), (2, 0, 2), (4, 0, 1)]
    a1 = [(1, 0, 1), (3, 0, 2), (5, 0, 1)]
    markers = ([a0] if include_a0 else []) + [a1] + [[(j, 1, 1)] for j in range(correction_count)]
    powers={(0,0):ring.one()}
    rows=[]
    for k in range(max_k+1):
        rows.append([
          sum(ring(weight)*v/(13-par+18*k-loss)
              for angle,par,weight in col
              for (step,loss),v in powers.items() if step==-angle and loss%2==par)
          for col in markers
        ])
        if k<max_k:
            new=defaultdict(ring.zero)
            for (step,loss),v in powers.items():
                for m,q,w in phase:
                    if step+m<=max_k-k-1:
                        new[(step+m,loss+q)] += v*w
            powers=new
    return rows

def exceptional_line():
    p=1009
    ff=GF(p)
    rr=PolynomialRing(ff,"b")
    b=rr.gen()
    bb=[1,1,0,0,1]
    rows = moment_rows(rr, 7, [b, *bb], 1, 1, 6)
    for n,D in enumerate((matrix(rr, [row[:5] for row in rows[:5]]).det(),matrix(rr, rows).det()),1):
        factors = D.factor()
        degrees = [factor.degree() for factor, _ in factors]
        assert D.degree() == (5, 12)[n-1]
        assert degrees == ([1, 1, 1, 2], [2, 3, 7])[n-1]
        assert D[0]!=0 and D.is_squarefree()
        if n==1:
            F=(155708*b**5+1407*b**4+8880*b**3+(3417+13362+1608)*b**2-(1608+402)*b+201)
            assert D.monic()==F.monic()
        dets=[]
        for v, _ in factors:
            ext=GF(p**v.degree(), name="s", modulus=v) if v.degree()>1 else ff
            x=ext.gen() if v.degree()>1 else v.roots()[0][0]
            assert D(x)==0
            det=evaluate_pair(ext,[x]+bb,1,1)
            dets.append((v.degree(),repr(det)))
        print("all exceptional line specializations good modulo",p,dets,flush=True)


def blind_branch(field):
    rr=PolynomialRing(field,"b")
    K=rr.fraction_field()
    b=K.gen()
    coeff=[b,K.one(),K.one(),K.zero(),K.zero(),K.zero()]
    rows=moment_rows(K, 3, coeff, 1, 1, 3, include_a0=False)
    C=matrix(K,[r[1:] for r in rows[1:]]).inverse()*matrix(K,[[r[0]] for r in rows[1:]])
    vector=[K.one()]+[-C[i,0] for i in range(3)]
    inner=lambda row: sum(x*y for x,y in zip(row,vector))
    pivots=[]
    for j in range(3,6):
        n=j+1
        zero=inner(moment_rows(K, n, coeff, 1, 1, 3, include_a0=False)[-1])
        coeff[j]=K.one()
        one=inner(moment_rows(K, n, coeff, 1, 1, 3, include_a0=False)[-1])
        pivot=one-zero
        assert pivot
        pivots.append(pivot)
        coeff[j]=-zero/pivot
    res=[inner(row) for row in moment_rows(K, 7, coeff, 1, 1, 3, include_a0=False)]
    assert not any(res[:7])
    return coeff,vector,pivots,res[7]

def first_eight_coherent_blindness():
    coeff,vector,pivots,rho=blind_branch(QQ)
    row_eight = moment_rows(rho.parent(), 8, coeff, 1, 1, 3, include_a0=False)[8]
    response_eight = sum(x*y for x, y in zip(row_eight, vector, strict=True))
    assert response_eight.denominator() == rho.denominator()
    f, g = rho.numerator(), response_eight.numerator()
    gcd, left, right = f.xgcd(g)
    assert gcd == 1 and left*f + right*g == 1
    p=1009
    field=GF(p)
    ring=PolynomialRing(field,"b")

    def reduce_poly(poly):
        assert all(c.denominator()%p for c in poly)
        return ring([field(c) for c in poly])

    reduced_coeff,reduced_vector,reduced_pivots,reduced_rho=blind_branch(field)
    assert rho.denominator()==rho.parent().gen()**6
    assert rho.numerator().degree()==12
    assert reduce_poly(rho.denominator())==reduced_rho.denominator()
    assert reduce_poly(rho.numerator())==reduced_rho.numerator()
    for original,reduced in zip(coeff+vector+pivots,reduced_coeff+reduced_vector+reduced_pivots):
        assert reduce_poly(original.numerator())==reduced.numerator()
        assert reduce_poly(original.denominator())==reduced.denominator()

    poly=reduce_poly(rho.numerator().monic())
    assert list(map(int,poly))==[109,847,249,565,381,59,859,657,853,49,620,0,1]
    b=field(519)
    assert not poly(b) and poly.derivative()(b)
    assert all(value.denominator()(b) for value in reduced_coeff+reduced_vector+reduced_pivots)
    lower=[value.numerator()(b)/value.denominator()(b) for value in reduced_coeff]
    inputs=[value.numerator()(b)/value.denominator()(b) for value in reduced_vector]
    measured=[sum(x*y for x,y in zip(row,inputs)) for row in moment_rows(field, 8, lower, 1, 1, 3, include_a0=False)]
    assert not any(measured[:8])
    assert measured[8] == 532
    print("Rows 0 through 8 give an exact unit certificate at b1=b2=c=e=1.")
    determinant=evaluate_pair(field,lower,1,1)
    assert determinant==869
    print("First-eight coherent blind branch at b=519 modulo 1009:",
          f"B={list(map(int, lower))}, (a1,C)={list(map(int,inputs))},",
          f"critical-value basis determinant {determinant}.")

if __name__ == "__main__":
    generic_seed()
    exceptional_line()
    first_eight_coherent_blindness()
