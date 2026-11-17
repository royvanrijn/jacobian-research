#!/usr/bin/env python3
"""Exact identities for the interval-cubic and constant-profile propositions.

Requires SymPy 1.14.0. No certificate file, search, floating point arithmetic,
Singular, or generated Groebner basis is used. The four-variable moment
polynomials are derived from the displayed multinomial sphere formula.
"""
from functools import lru_cache
from math import factorial,gcd,lcm
import sympy as S
from sympy.polys.domains import QQ
from sympy.polys.rings import ring

if not __debug__:
    raise RuntimeError('Run without -O: this verifier requires assertions')

assert S.__version__=='1.14.0'


def average(poly,t):
    return S.expand(sum(c/S.Integer(k+1) for (k,),c in S.Poly(S.expand(poly),t).terms() if k%2==0))


def interval_cubic():
    t,h,v,c,z=S.symbols('t h v c z')
    u=h*t**3+v*(t*t-S.Rational(1,3))+c*t
    M2,M3,M4=[average(u**m,t) for m in (2,3,4)]
    assert S.expand(M2-(105*c*c+126*c*h+45*h*h+28*v*v)/315)==0
    assert S.expand(S.rem(M3,M2,v)-8*v*(56*c*c+105*c*h+45*h*h)/2205)==0
    R=28028*c**4+84084*c**3*h+93275*c*c*h*h+46410*c*h**3+8955*h**4
    assert S.expand(S.rem(M4,M2,v)+16*R/735735)==0
    assert S.expand(S.rem(M3.subs(h,0),M2.subs(h,0),c)+S.Rational(256,4725)*v**3)==0
    R=R.subs({h:1,c:z});P=56*z*z+105*z+45;B=105*z*z+126*z+45
    assert S.rem(R,P,z)==-S.Rational(315,128)*(299*z+159)
    assert S.rem(R,B,z)==-S.Rational(84,125)*(117*z-110)
    assert P.subs(z,-S.Rational(159,299))==S.Rational(446976,89401)
    assert B.subs(z,S.Rational(110,117))==S.Rational(1169375,4563)


@lru_cache(None)
def df(n):
    return 1 if n<=0 else n*df(n-2)


def det_mod(matrix,prime):
    matrix=[[int(c)%prime for c in row] for row in matrix]
    determinant=1
    for column in range(len(matrix)):
        pivot=next((i for i in range(column,len(matrix)) if matrix[i][column]),None)
        if pivot is None:return 0
        if pivot!=column:
            matrix[column],matrix[pivot]=matrix[pivot],matrix[column]
            determinant=-determinant
        entry=matrix[column][column];determinant=determinant*entry%prime
        for row in range(column+1,len(matrix)):
            factor=matrix[row][column]*pow(entry,-1,prime)%prime
            matrix[row]=[(a-factor*b)%prime for a,b in zip(matrix[row],matrix[column])]
    return determinant%prime


def constant_profile():
    x,y,t,a,c,k=S.symbols('x y t a c k')
    v0,v1,d,e=S.symbols('v0 v1 d e');rr=S.symbols('r0:4')
    original=a*y+(v0+v1*x)*t*t+(c+d*x+e*x*x)*t+sum(rr[j]*x**j for j in range(4))
    rotation={x:x,t:t+k*x,y:y-k*t-k*k*x/2}
    assert S.expand((t*t+2*x*y).subs(rotation,simultaneous=True)-(t*t+2*x*y))==0
    transformed=S.Poly(S.expand(original.subs(rotation,simultaneous=True)),x,y,t)
    assert transformed.coeff_monomial(y)==a
    assert transformed.coeff_monomial(t)==c-a*k
    assert transformed.coeff_monomial(t*t)==v0 and transformed.coeff_monomial(x*t*t)==v1
    assert all(sum(powers)<=3 and powers[1]<=1 for powers,coefficient in transformed.terms())
    R,v0,v1,d,e,r0,r1,r2,r3=ring('v0,v1,d,e,r0,r1,r2,r3',QQ)
    variables=(v0,v1,d,e,r0,r1,r2,r3)
    weights=(1,2,2,3,1,2,3,4)
    phase=(0,1,1,2,0,1,2,3)
    tpower=(2,2,1,1,0,0,0,0)
    def vectors(n,index=0,prefix=()):
        if index==8:
            if n==0:yield prefix
            return
        for count in range(n//weights[index]+1):
            yield from vectors(n-count*weights[index],index+1,prefix+(count,))
    def moment(n):
        out={}
        for powers in vectors(n):
            k=sum(a*b for a,b in zip(powers,phase));j=sum(a*b for a,b in zip(powers,tpower))
            if j%2:continue
            denominator=df(2*k+j+1)
            for count in powers:denominator*=factorial(count)
            out[powers]=QQ(factorial(n)*df(j-1),denominator)
        return R.from_dict(out)
    assert moment(0)==1
    assert moment(1)==r0+v0/3
    assert moment(2)==r0*r0+2*r0*v0/3+v0*v0/5+2*r1/3+2*v1/15
    rr=(r0,r1,r2,r3)
    solutions=(-v0/3,-(2*v0*v0+3*v1)/15,-4*v0*(31*v0*v0+54*v1)/945,-(14175*d*d+4624*v0**4+10152*v0*v0*v1+3024*v1*v1)/28350)
    coefficients=(QQ(1),QQ(2,3),QQ(2,5),QQ(8,35))
    for n in range(1,5):
        value=moment(n)
        for variable,solution in zip(rr[:n-1],solutions[:n-1]):value=value.compose(variable,solution)
        assert value==coefficients[n-1]*(rr[n-1]-solutions[n-1])
    reduced=[];N=[]
    for n in range(5,9):
        value=moment(n)
        for variable,solution in zip(rr,solutions):value=value.compose(variable,solution)
        reduced.append(value)
        denominator=lcm(*(int(c.denominator) for c in value.values()))
        integers={powers:int(c*denominator) for powers,c in value.items()}
        content=gcd(*(abs(c) for c in integers.values()))
        if integers[max(integers)]<0:content=-content
        N.append(R.from_dict({powers:QQ(c//content) for powers,c in integers.items()}))
    assert reduced[3].compose(v0,R.zero)==-QQ(65536,41023125)*v1**4
    assert reduced[1].compose(v0,R.zero).compose(v1,R.zero)==QQ(8,231)*e*e
    BR,b,D,E,F=ring('b,D,E,F',QQ)
    rows=[]
    for n,polynomial in enumerate(N,5):
        row=BR.zero
        for powers,c in polynomial.items():
            i,j,h,k=powers[:4]
            assert not any(powers[4:]) and i+2*j+2*h+3*k==n
            row+=c*b**j*{(0,0):BR.one,(2,0):D,(1,1):E,(0,2):F}[(h,k)]
        rows.append(S.expand(row.as_expr()))
    bs,Ds,Es,Fs=S.symbols('b D E F')
    matrix=S.Matrix([[row.coeff(x) for x in (Ds,Es,Fs)]+[row.subs({Ds:0,Es:0,Fs:0})] for row in rows])
    expected=S.Matrix([
      [155925,155925,0,35504+94992*bs+52704*bs**2],
      [77395500-116093250*bs,77395500,-58046625,30907712+63343584*bs+3500496*bs**2-12573792*bs**3],
      [226800-297675*bs,226800+127575*bs,-212625,53504+17328*bs-211968*bs**2-108864*bs**3],
      [4819500-216877500*bs,4819500-65063250*bs,-75907125,7986496-27544704*bs-90351504*bs**2-16632864*bs**3+9144576*bs**4]])
    assert matrix==expected
    assert matrix[:,0]-matrix[:,1]-2*bs*matrix[:,2]==S.zeros(4,1)
    f=2485161*bs**3+6647562*bs**2+4908168*bs+1097984
    g=10216206*bs**4+24356619*bs**3-25372278*bs**2-41158800*bs-11711872
    assert S.Poly(matrix[[0,1,2],[0,2,3]].det(),bs).primitive()[1].as_expr()==f
    assert S.Poly(matrix[[0,1,3],[0,2,3]].det(),bs).primitive()[1].as_expr()==-g
    fp=S.Poly(f,bs,modulus=19);gp=S.Poly(g,bs,modulus=19)
    assert fp.degree()==3 and gp.degree()==4
    assert fp.as_expr()==-bs**3-6*bs**2-7*bs-7
    assert gp.as_expr()==bs**4+6*bs**3-bs**2+7*bs-6
    quadratic=gp.rem(fp);linear=fp.rem(quadratic);constant=quadratic.rem(linear)
    assert quadratic.as_expr()==-8*bs**2-6
    assert linear.as_expr()==8*bs+7 and constant.as_expr()==-5
    fc=list(map(int,S.Poly(f,bs).all_coeffs()));gc=list(map(int,S.Poly(g,bs).all_coeffs()))
    sylvester=[]
    for shift in range(4):sylvester.append([0]*shift+fc+[0]*(3-shift))
    for shift in range(3):sylvester.append([0]*shift+gc+[0]*(2-shift))
    assert det_mod(sylvester,19)==2
    corrupted=[row[:] for row in sylvester];corrupted[0]=[0]*7
    assert det_mod(corrupted,19)==0


def boundary_charts():
    x, y, t, u, z, v, s = S.symbols("x y t u z v s")
    q = 2*x*y+t*t
    affine = {x:(z*z-u*u)/(2*z), y:1/z, t:u/z}
    assert S.cancel(q.subs(affine, simultaneous=True)-1) == 0

    # All twenty coefficients of an ordinary representative of degree <=3.
    # The four redundant q-1 directions are retained, so this covers every
    # sixteen-dimensional sphere class before taking the quotient.
    indices = [(i,j,k) for i in range(4) for j in range(4-i)
               for k in range(4-i-j)]
    coeff = {alpha:S.Symbol("c%d%d%d" % alpha) for alpha in indices}
    representative = sum(coeff[i,j,k]*x**i*y**j*t**k for i,j,k in indices)

    def graph(poly, degree):
        result = S.cancel(z**degree*poly.subs(affine, simultaneous=True)-s*z**degree)
        S.Poly(result, u, z)
        return S.expand(result)

    cubic_graph = graph(representative, 3)
    origin = {u:0, z:0}
    assert cubic_graph.subs(origin) == coeff[0,3,0]
    assert S.diff(cubic_graph,u).subs(origin) == coeff[0,2,1]
    assert S.diff(cubic_graph,z).subs(origin) == coeff[0,2,0]
    assert S.diff(cubic_graph,s) == -z**3

    quadratic = sum(coeff[i,j,k]*x**i*y**j*t**k for i,j,k in indices if i+j+k<=2)
    quadratic_graph = graph(quadratic, 2)
    assert quadratic_graph.subs(origin) == coeff[0,2,0]
    assert S.diff(quadratic_graph,u).subs(origin) == coeff[0,1,1]
    assert S.diff(quadratic_graph,z).subs(origin) == coeff[0,1,0]
    assert S.diff(quadratic_graph,s) == -z**2

    def sphere_reduce(poly):
        result = 0
        for (i,j,k), c in S.Poly(S.expand(poly),x,y,t).terms():
            m = min(i,j)
            result += c*x**(i-m)*y**(j-m)*t**k*((1-t*t)/2)**m
        return S.expand(result)

    at_singular_boundary = representative.subs(
        {coeff[0,3,0]:0, coeff[0,2,1]:0, coeff[0,2,0]:0})
    correction = coeff[1,2,0]*y*(q-1)/2
    corrected = S.expand(at_singular_boundary-correction)
    assert S.diff(corrected,y,2) == 0
    assert sphere_reduce(corrected-at_singular_boundary) == 0
    normal = sphere_reduce(corrected)
    profile = S.diff(normal,y)
    assert S.diff(profile,x) == S.diff(profile,y) == 0
    assert S.Poly(profile,t).degree() <= 2
    assert S.diff(normal-y*profile,y) == 0
    assert S.diff(at_singular_boundary-2*correction,y,2) != 0
    quadratic_repeated = quadratic.subs({coeff[0,2,0]:0,coeff[0,1,1]:0})
    assert S.diff(sphere_reduce(quadratic_repeated),y) == coeff[0,1,0]

    # The complete affine-y family and every basis monomial of G(x,t).
    a, b, c = S.symbols("profile_a profile_b profile_c")
    gi = {(i,j):S.Symbol("g%d%d" % (i,j))
          for i in range(4) for j in range(4-i)}
    g = sum(gi[i,j]*x**i*t**j for i,j in gi)
    phase = y*(a+b*t+c*t*t)+g
    profile_graph = graph(phase,3)
    for i,j in gi:
        term = S.cancel(z**3*(x**i*t**j).subs(affine,simultaneous=True))
        expected = S.Rational(1,2**i)*(z*z-u*u)**i*u**j*z**(3-i-j)
        assert S.expand(term-expected) == 0
        assert min(sum(alpha) for alpha,_ in S.Poly(term,u,z).terms()) == 3+i
    jet2 = sum(value*u**i*z**j for (i,j),value
               in S.Poly(profile_graph,u,z).terms() if i+j<=2)
    assert S.expand(jet2-(c*u*u+b*u*z+a*z*z)) == 0

    # Both charts of the boundary-section blow-up.
    first = S.cancel(profile_graph.subs(u,z*v)/z**2)
    second = S.cancel(profile_graph.subs(z,u*v)/u**2)
    S.Poly(first,v,z)
    S.Poly(second,u,v)
    assert S.diff(first,s) == -z
    assert S.diff(second,s) == -u*v**3
    assert S.expand(first.subs(z,0)-(c*v*v+b*v+a)) == 0
    assert S.expand(second.subs(u,0)-(c+b*v+a*v*v)) == 0
    assert S.discriminant(c*v*v+b*v+a,v) == b*b-4*a*c
    assert S.diff(second,v).subs({c:0,u:0,v:0}) == b

    # Literal pullbacks of the residue volume, including the phase divisor.
    first_jacobian = S.det(S.Matrix([[S.diff(z*v,v),S.diff(z*v,z)],[0,1]]))
    second_jacobian = S.det(S.Matrix([[1,0],[S.diff(u*v,u),S.diff(u*v,v)]]))
    first_volume = S.cancel(first_jacobian/(2*z*z))
    second_volume = S.cancel(second_jacobian/(2*(u*v)**2))
    assert first_volume == 1/(2*z)
    assert second_volume == 1/(2*u*v*v)
    assert S.cancel(first_volume*S.diff(first,s)+S.Rational(1,2)) == 0
    assert S.cancel(second_volume*S.diff(second,s)+v/2) == 0
    assert S.cancel(S.diff(cubic_graph,s)/(2*z*z)+z/2) == 0
    assert S.cancel(S.diff(quadratic_graph,s)/(2*z*z)+S.Rational(1,2)) == 0

    # Repeated nonconstant profiles: full generic G, not a truncated model.
    r = S.symbols("profile_root")
    repeated = S.cancel(profile_graph.subs({a:r*r,b:-2*r,c:1})
                        .subs(u,z*(v+r))/z**2)
    S.Poly(repeated,v,z)
    assert S.diff(repeated,s) == -z
    assert repeated.subs(z,0) == v*v
    k = S.expand(S.diff(repeated,z).subs(z,0)+s)
    assert S.diff(k,s) == 0
    remainder = S.cancel((repeated-v*v-z*(k-s))/z**2)
    S.Poly(remainder,v,z)
    critical_value = k.subs(v,0)
    singular_point = {v:0,z:0,s:critical_value}
    assert all(S.diff(repeated,w).subs(singular_point)==0 for w in (v,z,s))
    assert S.det(S.hessian(repeated,(v,z,s)).subs(singular_point)) == -2

    # At the ordinary double point every ambient blow-up chart contributes
    # an exceptional Jacobian of order two, cancelling multiplicity two.
    sigma, exceptional, A, B = S.symbols("sigma exceptional A B")
    shifted = S.expand(repeated.subs(s,critical_value+sigma))
    point_variables = (v,z,sigma)
    quadratic_part = sum(value*v**i*z**j*sigma**k
                         for (i,j,k),value in S.Poly(shifted,*point_variables).terms()
                         if i+j+k==2)
    assert all(sum(alpha)>=2 for alpha,_ in S.Poly(shifted,*point_variables).terms())
    for mapping in (
        {v:exceptional,z:exceptional*A,sigma:exceptional*B},
        {v:exceptional*A,z:exceptional,sigma:exceptional*B},
        {v:exceptional*A,z:exceptional*B,sigma:exceptional},
    ):
        strict = S.cancel(shifted.subs(mapping,simultaneous=True)/exceptional**2)
        S.Poly(strict,exceptional,A,B)
        jacobian = S.det(S.Matrix([[S.diff(mapping[w],q)
                                    for q in (exceptional,A,B)] for w in point_variables]))
        assert S.cancel(jacobian/exceptional**2) in (1,-1)
        exceptional_quadric = S.cancel(quadratic_part.subs(mapping,simultaneous=True)
                                        /exceptional**2)
        assert S.expand(strict.subs(exceptional,0)-exceptional_quadric) == 0

    # New cusp case: constant nonzero profile and nonzero height cubic.
    cusp = S.expand(second.subs({b:0,c:0}))
    h = gi[0,3]
    leading_cusp = (a*v*v+u*(h+gi[0,2]*v+gi[0,1]*v*v+(gi[0,0]-s)*v**3))
    cusp_tail = S.cancel((cusp-leading_cusp)/u**2)
    S.Poly(cusp_tail,u,v)
    assert S.diff(cusp,s) == -u*v**3
    assert S.diff(cusp,u).subs({u:0,v:0}) == h
    assert cusp.subs(u,0) == a*v*v
    assert first.subs({b:0,c:0,z:0}) == a
    cusp_unit = S.cancel((cusp-a*v*v)/u)
    assert cusp_unit.subs({u:0,v:0}) == h
    implicit_ds = -S.diff(cusp,s)/S.diff(cusp,u)
    assert S.cancel(implicit_ds/(2*u*v*v)-v/(2*S.diff(cusp,u))) == 0

    # Orthogonal infinitesimal generator and its polynomial shear.
    metric = S.Matrix([[0,1,0],[1,0,0],[0,0,1]])
    delta_matrix = S.Matrix([[0,0,0],[0,0,-1],[1,0,0]])
    assert delta_matrix.T*metric+metric*delta_matrix == S.zeros(3)
    h0, h1 = S.symbols("shear_h0 shear_h1")
    shear_h = h0+h1*x

    def delta(poly):
        return x*S.diff(poly,t)-t*S.diff(poly,y)

    assert delta(q) == delta(x) == delta(shear_h) == 0

    def derivation(poly):
        return S.expand(-shear_h*delta(poly))

    forward = {x:x,t:t-x*shear_h,y:y+shear_h*t-x*shear_h**2/2}
    inverse = {x:x,t:t+x*shear_h,y:y-shear_h*t-x*shear_h**2/2}
    for w in (x,y,t):
        assert derivation(derivation(derivation(w))) == 0
        assert S.expand(forward[w]-(w+derivation(w)+derivation(derivation(w))/2)) == 0
        assert S.expand(forward[w].subs(inverse,simultaneous=True)-w) == 0
    assert S.expand(q.subs(forward,simultaneous=True)-q) == 0
    wrong = dict(forward)
    wrong[y] = y+shear_h*t-x*shear_h**2
    assert S.expand(q.subs(wrong,simultaneous=True)-q) != 0

    return {
        "outcome":"pass",
        "scope":"Exact graph, blow-up, residue, Hessian, cusp, and orthogonal-shear polynomial identities.",
        "coefficient_representative_dimension":len(indices),
        "sphere_class_dimension":16,
        "G_basis_monomials_checked":len(gi),
        "ordinary_double_point_blowup_charts":3,
        "cusp_hypotheses":"profile_a!=0 and g03!=0; checked local unit coefficient and complementary empty chart",
        "sphere_invariance_argument":"delta is in so(q), delta(h)=0, and O(q)-invariance gives L(h*delta(f))=L(delta(h*f))=0 for every f.",
        "limits":"Properness, resolution smoothness, and V-filtration theorems are proved or cited in the paper, not certified by these identities.",
    }



if __name__=='__main__':
    interval_cubic();constant_profile();boundary_charts()
    print('PASS interval cubics: four moments force zero')
    print('PASS constant-profile family with no t^3 term: all eight-moment identities, the three mod-19 Euclidean remainders, and the optional Sylvester determinant audit')
    print('PASS corrupted determinant rejected; no external polynomial certificate is needed')

    print('PASS exact boundary charts, cusp volume cancellation, ordinary-double-point Hessian, and sphere shear identities')
