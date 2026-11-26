# A Jacobi strip certificate for opposed unit factorial rows

For the ordinary exponent rows \((3,4),(2,1)\), the consecutive
moment polynomials are coprime in every degree. The conclusion also
holds after interchanging the columns, interchanging the rows, or
adding rows with equal exponents. The proof combines a real-part
bound for a Jacobi recurrence matrix with a linked polynomial of
degree four. It uses the same construction for any two opposed unit
rows, with a linked polynomial of degree bounded by the sum of the
two smaller exponents.

## The Jacobi strip

For nonnegative integers \(a,b\), take the rows
\((a,a+1),(b+1,b)\). Division by \((ma)!(mb)!\) gives the
moment polynomial
\[
 C_m(z)=\sum_{k=0}^m\binom mk (ma+1)_k(mb+1)_{m-k}z^k .
\]
Here \((v)_k\) denotes the rising factorial. Equal exponent rows
contribute a nonzero scalar to every moment, so they can be dropped.

Write \(\alpha=ma+1\), \(\beta=mb+1\), \(T=\alpha+\beta\) and
\(D=\alpha-\beta\). In the Cayley coordinate \(x=(z+1)/(z-1)\),
the monic polynomial
\[
 F_m(x)=\frac{(x-1)^m}{(T)_m}
        C_m\!\left(\frac{x+1}{x-1}\right)
       =\mathbb E(x+Y)^m,\qquad Y=2V-1,\quad
         V\sim\operatorname{Beta}(\alpha,\beta),
\]
has nonzero values at \(x=\pm1\). The representation follows by
expanding \((x+2V-1)^m=[(x+1)V+(x-1)(1-V)]^m\) and integrating.
Positivity gives \(C_m(1)>0\), so every original root corresponds to
a finite Cayley root. Integration by parts with the beta density
gives
\[
 (1-x^2)F_m''+[D+(T+2m-2)x]F_m'-m(T+m-1)F_m=0 .
\]
Indeed, if \(w\) is the density of \(Y\), then
\(((1-y^2)w)'=(D-Ty)w\), and \((1-y^2)w\) vanishes at the
endpoints. Apply this identity to \((x+y)^{m-1}\) and use
\(y=(x+y)-x\). The differential equation identifies \(F_m\) with
the monic Jacobi polynomial of degree \(m\) and parameters
\(A=-\alpha-m\), \(B=-\beta-m\). The leading factor in the
differential expression on a trial monomial of degree \(j<m\) is
\((j-m)(T+m-1-j)\ne0\), which also shows uniqueness at fixed leading
coefficient. The equation is regular at every root, so a multiple
root would contradict uniqueness for the ordinary differential
equation with zero value and derivative. [NIST DLMF, Jacobi differential
equation](https://dlmf.nist.gov/18.8.T1)

Set \(S=T+2m\). At these fixed Jacobi parameters the monic
three-term recurrence has coefficients
\[
 b_k=-\frac{DS}{(S-2k)(S-2k-2)},\quad 0\le k<m,\qquad
 c_k=\frac{4k(k+A)(k+B)(k+A+B)}
 {(2k+A+B)^2(2k+A+B+1)(2k+A+B-1)}<0,\quad 1\le k<m .
\]
The formulas follow from the Jacobi recurrence by dividing each
polynomial by its leading coefficient. The leading coefficients
through degree \(m\) are nonzero since their reciprocals in the
monic normalization are proportional to
\(\binom{2k+A+B}{k}^{-1}\), with \(2k+A+B<0\).
All displayed recurrence denominators are nonzero.
[NIST DLMF, Jacobi recurrence](https://dlmf.nist.gov/18.9.T1)

Thus \(F_m\) is the characteristic polynomial of a tridiagonal
matrix with real diagonal \(b_k\) and symmetric off-diagonal
entries \(i\sqrt{-c_k}\). For an eigenvector \(v\), the real part of
\(v^*Jv/(v^*v)\) is the weighted average of its diagonal entries
with weights \(|v_k|^2/(v^*v)\). Every eigenvector has nonzero first
and last coordinates, by the tridiagonal equations. If \(a>b\) and
\(m\ge2\), the diagonal is strictly decreasing, and every root
satisfies
\[
 L_m<\Re x<U_m,\qquad
 L_m=-\frac{DS}{(S-2m+2)(S-2m)},\qquad
 U_m=-\frac{D}{S-2}.
 \tag{1}
\]
When \(a=b\), the diagonal vanishes and all Cayley roots have real
part zero.

## The linked polynomial

Fix the upper order \(N=n+1\), and use \(\alpha,\beta,T,D\) at order
\(N\) in this section. Let
\[
 g(x)=(x-1)^{-\alpha}(x+1)^{-\beta},\qquad
 p(x)=(x-1)^a(x+1)^b.
\]
By Leibniz's rule for the two factors, \(F_N\) has the roots of
\(g^{(N)}\) off the endpoints, and \(F_n\) has the roots of
\((pg)^{(n)}\) there. At a root of \(F_N\), \(g^{(n)}\ne0\). To see
this, differentiate the identity
\((x^2-1)g'+(Tx+D)g=0\) a total of \(j\) times to get
\[
 (x^2-1)g^{(j+1)}+[(T+2j)x+D]g^{(j)}
               +j(T+j-1)g^{(j-1)}=0 .
\]
Two consecutive zero derivatives \(g^{(n)},g^{(n+1)}\) would give
zero derivatives down to \(g\), which is impossible off the
endpoints.

At \(F_N(x)=0\), the ratios
\(r_k=(n)_{\underline{k}}\,g^{(n-k)}/g^{(n)}\), where
\((n)_{\underline{k}}=n!/(n-k)!\), therefore satisfy the polynomial
recursion
\[
\begin{split}
 r_0&=1,\qquad r_1=-\frac{(T+2n)x+D}{T+n-1},\\
 r_k&=-\frac{(x^2-1)(j+1)r_{k-2}+[(T+2j)x+D]r_{k-1}}{T+j-1},
                 \qquad j=n-k+1,\quad 2\le k\le n .
\end{split}
\]
Leibniz's rule now makes every consecutive common root a zero of
\[
 H_{n;a,b}(x)=
 \sum_{k=0}^{\min(n,a+b)}\frac{p^{(k)}(x)}{k!}r_k(x).
 \tag{2}
\]
This polynomial has degree at most \(a+b\). For \(a>b\), equal
counts of zeros of (2) strictly to the right of \(L_N\) and \(U_N\)
exclude every common root by (1). In particular, the two half-plane
tests have size bounded in terms of the row exponents while the
moment order increases. The fixed-degree comparison is related to
the [linked beta remainder](LINKED_BETA_FIXED_DEGREE_REMAINDER.md),
which uses the lower polynomial as its reference.

## Every order for \(a=3,\ b=1\)

For \(n\ge4\), (2) is
\(-4(n+1)Q/[5(5n+2)(5n+3)(5n+4)]\), where
\[
\begin{split}
 Q(x)={}&(16n^2+20n+6)x^4-(52n^2+63n+18)x^3\\
        &+(34n^2+31n+6)x^2+(72n^2+115n+42)x\\
        &+186n^2+217n+60 .
\end{split}
\]
For each \(l\in\{L_N,U_N\}\), let
\(x^4+q_3x^3+q_2x^2+q_1x+q_0\) be the monic normalization of
\(Q(x+l)\). Both shifts satisfy
\[
 q_3<0,\qquad q_0>0,\qquad
 \Delta_2=q_3q_2-q_1<0,\qquad
 \Delta_3=\Delta_2q_1-q_3^2q_0<0 .
 \tag{3}
\]
The exact replay below verifies these signs by positive
coefficients after translating \(n\) by four.

Any real monic quartic satisfying (3) has two roots in the right
half-plane, counted with multiplicity. Substituting \(x=iy\) shows
that an imaginary-axis root would force \(q_0=0\) if \(y=0\), or
\(\Delta_3=0\) if \(y\ne0\). The coefficient set in (3) is
connected. Indeed write \(u=-q_3>0,\ v=q_0>0,\ w=-\Delta_2>0\) and
\(z=q_1+u^2v/w>0\). Then \(q_1=z-u^2v/w\) and
\(q_2=(w-q_1)/u\) express every such quartic in positive parameters.
The half-plane count is constant by continuity of roots. The
quartic \((x^2-x+1)(x^2+x/2+1)\) satisfies (3) and has two roots in
each half-plane. Therefore \(Q\) has no root in the strip (1), and
(2) excludes common roots for \(n\ge4\). The remaining three pairs
are verified exactly below.

~~~python
import sympy as sp

n, x, t = sp.symbols("n x t")
N = n + 1
T, D = 4*N + 2, 2*N
S = T + 2*N
lower = -D*S/((S-2*N+2)*(S-2*N))
upper = -D/(S-2)
p = (x-1)**3*(x+1)
r = [sp.S.One, -((T+2*n)*x+D)/(T+n-1)]
for k in range(2, 5):
    j = n-k+1
    r.append(sp.cancel(-((x*x-1)*(j+1)*r[k-2]
                         + ((T+2*j)*x+D)*r[k-1])/(T+j-1)))
H = sum(sp.diff(p, x, k)*r[k]/sp.factorial(k) for k in range(5))
Q = ((16*n*n+20*n+6)*x**4 - (52*n*n+63*n+18)*x**3
     + (34*n*n+31*n+6)*x*x + (72*n*n+115*n+42)*x
     + 186*n*n+217*n+60)
assert sp.cancel(H + 4*N*Q/(5*(5*n+2)*(5*n+3)*(5*n+4))) == 0
for endpoint in (lower, upper):
    monic = sp.cancel(Q.subs(x, x+endpoint)/(16*n*n+20*n+6))
    one, q3, q2, q1, q0 = sp.Poly(monic, x).all_coeffs()
    assert one == 1
    delta2 = q3*q2-q1
    delta3 = delta2*q1-q3*q3*q0
    for positive in (-q3, q0, -delta2, -delta3):
        numerator, denominator = sp.fraction(sp.cancel(positive))
        for polynomial in (numerator, denominator):
            assert all(c > 0 for c in
                       sp.Poly(polynomial.subs(n, t+4), t).all_coeffs())

def moment(m):
    return sum(sp.binomial(m, k)*sp.rf(3*m+1, k)*sp.rf(m+1, m-k)*x**k
               for k in range(m+1))

for m in range(1, 4):
    assert sp.degree(sp.gcd(moment(m), moment(m+1)), x) == 0
print("Two symbolic strip boundaries and three low-order pairs passed")
~~~
