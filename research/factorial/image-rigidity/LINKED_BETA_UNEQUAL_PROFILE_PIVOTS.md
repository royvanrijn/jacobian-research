# Unequal linked beta profiles and the first two Euclidean remainders

For every nonnegative integer pair \(r,u\) and every \(n\ge1\),
the first Euclidean remainder of \(J_{n+1}(z;r,u)\) by \(J_n(z;r,u)\)
has degree \(n-1\) and positive leading coefficient. For \(n\ge2\),
divide \(J_n\) by the monic normalization of that remainder.
The next remainder also has its expected degree, \(n-2\), and positive
leading coefficient. These are two steps of
`GAP-FACTORIAL-LINKED-BETA-REGULAR-CHAIN`. The proof applies to all
positive integer ratios in either orientation. The later remainders
remain an obligation at general starting orders.

Here
\[
 J_n(z;r,u)=\mathbb E[(z-T_n)^n],\qquad
 T_n\sim\operatorname{Beta}(nr+1,nu+1).
\]
This is the [actual phase reduction](canonical-density/moment_pair_obstructions.md)
for \(X^rY^u((z-1)X+zY)\). The two consecutive densities change with
the order on the same fixed exponent profile.

## 1. The first two Euclidean remainders

Let \(\mu_j,v_j\) be the mean and variance of \(T_j\). Comparing the
top three descending coefficients gives the leading coefficient of
the first remainder before monic normalization:
\[
 L=\frac{n(n+1)v_{n+1}-n(n-1)v_n
           +n(n+1)(\mu_{n+1}-\mu_n)^2}{2}.                         \tag{1}
\]
For example, if \(F_i,G_i\) denote the descending coefficients of
\(J_n,J_{n+1}\), start with \(F_0=G_0=1\) and set coefficients
below degree zero to zero. The quotient is \(z+q\) with
\(q=G_1-F_1\), and its remainder begins at coefficient
\(L=G_2-F_2-qF_1\). Substituting means and second moments gives (1).

Put \(s=r+u\). For a real variable \(x\ge1\), the required
monotonicity follows from the exact decomposition
\[
\begin{aligned}
 \frac{x(x-1)(rx+1)(ux+1)}{(sx+2)^2(sx+3)}
 ={}&
 \left(\frac{x}{sx+2}\right)^2\left(1-\frac1x\right)
                                     \frac{sx+1}{sx+3}\\
 &+ru\left(\frac{x}{sx+2}\right)^2\frac{x(x-1)}{sx+3}.
\end{aligned}
\]
Every factor is nonnegative and nondecreasing. The first product
increases strictly, since its first two factors increase strictly
and its third is positive. For the last factor of the second
product, write \(x(x-1)/(sx+3)=(x-1)x/(sx+3)\).
Thus \(x(x-1)v_x\) increases strictly for all real \(r,u\ge0\),
which proves \(L>0\) in (1).

Define \(M=G_3-F_3-qF_2\), \(K=G_4-F_4-qF_3\). On dividing
\(J_n\) by the first remainder normalized by \(L\), the coefficient
in degree \(n-2\) is
\[
                         F_2-K/L-F_1 M/L+(M/L)^2.                 \tag{2}
\]
Indeed, the divisor has descending coefficients \(1,M/L,K/L\),
and the quotient is \(z+F_1-M/L\). At \(n=2\), \(K=0\) because
the first remainder is linear, so (2) still applies.

To certify positivity of (2), let
\[
 D=24\prod_{i=2}^5(ns+i)((n+1)s+i).
\]
The scaled coefficients \(\widehat F_i=D F_i,\widehat G_i=D G_i\)
are polynomials over the integers for \(0\le i\le4\).
The scaled remainder coefficients \(\widehat L,\widehat M,\widehat K\)
are \(D^2\) times their unscaled versions. The numerator of (2)
over the positive denominator \(D\widehat L^2\) is
\[
 {\cal N}=\widehat F_2\widehat L^2-D\widehat K\widehat L
            -\widehat F_1\widehat M\widehat L+D\widehat M^2.        \tag{3}
\]
We may assume \(r\le u\). Reflection sends the monic
\(J_j(z;r,u)\) to \((-1)^jJ_j(1-z;r,u)=J_j(z;u,r)\).
It preserves the positive leading coefficient of each of the
first two remainders whenever their degrees are \(n-1,n-2\):
the alternating scalar from the dividend cancels the sign
from the variable reflection at the stated remainder degree.

In the sorted cases \((r,u)=(0,v)\) and
\((r,u)=(1+y,1+y+v)\), put \(n=x+2\). For each substitution,
the fully expanded numerator (3) has nonnegative coefficients in
\(x,y,v\) and strictly positive constant coefficient. The first
case has 767 nonzero coefficients and multidegrees \((36,0,34)\).
The second has 30,135 nonzero coefficients and multidegrees
\((40,37,34)\). The following exact polynomial-ring calculation
constructs the entire numerator and checks every coefficient.
It uses no sampled sign criterion.

```python
from math import factorial, prod
from sympy import ZZ
from sympy.polys.rings import ring

R, x, y, v = ring("x,y,v", ZZ)

def numerator(r, u):
    n, s = x + 2, r + u
    lower = [n*s+k for k in range(2, 6)]
    upper = [(n+1)*s+k for k in range(2, 6)]
    D = 24*prod(lower)*prod(upper)

    def scaled(order, first, second):
        return [
            (-1)**k*(24//factorial(k))*prod(first[k:])*prod(second)
            *prod((order-i)*(order*r+1+i) for i in range(k))
            for k in range(5)
        ]

    F, G = scaled(n, lower, upper), scaled(n+1, upper, lower)
    q = G[1]-F[1]
    L, M, K = [D*(G[i]-F[i])-q*F[i-1] for i in (2, 3, 4)]
    return F[2]*L**2-D*K*L-F[1]*M*L+D*M**2

for r, u, count, degrees, constant in [
    (R.zero, v, 767, (36, 0, 34), 11412608573767680000000000),
    (1+y, 1+y+v, 30135, (40, 37, 34),
     5081529250690564309014498958009958400000),
]:
    p = numerator(r, u)
    assert len(p) == count
    assert tuple(max(e[j] for e in p) for j in range(3)) == degrees
    assert all(value > 0 for value in p.values())
    assert p[0, 0, 0] == constant
```

This proves (2) for every integer \(r,u\ge0\), and also for
the real profiles with one exponent zero or both exponents at least
one. For orders \(n\le2\) the two formulas already prove the full
chain. For \(n=3\), the remaining constant follows from the
[known binary coprimality at this order](SFC_SHIFTED_RESULTANT_POSITIVITY_AUDIT.md),
`FACTORIAL-BINARY-SHIFTED-RESULTANT-POSITIVE-WINDOWS-THROUGH-FIVE`.
The all-order positive signs here use the remainder convention
\(G=(z+q)F+R\), with no minus sign inserted in \(R\).

## 2. Why the unit-circle argument does not extend

Mao–Wang–Zheng’s [Remark 3.4 and proof in Sections 2–3](https://arxiv.org/pdf/2609.34551v1)
apply to positive palindromic input polynomials after normalizing
their coefficients by \(\binom{n}{k}\) and checking log-concavity.
The symmetry gives a half-polynomial split with a reflected conjugate;
the normalized inequalities give the estimate on the first coefficient
of its reversed half and the monotone sequence used in Abel summation.
Fejér’s inequality bounds the half-polynomial’s argument between
zero and \(\theta/2\). The resulting finite Blaschke phase difference
then compares the two unit-circle root lists.

For the actual linked coefficient polynomial \(g_n\), the \(k\)-th
coefficient is
\[
      \frac{(nr+1)_k(nu+1)_{n-k}}{k!(n-k)!}.
\]
If \(r\ne u\), the two endpoint coefficients are unequal,
so even equality of their absolute values, necessary for all roots
to lie on the unit circle, fails. For \(r=u\), take the
coefficientwise reciprocal as the candidate input to the paper,
so the desired output is \(g_n\). Writing \(a=nr+1\), the first
normalized coefficient divided by the zeroth is
\((a+n-1)/(n^2a)\le1/n<1\) for \(n\ge2\).
Symmetry and log-concavity would require ratio at least one,
as used in the proof of Remark 3.4. Thus its proof step fails
even on the diagonal. The [direct Sturm proof](GAMMA_QUOTIENT_TRANSFER_AUDIT.md)
obtains diagonal interlacing by another argument.

A change to an affine line cannot extend that Sturm interlacing
to unequal profiles. For \(n\ge3\), the roots of an unequal
\(J_n\) are not collinear. A real line containing every root
must be invariant under conjugation. The real axis is impossible:
even \(J_n\) is strictly positive on it, and the derivative of
odd \(J_n\) is strictly positive there, giving one simple real root.
A conjugation-invariant vertical line through every root would
have real coordinate \(\mu_n\), the mean of \(T_n\), by the sum
of roots. Conjugate pairing of roots on that line would make
\(J_n(\mu_n+w)\) even or odd in \(w\). However, its coefficient
of \(w^{n-3}\) is
\(-\binom n3\mathbb E[(T_n-\mu_n)^3]\), and for a
\(\operatorname{Beta}(\alpha,\beta)\) variable the latter central
moment is
\[
 \frac{2\alpha\beta(\beta-\alpha)}
 {(\alpha+\beta)^3(\alpha+\beta+1)(\alpha+\beta+2)} .
\]
The identity follows by substituting the first three
rising-factorial moments. It is nonzero on the unequal profiles.
Simplicity of all roots, if needed to identify the invariant line,
follows independently below.

## 3. The exact Jacobi receiver and the individual discriminant

For arbitrary fixed positive \(\alpha,\beta\), put
\(s=\alpha+\beta\), \(P_n(z)=\mathbb E[(z-T)^n]\) under
\(T\sim\operatorname{Beta}(\alpha,\beta)\).
Integrate the derivative of
\(t^\alpha(1-t)^\beta(z-t)^{n-1}\) on \((0,1)\).
The endpoint terms vanish. On replacing \(t=z-(z-t)\), the result is
\[
 z(z-1)P_n''+[\alpha+n-1-(s+2n-2)z]P_n'
                          +n(n+s-1)P_n=0 .                        \tag{4}
\]
Both endpoints are nonzeros. At any other point, a repeated
root would give zero value and derivative for a regular
second-order ODE, contradicting uniqueness for a monic polynomial.

In fact (4) gives a product certificate for this simplicity.
Keep \(\alpha,\beta\) fixed and write \(D_n=\operatorname{Res}(P_n,P_n')\).
The monic derivative normalization is \(p=P_n'/n=P_{n-1}\).
At each root \(x\) of \(p\), (4) reads
\(P_n(x)=-x(x-1)p'(x)/(n+s-1)\). Multiply over those
roots, use \(\prod_x x(x-1)=p(0)p(1)\), and note
that swapping polynomials of degrees \(n,n-1\) has even
resultant sign. Therefore
\[
 D_1=1,\qquad
 D_n=\frac{n^n(\alpha)_{n-1}(\beta)_{n-1}}
           {(n+s-1)^{n-1}(s)_{n-1}^2}\,D_{n-1}>0\quad(n\ge2).
                                                                    \tag{5}
\]
The exact product applies to each physical order by taking
\(\alpha=nr+1,\ \beta=nu+1\) and fixing those values inside the
recursion. Its derivative step does not produce the lower-order
physical moment on the changing linked profile.

There is also an exact comparison to the primary
[Martínez–Finkelshtein–Orive Riemann–Hilbert paper](https://arxiv.org/pdf/math/0410320).
For the actual linked pair the same Möbius transformation takes
each polynomial to
\[
 J_n(z;r,u)=
 \frac{n!}{(n(r+u)+2)_n}(1-z)^n
 P_n^{(-n(r+1)-1,\ n(r+u)+1)}
                    \left(\frac{1+z}{1-z}\right),                 \tag{6}
\]
where the right side uses standard Jacobi normalization.
To see this from the finite Jacobi binomial formula, after
multiplication by \((1-z)^n\) its \(z^j\) coefficient is
\(\binom{-nr-1}{n-j}\binom{n(r+u+1)+1}{j}\). The prefactor
in (6) converts it to the coefficient in \(J_n\).
Neither \(z=0\) nor \(z=1\) is a root, so the transformation
preserves actual common roots.

For each **fixed** positive \(r,u\), the Jacobi parameters
divided by \(n\) tend to \(A=-(r+1)\), \(B=r+u\), with
\(A<-1\) and \(A+B=u-1>-1\). This is the paper’s case C.2
in (1.6)–(1.7). The negative integer first parameter is below
\(-n\); it is outside the endpoint degeneracies in (1.3)–(1.5).
Theorem 1, (2.15), gives a one-term expansion off the limiting
arc and its auxiliary cuts. Theorem 2(a), (2.17), gives two
terms on either side of the arc away from the endpoints.
After multiplying by the nonzero common factors, the zero
equation near a regular point has the form
\[
     H_n(\xi)^{2n}=\mp\frac{N_{12,n}(\xi)}{N_{11,n}(\xi)}
                                 \bigl(1+O(1/n)\bigr),
\]
with \(H_n=G_nw_n\) as in (4.3)–(4.4).
To compare consecutive roots one must use compatible local
branches for both degrees. On any shared patch where the
two-term formulas are uniform and the local phase and
amplitude ratio depend smoothly on the scaled parameters,
subtract the logarithms of the two root equations modulo
\(2\pi i\). The physical paths \(A_n=A-1/n,\ B_n=B+1/n\)
then force \(2\log H(\xi)\to0\) modulo \(2\pi i\) for a common
root sequence in that patch. Matching the sheets across moving
lenses and proving a uniform comparison near the Airy endpoints
require additional work. The expansion in Theorem 3, (2.23),
is an individual-order Airy approximation. Its storage and
scaling alone do not provide the adjacent-order error bound
for root separation at `GAP-FACTORIAL-BINARY-LINKED-BETA-COPRIMALITY`.

## Verification

The coefficient check above was run with SymPy’s exact polynomial
ring and was cross-checked with integer FLINT arithmetic.
Independent rational polynomial division checked the first two
remainders on 96 complete pairs at starting orders \(2,\ldots,17\).
Direct polynomial computation checked (4)–(6), including the full
value of the resultant in (5), on 48 individual moments of orders
\(1,\ldots,8\). Both checks used the profiles
\((0,0),(0,4),(1,2),(3,4),(7,2),(17,204)\).
