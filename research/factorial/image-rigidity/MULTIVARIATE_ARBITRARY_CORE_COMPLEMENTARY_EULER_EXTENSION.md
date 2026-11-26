# All-degree complementary extensions of an arbitrary binary core

Fix any two distinct ordinary monomials \(U,V\). Let \(D_+\) and
\(D_-\) be the sums of their positive and negative signed exponent
gaps, respectively, counted with absolute value. Take \(J\ge4\)
to be a computable
[absolute eventual ordinary cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL.md),
so every terminal \(N\ge J\) is regular. Set
\[
                M=\min\bigl(\lfloor(J-2)/2\rfloor,\max(D_+,D_-)\bigr).
\]
For \(a\) a positive integer, append \(M\) independent pairs of
coordinates with exponent rows
\[
             (a^s,a^s+1),\quad(Q-a^s,Q-a^s-1),
                          \qquad 1\le s\le M .                   \tag{1}
\]

**Theorem.** There is a computable \(a_0\) for the retained core
such that, for every \(a\ge a_0\) and every prime \(p\), there is a
computable \(L\) with the following property. Every positive integer
\(Q\) divisible by \(p^L\) gives an ordinary support (1) satisfying
SFC(2) at every degree. The modulus can ensure \(Q\ge2(a^M+1)\).
The additions preserve the original total-degree difference and
their total absolute curvature tends to zero uniformly in these
\(Q\) as \(a\) increases.

The argument separates the Euler comparisons even if the retained
core has multiple common roots. The number of additional pairs is
linear in the cutoff and can instead be bounded by the total gap
in one orientation.

## A degree-independent bound on multiplicities

Let \(P_j(z)=\sum_{k=0}^j c_{j,k}z^k\) be the order-\(j\)
core moment. For a row with increasing gap \(d\) and smaller
exponent \(b\), advancing \(k\) to \(k+1\) multiplies its factorial
by \(\prod_{h=1}^{d}(jb+dk+h)\). A decreasing row with the same
gap and smaller exponent divides by
\(\prod_{h=0}^{d-1}(jb+d(j-k)-h)\).
Hence there are real polynomials \(A,B\) of degrees \(D_+,D_-\)
such that
\[
    \frac{c_{j,k+1}}{c_{j,k}}
       =\frac{(j-k)A(k)}{(k+1)B(k)}\qquad(0\le k<j).
\]
The binomial factor accounts for \((j-k)/(k+1)\).
With \(D=z\,d/dz\), coefficient comparison gives the exact equation
\[
                    [D B(D-1)-z(j-D)A(D)]P_j=0.                \tag{2}
\]
The uppermost coefficient when this is written in ordinary
derivatives has the form \(z^R(\lambda+\mu z)\), where
\(R=1+\max(D_+,D_-)\), and \(\lambda,\mu\) are real and not both
zero. This follows from the coefficient of the highest power of
\(D\) and \(D^R=z^R(d/dz)^R+\) lower derivatives.
At any nonreal root (2) has a nonzero leading coefficient. If
\(P_j\) vanished to order at least \(R\) there, the term containing
the highest derivative would have the unique lowest Taylor order.
Thus nonreal root multiplicity is at most \(R-1\).

In every adjacent pair the even moment is strictly positive on
the real axis, by its integral representation as an even power
of a nonzero polynomial in independent exponential variables.
Every common root is nonreal. Its multiplicity in degree \(j\)
is consequently at most both \(\lfloor j/2\rfloor\) and \(R-1\),
where the first bound uses conjugation. At a terminal \(N<J\),
the smaller of its two multiplicities is at most \(M\).

## Separating all core collisions in the Euler comparison

For \(t=1/a\) consider the coefficientwise polynomial family
\[
              \mathcal E_j(t,z)
                  =\prod_{s=1}^{M}(1+t^s D/j)\,P_j(z).
                                                               \tag{3}
\]
At a root \(\zeta\ne0\) of \(P_j\), write
\(P_j(\zeta e^y)=b y^h+O(y^{h+1})\), where \(h\) is its
multiplicity and \(b\ne0\). In this coordinate \(D=d/dy\).
The coefficient of \(D^r\) in (3) is
\(j^{-r}e_r(t,t^2,\ldots,t^M)\), whose first term is
\(j^{-r}t^{r(r+1)/2}\).

At scale \(y=t^s x\), for \(1\le s\le\min(M,h)\), a term with
\(r\le h\) has lowest possible \(t\)-order
\[
                           r(r+1)/2+s(h-r).
\]
This is minimized exactly at \(r=s-1,s\), since the successive
difference is \(r-s\). Terms with \(r>h\) have order at least
\(r(r+1)/2\), greater than that minimum. All higher Taylor terms
in the coefficients or in \(P_j(\zeta e^y)\) also have strictly
higher order. After dividing by this minimum power, the limiting
polynomial in \(x\) is
\[
                     C x^{h-s}\left(x+\frac{h-s+1}{j}\right),
                             \qquad C\ne0.                     \tag{4}
\]
Its simple nonzero root gives exactly one branch at this scale,
with \(y/t^s\to-(h-s+1)/j\), by Rouché's theorem on a small
fixed disk around that root. These branches are distinct at
different scales. If \(h>M\), the derivative of order \(M\)
leaves the exact root \(\zeta\) with multiplicity \(h-M\) for
\(t\ne0\), since it has strictly lower \(y\)-order than all smaller
derivatives. Together the described roots exhaust the cluster
by continuity of the roots of the degree-\(j\) polynomial.

Now take \(N=n+1<J\). Clusters at distinct core roots stay
separate. At a shared core root, branches in different scales
cannot coincide. If two branches are in the same scale \(s\),
their respective velocities are
\(-(h_n-s+1)/n\) and \(-(h_N-s+1)/N\).
Each positive numerator is strictly less than its denominator,
by the multiplicity bound. Because \(n,N\) are relatively
prime, the velocities are different. Finally, at least one
of the orders has no stationary branch at this common root,
since \(\min(h_n,h_N)\le M\).
It follows that (3) is coprime at every terminal below \(J\)
for all sufficiently small positive \(t\).

More precisely, the vanishing order at \(t=0\) of the resultant
of (3) at a given pair is
\[
            \sum_{\zeta:\,P_n(\zeta)=P_N(\zeta)=0}
               \sum_{s=1}^{M}
                   (h_n(\zeta)-s+1)_+(h_N(\zeta)-s+1)_+ ,      \tag{5}
\]
where \(x_+=\max(x,0)\). Distinct cluster centers contribute
valuation zero. Within one common cluster the valuation of
a root difference is the smaller of its two scale exponents,
assigning exponent infinity to a stationary branch. At equal
scale the distinct velocities in (4) give the same valuation.
Summing the minimum counts the pairs of branches both at scale
\(s\) or smaller displacement for each \(s\). One of the
multiplicities is at most \(M\), so there are no terms beyond
\(M\). The leading coefficients have nonzero constant terms,
which proves (5) by the product formula.

## Exact ordinary lifts along a progression of totals

The integer Euler comparison for (1) is
\[
                    E_j(a,z)=\prod_{s=1}^M(ja^s+D)\,P_j(z).
\]
At \(a>0\) this differs from (3) with \(t=1/a\) by a nonzero
scalar. Each of its adjacent resultants below \(J\) is thus a
nonzero integer polynomial in \(a\).
To bound \(a_0\) without evaluating these resultants, let
\(u_i,v_i\) be the exponent columns of the core and put
\[
 C=2^J\prod_i(J\max(u_i,v_i))!,\qquad
 B_0=(2J)!\,[C(2J)^M]^{2J}.
\]
Every core coefficient is at most \(C\). The nonnegative
coefficients of each polynomial \([z^k]E_j(a,z)\) sum to at
most \(C(2J)^M\) by evaluation at \(a=1\).
Expanding the Sylvester determinant bounds the sum of absolute
coefficients of every Euler resultant by \(B_0\).
Its leading coefficient is a nonzero integer, so the Cauchy
root bound allows \(a_0=B_0+2\).
Equation (5) also gives its exact degree in \(a\) by subtracting the
stated vanishing order from \((n+N)M(M+1)/2\).

Fix \(a\ge a_0\) and a prime \(p\). The
[exact complementary-row identity](MULTIVARIATE_COMPLEMENTARY_ROW_EULER_REDUCTION.md#a-complementary-pair-over-a-retained-core)
compares each pair in (1) with \(ja^s+D\) after a scalar
normalization and \(z\mapsto-z\). All \(M\) signs combine into
the same change \(z\mapsto(-1)^M z\) in every degree.
If \(p^L\mid Q\), the normalized augmented moments are
coefficientwise congruent to \(E_j(a,z)\) with precision
\(L-\lfloor\log_p(J(a^M+1))\rfloor\), simultaneously for
\(j\le J\). The product of the coefficientwise comparisons
has that same precision since all factors are \(p\)-integral.
At fixed \(a\), the same determinant estimate bounds each
nonzero integer Euler resultant in absolute value by
\[
                        B_a=(2J)!\,[C(J(a^M+1))^M]^{2J}.
\]
Thus one may take
\[
      L=\lfloor\log_p(J(a^M+1))\rfloor+\lfloor\log_p B_a\rfloor+1.
\]
The resulting precision exceeds all the required valuations,
and \(p^L>J(a^M+1)\ge2(a^M+1)\).
Congruence of the Sylvester determinants excludes actual
common roots below \(J\). At every later terminal, (1) is
ordinary and the absolute eventual theorem applies.

Each added pair has gap \(+1\) and gap \(-1\). If
\(Q\ge2(a^M+1)\), its two curvature terms are
\(1/(a^s+1)\) and \(1/(Q-a^s)\), at most
\(1/(a^s+1)\) and \(1/(a^M+2)\), respectively. Their sum
over \(1\le s\le M\) tends to zero. The coprimality statements
transfer from \(\mathbb Q\) to every characteristic-zero field
by rational Bézout identities.
