# The exceptional polynomial for all orders at fixed signed gaps

Fix \(r\) unequal rows with positive integral gaps \(d_i\) and
orientations \(\epsilon_i\in\{-1,1\}\). Allow the smaller exponent
\(A_i\) of each row to range independently over the nonnegative integers.
Equal rows can be appended freely. Set \(D=\sum_i d_i\),
\(\Delta=\sum_i\epsilon_i d_i\), and
\[
                     D_* = D-\tfrac12\boldsymbol{1}_{\Delta=0}.
\]
Let \(J\ge3\) be any absolute terminal cutoff from the
[effective eventual theorem](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL.md),
and put \(M=J-1\).

**Theorem.** There is a nonzero polynomial \(E(A_1,\ldots,A_r)\) with
integer coefficients and degree \(D_* M(M^2-1)/3\) whose zeros on the
nonnegative integer grid are precisely the profiles that fail SFC(2)
at some order. It is nonnegative on an open real neighborhood of the
nonnegative orthant. The highest homogeneous part has positive
coefficients on its support and is given explicitly below.
In particular, for every integer \(H\ge0\) the number of exceptional
vectors in \(\{0,\ldots,H\}^r\) is at most
\[
                      \frac{D_* M(M^2-1)}6 (H+1)^{r-1}.        \tag{1}
\]
Along every ray with arbitrary nonnegative integer starting
bases and strictly positive integral increments in all the unequal
rows, at most \(D_* M(M^2-1)/6\) integer translations fail SFC(2)
at any order.

The result bounds genuine all-order exceptions at fixed signed gaps.
It also gives a density bound when the gaps vary.

## The highest homogeneous part of one resultant

Define
\[
 F_j(A,z)=\sum_{k=0}^j\binom jk z^k
       \prod_i\prod_{h=1}^{d_i\ell_{j,i}(k)}(jA_i+h),
 \qquad
 \ell_{j,i}(k)=
       \begin{cases}
          k &(\epsilon_i=1),\\
          j-k &(\epsilon_i=-1).
       \end{cases}
\]
For integer \(A_i\ge0\) this is the factorial moment polynomial
with the nonzero row scalars \((jA_i)!\) removed.
Let \(R_N(A)=\operatorname{Res}_z(F_n,F_N)\) with \(n=N-1\).
For \(s_i>0\), write \(D_-=\sum_{\epsilon_i=-1}d_i\) and
\[
     C_s=\sum_i d_i^2/s_i,\qquad
     U_s=\prod_i s_i^{\epsilon_i d_i},\qquad
     B_s=\prod_{\epsilon_i=-1}s_i^{d_i}.
\]
The top homogeneous part is of degree \(nND_*\) and equals
\[
 \operatorname{top}(R_N)(s)
 =(nN)^{nND_-}\left(\prod_i s_i^{nN d_i}\right)
 \begin{cases}
      (n^\Delta-N^\Delta)^{nN},&\Delta\ne0,\\[1mm]
      C_s^{\,nN/2}\,r_{n,N},&\Delta=0,
 \end{cases}                                                   \tag{2}
\]
where \(r_{n,N}>0\) is the resultant of the two monic polynomials
\[
               H_j^0(x)=j!\sum_{p=0}^{\lfloor j/2\rfloor}
                          \frac{(1/(2j))^p}{p!(j-2p)!}x^{j-2p},
                   \qquad j=n,N.
\]
In the balanced case the apparent denominators in (2) cancel.
Each power of \(1/s_i\) is at most \(nN/2\), while each
corresponding prefactor has exponent at least \(nN\).

Here is the leading calculation along an arbitrary direction,
extending the [common-translation formula](MULTIVARIATE_TRANSLATED_RESULTANT_DEGREE.md).
Take \(A=A^0+s/u\) and rescale \(z\) to \(u^\Delta z/U_s\).
Normalize order \(j\) by \(u^{jD_-}/C_j\), with
\(C_j=j^{jD_-}B_s^j\).
The resulting polynomial \(g_j(u,z)\) has \(z^k\)-coefficient
\[
 \binom jk j^{\Delta k}
        \prod_i\prod_{h=1}^{d_i\ell_{j,i}(k)}
                    (1+u(A_i^0/s_i+h/(j s_i))).
\]
Its pairwise resultant is
\[
 \operatorname{Res}(g_n,g_N)
     =u^{nND} U_s^{-nN} C_n^{-N} C_N^{-n} R_N(A^0+s/u).
                                                               \tag{3}
\]
If \(\Delta\ne0\), the constant term of the left side is
\((n^\Delta-N^\Delta)^{nN}\) from the two binomial powers.
This gives the unbalanced case of (2).

If \(\Delta=0\), the logarithm of the product in the \(k\)-th
coefficient has \(u^p\)-coefficient of \(k\)-degree at most
\(p+1\). At order one its quadratic coefficient is \(C_s/(2j)\).
Therefore the \(u^p\)-coefficient of the product has \(k\)-degree
at most \(2p\), with leading coefficient \((C_s/(2j))^p/p!\).
Take \(u=\lambda^2\), \(z=-1+\lambda x\).
The alternating sum with weights \((-1)^k\binom jk\) annihilates
all terms of \(\lambda\)-degree below \(j\), since they have
\(k\)-degree below \(j\). At degree \(j\) it sends \(k^j\)
to \((-1)^j j!\). Thus the exact constant term of
\(\lambda^{-j}g_j(\lambda^2,-1+\lambda x)\) is
\(C_s^{j/2}H_j^0(x/\sqrt{C_s})\).
The [scaled Hermite comparison](MULTIVARIATE_SCALED_HERMITE_INTERLACING.md)
excludes common zeros of \(H_n^0,H_N^0\). The resultant \(r_{n,N}\)
is positive since the even one has no real roots, and its conjugate
root pairs give a product of positive squared moduli.
Scaling the two Hermites gives resultant \(C_s^{nN/2}r_{n,N}\).
The \(x\)-substitution and division by the powers of \(\lambda\)
show that the left side of (3) has order \(nN/2\) in \(u\),
with this leading coefficient. This proves the balanced formula.

Both computations give degree \(nND_*\) along every positive
direction, independently of \(A^0\). They also determine the
multivariate top layer. Indeed a higher homogeneous coefficient
would have to vanish throughout the positive orthant, which is open.
Multiplying the constants in (3) gives (2).

## Partial translations and the smaller support

Translate just a subset \(I\) of the unequal rows in a positive
direction \(s\), leaving the other rows fixed. Let \(D_I\) be the
unsigned gap sum on \(I\), \(\Delta_I\) its signed sum, and
\(D_{I,-}\) its decreasing gap sum. Write \(F_j^{I^c}\) for the
moment polynomial with the rows in \(I\) deleted and its remaining
row factorials stripped as above. If all rows are deleted,
\(F_j^{I^c}=(1+z)^j\).

The coefficient of \(T^{nND_I}\) in \(R_N(A^0+sT)\), with the
direction zero outside \(I\), equals
\[
  (nN)^{nND_{I,-}}\left(\prod_{i\in I}s_i^{nN d_i}\right)
       \operatorname{Res}_z\!
       \left(F_n^{I^c}(n^{\Delta_I}z),\
             F_N^{I^c}(N^{\Delta_I}z)\right).                 \tag{3a}
\]
Indeed apply the normalization of (3) only to the translated
rows. Their \(u=0\) contribution in each coefficient is the same
factor \(j^{\Delta_I k}\) as before. All other factors are
independent of \(u\). The prefactors of the resultant therefore
give (3a). For proper \(I\), this is the generic leading degree:
in (2) the top monomial in the unbalanced case reaches
\(nND_I\); in the balanced case choose all inverse powers
from any one index outside \(I\) to reach the same degree.
There is never a larger degree, since each base exponent in
any resultant term is at most \(nN d_i\) by the coefficient
degrees and the resultant weight in \(k\).

Thus at fixed residual bases the leading degree persists exactly
when the indicated rescaled pair is coprime. For a balanced
inflated subset this is the ordinary consecutive resultant
of the smaller support. Angular separation of the smaller
support is sufficient for every signed sum of the inflated
subset. This gives an exact smaller-support test for the
leading partial-translation coefficient.

## From finitely many resultants to a grid bound

Set \(E=\prod_{N=2}^{M}R_N\). The effective eventual theorem excludes
all later terminals for ordinary exponent vectors, so the zero set of
\(E\) on the ordinary grid is exactly the failure set. Equation (2)
and \(\sum_{N=2}^M N(N-1)=M(M^2-1)/3\) give the degree
and positive coefficients of the top layer.

For real \(A_i>-1/(2J)\) and \(j\le M\), the same coefficient
form is an integral with Gamma factors at positive shapes.
The normalized even moment is strictly positive for real \(z\):
it integrates an even power of a real, nonzero sum of distinct
positive real monomial functions, divided by positive Gamma
scalars. Thus the even moment's roots occur in nonreal conjugate
pairs. Its leading coefficient is positive. By the product
formula for resultants, each \(R_N\) is nonnegative there. Swapping
which order is even does not change the sign, since \(nN\) is even.
In particular every real zero of \(E\) in this open neighborhood
has zero gradient, so its vanishing multiplicity is at least two.
The same is true of zeros along an ordinary positive ray, including
at its starting point. Equation (2) ensures that the polynomial
on each such ray is nonzero of the full degree, proving the
one-dimensional bound.

For completeness, if \(P\ne0\) has total degree \(d\) and \(X\) is a
finite real grid in one dimension, the sum of vanishing multiplicities
over \(X^r\) is at most \(d|X|^{r-1}\). For induction write
\(P(x,t)=P_k(x)t^k+\cdots\), with \(P_k\ne0\) and
\(\deg P_k\le d-k\). At a given \(x\in X^{r-1}\), pick a derivative of
total order \(\operatorname{mult}(P_k,x)\) that takes a nonzero value
there. Apply it to \(P\) in the \(x\)-variables before evaluating at
that \(x\). The result is a nonzero univariate polynomial of degree
\(k\). At every \(t\in X\), the original vanishing multiplicity is
at most \(\operatorname{mult}(P_k,x)\) plus the multiplicity of that
univariate polynomial at \(t\). Summing first over \(t\) and then
over \(x\), and using induction on \(P_k\), proves the claim.
Apply this bound to \(E\), with multiplicity at least two per
exceptional point, to get (1).

## Density with varying gaps

For \(r\) arbitrary rows, allow both exponents of every row to lie
in \(\{0,\ldots,H\}\). For \(H\ge1\) set
\(G=\min\{H,\lfloor\sqrt{27H}\rfloor\}\). The
[absolute curvature 27 theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md)
forces every exceptional profile to have \(|b_i-a_i|\le G\) for
each row. At a fixed signed gap vector, the smaller exponent
in a row is in \(\{0,\ldots,H\}\), so (1) applies. Equal rows can
be reinserted in counting by treating \(E\) as independent of
their variables. Bound \(D_*\) by \(D\). Summing \(D\) over all
signed gap vectors in \(\{-G,\ldots,G\}^r\) gives
\(r(2G+1)^{r-1}G(G+1)\). Hence the total number of all-order
exceptions is at most
\[
   \frac{M(M^2-1)}6\,
      r(2G+1)^{r-1}G(G+1)(H+1)^{r-1}.                        \tag{4}
\]
Let \(B=M(M^2-1)/6\). The exceptional fraction among the
\((H+1)^{2r}\) ordered exponent pairs is at most
\[
    B \frac{rG(G+1)}{(H+1)^2}
                  \left(\frac{2G+1}{H+1}\right)^{r-1}.
\]
For fixed \(r\) this is \(O_r(H^{-(r+1)/2})\).
For \(H\ge434\) it is at most \(28B/(H+1)\) uniformly in \(r\).
Indeed \(4\sqrt{27H}\le H-1\), \(G(G+1)\le28H\), and
\(r/2^{r-1}\le1\). The exact algebraic locus still requires
excluding its integer points to prove full SFC(2).
