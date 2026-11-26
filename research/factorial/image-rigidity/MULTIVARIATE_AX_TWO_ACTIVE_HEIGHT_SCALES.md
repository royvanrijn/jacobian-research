# All-degree exclusion at two active polynomial factorial height scales

Fix any ordinary two-monomial core, \(d\ge2\), integers
\(1\le a<b\) with \(v_2(a)\ne v_2(b)\), and an integer
\(K\ge0\). Put \(U=M^a,\ V=M^b+K\) and append the rows
\[
           (U,U+d),\quad (U+1,U),\quad
           (V,V+1),\quad (V+d,V)
\]
with counts \(r,\ dr,\ dr M^{b-a},\ rM^{b-a}\), respectively.
Here \(M,r\) are positive integers. There is a computable
height threshold, independent of \(r\), above which the
original ordinary binomial satisfies physical SFC(2) in
every degree. The core need not itself be regular. At
\(r/M^a\) tending to a positive limit both scales contribute
positive curvature. The
[polynomial divisor-rank theorem](MULTIVARIATE_POLYNOMIAL_SCALE_DIFFUSE_EXCLUSION.md)
covers positive-strength limits. The differential argument below
also covers zero strength, including a zero-strength core
collision.

The same rank calculation proves a simpler one-scale result:
appending \(r\) copies of \((M,M+d)\) and \(dr\) copies of
\((M+1,M)\) to any fixed ordinary core gives all-degree
coprimality above a computable height, uniformly over \(r>0\).
Here the positive and zero-strength conclusions follow from
[the affine Ax theorem](MULTIVARIATE_AX_AFFINE_COLLISION_BRANCHES.md).
Both statements extend over every characteristic-zero field.

## Each of the two pairs spans all divisor coordinates

Strip smaller-end factorials and subtract the count
at coefficient index zero at each order. For the first
pair, with multiplicities \(1,d\) per repetition, write
\(S_\ell(k)=1_{\{k\ge\ell\}}\) for \(1\le k\le j\).
Let \(B_\ell\) be the values of a putative annihilating
linear functional on this basis. At normalized offset
\(q=t/j\in(0,1)\), \(t=1,\ldots,j-1\), the divisor column
is \(S_{\lceil t/d\rceil}-dS_{j-t+1}\). The long row
starts sampling at \(\lceil t/d\rceil\), and the short
decreasing row stops at \(j-t\), giving the minus sign
after the index-zero subtraction.

At \(1<q=t/j<d\), only the long row contributes.
For each \(\ell=\lfloor j/d\rfloor+1,\ldots,j\)
one can choose a \(t>j\), nonzero modulo \(j\), with
\(\lceil t/d\rceil=\ell\). In the first block take \(t=j+1\);
in every later block of length \(d\ge2\) there is an
integer not divisible by \(j\). This produces \(S_\ell\).
In particular \(B_j=0\), and \(t=1\) in the previous
columns gives \(B_1=dB_j=0\). For \(k\ge2\) set
\(t=j-k+1\). The annihilation equation becomes
\[
                    B_{\lceil(j-k+1)/d\rceil}=dB_k .
\]
Iteration of the index map reaches \(1\) or cycles among the
indices \(2,\ldots,j\). A cycle of length \(\lambda\)
forces \(B_k=d^\lambda B_k\), which gives \(B_k=0\).
Thus every functional annihilating the columns is zero.
The first pair alone spans the full raw block at order \(j\).

The other pair is the column reversal of the same pair, with
short multiplicity \(d\) and long multiplicity one, and
therefore also spans. On normalized raw vectors the reversal
acts by \(v(k)\mapsto v(j-k)-v(j)\), an invertible linear
map. All selected offsets \(q\) are nonintegral, so they
are exclusive to their order across adjacent coprime orders.

In the two-scale setting their normalized factor polynomials
are \(M^a+q\) and \(M^b+K+q'\), respectively. These two
sets of divisors are disjoint: at a joint root, raising the
two equations to \(b/\gcd(a,b)\) and \(a/\gcd(a,b)\),
respectively, would equate real numbers of opposite signs.
Within a scale distinct offsets give polynomials differing
by a nonzero constant. Hence each scale alone gives full
raw row rank at the pair of orders. The
[linear normalization map](MULTIVARIATE_OVERLAPPING_AND_SEPARATED_PAIR_DIVISOR_RANK.md#raw-counts-and-the-normalization-map)
that sets the lower leading multiplier to one is surjective,
so the two matrices after normalization each have full row
rank separately.

## Excluding every nonzero repetition branch

Write \(u=1/M\), \(s=n+N-1\) and \(e=b-a>0\).
Let \(L=(L_1,\ldots,L_h)\) be the logarithms of the distinct
normalized irreducible polynomial divisors across both scales.
For a monic divisor \(F(M)\), use
\(L_i(u)=\log(u^{\deg F} F(1/u))\) with value zero at \(u=0\).
These are analytic germs and their exponentials are rational.
They are \(\mathbb Q\)-linearly independent modulo constants,
by distinct divisor supports, and
[Ax's power-series theorem](https://doi.org/10.2307/1970774)
gives their algebraic independence over \(\mathbb C(u)\)
as in equation (2) of the affine Ax proof. Rational
constant multipliers do not enter these normalized logarithms.

Let \(D_1,D_2\) be the two constant divisor matrices after
normalization, with zero columns in the other scale block.
Both have rank \(s\). Apart from the fixed core and its
binomial coefficients, the exact multiplier exponents are
the vector
\[
                   Y=r\big(D_1+u^{-e}D_2\big)L.             \tag{2}
\]
Here the matrices include the fixed integer count weights
within each pair and the common normalization factor.
Their column blocks are disjoint. The polynomial template
normalization makes \(Y\) analytic in \(u,\tau\), with
\(\tau=2(d^2+d)r u^a\). In fact every low-scale logarithm
starts at order at least \(a\), and every upper-scale
logarithm at order at least \(b\).
The signed gap is zero at each scale separately, so the
limiting common Gaussian clock has zero signed drift.

As in the [affine branch proof](MULTIVARIATE_AX_AFFINE_COLLISION_BRANCHES.md#a-contact-branch-has-constant-rational-repetition-parameter),
a sequence of accumulating common roots at a finite
limiting strength gives an analytic branch after finite
ramification, with \(u\) nonconstant and
\(r=\tau/(2(d^2+d)u^a)\) meromorphic. Both \(L\) and \(Y\)
are analytic on the parameter disk. The corank-two
coefficient locus gives
\(\operatorname{trdeg}_{\mathbb C}\mathbb C(\exp(Y))\le s-2\).

Suppose first \(r\notin\mathbb C(u)\) on this branch.
A nonzero rational linear relation modulo constants among
\(L,Y\) must have the form \(a_0L+b_0Y=c_0\) with
\(b_0\ne0\). In particular
\[
                   r=\frac{c_0-a_0L}
                         {b_0(D_1+u^{-e}D_2)L}.
\]
These numerator and denominator polynomials in \(L\) are
nonproportional over \(\mathbb C(u)\) and are coprime.
Two relations express the same ratio. Unique factorization
and algebraic independence over \(\mathbb C(u)\) imply
that both their numerator and denominator are proportional
with a scalar in \(\mathbb C(u)\). Full row rank then
makes their \(b_0\)-vectors proportional by that scalar,
which must be rational by comparing any nonzero coordinate.
Consequently there is at most one relation over \(\mathbb Q\).
With zero relations, Ax gives a lower transcendence bound
\(h+s+1\) for \(L,Y,\exp L,\exp Y\). The upper bound comes
from adjoining \(u,L,r,\exp Y\) and is
\(1+h+1+(s-2)=h+s\). If there is one relation, \(r\) is
already in the field generated by \(u,L\), reducing the
upper bound to \(h+s-1\), while Ax's lower bound is \(h+s\).
Both possibilities are impossible.

Now suppose \(r\in\mathbb C(u)^\times\).
A relation \(a_0L+b_0Y=c_0\) has \(c_0=0\) and
\(a_0+r b_0(D_1+u^{-e}D_2)=0\) by that same
algebraic independence. Both \(b_0 D_1\) and \(b_0 D_2\)
are nonzero for \(b_0\ne0\), since each matrix has full
row rank. Their supports are disjoint. Thus comparison in
the first block forces \(r\) to be constant, and comparison
in the second forces \(r u^{-e}\) to be constant. These
cannot both hold for nonzero \(r\) and nonconstant \(u\).
There are therefore no relations. Ax's lower bound is
\(h+s+1\), whereas adjoining \(u,L,\exp Y\) gives the
upper bound \(h+s-1\). This excludes all nonzero \(r\)
on a contact branch. The sole possible branch is \(r=0\).

## Effectivity of the exclusion

The [polynomial divisor theorem](MULTIVARIATE_POLYNOMIAL_SCALE_DIFFUSE_EXCLUSION.md)
supplies the effective neighborhoods at positive strengths
and the compact comparison away from the exceptional
points. At zero strength the only prospective common line
is \(\tau=0\). Compute the minimum vanishing order along
this line of the Sylvester resultant and the corank-two
minors. Testing a partial derivative there is exact: the
exponentials in (2) equal one, so the derivative is a
polynomial in \(L\) over \(\mathbb Q(u)\). Its vanishing is
equivalent to vanishing of its coefficients by algebraic
independence. Finite search terminates since the reference
resultant at \(u=0\) is not identically zero. Divide by
the corresponding common power of \(\tau\). The reduced
system has an isolated or empty common zero germ, since
every analytic common branch of the unreduced system was
on the removed line. The
[effective contour isolation argument](MULTIVARIATE_AX_AFFINE_COLLISION_BRANCHES.md#effective-isolation-at-zero-strength)
gives a computable punctured exclusion neighborhood.
The removed line has \(r=0\), while the ordinary inputs
have \(r>0\). The absolute curvature and terminal cutoffs
complete the uniform all-degree result.
