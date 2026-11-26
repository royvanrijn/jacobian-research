# Pole degrees exclude fast growth on opposed factorial boundary branches

Fix a finite ordinary two-monomial core, adjacent orders \(n,N=n+1\)
with \(n\ge4\), and two positive integer gaps \(d,e\). Add \(r_+\) copies
of \((L,L+d)\) and \(r_-\) copies of \((Y+e,Y)\) with fixed signed gap
\(h=dr_+-er_-\). Both heights are ordinary integers. Work near a positive
boundary contact of the [opposed diffuse branch theorem](MULTIVARIATE_OPPOSED_DIFFUSE_PUISEUX_HEIGHT_BRANCHES.md),
and exchange the labels and monomial columns if necessary so that the
\(L\) row has positive limiting curvature and the \(Y\) row has zero
limiting curvature. Put \(q=r/Y,\ u=1/r\) in the notation of that theorem.

**Theorem.** If an exact collision branch accumulates at the contact
with \(L,Y\to\infty\), let \(Y\sim cL^\beta\), with \(c>0\), be its
leading growth. Then the taller row has \(e\ge2\) and
\[
             1<\beta<2,\qquad
             \beta\le \frac{e}{e-\lceil e/n\rceil}\quad(e\ge3).       \tag{1}
\]
In particular the order of \(q\) in \(u\) is strictly between zero and
one. The rank \(\nu\) of the space of multiplicative relations is
between two and \(e\). If \(\beta=b/a\) in lowest terms, then
\[
                      \nu a\le 2ne,\qquad \nu b\le 2nd.           \tag{2}
\]
If \(e=1\), every positive boundary contact is excluded,
including multiple and stationary contacts. For every \(e\), a simple
nonreal root pair transverse in angle is excluded. The reduction holds at all
adjacent orders in the stated range and uses no finite terminal search.
It does not give a numeric height cutoff for the branches it excludes.

The two input restrictions are arithmetic and algebraic. Infinite
ordinary collisions require an algebraic height curve and two independent
factorial relations. Integral points on that curve force the local ratio
of pole orders to equal the global ratio of map degrees. The exact
factorial divisors bound the latter ratio.

## From local growth to a ratio of map degrees

Set \(g=\gcd(d,e)\), \(c_+=e/g,\ c_-=d/g\), and write
\(r_+=c_+r+s_+,\ r_-=c_-r+s_-\) with fixed integer offsets. Put \(p=r/L\).
By the local branch theorem \(p\) tends to a positive limit and \(q\)
vanishes to a positive rational order \(\tau\) in \(u\), unless the
branch has \(q=0\) identically and gives no finite \(Y\). Hence
\(Y\sim c L^\beta\) with \(\beta=1+\tau>1\).

The bias on this boundary is \(1/d\), so the
[auxiliary-polynomial theorem](MULTIVARIATE_OPPOSED_AUXILIARY_POLYNOMIAL_RATIONAL_BIAS.md)
forces any accumulating branch onto a rational algebraic height curve
\(C\). Use the irreducible closure of the branch over \(\mathbb Q\).
It is absolutely irreducible: if a rational irreducible equation
split geometrically, the distinct geometric components would be
conjugate, and its rational points could only lie in their finite
intersections. There are infinitely many integer height pairs on the
branch, and both coordinates are nonconstant.

On the smooth projective model of \(C\), a point at infinity is a pole
of at least one of the functions \(L,Y\). Apply
[Alvanos, Bilu and Poulakis, Theorem 1.2](https://arxiv.org/pdf/0907.2097)
with integer ring \(\mathbb Z\). Infinitely many integer pairs imply
that the model has genus zero and that its geometric points at infinity
are either a single point or two points conjugate over \(\mathbb Q\).
In the single-point case both functions have all their poles there.
In the two-point case the two points have equal pole orders for each
function, since each function is defined over \(\mathbb Q\). Along the
branch tending to \(L,Y=\infty\), the quotient of the pole orders is
\(\beta\). Therefore
\[
                          \frac{\deg Y}{\deg L}=\beta,             \tag{3}
\]
where the degrees are the degrees of the two maps from the normalized
projective curve to \(\mathbb P^1\). This also applies if the projection
to \(L\) is ramified at infinity.

## Divisors of the factorial relations

By the [boundary rank theorem](MULTIVARIATE_OPPOSED_AX_BRANCH_RANK.md#boundary-branches-need-at-least-two-character-relations)
the space of multiplicative relations between the normalized factorial
multipliers has dimension at least two over \(\mathbb Q\). Taking first
differences in the coefficient index of each parent is an invertible
integral change of generators. A nonzero integer vector in the resulting
relation space, with entries \(m_{j,k}\) for \(j=n,N,\ 1\le k\le j\),
therefore gives \(A(L)=B(Y)\) on \(C\), where
\[
 \begin{split}
 A(X)&=\prod_{j=n,N}\prod_{k=1}^j
      \left(\prod_{a=(k-1)d+1}^{kd}(1+a/(jX))\right)^{c_+m_{j,k}},\\
 B(X)&=\prod_{j=n,N}\prod_{k=1}^j
      \left(\prod_{a=(j-k)e+1}^{(j-k+1)e}(1+a/(jX))\right)^{c_-m_{j,k}}.
 \end{split}                                                     \tag{4}
\]
The constant in the character relation is one because every multiplier
tends to one. The rational functions \(A,B\) are nonconstant: if either
were constant, the other would also be constant by their equality along
the branch, contradicting formal independence in two independent height
variables. This independence is the [factorial-ramp divisor argument](MULTIVARIATE_OPPOSED_AX_BRANCH_RANK.md#characters-of-the-factorial-ramps).
The equality of the nonconstant composite maps and (3) now give
\[
                              \beta=\frac{\deg A}{\deg B}.        \tag{5}
\]

Put \(E=\sum_{j,k}|m_{j,k}|\) and \(S=|\sum_{j,k}m_{j,k}|\).
Write \(C_L,C_Y\) for the loss of absolute multiplicities at the
nonzero finite points of the two divisors when equal offsets are
combined, before multiplying by \(c_+,c_-\), respectively. There is no
overlap within one order, since its blocks partition the factorial
ramp. Between orders, equality \(a/n=b/N\) forces the common offset to
be an integer because \(n,N\) are coprime. The cancellation at each
shared point is nonnegative, and is at most the sum of the absolute
weights of the factors at that point. The point \(X=0\) contributes
absolute multiplicity \(c_+dS\) or \(c_-eS\), and both rational
functions have value one at infinity. Consequently
\[
 \deg A=\frac{c_+}{2}\bigl(d(E+S)-C_L\bigr),\qquad
 \deg B=\frac{c_-}{2}\bigl(e(E+S)-C_Y\bigr).                        \tag{6}
\]
A block of \(e\) consecutive integers contains at most
\(\lceil e/j\rceil\) multiples of \(j\). Therefore
\(C_Y\le\lceil e/n\rceil E\). Since \(c_+d=c_-e\), equations (5) and
(6) prove the upper bound in (1) when \(e\ge3\). In this case
\(\lceil e/n\rceil<e/2\): \(e=3\) is immediate, and for \(e\ge4\) use
\(\lceil e/4\rceil\le(e+3)/4<e/2\).

There are \(e\) shared integer offsets in the taller ramp. If the
relation space has rank greater than \(e\), impose \(e\) linear
conditions by setting one coefficient incident to each shared
offset equal to zero. They leave a nonzero relation with
\(C_Y=0\), which gives \(\beta\le1\) in (5)–(6). Thus the rank is
at most \(e\). Since it is at least two, this also excludes \(e=1\).

## The equality case when the taller gap is two

For \(e=2\), (6) gives \(\beta\le2\). If equality holds, it must
hold for every nonzero relation vector. In particular, equality
forces \(S=C_L=0\) and \(C_Y=E\). The only common offsets of the
taller ramps are \(1,2\), and a block contains at most one of them.
The loss at a shared offset with weights \(v,w\) is
\(|v|+|w|-|v+w|\). Thus \(C_Y=E\) requires that every nonzero
weight lies in a block containing a shared offset and that the two
weights at each shared offset sum to zero. All relation vectors
then lie in the two-dimensional space spanned by the opposite
pairs for offsets \(1,2\). Since their space has dimension at
least two, each individual opposite pair gives an identity (4).

Use the pair for the shared offset two, in the first coefficient
block of both orders. Take the lower-order weight to be \(1\)
and the upper-order weight to be \(-1\). The active function has
\[
  A(L)=\prod_{a=1}^d
     \left(\frac{1+a/(nL)}{1+a/(NL)}\right)^{c_+}
      = 1+\frac{c_+d(d+1)}{2L}
                \left(\frac1n-\frac1N\right)+O(L^{-2}).
\]
The coefficient of \(1/L\) is positive. But \(B(Y)=1+O(Y^{-1})\),
which is \(1+o(L^{-1})\) since \(\beta>1\), contradicting
\(A(L)=B(Y)\). This excludes equality.

For example, at \(n=14\) the two nonintegral pairs remaining after
the common integer offsets cancel are \((13/14,16/15)\) and
\((27/14,29/15)\). The latter is the first coefficient block,
and its active response above has the clock difference
\(1/14-1/15=1/210\). This checks the indexing of the equality
argument at one terminal.

## A finite list of boundary powers

For any nonzero relation vector, put
\(k_L=\operatorname{ord}_\infty\log A\) and
\(k_Y=\operatorname{ord}_\infty\log B\), using the logarithm with
constant term zero at infinity. Comparing the first nonzero terms
in \(A(L)=B(Y)\) gives \(k_L=\beta k_Y\). Thus when
\(\beta=b/a\) is reduced, every nonzero \(\log B\) in the
relation space has leading index divisible by \(a\), and every
nonzero \(\log A\) has leading index divisible by \(b\).

There are \(e(n+N)-e=2ne\) distinct nonzero taller offsets,
because each order has a complete ramp and the overlap is exactly
the integer offsets \(1,\ldots,e\). The shorter ramp similarly
has \(2nd\) offsets. At the taller offsets \(x_i\), write
\(B(X)=\prod_i(1+x_i/X)^{t_i}\). The coefficient of \(X^{-k}\)
in \(\log B\) is \((-1)^{k+1}\sum_i t_ix_i^k/k\).
If the first \(2ne\) coefficients vanish, the Vandermonde matrix,
with its columns scaled by the nonzero \(x_i\), gives \(t_i=0\)
for every \(i\). Then \(B=1\), which is impossible for a nonzero
relation. Hence the map from rational relation vectors to their
first \(2ne\) tall log coefficients is injective. In its image,
the dimension of the subspace of vectors vanishing through
index \(k-1\) drops by at most one at \(k\), and it can drop only
at multiples of \(a\), the admissible leading indices. This gives
\(\nu\le\lfloor 2ne/a\rfloor\). The same argument for \(A\)
gives \(\nu\le\lfloor 2nd/b\rfloor\), proving (2).

These bounds produce a finite list of possible reduced powers
using only the gaps and an upper bound on \(n\). They also
restrict the analytic contact data: the ramification needed
for \(q\) to have order \((b-a)/a\) is at least \(a\).

## Transverse boundary contacts

At a specified simple nonreal common root, continue the lower and
upper roots analytically in \((u,T,\lambda)\), using the exact
normalized family of the branch theorem. Their ratio tends to one
at the physical contact. The real and imaginary parts of its local
logarithm are real analytic on the real parameter space. At \(u=0\),
the imaginary part is independent of \(\lambda\) because \(\lambda\)
only makes a positive dilation in the upper parent. The
\(\lambda\)-derivative of the real part is nonzero, since that
dilation factor is \(v=\exp(-\lambda/(2nN))\). Transversality means
the \(T\)-derivative of the imaginary part is nonzero. These two
derivatives give an invertible \((T,\lambda)\) Jacobian, so the
physical equality branch is analytic in \(u\) by the implicit
function theorem. Then \(q\) is analytic in \(u\) and either
vanishes identically or has integer order at least one. The
identically zero branch has no ordinary profiles with finite \(Y\),
and (1) excludes the other branches. This excludes a contact at
which every possible nonreal common root is simple and transverse
and neither parent has a multiple common real root. Indeed a simple
common real root stays real in the even parent under small real
perturbations, whereas the exact even factorial parent is strictly
positive on the real axis.

For a simple nonreal pair let \(k\) be the order of vanishing
of the angular difference as a function of \(T\) at \(u=0\).
This order is finite because a shared root argument on an
interval would give an identically zero off-diagonal
[quotient resultant](MULTIVARIATE_OPPOSED_DIFFUSE_PUISEUX_HEIGHT_BRANCHES.md#the-exact-complex-one-parameter-germ).
First solve the real log-ratio equation for \(\lambda\), using
its nonzero partial derivative as above. The imaginary equation
in \((u,T)\) then has Weierstrass degree \(k\) in \(T-T_*\).
The holomorphic continuations of its distinct roots permute
around \(u=0\); each orbit has size \(M\le k\), giving convergent
expansions in \(u^{1/M}\) for \(T,\lambda\) and hence \(q\).
A nonzero vanishing \(q\) has order at least \(1/M\). If there
is accumulation, (1) now requires
\[
          k\ge 2,\qquad
          k\ge
          \left\lceil\frac{e-\lceil e/n\rceil}{\lceil e/n\rceil}
          \right\rceil\quad(e\ge3).
\]
These are necessary contact conditions from exact algebraic
reference polynomials. In the holomorphic complexification of
the real log-ratio equations, each part is continued by
combining the selected roots and their conjugate branches.
At \(u=0\), derivatives in \(T\) can be taken by differentiating
the reference polynomials in the algebraic clock \(w\);
\(d/dT=-(w/(2nN))d/dw\). Hence the simplicity and finite-order
tests can be checked by exact algebraic operations.

For fixed core, gaps and signed gap the
[absolute terminal cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
and the finite local branch theorem apply the reduction across all
relevant orders. The remaining two-type boundary profiles must be
on genuinely ramified branches of growth (1), with at least two
independent factorial relations, the stated pole configuration,
and the additional [prime-grid restrictions](MULTIVARIATE_OPPOSED_ALGEBRAIC_CURVE_PRIME_GRID.md)
for their algebraic height curves. This argument does not exclude
all such ramified branches or treat tails with an unrestricted
number of distinct row types.
