# Exact deletion of complete exponent shift blocks

For an ordinary two-monomial profile, a block of exponent rows can
be deleted when its contributions to both adjacent coefficient
arrays are congruent to a constant at sufficient precision.
Complete arithmetic progressions of common exponent shifts give
an explicit modulus for this deletion. The modulus can be composite,
and the retained pair may have zero resultant at its first residue
step. The argument works at arbitrary total curvature, with the
original rows and the original physical variable.

## A boundary product for a block of rows

Fix a length \(P\ge1\), a positive integer stride \(q\), and a list
of starting pairs \((A_g,B_g)\) of nonnegative integers. The block
contains the \(P\) rows
\[
                     (A_g+iq,B_g+iq),\qquad 0\le i<P,
\]
for each \(g\). All occurrences are counted, even if two rows
coincide. Fix an order \(j\), a coefficient index \(0\le k\le j\),
and write
\[
          t=jq,\qquad d_g=B_g-A_g,\qquad
          J_g(k)=jA_g+kd_g .
\]
Define the finitely supported signed interval counts on positive
integers by
\[
      m_{j,k}(u)=\sum_g
        \left(\mathbf 1_{\,u\le J_g(k)}
                         -\mathbf 1_{\,u\le J_g(0)}\right)
\]
and let \(M_c\) be their sum on the residue class
\(c\pmod t\), \(1\le c\le t\). For \(v\ge1\) put
\[
             e_{j,k}(v)=\sum_{\substack{u>v\\u\equiv v\bmod t}}
                                               m_{j,k}(u).
\]
The contribution \(T_{j,k}\) of the block to the coefficient,
divided by its contribution \(T_{j,0}\) at index zero, is exactly
\[
 \frac{T_{j,k}}{T_{j,0}}=
  \prod_{c=1}^{t}\left(\prod_{i=0}^{P-1}(t i+c)\right)^{M_c}
             \prod_{v\ge1}(1+tP/v)^{e_{j,k}(v)} .              \tag{1}
\]
To see this, put \(W(u)=\prod_{i=0}^{P-1}(t i+u)\).
Changing a factorial argument from \(ti+J_g(0)\) to
\(ti+J_g(k)\) gives exactly the interval product of \(W(u)\)
with its signed exponent. The identity
\(W(u+t)=W(u)(1+tP/u)\) telescopes from the positive
representative \(c\) of the residue class. Collecting its
contributions gives (1) for positive and negative gaps alike.

The condition \(M_c=0\) for all \(c\) is equivalent to
\[
 \sum_g J_g(k)=\sum_g J_g(0),\qquad
 \{J_g(k)\bmod t\}_g=\{J_g(0)\bmod t\}_g                       \tag{2}
\]
as multisets. Indeed, the count at \(c\) up to an integer
\(J\ge0\) is \(\lfloor J/t\rfloor+\mathbf 1_{c\le J\bmod t}\).
Taking consecutive differences in \(c\) recovers the counts at
each nonzero residue, and the number of starts fixes the count
at zero. Summing the interval counts over all \(c\) recovers
the first equality in (2). This also proves the converse.

For stride \(q=1\), condition (2) holds at all indices of
a specified list of orders if and only if the signed gaps
sum to zero and each gap is divisible by the least common
multiple of those orders. Necessity follows by using \(k=1\)
and the fact that all \(J_g(0)\) are multiples of \(j\).
For general stride \(q\), divisibility of every gap by
\(q\) times that least common multiple and zero signed sum
is sufficient. Condition (2) also permits the residues
to be permuted between starting pairs, as opposed to being
fixed within each pair.

## A simultaneous resultant transfer

Suppose (2) holds at every index of orders \(n,N=n+1\).
For each active boundary factor in (1), meaning an index
with \(e_{j,k}(v)\ne0\), form the reduced numerator of
\(jqP/v\). Let \(D\) be the gcd of those numerators across
both orders. If the product is empty in both orders the
deletion gives equality, and no modulus is needed.

At every prime dividing \(D\), each active factor in (1) is
a unit congruent to one modulo the prime's full power in
\(D\), regardless of the sign of its integer exponent.
Therefore \(T_{j,k}/T_{j,0}\equiv1\pmod D\), in the sense
of rational reduction on the prime support of \(D\).
This test only uses active offsets. The size of a large
inactive offset has no effect on \(D\). A simpler sufficient
condition for a prescribed modulus \(D_0\) is
\(D_0 v\mid jqP\) for every active \(v\) at order \(j\).
For example, if all starting exponents are at most \(U\),
any \(P\) divisible by \(D_0\,\operatorname{lcm}(1,\ldots,NU)\)
works whenever (2) holds. The actual modulus \(D\) often
allows smaller periods.

Take any integral scalar normalizations \(F_n,F_N\) of the
retained core polynomials. After appending the block and
dividing its constant contribution out separately in each
order, the new polynomial is obtained by multiplying
coefficient \(k\) of \(F_j\) by \(T_{j,k}/T_{j,0}\).
Both arrays remain integral at the primes dividing \(D\)
and have the original formal degrees. Their Sylvester
determinant is congruent modulo \(D\) to
\(\operatorname{Res}(F_n,F_N)\). In particular, if
\[
                      D\nmid \operatorname{Res}(F_n,F_N),     \tag{3}
\]
the original enlarged pair is coprime. One need not factor
\(D\), choose a prime in advance, or require a unit
resultant at the first residue step. A rational coefficient
dilation common to both orders is also allowed if the
retained pair is integral at the tracked primes after
suitable scalar normalizations. For several disjoint
deletable blocks, use the gcd of their moduli. For the
same block repeated \(D^r\) times, one can replace \(D\)
by \(D^{r+1}\), since \(u\equiv1\pmod D\) implies
\(u^{D^r}\equiv1\pmod{D^{r+1}}\) by the binomial theorem
and induction, also for rational units on this support.

The finite comparisons may use a different deletion
partition and modulus at each adjacent pair. Let \(J\) be an
[effective absolute terminal cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
If all pairs with \(3\le N<J\) pass (3) or its exact
equality case, the enlarged profile is coprime in every
degree. The even moment excludes the order-one real root,
and the absolute theorem covers \(N\ge J\). In particular,
for a fixed all-degree regular core one can use a single
\(D_0\) exceeding the absolute values of all its nonzero
low-order integral resultants. With balanced starting gaps
divisible by all the low orders, take a length \(P\) as
in the sufficient condition for \(D_0\), independently of
the size or sign of the individual gaps.

## Independent lifts and arbitrary diffuse strength

Completeness of a block can be taken modulo a period rather
than imposed on its actual row heights. For a known modulus
\(D_0\) dividing \(D\), let \(L_0\) be a common multiple of
all nonempty literal ramp factors at the queried orders in
the reference block. These are the factors \(j\ell+u\)
with \(1\le u\le j|b-a|\) for a reference row \((a,b)\)
of smaller exponent \(\ell\). Independently add any
nonnegative multiple of \(D_0L_0\) to both entries of each
actual row. Every affected factor \(f\) is replaced by
\(f+jD_0L_0 w\), whose ratio to \(f\) is one modulo
\(D_0\) as an integer. The normalized factorial ratios,
including those for decreasing rows, consequently retain
their congruences. This permits arbitrary extra rows already
in the core and independent growth of every lifted height.

Use a nonempty balanced block with a nonzero gap. A common
lift tending to infinity makes its positive strength
\(\sum |b-a|^2/\max(a,b)\) and maximum relative gap tend
to zero. Repeating such a lifted block the requisite number
of times makes its strength converge to any prescribed
positive real number, while every individual relative gap
tends to zero. All copies remain deletable in (3), and the
tail's signed gap stays zero. The conclusion therefore does
not require a weak total tail. [Complete base-\(p\) exponent orbits](MULTIVARIATE_FACTORIAL_DIGIT_ORBITS.md) is the related first-residue construction with a nontrivial residual coefficient kernel for unbalanced gap classes.

## An ordinary example requiring the second precision

Use the core \((6,7)\) at orders \(2,3\). Its
constant-normalized coefficient arrays are
\[
            (1,26,182),\qquad (1,57,1140,7980).
\]
The absolute resultant is \(431396=31\cdot13916\) with
\(31\nmid13916\). Both leading coefficients are units at
\(31\), so the first reductions have a common factor.

Append two progressions of stride one and length \(P=31^2\)
from starting profiles \((1,7)\) and \((7,1)\), respectively.
The signed gaps \(6,-6\) are balanced and divisible by
both orders. All possible active offsets are at most \(21\),
hence their factors in (1) are one modulo \(31^2\).
The same two rows at each position are mirror pairs,
but independent lifts by a common multiple of the stated
ramp modulus can remove all such reflections.
For every allowed lift or repetition the augmented
resultant remains nonzero at precision \(31^2\), even
though its two first reductions still have a common
factor. This is a pair-level illustration of (3). The
all-degree conclusion for a chosen core uses a successful
test at each of its remaining orders.
