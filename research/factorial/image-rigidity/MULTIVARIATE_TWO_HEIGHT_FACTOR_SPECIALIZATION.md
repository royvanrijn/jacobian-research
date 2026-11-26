# Exact finite height reduction for two independent unit-gap rows

Fix any ordinary core with factorial moment coefficients
\(c_{j,k}=[z^k]Q_j^C(z)>0\). Append \(r\) copies of
\((A,A+1)\) and \(s\) copies of \((B,B+1)\), with \(r,s\)
positive integers and \(A,B\) independent nonnegative integers.
There is no bound on either height. The two columns can be
exchanged before applying the result.

Label two distinct coprime orders by \(n,N\ge2\).
The order \(N\) is the one that degenerates in the proof.
Either order may be given this label.
For \(1\le h<n\) let \(w=n-h+1\) and put
\[
 \begin{split}
 \mathcal L_h&=\frac1w\left[
 (-1)^s\frac{c_{n,h-2}}{c_{n,h-1}(h-1-n/N)^r}
 -\frac{c_{n,n-1}}{c_{n,n}[n(N-1)/N]^r(w-1)^s}\right],\\
 \mathcal U_h&=-\frac{c_{N,N-1}}
                {N c_{N,N}(N-1)^r\,[N(w-1)/n]^s}.
 \end{split}                                                 \tag{1}
\]
For \(h=1\) omit the first term in brackets.

**Theorem.** If \(\mathcal L_h\ne\mathcal U_h\) whenever
\(\gcd(w,N)>1\), the degree-at-least-two common-factor locus in the full-degree
chart is finite. For an ordinary pair with at least one even
order, a common root therefore occurs at only finitely many
\((A,B)\). These pairs can be
computed by exact algebraic elimination.

The hypothesis is empty when \(N>n\) and \(N\) is prime.
It also holds automatically when both \(r,s\) are even
and \(N\) is odd, with arbitrary fixed core. When
\(N>n\), all fixed rows are increasing or equal, and at
least one multiplicity is odd, the test is automatic
after exchanging the two appended types if necessary.
Thus for a comparable fixed core and at least one even
order, the finite reduction holds for every choice of
the two multiplicities.

For consecutive orders this gives an all-degree finite decision
for the two-height family when the tests hold for all terminals below
the [absolute cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
The finite coefficient tests and the elimination involve just
the core and the two multiplicities. They require no sampled
upper bound on either height.
 
## Algebraic parameter space

Strip the factorials of the smaller endpoints of the added
rows. At ordinary \(A,B\) the resulting polynomials are
\[
 F_j(A,B,z)=\sum_{k=0}^j c_{j,k}(jA+1)_k^r (jB+1)_k^s z^k,
 \qquad j=n,N,                                               \tag{2}
\]
where \((U)_k=\prod_{i=0}^{k-1}(U+i)\).
Use (2) throughout the complex parameter plane. On the
open chart where both leading coefficients are nonzero,
the degree-at-least-two condition is the vanishing of the
Sylvester minors of size \(n+N-1\). If one order is even, its ordinary
moment is positive on the real line, as the expectation
of an even power of a nonzero real binomial under product
exponential measure. A common ordinary root therefore
comes with its distinct complex conjugate and belongs
to this locus.

A positive-dimensional component on this chart cannot have
\(A\) constant. Otherwise specialize \(B=-1/N\).
The upper polynomial becomes its nonzero constant and the
lower retains its leading coefficient. In the homogeneous
variable pair there is no common root, including at infinity.
The resultant restriction to that vertical line would thus
be nonzero. The case of a two-dimensional component is
excluded by the same specialization with \(A\) generic.

Suppose instead that there is a curve component on which
\(A\) is nonconstant. There are points of the curve in
the full-degree chart with \(A\to A_0=-1/N\).
Indeed a plane irreducible curve with nonconstant \(A\)
projects to all \(A\) outside a finite set, by the
leading coefficient of its defining polynomial in \(B\).
Intersecting with the complement of the chart removes
only finitely many points of this component. Pass to a
subsequence with \(B\to B_0\in\mathbb P^1(\mathbb C)\).

If \(B_0=\infty\), rescale the variable by \(z=B^{-s}x\).
At each order \(B^{-sk}(jB+1)_k^s\to j^{sk}\).
Thus the lower leading coefficient is nonzero and the
upper limit is its nonzero constant. This is impossible
by homogeneous resultant continuity. If \(B_0\) is finite,
the upper again limits to its nonzero constant, so the
lower leading coefficient must vanish at \(B_0\).
It follows that \(B_0=-h/n\) for \(1\le h\le n\).

At \(h=n\) the lower limit has degree \(n-1\).
It has exactly one root escaping to infinity, counted
with multiplicity, by continuity of projective roots.
Since all upper roots escape, at most one is common.
Therefore \(1\le h<n\). Now \(N B_0=-Nh/n\) is nonintegral.
Every upper coefficient at \(k>0\) vanishes at order
exactly \(r\) in \(A-A_0\), with multiplier a unit at the
limiting \(B_0\). In the lower polynomial, the coefficients
for \(k\ge h\) vanish at order exactly \(s\) in \(B-B_0\),
and the lower limit has degree \(h-1\).

## The escaping root ratios

Take a choice of \(t\) with \(t^w=(B-B_0)^s\).
Rescaling by \(z=x/t\) gives uniform polynomial expansions
on bounded \(x\)-sets,
\[
 t^{h-1} F_n(A,B,x/t)
 = x^{h-1}(\ell_0+\ell_1 x^w)
        +t(\ell_-x^{h-2}+\ell_+x^{n-1})+O(t^2).
\]
Here the four coefficients are the exact functions of \(A,B\)
from the indicated degrees in (2). Those at \(n,n-1\)
are divided by \((B-B_0)^s\) and all are continuous at
\((A_0,B_0)\), with \(\ell_0,\ell_1\) nonzero there.
Omit the term involving \(\ell_-\) when \(h=1\).
The error is uniform since all other degree differences
in this scaled polynomial are at least two.
Each of the \(w\) nonzero binomial roots is simple and
gives one escaping root. The remaining \(h-1\) roots
are bounded, as follows from the unscaled lower limit.

The analogous upper scaling uses a choice of
\(v\) satisfying \(v^N=(A-A_0)^r\), with all \(N\)
roots at scales bounded above and below times \(1/|v|\).
Since at least one root is common, \(v/t\) stays bounded
above and away from zero. After another subsequence
\(V=(A-A_0)^r/t^N\) tends to a nonzero finite number.
In the lower scale,
\[
 F_N(A,B,x/t)=c_{N,0}+V u_1 x^N
                              +t V u_+x^{N-1}+O(t^2).
\]
The exact functions \(u_1,u_+\) arise by division of the
two corresponding coefficients by \((A-A_0)^r\).
They are nonzero at the limiting point. Their continuous
variation and that of the lower coefficients need not
be \(O(t)\).

At a nonzero root of the lower leading binomial, the
first correction is the translation
\[
              t\,\frac{\ell_-/\ell_0-\ell_+/\ell_1}{w}.
\]
At a root of the upper leading binomial the correction
is the translation \(-t\,u_+/(Nu_1)\).
Both statements follow by evaluating the order-one term
and dividing by the derivative at the simple leading root.
Using the exact varying binomial coefficients first, the
errors are \(O(t^2)\) on every isolated root disk.
Only the limits of the translation coefficients matter
after division by \(t\). Evaluating the rising factors
at \(A_0=-1/N,B_0=-h/n\) gives exactly the limits (1).
For the lower left coefficient ratio the final factor
from \(B_0\) is \(-1\); for the lower right one it is
\(n-h=w-1\). In the upper right ratio omitting the last rising \(A\)
factor gives \(1/(N-1)^r\) and the
last \(B\) factor is \(N(w-1)/n\).

If two distinct roots are common, their limits in \(x\)
are distinct nonzero roots \(\lambda,\mu\) of both
leading binomials. The ratio of the two moving leading
roots is exactly \(\lambda/\mu\) in each order,
because changes of the leading coefficients apply one
common scaling within each binomial. Thus
\((\lambda/\mu)^w=(\lambda/\mu)^N=1\), with
\(\lambda/\mu\ne1\). This requires \(\gcd(w,N)>1\).
For the ratio of the two *actual* roots, its difference
from \(\lambda/\mu\), divided by \(t\), tends to
\[
       \mathcal L_h(\mu-\lambda)/\mu^2
       \quad\text{and}\quad
       \mathcal U_h(\mu-\lambda)/\mu^2
\]
in the respective orders. Since their ratios agree,
\(\mathcal L_h=\mathcal U_h\). Every candidate root along
the chosen sequence is simple by the isolated binomial
limits. The hypothesis excludes the curve.

Hence the algebraic locus has dimension zero in the
full-degree chart. Saturating the Sylvester minor ideal
by the product of the leading coefficients and solving
the resulting zero-dimensional system gives the finite
list. Alternatively eliminate all but one height at a
time. Each elimination ideal contains a nonzero univariate
polynomial, and its Cauchy bound gives a finite upper
cutoff in that height. Checking the resulting integer box
requires only the original factorial arrays.

## Automatic tests and unbounded counts

Write \(\rho_j=c_{j,j-1}/(j c_{j,j})\).
If all fixed core rows are increasing or equal, \(\rho_N\le
\rho_n\). For each increasing fixed row \((a,b)\) its
contribution is
\[
 \frac{((j-1)b+a)!}{(jb)!}
            =\prod_{i=0}^{b-a-1}(jb-i)^{-1},
\]
which decreases in \(j\). If \(s\) is odd, the two terms
of \(\mathcal L_h\) have the same negative sign, with
only the high term present at \(h=1\). The ratio of
the magnitude of the high term to \(|\mathcal U_h|\) is
\[
                    \frac{\rho_n}{\rho_N}\frac n w
                                      \left(\frac N n\right)^{r+s}>1 .
\]
If \(r\) is odd, exchange \(A,B\) and \(r,s\) to obtain
the same conclusion.

When both multiplicities are even, label the odd order by
\(N\). The cases \(w=n\) and \(w=2\) have \(\gcd(w,N)=1\).
All tests still required have \(k=h-1\ge1,\ w=n-k\ge3\).
Write the terms in the brackets of \(\mathcal L_h\) as
\(P-Q\), both \(P,Q>0\). The order-\(n\) core factorial
coefficients after removing the binomials are log-convex,
by the log-convexity of each row factorial. Hence
\[
  \frac{c_{n,k-1}/c_{n,k}}{c_{n,n-1}/c_{n,n}}
          \ \ge\ \frac{k}{n(w+1)},\qquad
  \frac P Q\ \ge\ \frac{k}{n}\frac{(w-1)^s}{w+1}
               \left(\frac{n(N-1)}{|Nk-n|}\right)^r .        \tag{3}
\]
If \(Nk>n\), the ratio in parentheses exceeds \(n/k\)
since \(k<n\). Also \((w-1)^s\ge w+1\) for \(s\ge2,w\ge3\).
Thus \(P/Q>1\).

If \(Nk<n\), then \(w\ge 2n/3\) and the ratio in
parentheses is strictly greater than \(N-1\ge2\).
The bound \((w-1)^2/(w+1)\ge n/4\) follows from
\((w-1)^2/(w+1)\ge3w/8\) for \(w\ge4\). For \(w=3\)
one has \(n=4\), and the bound is equality. Now (3)
again gives \(P/Q>1\), since \(k\ge1,r\ge2\).
Equality \(Nk=n\) is impossible by coprimality.
In both cases \(\mathcal L_h>0>\mathcal U_h\).

For a fixed comparable core the multiplicities need not
be prescribed in advance. The same applies to a fixed
mixed core if both multiplicities are even.
At a possible ordinary collision the
[absolute degree-imbalance bound](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md)
requires \(|\Delta_C+r+s|<146J-73\), where \(\Delta_C\)
is the signed gap of the core and \(J\) is a certified
integer absolute cutoff. Only finitely many pairs
\(r,s\) need to be used in the preceding elimination.
Thus all-degree SFC(2) is decidable for the entire
two-height comparable family, and for the mixed-core
family with both multiplicities even. Every potential
failure has a computable finite height cutoff. This does not assume bounds on the
heights before the calculation.

The core can also vary in either family, over all
allowed cores of height at most a prescribed integer \(B_0\ge1\).
Discard its equal rows. At a failure the
[absolute curvature bound](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
forces fewer than \(12B_0\) unequal core rows, since each has
curvature at least \(1/B_0\). There are finitely many such cores,
so taking the maximum of the calculated cutoffs gives a
finite decision for this enlarged family. Both appended
types still have independently unbounded heights and
unrestricted multiplicities in the comparable case, or
arbitrary even multiplicities in the mixed case.

The [exact replay](verify_multivariate_two_height_factor_specialization.py)
checks the two limiting translation coefficients directly
from the rising products, including even multiplicities.
