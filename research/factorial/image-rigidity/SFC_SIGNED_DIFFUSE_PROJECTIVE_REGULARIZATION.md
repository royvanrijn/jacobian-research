# Projective regularization by signed diffuse factorial rows

Let \(U_1,\ldots,U_s\) be a fixed ordinary monomial support, with
\(s\ge2\). Write its order-\(j\) factorial moment form as
\[
       F_j(x)=\sum_{|\alpha|=j}c_{j,\alpha}x^\alpha .
\]
Fix integers \(v_i>0\) and pairwise distinct integers \(w_i\) with
\(|w_i|\le v_i\). For \(\tau>0\) define the signed diffuse moment forms
\[
 F_j^\tau(x)=\sum_{|\alpha|=j}
       c_{j,\alpha}\exp\!\left(
            \frac{\tau}{2j}\sum_i(v_i\alpha_i^2+w_i\alpha_i)
                       \right)x^\alpha.                       \tag{1}
\]
At each order these are moments of the original monomials times
independent positive lognormal variables of respective log variances
\(\tau v_i/j\) and log means \(\tau w_i/(2j)\).

**Theorem.** If \(e^{-\tau/2}\) is transcendental, then any \(s\)
distinct positive orders of (1) have no common projective zero.
This includes all \(s\)-term consecutive windows and every positive
algebraic \(\tau\). The original support and its total degree
differences are arbitrary.

For any finite collection of sets of \(s\) distinct orders, the
same simultaneous projective regularity can be realized by appending
only fully shared ordinary unit-gap rows to the original support.
The extension leaves each total degree difference unchanged, and its
added pairwise curvatures can be made arbitrarily small. It works for
all coefficients over every characteristic-zero field. The threshold
for the ordinary realization depends on the chosen finite collection.

## The signed drift separates pure-term ties

In fact the first claim applies to any homogeneous forms over the
algebraic numbers with all pure coefficients nonzero. Fix a set
\(I\) of \(s\) distinct positive orders, let \(H\) be their least
common multiple, and put \(q=e^{-\tau/(2H)}\). The coefficients of
(1) are Laurent polynomials in \(q\) with algebraic coefficients.
Multiplying each form by an appropriate monomial in \(q\) makes
them polynomial in \(q\) without changing its projective zero
set for \(q\ne0\).

Their homogeneous resultant is a nonzero polynomial in \(q\).
Here is a quantitative proof. Put
\[
 \delta=\min_{i\ne h}|w_i-w_h|,\qquad
 \sigma=\min_{\substack{j,k\in I\\j\ne k}}|1/j-1/k|,\qquad
 M=\max_{j\in I}\log
           \frac{2\sum_{|\alpha|=j}|c_{j,\alpha}|}
                {\min_i |c_{j,j\mathbf e_i}|}.
\]
All three numbers are positive. We will exclude common zeros
when \(\tau>4M/(\delta\sigma)\).

Fix a nonzero coefficient vector \(x\). Retain only the nonzero
coordinates and associate to coordinate \(i\) the real line
\[
              \ell_i(t)=\frac{2}{\tau}\log|x_i|+v_i+w_i t .
\]
At \(t=1/j\), the log magnitude of a term of (1) with nonzero
coefficient and supported on these coordinates is exactly
\[
 \log|c_{j,\alpha}|+
 \frac{\tau}{2}\sum_i\alpha_i\ell_i(1/j)
 -\frac{\tau}{2j}\sum_i v_i\alpha_i(j-\alpha_i).
\]
The last term is nonpositive, and it is zero for the pure terms.
If the largest line value exceeds the second largest by at least
\(2M/\tau\), any term other than the winning pure term uses at
least one unit of a different coordinate. The total modulus of
those terms is at most
\[
    e^{(\tau/2)j\ell_{\max}-M}\sum_{|\alpha|=j}|c_{j,\alpha}|,
\]
which is at most half the modulus of that pure term.
When there is only one retained coordinate the form is already
nonzero. Thus a zero requires a tie or a gap less than
\(2M/\tau\) between the top two lines at \(t=1/j\).

Two such near-tied lines intersect within
\(2M/(\tau\delta)\) of \(1/j\). On the segment from \(1/j\)
to the intersection there must be a breakpoint of the maximum:
either the leading line ceases to lead earlier or it ties at
the intersection. The maximum of at most \(s\) lines with
distinct slopes has at most \(s-1\) breakpoints, since its
successive active slopes strictly increase. The \(s\) intervals
of radius \(2M/(\tau\delta)\) about the distinct points \(1/j\)
are disjoint by our bound on \(\tau\). Thus the \(s\) orders
cannot all vanish. This gives the claimed nonzero resultant.

If \(e^{-\tau/2}\) is transcendental, so is \(q\), since
\(q^H=e^{-\tau/2}\). A nonzero polynomial over the algebraic
numbers cannot vanish at \(q\). The
[Hermite–Lindemann theorem](https://isa-afp.org/entries/Hermite_Lindemann.html)
covers every positive algebraic \(\tau\). The argument uses one
set \(I\) at a time and applies to every such set for the same
\(\tau\).

One can also take \(v_i,w_i\) real, with the same strict signs,
distinct slopes and magnitude bounds. At each window the
resultant is then a finite linear combination of real exponentials
in \(\tau\). The same envelope argument for large \(\tau\)
shows it is nonzero as a function. If all the parameters
\(v_i,w_i,\tau\) are algebraic, group equal exponents and apply the
[Lindemann–Weierstrass linear-independence theorem](https://isa-afp.org/entries/Hermite_Lindemann.html)
to conclude that its value is nonzero. For arbitrary real
parameters, each finite window has only finitely many real
exceptional \(\tau\). A nonzero real linear combination of \(r\)
distinct real exponentials has at most \(r-1\) real zeros, by
dividing by one exponential, differentiating, and using Rolle's
theorem inductively. For complex coefficients one can apply
this argument to a nonzero real or imaginary part. Consequently
for fixed real \(v_i,w_i\) only a countable set of \(\tau>0\)
is excluded over all order sets.

## Ordinary realization with exactly balanced total degrees

For an integer \(A\ge2\), a positive tag for column \(i\) is
a new row with value \(A+1\) in column \(i\) and value \(A\)
elsewhere. A negative tag uses \(A-1\) in column \(i\).
Take respectively
\[
 K_{i,+}=\left\lfloor\frac{\tau A}{2}(v_i+w_i)\right\rfloor,
 \qquad
 K_{i,-}=\left\lfloor\frac{\tau A}{2}(v_i-w_i)\right\rfloor
\]
copies and put \(D_i=K_{i,+}-K_{i,-}\). At base \(B=A^3\),
add \(|D_i|\) tags of the opposite orientation to \(D_i\)
with the same column \(i\) tagged. Every new row is fully
shared. The signed increment above its row baseline sums to
zero in every column, so the new rows increase every total
degree by the same integer.

Divide out the row scalar \((jA)!\) or \((jB)!\) at order
\(j\). Let \(k=\alpha_i\). The contribution of a positive
tag at base \(A\) to the log coefficient is
\[
          k\log(jA)+\frac{k(k+1)}{2jA}+O(j/A^2).
\]
For a negative tag it is
\[
         -k\log(jA)+\frac{k(k-1)}{2jA}+O(j/A^2).
\]
These follow by expanding \(\log(1+h/(jA))\) for
\(1\le h\le k\), or \(-\log(1-h/(jA))\) for \(0\le h<k\).
The error bounds are uniform for \(0\le k\le j\) since
\(|h|/(jA)\le1/2\). The same formulas hold at \(B\).
The two sets of tags sum to the linear term
\(D_i k\log(A/B)\). The common positive diagonal change
\(x_i\mapsto(B/A)^{D_i}x_i\) removes it at every order.
The quadratic terms from the balancing rows at \(B\) are
\(O(jA/B)\), and the other remainders and the count rounding
are \(O(j/A)\) for fixed \(s,\tau,v_i,w_i\).
The resulting ordinary coefficients converge to those of (1)
in each order as \(A\to\infty\).

At a fixed finite collection of order sets take any positive
rational \(\tau\). The resultants of all the normalized limit
systems are nonzero. The resultants are polynomials in the
coefficients, so the actual ordinary resultants are nonzero
for every sufficiently large \(A\). The common diagonal change
and row scalars are invertible. Before normalization the
coefficients are integers, which gives the conclusion over
every characteristic-zero field. The threshold is computable
from rational enclosures of the finitely many nonzero limit
resultants and the explicit log expansion bounds.

For columns \(i\) and \(h\), only tags for \(i\) and \(h\)
contribute to their added curvature
\(\sum (\text{gap})^2/\max(\text{two exponents})\).
It tends to \(\tau(v_i+v_h)\). The balancing rows contribute
\(O(A/B)\). Thus one can achieve any prescribed positive bound
on all added pairwise curvatures by first choosing a small
positive rational \(\tau\) and then taking \(A\) large.

The signed drift in (1) gives the slopes \(w_i\) in the pure
costs. Its ordinary realization uses precisely the simultaneous
quadratic and signed affine corrections underlying the
[joint binary diffuse profile closure](MULTIVARIATE_JOINT_CURVATURE_PROFILE_CLOSURE.md).
The limit theorem treats every order at the same \(\tau\).
For the actual ordinary extension the coefficient convergence
proves simultaneous regularity on each prescribed finite
collection and supplies no order-independent threshold when
\(s>2\).
