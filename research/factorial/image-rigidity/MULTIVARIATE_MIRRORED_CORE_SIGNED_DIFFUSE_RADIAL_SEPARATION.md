# All-order radial separation for mirrored cores with signed diffuse tails

Split an ordinary two-monomial profile into a mirrored core and a
nonempty tail of unequal rows. The core is any multiset of pairs
\((a,b),(b,a)\), including private pairs and repeated pairs.
Equal rows may be removed. For each tail row put
\(d_i=b_i-a_i\) and \(M_i=\max(a_i,b_i)\), and set
\[
 \tau=\sum_{\rm tail}\frac{d_i^2}{M_i},\qquad
 \eta=\max_{\rm tail}\frac{|d_i|}{M_i},\qquad
 \Delta=\sum_{\rm tail}d_i,\qquad
 L=\sum_{\rm tail}\log(b_i/a_i).
\]
Use \(\eta\le1/2\), so the logarithms are defined.
The integer \(\Delta\) is also the total degree difference
of the full profile.

**Theorem.** For every rational \(\kappa>0\) there is a computable
absolute \(\eta_\kappa>0\), independent of the core and tail,
with the following property. At every terminal \(N=n+1\), \(n\ge2\),
for which
\[
            |2nN\Delta\log(N/n)-L|\ \ge\ \kappa\tau ,
                                                               \tag{1}
\]
the original ordinary moment pair is coprime if
\(\eta\le\eta_\kappa\). Terminal two is always separated.
In particular, with \(\Delta=0\) and \(|L|\ge\kappa\tau\),
the original ordinary profile satisfies SFC(2) in every degree.

The test is by disjoint physical root radii at the remaining
low terminals. The tolerance uses an unevaluated absolute
terminal cutoff. It imposes no upper bound on the number of
rows or their exponents and allows arbitrary curvature of the
core. It applies to the actual coefficient ratio shared by the
two moments.

## A comparison at exactly the same signed center

Let \(P_j\) be the order-\(j\) moment polynomial after making
the two endpoint coefficients one by positive rescaling.
The variable of \(P_j\) is \(z=G_j t\), where \(t\) is the
original coefficient ratio and \(G_j>0\) is the endpoint center.
Let \(Q_j\) be the same normalized polynomial for the core
alone, and define
\[
       F_j(z)=\sum_{k=0}^j [z^k] Q_j(z)\,
                  e^{-\tau k(j-k)/(2j)} z^k .
                                                               \tag{2}
\]
This is the core with the Gaussian tail symbol. Its roots
are simple and on the unit circle. In fact, the
[mirrored circle theorem](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md)
writes \(Q_j\) as a Schur–Szegő composition of simple
circle-rooted Gegenbauer factors. The Gaussian factor in
(2) is also circle-rooted. For these palindromic factors
composition preserves the circle location and simplicity,
and it preserves at least the angular mesh of each factor by
the [fixed-degree mesh argument](MULTIVARIATE_MATCHED_LAPLACE_NONLINEAR_TRANSFER.md#actual-mirrored-reference-profiles).

Normalize the core and tail separately and multiply their
binomial-stripped coefficients. The exact core cancels when
comparing (2) to \(P_j\). For the tail the
[joint curvature bound](MULTIVARIATE_JOINT_CURVATURE_PROFILE_CLOSURE.md#all-order-curvature-error-for-an-escaping-row)
bounds the second differences of its binomial-stripped log
coefficient array at order \(j\ge2\) between
\[
          \frac{\tau}{j(1+\eta/2)}\quad\text{and}\quad
          \frac{\tau}{j(1-\eta)} .
\]
The lower bound uses \(M_i^{-1}\le\eta\), since every unequal
integer row has \(|d_i|\ge1\). The difference from the Gaussian
second difference \(\tau/j\) is at most \(2\tau\eta/j\).
The log arrays agree at both endpoints. Comparison with
the quadratics of constant second difference gives a log
coefficient error at most \(j\tau\eta/4\). The normalized
binomial-stripped log coefficients are nonpositive on both
sides, so exponentiation gives
\[
                     |P_j(z)-F_j(z)|
                         \le(j\tau\eta/4)(1+|z|)^j .          \tag{3}
\]

For the relative physical center write \(C=\log(G_N/G_n)\).
The [exact center formula](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#normalized-polynomials-and-row-data)
has canceling contributions from every mirrored pair.
The [center remainder inequality](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#uniform-estimates-for-rows-of-large-size)
applied rowwise on the tail gives
\[
 \left|C-\Delta\log(N/n)+\frac{L}{2nN}\right|
        \le\sum_{\rm tail}\frac{|d_i|}{4n^3\min(a_i,b_i)^2}
        \le\frac{\eta\tau}{4n^3(1-\eta)^2}
        \le\frac{\eta\tau}{n^3}.                             \tag{4}
\]
For the middle bound use \(|d_i|/M_i^2\le
\eta d_i^2/M_i\), where \(|d_i|\ge1\).

## A mesh bound uniform over the core

Take a computable integer \(J\ge4\) at least the
[absolute ordinary eventual angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Only \(3\le N<J\) need comparison. At terminal two the
even moment is positive along the real axis, since it is
the expectation of the square of a nonzero real linear
combination of two distinct positive monomials.
Also, by the [absolute large-curvature theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md)
we can assume the total curvature is less than \(27\)
and hence \(\tau<27\).

Take a rational \(\eta_0\le1/2\) from the
[absolute diffuse theorem](MULTIVARIATE_DIFFUSE_ALL_ORDER_RESULTANT_TRANSFER.md).
If all core rows also have \(|b-a|/\max(a,b)\le\eta_0\),
that theorem handles the whole profile, provided
\(\eta\le\eta_0\). Otherwise, fix one mirrored core pair
with \(d/\max(a,b)>\eta_0\), where \(d=|b-a|\), and
write \(\ell=\min(a,b)\). Its last Gauss shift is \(h=d\)
and gives a Gegenbauer factor of parameter
\(\alpha=j\ell/d+1\). In particular
\[
                       1\le\alpha\le j/\eta_0+1.
\]
The [Gegenbauer equation mesh bound](MULTIVARIATE_MATCHED_LAPLACE_NONLINEAR_TRANSFER.md#actual-mirrored-reference-profiles)
gives this factor a circular mesh at least
\(2\pi/(j+\alpha)\), including the wraparound gap.
Composition carries the mesh to \(F_j\). Thus for
\(2\le j<J\) the distinct roots of \(F_j\) have chord gap
at least the positive rational number
\[
                        g=\min\!\left(1,
                              \frac{4}{J(2+1/\eta_0)}\right).
\]
The conversion uses the lesser circular distance and the
inequality \(\sin y\ge(2/\pi)y\) for \(0\le y\le\pi/2\).

## Linear radial stability and the all-degree criterion

Put \(B=J3^J(2/g)^{J-1}\). Take
\[
          \eta_\kappa=\min\!\left(\eta_0,\ \frac{g}{54B},
                                         \ \frac{\kappa}{32BJ^2}\right).
\]
Then \(r=B\tau\eta\le g/2\). On any circle of radius \(r\)
about a root \(\zeta\) of \(F_j\), the monic product satisfies
\[
                     |F_j(z)|\ge r(g/2)^{j-1}.
\]
The bound (3) is strictly smaller there, since \(|z|\le3/2\)
and \(B\) was chosen with room for the factor \(j/4\).
Rouché's theorem places exactly one root of \(P_j\) in each of
the disjoint open disks of radius \(r\) about the \(F_j\) roots.
Hence every root of \(P_j\) has absolute log radius at most
\(2r\). This follows from \(r\le1/2\) and
\(-\log(1-r)\le2r\).

If the original moments share a nonzero physical root,
its two normalized radii differ in logarithm by exactly \(C\).
This would require \(|C|\le4r\).
But (1) and (4) imply
\[
                    |C|\ge \frac{\kappa\tau}{4J^2},
                  \qquad
                    4r\le\frac{\kappa\tau}{8J^2},
\]
since \(B\ge1\) and \(n,N<J\). This is a contradiction.
The higher orders were covered by the absolute eventual
comparison. Rational Bézout gives the characteristic-zero
claim for ordinary factorial moments.

For an example with \(\Delta=0\), append the same number of rows
\((A,A+1)\) and \((3A+1,3A)\) to any fixed mirrored core.
The ratio \(L/\tau\) tends to \(1/2\) as \(A\to\infty\),
and the relative gaps tend to zero. By choosing the number of
copies proportional to \(A\), the total added curvature can
converge to any prescribed positive value. The shared degree
difference remains exactly zero.

At bounded total curvature the remaining nearly mirrored
profiles must therefore approach a signed center resonance
at the low terminal in question. When the total degrees agree,
that condition means \(L=o(\tau)\) along any sequence with
\(\eta\to0\) of purported collisions. The angular comparison
at those profiles still requires information beyond the radial
argument.
