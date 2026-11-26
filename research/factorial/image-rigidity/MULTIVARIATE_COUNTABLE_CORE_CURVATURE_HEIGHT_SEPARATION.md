# Curvature and height bounds for countable ordinary cores

Fix adjacent orders \(n,N=n+1\), \(n\ge2\), and a bounded
curvature \(S=\sum_w |b-a|^2/\max(a,b)\le K\). An ordinary
retained core \(R\) can have finitely or countably many rows,
with integer multiplicities. Equal rows are discarded.
At cutoff \(B\), let \(F_B\) contain the rows with
\(r(w)=\max(a,b)\le B\), and let \(s_B\) be the curvature
of the omitted rows. For a finite core \(F\) put
\[
           \mathcal H(F)=
               \prod_{w=(a,b)\in F}(r(w)+1)^{|b-a|}.
\]
This cost counts gap length and the logarithm of a row's height.
Multiplicities enter the product. Finite prefixes of a bounded
curvature core exist even if its row count is infinite.

Consider the endpoint-normalized Gaussian reference of the
[exact compactification](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#exact-compactified-data)
at strength \(T>0\). Disjoint root arguments imply physical
coprimality for every signed bias and total degree difference.

**Theorem.** For each positive algebraic \(T\) there are effective
constants \(c,E>0\), depending on \(n,T,K\), such that a
common root argument of the two reference polynomials forces
\[
                    s_B\ \ge\ c\,\mathcal H(F_B)^{-E}
                      \qquad(B>4K).                            \tag{1}
\]
Thus a countable core with \(\liminf_B s_B\mathcal H(F_B)^E=0\)
has disjoint reference arguments at that strength. This also
excludes multiple reference contacts.

For a fixed pair there is a computable integer \(D\ge1\) with
the following metric version. For any \(\nu>0\), the set of
positive strengths at which **some** retained core of bounded
curvature with
\[
                    \liminf_B s_B\mathcal H(F_B)^\nu=0          \tag{2}
\]
has a common reference argument has Hausdorff dimension at most
\(D/\nu\). This dimension estimate includes all curvature
bounds \(K<\infty\). In particular, when (2) holds for every
\(\nu>0\), the exceptional strengths, over all such cores
and all adjacent orders, form a set of Hausdorff dimension
zero. This set contains no positive algebraic strength.

There is also a uniform ordinary transfer at a prescribed
positive algebraic strength \(T\). Fix a curvature bound \(K\)
and an integer \(J\) at least the
[absolute ordinary terminal cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Take \(E\) at least the exponents in (1) for \(3\le N<J\).
Suppose a family of finite ordinary cores has computable cutoffs
\(B_m\to\infty\), a common curvature cap \(K\), and computable
bounds \(\mathcal H(F_{B_m})\le H_m\), \(s_{B_m}\le\epsilon_m\)
such that \(\liminf_m H_m^E\epsilon_m=0\). Then computable
\(\delta,\eta_0>0\), uniform over the family, have the following
property. Append arbitrary ordinary unequal shared rows whose
curvature \(\tau\) satisfies \(|\tau-T|\le\delta\) and whose
largest relative gap is at most \(\eta_0\). The resulting
ordinary supports satisfy SFC(2) in every degree. There is no
retained row-count bound, and neither a bound on the number of
added rows nor a signed bias restriction is needed.

## A quotient resultant with controlled integer height

Write \(q=e^{-T/(2nN)}\) and \(j=n,N\).
For a finite core, first divide its raw moment by the constant
coefficient. In any increasing row, the resulting multiplier
at index \(k\) is a product of integer factors from the span
\[
      j\min(a,b)+1,\ldots,j\max(a,b).
\]
In a decreasing row it is the reciprocal of such a product.
Let \(U_j^-\) be the product of the entire displayed spans
over decreasing rows, with multiplicities. Multiplying the
constant-normalized Gaussian reference by \(U_j^-\) gives an
integer polynomial \(B_j(q,z)\). Its coefficient at \(z^k\)
is one integer times \(q^{(nN/j)k(j-k)}\). The integer is at
most \(2^j\mathcal H(F)^{j^2}\), because each row uses at most
\(j|b-a|\) factors no greater than \(jr(w)\le(r(w)+1)^j\).
The constant and leading coefficients \(C_j,L_j>0\) are
independent of \(q\).

Let \(\mathscr Q_j(q,X)\) be the monic polynomial of ordered
off-diagonal root quotients \(r_b/r_a\) of \(B_j\), with
the roots indexed with multiplicity. The exact integer identity is
\[
     A_j(q,X)=(-1)^j
         \frac{\operatorname{Res}_z(B_j(q,z),B_j(q,Xz))}
                {(X-1)^j}
          =(L_jC_j)^j\mathscr Q_j(q,X).                          \tag{3}
\]
Indeed the resultant is
\(L_j^{2j}\prod_{a,b}(Xr_a-r_b)\).
The diagonal factors give \((X-1)^j\); the leading
coefficient of the whole product is
\((-1)^j(L_jC_j)^j\), as claimed. Division by the monic
polynomial \((X-1)^j\) takes place over \(\mathbb Z[q]\).
If \(j\) is odd, adjoin the factor \(X-1\) to both
\(A_j\) and \(\mathscr Q_j\), to account for the odd real
root rays in the angular test. Keep the same names for
these augmented polynomials and put \(\ell_j=j(j-1)+[j\text{ odd}]\).

The integer polynomial \(P_F(q)=\operatorname{Res}_X(A_n,A_N)\)
is a positive constant multiple of
\(\mathscr R_F(q)=\operatorname{Res}_X(\mathscr Q_n,\mathscr Q_N)\).
The constant is
\[
                (L_nC_n)^{n\ell_N}(L_NC_N)^{N\ell_n}.
\]
By the two Sylvester determinants, their degree in \(q\) is
at most
\[
     D=2(n\ell_N v_n+N\ell_n v_N),\qquad
                v_j=(nN/j)\lfloor j^2/4\rfloor .
\]
For \(A=2(n^3\ell_N+N^3\ell_n)\), both the height of \(P_F\)
and the positive constant are bounded above by
\(c_n\mathcal H(F)^A\) with an effective \(c_n\).
In fact the first determinant is homogeneous of total degree
\(2j\) in the original coefficients. Division in \(X\) in
(3) adds only a fixed finite integer linear combination of
its coefficients. The outer determinant has respective degrees
\(\ell_N,\ell_n\) in the two quotient arrays. These operations,
and the displayed bound on the original coefficients, give
the asserted degree and height estimates.

At \(q=0\) each root set is a positive scaling of the roots
of \(1+z^j\). The two nontrivial root-of-unity quotient sets
are disjoint by \(\gcd(n,N)=1\), and the additional root \(1\)
appears in just one of the two polynomials. Therefore
\(\mathscr R_F(0)=R_0\ne0\) with \(R_0\) depending only on
the orders. In particular \(P_F\) is nonzero for every core.

## From the omitted curvature to the resultant

Positive argument rescaling does not affect ordered root
quotients. In the endpoint normalization both endpoint coefficients
are one, so (3) shows that the quotient resultant is a universal
polynomial in the two input coefficient arrays. These arrays are
in the fixed box with coordinates between \(0\) and \(\binom jk\).
The resultant is consequently Lipschitz there with an effective
constant obtained from its partial derivatives.

For \(B>4K\) all omitted rows are shared and
\(|b-a|/r(w)<1/2\). Along the factorial chord the trigamma
is bounded by \(2/(jr(w))\), and the log contribution at
index \(k\) satisfies
\[
                 |J_{j,k}(w)|\le
                      \frac{k(j-k)}j\frac{(b-a)^2}{r(w)} .
\]
This is the second-derivative quadratic comparison from the
compactification. Summing the bounds gives \(j s_B/4\) for
the log coefficient error. On the fixed coefficient box,
the normalized coefficient errors are bounded by a constant
times \(s_B\). Thus the full reference quotient resultant and
the finite one differ by at most \(L s_B\) with effective \(L\).
This estimate uses the actual normalized arrays and does not
require simple roots.

At positive \(T\), the even reference polynomial has no real
root. Its undeformed core polynomial is nonnegative on the
real axis, even for a countably retained core by finite
truncation. The Gaussian multiplier averages positive lognormal
rescalings, and at each nonzero real argument the continuously
distributed scale misses the finite zero set almost surely.
At zero its constant coefficient is one. If nonreal roots
in the two orders have the same argument, their conjugate
quotients coincide. Hence a shared argument makes the full
reference quotient resultant zero and gives
\[
                         |\mathscr R_{F_B}(q)|\le Ls_B.         \tag{4}
\]

For algebraic \(T>0\), use the quantitative
Lindemann–Weierstrass bound for the fixed nonzero algebraic
\(\alpha=-T/(2nN)\), with \(q=e^\alpha\). Theorem 1 and the
explicit bound of Proposition 1 of
[Fischler and Rivoal, 2026](https://pmb.centre-mersenne.org/articles/10.5802/pmb.67/) ([author copy](https://rivoal.perso.math.cnrs.fr/articles/mesureexp.pdf))
give effective \(c_1,\gamma>0\) such that
\(|P(q)|\ge c_1 H(P)^{-\gamma}\) for every nonzero integer
polynomial of degree at most \(D\). The number field in that
theorem can be \(\mathbb Q\); its parameter \(d\) is the
degree of \(\alpha\), and the bound includes \(d=1\).
Combining the integer height and the scaling bound gives
\[
              |\mathscr R_F(q)|\ge c_2\mathcal H(F)^{-E},
                       \qquad E=A(\gamma+1).
\]
Together with (4) this proves (1).

## Size of the exceptional strength set

For any \(\sigma>1\) the ordinary finite multisets have
\[
       \sum_F\mathcal H(F)^{-\sigma}
        = \prod_{r\ge1}\prod_{d=1}^r
                      \bigl(1-(r+1)^{-\sigma d}\bigr)^{-2}
        <\infty .                                               \tag{5}
\]
At each maximum \(r\) and absolute gap \(d\) there are exactly
two orientations. The convergence follows by summing the
geometric bound for \(\sum_{d=1}^r(r+1)^{-\sigma d}\) and using
\(\sum_r(r+1)^{-\sigma}<\infty\).

Since \(\mathscr R_F(0)=R_0\ne0\), the small-value condition
\(|\mathscr R_F(q)|\le e\) on \(0\le q\le1\), for small \(e>0\),
is covered by at most \(D\) real intervals of lengths
\(O_n(e^{1/D})\), uniformly in \(F\). To verify this, factor
the polynomial as \(R_0\prod_{j=1}^m(1-q/\xi_j)\), \(m\le D\).
At least one factor has modulus at most \((e/|R_0|)^{1/D}\).
If this last quantity is at most \(1/2\), that root has
\(|\xi_j|\le2\), and \(q\) is within twice that quantity of
\(\xi_j\). The interval can be centered at \(\Re \xi_j\).
A constant polynomial has no small-value set here.

Consider the limsup of the sets
\[
                  |\mathscr R_F(q)|\le \mathcal H(F)^{-\nu}
\]
over finite cores with \(\mathcal H(F)\to\infty\).
For any \(s>D/\nu\), the coverings have summable \(s\)-powers
of diameters by (5), with exponent \(\sigma=\nu s/D>1\).
The limsup therefore has Hausdorff dimension at most \(D/\nu\).
If an infinite core has a contact and satisfies (2), (4) places
its \(q\) in this limsup. If the core is finite, its possible
\(q\) belong to the countable union of the zero sets of the
nonzero polynomials \(\mathscr R_F\). The same dimension bound
holds on adding this set. The change \(q=e^{-T/(2nN)}\) and
its inverse are Lipschitz on compact positive strength intervals,
which suffice by a countable exhaustion. The all-order
zero-dimensional assertion follows by taking arbitrarily large
\(\nu\) at each fixed pair, then the countable union over pairs.

## Uniform transfer for count-unbounded families

For every pair below \(J\), choose \(m\) with \(B_m>4K\) so
that \(L\epsilon_m<\tfrac12 c_2 H_m^{-E}\), making a common
choice for the finitely many pairs. The liminf hypothesis
ensures this search terminates. All the full-core Gaussian
references then have a common positive lower bound on the
moduli of their quotient resultants. The
[uniform diffuse coefficient transfer](MULTIVARIATE_UNIFORM_TIGHT_CORE_DIFFUSE_STABILITY.md#transfer-to-the-appended-ordinary-rows)
bounds the normalized log coefficient errors after adding
ordinary rows by \(j\tau\eta/4+j|\tau-T|/8\), where \(\eta\)
is the largest added relative gap. The same polynomial
Lipschitz bound gives effective \(\eta_0,\delta\) keeping
the quotient resultant nonzero under this transfer.
Orders \(1,2\) separate by real positivity and the
absolute ordinary cutoff handles \(N\ge J\).
The cutoff value is currently unevaluated.

For a concrete count-unbounded family, allow one unit-gap
row at each of any initial segment of levels, with arbitrary
orientation and maximum \(r_i\in[R_i,2R_i]\). Choose
\(R_{i+1}>2R_i\) and
\(R_{i+1}\ge 2\bigl(\prod_{k\le i}(2R_k+1)\bigr)^{i+1}\).
Then the cutoffs \(B_i=2R_i\) have
\(H_i=\prod_{k\le i}(2R_k+1)\) and
\(s_{B_i}\le2/R_{i+1}\le H_i^{-i-1}\).
The same family can contain the countably infinite complete
lists. The metric and algebraic tests here use the rapid
curvature decay relative to retained factorial height.
They leave open the corresponding claim without a decay
condition and at specific nonalgebraic exceptional strengths.
