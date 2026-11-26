# Exclusion of positive opposed factorial boundary contacts

Fix an ordinary two-monomial core, two positive integer gaps \(d,e\),
and adjacent orders \(n,N=n+1\) with \(n\ge4\). Append \(r_+\)
copies of \((L,L+d)\) and \(r_-\) copies of \((Y+e,Y)\), where
\(L,Y\) are ordinary heights and the signed gap \(h=dr_+-er_-\)
is fixed. Use the
[positive Gaussian contacts](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md#the-finite-set-of-physical-reference-contacts)
of the core at twist \(h\).

**Theorem.** Each positive boundary contact, at which one of these
two appended types has zero limiting curvature share, has a
neighborhood free of exact ordinary collisions for all sufficiently
large \(L,Y\). This holds at multiple and ramified contacts. The
cutoff is uniform over the integer heights and counts in that
neighborhood, with the core, gaps, and signed gap fixed. The proof
does not compute the cutoff.

The [pole-degree theorem](MULTIVARIATE_OPPOSED_BOUNDARY_POLE_DEGREES.md)
shows that any surviving branch would have distinct height scales
\(L\) and \(Y\) with \(L<Y\ll L^2\) asymptotically. On the
algebraic height curve of such a branch, the first correction from the
taller row acts on the two Gaussian clocks before any finite-height
correction from the shorter row. Treating its transcendental count
constant as an independent analytic parameter then contradicts
isolation of the physical Gaussian contact.

## The algebraic small parameter

Assume that ordinary collisions accumulate at a positive boundary
contact as \(L,Y\to\infty\). Exchange the columns and row labels
if needed so that the increasing \(L\)-row contributes the positive
limiting share, and the decreasing \(Y\)-row has zero limiting
share. Write \(g=\gcd(d,e)\), \(c_+=e/g,\ c_-=d/g\), and
\(r_\pm=c_\pm r+s_\pm\) for fixed integer offsets \(s_\pm\);
if \(g\nmid h\) there are no ordinary points. Let
\[
 \sigma=r/L,\qquad
 \sigma_*=\frac{T_*}{c_+d^2}>0 .
\]
The boundary bias is \(\lambda_*/T_*=1/d\), so
\(\sigma_*\) is a nonzero rational multiple of
\(\log w_*\), where \(w_*=\exp(-T_* /(2nN))\) is positive
algebraic and different from one. Thus \(\sigma_*\) is
transcendental by
[Hermite–Lindemann](https://isa-afp.org/entries/Hermite_Lindemann.html).

By the
[opposed auxiliary-polynomial theorem](MULTIVARIATE_OPPOSED_AUXILIARY_POLYNOMIAL_RATIONAL_BIAS.md),
an accumulating sequence at rational bias has a subsequence on a local
branch whose heights lie on a rational algebraic curve. The
[pole-degree theorem](MULTIVARIATE_OPPOSED_BOUNDARY_POLE_DEGREES.md)
then gives \(Y\sim C L^\beta\) with \(C>0\) and \(1<\beta<2\).
Put \(t=1/L\). The ratio \(\xi=L/Y\) is algebraic over
\(\mathbb Q(t)\) on this curve. Choose the branch at \(t=0\)
that contains the sequence. By Puiseux expansion, after taking
\(t=s^M\) for a positive integer \(M\), its ratio is holomorphic
in \(s\) with algebraic Taylor coefficients. The growth restriction
gives
\[
 \xi(s)=\kappa s^b+O(s^{b+1}),\qquad
 \kappa>0,\qquad b=M(\beta-1)\in\{1,\ldots,M-1\}.
                                                               \tag{1}
\]
On ordinary points use the positive real \(s=L^{-1/M}\).
Only this expansion of the height ratio is needed. No algebraicity
assumption on the count along the branch is made.

## The two independent Gaussian directions

Let \(a_{j,k}=[z^k]Q_j^C(z)/Q_j^C(0)\) be the fixed
core coefficients after constant normalization, \(j=n,N\).
For \(0\le k\le j\) define the finite log sums
\[
 U_{j,k}(x)=\sum_{a=1}^{kd}\log(1+(a/j)x),\qquad
 V_{j,k}(x)=-
       \sum_{a=(j-k)e+1}^{je}\log(1+(a/j)x).
\]
Take all logarithms at one. They have rational Taylor coefficients.
Strip the factorials at index zero in each row. The remaining
increasing factors have principal powers \((jL)^{kd}\)
and the decreasing factors have powers \((jY)^{-ke}\);
their residual logarithms are \(U_{j,k}(1/L)\) and
\(V_{j,k}(1/Y)\). The common argument change by
\(L^{-dr_+}Y^{er_-}\) leaves exactly \(j^{hk}\) from those
principal powers. Thus the original moment pair is represented
in one common variable by the coefficients
\[
             a_{j,k}j^{hk}
                    \exp\big(\sigma F_{j,k}(s)+G_{j,k}(s)\big),
                                                               \tag{2}
\]
where
\[
 \begin{split}
 F_{j,k}(s)
   &=\frac{c_+ U_{j,k}(t)+c_- V_{j,k}(t\,\xi(s))}{t},\\
 G_{j,k}(s)
   &=s_+ U_{j,k}(t)+s_- V_{j,k}(t\,\xi(s)),\qquad t=s^M .
 \end{split}
\]
The apparent division is removable, and \(G(0)=0\).
By (1) and \(b<M\),
\[
 \begin{split}
 F_{j,k}(0)&=c_+ U_{j,k}'(0),\\
 B_{j,k}(s)&=\frac{F_{j,k}(s)-F_{j,k}(0)}{s^b},\\
 B_{j,k}(0)&=\kappa c_- V_{j,k}'(0).
 \end{split}                                                     \tag{3}
\]
Here \(B\) is holomorphic with algebraic Taylor coefficients.
The active row's remainder in \(F\) is \(O(t)=O(s^M)\).
The inactive row contributes its stated linear term at order \(b\)
with error \(O(t\,\xi^2)\) beyond its linear function of \(\xi\).
This proves (3) even when \(b>1\).

Now introduce independent complex variables \(\eta,\zeta\) and
replace (2) by the coefficient family
\[
 a_{j,k}j^{hk}\exp\big(\sigma_*F_{j,k}(0)\big)
       \exp\big(G_{j,k}(s)+\eta F_{j,k}(s)+\zeta B_{j,k}(s)\big).
                                                               \tag{4}
\]
The prefactor involving \(\sigma_*\) is algebraic, since
\(F_{j,k}(0)\) is rational and \(\sigma_*\) is a rational
multiple of \(\log w_*\). Therefore (4) is holomorphic with
algebraic Taylor coefficients in \((s,\eta,\zeta)\). At
\[
                       \eta=\sigma-\sigma_*,\qquad
                       \zeta=\sigma_* s^b
\]
it is the exact ordinary family (2).

At \(s=0\), (4) is precisely the Gaussian family of the original
core, with parameters \(p=\sigma_*+\eta\) and \(q=\kappa\zeta\)
in the notation of the
[exact normalization](MULTIVARIATE_OPPOSED_DIFFUSE_PUISEUX_HEIGHT_BRANCHES.md#the-exact-complex-one-parameter-germ).
To check this directly, the linear and quadratic terms in \(k\)
of the two derivatives in (3) are obtained by summing
\(a/j\) over their indicated blocks. For arbitrary \(p,q\)
these terms give shape curvature
\(T=c_+d^2p+c_-e^2q\) and relative endpoint drift
\(-\lambda/(2nN)\) with \(\lambda=c_+dp-c_-eq\).
Removing the lower-order linear term by a common argument change
therefore gives the clocks \(w=\exp[-T/(2nN)]\) and
\(v=\exp[-\lambda/(2nN)]\) at twist \(h\).
The linear map from \(\eta,\zeta\) to \(T,\lambda\) has
nonzero determinant \(-\kappa c_+c_-de(d+e)\), so these
are independent complex clock variables near the contact.

## Isolation contradicts the fixed transcendental slope

For the polynomials in (4), let \(A(s,\eta,\zeta)\) be
the resultant of their monic off-diagonal root-quotient
polynomials, and \(D(s,\eta,\zeta)\) their physical
resultant. The endpoint coefficients are nonzero near
the origin. The quotient coefficients are symmetric rational
functions of the parent coefficients with only endpoint
denominators, so \(A,D\) are holomorphic with algebraic
Taylor coefficients. An exact ordinary collision makes
both resultants zero. Indeed even-order real positivity
makes any common root nonreal, and its conjugate is another
common root, giving a shared off-diagonal quotient.

At \(s=0\) the angular eliminant is a nonzero polynomial
in \(w\) alone. At \(w=0\) the endpoint binomials have
disjoint nontrivial root quotients because \(n,N\) are
coprime. At any fixed nonzero \(w\), the physical
resultant is nonzero as a function of \(v\), since only
finitely many nonzero ratios of the fixed Gaussian roots
are possible. Common argument and scalar normalizations
multiply the physical resultant by a nonzero analytic factor
and preserve these statements about its zeros. Thus the
central simultaneous zero at \(\eta=\zeta=0\) is isolated.
Also \(A(0,\eta,0)\) is not identically zero, because
\(w\) changes nontrivially along this active direction.

Apply Weierstrass preparation to \(A\) in \(\eta\), and
divide \(D\) by the resulting monic polynomial. Take the
norm of the remainder, namely the determinant of multiplication
on the free module of that polynomial's quotient, to obtain
a holomorphic function \(E(s,\zeta)\). Vanishing of \(E\)
is equivalent to a common zero with \(A=D=0\) among the
prepared roots. Formal preparation and division over
\(\overline{\mathbb Q}\), followed by uniqueness of the
analytic preparations, show that its Taylor coefficients
are algebraic. After shrinking the discs all the prepared
roots tend to \(\eta=0\) as \((s,\zeta)\to(0,0)\).
The isolated central simultaneous zero implies that
\(E(0,\zeta)\) is not identically zero.

The accumulating ordinary sequence gives zeros of
\(E(s,\sigma_* s^b)\) with \(s\to0\), hence this composition
vanishes identically. For each \(k\) the coefficient of
\(s^k\) in the composition is the finite sum
\[
               \sum_{i+bj=k}[s^i\zeta^j]E(s,\zeta)\,\sigma_*^j.
\]
Its coefficients are algebraic and \(\sigma_*\) is
transcendental, so every coefficient in this sum must
vanish. This forces \(E\) to vanish identically, contrary
to the central isolation. The contradiction proves the
claimed neighborhood exclusion.

The result, together with the
[interior rational-bias exclusion](MULTIVARIATE_OPPOSED_AX_BRANCH_RANK.md),
leaves only interior irrational-bias contacts among positive
two-type limits. The orders below this note's range are already
covered by [terminal fourteen](MULTIVARIATE_TERMINAL_FOURTEEN_ALL_PROFILES.md),
and the [absolute cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
leaves finitely many relevant orders. No conclusion for zero
limiting strength or arbitrary growing lists of distinct tail types
is used here.
