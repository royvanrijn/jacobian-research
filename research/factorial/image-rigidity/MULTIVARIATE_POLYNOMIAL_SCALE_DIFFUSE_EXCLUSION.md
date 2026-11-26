# Exact diffuse exclusion at several polynomial exponent scales

Fix a finite ordinary two-monomial core and a finite list of
nonzero signed integer gaps \(d_i\). Use ordinary row types
\[
                  \bigl(a_i(M),a_i(M)+d_i\bigr),
\]
where each \(a_i\in\mathbb Q[M]\) has degree \(e_i\ge1\) and
positive leading coefficient \(\beta_i\). Give such a row
multiplicity
\[
                             r C_i(M)+S_i(M)
\]
where \(C_i,S_i\in\mathbb Z[M]\), \(C_i\) is zero or has
positive leading coefficient, some \(C_i\ne0\), and \(r\) is
a nonnegative integer. Admissible \(M\) are positive integers
giving ordinary endpoints and multiplicities. Impose the
polynomial signed-gap constraints
\[
           \sum_i d_i C_i(M)=0,\qquad \sum_i d_i S_i(M)=h,
\]
for a fixed integer \(h\). For nonzero \(C_i\) put
\(b_i=\deg C_i\) and \(c_i=[M^{b_i}]C_i\). Assume
\[
             E=\min_{C_i\ne0}(e_i-b_i)>0,\qquad
        \deg S_i<e_i\quad\text{whenever }S_i\ne0.
\]
Thus the common repetition parameter can contribute positive
curvature at several different exponent scales. The counts
within a single repetition can grow at different polynomial
rates. Define \(I=\{i:C_i\ne0,\ e_i-b_i=E\}\) and
\[
         A=\sum_{i\in I}\frac{c_i d_i^2}{\beta_i},\qquad
         B=\sum_{i\in I}\frac{c_i d_i}{\beta_i},\qquad
         \tau=\frac{Ar}{M^E}.
\]

**Theorem.** At each adjacent pair \(n,N=n+1\) with \(n\ge2\),
a computable divisor-rank test excludes accumulating exact
ordinary common roots at every positive Gaussian reference
strength on the signed-bias slope \(B/A\). The tests apply with
no restriction on the multiplicity of the reference roots.
If the tests pass at every paired order below the absolute
cutoff and the reference pairs at zero strength are coprime,
then for every sufficiently large admissible \(M\) the profile
satisfies SFC(2) in every degree, uniformly in \(r\).
The height threshold is computable. The threshold also holds
when the zero-strength end is covered directly by an
all-degree weak-tail theorem.

## The polynomial divisor criterion

Put \(u=M^{-1}\) and \(\Lambda=N/n\). For \(j=n,N\) and
\(0\le k\le j\) define the finite rational functions \(V_{i,j,k}\)
by
\[
 R_{i,j,k}(u)=\sum_{\ell=1}^{k d_i}
         \log\!\left(1+\frac{\ell}{j a_i(1/u)}\right),\qquad
 V_{i,j,k}(M)=
   \exp\!\left(nR_{i,j,k}(1/M)-kR_{i,n,n}(1/M)\right).
                                                               \tag{1}
\]
For negative upper endpoint \(m\) the sum means
\(-\sum_{\ell=m+1}^0\); the logarithms use the branch tending
to zero at \(u=0\). Thus \(V_{i,j,k}\in\mathbb Q(M)^\times\).
The coordinates at the two constants and the lower leading
coefficient are identically one. Use the remaining
\(q=n+N-1\) coefficient coordinates in the following matrix.
For every distinct monic irreducible polynomial divisor \(F\)
of the functions \(V_{i,j,k}\), make a column
\[
          E_F(M)_{j,k}=\sum_i C_i(M)\operatorname{ord}_F V_{i,j,k}
                                       \in\mathbb Z[M].       \tag{2}
\]
The divisor test is full row rank over \(\mathbb Q(M)\).
It only uses the template, including the polynomial main
counts \(C_i\).

There is a relative version. Use the projected corank-two
locus \(Z\) of the
[affine divisor-rank proof](MULTIVARIATE_AFFINE_DIVISOR_RANK_DIFFUSE_EXCLUSION.md#algebraic-projection-of-a-common-factor),
with the same core, twist \(h\) and Gaussian direction
\[
        f^0_{j,k}=\frac{k(k-j)}{2j}
                +\frac BA\frac{k}{2}\left(\frac1j-\frac1n\right).
\]
For any Laurent polynomial \(P(Y)=\sum_\alpha p_\alpha Y^\alpha\)
vanishing on \(Z\), group its monomials according to the
polynomial vectors \((\langle\alpha,E_F(M)\rangle)_F\).
The relative test passes if some group has nonzero sum
of the coefficients \(p_\alpha\). Computing an elimination
ideal for \(Z\) and checking its generators gives a finite
test. Indeed the proposed grouping is the grouping on the
column subtori for every sufficiently large \(M\), since
there are only finitely many nonidentical polynomial
differences between exponents of the generators. Full row
rank implies the relative test, since \(Z\) is proper.

## The analytic normalization

Set \(g_{i,j,k}(u)=R_{i,j,k}(u)-(k/n)R_{i,n,n}(u)\).
The exact normalization by the common lower endpoint scale
gives the moment pair with coefficients
\[
 ([z^k]Q_j^C)\,\Lambda^{hk\boldsymbol 1_{\{j=N\}}}
                  \exp\bigl(\tau f_{j,k}(u)+G_{j,k}(u)\bigr),
\]
where
\[
 f_{j,k}(u)=\frac{\sum_i C_i(1/u)\,g_{i,j,k}(u)}{Au^E},
 \qquad G_{j,k}(u)=\sum_i S_i(1/u)\,g_{i,j,k}(u).
\]
In fact stripping the order-\(j\) constant row factorial
contributes \((j a_i(M))^{kd_i}\exp R_{i,j,k}(u)\).
The powers of \(a_i(M)\) cancel under the common normalization
and the residual order ratio is exactly \(\Lambda^{hk}\)
in the upper order. Each \(g_i\) starts at degree at least
\(e_i\) in \(u\), so \(f,G\) are analytic with rational Taylor
coefficients, and \(G(0)=0\). For \(\ell=1,\ldots,q\) with
signed upper limit \(q\) use \(\sum_\ell\ell=q(q+1)/2\)
in (1). The active indices \(I\) then give \(f(0)=f^0\).
The exact tail curvature equals \(\tau(1+O(u))+O(u)\);
the signed bias has limit \((B/A)\tau\) on bounded
strength sets. These follow by dividing each row's
\(d_i^2\), respectively \(d_i\), by its larger exponent;
the leading terms are exactly those picked out by \(I\).

## Simultaneous polynomial specialization

The columns (2), evaluated at suitable large \(M\), can be
read from separate primes. For each divisor \(F\), choose a
different prime \(p_F\) at which its primitive integer form
has a simple root and at which the other divisors cannot
vanish at that root. Exclude the finite primes dividing
the leading coefficients, discriminants or pairwise
resultants, as well as the rational constant factors of
the functions in (1). The elementary prime divisor theorem
for a nonconstant integer polynomial provides infinitely
many choices for each \(F\). To see the latter, if a finite
prime set covered all its values, let \(c=F(m)\ne0\), and
let \(P\) be the product of those primes. Then
\(F(m+cPt)/c\) is a nonconstant integer polynomial in \(t\)
with all its values congruent to one modulo \(P\). Such
values eventually have absolute value greater than one,
contradicting that prime cover.

Lift the simple root of \(F\) modulo \(p_F\) to a base
residue modulo \(p_F^2\) at which the valuation is one.
Use the Chinese remainder theorem to obtain arbitrarily
large positive \(M\) in these residue classes. At each such
\(M\), the vectors of prime valuations of the actual
rational numbers
\[
                 \prod_i V_{i,j,k}(M)^{C_i(M)}
\]
at the chosen primes are exactly \(E_F(M)\). These samples
are evaluations of rational functions and need not themselves
give ordinary heights. They are used only inside the
analytic argument.

Fix a Laurent polynomial \(P\) passing the relative test.
For sufficiently large sampled \(M\), monomials with the
same exponential frequency on \(y\mapsto\exp(y f(1/M))\)
pair equally with every \(E_F(M)\). This follows by taking
the chosen prime valuations of their equality of
monomial values in the above positive rationals. They
therefore belong to the same one of the polynomial
groups in the relative test. Multiplicative translation
by \(\exp G(1/M)\) multiplies every monomial coefficient
by a factor tending to one. The nonzero group sum
persists, so the translated Laurent polynomial cannot
vanish identically on the exponential curve, by linear
independence of distinct exponentials.

## Exact exclusion at the reference contacts

At a positive reference contact \(\tau_*\), the common
ordinary Gaussian clock is algebraic by the nonzero
Laurent reference resultant and \(\tau_*\) is transcendental.
The [affine divisor-rank proof](MULTIVARIATE_AFFINE_DIVISOR_RANK_DIFFUSE_EXCLUSION.md#exact-elimination-at-positive-strengths)
constructs a contour eliminant for the analytic multiplier
arrays
\[
    \exp\left(\tau_* f^0+G(u)+v f(u)
                             +t(f(u)-f^0)/u\right)
\]
entry by entry. The hypotheses of that construction here
are the rational Taylor coefficients of \(f,G\), the
algebraic reference coefficient values, and the nonzero
one-variable reference resultant. At an actual collision
\(v=\tau-\tau_*,\,t=\tau_*u\), even positivity forces a
common factor of degree at least two. An accumulating
sequence forces the contour eliminant to vanish identically,
since each of its coefficients along \(t=\tau_*u\) is a
polynomial in the transcendental \(\tau_*\) over the
algebraic numbers.

Fix one sufficiently small sampled \(u\). Such an identically
vanishing eliminant supplies a contact for every real \(t\)
in a small interval, with \(v\) in a bounded disk.
Set \(y=v+t/u\) and \(X=\exp((\tau_*-t/u)/L)\), where
\(L\) clears the denominators of \(f^0\). The multiplier
tuple lies on (2) of the affine divisor-rank proof, with
coefficient multiplier \(\exp G(u)\exp(y f(u))\).
All possible \(y\)'s are bounded. At a contact the
Laurent polynomial \(P\) above vanishes, whereas its
substitution in this exponential curve is a nonzero
entire function by the specialization. Thus there are
only finitely many possible \(y\)'s. Some fixed tuple of
coefficient multipliers would admit infinitely many
different \(X\)'s. This is impossible by the endpoint
modulus separation for nonzero signed-bias slope and
the distinct even/odd large-clock root scales for zero
bias, proved there even for nonzero complex coefficients.
This gives the exclusion. Taylor coefficient search and
analytic remainder bounds in the same contour construction
make it effective.

At an ordinary failure the tail curvature is bounded by
the absolute curvature theorem, hence \(\tau\) is bounded.
The finitely many positive exceptional strengths at
the finitely many low orders can be isolated algebraically
in their exponential coordinate. Away from their certified
neighborhoods and from zero strength, exact coefficient
continuity gives uniform exclusion on compact strength
sets. The zero-strength tests or an applicable weak-tail
theorem handle the neighborhood of zero. The absolute
angular cutoff completes the all-degree claim.

## An exact example with two distinct positive curvature scales

Take the four common-slope balanced types
\[
       (z,z+d),\quad(z+1,z),\quad
       (z+d+1+d,z+d+1),\quad(z+d+1,z+d+2),\qquad d\ge2
\]
at \(z=M\), each with count \(r\). Add the same four offset
types at \(z=M^2+K\), each with count \(rM\), where \(K\ge2\)
is fixed. Both groups have signed gap zero and contribute
positive limiting curvature for \(r/M\) tending to a
positive constant. The quadratic-scale linear-factor
arguments become quadratic polynomials in \(M\) with
positive constant term, and have no real zero. Thus
they do not occur among the linear divisors of the first
group. The [reflected-window rank proof](MULTIVARIATE_AFFINE_DIVISOR_RANK_DIFFUSE_EXCLUSION.md#an-infinite-family-beyond-the-single-mark-test)
for the four affine types shows that the first group alone
supplies full rank at every paired order. Its witness
columns are unchanged by the quadratic group. This yields
all-degree coprimality for every sufficiently large \(M\)
and every positive \(r\), by itself or with any fixed
ordinary core already satisfying SFC(2). For the empty
core the zero-strength limit uses the absolute weak
curvature cap.
