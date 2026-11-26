# Newton residual degrees without prime factorization

Let \(P_n,P_N\), \(N=n+1\), be ordinary two-monomial factorial
moments of an arbitrary finite exponent profile, with \(n\ge2\).
Normalize each constant coefficient to one and write \(c_{j,k}\)
for the resulting rational coefficients. Every nonconstant rational common
factor has even degree because the even factorial moment is strictly
positive on the real line. This gives finite exact necessary tests
for each such degree on the original profile. The arithmetic uses
gcds, exact division, modular polynomial division and integer root
tests. Its coefficient support is profile-dependent. There is no
bound on row heights or multiplicities in the statements below.

## Coprime bases and the common residual widths

Choose pairwise coprime integers \(B>1\), with
\(c_{j,k}=\prod_B B^{e_{B,j,k}}\), \(e_{B,j,k}\in\mathbb Z\).
To obtain them, start with the distinct nonunit numerators and
denominators of the reduced coefficients. For distinct overlapping
bases \(A,B\) with \(g=\gcd(A,B)>1\), replace \(A,B\) by
\(g,A/g,B/g\), recording the new powers and merging duplicates and
units. The product of the current distinct bases decreases at every
step, so the procedure terminates. This is the elementary version
of the coprime-base construction. [Bernstein, Factoring into coprimes in essentially linear time](https://doi.org/10.1016/j.jalgor.2004.04.009)

For a fixed \(B\), use the lower polygons of the integer ordinates
\(e_{B,j,k}\). For each common slope \(s=a/b\) in lowest terms
(\(b>0\)), let \(l_j\) be the left endpoint of its edge in order \(j\).
Define
\[
 R^B_{j,s}(Y)=\sum_{t:\ l_j+bt\text{ on the edge}}
       \frac{c_{j,l_j+bt}}{B^{e_{B,j,l_j+bt}}}Y^t .
\]
An index is on the edge if its ordinate lies on the supporting
line. The displayed rational coefficients are units modulo \(B\).
For every \(p\mid B\) put \(v=v_p(B)\),
\(g=\gcd(v,b)\) and \(b_p=b/g\). Then the first stepped residual
gcd at actual slope \(vs\) has degree
\[
     g\deg\gcd\bigl(R^B_{n,s}\bmod p,R^B_{N,s}\bmod p\bigr).      \tag{1}
\]
Consequently its common width in the degree of the original
variable is \(b\) times the gcd degree on the right.

Indeed, the coefficient valuations are \(v e_{B,j,k}\).
Write \(B=p^v U\) and substitute \(z=p^{-va/b}x\) in a valued
extension with a compatible \(b\)-th root of \(p\). After removing
the common minimum power of \(p\), edge terms at \(k=l_j+bt\)
reduce to
\[
                   \overline U^{\,e_{B,j,l_j}}\,
                   x^{l_j} R^B_{j,s}(\overline U^{\,a}x^b).
\]
Use \(X=x^{b_p}\) and strip the unit and left monomial. The common
gcd of the two polynomials \(R^B_{j,s}(\overline U^{\,a}X^g)\)
is the same substitution into their raw gcd, up to a unit.
To see this, apply the substitution to a Bézout identity for
the raw cofactors. This proof includes \(p\mid g\).
It proves (1) and preserves multiplicity.

For clarity about the local bound, a monic factor formed from
the roots at actual slope \(vs\) after the substitution has integral
coefficients and unit constant. Its reduction divides both
unstepped edge polynomials after the initial powers of \(x\) are
stripped. Thus their gcd bounds its degree. The polygon product
rule also makes that degree a multiple of \(b_p\), since the
original rational factor has integral valuation ordinates and
\(b_p\) is the reduced denominator of \(vs\). A shared edge is
necessary at every nonzero multiplicity. [Filaseta and Rocks, Theorem 2.1](https://doi.org/10.4064/cm6761-3-2016) [Physical comparisons of truncated Newton faces](MULTIVARIATE_PHYSICAL_TRUNCATED_NEWTON_FACES.md)

## Exact residue gcd degrees on each coprime block

For a modulus \(M\mid B\) containing whole prime powers of \(B\),
one can run the polynomial Euclidean algorithm simultaneously
over all the residue fields \(\mathbb F_p\) with \(p\mid M\).
For any current leading coefficient represented by \(c\bmod M\),
put
\[
       T=\gcd\bigl(M,c^{\,1+\lfloor\log_2 M\rfloor}\bigr).
\]
This is the product of the full prime powers in \(M\) supported
on primes dividing \(c\). If \(T=1\), the coefficient is
invertible modulo \(M\) and one can take the next remainder.
If \(T=M\), drop the coefficient, which is zero at every
relevant residue field. Otherwise split the calculation into
the coprime moduli \(T,M/T\), retaining a copy of the input
at the split. Handle a zero polynomial by the usual Euclidean
termination rule. All operations now follow the corresponding
calculation at every prime in the retained modulus. Divisions
decrease a degree, leading zero removal decreases the current
degree, and proper modulus splits are finite. The algorithm can
record a monic representative of the final gcd at each block.

Refine the blocks for every common slope. The result is a
pairwise coprime partition of \(B\) into moduli \(M\) containing
whole prime powers and a vector \(w_s(M)\), constant for all
\(p\mid M\), of exact widths from (1). This is exact at the
first residual step and at every prime in the coefficient support.
A supplemental modulus supported outside the coefficient support
can be tested by the same polynomial algorithm with the original
coefficient arrays, which have horizontal polygons at such primes.

## Degree and denominator restrictions

For each block \(M\), fix a proposed degree \(d>0\). Take all
integer allocations \(\mathbf r=(r_s)\) with
\(0\le r_s\le w_s(M)\) and \(\sum_s r_s=d\). Use the raw slope
denominators \(b_s\) and set
\[
        \ell(\mathbf r)=\operatorname{lcm}_s
                   \frac{b_s}{\gcd(b_s,r_s)} .
\]
The convention \(\gcd(b_s,0)=b_s\) removes unoccupied slopes.
For every \(p\mid M\), a degree-\(d\) common factor requires
some allocation with \(\ell(\mathbf r)\mid v_p(M)\). This is
exactly the denominator-multiple restriction
\(b_s/\gcd(b_s,v_p(M))\mid r_s\), by prime-power divisibility.
An empty menu excludes degree \(d\). Otherwise, for
\(q=\gcd_{\mathbf r}\ell(\mathbf r)\), the integer \(M\)
must be an exact \(q\)-th power. This last test is necessary
even for a composite block of unknown factorization. Multiple
different menu entries can make it weaker than the individual
prime conditions. If the menu has just one distinct value of
\(\ell\), the perfect-power condition is equivalent to all
the denominator-multiple restrictions on that block.

There is a simultaneous projective denominator bound. Let \(L_H\)
be the least integral multiplier of the projective polynomial
of a monic degree-\(d\) common factor \(H\).
For each allocation, sort the list of \(d\) raw root valuations
\(t_i\), with \(-s\) repeated \(r_s\) times, and set
\[
        W(\mathbf r)=\max_{0\le k\le d}
                  \left(k\sum_i t_i-d\sum_{i\le k}t_i\right).
\]
Take \(W_M\) to be the maximum over the block's allocations.
The exact projective denominator formula is homogeneous in
the root valuations, giving \(v_p(L_H)\le v_p(M)W_M\).
Primes outside the coefficient support have all root valuations
zero. Thus for any common denominator \(h>0\) of the \(W_M\),
\[
       L_H^h\ \mid\ \prod_M M^{h W_M},
       \qquad
       \log L_H \ \le\ \sum_M W_M\log M.                         \tag{2}
\]
This includes arbitrary root multiplicities, and avoids any
need to factor the coefficient numerators or denominators.
The raw allocation menus are relaxations of the possible
allocations. Additional local exclusions can narrow the
maxima when supplied. [Projective common-factor sieve](MULTIVARIATE_DIFFUSE_PROJECTIVE_FACTOR_SIEVE.md#exact-projective-denominators)

For a proposed common quadratic one can use the same bound
on its trace denominator, since \(W(\mathbf r)\) is the
absolute difference of the two raw root valuations in that case.
If the auxiliary analytic argument for a profile has excluded
every pair above a specified terminal threshold, exclusion of
every even degree at each remaining pair by these arithmetic tests
certifies all consecutive orders of that profile. A bound that
would force a successful test across all unrestricted profiles
is still needed. The basis can grow with the profile.

## An ordinary factorial profile where the denominator condition matters

Take the two rows \((42,44),(39,40)\) at orders \(4,5\) and
\(p=43\). The factorial ramp of the first row at order \(4\)
contains the sole multiple \(172=4p\) at offset four. Its order
\(5\) ramp contains \(215=5p\) at offset five. Both ramps of
the second row are units. With constant normalization the
valuation arrays are \((0,0,1,1,1)\) and \((0,0,0,1,1,1)\).
The shared slopes are zero and \(1/3\). For the zero slope,
the face polynomials modulo \(43\) are \(1+27Y\) and
\(1+21Y+36Y^2\), which are coprime. For slope \(1/3\) they
are \(27+16Y\) and \(36+7Y\), each a unit times \(1-Y\),
giving common width three. Width alone permits two roots.
The denominator \(3\) requires a multiple of three at that
slope and excludes degree two. A degree-four common factor is
excluded by the independent nondivisibility theorem, so this
pair is coprime. [Two-monomial nondivisibility](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md)

## An optional Jacobi condition

For some common slopes, the monic raw residual gcd on a block
can have degree two. Choose any collection of these slopes
and refine the block so it is odd and their discriminants
\(\Delta_s\) are units. This uses the same saturated gcd
operation as above. For a prime \(p\) of the refined block at
which \((\Delta_s/p)=-1\), every irreducible local factor
at that slope has degree divisible by \(2b_{s,p}\).
Indeed, take one root \(z\) of that local factor.
The element
\[
             X=z^{b_{s,p}}p^{v_p(M)a_s/\gcd(b_s,v_p(M))}
\]
is a unit of \(\mathbb Q_p(z)\). Its residue solves the
common stepped residual equation. Consequently
\(\overline U^{a_s}\overline X^{\gcd(b_s,v_p(M))}\) is a root
of the nonsplit quadratic raw gcd. The residue degree of the
local field is even. Its ramification degree is a multiple
of \(b_{s,p}\) by the valuation of \(z\). Multiplying these
degrees gives the claim, also when a local factor is repeated.

For a vector of hypothetical Legendre signs \(\varepsilon\),
set \(\tau_s=2\) at the selected negative signs and
\(\tau_s=1\) elsewhere. Keep only allocations with
\(\tau_s\mid r_s\) for all \(s\), and replace \(\ell\) by
\(\operatorname{lcm}_s b_s/\gcd(b_s,r_s/\tau_s)\).
At a prime whose valuation in \(M\) is odd, its actual
sign vector must admit such an allocation with odd modified
\(\ell\). Let \(E\) be the finite set of sign vectors so
admitted. The vector \(((\Delta_s/M))_s\) of Jacobi symbols
must be in the subgroup generated by \(E\), since primes
with even valuation contribute the identity. If \(E\) is
empty, \(M\) must be a square. These tests are necessary.
They can be calculated without decomposing \(M\) into primes.
