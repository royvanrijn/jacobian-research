# Exact projective regularization by triangular ordinary row additions

Let \(U_1,\ldots,U_s\) be any fixed list of ordinary monomials,
with \(s\ge2\), and fix positive integers \(d_1,\ldots,d_{s-1}\).
For \(r=1,\ldots,s-1\), append the exponent row
\[
              (L_r+d_r\mathbf1_{i>r})_{i=1}^s .
\]
The bases \(L_r\) are independent nonnegative integers.
Put \(D=\sum_r d_r\). One may also append the fixed row
\[
                  \left(C+D-\sum_{r<i}d_r\right)_{i=1}^s
\]
for any nonnegative integer \(C\). With this last row every total
degree increases by the same amount, so the original total-degree
differences are preserved.

**Theorem.** For every set \(I\) of \(s\) distinct positive orders,
the homogeneous resultant of the normalized ordinary factorial
moment forms on these rows is a nonzero integer polynomial in the
bases. Its total degree is at most \(sD\prod_{j\in I}j\).
For any finite family of such sets, simultaneously:

- The joint exception set on the ordinary grid is the zero set of
  the product of these resultants. From any given baseline there
  is a strictly positive integral direction with only finitely many
  exceptional integer translates.
- There is an explicit product of residue classes of the bases
  which forces projective regularity at all prescribed windows.
  Each coordinate residue can be increased arbitrarily far.
- The added rows can be fully shared, and their pairwise added
  curvatures can be made arbitrarily small. The compensating row
  preserves total-degree differences in these conclusions.

The projective conclusions hold for coefficients over every
characteristic-zero field. In the binary case one can take
all windows required by the
[absolute ordinary cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL.md),
so both conclusions give ordinary extensions regular in every degree.
For higher support size the conclusion concerns the prescribed
finite collection of windows.

## The triangular specialization

Absorb the optional fixed row into the original core. At degree
\(j\), write \(b_{j,\alpha}\) for its ordinary factorial coefficient
at the monomial \(x^\alpha\) in the coefficient variables \(x_1,\ldots,x_s\).
For integer bases strip the nonzero scalars \((jL_r)!\) to obtain
the polynomial family
\[
 F_j(L,x)=\sum_{|\alpha|=j}b_{j,\alpha} x^\alpha
           \prod_{r=1}^{s-1}(jL_r+1)_{d_r\sum_{i>r}\alpha_i}.
\]
The rising factorial is one at length zero. Each \(b_{j,\alpha}\)
is a positive integer.

Order the elements of \(I\) as \(j_1>j_2>\cdots>j_s\). Specialize
over the rationals to \(L_r=-1/j_r\) for all \(1\le r<s\).
At order \(j_r\), the factor on row \(r\) kills every term
involving a variable \(x_i\) with \(i>r\). All coefficients
of the pure term \(x_r^{j_r}\) survive. The rows \(t\ge r\)
contribute one to that coefficient, while the rows \(t<r\)
use factors \(h-j_r/j_t>0\) for integers \(h\ge1\).
At order \(j_s\), the same reasoning shows that its final
pure coefficient \(x_s^{j_s}\) survives.

Thus \(F_{j_1}=0\) forces \(x_1=0\). Given
\(x_1=\cdots=x_{r-1}=0\), the equation \(F_{j_r}=0\) forces
\(x_r=0\). The final equation forces \(x_s=0\). The forms have
no common projective zero.

By the [classical homogeneous resultant and its degree formula](https://www.math.ucdavis.edu/~deloera/MISC/LA-BIBLIO/trunk/Sturmfels/cbms.pdf),
the resultant vanishes exactly at a projective common zero
and is homogeneous of degree \(\prod_{h\in I\setminus\{j\}} h\)
in the coefficients of \(F_j\). Each coefficient has total
base degree at most \(jD\), giving the degree bound. In a finite
family the product of the resultants has a nonzero top homogeneous
part. Choose a strictly positive integral direction off its zero
set, using rational density in the positive orthant and homogeneity.
The highest coefficient on every line in this direction is
nonzero, independently of the starting baseline. Only finitely
many translates on each such line can be exceptional.

## An explicit ordinary congruence construction

For one window \(I\), let \(J_I=\max I\), let \(a_*\) be the
largest fixed exponent after including the optional row, and let
\(d_*=\max_r d_r\). Choose a prime
\[
                     p_I>\max(J_I^2d_*,J_I a_*).
\]
Impose \(L_r\equiv -j_r^{-1}\pmod {p_I}\) for \(r<s\), where
the \(j_r\) are ordered separately for each window. In this
characteristic, the triangular proof is unchanged. The pure
core coefficients are units. Every numerator \(j_t h-j_r\)
in a required surviving rising factor has positive magnitude
less than \(p_I\), and each denominator \(j_t\) is a unit.
The specialized forms are therefore projectively regular over
the algebraic closure of \(\mathbb F_{p_I}\). A characteristic-zero
common projective point, if present, could be taken over the
algebraic numbers, embedded in a finite \(p_I\)-adic extension,
and scaled so its coordinates are integral with at least one
unit. Its reduction would then give a common projective point,
a contradiction.

Choose distinct primes for the finite collection and impose all
these conditions with the Chinese remainder theorem. Each
\(L_r\) has one residue class modulo the product of the primes.
Every nonnegative integer baseline in the resulting product has
a nonzero resultant at each required window. In particular
adding the same multiple of the product to every base gives
an ordinary ray with no exceptions. The identity of the monomials
is ensured by the distinct column tags. Total-degree balancing
uses only the fixed compensating row. Its baseline and all
variable baselines can be made arbitrarily large with their
gaps fixed, which makes every pairwise added curvature
\(\sum_{\rm added}(a_i-b_i)^2/\max(a_i,b_i)\) arbitrarily small.

The exact construction uses \(s-1\) variable rows and at most
one fixed row. The
[signed diffuse regularization theorem](SFC_SIGNED_DIFFUSE_PROJECTIVE_REGULARIZATION.md)
also gives simultaneous finite-window ordinary regularization
with small balanced curvature, and treats all windows in a
common transcendental limit. Here integer residues force the
finite-window conclusion with the exact factorial coefficients.
