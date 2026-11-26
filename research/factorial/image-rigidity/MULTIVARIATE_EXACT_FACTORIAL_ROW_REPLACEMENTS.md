# Exact factorial row replacements and all-degree separation

A row is its ordered exponent pair \((a,b)\). A profile is a finite
multiset of rows. At order \(m\), the contribution of a row to the
coefficient of index \(k\), apart from the binomial coefficient, is
\[
                     ((m-k)a+kb)!.
\]
For a signed integer row multiset \(\omega\) and a nonnegative
integer direction \((p,q)\ne(0,0)\), write \(\omega_{p,q}(h)\)
for the total signed multiplicity of rows with \(pa+qb=h\).

**Theorem.** Fix \(B\ge1\) and two ordinary profiles \(A,A'\).
Translate both entries of every row in each profile by the same
nonnegative integer \(L\). The following are equivalent.

1. The two translated profiles have identical raw factorial
   moment polynomials in all orders \(1,\ldots,B\) for infinitely
   many \(L\).
2. For the signed difference \(\omega=A-A'\), every projected
   count \(\omega_{p,q}(h)\) is zero when \(1\le p+q\le B\).

When they hold, the polynomial equality holds for every
nonnegative integer \(L\). Equality is preserved after appending
these rows to any common fixed core. Choose \(B\) at least one
less than the [absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Then this replacement preserves all-degree angular separation and
physical SFC(2) in both directions, for every such core and every
\(L\) for which the two columns are distinct. The numerical value
of that cutoff is unevaluated.

There are nontrivial replacements for every \(B\). In particular,
one can replace an exactly mirrored fully shared block by a
nonmirrored fully shared block with **no reflected row pairs**.
The replacement preserves the multiset of unordered row pairs,
both column degrees, and every curvature statistic depending
only on unordered rows. These examples transfer all-degree
angular separation of a private unit pair with an arbitrary
mirrored core to ordinary profiles with no reflected shared pairs
in the replaced block. The integer multiplicities can be held
fixed while the smaller exponent in every moving row tends to
infinity.

## Classification along integer translations

Fix \((p,q)\) and \(m=p+q\). The ratio of the translated factorial
contributions at \(k=q\) is exactly
\[
                      \prod_{h\ge0}(mL+h)!^{\,\omega_{p,q}(h)}.
\]
If all projected counts vanish, every factor cancels. Conversely
suppose this ratio equals one for infinitely many nonnegative
\(L\). Write \(s=\sum_h\omega_{p,q}(h)\). On factoring out
\((mL)!^s\), the remaining ratio is a rational function of \(L\)
formed from factors \(mL+h\) with \(h\ge1\). A nonzero power of
the factorial grows or decays faster than any rational function
along \(L\to\infty\), so \(s=0\). The rational factor equals one
identically, because it equals one at infinitely many \(L\).
The exponent of each distinct linear factor \(mL+h\) is
\(\sum_{j\ge h}\omega_{p,q}(j)\). Unique factorization gives
zero for all these exponents, so all projected counts at \(h>0\)
are zero. Their total \(s=0\) gives the count at zero as well.
This proves the equivalence coefficient by coefficient.

There is a corresponding normalized version. If the two degree-\(m\)
polynomials need only agree up to a positive scalar and positive
argument change, the criterion for agreement at infinitely many
translations is
\[
       m\,\omega_{m-k,k}
         =(m-k)\,\omega_{m,0}+k\,\omega_{0,m}
                           \qquad (0\le k\le m)
\]
as vectors of projected multiplicities indexed by \(h\).
To verify this, raise the quotient of the \(k\)-th coefficients
to the power \(m\) and divide by the two endpoint quotients
to the powers \(m-k\) and \(k\). Equality to one is precisely
normalized equality, since all coefficients are positive. The
preceding factorial argument then applies to this signed
combination of histograms.

Appending any fixed rows multiplies corresponding coefficients
with the binomial factor removed by the same positive integer.
Thus raw equality remains exact. For \(B\) chosen as in the
theorem, both orders of every terminal below the absolute cutoff
are covered by this equality. At and above the cutoff both
profiles have angular separation independently. The consequence
for coprimality transfers from the rationals to characteristic
zero by a rational Bezout identity.

## A reciprocity criterion using the symmetric average

Suppose the total multiplicity of each unordered pair in a block
is even. Replace its oriented counts by an equal split between
the two directions. In every degree the binomial-normalized
endpoint-normalized log coefficient of the new block is the
arithmetic average of the old log coefficient and its reflection.
This follows by additivity of log factorial coefficients and the
two endpoint corrections. Therefore, if the old block has
reciprocal endpoint-normalized polynomials through order \(B\),
the symmetric replacement has exactly the same normalized
polynomials in those orders. This equality remains true after
joining the same fixed core to both blocks, by additivity again.

In particular, once \(B\) reaches the absolute cutoff minus one,
a block satisfying this finite reciprocity criterion can be joined
to a private unit pair and any additional mirrored core. The
resulting original profile has all-degree angular separation.
The two physical centers can differ. They do not affect the
argument comparison.

## Reorienting shared rows without changing the bounded orders

Let \(\mathcal D_B\) consist of the primitive nonnegative integer
pairs \((p,q)\) with \(1\le p+q\le B\). In the Laurent polynomial
ring over the integers, put
\[
                 F(X,Y)=\prod_{(p,q)\in\mathcal D_B}
                                             (1-X^qY^{-p}).
\]
A signed row multiset with generating polynomial \(F\) has zero
projection in each of these directions. Indeed, substitution
\(X=t^p,\ Y=t^q\) kills the factor for \((p,q)\). Multiples of
the primitive directions have zero projections as well. This
is the generating-polynomial construction used for zero line
sums in discrete tomography. [Hajdu and Tijdeman, Algebraic aspects of discrete tomography](https://math.unideb.hu/sites/default/files/inline-files/disctom8.pdf)

Take a monomial multiple so the support lies in \(1\le a<b\).
For example choose \(X Y^{2+\sum_{\mathcal D_B}(p+q)}\).
Call the signed multiplicities of this shifted polynomial
\(\nu(a,b)\). They are finitely supported integers and they
are not all zero, since the Laurent polynomial ring is a domain.
Translation and the positive shift keep every vanishing
projection. Each of the two coordinate-axis projections also
vanishes, so \(\sum a\nu(a,b)=\sum b\nu(a,b)=0\).

Make the reference block with \(|\nu(a,b)|\) copies of \((a,b)\)
and \(|\nu(a,b)|\) of \((b,a)\) for each \(a<b\). Make its
replacement with \(|\nu(a,b)|+\nu(a,b)\) copies of \((a,b)\)
and \(|\nu(a,b)|-\nu(a,b)\) of \((b,a)\). Exactly one direction
is used for every unordered pair in the replacement. Its signed
difference from the reference is \(\nu-\nu^{\mathsf T}\). Transposition
interchanges the directions of projection, and \(\mathcal D_B\)
is invariant under this interchange. The two blocks therefore
have identical raw moments in the bounded orders, including
after any common translation. They have the same unordered
row multiplicities and the same column degree totals, using
the two vanishing coordinate sums. Their projected shapes
are exactly equal at the bounded orders even when raw leading
and constant coefficients and their relative physical scale
are retained.

For an all-degree separation family, append the distinguished
private units \((1,0),(0,1)\), with or without additional fixed
mirrored shared rows. The reference pair separates in every
degree by the
[private-unit mirrored interlacing theorem](MULTIVARIATE_ONE_PRIVATE_MIRRORED_LINKED_PHASE_INTERLACING.md).
For \(B\) at least the absolute cutoff minus one, the replacement
therefore separates in every degree as well. With only the
distinguished private units in the fixed core, none of its shared
rows have reflected partners. In particular the transfer permits
strict asymmetry at the level of native integer rows even though
all bounded-order polynomials are exactly reciprocal.

For fixed multiplicities the shared curvature of both translated
blocks goes to zero. To keep their curvature in any prescribed
positive band containing \(s>0\) in its interior, replicate each
block \(r_L\) times with \(r_L/L\to s/C\), where
\(C=2\sum_{a<b}|\nu(a,b)|(b-a)^2>0\).
Their common shared absolute curvature then tends to \(s\).
The equality and the regularity transfer still hold separately
for each integer \(L\).

For a finite example, the displayed product with \(B=5\) has
122 nonzero signed cells after collecting powers. Both blocks
have 328 rows. Exact integer factorial products agree through
degree five and differ in degree six, for both translations
\(L=0\) and \(L=7\).

Thus raw moment polynomials through any fixed cutoff, even with
all unordered row counts given, can fail to distinguish a mirrored
integer profile from a completely unpaired orientation. The
transfer only uses equality at the indicated orders.
[Single sufficiently high order rigidity](MULTIVARIATE_SINGLE_ORDER_ROW_RIGIDITY.md)
