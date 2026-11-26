# Ordinary amplification of the paired and evaluated norm budgets

Fix adjacent orders \(n,N=n+1\) with \(n\ge4\), any finite
ordinary core \(C\), and a signed tail gap \(h\).
For the evaluated budget fix a nonzero algebraic integer
\(\alpha\) in the [common rational chart](MULTIVARIATE_EVALUATED_FACTOR_NORM_HEIGHT.md).
It can be the reference root used in that bound.
The common argument change from the original variable is
\(z\mapsto z/(q_0 R_{n,1}^{\rm tail})\), where \(q_0>0\) is fixed
and rational. Here \(R_{n,1}^{\rm tail}\) is the first
binomial-stripped factorial ratio of the tail.

**Theorem.** Start from any ordinary finite profile with core
\(C\) and a nonempty shared unequal tail of signed gap \(h\).
Preserve all its rows and the exact maximum relative tail gap
\(\eta\). In every real coefficient neighborhood of the starting
profile, there are ordinary extensions at the same signed gap
with all-degree angular separation and the following two
properties at a fixed sufficiently large prime \(p\):

1. For every \(2\le r\le n\), the quantity
   \((-v_p(S_+-S_-))_+\) from the
   [paired factorial norm budget](MULTIVARIATE_PAIRED_FACTORIAL_DIFFERENCE_NORM_HEIGHT.md#rational-scales-and-the-local-test)
   is arbitrarily large, using the entire tail outside \(C\).
   Both parent polygons at \(p\) have only unit roots. Their
   common quadratic allowance survives any prescribed finite
   depth of the truncated face tests. Consequently the indicated
   local exponent counts in the paired denominator cap at that
   precision, as the candidate norm has valuation zero.
2. In the common rational chart the primitive parents \(F,G\)
   satisfy
   \[
      \min\bigl(v_{\mathfrak p}(F(\alpha)),
                v_{\mathfrak p}(G(\alpha))\bigr)\ge L
   \]
   at the same degree-one prime \(\mathfrak p\) over \(p\), for any
   prescribed \(L\). This is an actual divisor of their common
   evaluation ideal. It also counts in the evaluated norm cap as
   long as the finite tests retain their two-unit-root allowance.

Here the added total curvature can be made arbitrarily small.
The prime can be chosen outside any fixed finite set. The two
depth parameters are independent. First fix the depth used for
the local quadratic and the evaluation, then grow the paired
denominator exponents with a fixed ordinary calibration.

At a simple interior [sign-changing diffuse reference](MULTIVARIATE_DIFFUSE_INTERIOR_CONTACT_PERSISTENCE.md),
the starting profile can retain a fixed sufficiently weak
anchor and have arbitrarily close roots in the common physical
normalization. Thus the two budgets can be inflated
simultaneously in any chosen neighborhood of this reference,
with the anchor fixing \(\eta\) as the local depths grow.
The construction is also valid for a sequence of such anchors
with \(\eta\to0\).

## The two literal factorial marks

Choose \(p>nN\) outside the primes of the fixed rational data.
For the evaluated assertion also take \(p\nmid
T(0)\operatorname{disc}(T)\), with \(T\) the monic minimal
polynomial of \(\alpha\), and require \(T\bmod p\) to have a
linear factor. There are arbitrarily large such primes.
Indeed a nonconstant integer polynomial has values divisible
by infinitely many primes. To see this, if their prime set
were finite with product \(P\), choose \(A\) with \(T(A)\ne0\).
The nonconstant integral polynomial
\(T(A+T(A)Px)/T(A)\) is congruent to one modulo \(P\).
Its nonzero values would have no prime divisors, contradicting
their unbounded size. Removing the fixed bad primes gives the
claimed simple unit root.

Let \(K\) be a fixed precision at least as large as any
tested \(p\)-depth. Choose an increasing and a decreasing
unit-gap row, with smaller exponents \(a,a'\), so that
\[
                   Na+1\equiv p^B\pmod{p^{B+K}},\qquad
             Na'+N-1\equiv p^B\pmod{p^{B+K}},
             \qquad B\ge K .
\]
These congruences have integer solutions. At the upper order
the increasing row contributes \(+B\) to the valuation of
the first adjacent factorial ratio and zero to the rest.
The decreasing row contributes \(-B\) at the second adjacent
ratio and zero elsewhere. Indeed its inverse suffix first
contains the marked offset \(N-1\) at coefficient index two.
At the lower order all factors of both ramps are units:
their offsets after multiplication by \(N\) are nonzero
integers of absolute value less than \(nN<p\).
The total coefficient valuation is therefore \(B\) at upper
index one and zero elsewhere in both orders. All unit residues
of the row ratios modulo \(p^K\) are independent of \(B\),
because the marked units are one and the bases modulo \(p^K\)
are fixed.

Keep the two depths while adding the same arbitrarily large
multiple of \(p^{B+K}\) to both bases. This preserves the
data. Their signed gaps cancel, their heights become
asymptotically equal, and the product of their factorial
ratio multipliers tends to one in each of the two orders.
This last assertion follows termwise from the leading powers
of the common translated height, which cancel between the
positive and negative gaps. The same comparison holds at
every other fixed order. Both rows become more diffuse than
the original row attaining \(\eta\).

## A fixed calibration with a shared local quadratic

In \(\mathbb Z_p\) lift a simple residual root of \(T\) to a
unit \(u_1\), and choose a second unit \(u_2\) with different
residue. Put
\[
              J(z)=(1-z/u_1)(1-z/u_2)=1+j_1z+j_2z^2 .
\]
Choose \(u_2\) also so that \(j_1\) is a unit.
The restriction removes at most one more residue.
In the lower raw constant-normalized polynomial prescribe
the first coefficient to be \(n r_{n,1}^C/q_0\), a unit by
the choice of \(p\). Take the lower target \(JA\) with
\(A(0)=1\), degree \(n-2\), and its linear coefficient
chosen accordingly. Take the upper target \(JD_B\) with
\(D_B(0)=1\), degree \(N-2\), and linear coefficient
\(-j_1+p^B\). Their free cofactor coefficients may be
chosen so that all lower coefficients and all upper
coefficients except the linear one are units. One can also
require the two cofactors to be coprime to \(J\) and to
each other modulo \(p\). Indeed each forbidden coefficient
vanishing or evaluation at a previously fixed nonzero
root excludes a proper affine condition on the free
coefficients. For the coefficient constraints, the
cofactor coefficient at the same index is free when the
index is at most the cofactor degree. At the last two
indices its leading coefficient enters with the respective
nonzero multipliers \(j_1,j_2\). For an evaluation at a
nonzero root the free quadratic coefficient has nonzero
multiplier. The counts of constraints depend only on
\(n,N\), so a sufficiently large \(p\) permits these choices.
The lower cofactor has a free quadratic coefficient since
\(n\ge4\).

All target unit data modulo \(p^K\) are independent of \(B\),
including the upper linear unit after its \(p^B\) is removed.
Divide these unit data by the corresponding data of the
two marked rows. The
[real and finite-prime approximation theorem](MULTIVARIATE_REAL_LOCAL_COEFFICIENT_APPROXIMATION.md)
now gives a fixed ordinary calibration with zero local
coefficient valuations and those divided unit data.
The calibration preserves the original rows and signed gap,
and can be chosen all-degree angular separated. It may also
match additional prescribed prime data at arbitrary finite
precisions. Adding the marked rows recovers the targets at
the prescribed precision and gives the stated coefficient
valuation spike. The all-degree separation persists when
both added rows are sufficiently weak. At the finitely
many orders below the absolute cutoff this follows from
openness of disjoint root arguments. At higher orders it
follows from the [absolute cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).

The lower first coefficient makes the common rational
argument change the identity modulo \(p^K\). Passing to
primitive parents only multiplies the local arrays by
units, because their coefficients are integral at \(p\)
with unit constant coefficients. The two targets have
unit endpoints and share exactly the two distinct residual
roots of \(J\). They share \(J\) also modulo \(p^K\), so
the truncated Sylvester maps retain width at least two
through this precision. At the degree-one prime given by
\(u_1\), both actual parent evaluations have valuation at
least \(K\), since both locally agree with multiples of
\(J\). This proves the assertion on the evaluated budget.

## All paired orders from this same spike

Before applying the orientation sign, the upper-grid
multipliers on the first and second adjacent logs are
\[
       (-1)^{r-1}N^{r-1},
       \qquad (-1)^{r-2}N^{r-1}
                         \bigl(r-1-\mathbf 1_{\{r=n\}}\bigr).
\]
The last extreme window meets the second slot precisely
when \(r=n\), and its first alternating weight accounts
for the correction of one. At the lower order both
marked rows have valuation zero. The fixed core, the
row prefactors from the lower first ratio, and the
\(N/n\) twist give no additional \(B\)-dependent
contribution at \(p\).

For even \(r\), the increasing row places its \(+B\)
contribution in \(S_-\), and the decreasing row places
its \(-B\) contribution in \(S_+\). The coefficients
of \(B\) in \(v_p(S_+),v_p(S_-)\) are therefore
\[
  -N^{r-1}(r-1-\mathbf 1_{\{r=n\}}),\qquad +N^{r-1}.
\]
For odd \(r\) the increasing row flips the signs, while
the decreasing row uses the unflipped upper weights.
Both contributions are then in \(S_-\), and the
coefficients of \(B\) are
\[
                    0,\qquad
       -N^{r-1}(r-2-\mathbf 1_{\{r=n\}}).
\]
The negative weights are strictly negative for \(n\ge4\).
Contributions from the fixed ordinary calibration are
constant in \(B\). Thus for all sufficiently large \(B\)
the two valuations differ, \(v_p(S_+-S_-)\) is their
minimum, and the asserted paired cap grows linearly in
\(B\) at \(p\) for every \(r\).

The evaluated common ideal is an exact arithmetic object
of the constructed ordinary arrays. The paired cap and
the two-slot evaluated cap use the hypothetical
common-factor allowances from the specified finite tests.
These allowances can disappear at deeper precision or at
other primes. In particular the constructed separated
parents have true exact common-root width zero at \(p\).
The result therefore obstructs uniform bounds on these
concrete budgets from fixed real tail controls and
preassigned finite precision. An upper bound specific to
an exact common factor, adaptive local depth or prime
choice, or a method selecting different tail subsets,
requires additional information.
