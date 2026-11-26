# Normalized joint modular words and exact recovery

Two products of genuine abc letters recover the first two
symmetric functions of a normalized pair on a full-divisor
chart. Their combined precision can force exact recovery
across several locks. The same records have an upper count
that charges repeated prime incidences across the three
original roles.

## Exact pair recovery from two modular products

**NJ1.** Fix positive primitive rows \(A+B=C\le X\).
Choose a nonempty subset of the three roles \(W\).
For each chosen role fix \(D_W\mid W\) on every row and
an integer precision \(\ell_W\ge1\). The fixed divisors
on different roles are pairwise coprime by original
primitivity. In the orientation with \(W\) as target,
write
\[
                    U+V=W,\qquad
          \theta_W= U/(W/D_W) \pmod {D_W^{\ell_W}},
\]
where \(U,V\) are signed versions of the other two
positive letters with prescribed signs. Assume the
\(D_W\)-supported cofactor of \(W\) has been removed
in full, so \(W/D_W\) is prime to \(D_W\).
No new values are used in the second summand of this
signed expression. For a pair of original rows record
\[
  P_W=\theta_{W,1}\theta_{W,2}\pmod{D_W^{\ell_W}},
  \qquad
  Q_W=\prod_{i=1}^2(-V_i/U_i)\pmod{D_W^{\ell_W+1}} .
\]
When
\[
                  \prod_W D_W^{\,4+2\ell_W}>X^8 ,
\]
the simultaneous records recover the unordered pair,
including repeated rows. If all counted original
\(A\)'s satisfy \(A\le X^a\) for fixed \(0<a\le1\),
replace \(X^8\) by \(X^{4(1+a)}\). The same improvement
uses \(B\) in place of \(A\). The condition refers
to the numerical values of the divisors. For example,
at precision two on all chosen roles, a combined
weight strictly exceeding one suffices at full height.

To prove this, note the exact local identity
\(-V/U=1-D_W/\theta_W\), interpreting the ratios as
rational numbers with unit denominators at \(D_W\).
Write \(S_W=\theta_{W,1}^{-1}+\theta_{W,2}^{-1}\).
Then
\[
               Q_W=1-D_W S_W+D_W^2/P_W .
\]
The two records determine \(S_W\bmod D_W^{\ell_W}\),
and hence determine the monic degree-two polynomial
of the \(\theta_W\)'s, whose linear coefficient is
\(-P_W S_W\). They also determine the polynomial of
the inverses, since those are units.

More generally, suppose that for an \(m\)-tuple of
rows the monic polynomial of its normalized residues,
of degree \(m\), is given modulo \(D_W^{\ell_W}\) at
every chosen role. If two such tuples have equal
registers, their mutual cross-determinant product
\[
             T=\prod_{i,j=1}^m
                  \bigl(A_i C'_j-A'_j C_i\bigr)
\]
is divisible by \(\prod_W D_W^{m^2+m\ell_W}\).
Here is the local check for \(p^e\parallel D_W\).
Each determinant has a factor \(p^e\) coming
from the role-\(W\) lock. The normalized determinant
after division by \(D_W\) is a unit multiple of the
difference of normalized slopes at this prime.
For \(W=C\), choose \(U=A\) and use
\((A_i(C'_j/D_W)-A'_j(C_i/D_W))\); choosing \(B\)
instead only changes a sign. For \(W=A\), the
slope after this division is \((A_i/D_W)/C_i\).
Its residue equals \(\theta_{W,i}^{-1}\) up to a
fixed sign, since the free numerator is congruent
to either \(C_i\) or \(-C_i\) at the lock. To
obtain this with precision \(\ell_W\) when
\(\ell_W>1\), choose \(C\) itself as the free
numerator for the locks on \(A\) and \(B\).
For \(W=B\), use \(B/C\) with the sign change of
the original determinant. With these prescribed
numerators the slope identities hold exactly as
rational identities, so there is no loss of
precision. For \(W=C\) use \(A\) as above.
These are the orientations used in the higher
precision statement.

The product of the differences of the slopes is
the resultant of their two monic polynomials.
If these polynomials are congruent modulo
\(p^{e\ell_W}\), then, at every root of the
first polynomial, the second has valuation at
least \(e\ell_W\). All roots here are actual
local units. Hence the resultant has valuation
at least \(me\ell_W\) in addition to the
\(m^2 e\) baseline. The divisors at the three
roles are coprime, giving the stated product.
If this product exceeds \(X^{2m^2}\), the
cross-determinant product vanishes, so at
least one original primitive ratio \(A/C\)
is common to the two tuples. Remove its
monic factor at each lock. Cancellation
modulo prime powers is valid because the
factor is monic. The same numerical test
for the smaller degree follows from the
original one: the quotient of the log
divisor bound by the square of the degree
is \(\sum_W(1+\ell_W/m)\log D_W\) and
increases as \(m\) decreases. Induction
recovers the unordered tuple. An \(A\) or \(B\)
height cap bounds each determinant by
\(X^{1+a}\), proving the refined threshold.

For NJ1 the initial numerator can be chosen
as either of the free letters. At precision
greater than one, for the two locks on
summands use the prescribed target numerator
\(C\) or convert the polynomial by the exact
translation \(\theta\mapsto\pm\theta\pm D_W\)
from \(W=U+V\). The conversion is monic
polynomial substitution modulo the specified
power. The criterion is unconditional on
radical sizes and all three caps may be
imposed on its rows. In particular if a
chosen collection of pairs has exponent
\(2d-\rho\), its joint product catalogue has
at least that exponent. A derived smaller
upper catalogue at the required precision
would thus contradict the extraction. Separate
marginal bounds at a shallower precision give
no such comparison.

## Repeated primes in normalized signed words

**NJ2.** Fix signed nonzero integers \(e_1,\ldots,e_m\).
Use genuine primitive rows of polynomial height at
\(n=X^L\), with all three original radical caps. For
one fixed full lock \(D\mid W_i\) of fixed positive
numerical weight, put \(T_i=W_i/D\) and choose either
signed free numerator \(U_i\) in each slot. The
other free letter is congruent to \(-U_i\) modulo
\(D\), so changing this choice changes the
normalized residue word only by a fixed sign.
Off any fixed polynomial-height prime pool and
off \(D\), write
\[
 b_p=\#\{ \text{occurrences of }p
                         \text{ among the selected }U_i,T_i\},
   \qquad \lambda_p=\sum_i e_i(v_p(U_i)-v_p(T_i)).
\]
Set \(H=\sum_{b_p>0}b_p\log_n p\) and
\(\Psi=\sum_{b_p>0}(b_p-1-\mathbf1_{\lambda_p=0})\log_n p\).
A prescribed list of \(R\) outputs for
\(\prod_i(U_i/T_i)^{e_i}\bmod D\) gives count
\[
             O\bigl(n^{(H-\Psi)/2+\eta}
                          R^{1/2} X^\varepsilon\bigr)
\]
on any sufficiently fine finite mesh of the shown
masses. The right side can use the relative
expansion price of RX1 at polynomial auxiliary
size. The selected complete values \(U_i,T_i\)
determine the original row, because the divisor
\(D\) is prescribed and \(W_i=D T_i\).
Each prime in a canceled slot has at least two
incidences by primitivity and nonzero row
weights. Apply RX7 to the indicated \(2m\)
signed slots, using its joint valuation
accounting. This prices a visible prime at
half once and pays a canceled prime once.
Their raw radical incidence bound comes from
the original selected roles. All original
caps can still be tested during recovery.

There is an additive test at one extra digit.
For one chosen numerator per slot,
\[
     \prod_i(-V_i/U_i)^{e_i}\equiv
         1-D\sum_i e_i\,\theta_i^{-1}\pmod {D^2}.
\]
Binomial expansion modulo \(D^2\) gives this
also for negative weights. Thus a list for the
reciprocal sum modulo \(D\) gives a list of
the same size for a true signed word on the
two free original letters modulo \(D^2\).
Apply RX7 to this latter word. The chosen
free operands recover each primitive row,
and both of their original radical caps
enter the estimate. For a pair this additive
test and the normalized product are the two
records in NJ1 at precision one. Its identity
at arbitrary fixed precision uses the
quadratic relation in NJ1 rather than
truncation.

[Signed word count and repeated incidences (RX7)](ABC_RESIDUE_EXPANSION_AND_JOINT_RAY_COUNTS.md#signed-residue-words-and-diagonal-energy-on-the-cells)

### The three-role overlap test

**NJ3.** Suppose all three roles have such fixed full
locks. For one chosen family of \(m\)-tuples, suppose
the normalized signed word on each role lies in
a list of exponent at most \(b\). Define \(\Omega\)
to be the normalized radical mass, counted once
per prime outside the prescribed pool and all the
locks, of primes occurring on at least two rows.
Omit from \(\Omega\) a prime occurring in exactly
two rows on the same original role if its two
weighted full valuations cancel exactly. If the
original radical caps are \(n^{1/3+\tau}\) per
letter, then for every \(\eta>0\) the chosen word
sector has count
\[
           O\bigl(n^{m/3+m\tau-\Omega/24+
                              b/2+\eta} X^\varepsilon\bigr)
\]
with \(\Omega\) replaced by any lower bound for
it. In particular for positive weights every
repeated outside prime is included.

For the proof choose the target role uniformly
from the three and each free numerator
independently uniformly from the two available
letters. A prime occurring in at least three
rows is included at least three times with
probability at least \(1/8\) and contributes
at least one unit to \(\Psi\) in that event.
For two occurrences in the same role, their
noncanceling weighted sum contributes whenever
both are included, an event of probability
at least \(1/4\). For exactly two occurrences
on different roles, write their nonzero
weighted orders \(x,y\). If \(x+y\ne0\),
choose the third role as target and include
both occurrences as numerators, an event of
probability \(1/12\). When \(x+y=0\), choose
one of these roles as target and include the
other occurrence as numerator, giving the
nonzero difference \(x-y\) or \(y-x\) with
at least that probability. Thus
\(\mathbb E\Psi\ge\Omega/12\).
There are finitely many choices at fixed \(m\).
Assign tuples to a maximizing choice and
apply NJ2 in each, with raw incidence
\(H\le2m(1/3+\tau)\). The rotations use
all three original role caps; a different
pair of roles can win for different tuples.

If only two fixed deep locks are available,
the same proof works after also omitting
two-row occurrences on the two locked roles
when the two weighted orders are equal.
Those two target choices both read the
difference of the weighted orders; the
omitted third lock would read their sum.
A pair in any other two-row case is caught by
at least one of the available choices with
probability at least \(1/12\). With the two
locks, the bound therefore still applies
at occurrence if the stated adjusted overlap
is present. For either number of locks, a
near-full sector with loss \(\rho\) from
\(m/3\) and output maximum price \(b\)
requires \(\Omega\le24\rho+12b+o(1)\).
This identifies which overlaps can give a
strict gain and which modular tests miss them.


## Exact complementary words at high fixed precision

**NJ4.** Suppose two different roles carry fixed full
divisors \(D_1,D_2\) of fixed positive numerical
weight. Divide the two locked letters by the remaining
complete letter, absorbing the actual signs so that the
two ratios \(x_i,y_i\) satisfy \(x_i+y_i=1\).
The ratios may be signed when the two locked roles
are on different sides of the original equation.
At \(D_1\) use \(x_i/D_1\), and at \(D_2\) use
\(y_i/D_2\), both local units. For fixed nonzero
weights \(e_i\), prescribe the two respective modular
products at precisions \(\ell_1,\ell_2\).

Let \(E=\sum_i|e_i|\). If
\(D_1^{\ell_1}>X^{2E}\) and
\(D_2^{\ell_2}>X^{2E}\), equality of these records
for two tuples forces equality of both exact
complementary products \(\prod x_i^{e_i}\) and
\(\prod y_i^{e_i}\). Indeed the common power
of the fixed divisor cancels in the ratio of two
outputs, which has rational height at most
\(X^{2E}\). The corresponding congruence to one
then gives equality by rational cross multiplication.
In the ordinary \(A,B\) lock orientation, both
ratios lie in \((0,1)\). In the other orientations
one recovers the two positive complementary words
for \(A/C,B/C\) by taking quotients and inverting
the two exact outputs, with the known fixed signs.

Consequently, for any fixed \(j\), the joint
normalized product moment of order \(j\) (equal
positive slot weights on each side) on this chart
is diagonal at the exponent scale at saturation.
Take both precisions large enough for the doubled
word comparison before choosing the radical
budget threshold. The exact-word energy bound
(SWP1) on genuine original rows applies to the
nontrivial active patterns. The remaining patterns
are the diagonal pairs of multisets. This
holds even at the minimum complete dictionary,
where the ordinary residue surplus is zero.
For an occurrence chart it gives an exact
alternative to the shallower joint records of
NJ1. The required fine-precision joint
catalogue still has to be extracted.
[Exact complementary words (SWP1)](ABC_CROSS_TARGET_COMPLEMENTARY_WORDS.md#exact-signed-complementary-words-with-all-three-caps)

## Tight full locks and the third register

**NJ5.** Fix a chart of original primitive rows
with full \(D_B\mid B\) and \(D_C\mid C\). Assume
\(D_C\) has fixed positive numerical height and
fullness at least \(\lfloor\kappa/L\rfloor\)
for a fixed \(\kappa>0\). Let \(r_A,r_B\) bound
the normalized radical of \(A\) and of \(B/D_B\)
in base \(n\). Let \(\alpha_B\) cap the exponent
of the complete \(B\) dictionary.
If there is a third full lock \(D_A\mid A\),
put \(Q=D_A D_B D_C\). Write \(S_A\) for
the actual \(B\bmod D_A\) register list and
\(q_A=\log_n|S_A-S_A|\). Let
\[
 \Theta=\min\left\{\log_n(1+X/(D_BD_C)),\
                  q_A+\log_n(1+X/Q)\right\}.
\]
In the absence of the third lock use the first
entry only. Then the row-count exponent on
the chart is at most
\[
             \frac{r_A}{2}+
                \frac12\max\{\alpha_B,\Theta+r_B/2\}
                     +\eta ,
\]
with a factor \(X^\varepsilon\) and arbitrary
fixed \(\eta,\varepsilon>0\) at sufficiently
small fixed budgets. The first radical cap
may omit prescribed support as long as it
can be fixed at subpower cost.

Let \(a_r,b_r\) be the dictionary sizes in the
shared classes \(r=A\equiv-B\bmod D_C\).
The absolute residue count for any \(k\) classes
bounds the sum of the \(k\) largest \(a_r\)'s at
\(n^{r_A/2+\eta'}\sqrt{k}X^{\varepsilon'}\).
Hence the \(k\)-th largest is bounded by
\(n^{r_A/2+\eta'}X^{\varepsilon'}/\sqrt{k}\)
and summing squares costs \(n^{r_A+o(1)}\)
with arbitrary fixed allowances. For distinct
\(B\)-values in a shared class,
\[
             t=(B_1-B_2)/(D_BD_C)
\]
is a nonzero integer of height at most
\(X/(D_BD_C)\). Its residue at a third lock
is in the fixed unit multiple
\((D_BD_C)^{-1}(S_A-S_A)\bmod D_A\).
Thus there are \(O(n^\Theta)\) possible
integer shifts across the chart, counting
actual differences.

For each fixed \(t\), write
\(g=\gcd(B_1/D_B,B_2/D_B)\).
It is prime to \(D_C\) by primitivity,
hence divides \(t\), with subpower choices.
The two reduced cofactors differ by
\(D_C t/g\). Apply FS1 with its deep divisor
\(D_C\) to the genuine coprime reduced terms,
whose radical caps remain bounded by
\(n^{r_B}\). This costs
\(n^{r_B/2+\eta'}X^{\varepsilon'}\)
per shift. Add the diagonal
\(n^{\alpha_B}\) to the second-moment count.
By original addition each pair of complete
\(A,B\) gives at most one row, so
\(\#\mathcal E\le\sum_r a_r b_r\).
Cauchy–Schwarz gives the asserted exponent.

For instance at balanced raw caps
\(r_A,r_B\le1/3+o(1)\), if the combined
actual lock catalogue costs \(t_0\), the
chart estimate gives a strict saving whenever
\(t_0<u/2\) and \(t_0+\Theta/2<1/12\)
by fixed margins, with \(u=1/3-\alpha_B\).
In a near-full chart with \(u>0\), the
same-residue distinct pairs alone have exponent
\(s+2u=1/3\) by RB4, so necessarily
\(\Theta+r_B/2\ge1/3-o(1)\).
If the full three-lock product is at least
\(X\), this forces \(q_A\ge1/6-o(1)\) at
the default cap. In particular a third list
with minimally growing difference set,
\(q_A=s+o(1)\), is excluded when
\(\alpha<1/4\) with fixed margin, after the
actual chart price is included. At occurrence
the diagonal branch of the chart estimate
reaches \(1/3\).
[Fixed-shift half count (FS1)](ABC_FULL_SHIFT_DICTIONARIES_AND_CODEGREES.md#a-fixed-shift-count-without-a-radical-cap-on-the-shift)
· [Actual residue-cell sizes (RB4)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md#the-actual-residue-cell-sizes)
