# Balanced complete dictionaries and common residue cells

The square-root cost of a list of absolute residues gives a
symmetry constraint on an actual balanced additive graph.
At the one-third exponent, all three complete dictionaries
have the same exponent. After fixing a cheap full divisor on one
letter, almost all rows lie in residue cells whose two moving
dictionaries have the size of a single complete-letter fibre.
The cells have the full pair-count exponent at that smaller scale.

## Regular balanced cohorts and absolute residue pricing

Put \(n=X^L\) and use original positive primitive triples
\(A+B=C\), \(X/2<C\le X\). The raw radical cap on each letter
is \(n^{1/3+o(1)}\) in the outer small-budget limit. Work with
a family of size \(n^{1/3-o(1)}\) in the height-first convention.
Cap losses from box refinements and inner \(X^\varepsilon\)
allowances are taken arbitrarily small at each fixed budget.
The same conventions apply to all normalized logarithms and
entropies in this note. Eligibility conditions, including an
exact least-repair cutoff, may be imposed on the original rows.

On a regular cohort, each complete value of \(T\) has degree
\(n^{u_T+o(1)}\), and its dictionary has size
\(n^{\alpha_T+o(1)}\). Here and below an asserted uniform
degree exponent on a cohort means upper and lower bounds
with arbitrarily small fixed tolerances. One has
\[
                 \alpha_T+u_T=\frac13,
                 \qquad 0\le u_T\le\frac16 .                  \tag{RB1}
\]
The second bound follows from (FL1). One may reach a regular
cohort by binning original complete degrees in all three roles.
With any fixed positive mesh, the number of bins is bounded
before the height limit. In a retained joint bin, the original
lower degree in a role bounds the size of its dictionary by the
total original row count divided by that lower degree. Hence
iterative deletion of vertices whose current degrees have
dropped below that exponent by a slightly larger tolerance
removes only a power-smaller number of rows. Each deleted
vertex is charged once. This proves the regularization
claim, with mesh and deletion tolerances taken to zero in
the outer sense. For a common valuation profile one can
further select one anatomic box with vanishing exponent
loss. Its inherited dictionary and degree upper bounds
and the same degree pruning recover the original regular
exponents up to arbitrarily small allowances. Pass to
the weak degree measures on that box before choosing
the full divisors.
[Fixed-letter bound (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

We use the degree-one absolute residue estimate underlying
(AR1). If \(D\) has fixed positive numerical weight, polynomial
height, and \(\mathcal R\) is a list of \(R\) unit residues
modulo \(D\), then the number of positive complete integers
of polynomial height and radical at most \(n^{r+o(1)}\)
in those classes is at most
\[
                         R^{1/2}n^{r/2+o(1)}.                \tag{RB2}
\]
All uses below have \(r=1/3\) up to the balanced tolerances.
The estimate is uniform in the residue list and the modulus
with a fixed positive lower bound for its numerical weight.
The sign in the relation between two free letters permutes
their unit residues and does not change \(R\). The integer
being counted in (RB2) can be any original letter, because
the degree-one estimate only uses its own height cap and
radical cap. This also explains the free-target use of (AR1)
when the selected role is a summand.
[Absolute residue bound and its elementary degree-one proof (AR1)](ABC_ABSOLUTE_RESIDUE_TARGET_TAILS.md#one-full-divisor-and-the-multiplicity-threshold)

If a letter \(W\) has positive nonatomic mass, choose a
small compact positive-degree block and fix its full
product \(D\mid W\). For an arbitrarily small fixed \(f>0\)
the record dictionary can have exponent at most \(f\),
and the numerical weight stays bounded below by a positive
constant after that choice. Thus a single chart retains
at least \(n^{1/3-f-o(1)}\) rows. Its two free letters
\(U,V\) satisfy
\[
                   U\equiv \epsilon V\pmod D,\qquad
                   \epsilon\in\{1,-1\}.
\]
Both residues are units. The lower bound on the modulus
and the block are chosen before applying (RB2) at
sufficiently late budgets. In each use of this chart,
choose \(f\) smaller than the desired cell or degree
tolerance first. Multiple locks use their total record
cost in place of \(f\).
[Cheap full divisors from a positive nonatomic band (MD3)](ABC_MODULAR_DIGIT_SUPPORT_AND_ATOMIC_COMPLETION.md#completing-all-balanced-profiles-after-small-modulus-support-recovery)

## Equality of the three complete dictionary exponents

**RB3.** Every regular full-exponent balanced cohort satisfies
\[
 \alpha_A=\alpha_B=\alpha_C=\alpha,\qquad
 u_A=u_B=u_C=u=\frac13-\alpha,\qquad
               \frac16\le\alpha\le\frac13 .                  \tag{RB3}
\]

Fix a cheap \(W\)-lock with free roles \(U,V\) as above.
The complete dictionary caps on the chart are inherited from
the original regular cohort. Delete from this chart the rows
on \(U\)-values of chart degree less than
\(n^{u_U-f-\omega}\), with a fixed positive tolerance \(\omega\).
This deletes at most
\[
            n^{\alpha_U+u_U-f-\omega+o(1)}
                  =n^{1/3-f-\omega+o(1)}
\]
rows. If the threshold is below one, there are no such
occurring values. Every surviving \(U\) has that many
distinct \(V\) witnesses in its shared residue class.
Thus at most \(n^{\alpha_V-u_U+f+\omega+o(1)}\) classes
are needed. The degree upper bound \(n^{u_U+o(1)}\)
gives at least \(n^{\alpha_U-f-o(1)}\) distinct surviving
complete \(U\)-values. Applying (RB2) to those \(U\) gives
\[
 \alpha_U-f
       \le \frac{1/3+\alpha_V-u_U+f+\omega}{2}+o(1).
\]
Using (RB1) this implies
\(\alpha_U\le\alpha_V+3f+\omega+o(1)\).
Exchange \(U,V\) with the same congruence and let the fixed
allowances tend to zero to obtain \(\alpha_U=\alpha_V\).
At least two different roles have positive nonatomic mass
at a balanced endpoint, hence cheap locks on those two
roles give equality of all three dictionary exponents.
The same argument gives a quantitative inequality before
taking the balanced limit. Suppose the original row
exponent is \(d\), a regular complete degree exponent is
\(u_U\), and \(\alpha_U=d-u_U\). With the raw normalized
cap \(r_U\) on \(U\) and a record of cost at most \(f\)
on \(W\), it yields
\[
       d\le r_U-(\alpha_U-\alpha_V)+3f+O(\tau)+\eta ,          \tag{RB8}
\]
where \(\tau\) is the regularization mesh and \(\eta>0\)
is an arbitrary fixed counting allowance at sufficiently
small fixed budgets. Indeed the two sides of the count
above become \(\alpha_U-f-O(\tau)\) and
\((r_U+\alpha_V-u_U+f)/2+O(\tau)+\eta\), choosing the
degree-discard margin and the inner allowances inside
\(\eta\). If at least two marked roles have economical
locks, the resulting two opposite-role comparisons
connect the three dictionaries. At least one available
comparison then sees at least half the difference
between the largest and smallest dictionary exponents.
Use its larger dictionary as the counted letter \(U\).
A fixed spread gives a strict saving when its half
exceeds the raw-cap tolerance and the lock and counting
allowances in (RB8).
[Pairwise nonatomic lower bounds at the balanced endpoint](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)

## The actual residue-cell sizes

Set
\[
               s=2\alpha-\frac13=\alpha-u\ge0 .
\]
Fix one cheap \(W\)-lock and use the common residue
\(\Theta=U\bmod D=\epsilon V\bmod D\) on its rows.
The complete dictionaries in the chart have exponent
\(\alpha\) up to the arbitrarily small outer losses.
Write \(\mathcal U_a,\mathcal V_a\) for their complete
values belonging to a residue \(a\), with the sign for
\(\mathcal V_a\) absorbed in its label.

**RB4.** Given a fixed cell tolerance \(\delta>0\), choose
the chart price \(f\) small in terms of \(\delta\). Apart
from a power-smaller proportion of the chart rows, their
occupied cells have
\[
 |\mathcal U_a|,\ |\mathcal V_a|=n^{u+O(\delta)+o(1)},
 \qquad
 \#\{\text{original rows in cell }a\}
                           =n^{2u+O(\delta)+o(1)} .             \tag{RB4}
\]
There are \(n^{s+O(\delta)+o(1)}\) such cells. The formulas
describe exponent bounds, including lower exponent zero
when \(u=0\) or \(s=0\), with existing nonempty cells
giving the corresponding lower bound of one.

For the upper dictionary bounds, use the original chart
dictionaries before deleting any rows. There are at most
\(n^{\alpha-u-\delta+o(1)}=n^{s-\delta+o(1)}\) residues with
\(|\mathcal V_a|>n^{u+\delta}\). If the exponent is negative
there are none at late heights. Otherwise (RB2) bounds
all complete \(U\) over these residues by
\[
                n^{1/6+(s-\delta)/2+o(1)}
                     =n^{\alpha-\delta/2+o(1)} .
\]
Their total row count is at most \(n^{1/3-\delta/2+o(1)}\)
by the upper original degree of \(U\). Relative to the
lower bound \(n^{1/3-f-o(1)}\) for the chart, this is
power-smaller when \(f<\delta/4\). The same argument
handles the cells with too many \(U\)-values.

For a cell with fewer than \(n^{u-\delta}\) values on
the \(V\) side, every \(U\) in that cell accounts for at
most this many rows, since an exact pair determines a
row. Summing the \(U\) dictionary gives at most
\(n^{1/3-\delta+o(1)}\) rows on these cells. The symmetric
argument handles small \(U\) cells.

For an upper bound on the number of active cells,
first apply the degree deletion in the proof of (RB3)
with \(\omega<\delta/8\). Every remaining active residue
has at least \(n^{u-f-\omega}\) distinct \(V\) witnesses
in the original chart dictionary. So at most
\(n^{s+f+\omega+o(1)}\) cells remain active.
Cells having fewer than \(n^{2u-\delta}\) actual rows
then contain in total at most
\(n^{1/3+f+\omega-\delta+o(1)}\) rows, again power-smaller
relative to the chart if \(2f+\omega<\delta\). Each of the
other cells has at most \(n^{2u+2\delta}\) rows by its
two complete dictionary caps. These bounds prove (RB4)
and give the lower bound on the count of cells by
division of the retained row count.

In particular a typical cell is dense at the exponent
scale on the product of its two complete dictionaries.
The third letter remains the actual signed linear
combination in each row, with the original radical cap.
The argument gives no arithmetic bound on the absolute
density of an individual cell.

The complete residue lists also have a multiplicative
restriction. At \(s>0\), every retained list of size
\(n^{s-o(1)}\) grows by \(n^{s-o(1)}\) against every
auxiliary list whose product has polynomial size in \(n\).
The same is true after any fixed nonzero pointwise
power. This follows from charging the original complete
integers by a relative expansion factor in (RB2).
[Relative expansion and bounded products (RX1–RX2)](ABC_RESIDUE_EXPANSION_AND_JOINT_RAY_COUNTS.md#every-large-balanced-residue-list-has-maximal-finite-order-growth)

## Conditional entropy identities and three registers

Use the uniform law on the cheap chart. All entropy
equalities below are in units of \(\log n\), with the
fixed outer tolerances and total chart costs made
arbitrarily small. For fixed cell tolerance \(\delta\),
the exceptional mass in (RB4) is power-small if the
chart cost is selected as above. Entropies of deleting
such rows change by \(o(\log n)\), since the original
dictionary and row counts have bounded normalized
logarithms. Use (RB4), then decrease \(\delta\), to obtain
\[
 \begin{split}
  H(\Theta)&=(s+o(1))\log n,\
  H(U\mid\Theta)&=(u+o(1))\log n,\qquad
  H(V\mid\Theta)=(u+o(1))\log n,\
  H(U,V\mid\Theta)&=(2u+o(1))\log n,\qquad
  I(U;V\mid\Theta)=o(\log n).                                \tag{RB5}
 \end{split}
\]
For the first equality, the upper count of active cells
bounds its entropy from above. The upper row count per
cell bounds \(H(U,V\mid\Theta)\) from above. Since the
pair determines the row and \(H(\text{row})=(1/3+o(1))\log n\)
with the arbitrarily small chart price absorbed, this
also gives the matching lower bound on \(H(\Theta)\).
The conditional complete-letter entropies then follow
from \(H(U)=H(V)=(\alpha+o(1))\log n\), as follows
from the regular degree and dictionary caps, and from
\(\Theta\) being a function of either letter.

A bounded list of cheap locks on the same marked role
has one joint residue register modulo their least
common multiple, again with fixed positive numerical
weight and polynomial height. Applying (RB5) to the
list and to its nonempty sublists shows that all of
these joint registers have entropy \((s+o(1))\log n\).
The residue at any one of these locks therefore gives
the other residues with conditional entropy \(o(\log n)\).

Suppose that all three roles supply cheap full locks,
which can be fixed on one chart at arbitrarily small
combined price. This holds, for example, when \(u>0\):
pigeonhole one complete letter in the regular cohort,
obtaining a fibre of exponent \(u-o(1)\); (FL2) forces
positive nonatomic mass on both other letters, and
the orientations can be exchanged. Let \(J_A,J_B,J_C\)
be the shared residue register of the two letters
opposite each locked role, using a consistent sign
for each register.
[Positive nonatomic mass in positive complete fibres (FL2)](ABC_FIXED_LETTER_RADICAL_CAPS.md#strict-positive-scale-inequalities-on-the-two-free-letters)

Then \(H(J_T)=(s+o(1))\log n\) on each role. Also
\[
 I(J_B;J_C\mid J_A),\quad
 I(J_A;J_C\mid J_B),\quad
 I(J_A;J_B\mid J_C)=o(\log n).                               \tag{RB6}
\]
For example, given \(J_A\) the two complete letters
\(B,C\) have conditional mutual information \(o(\log n)\)
by (RB5); \(J_C\) is a function of \(B\), and \(J_B\)
a function of \(C\). Data processing gives the first
bound. The other two are obtained by exchanging roles.

The normalized pairwise mutual informations of
\(J_A,J_B,J_C\) therefore agree up to \(o(1)\). Write
\(z\) for a common limiting value on a subsequence.
To see the equality, subtract the conditional mutual
informations in the standard chain rule
\(I(J_A;J_B)-I(J_A;J_B\mid J_C)
 =I(J_A;J_C)-I(J_A;J_C\mid J_B)\).
The pair entropies are \(2s-z+o(1)\), and the triple
entropy is \(3s-2z+o(1)\), in units of \(\log n\).
Each register pair is a function of the remaining
complete letter, which gives
\[
              \max(0,3\alpha-2/3)\le z\le s .                 \tag{RB7}
\]
In particular the three residue registers have positive
pairwise mutual information on the exponent scale
when \(\alpha>2/9\). Their pairwise conditional
independence in (RB6) is an entropy assertion at that
scale.

## Consequence for shared CRT charges

If one works in the setting of (SW2) with the common
cheap \(W\)-lock and an original near-saturating regular
cohort, the typical-cell pruning above allows the
opposite local dictionary exponent \(u+\delta\),
in place of the general one-letter bound \(1/6\).
Adding a shared numerical interval only decreases
the cell sizes. Thus the repeated opposite-dictionary
charge in (SW3) is at most
\[
                         \alpha+(\ell-1)(u+\delta)
\]
on the retained family, with the chart price and all
remaining allowances paid separately. The multiplicity
and packet hypotheses of (SW2) are imposed on the
retained witnesses after trimming. For any fixed
\(\ell\) one can take \(\delta\) arbitrarily small
and then take the chart price smaller. This is an
improvement over summing \(\ell\alpha\) when
\(\ell>1\) and \(\alpha>1/6\). The conditional packet extraction,
the overlap bounds and the numerical least-common-multiple
threshold are still required by (SW2). In the limiting
formula a strict target saving requires
\[
              \sum_i\bigl(\min(u,t_i)+\gamma_i\bigr)<u ,
\]
since the target exponent to beat is \(1/3-u=\alpha\).
Consequently every \(t_i<u\) and \(\sum_i t_i<u\) if
\(u>0\); at \(u=0\) this formula gives no strict saving.
In a tall chart, uniformly bounded scaled orders
\(S_{\max}<6\) still leave a fixed numerical deficit
for the strict least-common-multiple threshold when the
common lock is chosen at arbitrarily small radical cost
from the same bounded spectrum. The extracted numerical
weight is at most \(S_{\max}\sum_i t_i\) before adding
that common lock's vanishing weight.
[Shared opposite dictionary charges (SW2–SW3)](ABC_SHARED_CRT_WITNESS_DICTIONARIES.md#data-and-shared-dictionary-bound)

The equality and residue-cell identities do not force
\(\alpha=1/6\), a low target radical, or any
individual numerical packet inequality. The inequalities leave a range of possible common
dictionary exponents. Closing the unrestricted endpoint
requires an arithmetic exclusion or a weighted target-image
bound on these regular residue components.
