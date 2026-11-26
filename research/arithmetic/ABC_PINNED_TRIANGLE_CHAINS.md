# Affine lifts and inherited moduli in pinned triangle chains

A sequence of critical pinned triangle transfers retains more than
the primitive output counts. A terminal complete value gives a
family of actual shifted values at every earlier level. Fixed
high-order parts in two different original legs pass to residue
classes on two output roles. The next fixed-gap transfer makes
these two classes constrain one common scaling factor.

## Chain hypotheses and the height convention

Put \(\alpha=1/6\), let \(X_0\) be the first shell height, and
write \(n=X_0^L\). Take the height limit at fixed \(L>0\)
before the small-budget limit. A level \(j\) has a family
\(\mathcal E_j\) of exponent \(2\alpha-o(1)\) in base \(n\),
consisting of positive primitive triples \(a_j+b_j=c_j\)
in a shell \(X_j/2<c_j\le X_j\), chosen with \(X_j\le X_0\).
Its three complete
dictionaries have sizes at most \(n^{\alpha+o(1)}\), and its
individual radical caps are \(n^{2\alpha+o(1)}\). The chain
length \(k\) is fixed before both limits.

A transfer at level \(j\) uses a selected tagged relation of
\(n^{4\alpha-o(1)}\) triangles on fixed complete \(b_j\)-fibres.
The three increasing values are in the \(a_j\) role, with
normalized primitive difference radical at most
\(n^{6\alpha+o(1)}\). Use the exact scaling-and-pin selection
to obtain one gcd \(G_j\), one first \(a_j\)-value \(p_j\),
and one common fixed \(b_j\)-value \(T_j\). Write the
primitive outputs in the standard order as
\[
   a_{j+1}=(x_2-p_j)/G_j,\qquad
   b_{j+1}=(x_3-x_2)/G_j,\qquad
   c_{j+1}=(x_3-p_j)/G_j.
\]
The next input family is the selected primitive output in a
dyadic height shell. The hypothesis at every stage concerns
a full relation on this actual retained family. If further
refinements are made, availability of the next full relation
is assumed on their retained output. The existing hereditary
transfer gives the cap and dictionary conclusions within each
such stage. Every fixed tolerance below is chosen before
taking the corresponding inner height limit.
[Exact scaling, pinning and hereditary dictionary bounds](ABC_HEREDITARY_MINIMUM_TRIANGLE_TRANSFER.md#pinned-transfer-at-the-critical-radical-cap)

## Affine propagation and the cumulative numerical floor

For either \(x=a_{j+1}\) or \(x=c_{j+1}\), the parent family
contains the original row
\[
             (p_j+G_jx,\ T_j,\ p_j+T_j+G_jx).                 \tag{PT1}
\]
Thus one value in either of these two output roles lifts to
both moving parent roles. Begin with a list containing the
identity function at level zero. At step \(j\) compose each
function in the list with both maps
\(x\mapsto p_j+G_jx\) and \(x\mapsto p_j+T_j+G_jx\), then
append the identity. On any eventual terminal value in the
\(a_k\) or \(c_k\) role, the nonidentity functions evaluate to
actual earlier \(a\)- or \(c\)-values. Every evaluation comes
from the corresponding original fixed-\(T_j\) fibre and
retains its individual cap. In particular from each level
\(j<k\) there are \(2^{k-j-1}\) forms in each of its two
moving roles.

The pins can be chosen so that the forms in the list always
have distinct rational zeros. Initially the only zero is
zero, and all zeros remain nonpositive. If the old zeros
are \(z_i\), the two copies have zeros
\((z_i-p_j)/G_j\) and \((z_i-p_j-T_j)/G_j\).
A collision between these copies requires \(T_j\) to equal
a difference between two old zeros. The finitely many
forbidden \(T_j\)'s account for \(O_k(n^{3\alpha+o(1)})\)
tagged triangles after the exact-cell selection. Discard
them before choosing the pins by averaging on the full
\(n^{4\alpha-o(1)}\) retained scale. Each copied zero is strictly
negative, so the appended zero cannot coincide with one of
them. All form coefficients have polynomial height in \(X_0\).
For example their slopes and intercepts are positive
products and sums of the actual scaling factors and pins.
Evaluating a nonidentity form on one positive terminal
value bounds both its slope and its intercept by \(X_0\).

Write
\[
                      P_k=\prod_{j=0}^{k-1}G_j .
\]
One of the level-zero affine lifts has slope \(P_k\).
The complete \(a_k\)-dictionary has size \(n^{\alpha-o(1)}\)
because its row count is \(n^{2\alpha-o(1)}\) and one other
dictionary has size at most \(n^{\alpha+o(1)}\). The integer
parameter \(x=a_k\) lies in an interval of length at most
\(\min(X_k,X_0/P_k)\), using the bound on its original lift
for the second term. Apply the rank-one fixed-fibre bound
to that lift in the original \(T_0\)-fibre. The exponent
\(r=\alpha/2\) suffices in that bound. There is an absolute
\(\theta>0\), independent of the fixed chain length, such
that on a full-exponent chain
\[
               \log_{X_0}X_k\ge\theta,\qquad
               \log_{X_0}P_k\le 1-\theta .                     \tag{PT2}
\]
Use any smaller fixed \(\theta\) than the width threshold
in the rank bound to absorb constant factors. The same
bounds apply to every prefix. In particular the shell
height change in applying (FL1) at an earlier level has
a bounded factor \(\log X_0/\log X_j\).
[Small affine features in a fixed fibre](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#a-small-feature-bound-on-a-single-fibre)

## High-order parts pass to the next exact scale

At any stage \(j\), fix a finite cutoff \(R>0\) in the
common scaled degree \(Lv_p\), and on each of the three
tagged original legs fix the whole parts of both \(a_j\)
and \(c_j\) at orders \(Lv_p\ge R\). Denote the two
fixed parts on leg \(i\) by \(F_{a,i},F_{c,i}\), and put
\(M_i=F_{a,i}F_{c,i}\). The normalized radical cost is
at most \(1/R\) for each part because its numerical height
in base \(X_0\) is at most one. Uniform fixed-support
enumeration bounds the loss to fix all six parts by
\(n^{6/R+o_X(1)}\). This fixes records separately on the
three legs, which keeps a large subclass even for an
arbitrary selected triangle relation. In applying the
exact-cell and pin selections, take the tolerances larger
than this loss, or make power-saving exclusions before
fixing the records. Then let \(R\) increase as the losses
tend to zero.
[Uniform fixed-support enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

At the fixed \(T_j\), the \(a_j\)-value \(x_i\) on leg
\(i\) lies in one residue class modulo \(M_i\).
Indeed \(F_{a,i}\mid x_i\), \(F_{c,i}\mid x_i+T_j\),
and the two factors are coprime by primitivity. From
\(x_2=p_j+G_j a_{j+1}\) and \(x_3=p_j+G_j c_{j+1}\)
we obtain prescribed classes for \(a_{j+1},c_{j+1}\)
modulo
\[
        Q_2=M_2/\gcd(M_2,G_j),\qquad
        Q_3=M_3/\gcd(M_3,G_j),
\]
respectively. Division of the congruence by its gcd is
valid because the original rows give a solution. These
classes persist on every later restricted subfamily of
the output.

On the next fixed-\(b_{j+1}\) fibre, each difference
between \(a_{j+1}\)-values equals the corresponding
difference between \(c_{j+1}\)-values. It is therefore a
multiple of both \(Q_2\) and \(Q_3\). The next exact gcd
is divisible by their lcm. In particular
\[
                    \operatorname{lcm}(M_2,M_3)
                              \mid G_jG_{j+1}.                \tag{PT3}
\]
At a prime with valuations \(e_2,e_3\) in the two
products and valuation \(g\) in \(G_j\), the statement
is \(g+\max\{(e_2-g)_+,(e_3-g)_+\}\ge\max\{e_2,e_3\}\).
This accounts for all shared powers and for overlap
between the incoming residues and the next records.

One can also use this inherited modulus without a next
triangle transfer. Fix \(b_{j+1}\) in the full output so
its fibre has exponent \(\alpha-o(1)\). The complete
\(a_{j+1}\)-values are in a class modulo
\(Q_*=\operatorname{lcm}(Q_2,Q_3)\). Apply the rank-one
fixed-fibre bound to their level-zero lift, parameterized
by this progression. This gives
\(\log_{X_0}(X_{j+1}/Q_*)\ge\theta\), reducing \(\theta\)
if necessary. Thus the inherited classes shorten the
parameter interval even before a second gcd is fixed.

For the moment interpretation, let \(e_{W,j}\) denote
escaped numerical weight measured in base \(X_0\).
It is the limit of the complete numerical height
\(\log_{X_0} W\) minus the first moment of the weak degree
measure in the common units \(Lv_p\) and \(\log_n p\).
At a fixed continuity cutoff, the fixed tail weight
tends to the complete numerical height minus its
truncated moment. With a slow-cutoff choice that
captures the escaping moment along the selected chain,
(PT3) gives
\[
          \log_{X_0}(G_jG_{j+1})
                      \ge e_{a,j}+e_{c,j}-o(1).             \tag{PT4}
\]
The finite-cutoff inequality (PT3) holds for the actual
tail records. Formula (PT4) uses a diagonal choice
tending to infinity slowly enough for the truncated
moments of the retained profiles, with the next full
relation available on the selected outputs. Summing
over disjoint consecutive step pairs and using (PT2)
bounds the total escaped moving weight charged across
those pairs by \(1-\theta+o(1)\). Thus a uniformly
positive base-\(X_0\) escaped weight in the two moving
roles at every stage bounds the possible chain length.

## CRT packing on the ancestral forms

Use the distinct-zero selection in (PT1), and take the
terminal \(a_k\)-dictionary as the query set. Its
exponent is \(\alpha-o(1)\). Use any fixed collection
of nonidentity ancestor forms having prescribed compact
prime bands from their actual complete roles. Put in a
pool all primes dividing a slope or a nonzero pair
determinant of those forms. This has polynomial height.
Outside it the prime supports of different forms at
the same query are disjoint, and every slope is a
unit at those primes.

Each form is injective on the query set and its value
participates in a single fixed \(T_j\)-fibre at its
ancestral level. If it shares normalized radical mass
at least \(\epsilon>0\) with the pool, (FL1) on that
fibre gives an exponent at most \(\alpha-\epsilon/2+o(1)\)
in the original base \(n\), allowing an arbitrarily
small fixed exponent loss. The caps in that application
are the actual level caps converted to height \(X_j\).
By (PT2) the conversion and the pool height exponent
are bounded uniformly at the fixed chain length.
Remove this power-smaller class for each form.
[Fixed-letter count with prescribed prime pools (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

On the surviving queries, take full marked powers of
individual costs in fixed \([t_-,t_+]\) and scaled
orders at least \(s>0\) on their chosen forms. A mark
survives if \(t_->\epsilon\). One can use several
marked types on distinct forms. For a filler type with
at least \(\rho>2\epsilon\) mass per form in a fixed
compact band of scaled orders at least \(u>0\), and
maximum individual cost \(b\), use enough distinct
copies that their aggregate available mass exceeds a
prescribed \(q>0\). Greedy selection costs at most
\(q+b\) and supplies numerical weight at least \(uq\).
The same argument works with separate filler types
on distinct lists of forms.

The product of all selected prime powers determines
one CRT residue for the terminal query. A lower
numerical weight exceeding \(\log_{X_0}X_k\) by a
fixed margin leaves at most one query in the interval.
If the combined radical cost, including all filler
overshoots, is below \(\alpha\) by a fixed margin, enumerate
the supported powers and their form assignments at that
cost plus a height-subpower loss. This contradicts the
terminal dictionary exponent. Thus such simultaneous
weight and cost bounds are impossible along a full
chain. The packing mixes actual prime powers across
levels and roles without a common-column extraction.
[Mixed compact packing in one graph](ABC_CROSS_ROLE_COMPACT_PACKING.md)

For example, a compact band of positive nonatomic
mass with lower scaled order \(u\) on one ancestor
role can supply arbitrarily small individual filler
costs after removing neighborhoods of the finitely
many larger atoms. If there are sufficiently many
copies of that role, the filler-only test gives
\(\alpha u\le\liminf\log_{X_0}X_k\). A uniform choice
of the compact band and its positive mass is made
before taking the chain length and the small-budget
limit. Orders measured in the relative budget of
level \(j\) contribute their weight multiplied by
\(\log_{X_0} X_j\) when converted to these common
units. Thus the numerical constraint concerns the
whole chain width. It allows chains with small
exact scalings and no forced collapse of the
terminal height. The input still requires critical
triangle relations on the retained families.
