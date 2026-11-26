# Mixed compact prime packing across the letters

At the minimum complete-dictionary scale one can sample a single
common column using fixed anchors in both other roles. Full prime
powers on the two resulting lists of affine forms share one
Chinese-remainder reconstruction bound. In particular the highest
diffuse scaled orders on any letter constrain the compact high-order
prime powers on every other letter.

## A mixed packing inequality in actual additive graphs

Put \(c_0=1/6\), \(n=X^L\), and use primitive positive triples
\(a+b=c\) in \(X/2<c\le X\). Assume that each individual radical
is at most \(n^{2c_0+\kappa}\) and each complete dictionary has
at most \(n^{c_0+\kappa}\) values. Write the additive equation
as a signed sum of the complete values in their respective roles.

Fix two roles \(A,B\), which can be the same. On \(A\), suppose
every retained row has a marked prime power with normalized
prime cost in \([t_-,t_+]\subset(0,\infty)\) and scaled
order in \([s,S]\subset(0,\infty)\). On \(B\), suppose there
are prime powers in a compact scaled band \([u,U]\) whose
individual normalized costs are at most \(b_0\) and whose
combined normalized cost is at least \(\rho>0\) on every
row. The normalized prime cost is \(\log_n p\) for a prime
\(p\), and its scaled order is \(Lv_p(W)\) on its complete
letter \(W\). The marked and filler role requirements
can be used on different copies of the same role.
Assume all values of the free query role are in one
residue class modulo \(G\ge X^\gamma\), with \(1\le G\le X\).
For distinct \(A,B\) use the third role as the free query
role. When \(A=B\), use either different role. Take \(G=1\)
and \(\gamma=0\) without a congruence hypothesis.

Fix an integer \(k\ge0\) and \(q\ge0\), with at least one
positive. The corresponding marked or filler conditions
can be omitted when its choice is zero. If
\[
 \begin{split}
       k t_+ + q + b_0{\bf1}_{q>0}&<c_0,\\
                         k s t_- +u q&>1-\gamma ,
 \end{split}
 \tag{MPK1}
\]
there are positive \(\kappa_0,L_0,\delta\), depending
only on the fixed parameters and the strict margins,
such that the retained graph has
\(O_L(n^{2c_0-\delta})\) rows for \(0<L<L_0\)
and \(0\le\kappa\le\kappa_0\). All constants, including
\(\gamma\), are held fixed before the height limit.
In particular the hypothesis cannot hold along a full
minimum-dictionary endpoint sequence with fixed margins.

Choose a fixed number \(\ell\) of filler copies with
\(\ell\rho/2>q\), or \(\ell=0\) when \(q=0\). Set
\(m=k+\ell\). For a free query value with degree \(d_v\),
each role different from its own has \(d_v\) distinct
complete values in its neighborhood. This follows because
fixing two roles determines the third value on each
original row. Choose \(k\) anchors in the role opposite
\(A\) and the query, then \(\ell\) anchors in the role
opposite \(B\) and the query. Each anchor gives the
desired marked or filler value as an affine form in
the complete query. Before changing the query parameter
its slope is \(+1\) or \(-1\) according to the signs.
At each choice, an earlier affine zero forbids at
most one of the available anchors. This is true also
across the two kinds of anchors; the signed zero
uniquely determines a proposed anchor of the given role.
Thus the forms can be chosen with distinct zeros and
there are at least \((d_v-m)^m\) valid choices if \(d_v>m\).
No independence of the neighborhoods in the two roles
is required.

If \(E\ge n^{2c_0-\tau}\), the low-degree query vertices
have total degree at most \(2m n^{c_0+\kappa}\), which
is negligible at small fixed \(\tau,\kappa\). By the
power mean, the total number of ordered stars is
\(\gg_m E^m/n^{(m-1)(c_0+\kappa)}\). There are at most
\(n^{m(c_0+\kappa)}\) anchor lists. A fixed valid list
therefore has at least a constant times
\[
                   n^{c_0-m\tau-(2m-1)\kappa}
\]
common query values. The same proof covers anchors in
a single role. The mixed choice matters when the two
types of packets come from different original letters.

Write the query as \(x_0+G z\). The integer parameter
has interval length \(O(X^{1-\gamma})\). The target
forms are now nonproportional affine forms
\(F_i(z)=\sigma_i G z+a_i\), where \(\sigma_i\in\{-1,1\}\)
and all coefficients have polynomial height.
Put into one prescribed pool all primes of \(G\) and of
the pairwise determinants of these forms. The latter
integers are nonzero because the zeros were distinct,
and the pool has polynomial height for the fixed list
size. Off this pool the prime supports of the forms
on a query are disjoint, and all slopes are units.
At every fixed anchor the original fixed-letter count
(FL1) discards the queries carrying pool radical
of normalized mass at least \(\epsilon_p>0\) on its
target \(F_i(z)\) with normalized exponent at most
\(c_0+\kappa/2-\epsilon_p/2\), up to the arbitrarily
small coefficient allowance. Apply this separately
on all original fixed-anchor fibres. Choose
\(\epsilon_p<\min(t_-,\rho/2)\) when both are present,
and choose \(\tau,\kappa\) and the allowances small
after fixing the list. Thus on the surviving query
class no marked prime was lost, and each filler form
has at least \(\rho/2\) mass off the pool.
[Fixed-letter cap with prescribed prime pools (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

On each query select the marked full prime power
on each marked form and fill greedily from the filler
forms to cost at least \(q\), with overshoot at most \(b_0\).
The prime powers are pairwise coprime by the pool
construction. Their product has numerical exponent
at least \(k s t_-+u q\); each selected congruence
for \(z\) has one residue at its prime power.
By (MPK1) their Chinese remainder solution determines
at most one integer \(z\) in the parameter interval.
The squarefree support of the product has at most
\(n^{k t_++q+b_0{\bf1}_{q>0}}\) possibilities.
Its height is at most \(X^m\), so after choosing this
support there are subpower many assignments of the
prime powers and their roles by fixed-support enumeration
and the divisor bound. The resulting total exponent is
smaller than \(c_0\) by a fixed amount. This contradicts
the retained star exponent for sufficiently small fixed
\(\tau,\kappa\). Taking \(\delta=\tau\) proves (MPK1).
[Single-role compact packing](ABC_MINIMUM_DICTIONARY_COMPACT_PRIME_PACKING.md)

## Coupling compact prime masses to the diffuse support on any role

Consider a limiting full minimum-dictionary profile.
Take \(\gamma\) as a lower limit for the height exponent
of a common original congruence modulus, using all
three original residue classes as in the exact triangle
transfer. Put \(H=1-\gamma\). Let \(b\) be the largest
upper support point of the three nonatomic scaled-degree
measures. By the one-role progression-width bound,
\(0\le b\le6H\).

Suppose a compact prime pattern survives on a
full-exponent cohort on some original role: on each
of its edges there is a prime power on that role
whose scaled order approaches \(s>0\) and whose
normalized prime cost approaches \(t>0\).
Keep the limiting degree profiles on this cohort.
Then for every positive integer \(k\) with \(kt<c_0\),
\[
                    kst+b(c_0-kt)\le H .
 \tag{MPK2}
\]
For \(b=0\), use only the repeated marked powers and
the Chinese-remainder count above. Otherwise use as
filler the role whose nonatomic measure has upper
support \(b\). For any fixed \(u<b\) sufficiently
close to \(b\), choose a compact positive band
above \(u\) carrying positive nonatomic mass.
Discard small neighborhoods of its finitely many
atoms of mass above a small cutoff. The remaining
band retains positive measure, and the individual
prime costs can be made at most any prescribed
positive \(b_0\) in late boxes. This uses exactly
the compactness argument for the affordable full
powers in (BC0). It is valid on the same cohort
as the marked prime. The filler may come from
an original role different from the marked role.
Take \(q\) just below \(c_0-kt\).
If (MPK2) fails there is a strict margin in both
conditions of (MPK1) after making the parameter
boxes and the slack sufficiently small.
Use a smaller fixed lower height exponent for \(G\)
when taking its limiting value. The count contradicts
the full cohort exponent. This proves (MPK2).
[Affordable compact full divisors (BC0)](ABC_BOUNDED_DEGREE_CURVES_AND_SECANT_GCD.md#uniformly-affordable-full-divisors)
· [Exact tagged triangle normalization](ABC_HEREDITARY_MINIMUM_TRIANGLE_TRANSFER.md#fixing-the-exact-scaling-and-the-three-residue-classes)

In particular a finite scaled degree \(s>6\)
of a persistent prime pattern has \(t\le1/s<c_0\)
from the complete-letter height. Taking \(k=1\)
first gives \(st\le H\), also when the filler is
omitted. Take \(k=\lceil c_0/t\rceil-1\) so
\(c_0-t\le kt<c_0\). As \(s>b\), (MPK2) now gives
\[
 H\ \ge\ bc_0+(s-b)kt
     \ \ge\ sc_0-(s-b)H/s .
\]
Consequently
\[
                 s\le 6H+\sqrt{36H^2-6Hb}.
\]
It follows that all three finite scaled-degree measures
have support in
\[
 \left[0,\ \max\{6,\ 6H+\sqrt{36H^2-6Hb}\}\right].
 \tag{MPK3}
\]
To check the passage from marked primes to the full
measures, take a compact band \([s_*,S_*]\) strictly
beyond the displayed cutoff with positive mass in
the limiting profile. Work in refined boxes so this
mass and the filler profile hold uniformly on the
rows of any retained cohort. Choose a small fixed
cost cutoff \(b_0<c_0-H/s_*\). For costs at least
\(b_0\), the possible pattern parameters lie in the
compact set of \((s,t)\) with \(s\in[s_*,S_*]\),
\(t\in[b_0,1/s_*]\) and \(st\le1\).
At each of these parameters the preceding calculation
forces a violation of a specific repeated-packet test:
either \(st>H\), or one can find \(k\) with \(kt<c_0\)
and \(kst+b(c_0-kt)>H\). Each strict violation gives
a neighborhood with a fixed choice of the parameters
in (MPK1). Finitely many such neighborhoods cover
the compact set, so the union of their excluded
cohorts has a fixed power saving. In the complement
the high-band primes are all small-cost filler.
Since \(b_0+H/s_*<c_0\), the filler-only (MPK1) excludes
positive high-band mass after taking enough copies.
Continuity windows on the original profiles justify
retaining a positive uniform mass in this band.

At \(\gamma=0\) the new upper cutoff is
\(6+\sqrt{36-6b}\), improving \(12\) whenever
the maximal diffuse cutoff \(b\) is positive.
If \(b=6\), even compact atomic prime weights
above scaled order six are excluded on every role.
The statement concerns finite scaled-degree bands.
The common residue width can itself come from
the high-order tail alignment in the compact
packing theorem, as well as from exact triangle
scaling. In the tail alignment, write \(e_{\rm tot}\)
for the sum of the three escaped moments. One may use
\(H=1-e_{\rm tot}\) in (MPK2) and (MPK3), by first
using any fixed numerical exponent below \(e_{\rm tot}\)
and then taking the limit. If the aligned modulus
initially has an upper height bound larger than \(X\),
fullness of the residue-aligned graph forces \(G\le X\)
at large heights, since \(G>X\) leaves at most one
complete value in each prescribed residue. No abundance assertion for bounded-radical
triangles follows from this packing constraint.

### The indivisible atom above the mixed cutoff

Let \(R_*=6H+\sqrt{36H^2-6Hb}\). On each role any
finite scaled-degree mass above \(R_*\) is confined
to a single possible atom of mass \(c_0\). When
this atom occurs the remaining mass of that role
is entirely nonatomic and has mass \(c_0\).
In particular, in the absence of an atom of mass
\(c_0\), the full finite support of the role is contained
in \([0,R_*]\).
\[
  \nu_W((R_*,\infty))=0
  \quad\text{or}\quad
  \nu_W=\nu_{W,\mathrm{na}}+c_0\delta_{s_W},
  \quad c_W=c_0,\ s_W>R_* .
 \tag{MPK4}
\]

To see this, the argument from (MPK2) giving the
quadratic cutoff works whenever \(s>b\) and \(t<c_0\),
whether or not \(s>6\). A compact band above \(R_*\)
with positive mass and no atom of mass at least \(c_0\)
has every individual prime cost uniformly below \(c_0\)
on sufficiently late refined boxes. Otherwise an
exact valuation position of cost tending to \(c_0\)
would give an atom of at least that mass in the
compact band. Repeat the finite-cover proof above,
with the compact cost interval bounded away from
\(c_0\), and the filler-only lower-cost branch.
It excludes the band. Hence any remaining mass is
on an atom of size at least \(c_0\). Pigeonholing a
fixed complete opposite letter gives a fibre of exponent
\(c_0-o(1)\) in the same refined profiles. By (FL2)
its free-letter nonatomic mass is at least \(c_0\).
The total mass of a letter is \(2c_0\), so the atom
must exhaust its atomic mass. Its numerical moment
satisfies \(c_0s_W+\int r\,d\nu_{W,\mathrm{na}}(r)\le1\).
[Nonatomic lower bound in a fixed fibre (FL2)](ABC_FIXED_LETTER_RADICAL_CAPS.md#strict-positive-scale-inequalities-on-the-two-free-letters)

The packing inequalities still allow balanced numerical
profiles. For illustration, for \(H=1\) take on each role
nonatomic density \(1/12\) on \([2,4]\) and an atom
of mass \(1/6\) at \(3\). Each measure has total mass
\(1/3\) and numerical moment one. Every full prime
packet in these bands weighs at most four times its
normalized radical cost, so a packet of cost less
than \(1/6\) lies below the reconstruction threshold.
This is a compatibility check for the new packing
inequalities on the proposed profiles, with no
arithmetic realization claimed.
