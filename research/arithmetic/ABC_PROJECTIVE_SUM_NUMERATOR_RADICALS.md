# Fresh numerator radicals in sums of projective abc coordinates

In a tuple of primitive rows, a prime used on one original
row alone survives in the reduced denominator of the sum of
their projective coordinates. The full surviving prime powers
also determine residues of the product of the other rows'
denominators once the output numerator is known. Exposing a
numerically larger block recovers the whole remaining block.
A second estimate applies the two-letter abc count to two
shifted numerators of the same sum.

## Statements and the height-first setting

Put \(n=X^L\), and take the height limit at each fixed \(L>0\)
before sending \(L\) to zero. In every slot use positive
primitive rows \(a+b=c\) with \(X/2<c\le X\) and
\(\operatorname{rad}(abc)\le n\), allowing arbitrary subfamilies.
Fix \(k\ge2\). Denote the row in slot \(i\) by \(x_i\).
The slots may use different rational functions or assignments
\(f_i(x_i)=A_i/D_i\) in reduced terms. Assume that their
reduced heights are bounded by \(X^K\) for a fixed \(K\).
The required denominator hypothesis is as follows.
There is a prescribed prime pool \(P_*\) of polynomial
height such that for every \(p\notin P_*\),
\[
 v_p(D_i)=d_{i,W}v_p(W)
      \quad\text{if \(p\) divides the original role \(W\)
                 of \(x_i\)},                               \tag{NR1}
\]
and \(v_p(D_i)=0\) if \(p\nmid a_i b_i c_i\).
All \(d_{i,W}\) are fixed positive integers. Here and below
a prime pool of polynomial height is supported on an
integer of size at most \(X^{K_1}\) for a fixed \(K_1\).
The notation \(\operatorname{rad}_{P^c}(v)\) denotes the
product of the prime divisors of \(v\) outside \(P\).

For \(\mathbf x=(x_1,\ldots,x_k)\), write
\[
         t=\sum_{i=1}^k w_i f_i(x_i),\qquad
         t+s=\frac{N_s}{D_s}
\]
in reduced terms with \(D_s>0\), for nonzero rational weights
\(w_i\) and a rational shift \(s\). We allow the weights and
shift in the following first estimate to vary with \(X\) and
\(L\) subject to a fixed polynomial bound on their reduced
heights. For any prescribed polynomial-height pool \(P\)
and any \(\eta>0\), if \(k\ge3\) then there are
\(\delta,L_0>0\) such that
\[
 \#\left\{\mathbf x:
    N_s=0\ \text{or}\
    \operatorname{rad}_{P^c}(|N_s|)
                \le n^{\lfloor(k-1)/2\rfloor/3-\eta}\right\}
       \ll X^{(k/3-\delta)L+\varepsilon}\qquad(0<L<L_0).
                                                               \tag{NR2}
\]
Constants are uniform in all the assignments, weights,
shift and prescribed pool under the displayed bounds.
In particular this estimate holds for the radical of the
numerator on primes outside all input letters, with a
possibly smaller fixed \(\delta\).

For the second estimate, assume \(P_*\) is supported on
a fixed finite set of rational primes of bounded size.
Take the weights and shifts from fixed finite sets of
rational numbers. Let \(s_1,\ldots,s_b\) be distinct such
shifts, where \(b\ge2\) is fixed, and let
\(\mathcal I(\mathbf x)\) be the support of all input
letters. For each \(\eta>0\), outside
\(O(X^{(k/3-\delta)L+\varepsilon})\) tuples, the following
hold simultaneously for all \(j\ne l\):
\[
 \begin{split}
     \log_n\operatorname{rad}(D_{s_j})&\ge k-\eta,\\
     \log_n\operatorname{rad}_{\mathcal I(\mathbf x)^c}(|N_{s_j}|)
      +\log_n\operatorname{rad}_{\mathcal I(\mathbf x)^c}(|N_{s_l}|)
                   &\ge 2k/3-\eta .
 \end{split}                                                   \tag{NR3}
\]
All these numerators are nonzero on the retained tuples.
The finite-shift restriction controls the raw radical
introduced in a second additive triple below.

In particular, for \(s_j=0,-1\), the reduced output
\(t=N/D\) defines the positive primitive triple obtained
by reordering the absolute values of \(N,N-D,D\).
On all but a power-smaller set of input tuples it satisfies
\[
                      \operatorname{rad}\bigl(|ND(N-D)|\bigr)
                                \ge n^{5k/3-2\eta}.            \tag{NR4}
\]
The assignments include the coordinate \(F(a/c)\),
its fixed positive powers, and the Frey \(j\)-invariant,
whose reduced denominator away from two is \((abc)^2\)
at the primes dividing \(abc\).
[Projective reciprocal numerator and denominator](ABC_THREE_CAP_RATIO_SUMSETS.md#a-single-ratio-and-an-observable-separable-by-letters)
· [Frey invariants](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md)

## Private denominator primes and support records

Use the limiting one-third estimate and (EB5) for prescribed
pools. For each fixed \(\rho>0\), tuples whose shared input
radical (the product of the primes on at least two rows)
exceeds \(n^\rho\) form a power-smaller class. Indeed some
pair has shared normalized mass at least
\(\rho/\binom{k}{2}\). Fix that row of the pair first and
discount its support on the other row by (EB5).
Count all other slots using the one-third bound with
arbitrarily small allowances. The same argument deletes
tuples having mass more than \(\rho\) on any one row in a
fixed polynomial-height pool. Include in this pool \(P_*\),
all coefficient denominator and numerator primes, shift
denominators, and the pool \(P\) from (NR2). Call the
enlarged pool \(P_0\). Its product has polynomial height.
Also discard rows of raw radical mass below \(1-\rho\) when
needed, using the one-third bound at a reduced budget.
[One-third effective counting and fixed prime pools](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

For the remaining tuples let \(J\) be the shared input
radical. Its support can be recorded for at most \(n^\rho\)
possibilities. Outside \(P_0\) and \(J\), a prime on the
input can belong to just one slot. In each shifted sum
it occurs to exactly the exponent contributed by that
one \(D_i\), since the weight is a unit and all other
terms, including the shift, are integral. In particular
it never divides a reduced shifted numerator. For any
fixed shifted reduced denominator, allocate these private
prime powers to their slots at divisor cost \(X^{o_X(1)}\).
The full parts of every original complete letter on
\(P_0\cup\operatorname{supp}(J)\) are fixed-support choices
at height-subpower cost. Each outside denominator power,
once its role is specified, gives the original prime
valuation by (NR1); role allocations are also subpower.
Each possible denominator thus determines at most
\(n^\rho X^{o_X(1)}\) input tuples. All counts are uniform
in the assigned rational numerators whenever (NR1) holds.
[Uniform enumeration on fixed supports](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)


## Reconstruction from one shifted numerator

For (NR2) use the uniform compact-band estimate at the
one-third endpoint. Outside a power-smaller one-row class,
each row has at least a fixed mass \(\gamma>0\) between
two fixed positive scaled orders \(r_0,R_0\). It suffices
here to have this band on any one role per row. The
constants and the loss are uniform at sufficiently small
fixed budgets.
[Uniform compact valuation mass](ABC_RECIPROCAL_POWERS_VARIABLE_TARGET.md#compact-positive-valuation-mass)

Fix a large \(R>R_0\). For each slot record the whole
product of original prime powers with \(Lv_p(W)>R\).
Each of these tail values has radical at most \(n^{3/R}\)
since \(abc\le X^3\). Recording all tails costs at most
\(n^{3k/R}X^{o_X(1)}\). Add their prime supports to the
pool \(P_0\cup\operatorname{supp}(J)\), giving \(P_1\).
Let \(d_i\) be the part of \(D_i\) outside \(P_1\).
The \(d_i\) are pairwise coprime and, for the shifted
denominator, the identity is
\[
                         D_s=D_0\prod_i d_i,
\]
where \(D_0\) is supported on \(P_1\). After fixing the
tails and \(J\), one can record \(D_0\) at subpower cost,
using the polynomial bound on the whole output height.

Let \(\sigma_i\) be the numerical weight in \(D_i\) of the
primes outside \(P_0\) of original order at most \(R\),
in base \(X\). By (NR1) and the compact band, each
\(\sigma_i\ge h_0:=r_0\gamma/2>0\) once \(\rho\) is small.
Every low-order prime deleted by a tail record from another
slot was already on \(J\). Consequently, if \(d_*=\max d_{i,W}\),
\[
                 \sigma_i-d_*R\rho\ \le\ \log_X d_i
                                            \ \le\ \sigma_i.
                                                               \tag{NR5}
\]
Use a fixed finite mesh to record the \(\sigma_i\). Put the
\(m=\lfloor(k-1)/2\rfloor\) smallest \(\sigma_i\) at the end
of the order, breaking ties by the slot index, and expose
all \(h=k-m\) preceding rows. The cases for this choice of
slots are finite. The sum of their \(\sigma_i\) is at least
the sum on the remaining slots plus \(h_0\). Choose \(\rho\)
and the mesh, after \(R\), so that the loss from (NR5) and
from the upper ends of these fixed meshes is smaller than
\(h_0/2\). The modulus
\[
                         Q=\prod_{i\le h}d_i
\]
then exceeds the upper bound for the unknown product
\(V=\prod_{j>h}d_j\), uniformly among candidate tuples
with the same records. The positive margin is measured
in base \(X\). This comparison also gives a conditional
version with any prescribed block of exposed slots whose
\(\sigma_i\)-sum beats that of the complementary block
by a fixed positive numerical margin.

The numerator \(N_s\) is nonzero since \(D_s\) has private
prime factors. Record \(N_s\) under the asserted radical cap
using support enumeration, at cost at most
\(n^{m/3-\eta+o_X(1)}\). At any prime \(p^e\parallel d_i\)
of an exposed row, write \(f_i=A_i/D_i\) in reduced terms.
Multiplying the shifted sum by \(D_s\) gives the local
congruence
\[
             N_s\equiv w_i A_i\,\frac{D_s}{D_i}\pmod{p^e}.
\]
The other terms vanish, since their expanded denominators
are \(p\)-units. In this equation the coefficient of \(V\)
in the right-hand side is a known \(p\)-unit. In fact
\(D_s/D_i=(D_0/(D_i/d_i))(\prod_{\ell\le h,\ell\ne i}d_\ell)V\),
where each factor besides \(V\) is fixed and \(p\)-integral
with \(p\)-integral inverse. All these congruences determine
\(V\) modulo \(Q\). Since \(0<V<Q\), \(V\) is recovered
exactly. Allocate the prime powers of \(V\), and then
fill the supported parts, as in the preceding denominator
argument. This gives subpower ambiguity after the known
records and the exposed rows.

The exposed rows cost \(n^{h/3}\) up to arbitrary fixed
allowances. Adding the numerator, \(J\), and tail costs
gives an exponent at most
\[
                     k/3-\eta+\rho+3k/R
\]
before arbitrary small allowances. Choose \(R\) first
large in terms of the desired loss, then \(\rho\) and
the mesh in terms of \(R,h_0,\eta\). The discarded bands
and overlap classes have fixed positive losses at
these choices. This proves (NR2), taking into account
the removed classes. On the screened tuples the input
support of \(N_s\) outside \(P_0\) is contained in \(J\).
Thus adding the input support to the excluded primes
in (NR2) costs at most \(\rho\) in the numerator record
and gives the asserted fresh version.

## Two shifted numerators and the output triple

In the setting of (NR3) there is only a fixed finite
pool \(P_0\). Each slot contributes its outside-private
radical to \(D_s\). After discarding rows of low total raw
mass, the input denominator mass is at least
\(k-k\rho\) in base \(n\), prior to removing \(P_0\) and
the shared radical \(J\). These removals cost
\(O_k(\rho)+o_X(1)\). This proves the denominator
part of (NR3). Also \(D_s\ge X^{h_0'}\) for some fixed
\(h_0'>0\) from a compact band on some row after removing
the fixed primes and shared support. No high tails
need be recorded for this part. Take \(\rho\) sufficiently
small in terms of the band mass.

Write \(t=N/D\) in reduced terms and clear the denominators
of two fixed shifts \(s_j,s_l\) with a positive constant
integer \(B\). Use the two forms \(U_j=BN+Bs_jD\),
\(U_l=BN+Bs_lD\) and their difference
\(\Delta D\), where \(\Delta=B(s_j-s_l)\) is a fixed
nonzero integer. Divide by their common gcd and reorder
the absolute values to get a positive primitive triple.
No value vanishes, since each shifted denominator contains
private powers. The gcd divides \(B\Delta\), by
\(\gcd(N,D)=1\), and the output height is at least a fixed
multiple of \(D\), hence exceeds a fixed power of \(X\).
It is bounded above by a fixed power of \(X\) by the input
coefficient-height bounds.

The primes in the unreduced \(U_j\) absent from the reduced
numerator \(N_{s_j}\) lie in the fixed support of \(B\).
Indeed at a prime of \(D\) outside \(B\), \(U_j\) is a
unit. The same holds for \(l\). The reduced numerators
have input primes only in \(J\) and the fixed pool, so
their total raw mass is at most their total fresh mass
plus \(2\rho+o_X(1)\). The output role from \(\Delta D\)
has a bounded normalized raw cap because it is supported
on the input letters and a fixed prime set. All output
radical caps are thus \(O(L)\) in height units in the
sector where the hypothesized two fresh masses have
bounded sum.

Prescribe dyadic height shells for the output. At each
fixed initial budget, the relative budget can be rounded
up on a finite mesh, since the output height is between
two fixed powers of \(X\). Apply (GD3), using the two
shifted numerator roles for the two charged caps. The
number of primitive outputs whose two fresh quantities
in (NR3) sum to at most \(b_0\) is at most
\(n^{b_0/2+\rho+o(1)}X^\varepsilon\) up to arbitrarily
small fixed mesh and budget allowances. From the output
and a bounded choice of signs and roles, the ratio
\(U_j/U_l=(t+s_j)/(t+s_l)\) is known. Inverting this
fractional linear function gives \(t\). Knowing its
denominator gives at most \(n^\rho X^{o_X(1)}\) tuples
by the private-prime recovery. The original tuple
count is at most
\[
                              n^{b_0/2+2\rho+o(1)}X^\varepsilon.
\]
Take \(b_0=2k/3-\eta\) and choose \(\rho\) small enough
to get the strict loss. There are only a fixed number
of pairs of shifts, so the cuts hold simultaneously.
For \(0,-1\) the three original output parts \(N,D,N-D\)
are pairwise coprime. Their bounds give (NR4).
[The two-cheapest-letter inequality](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#group-dictionary-count)

For finitely many distinct shifts, their numerator
radicals overlap only on a fixed set of primes.
Outside denominators of the coefficients, a common
numerator prime would force \(s_j-s_l\equiv0\) at that
prime. Summing (NR3) over the pairs then gives a lower
bound \(bk/3-O_b(\eta)\) for the normalized radical of
their union outside the input support. The numerator
recovery (NR2) also controls each individual shift.
These are exclusions for small-radical **sum outputs**
of regular input tuples. They require no assumption
on the size or distribution of the original complete
dictionaries. To make an unrestricted abc saving one
would need a lower bound for the frequency of such
small-radical outputs or a different way of using the
large output radicals. The graph structure alone does
not supply that lower bound.
[Projective reciprocal identities in dense primitive graphs](ABC_PROJECTIVE_RECIPROCAL_GRAPH_OBSTRUCTION.md)


## Fixed complete fibres and powerful targets

The one-row count used in the reconstruction can be the
fixed-letter count instead of the unrestricted one-third
count. Fix a complete letter \(T_0\) of radical at most
\(n^C\) and put \(B=\ell_U+\ell_V>0\), where the two raw
free-letter caps are \(n^{\ell_U},n^{\ell_V}\).
Use \(e=B/4\). The results below also hold on a fixed
actual powerful target \(A+B'=T_0\), \(A\le B'\),
with orders \(p,q\) on \(A,B'\) and no gcd assumption.
In the latter case use \(X=T_0\), \(L=1/p+1/q\),
\(e=1/4\), \(B=1\), and assume
\(\operatorname{rad}(T_0)\le n^C\). The rational outputs
are formed on the primitive quotients; in particular
\(F(A/T_0)\) has this form.

In both settings take the prescribed pool for freshness
to contain \(\operatorname{supp}(T_0)\). Permit an arbitrary
additional polynomial-height pool in the single-numerator
result. For \(k\ge3\), the analogue of (NR2) has threshold
\(\lfloor(k-1)/2\rfloor e-\eta\) and exceptional count
\(O(X^{(ke-\delta)L+\varepsilon})\). For bounded rational
weights and two distinct fixed rational shifts, the
analogue of (NR3) is
\[
 \begin{split}
 \log_n\operatorname{rad}_{\operatorname{supp}(T_0)^c}(D_s)
                            &\ge kB-\eta,\\
 \sum_{j=1}^2
   \log_n\operatorname{rad}_{(\operatorname{supp}(T_0)
                     \cup\mathcal I(\mathbf x))^c}(|N_{s_j}|)
                            &\ge 2ke-C-\eta .
 \end{split}                                                   \tag{NR6}
\]
In the second line the right side is useful when
\(2ke>C+\eta\). If the two numerator radicals are taken outside
\(\mathcal I(\mathbf x)\setminus\operatorname{supp}(T_0)\),
so that target primes are retained, the \(C\) term can be omitted.
The denominator hypothesis (NR1) remains in force and,
for the two shifts, the exceptional and coefficient primes
beyond the target support are in a fixed finite set.
The constants allow dependence on the fixed caps and order.
A free letter belonging to the original role opposite to
the prescribed \(T_0\) is still a free role in this count.

To justify the change of input price, (FL1) bounds the
normalized free-letter exponent by \(e\). An overlap of
mass at least \(\rho\) on a prescribed pool in the two
free roles replaces this bound by \(e-\rho/4\), by summing
the two outside caps. Thus the overlap screening works on
free letters outside \(\operatorname{supp}(T_0)\). The uniform
compact-band estimate for those parts follows from (FT1)
with that prescribed pool. This is the same height-first
failure diagonal as for (VP3): a failure along bands
expanding to \((0,\infty)\) would produce a zero-order
profile at positive exponent. The numerator reconstruction
therefore gives \(he+\ell+\rho+3k/R\) for numerator cap
\(n^\ell\), which proves the claimed loss.

On the actual powerful target, first prescribe the gcd
\(g_i\mid T_0\) in each slot, at divisor cost. The complete
quotient target is \(W_i=T_0/g_i\). Free prime powers outside
\(T_0\) retain their original powerful orders. If
\(W_i<X^{1/8}\), enumerate just the outside radical of the
quotient free letter with smaller \(1/p,1/q\), and then
its exact supported values. Fixing \(W_i\) and that letter
recovers the row. Including the original target pool in
this fixed-support enumeration costs subpower. This class
has exponent at most \(1/16\) in base \(n\) and can be
deleted. For all other quotient targets use the budget
\(8L\) in the shell of \(W_i\). Converting the two outside
caps back to base \(n\) gives the same bound \(e-\rho/4\).
At exponent tending to \(e\) in base \(n\) the exponent in
base \(W_i^{8L}\) is bounded below, so (FT1) on the outside
parts again supplies the uniform compact band. The
conversion of scaled orders is by the fixed factor eight.
Recovered quotient rows give their original representations
at the prescribed actual target.
[Fixed-letter count after prescribed prime pools](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)
· [Outside-pool nonatomic mass](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#nonatomic-mass-outside-prescribed-prime-pools)
· [Gcd splitting and the short quotient bound](ABC_PROJECTIVE_RECIPROCAL_GENUS_TWO_RECOVERY.md#scalar-relations-on-powerful-representations-of-a-fixed-target)

For (NR6), delete also outside-target free support below
\(B-\rho\), which has a uniform count loss by (FL1).
The private denominator then contains \(kB-O_k(\rho)\)
outside-target mass. The two-shift argument uses (GD3)
on the actual output with total raw cap \(O((kB+C)L)\).
Its proof of the lower bound \(2ke-\eta\) on the sum of
the two numerator masses outside the input letters but
including the target pool may be read as follows. The
raw mass on primes in the input free letters outside
the target pool contributes only \(O_k(\rho)\), because
private numerator primes are excluded on those inputs.
Thus the two-numerator price at most half their remaining
combined mass forces the lower bound. Across the two
shifted numerators, primes of \(T_0\) can contribute
at most \(C+o_X(1)\) in total. A prime outside the fixed
shift exceptional set cannot divide both, since their
unreduced forms differ by a fixed nonzero multiple of the
denominator and the original fraction is reduced.
Removing the target support therefore gives the second
line of (NR6). This keeps the cost of the entire target
pool independent of the number of released slots.


## Half price for fresh radicals of fixed-degree norm forms

A one-row analogue is available when the numerator comes
from an irreducible binary norm form. Fix a finite collection
of irreducible homogeneous forms \(G_1,\ldots,G_v\) in
\(\mathbb Z[A,B]\), each of fixed degree at least two, with
no two proportional over \(\mathbb Q\). Their coefficients
may vary subject to a fixed polynomial-height bound, and
are fixed during each count. Let \(K_j\) be the field of a
root of \(G_j(Z,1)\), and suppose
\[
                     |\operatorname{Disc}(K_j)|
                                        \ll X^{g_j L}
\]
with each \(g_j\) a fixed nonnegative parameter. Allow an
arbitrary polynomial-height prescribed pool \(P\), and put
\(\mathcal I=\operatorname{supp}(abc)\). On the same
unrestricted input cap \(\operatorname{rad}(abc)\le n\),
for any fixed \(\eta>0\), all but
\(O(X^{(1/3-\delta)L+\varepsilon})\) rows at sufficiently
small fixed \(L\) satisfy
\[
 \log_n\operatorname{rad}_{(P\cup\mathcal I)^c}
            \left(\prod_j |G_j(a,b)|\right)
                \ \ge\ \sum_{j=1}^v(2/3-g_j/2)_+-\eta .
                                                               \tag{NR7}
\]
The gaps and thresholds can depend on the fixed degree,
height and discriminant-cap parameters and on \(\eta\).
For fixed forms one may take \(g_j=0\).

To see this, enlarge the pool by all form coefficients
that are excluded in (PR1), all pairwise resultants
between the \(G_j\), and the resultants with each of
\(A,B,A+B\). These integers are nonzero and of
polynomial height. Outside this fixed pool every
evaluated form support is disjoint from the other
forms and from the three original letters, by
primitivity and homogeneous resultant identities.
One may prove the asserted lower bound with this
larger excluded pool. It suffices to show the
one-form bound with an arbitrarily small fixed
margin below \(2/3-g_j/2\) when this is positive.

Suppose the corresponding small-radical class has a
count exponent tending to \(1/3\) at budgets tending
to zero. Use exact valuation profiles with the height
limit first. Effective prime pricing discards classes
placing fixed positive raw overlap in the enlarged
pool, so one may take the limiting profiles outside
it. By endpoint pairwise diffusion, some original
letter contains positive nonatomic mass there.
For any small fixed \(\mu>0\), a finite compact
positive-order band can supply a subpacket of cost
at most \(\mu\) and numerical weight bounded below by
a positive constant \(\nu\). To see the rounding
claim, choose a positive nonatomic subinterval away
from large atoms, then choose whole prime powers
greedily up to a smaller positive mass inside that
subinterval. On a compact continuity region the
largest indivisible cost is bounded by its largest
remaining atom, which can be made arbitrarily small
before making the choice. Prescribed-pool mass is
removed before this selection. All cutoffs, the lower
bound \(\nu\) and the packet price \(\mu\) are fixed
before taking the inner height subsequences. The
selected powers supply a full divisor \(M\) of an
original letter \(A,B\) or \(A+B\), with no excluded
primes. The projective norm moment (XL1) applies to
\(G_j\) with this modulus. If the outside-pool radical
cap of the form is \(n^\ell\), it gives exponent at
most \(\mu+\ell/2+g_j/4\) with arbitrarily small
allowances. Choose \(\mu\) below the fixed gap arising
from \(\ell<2/3-g_j/2\). This contradicts the assumed
exponent \(1/3\). The finite union over the forms and
their disjoint fresh supports proves (NR7).
[Projective ray class count with a divisor of another form](ABC_CROSS_FORM_PROJECTIVE_RAY_LOCKS.md#a-divisor-of-a-different-form)
· [Endpoint diffusion and fixed root masses](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)

For example \(a^2+ab+b^2\), the norm factor in the
Frey \(j\)-invariant, is irreducible with fixed
Eisenstein field. The numerator form
\(a^2b+(a-b)(a+b)^2\) of \(F(a/c)\) is also
irreducible. Its monic dehomogenization is
\(Z^3+2Z^2-Z-1\), with neither possible rational
root \(1,-1\). Its leading coefficients and
discriminant are fixed. They are coprime as forms
over \(\mathbb Q\) because their degrees differ
and both are irreducible. Their two radical
contributions together have normalized fresh mass
at least \(4/3-\eta\) at the unrestricted endpoint.

