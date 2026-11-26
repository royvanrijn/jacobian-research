# Fixed-letter counts using the remaining radical caps

Fix one complete letter of a primitive abc triple. The two remaining
radical caps give a pointwise bound, including after recording all
powers supported on prescribed prime pools. At the balanced endpoint
this yields a separate positive-scale requirement on each letter in
the sector close to the minimum nonatomic mass.

## A pointwise count with prescribed prime pools

Let \(a+b=c\) be positive and primitive, with \(X/2<c\le X\).
Fix the role and the complete value \(T\) of one letter and write
\(U,V\) for the other two. Let \(L>0\) be fixed before the height
limit. Assume \(\operatorname{rad}(T)\le X^{CL}\) and let \(R_U,R_V\)
be prescribed positive integers at most \(X^B\). For \(W=U,V\)
require
\[
             \prod_{\substack{p\mid W\\p\nmid R_W}}p
                        \le X^{L\ell_W},\qquad 0\le\ell_W\le C .
\]
For each fixed \(B,C,\eta>0\) and sufficiently small fixed \(L\),
the number of triples is at most
\[
     O_{B,C,\eta,L,\varepsilon}
           \bigl(X^{L\min(\ell_U,\ell_V)/2+\eta L+\varepsilon}\bigr).
                                                            \tag{FL1}
\]
The constant is uniform in \(T,R_U,R_V\), in the caps, and in the
choice of the fixed role. The two pools may overlap and may depend
on \(T\).
With empty pools the bound agrees with (GD5). The extra input
here is the separate outside radical cap on each of the two
free letters, uniformly after recording their prescribed powers.
[Fixed-value count (GD5)](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#fixed-complete-values-and-prescribed-least-repaired-targets)

The singleton dimension in (FS1) is fixed. The greedy concentration
step there bounds the number of retained heavy groups independently
of the length of the initial valuation list. Constants for that list
are taken at the fixed-budget height limit.

For the proof, split off the complete prime powers on each prescribed
pool and fix their values by fixed-support counting. This costs
\(X^{\varepsilon'}\) for each letter, for arbitrarily small
\(\varepsilon'>0\), uniformly at polynomial height. The fixed
parts enter the coefficients of the signed shift equation with
shift \(T\). The original two terms remain coprime. Write the
remaining terms in exact valuation layers and use dyadic root
boxes with singleton groups. Their total raw caps on the two
letters are \(L\ell_U,L\ell_V\), up to arbitrary truncation and
height allowances. Concretely, use a fixed anatomic truncation at
each fixed \(L\) and put the powers above it into the packed root
or count their supports before fixing the corrections. The latter
cost is at most \(X^{2/D+\varepsilon'}\) for a degree cutoff \(D\),
since a full factor of order greater than \(D\) has radical at
most \(X^{1/D}\). Choose \(D\) with \(2/D\) smaller than the allotted
positive multiple of \(L\). The corrected terms have fixed
coefficients of height at most \(X\).
[Fixed-support counting](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)
· [Exact layers and truncation](ABC_SECANT_ENDPOINT_TRANSFER.md#the-height-first-setting)

If \(T<X^{1/2}\), it is shorter than the maximum by \(X^{1/3}\)
for all large \(X\). Otherwise choose \(\kappa>0\) with
\(\kappa C<1/8\). The part of \(T\) on degrees less than
\(\lfloor\kappa/L\rfloor\) is at most \(X^{\kappa C}\). The complementary
full divisor is at least \(X^{3/8}\) and has order at least
\(\lfloor\kappa/L\rfloor\). In either case (FS1) applies with fixed
margins. Choose its coefficient allowance and the truncation costs
to sum to less than \(\eta L\) and absorb the remaining inner
allowances in \(\varepsilon\). For uniformity in the raw caps one
may round both to a finite grid of sufficiently small mesh
depending on \(\eta\).
[Fixed shift with a full divisor or height gap (FS1)](ABC_FULL_SHIFT_DICTIONARIES_AND_CODEGREES.md#a-fixed-shift-count-without-a-radical-cap-on-the-shift)

In particular, with total radical cap \(X^L\) and a fixed value
\(T\), the exponent in (FL1) is at most
\((L-\log_X\operatorname{rad} T)/4+\eta L+\varepsilon\).
With individual caps \(X^{L(1/3+o(1))}\), the normalized bound is
\(1/6+o(1)\). A dictionary of complete \(T\)'s can be counted by
summing (FL1) with its actual fibre caps.

There is also a support consequence at half-budget saturation.
Suppose fixed-letter fibres have limiting exponent \(d>0\) in
base \(X^L\), and a free letter \(W\) has limiting raw radical
cap \(2d\). For any fixed \(\xi>0\), the subfamily with a
normalized radical mass at least \(\xi\) on a prescribed pool
has exponent at most \(d-\xi/2+o(1)\). This follows by applying
(FL1) with the outside cap reduced by \(\xi\). It is uniform in
pools of bounded polynomial height and holds in every subfibre
with the same limiting exponent.

## Strict positive-scale inequalities on the two free letters

Work with exact valuation measures
\[
 \nu_{W,X,L}=\sum_{p^v\parallel W}
                \frac{\log p}{L\log X}\delta_{Lv},\qquad W=U,V .
\]
Consider fixed-letter fibres in height shells as above with all
three radical caps \(X^{O(L)}\) and normalized count exponent
limiting to \(d>0\). As usual, first take arbitrarily large
heights at each fixed budget and then let \(L\to0\). Refine by
anatomic boxes and suppose the two free-letter measures converge
weakly to \(\nu_U,\nu_V\). Write \(c_W\) for the nonatomic mass
and \(r_W=\nu_W((0,\infty))\). Then
\[
                   d\le\min(c_U,c_V),\qquad
                   d<\min(r_U,r_V).                             \tag{FL2}
\]
This controls each free letter separately.

The full-factor operations can be performed on exact complete
prime powers. One can use the truncated anatomic representations
with normalized error tending to zero, fixing the omitted
high-order support at vanishing normalized cost. At each fixed
budget, root and letter dyadic bins cost a height subpower.
Thus passing to these boxes does not change the exponents or
weak limits in (FL2). All subsequent choices of windows and
positive numerical margins are made at the limiting profile and
held fixed in the height limit.

First suppose just one free letter has positive nonatomic mass.
It supplies a positive numerical full modulus at arbitrarily
small root cost, by selecting a sufficiently small compact
positive continuity band and fixing its powered factors. The
opposite free letter lies in a prescribed unit residue class,
since \(T\) is fixed and the solutions are primitive. Free
finitely many positive atomic windows and the zero window on
that other letter and freeze all remaining positions at
arbitrarily small cost. The grouped roots and the cofactor
with vanishing numerical weight are counted by (GW), with a
single exact complete value determining the triple. This gives
\(d=0\), impossible.

If both nonatomic masses vanish, use the same atom grouping
on the two free letters and use the full divisor of \(T\)
from the proof of (FL1) when \(T\ge X^{1/2}\). The unit ratio
congruence modulo that divisor costs an arbitrarily small
normalized exponent by (GW). When \(T<X^{1/2}\), both free
letters are the tall terms and their ratio is in a real strip
of fixed positive relative decay. The grouped ratio costs an
arbitrarily small exponent by (RUW1). Each exact ratio
determines the triple by primitivity. Again \(d=0\).
Finally, if both nonatomic masses are positive, use an
arbitrarily cheap positive-weight lock on \(U\), free the
atoms on \(V\), and freeze the remaining positions on \(V\)
at cost \(c_V+o(1)\) in the limit. The symmetric argument
gives \(d\le\min(c_U,c_V)\).
[Atom grouping with a cheap full modulus (GE1)](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)
· [Real strip count (RUW1)](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-real-strip-with-unequal-root-heights)

For strictness suppose \(d=r_V\), so \(c_V=r_V=d\). The
entire positive-scale mass of \(V\) is nonatomic. Choose a
band \([b,2b]\), \(b>0\), carrying positive nonatomic mass of
\(U\). On \(V\) choose a positive compact continuity
interval of width at most \(\sigma\) carrying small mass
\(t>0\). Such intervals with arbitrarily small \(t\) exist
with arbitrarily small \(\sigma\). Leave that interval and
a sufficiently small zero window free and fix everything
else on \(V\), at normalized cost at most \(d-t+o(1)\).
The remaining root for the interval has normalized height
at most \(t\); the wider rational cofactor has numerical
height at most \(\sigma t+o(1)\) after shrinking the zero
window and grouping at a representative scaled degree.

In the band on \(U\) select continuity intervals of total
mass \(\mu\) between \(t/64\) and \(t/32\) and freeze their
complete powers. The choice works for sufficiently small
\(t\), by regularity of the measure and by restricting the
nonatomic part away from the countable atoms. The modulus
has weight \(w\ge b t/64\) and fullness proportional to
\(b/L\). The fixed \(T\) makes \(V\) a prescribed unit
residue. Apply the one-root patch estimate (UW) with narrow
height bound a fixed multiple of \(t\). The ratio of the
wider height to the modulus height can be bounded by
\(64\sigma/b+o(1)\). The one-root estimate and its
amplification give cost at most \(t/4+o(1)\) for sufficiently
small \(\sigma\) in terms of the chosen \(U\)-band. This is
the one-patch estimate used in the proof of (IP1), with the
fixed \(T\) in place of the frozen opposite variable. Its
count applies to the grouped root and cofactor; their exact
value determines \(V\) and thus the whole triple.
The total exponent is at most \(d-t+t/32+t/4+o(1)<d\).
Interchanging \(U,V\) proves the second inequality of (FL2).
[One-root wide-cofactor estimate (UW)](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#a-product-bound-for-the-admissible-wider-height)
· [Choice of patch and modulus, including the cost constants](ABC_IMBALANCED_PRODUCT_RADICAL_GROWTH.md#freeing-atoms-with-an-untracked-cofactor)

## Consequence for each letter close to minimum nonatomic mass

There exist absolute \(s,\gamma_0,\delta>0\) such that every
balanced endpoint with \(0\le\Gamma-1/2\le\gamma_0\) has
\[
                   \nu_W((s,\infty))>1/6+\delta
                            \qquad(W=A,B,C).                    \tag{FL3}
\]
Here balanced means limiting count exponent \(1/3\) and
each individual radical mass \(1/3\), and \(\Gamma\) is
the sum of the three nonatomic masses. In particular
the mass at zero on each letter is at most \(1/6-\delta\).

To prove the uniform assertion, suppose it fails and
choose profiles with \(\Gamma-1/2\to0\), some \(s_j\to0\),
and \(\nu_W((s_j,\infty))\le1/6+o(1)\) on a fixed
letter \(W\). Choose the profiles sufficiently close to
minimum for the complete-label reduction with tolerances
tending to zero. Each profile is realized in arbitrarily
small fixed budgets. In the \(j\)-th such profile take a
sufficiently small fixed budget, using a continuity cutoff
between \(s_j\) and \(2s_j\) for the small positive-scale
mass bound. Use the complete-label reduction to retain
exponent \(1/3-o(1)\) with each complete part of exponent
at most \(1/6+o(1)\), including any opposite letter \(T\).
Pigeonhole its complete value. This gives fixed-letter
fibres of exponent at least \(1/6-o(1)\) with \(W\) free.
The prescribed valuation measures are unchanged up to the
box tolerances because all selected rows are in the original
boxes. Compactness gives a limiting fibre with \(r_W\le1/6\):
at each fixed positive continuity cutoff its mass above
the cutoff is at most \(1/6\), by the bounds with cutoffs
tending to zero. Its exponent is \(d=1/6\) by (FL1),
contrary to (FL2). This proof chooses the height allowances
at each fixed budget before taking the diagonal over
budgets. It requires no bound on the number of atoms.
[Complete labels with vanishing count loss near minimum mass](ABC_DIFFUSE_FACTOR_GRAPH.md#the-complete-label-reduction-close-to-minimum-mass)

The same fixed-fibre argument has a dictionary formulation
away from the near-minimum sector. In any limiting family of
exponent \(D\) with a complete \(T\)-dictionary of exponent
at most \(K<D\), pigeonholing \(T\) gives
\(D-K\le\min(c_U,c_V)\) and \(D-K<\min(r_U,r_V)\).
Obtaining useful complete-letter dictionary caps at
arbitrary diffuse mass remains an additional obligation.

## A height gap or escaped weight on a free letter

Continue with (FL2), and take the limiting numerical heights
\(h_U,h_V,h_T\) in base \(X\). Let
\(e_U=h_U-\int r\,d\nu_U(r)\) be the escaped weight on a free
letter \(U\). For an absolute \(\kappa>0\) when the raw caps
are bounded by a fixed absolute constant,
\[
 d\le c_V(1-\kappa e_U).
 \tag{FL4}
\]
Interchange \(U,V\) for the analogous bound. If \(T,V\) are
the two tall letters and \(h_U<1\), there is also
\(d\le c_V(1-\kappa(1-h_U))\).

For both statements use the variable-mesh construction (VM3)
on the single letter \(V\). Choose a compact interval carrying
at least \(c_V/2\) of its nonatomic mass. The proof of (VM3),
with \(M=\lceil4/\theta\rceil\) and just this letter, gives
\(t\ge c_V\theta/10\), a freeze of cost at most
\(c_V-t+t/16\), and a remaining factor \(v z^q\) with
\[
 H(z)\le X^{(C/b+o(1))L},\qquad
 H(v)\le X^{C\theta bt+o(1)} .
\]
Here \(b>0\) is the chosen mesh, \(0<\theta\le1/8\), and
\(C\) is absolute after bounding the initial caps.
The simultaneous approximation of finitely many atom scales
and the finite cover by mesh windows use just the measure of \(V\).

If \(e_U>0\), choose \(\theta\) a sufficiently small absolute
multiple of \(e_U\). After fixing the resulting \(t,b\), lock
the complete high tail of \(U\) above a fixed continuity cutoff
\(R>128/t\). Its radical cost is at most \(1/R\), and its numerical
weight \(w\) is at least \(e_U\). It is \(k\)-full for a choice
\(kL\to8w/t\), since \(w\le1\). With \(T\) fixed, this modulus
makes \(V\) a prescribed unit residue. In (UW) the normalized
narrow height constant is \(O(1/(bt))\) and the relative wider
height is \(O(\theta bt/w)\). Their product is sufficiently
small because \(\theta\ll e_U\le w\). The exponent of the count
is at most \(t/4+o(1)\) in base \(X^L\). Together with the two
freezes this gives \(d\le c_V-t/2\), and hence (FL4).
When \(e_U=0\), use (FL2).

For the real bound, \(V/T=1+O(X^{-\tau})\) for each fixed
\(0<\tau<1-h_U\). Divide \(T\) by the frozen factor and apply
(RUW1) with \(s=t/8\) and \(\theta\) a sufficiently small
absolute multiple of \(\tau\). The narrow constant \(O(1/b)\)
and wider height \(O(\theta bt)\) meet its wider-height
condition, and its count costs \(t/4+o(1)\). This gives
\(d\le c_V-t/2\) again. Let \(\tau\) increase to \(1-h_U\).
In both arguments one exact \(V\) and the fixed \(T\) recover
the original triple.
[Variable-mesh construction (VM3)](ABC_VARIABLE_MESH_DEGREE_GROUPING.md#one-grouped-root-on-the-free-pair)
· [One-root wide-cofactor estimate (UW)](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#a-product-bound-for-the-admissible-wider-height)
· [Real strip count (RUW1)](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-real-strip-with-unequal-root-heights)

## Size of a minimally repaired target-radical fibre

Fix \(0<a_-\le a_+\) and consider positive primitive
\(a+b=c\) with \(\operatorname{rad}(W)\in[c^{a_-\lambda},c^{a_+\lambda}]\)
for each letter \(W\). For a fixed order \(M\ge2\), let
\(H_M(c)=\prod_{p\mid c}p^{\max(M,v_p(c))}\), and consider
the least-target-repaired triples of height \(H_M(c)\le Y\).
Let \(m_Y(R)\) count their primitive quotients with
\(\operatorname{rad}(c)=R\). For each \(\zeta>0\), all sufficiently
small fixed \(\lambda\), and every \(\varepsilon>0\), one has
\[
 \begin{split}
 m_Y(R)&\ll_{\lambda,a_\pm,\zeta,\varepsilon}
             R^{a_+/(2a_-)+\zeta+\varepsilon},\\
 \max_R m_Y(R)&\ll_{\lambda,a_\pm,\zeta,\varepsilon}
             Y^{(a_+/(2a_-)+\zeta)/M+\varepsilon}.
 \end{split}                                                     \tag{FL5}
\]
These estimates are uniform in \(R,Y,M\). Indeed all such
\(c\) have support \(R\) and \(c\le R^{1/(a_-\lambda)}\).
There are \(R^\varepsilon\) complete choices uniformly, by
fixed-support enumeration, with the inner allowance chosen
depending on \(\lambda,a_-\). For each \(c\), take \(X=c\)
in (FL1) with \(T=c\), \(L=\lambda\) and both free caps
\(\ell=a_+\). Choose the loss below \(\zeta a_-\lambda\)
in the primitive height exponent. This proves the first bound.
Then \(R^M\le H_M(c)\) proves the second. The argument also
works when the selected repair letter is a summand, by taking
the primitive maximum as \(X\) for each solution and using
height shells and the small-\(T\) alternative of (FL1).
This counts each primitive quotient once, as least repair
is prescribed by that quotient.

In particular, at count exponent approaching \(1/M\) in
repaired height and as \(a_-,a_+\to1/3\), the occurring target
radicals have exponent at least \((1/2-o(1))/M\). A target
radical dictionary of exponent bounded below \(1/(2M)\)
by a fixed positive proportion therefore gives a strict
saving in that balanced repair sector. This sharpens the
positive marginal exponent previously obtained by the
fixed-shift gap. It does not estimate total occurrence when
both the target marginal and its representation fibres are large.
[The weighted target-support count and the earlier marginal restriction](ABC_MINIMAL_REPAIR_SUPPORT_OCCURRENCE.md#marginal-complexity-forced-by-the-fixed-shift-gap)

### The remaining multiplicity estimate

Put \(Z=Y^{1/M}\) and \(b=a_+/(2a_-)+\zeta\). In this
normalization (FL5) gives \(\max_R m_Y(R)\le Z^{b+o_Y(1)}\)
and \(R\le Z\), with \(M,\lambda,a_\pm\) fixed first.
The joint estimate sufficient for a saving is, for some
absolute \(\sigma>0\) and all \(u\in[0,b+\sigma]\),
\[
         \#\{R:m_Y(R)\ge Z^u\}
                     \ll Z^{1-u-\sigma+\varepsilon}.            \tag{FL8}
\]
It suffices to prove this for a mesh of spacing less than
\(\sigma/4\) in \(u\), including \(u=0\), taking errors
smaller than the saving. The weighted count is bounded by
summing over these finitely many bins, using the next mesh
point for the multiplicity cap in each bin. Conversely a
count bound \(Z^{1-\sigma+\varepsilon}\) implies (FL8)
directly. The earlier weighted-support equivalence used the
same bins with the full multiplicity range; (FL5) limits
the range needing a joint estimate to \(u\le1/2+o(1)\)
as the balance tolerance is decreased. In particular
(FL8) still needs estimates at intermediate \(u\), even
if the number of distinct target radicals has a positive
power saving.
[Weighted support equivalence](ABC_MINIMAL_REPAIR_SUPPORT_OCCURRENCE.md#exact-weighted-form-on-the-repaired-target-radical-alone)

## Several least repairs in a single multiplicative class

On an original primitive triple, put
\(\mu_M(W)=\prod_{0<v_p(W)<M}p^{M-v_p(W)}\).
Fix arbitrary integer weights \(e_W\), and prescribe the
class of \(\prod_{W}\mu_M(W)^{e_W}\) modulo
\((\mathbb Q_{>0}^{\times})^M\). Use height \(X\) and raw
radical caps \(\ell_W L\) as in (FL1) on all three letters.
For each letter define \(g_W=\gcd(M,e_W)\) (including
\(\gcd(M,0)=M\)) and put \(b_W=\min(\ell_W L,g_W/M)\).
For each fixed allowable budget and allowance, the count in the
prescribed class is bounded by
\[
 O\!\left(
 X^{\min_{T}\,[\,b_T+\frac12\min(b_U,b_V)\,]+\eta L+\varepsilon}
 \right),                                                        \tag{FL6}
\]
uniformly in the repair order, weights and prescribed class.

To prove this, take one original triple in the class, if it
is nonempty, and let \(S\) be the pool of its primes. At a
prime outside \(S\), primitivity gives at most one letter
\(W\) with a nonzero valuation. If that valuation \(v\) is
below \(M\), the class constraint gives \(M\mid e_W v\),
so \(v\ge M/g_W\). If \(v\ge M\), the same lower bound is
automatic. Thus each letter has outside radical cap \(X^{b_W}\).
Its complete dictionary has size at most \(X^{b_W+\varepsilon'}\)
by enumerating the outside radical and then the powers on
the union with \(S\), whose radical is at most \(X^4\).
For each \(T\) in that dictionary, use (FL1) with the original
signed equation and this same prime pool to count the other
two. Optimize the choice of \(T\). This also proves that a
single class on any fixed complete \(T\) costs at most
\(X^{\min(b_U,b_V)/2+\eta L+\varepsilon}\).

In particular, if two weights are units modulo \(M\) and
\(ML\to\infty\), a balanced endpoint requires at least
\(X^{(1/3-o(1))L}\) classes for their weighted product,
even if the third weight is zero. A uniform bound below
that class exponent by a fixed margin would give a saving
at all balanced profiles. Raw radical bounds alone do not
supply such a bound.

## Least common multiples of simultaneous valuation classes

The same class method also works with effective orders on
the fine-grid scale, where (FL6) need not give a gain.
For \(j\) in a finite list choose a modulus \(M_j\ge2\)
and integer weights \(e_{j,W}\). Prescribe one class
modulo rational \(M_j\)-th powers for each of the products
\(\prod_W F_{j}(W)^{e_{j,W}}\), where \(F_{j}(W)=W\)
or \(F_j(W)=\mu_{M_j}(W)\). Use the same formula for all
three letters at a given \(j\). For each letter choose any cutoff \(B_W\in[2,\infty]\).
Retain all complete-value constraints, and retain a repair
constraint for \(W\) only when \(M_j\ge B_W\).
Define \(k_W\) as the least common multiple of
\(M_j/\gcd(M_j,e_{j,W})\) over the retained constraints.
Suppose \(k_W\ge K\ge2\) on all three letters. Among positive
primitive abc triples with an individual radical cap \(X^{f_W}\)
on letter \(W\), the joint class has at most
\[
           O_{K,\varepsilon}\left(
               X^{\,\sum_W \min(f_W,1/B_W)+3K^{-3/2}+\varepsilon}
                               \right) .                         \tag{FL7}
\]
Here \(1/\infty=0\), and all heights are at most \(X\). The bound
is uniform in the list, moduli, weights and prescribed classes.
The caps may be omitted by taking \(f_W=1\). There is no
assumption about the limiting valuation measures.

Use one original triple in the joint class as a reference.
Its support \(S\) has radical at most \(X^3\), independently
of the size of the list. One can optimize the cutoffs separately
on the three letters. For each letter freeze the
\(S\)-supported part at height subpower cost. For primes
outside \(S\), primitivity ensures that a prime dividing \(W\)
contributes only on \(W\) to every prescribed product.
If \(v_p(W)<B_W\), a retained class condition from a least repair has
the same divisibility consequence on \(v_p(W)\) as one
from the complete original value. In fact the exponent
of the repair at \(p\) is \(M_j-v_p(W)\) for each active
repair with \(v_p(W)<B_W\le M_j\). Discarded or inactive terms impose
no condition in the estimate. For a complete original value the class
always gives this condition, at every depth. Thus
\(k_W\mid v_p(W)\) for all the outside primes below the
cutoff on that letter.

Freeze the remaining outside prime powers on \(W\) of
valuation at least \(B_W\). Their total radical is at most
\(X^{\min(f_W,1/B_W)}\), and there are at most that many
choices for their supports. For each support, the full
part has height subpower many possibilities. The combined
frozen parts are the three fixed coefficients of a
generalized Fermat equation in powers of orders
\(k_A,k_B,k_C\). For fixed coefficients the count is at
most \(O_{K,\varepsilon}(X^{3K^{-3/2}+\varepsilon})\), as
in the exact fine-grid sector estimate. That proof is
uniform even when orders are large or vary with height:
bounded orders are treated by Browning and Verzobio's
coefficient-uniform theorem, while higher orders are
enumerated by their actual values. Multiplying by the
coefficient choices proves (FL7).
[Uniform generalized Fermat count at arbitrary orders](ABC_FINE_VALUATION_DIVISIBILITY.md#exact-sector-bound)
· [Browning and Verzobio, Theorem 1.2](https://arxiv.org/html/2608.24512v1#S1)

For example take a fixed \(\alpha>2/3\), \(K=\lceil L^{-\alpha}\rceil\)
and suppose \(LB_W\to\infty\) for each applicable repair cutoff.
Then each joint class has normalized exponent \(o(1)\) as
the fixed budgets tend to zero. Many individually small
effective orders can achieve this \(K\) by their least
common multiple. One could use a list of complete-value
classes only, in which case there is no cutoff cost at
all. It follows that a near-full balanced endpoint has
joint class dictionary exponent at least \(1/3-o(1)\)
for every such list. This conclusion concerns the
primitive quotients. The joint class record of a uniformly chosen triple must
also have Shannon entropy at least \((1/3-o(1))L\log X\),
since the maximum class probability is at most
\(X^{-(1/3-o(1))L}\). Consequently any upper bound below
this threshold for the sum of the individual class
entropies would imply a count saving. The common multiplier
of the minimally repaired triple plays no part in these
class constraints.
