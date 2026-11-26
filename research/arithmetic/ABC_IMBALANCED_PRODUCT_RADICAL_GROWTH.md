# Radical growth and low-order valuations in imbalanced products

Products of ratios from balanced exceptional families gain fresh
radical when the weights are sufficiently unequal. The lower bound
can exceed any prescribed multiple of the original radical budget.
The conclusion also controls the numerical cofactor left after
recording a divisor of bounded radical. The count uses the
primitive output triple and permits the weights to vary with the
budget.

## A degree-measure bound before balance

Consider a sequence of families of positive primitive abc triples
of height in \([Y/2,Y]\). Set \(\ell>0\) before \(Y\to\infty\), then
let \(\ell\to0\). The count has normalized exponent \(d>0\) in
base \(Y^\ell\). On two labelled letters \(P,Q\), track the full
exact valuation measures
\[
     \sum_{p^v\parallel R}\frac{\log p}{\ell\log Y}\delta_{\ell v},
                          \qquad R=P,Q .
\]
On the third letter \(C\), track just a full divisor \(Z\) of \(C\),
where \(W=C/Z\) has numerical height \(Y^{o(1)}\) in this iterated
limit. A full divisor is a product of complete prime powers of
\(C\). Suppose the radicals of the three tracked integers are
\(Y^{O(\ell)}\), uniformly in the small-budget limit. Take weak
limits \(\mu_P,\mu_Q,\mu_Z\) in the corresponding anatomic boxes.
Put
\[
   r_R=\mu_R((0,\infty)),\qquad
   c_R=\text{nonatomic mass of }\mu_R .
\]
Then
\[
                              d<r_P+r_Q.                         \tag{IP1}
\]
The letter \(C\) may have any numerical height. No radical bound
on \(W\) is needed.

To pass from exact valuations to boxes, use the anatomic
truncation separately on \(P,Q,Z\) with parameter \(o(\ell)\).
The initial coefficient lies in the vanishing numerical
correction. The packed coordinate and the omitted exact layers
have total normalized root mass \(o(1)\), and the numerical
heights of the tracked integers are at most \(Y\).
At fixed budget, dyadic root-coordinate boxes and height bins
cost \(Y^{o_Y(1)}\). A choice of \(Z\) per row suffices.
There are only divisor-many possible choices per row if all
presentations are retained. These reductions give bounded
total mass and tightness by the first moment bound.
[Exact valuation layers and truncation](ABC_SECANT_ENDPOINT_TRANSFER.md#the-height-first-setting)

### Freeing atoms with an untracked cofactor

Use the argument of (GE1). If one tracked letter has a fixed
full modulus of positive numerical weight, bought at cost
\(\mu\) in units of \(\ell\), free finitely many shrinking
positive atom windows and the zero window on the other
two tracked integers. Freeze all remaining root positions on
those two at cost tending to the sum of their nonatomic masses.
Grouping a positive window uses one short root with a
representative exponent and a signed power correction of
numerical height \(Y^{o(1)}\). The zero windows have the
same numerical height bound. If the two free letters include
\(C\), include its \(W\) in this wider rational correction.
This does not change the applicable wide-cofactor estimate.
A fixed full divisor of any one letter makes the ratio of
the other two a fixed unit residue. Knowing the exact ratio
recovers the primitive triple. The result is
\[
                d\le \mu+c_{R_1}+c_{R_2},                        \tag{IP2}
\]
where the notation for \(C\) means its tracked divisor \(Z\).
A tracked letter of positive nonatomic mass supplies such
moduli at arbitrarily small positive cost. Choose a compact
positive band meeting its nonatomic part and use small
continuity subintervals avoiding its large atoms. A tracked
letter with positive mass on \((0,\infty)\) also supplies a
positive-weight modulus at cost at most its entire positive
mass, by fixing the root positions of a compact positive
continuity band. The chosen factors may be taken to be
actual complete prime powers. If a packed coordinate is
present, ignoring it only decreases the modulus cost and
does not alter the positive mass in the band.
[Grouping with a cheap full modulus (GE1)](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)

When \(c_Z=0\) and both \(r_P,r_Q>0\), applying (IP2) with
locks on \(P\) and on \(Q\) gives
\[
              d\le r_P+c_Q,\qquad d\le r_Q+c_P.
\]
If \(d\ge r_P+r_Q\), then equality forces \(c_P=r_P>0\)
and \(c_Q=r_Q>0\). Either pure letter now supplies an
arbitrarily cheap modulus. Applying (IP2) with that letter
gives \(d\le c_Q\) or \(d\le c_P\), a contradiction.

For the remaining \(c_Z=0\) case exchange \(P,Q\) if needed
so \(r_Q=0\). If \(Q\) has positive limiting numerical height,
that moment has escaped to high orders. The full part on
orders \(v\ge 1/\ell\) then has numerical weight bounded below
and radical exponent \(o(1)\). Its exact value can be fixed
at normalized cost \(o(1)\) by radical and fixed-support
enumeration. Using it as a modulus in (IP2) gives \(d\le c_P\).
If \(Q\) has zero limiting numerical height, \(P,C\) are the
tall terms and \(P/C=1+O(Y^{-1/2})\). On \(P,Z\) free the
zero windows and finitely many positive atomic windows as
in (IP2), and include \(W\) in the wide correction. Freeze
the remaining positions at cost \(c_P+c_Z+o(1)\). After the
frozen factors are removed, the grouped ratio lies in a
fixed real strip of relative width \(O(Y^{-1/2})\), with
boundedly many short powered roots and a wide factor of
height \(Y^{o(1)}\). By (RUW1) that strip costs \(o(1)\)
in these units. Its exact ratio recovers the primitive
triple, so one may choose one presentation per exact ratio.
If there are no grouped roots, rational spacing in the
strip gives the same bound. Thus again \(d\le c_P\).
Now \(d\ge r_P+r_Q=r_P\) would give \(c_P=r_P=d>0\).
A cheap \(P\) modulus in (IP2) then gives \(d\le c_Q+c_Z=0\).
[Real strip count (RUW1)](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-real-strip-with-unequal-root-heights)

When \(c_Z>0\), arbitrarily cheap locks on \(Z\) give
\(d\le c_P+c_Q\le r_P+r_Q\). Equality forces both positive
parts to be nonatomic, including the case of one zero
part. Since \(d>0\), one of the two parts has positive
mass. On that letter take a compact positive continuity interval
of arbitrarily small width \(\sigma\) and mass \(t>0\)
tending to zero. Retain its positions and both zero windows
on \(P,Q\), and freeze the rest at cost \(d-t+o(1)\).
The patch is a single short root with a rational correction
of numerical logarithm at most \(\sigma t+o(1)\).

Choose a compact positive band \([b,2b]\) of positive
nonatomic mass on \(Z\). Inside it one can choose a finite
union of continuity intervals of total mass between \(t/64\)
and \(t/32\) for sufficiently small \(t\). For example choose
a compact subset of the nonatomic part of mass strictly
between the two thresholds, avoiding the countable atoms
of the full measure, and cover it inside an open set with
mass still below \(t/32\). A finite family of continuity
intervals in that open set suffices. Freeze its full factors
as a modulus of numerical weight \(w\ge bt/64+o(1)\)
and fullness at least \(\lfloor b/(2\ell)\rfloor\).
The one-root wide-cofactor count (UW), with narrow cap
constant \(64\) or any fixed enlargement of it, counts the
remaining ratio at cost at most \(t/4+o(1)\) by taking
\(\sigma\) sufficiently small compared with \(b\).
In fact the ratio of the wider logarithmic height to \(w\)
is at most \(64\sigma/b+o(1)\), and the unamplified saving
can be made this small. This is the patch argument for
(VM7), with \(d\) in place of its balanced value \(1/3\).
The combined cost
\[
                       d-t+t/32+t/4+o(1)<d
\]
is impossible. This proves (IP1). The argument needs the
untracked cofactor only in the case where \(C\) is one of the
two free letters.
[One-root count (UW)](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#a-product-bound-for-the-admissible-wider-height)
· [Zero-atom patch argument in (VM7)](ABC_VARIABLE_MESH_DEGREE_GROUPING.md#atomic-and-nonatomic-weight-near-the-minimum-mass)

## Imbalanced products of independent input directions

Fix an integer \(m\ge2\) and \(B>0\). There exist
\(\rho,\eta,\delta,L_0>0\) with the following property.
For every fixed \(0<L\le L_0\), let \(n=X^L\). In each of \(m\)
slots use the primitive positive abc triples of height in
\([X/2,X]\) with each individual letter radical at most
\(n^{1/3+\delta}\). One may impose additional restrictions in
any slot. Choose an ordered pair \(U_i,V_i\) of distinct letters in
slot \(i\) and positive integer weights \(u_i,v_i\) respectively.
Write \(K=\max_i\max(u_i,v_i)\). Suppose that at least \(m\) of
the \(2m\) weights are at most \(\rho K\).
Form the primitive pair
\[
     A:B=\prod_i U_i^{u_i}:\prod_i V_i^{v_i},\qquad
                       \gcd(A,B)=1.
\]
For the synchronized choice \(u_i=v_i\), the hypothesis holds
if \(\lceil m/2\rceil\) slots carry small weights. Signed
ratio exponents can be implemented by swapping the selected
letters.
Fix nonzero integers \(a,b\). When \(aA+bB\ne0\), divide
\((|aA|,|bB|,|aA+bB|)\) by \(\gcd(aA,bB)\).
Denote the three labelled positive output letters by \(P,Q,C\)
and their maximum by \(Y\). Count only tuples on which \(C\)
has a full divisor \(Z\) with
\[
                   \operatorname{rad}(Z)\le n^B,\qquad
                           C/Z\le Y^\eta .
\]
Then the tuple count is
\[
                       \ll X^{(m/3-\delta)L+\varepsilon}          \tag{IP3}
\]
for every \(\varepsilon>0\). Choices of orientations, weights
and the integers \(a,b\) are held fixed before \(X\to\infty\).
They may vary with \(L\). The exponent margins and budget
threshold are independent of these choices. Implied constants
may depend on them, the fixed budget and \(\varepsilon\).

### Proof and fixed-budget quantifiers

If the assertion fails, choose a sequence of violations with
\(L,\delta,\eta,\rho\to0\). At each of those budgets there are
arbitrarily large heights with near-full tuple exponent
\(m/3-o(1)\). Choose the height allowances there below the
corresponding fixed-budget margins. Each slot has exponent at
most \(1/3+o(1)\) by the global small-budget count, and hence
each projection of the eligible tuple family has exponent
\(1/3-o(1)\). Refining before the height limit by dyadic
anatomic root boxes in all slots costs inner subpower. The
input budget can be padded by \(o(L)\) in this limiting
argument, so each resulting marginal is a genuine balanced
endpoint. Choices of letters and of the indices of the small
weights can be fixed on a subsequence.
[The limiting small-budget count](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)

For each fixed overlap tolerance, (DRR4) discards a
power-smaller set of tuples with excessive common prime
mass. Decrease this tolerance along the budgets. Any prime
exclusive to an input row contributes its selected full
valuation, multiplied by its weight \(u_i\) or \(v_i\), to the corresponding
pure output letter. Shared supports have total normalized
radical \(o(1)\). The reconstruction from a fixed primitive
output pair costs only \(n^{o(1)}\), by recording the full
parts on shared supports and allocating the remaining primes
to slots before dividing their valuations by their weights.
Dividing by constants supported on \(ab\) and recovering the
primitive output pair from the labelled pure output ratio
give only inner losses. These are the support and
multiplicity arguments for (ST2). With varying weights all
integers involved still have polynomial height at each
fixed budget. The allowed polynomial exponent can depend on
that budget. The fixed-prime enumeration and the bound
on the number of allocations are taken in the height limit
there, before the weights vary.
[Small intersections and output multiplicities](ABC_SECANT_ENDPOINT_TRANSFER.md#small-intersections-and-output-multiplicities)

Let \(S=\sum_i(u_i+v_i)\). A dominant weighted factor
contributes a uniformly positive compact positive-scale
packet to its pure letter, by the input endpoint theorem
(SH1). Discarding \(o(1)\) radical at bounded scales removes
only \(o(1)\) of the input numerical weight in this packet.
Its output numerical weight in base \(X\) is therefore a
positive fixed fraction of the largest weight. In particular
if \(H=\log_X Y\) then \(H/S\) is bounded below by a fixed
positive number and above by \(1+o_X(1)\). Choose finite
bins for \(H/S\) at each fixed budget, refining the mesh
along the budgets, and a dyadic output height shell. In a
chosen bin write \(H\) for the numerical exponent of its
shell, with the vanishing relative bin allowance absorbed
in the normalization. Set \(\ell=L/H\).
This is fixed before the inner height limit and tends to
zero, regardless of whether the weights grow. The three
tracked radicals in (IP1) are of exponent \(O(\ell)\) in
base \(Y\), since \(Y^\ell=n\) up to the bin allowances.
There are \(n^{m/3-o(1)}\) distinct primitive labelled
outputs, with \(d\ge m/3\) in these units.
[Compact positive-scale input mass (SH1)](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#positive-scales-in-every-endpoint-letter)

Pass to limits of \(u_i/H,v_i/H\) in their bounded range.
The exact weak law of a pure output is the sum of the
input laws pushed forward by the corresponding dilations,
including dilation by zero. This form of (ST2) follows
by testing the exclusive valuations. Tightness of the
input measures justifies dilation by factors tending to
zero. At most \(m\) among the \(2m\) dilations have
nonzero limits, each contributing at most \(1/3\) of
positive-scale mass. Thus \(r_P+r_Q\le m/3\).
The remaining cofactor \(W=C/Z\) has vanishing numerical
height because the violating sequence has \(\eta\to0\).
Apply (IP1) to the tracked output measures. This gives
\(m/3\le d<r_P+r_Q\le m/3\), a contradiction.
The selection used a finite number of types at each
fixed budget, with an inner subpower allowance, so the
contradiction proves the stated order of quantifiers.

### Fresh radical and numerical weight at low orders

First, the fresh radical means the radical outside all three
complete letters of every input row. For every fixed \(B\),
(under the corresponding constants in (IP3)) the tuples
with fresh radical of \(aA+bB\) at most \(n^B\) are
power-smaller. More generally the theorem permits unrecorded
fresh factors of vanishing numerical height. The radical of the
nonfresh part costs at most \(m+o(1)\) in base \(n\), from all three
original letters of every slot. Primes of the fixed coefficients
have inner exponent zero. Take \(Z\) to be the full nonfresh part
times the full selected fresh divisor, and apply (IP3) with budget
\(B+m+1\). This argument accommodates independent weights on the
two pure products. The omitted-letter congruence in the
synchronized case can lower the nonfresh allowance to \(o(1)\),
but an upper bound by a fixed constant suffices here.

Second, for each fixed \(u>0\) there are corresponding
\(\rho,\eta,\delta,L_0>0\) such that outside a power-smaller
tuple set,
\[
       \sum_{\substack{p^v\parallel C\\vL/H<u}}
                         v\log_Y p>\eta,\qquad H=\log_X Y.        \tag{IP4}
\]
Use \(Z\) equal to the full part of \(C\) on orders \(v\ge uH/L\).
Then
\(\log_n\operatorname{rad}(Z)\le(\log_Y C)/u\le1/u\).
If (IP4) fails, the complementary part is at most \(Y^\eta\).
Apply (IP3) with \(B=1/u\). The same conclusion holds if
the numerical weight in (IP4) is restricted to fresh primes:
absorb the nonfresh full part in \(Z\) at a bounded extra
normalized radical cost as above.
In particular a product of two ratios with weights \((K,1)\)
has a positive amount of low-order fresh numerical weight
for sufficiently large \(K\), in the indicated order of
parameters, with no prior upper bound on its fresh radical.

### An exclusion for original balanced profiles

For original triples in the balanced sector there are absolute
\(\sigma,\theta,\delta,L_0>0\) as follows. In the height shell
\([X/2,X]\), restrict each letter radical to
\(n^{1/3+\delta}\), where \(n=X^L\) and \(0<L\le L_0\)
is fixed before the height limit. Suppose two letters
together carry at least \(1/3-\theta\) normalized radical
mass on primes of orders \(vL\le \sigma\). Then the triple
count is \(O(X^{(1/3-\delta)L+\varepsilon})\).
Indeed, failure would give a limiting endpoint along
\(\theta,\sigma,\delta,L\to0\). By weak convergence the
combined zero atom on those two letters has mass at least
\(1/3\). Their total mass is \(2/3\), so their positive
mass is at most \(1/3\). Apply (IP1) with the full remaining
letter tracked to contradict \(d=1/3\). This applies at any
value of the total nonatomic mass. Both low-order and
positive-scale masses here refer to radical weight.

For synchronized weights there is a further profile-dependent
criterion. Let \(I\) be the slots whose weight has a positive
limiting ratio to \(H\). Along any such output sequence, if
\[
                   \sum_{i\in I}(c_{U_i}+c_{V_i})<m/3 ,
\]
then a bounded-radical part capturing all but vanishing numerical
weight of \(C\) is impossible on a near-full tuple family.
For a diffuse \(Z\), a cheap modulus there gives
\(d\le c_P+c_Q<m/3\). For atomic \(Z\), at least one pure
letter has positive nonatomic mass because \(I\) is nonempty
and every input selected pair has nonatomic total at least
\(1/3\). A cheap modulus on that pure letter then gives the
same contradiction, via (IP2). In particular near the
minimum input nonatomic total, even one vanishing slot weight
is enough for fixed \(m\). Its live pair sum is at most
\((m-1)(1/3+2\gamma)\) if every input total is at most
\(1/2+\gamma\).

### Packing high powers on small prime roots

There is a separate fixed-coefficient estimate with a numerical
threshold independent of the dominant exponent. Let \(K\ge2\)
and \(L>0\) be fixed before the height limit, and let
\(\alpha,\beta\ne0\) be prescribed integers of polynomial height.
For coprime \(1\le u,v\le X\) put \(D=\alpha u^K+\beta v^K\ne0\).
Fix \(\kappa,\tau>0\), and use the full part of \(D\) on primes
\(p\le n^\tau\) with \(p\nmid K\alpha\beta\) and order
\(v_p(D)\ge \kappa K/L\). The number of primitive pairs for
which this full part exceeds \(2X^2\) is
\[
                       \ll n^{\tau+2/(\kappa K)}X^\varepsilon .  \tag{IP5}
\]
Multiply its prime powers until their product \(M\) first exceeds
\(2X^2\). The previous product costs at most
\(2/(\kappa K)+o_X(1)\) in radical units, and the last prime
costs at most \(\tau\). Enumerate the radical and then the
prime powers at inner subpower cost. At every prime dividing
\(M\), primitivity and \(p\nmid\alpha\beta\) imply \(p\nmid uv\).
The congruence \((u/v)^K=-\beta/\alpha\pmod M\) has at most
\(K^{\omega(M)}=X^{o_X(1)}\) solutions in unit residues, by
simple-root Hensel lifting and the Chinese remainder theorem.
Each residue allows just one primitive pair: the determinant
of two such pairs would be divisible by \(M>2X^2\), which
forces its vanishing. This proves (IP5), uniformly in the
prescribed coefficients at fixed polynomial height. In a
two-slot imbalanced product, the minority row can enter these
coefficients. The estimate then shows that sufficiently high
orders on small primes are sparse on the dominant-row fibre.
Primes above the selected cutoff are not part of this estimate.

The ordinary two-row secant uses equal synchronized weights.
An unrestricted endpoint can avoid the small-degree input
sector above. A transfer back to the ordinary secant would
need an additional arithmetic implication.
