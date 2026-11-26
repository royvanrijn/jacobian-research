# Laurent-monomial output rigidity on actual abc boxes

On almost all tuples of synchronized near-full shared boxes, an
additive output of polynomial primitive height given by three
Laurent monomials in the original labels is an original row or
a constant output. The coefficients may come from a dictionary
chosen separately for each ordered residue list. The conclusion
holds with unrestricted integer exponents. A private full prime
power part on each original edge gives a uniform exponent bound
even after removing the numerical gcd of the output.

## Joint private numerical height

Use the ordered \(m\)-boxes of
[LB1](ABC_SHARED_VALUE_RESIDUE_BOXES.md#a-finite-difference-array-theorem),
with \(m\ge2\), residue alphabet bounds \(R_i\), product
\(R=\prod_i R_i\), and pair-class bound \(K\). Every triangle
is an original positive primitive row in \(X/2<C'\le X\).
Every complete edge label has radical at most \(n^C\) for a
fixed \(C>0\), where \(n=X^L\). Let \(h=3m\).
Take an ordered tuple of \(j\ge1\) boxes with the same three
ordered residue lists. The integer \(j\) is fixed. Put
\(u=j(h-1)\) and \(q=3jm^2\), the numbers of independent tree
coordinates and of complete edge positions, respectively.

Fix an external prime pool \(P\) whose prime product is at most
\(X^B\), for a fixed \(B>0\). The pool may be chosen separately
for each ordered residue list. For an edge position \(e\) of
the tuple define its private full prime power part by
\[
 W_{e,P}=\prod_{\substack{p\mid L_e,\ p\notin P\\
                   p\nmid L_f\ \text{for every } f\ne e}}
                   p^{v_p(L_e)} .
\]
The test against the other positions ranges over all \(j\)
boxes, including positions in other roles.

**LP1.** For every fixed \(r>0\) there are
\(\theta,L_0>0\) such that for every fixed \(0<L<L_0\) the
number of ordered synchronized tuples with \(W_{e,P}\le
X^\theta\) at any position is
\[
             O_{m,j,B,C,r,L}\bigl(R^mK^{u-1} n^r\bigr).
                                                               \tag{LP1}
\]
The constants are uniform in the original family, the alphabet
and pair-class bounds, and the pools.

Fix the lists, the other \(j-1\) boxes, and all potentials
of the tested box except one endpoint \(v\) of \(e\). The
other boxes cost at most \(K^{(j-1)(h-1)}\) tree labels.
The induced graph of the remaining \(h-1\) potentials is
connected, so \(h-2\) tree labels determine those anchors
up to translation. Thus there are at most \(R^m K^{u-1}\)
such records. The residue of \(v\) is already fixed.

Enlarge the pool by the complete labels of all other boxes,
the complete labels between two anchors, and all nonzero
differences between the neighbors of \(v\). All these
integers have polynomial height uniformly in the recorded
data. Different neighbors give different potentials, by
distinct residues or positivity of the cross labels. If a
prime divides \(L_e\) and another label incident to \(v\),
it divides the difference between the two corresponding
neighbors. Consequently the part of \(L_e\) outside the
enlarged pool divides \(W_{e,P}\). This is the prime
exclusion in (PV1) with the other boxes included in the
pool.

Complete \(e\) to a triangle using an anchor on the third
axis. The side between its two anchors is a fixed complete
letter of the original family. Freeze the entire part
\(U_0\) of \(L_e\) supported on the enlarged pool. There
are \(O_\varepsilon(X^\varepsilon)\) choices uniformly at
fixed polynomial pool height. When \(W_{e,P}\le X^\theta\),
the remaining factor \(z\) is a positive integer at most
\(X^\theta\), and the moving complete letter is \(U_0z\).
The rank-one small-feature bound at exponent \(r\)
gives \(O(n^{r-\delta})\) possibilities for some
\(\delta>0\), after choosing its \(\theta,L_0\).
The bound allows the affine step \(U_0\) to be arbitrarily
large. Use the supported-part allowance
\(\varepsilon<\delta L\) at each fixed \(L\).
The value on \(e\) recovers the omitted potential.
Sum the bounds over the \(q\) positions. The proof always
applies the fixed-letter bound to original primitive rows
in their original height shell and retains all three
individual radical caps.
[Omission and prime exclusion (PV1)](ABC_SHARED_BOX_PRIVATE_PRIME_SUPPORT.md#a-bound-for-one-position)
· [Rank-one small-feature bound](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#a-small-feature-bound-on-a-single-fibre)
· [Supported-value enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)
· [Private numerical reserves on one box (NB3)](ABC_SHARED_BOX_PRIVATE_VALUATION_BANDS.md#uniform-application-to-all-edges-of-a-box)

## Exponents after primitive normalization

For each ordered residue list, prescribe a dictionary
of at most \(D\) coefficient vectors
\(\alpha=(\alpha_1,\alpha_2,\alpha_3)\) of nonzero integers
with \(\max_i|\alpha_i|\le X^{H_c}\). These dictionaries
may depend on the height, the original family, and all
residue lists. Each dictionary is fixed before choosing
the actual tree labels of its tuple. Rational coefficients
of componentwise reduced height at most \(X^{H_c}\) can
also be used, by multiplying each vector by the product of
its three denominators and replacing \(H_c\) by \(3H_c\).
Fix a numerical output-height exponent \(H>0\).

A candidate output is obtained from the three nonzero
rational numbers \(\alpha_i\prod_e L_e^{a_{i,e}}\) with
\(a_{i,e}\in\mathbb Z\). By primitive normalization we
mean multiplication by a common rational scalar to obtain
a signed integral triple with gcd one, determined up to
overall sign. For each edge subtract \(\min_i a_{i,e}\)
from its three exponents. This common Laurent factor can
be canceled before evaluation. The resulting nonnegative
exponents will again be denoted \(a_{i,e}\).
Let
\[
       M_i=\alpha_i\prod_e L_e^{a_{i,e}},\qquad
       g=\gcd(|M_1|,|M_2|,|M_3|),\qquad T_i=M_i/g.
\]
Consider candidates satisfying \(\sum_i T_i=0\) and
\(\max_i|T_i|\le X^H\).

**LP2.** If \(W_{e,P}>X^\theta\) holds on all positions,
then for each \(i\)
\[
  \prod_e W_{e,P}^{a_{i,e}}\ \bigm|\
            |T_i|\!\prod_{k\ne i}|\alpha_k|,
  \qquad
  \sum_e a_{i,e}\le \frac{H+2H_c}{\theta}.                  \tag{LP2}
\]
Also \(a_{i,e}\le (H+H_c)/\theta\) individually.
For an edge with \(a_{i,e}>0\), choose \(k\ne i\) with
\(a_{k,e}=0\). For \(p\mid W_{e,P}\), none of the other
labels is divisible by \(p\). Hence \(v_p(g)\le v_p(M_k)
=v_p(\alpha_k)\) and \(v_p(|\alpha_kT_i|)\ge a_{i,e}v_p(L_e)\).
This gives the individual bound. The \(W_{e,P}\) are
pairwise coprime, so the product of the required prime
powers divides \(|T_i|\) times the product of the two
other coefficients, even when the choice of \(k\) differs
between edges. The same divisibility works with any fixed
number of output terms. It allows the coefficients to
contain arbitrary private primes and uses the height
after taking the actual numerical gcd.

## Exclusion of nonidentity specializations

**LP3.** Under the hypotheses of (LP1), with \(H,H_c\)
fixed, outside
\[
    O_{m,j,B,C,r,H,H_c,L}
        \bigl(R^m K^{u-1}(n^r+D)\bigr)                       \tag{LP3}
\]
ordered synchronized tuples, every candidate is the
evaluation of a formal polynomial identity in the
independent potentials of the \(j\) boxes, after common
Laurent factors are removed. Its primitive normalized
values are either the signed coefficient vector divided
by its numerical gcd, or an original triangle of one of the
boxes up to ordering and overall sign.

First omit the exceptional tuples in (LP1). By (LP2) there
is a fixed bound \(A=\lceil(H+H_c)/\theta\rceil\) on all
reduced exponents, independent of \(L\) as \(L\to0\).
For each fixed ordered residue list, assign arbitrary fixed
spanning trees in the \(j\) boxes. Every tuple lies in a
product of \(u\) coordinate sets of size at most \(K\),
one for each tree label. The tree coordinates parameterize
all potential differences by invertible linear coordinate
changes after translation is removed on each box.

For a coefficient choice and an exponent array the
candidate requires the polynomial
\[
                       \sum_{i=1}^3
                        \alpha_i\prod_e L_e^{a_{i,e}}
\]
to vanish. It has total degree at most \(qA\) in the tree
coordinates. A nonzero polynomial of degree \(d\) has at
most \(dK^{u-1}\) zeros on a grid of this kind. One can
enlarge each coordinate set to size \(K\) in \(\mathbb Q\)
and induct on \(u\). For a polynomial of degree \(b\)
in the last variable its leading coefficient has total
degree at most \(d-b\). It vanishes on at most
\((d-b)K^{u-2}\) prefixes, and each other prefix gives at
most \(b\) zeros. The one-variable statement is the root
bound. There are at most \(D(A+1)^{3q}\) polynomial
choices per ordered residue list. This bounds the
exceptions from all nonidentities by a constant times
\(D R^m K^{u-1}\).

The surviving polynomial identities are classified by
(GP1), applied to the disjoint union of \(j\) copies of
\(K_{m,m,m}\). Removal of a translation in each component
does not change polynomial identities in its differences.
Every triangle is within one component and there is no
four-clique. Thus a nonconstant identity after removing
the common monomial is one triangle with a common scalar
coefficient and the signs dictated by the orientations.
The underlying triangle is one original primitive row.
Common scalar and monomial factors disappear under
primitive normalization. In the constant case all
reduced exponents vanish and the coefficient vector sums
to zero.
[Triangle and four-clique identities (GP1)](ABC_SHARED_BOX_POLYNOMIAL_TRANSFERS.md#three-term-monomial-identities-on-a-difference-graph)
· [Aprodu–Farkas–Raicu–Suciu, graphic arrangement multinets, §9.8](https://arxiv.org/pdf/2512.10160v1)

### Near-full strata and the dictionary threshold

Use a simultaneous packet stratum of (LB4) with
\(d=s_1+s_2+s_3\), \(t=1/6-d\), and fixed positive
margins for \(t\) and the \(s_i\). Assume \(t\ge t_*>0\).
Let \(\mathcal B\) be its ordered box family and
\(B_\ell\) its box count on residue list \(\ell\).
Hölder gives
\[
       \sum_\ell B_\ell^j\ \ge\
              |\mathcal B|^j/(R^m)^{j-1}.
\]
Together with (LB4) this is at least \(R^mK^u\) up to
arbitrarily small normalized power losses at fixed \(m,j\).
Here \(K=n^{t+O(\eta)}\) in the notation of the packet
stratum. Use \(r=t_*/2\) in (LP3). For every fixed
\(\sigma>0\), if \(D\le n^{t-\sigma}\), the exceptional
tuples form a power-small fraction, with sufficiently small
fixed tolerances, sufficiently small fixed \(L\), and then
sufficiently large \(X\). This applies uniformly across
the surviving balanced profiles with the indicated price
margins. It excludes a full-exponent subclass of these
tuples producing a new primitive additive output by
Laurent monomials in the original labels with this
coefficient cost. No bound on the radical or on the
unreduced height of the candidate output is assumed.
[Forced box exponent (LB4)](ABC_SHARED_VALUE_RESIDUE_BOXES.md#consequences-for-simultaneous-packet-strata)

The order \(K\) dictionary threshold is sharp for arbitrary
dictionaries and the stated output-height test. Fix one
edge position and, for every possible value \(a\) in its
pair class, include the coefficient vector \((1,1,-a-1)\).
There are at most \(K\) records per list. Use the monomials
\((1,L_e,1)\). On every tuple with \(L_e=a\), the resulting polynomial
is nonzero but gives
the primitive output \((1,a,-a-1)\). On the tuples
with all \(W_{f,P}>X^\theta\), no original label is \(1\),
so this output is new. The coefficient and primitive
output heights are at most \(X+1\).

## Growing synchronized tuples from one dense residue list

Fix \(m,C,B,H,H_c,\tau,\sigma>0\), with integer \(m\ge2\).
Use a fixed ordered residue list on which the number of
actual boxes is \(M\ge K^{h-1}n^{-\eta}\), and assume
\(K\ge n^\tau\). Fix the prime pool \(P\) and a coefficient
dictionary for this list as above. Put \(r=\tau/2\).
Choose \(\theta,L_0,\delta>0\) from the rank-one
small-feature bound used in (LP1) at exponent \(r\), and
put \(A_*=\lceil(H+2H_c)/\theta\rceil\).
Take
\[
          0\le\eta<\min\{\tau/8,\ \sigma/(12A_*)\},
               \qquad D\le K n^{-\sigma}.
\]
These inequalities are at the original scale \(n=X^L\).

**LP4.** At each sufficiently small fixed \(L>0\), and with
\(1\le j\le(\log X)^{\delta L/8}\), an
\(O(n^{-\min(\tau/4,\sigma/2)})\) fraction of the ordered
\(j\)-tuples on this list can fail the joint private
numerical height test or the conclusion of (LP3).
The constants and height threshold are uniform in \(j\)
in the stated range, the residue list, the coefficient
dictionary and the original family.

To test one edge on one member of the tuple, choose the
other \(j-1\) boxes from the \(M^{j-1}\) actual choices.
Omit a potential in the tested box as in (LP1) and record
the \(h-2\) anchor labels from their pair classes.
The new anchor pool is supported on the product of
\(P\), all labels of the other boxes, the anchor labels
between any two anchor vertices, and the nonzero
differences of neighbors of the omitted vertex.
The prime product is at most \(X^{O_{m,B}(j+1)}\).
It has \(k=O_{m,B}((j+1)\log X)\) distinct primes.

The supported-part count is still subpower with the
needed explicit allowance. Put \(s=\delta L/4<1\),
decreasing \(L_0\) if needed. The number of integers at
most \(X\) supported on the pool is at most
\[
 X^s\prod_{p\ {\rm in\ pool}}(1-p^{-s})^{-1}
       \le X^s\exp\!\bigl(O_s(k^{1-s})\bigr).
\]
The last bound follows by ordering the primes, lower
bounding the \(i\)-th of them by \(i+1\), and summing
\((i+1)^{-s}\). The assumed \(j\) gives
\(k^{1-s}=O_{m,B}((\log X)^{(1+s/2)(1-s)})=o(\log X)\)
at fixed \(L\). There are therefore at most \(X^{2s}\)
choices for the supported part at all sufficiently large
heights, uniformly in the pool. The rank-one
small-feature count is unchanged. Its \(n^{r-\delta}\)
bound after a supported part is fixed can absorb the
factor \(X^{2s}=n^{\delta/2}\).
Summing over the tested positions, the bad fraction is
\(O_m(j K^{h-2}n^r/M)\), which is
\(O_m(j n^{r+\eta}/K)\).

On the remaining tuples the total exponent bound (LP2)
gives at most \(A_*\) factors counted with multiplicity
in each reduced term. Any candidate thus uses at most
\(3A_*\) boxes, and there are \(O_{m,A_*}((1+j)^{3A_*})\)
reduced exponent arrays over all positions. For one
candidate with \(b\) active boxes, where \(1\le b\le3A_*\),
apply the polynomial grid bound to the \(b(h-1)\)
tree coordinates of those boxes. A nonzero polynomial
has degree at most \(A_*\). The other boxes have
\(M^{j-b}\) choices. Its zero proportion is at most
\[
             \frac{A_*}{K}\left(\frac{K^{h-1}}{M}\right)^b
               \le \frac{A_*n^{3A_*\eta}}{K}.
\]
When no box is active the polynomial is constant and
a nonidentity has no zeros. The dictionary may also depend
on the active set and on the actual boxes on its complement.
For a prescribed exponent array and each fixed complement
the same grid bound applies to at most \(D\) coefficient
vectors chosen before the active boxes. Sum over the coefficients
and exponent arrays. The exceptional fraction is
\[
 O\left(\frac{j n^{r+\eta}}{K}
       +(1+j)^{3A_*}\frac{D\,n^{3A_*\eta}}{K}\right),
\]
which gives (LP4) at each fixed \(L\). This range in
particular includes any choice \(j=(\log X)^{o(1)}\)
at that fixed budget.

Suppose more generally that the total number of boxes
is at least \(R^m K^{h-1}n^{-\eta/2}\) for a fixed
\(\eta>0\) as above, and apply the same dictionary cap
on all lists. Hölder gives at least
\(R^m K^{j(h-1)}n^{-\eta j/2}\) synchronized tuples
over all lists. Lists with \(M<K^{h-1}n^{-\eta}\)
contribute at most \(R^mK^{j(h-1)}n^{-\eta j}\).
Their relative contribution is at most \(n^{-\eta j/2}\).
Applying (LP4) on the other lists gives the global
exceptional fraction
\(O(n^{-\min(\tau/4,\sigma/2)}+n^{-\eta j/2})\).
This total box count hypothesis follows from (LB4) for
every fixed positive \(\eta\) on a surviving near-full
stratum with sufficiently small fixed tolerances.
The count permits repetitions when drawing tuples.
Repeated boxes fail the joint private part test and
thus are absent from every retained tuple.

The graphic identity classification is already (GP1),
with the corresponding triangle and four-clique patterns
in graphic arrangement theory. The additional input here
is private numerical height after primitive reduction,
which removes the degree restriction even for mixed
Laurent monomials of several boxes. The result controls
one proposed source of an economical arithmetic transfer.
It leaves open whether the original caps force a transfer
using additional factors or a direct loss in the original
box count.
