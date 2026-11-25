# Fixed shifts and energy without primitivity

The arbitrary-shift count extends to signed pairs with common divisors.
The common divisor is supported on the prescribed shift. Stripping
that support can leave quotients of vanishing numerical height.
In that case the low-degree root dictionaries bound the count.
The resulting estimates are uniform in the prescribed shift and give
an energy saving for two complete-value dictionaries with diffuse
degree profiles, even when their values overlap in prime support.

## Two dictionary bounds for any nonzero shift

Let \(X\ge 2\), fix \(K,C_0,a>0\), and use signed equations
\[
                \sigma_U U+\sigma_V V=t\ne0,\qquad 0<U,V\le X^K.
\]
Here \(\sigma_U,\sigma_V\in\{-1,1\}\), and \(t\) is fixed in each
count. The complete values are each a fixed positive integer
coefficient times a product of positive integral roots. The
positive-index coordinates of each term are pairwise coprime.
The two terms can have an arbitrary common divisor. Fix the finite
position and degree lists before taking the height limit at each
fixed small budget \(\lambda\). The sums of the raw root height
exponents on the two sides have total at most \(C_0\lambda\).

Give position \(i\) an individual dictionary of at most
\(X^{W_i+o_X(1)}\) elements, with \(W_i\) at most its raw cap. Write
\(W_U,W_V\) for the respective sums of dictionary costs and
\(B_U,B_V\) for their sums over degrees less than \(a/\lambda\).
For every \(\eta,\varepsilon>0\) and sufficiently small fixed
\(\lambda\), the count of distinct complete pairs is
\[
 \ll X^{Q+\eta\lambda+\varepsilon},
 \quad
 Q=\min\left\{
 \begin{aligned}
  &W_U,\ W_V,\\
  &\max\{(W_U+W_V)/3,\min(B_U,B_V)\},\\
  &\max\{\min(W_U,W_V)/2,B_U+B_V\}
 \end{aligned}\right\}.                                           \tag{NS1}
\]
The budget threshold may depend on \(K,C_0,a,\eta\).
At each budget the constants can also depend on the finite degree
lists and \(\varepsilon\). They are uniform in the shift and both
coefficients under the stated height cap. Use \(\log n=\lambda\log X\)
for normalized exponents.

Each complete value determines its partner, giving the first two
trivial caps in (NS1). To prove the others, put \(S=\operatorname{supp}(t)\),
and fix the entire \(S\)-supported part of every original root.
At fixed truncation that takes a subpower number of choices, uniformly
in \(t\), since \(|t|\le2X^K\). The fixed-support bound used here is
the Rankin count at arbitrarily small fixed parameter: an integer
supported on a set of at most \(O_K(\log X/\log\log X)\) primes has
\(X^{o_X(1)}\) choices at polynomial height. For instance, the Euler
product for Rankin's trick at any fixed positive parameter has logarithm
\(o(\log X)\) by comparing the primes with the first integers of the
same cardinality.

Choose also \(g=\gcd(U,V)\mid t\), at divisor cost. The two quotients
\(\widetilde U=U/g,\widetilde V=V/g\) are coprime, and their equation
has shift \(t/g\ne0\). All primes of \(g\) lie in \(S\), so the
quotients are products of the residual original roots and fixed
integer coefficients. The raw caps and the dictionary costs do
not increase, since the operation on every dictionary is just
removing its supported part. Rebox the residual roots dyadically
at subpower cost.

The final cap in (NS1) now follows from the primitive arbitrary-shift
estimate (FS5) applied to singleton groups. It allows arbitrary
fixed coefficients and every numerical height at most \(X^K\).
The low-degree bound on the stripped dictionaries is at most the
original \(B_U+B_V\).
[Arbitrary shifts with high free degrees (FS5)](ABC_FULL_SHIFT_DICTIONARIES_AND_CODEGREES.md#arbitrary-shifts-with-high-free-degrees)

### Positive-height quotients with arbitrary coefficients

For the middle cap in (NS1), use successive limits. At a candidate
contradiction take arbitrarily large heights at each fixed budget,
and then a sequence of budgets tending to zero. Pass to subsequences
on which the normalized costs, the count exponent, and the numerical
heights converge. Denote normalized costs by \(w_U,w_V,b_U,b_V\), and
the normalized primitive pair count exponent by \(d\). In this
argument all costs are those of the stripped dictionaries, bounded
by the original costs. If the maximum of the numerical quotient
heights tends to zero, the high-degree coordinates on each side
have total raw cap \(o(\lambda)\): their degrees are at least
\(a/\lambda\), and the total numerical height of all their powers
is \(o(1)\). On either side the remaining dictionary costs are
therefore at most \(B_U+o(\lambda)\) or \(B_V+o(\lambda)\),
respectively. An exact value fixes the other term. This gives
\(d\le\min(b_U,b_V)\).

Suppose instead that at least one quotient has a positive limiting
numerical height. To prove \(3d\le w_U+w_V\), assume the reverse
strict inequality. Trivial enumeration gives \(d\le\min(w_U,w_V)\),
so \(w_U,w_V<2d\). The one-residue concentration and its symmetric
application (FT3)–(FT4) apply directly to the stripped primitive pairs.
They allow arbitrary fixed coefficients. Each variable part has just
boundedly many original roots of positive limiting numerical weight,
and each such root costs at least \(d\). Since the whole cost on
either side is below \(2d\), there is at most one positive-weight
root per side.
[Fixed-shift concentration (FT3)–(FT4)](ABC_FIXED_SHIFT_TWO_ROOT_BOUND.md#one-fixed-residue-forces-concentration)

If there is one positive-weight root on each side, fix both
complements at normalized cost at most \(w_U+w_V-2d\). The
two remaining positive roots of height \(X^{O(\lambda)}\)
satisfy \(A x^p+B y^q=t/g\) with signed nonzero fixed coefficients,
nonzero shift, and \(p,q\gg 1/\lambda\). This affine plane curve is
absolutely irreducible by the simple zero of \(t/g-Ax^p\), and
its degree is at least a fixed positive multiple of \(1/\lambda\).
The coefficient-uniform plane-curve bound contributes
\(O(X^{O(\lambda^2)+\varepsilon})\) solutions. This gives the
contradiction \(3d\le w_U+w_V\).
[Two-root release (FT1)](ABC_FIXED_SHIFT_TWO_ROOT_BOUND.md#releasing-two-original-roots)

Suppose one whole variable part, say the one on \(\widetilde U\),
has vanishing numerical weight. If its fixed coefficient
\(A'_0\) has positive weight and \(\widetilde V\) has a positive
variable root \(z\), fix its complement at cost \(w_V-d<d\).
The remaining equation is a unit congruence
\[
                         z^e \equiv u\pmod{A'_0}
\]
for a fixed unit \(u\). The coefficient of \(z^e\) is a unit
modulo \(A'_0\) by primitivity. At a fixed degree \(e\), the
kernel of the \(e\)-power map on the units of any prime power
has size at most \(2e\): for odd primes the unit group is cyclic,
and at \(2\) it is a product of a cyclic group and a cyclic group
of order at most two. The modulus therefore has at most
\((2e)^{\omega(A'_0)}=X^{o_X(1)}\) solution classes. Since
\(z\le X^{O(\lambda)}<A'_0\), each class specifies at most
one \(z\). This contradicts the assumed count. If the variable
parts on both sides have vanishing weight, a positive-weight
coefficient \(A'_0\) instead determines the whole opposite variable
part directly modulo \(A'_0\). The opposite fixed coefficient
is a unit and the variable part has height \(X^{o(1)}<A'_0\).
At least one coefficient has positive weight in this case.

In the remaining case \(\widetilde U\) also has vanishing complete
height, so \(\widetilde V\) is tall and has ratio to the fixed
\(|t/g|\) within \(O(X^{-\tau})\) of one for some fixed
\(\tau>0\). Use (RUW1) on the narrow roots and a wider cofactor
of the tall term, with its arbitrary fixed coefficient absorbed
into the centre of the ratio. Each exact value of the tall term
determines both complete values, with subpower representation
multiplicity. Thus the count exponent is \(o(\lambda)\).
This proves the middle cap in (NS1). The use of a positive limiting
height in this step does not impose a new uniform margin: any
violating sequence whose height limit is zero has already been
counted by the low-degree limb.
[Unequal-height real strips (RUW1)](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-real-strip-with-unequal-root-heights)

## Packet size and varying-shift codegrees

In the setting of (NS1), suppose each pair in a selected family
has, on side \(U\), a divisor record \((g_i)\) in at most \(M\)
original root positions, all of degree at least \(a/\lambda\).
The choices can vary with the pair. Require \(g_i\mid x_i\) and
\[
                           \prod_i g_i\ge X^{\beta\lambda}
\]
for a fixed \(\beta>0\). Let \(\ell_U\) be the normalized sum of
the raw root caps on \(U\), and \(b_U,b_V\) the normalized
low-degree dictionary costs as above. The normalized pair
count exponent \(d\) satisfies, in the same limiting convention,
\[
                    d\le\max\{b_U+b_V,\ell_U-\beta\}.               \tag{NS3}
\]
The bound is uniform in the shift and both coefficients at
polynomial height. It applies to the entire selected family
of pairs satisfying the packet condition.

Fix the positions of the packet by pigeonholing. Their lists
are finite before the height limit. Choose one packet per
actual complete pair. Use the support and gcd refinements
in the proof of (NS1), where \(S=\operatorname{supp}(t)\).
If \(h_i\) is the now fixed \(S\)-supported part of \(x_i\),
write \(x'_i=x_i/h_i\) and let \(g'_i\) be the part of \(g_i\)
supported outside \(S\). Then \(g'_i\mid x'_i\), and
\(h_i g'_i\ge g_i\). For a given transformed packet vector
\((g'_i)\), the number of residual root tuples in these
positions is at most
\[
                X^{\sum_{i\text{ in packet}}\alpha_i
                            -\beta\lambda+o_X(1)} ,
\]
where \(\alpha_i\) are their original raw caps. To see this,
count multiples of \(g'_i\) up to
\(X^{\alpha_i+o_X(1)}/h_i\). This bound depends only on the
fixed supported parts and the chosen transformed packet.
Their product has the required lower bound because the vector
comes from an actual pair.

There are at most \(X^{d\lambda+o_X(1)}\) transformed vectors
on a family with normalized exponent \(d\), since the original
packets were chosen one per pair. Thus the selected residual
positions form one group of size at most \(M\) with cost at most
their raw cap minus \(\beta\lambda\) plus \(d\lambda\), up to the
arbitrarily small budget allowance. Count every other position
on that term by its raw cap, and use at most the original
low-degree dictionary costs for the low positions. The normalized
total dictionary cost on \(U\) is no more than
\(\ell_U-\beta+d+o(1)\). Apply (FS5), using this one bounded
high-degree group and the remaining singleton groups. If
\(d>b_U+b_V\), the bound gives \(2d\le\ell_U-\beta+d\),
which proves (NS3).
[Arbitrary shifts with high free degrees (FS5)](ABC_FULL_SHIFT_DICTIONARIES_AND_CODEGREES.md#arbitrary-shifts-with-high-free-degrees)

For example, take two fixed distinct columns of a varying-shift
graph, with original targets \(T_i=Fxy_i+\sigma Gk_i\). Their
common rows have the nonzero fixed weighted difference
\[
                         y_jT_i-y_iT_j=D_{ij}.
\]
Count the pairs of weighted terms in (NS3), absorbing the fixed
\(y_i,y_j\) into their respective coefficients. Their original
variable positions are the target root positions. If the
packet condition is imposed on \(T_i\), its original mass
lower bound is still valid in the dictionary estimate after
the supported parts and gcd are fixed. A prime lying in the
support of \(D_{ij}\) is charged by the reduction of the range
by \(h_i\), even if its full power disappears in primitive
normalization. Recovering the weighted terms recovers a
common row. Thus (NS3) gives a codegree bound with no radical
hypothesis on the difference. In particular when the target
raw cap tends to \(1/3\), the combined low-degree cost is
strictly below \(1/6\), and a packet has mass bounded away
above \(1/6\), the common-row exponent is strictly below \(1/6\).
[Weighted differences in a varying-shift graph](ABC_FULL_SHIFT_DICTIONARIES_AND_CODEGREES.md#weighted-differences-in-a-varying-shift-graph)

## Energy of complete-value dictionaries

Let \(P,Q\) be two complete-value dictionaries with the stated
coordinate representations, with \(n^{s_P+o(1)}\) and
\(n^{s_Q+o(1)}\) distinct values, respectively. Use normalized raw
or charged masses \(w_P,w_Q\). Suppose for a fixed cutoff
\(a>0\) the sum of their low-degree charged masses satisfies
\(b_P+b_Q\le\min(w_P,w_Q)/2\). Then every nonzero sum or difference
fibre has size
\[
                           n^{\min(w_P,w_Q)/2+o(1)}.
\]
The error uses the height-first budget limit. Hence
\[
 \begin{split}
 E(P,Q)&=\sum_t r_{P+Q}(t)^2
       \le n^{s_P+s_Q+\min(w_P,w_Q)/2+o(1)},\\
 |P+Q| &\ge n^{s_P+s_Q-\min(w_P,w_Q)/2-o(1)}.
 \end{split}                                                     \tag{NS2}
\]
Indeed \(\sum_t r_{P+Q}(t)=|P||Q|\) and the sums are positive.
The mixed difference energy equals the sum energy, by exchanging
one pair of variables. Its zero fibre has size \(|P\cap Q|\), whose
square is at most \(|P||Q|\); it does not affect the exponent in
(NS2). If both dictionaries have positive exponents, the same
sumset bound applies to differences by discarding the zero fibre
before dividing the pair count by the largest nonzero multiplicity.

For example if \(S\) has raw mass \(\ell>0\), cardinality
\(n^{\ell+o(1)}\), and no limiting atom at scaled degree zero, choose
a fixed positive cutoff with raw mass below \(\ell/4\) on each
copy. The bound applies even if some first moment escapes to
infinite scaled degrees. It gives
\[
                       E(S)\le n^{5\ell/2+o(1)},\qquad
                       |S\pm S|\ge n^{3\ell/2-o(1)}.
\]
For a wholly diffuse dictionary these are fixed exponents
instead of a profile-dependent positive saving. The proof
needs no coprimality between the two complete values in a
pair. The all-three-letter incidence at the balanced endpoint
can still have complete-value dictionaries of size
\(n^{1/6+o(1)}\) with root masses \(1/3+o(1)\). In that range the
shift exponent supplied here is \(1/6\); an additional power
saving in the simultaneous incidence count remains necessary.
[Earlier wholly diffuse shift and energy estimates (DS1)](ABC_DIFFUSE_FACTOR_GRAPH.md#additive-energy-of-wholly-diffuse-complete-value-dictionaries)
