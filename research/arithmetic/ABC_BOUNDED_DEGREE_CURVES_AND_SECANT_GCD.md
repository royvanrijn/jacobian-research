# Bounded-degree curves and secant gcds in balanced abc families

A moving plane curve of fixed degree contains a power-smaller part of
a balanced small-budget abc shell. One can use this to control the
common scaling gcd of affine secants. For almost every query, a bounded
number of its secants supplies a prime pool containing all but a short
factor of the gcd on almost every other secant from that query.

## Residues from a polynomial on one locked letter

Fix \(D\ge1\). If \(G\in\mathbb Z[Y]\) has degree at most \(D\)
and has a coefficient prime to \(p\), its roots modulo \(p^e\)
occupy at most \(D\) classes modulo \(p^{\lfloor e/D\rfloor}\)
when \(e\ge D\). Indeed, otherwise choose \(D+1\) integer
representatives \(v_i\) of roots in distinct classes at that
precision. In the Lagrange interpolation formula
\[
    G(Y)=\sum_i G(v_i)
               \frac{\prod_{j\ne i}(Y-v_j)}
                    {\prod_{j\ne i}(v_i-v_j)},
\]
each denominator has \(p\)-adic valuation at most
\(D(\lfloor e/D\rfloor-1)<e\). As each \(G(v_i)\) is divisible
by \(p^e\), every term has all coefficients in \(p\mathbb Z_p\).
This contradicts the unit coefficient of \(G\).

Let \(M\le X\) be a divisor of a complete letter \(W\), prime to one
nonzero coefficient of \(G\), and assume all \(e=v_p(M)\) are at
least \(2D\). Define
\[
          M_0=\prod_{p^e\parallel M}p^{\lfloor e/D\rfloor}.
\]
Then \(M_0\ge M^{1/(2D)}\). By the Chinese remainder theorem,
\(G(V)\equiv0\bmod M\) constrains \(V\) to at most
\(D^{\omega(M)}=X^{o_X(1)}\) classes modulo \(M_0\). If \(V\) is a
different complete letter of a primitive triple, each relevant class
is a unit. The number of classes is uniform in \(G\) and \(M\).

### One complete letter in a unit residue

For fixed \(A,\omega,\xi>0\), sufficiently small fixed \(L>0\), and
\(0\le a\le A\), uniformly in moduli \(X^\omega\le m\le X\) and
unit classes \(r\bmod m\), one has
\[
 \#\{1\le V\le X:\operatorname{rad}(V)\le X^{aL},\
                               V\equiv r\bmod m\}
       \ll_{A,\omega,\xi,L,\varepsilon}
                       X^{(a/2+\xi)L+\varepsilon}.
 \tag{BC1}
\]
Here is the single-letter form of the multiplicative moment underlying
(EB1). Choose an exact valuation cutoff \(J\) at each fixed budget.
The full product of the powers above order \(J\) has radical at most
\(X^{1/J}\), and fixing it costs \(X^{1/J+o_X(1)}\), uniformly.
Write the remaining part as \(\prod_{j\le J}x_j^j\), with \(x_j\)
the squarefree exact valuation roots. Partition the root sizes into
dyadic ranges and enlarge to the Cartesian product of their allowed
unit values. The dyadic upper exponents \(\alpha_j\) have sum
\(H\le aL+o_X(1)\), since the original exact roots are disjoint.

For one root set \(\mathcal D_j\), put
\(S_j(\chi)=\sum_{x\in\mathcal D_j}\chi(x)^j\).
Round the \(J\) upper exponents to a positive mesh small enough
that the total rounding error is at most \(\xi L/2\). At fixed
\(L,J\) this gives finitely many possible cap lists and
\(H\le(a+\xi/2)L+o_X(1)\). With the rounded caps set
\(h_j=\lceil H/\alpha_j\rceil\), so
\(h_j\alpha_j\le2H\). Products in the \(2h_j\)-th moments have
height \(X^{2H+o_X(1)}<m\) if \(2(A+\xi)L<\omega\).
All \(h_j\), the cutoff, and the powers \(j\) belong to fixed
finite lists at the fixed budget.

Given the first \(h_j\) factors of a moment pair, the product of the
other \(h_j\) factors lies in one coset of the kernel of the \(j\)-th
power map on \((\mathbb Z/m\mathbb Z)^\times\). This kernel has size
at most \(j^{2\omega(m)}\), by the unit groups of prime powers, and
hence \(X^{o_X(1)}\). Each such residue has at most one integer
product in the allowed height range. Factoring that product costs
\(X^{o_X(1)}\) at fixed \(h_j\). Consequently, with normalized
character averages,
\[
 \mathbb E_\chi |S_j|^{2h_j}
                \le X^{h_j\alpha_j+o_X(1)}.
\]
Apply Hölder with \(\sum_j h_j^{-1}\le1\), using the constant
function for the remaining reciprocal exponent, to obtain
\(\mathbb E_\chi|\prod_j S_j|^2\le X^{H+o_X(1)}\).
Orthogonality and Cauchy--Schwarz now bound the count in a fixed unit
class by \(X^{H/2+o_X(1)}\). The fixed tail contributes a unit twist.
Take \(1/J\le\xi L/4\) and absorb the inner losses in
\(\varepsilon\). The finite cap and moment lists make the bound
uniform in the dyadic exponents and moving moduli.
[Multiplicative residue moments and their numerical height condition](ABC_EFFECTIVE_RADICAL_BUDGET.md#a-diagonal-moment-with-two-scales)
· [Uniform fixed-support enumeration](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

## Uniformly affordable full divisors

Fix absolute \(p_0\in(0,1/48)\) and \(\gamma=1/100\). There
are absolute \(0<s<S<\infty\) and \(\kappa_0,\delta_0,L_0>0\)
with the following property. In the shell \(X/2<c\le X\) with
\(\operatorname{rad}(abc)\le c^L\) and
\(\operatorname{rad}(W)\le n^{1/3+\kappa_0}\) on each letter,
all but \(O_{L,\varepsilon}(X^{(1/3-\delta_0)L+\varepsilon})\)
rows have some role \(W\) for which
\[
  \sum_{\substack{p^v\parallel W\\
                   s\le Lv\le S,\ \log_n p\le p_0/2}}
                 \log_n p\ \ge 3\gamma .
 \tag{BC0}
\]
This holds at every fixed \(0<L<L_0\).

Otherwise take heights arbitrarily large at fixed budgets tending
to zero, and take \(\kappa_0,\delta_0\) tending to zero and tested
bands \([s,S]\) increasing to \((0,\infty)\). Refine the violating
rows into exact valuation boxes before taking the budget limit.
The unrestricted limiting count gives exponent \(1/3\) in base
\(n\). The cap on each role and the exclusion of a fixed total
radical deficit give limiting mass \(1/3\) on each role. The sum
of the limiting nonatomic masses is at least \(1/2\), so some role
has nonatomic mass at least \(1/6\).

Choose a compact positive-degree band on this role with nonatomic
mass above \(1/12\). Within that band remove small neighborhoods
of the finitely many atoms of mass at least \(p_0/4\), losing less
than \(1/48\) of nonatomic mass, and take the band and deleted
neighborhood boundaries to be continuity points. The compact
remainder retains mass greater than \(1/16\). Every exact-degree
root with scaled degree there has normalized mass at most \(p_0/2\)
in sufficiently late boxes. A sequence of roots violating this
would give an atom of mass at least \(p_0/2\) at a limit point of
the compact remainder. The full normalized mass of each such
root is the sum of its prime masses, so every prime in it meets
the prime-size bound in (BC0). Weak convergence on the compact
remainder and \(1/16>3\gamma\) give the contradiction. At each
fixed budget the layers in any bounded scaled-degree band are
finite. Their root heights can be binned dyadically, making
their masses uniform in each refined box to \(o_X(1)\). The
radical mass above scaled degree \(R\) is at most \(1/R\), so
the finite refinements and eventual compact-band tests also
apply to exact valuation measures.
[Limiting small-budget count](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)
· [Nonatomic mass inequalities](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)
· [Exact valuation measures and the height-first convention](ABC_SECANT_ENDPOINT_TRANSFER.md#the-height-first-setting)

## A uniform saving on moving curves

There are absolute \(\delta,\kappa>0\) such that for every integer
\(D\ge1\) there is \(L_D>0\) with the following property. Fix
\(0<L<L_D\), put \(n=X^L\), and take positive primitive \(a+b=c\)
with \(X/2<c\le X\), \(\operatorname{rad}(abc)\le c^L\), and
\(\operatorname{rad}(W)\le n^{1/3+\kappa}\) for \(W=a,b,c\).
For every nonzero real polynomial \(F(a,b)\) of degree at most \(D\),
the number of these rows on \(F(a,b)=0\) is
\[
                  O_{D,L,\varepsilon}
                         (X^{(1/3-\delta)L+\varepsilon}).
 \tag{BC2}
\]
The implied constant is uniform in \(F\), which may move with \(X\).
The exponent saving and the cap tolerance are independent of the
fixed degree. The small-budget threshold may depend on \(D\).

On the integer points of the curve one may replace \(F\) by a
nonzero integer polynomial of degree at most \(D\) and coefficient
height \(X^{O_D(1)}\). The degree-\(D\) monomial evaluation matrix has
integer entries of height at most \(X^D\) and a nontrivial kernel.
Choose a nonsingular \(r\times r\) minor at its rank \(r\), and any
nonpivot column. Cramer's rule expresses the latter in terms of the
pivot columns, with integral coefficients after multiplication
by the minor. This gives the required integer polynomial and
coefficient bound. In rank zero use any standard basis vector.

Discard the exceptions to (BC0) and split the other rows into
three classes according to a qualifying role \(W\). On each class
fix any other role \(V\). The change from \((a,b)\) to \((W,V)\)
is integral and invertible over the integers. Express the integer
polynomial in \((W,V)\) and divide off the largest monomial power
of \(W\) dividing it. This preserves vanishing on the positive
sample and ensures that \(G(V)=F(0,V)\) is a nonzero polynomial.
Choose a nonzero coefficient \(C\) of \(G\), of height \(X^{O_D(1)}\).

Discard rows on which the normalized radical mass of \(W\) on the
primes of \(C\) is at least \(\gamma\). By (EB5) these rows number
at most \(O_{D,L,\varepsilon}(X^{(1/3-\gamma/4)L+\varepsilon})\)
for sufficiently small fixed \(L\). The prime pool of \(C\) has
\(O_D(\log X/\log\log X)\) elements, so it admits zero prices.
This discard is on the original triples and is uniform in \(C\).
[Prime-pool discount (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

At least \(2\gamma\) normalized mass from (BC0) remains outside
the primes of \(C\). Select its full prime powers greedily up to
a normalized mass between \(p_0/2\) and \(p_0\), using the bound
on the size of each prime. Their product \(M\) satisfies
\(M\ge X^{sp_0/2}\), \(\operatorname{rad}(M)\le n^{p_0}\) and
\(v_p(M)\ge s/L\) at every \(p\mid M\). There are
\(n^{p_0+o_X(1)}\) such integers \(M\), by radical and
fixed-support enumeration. It suffices to count each lock class.

On a locked class \(G(V)\equiv0\bmod M\). For \(L<s/(2D)\)
the first reduction gives \(X^{o_X(1)}\) unit classes modulo
\(M_0\ge X^{sp_0/(4D)}\). Apply (BC1) to \(V\) with
\(a=1/3+\kappa\) and fixed \(\xi=1/48\). Except for at most \(D\)
complete values of \(V\), the polynomial \(F(W,V)\) with \(V\)
fixed is nonzero in \(W\) and has at most \(D\) possible \(W\).
Indeed, exceptional \(V\) must be a root of any nonzero
coefficient polynomial in \(V\) of \(F\) as a polynomial in \(W\).
For each exceptional \(V\), (FL1) bounds its original fibre
by exponent \(1/6+\kappa/2+1/48\) in base \(n\) with inner
allowances, using both other radical caps. Summing over locks
gives normalized exponent at most
\(p_0+1/6+\kappa/2+1/48\) with inner allowances. For
\(\kappa\le\min(\kappa_0,1/48)\) this is below \(1/4\).
Together with the two preceding discards this proves (BC2)
with any sufficiently small absolute
\(\delta<\min(\delta_0,\gamma/4,1/12)\). The weakening of the
modulus by a factor depending on \(D\) enters only the choice of
\(L_D\) in (BC1).
[Fixed-complete-letter cap (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

## Local prime pools for common secant gcds

For two distinct rows \(x=(a,b)\), \(y=(a',b')\), put
\[
                g(x,y)=\gcd(|a-a'|,|b-b'|).
\]
The gcd is positive and at most \(X\), and also divides \(c-c'\).
Fix \(\eta>0\). There exist an integer \(k=k(\eta)\ge1\)
and positive \(\kappa,L_0,\sigma\), depending only on \(\eta\),
with the following property. At every fixed \(0<L<L_0\) take
any family \(\mathcal E\) in the shell of (BC2), with \(\kappa\)
no larger than its absolute cap tolerance. There exist integers
\(1\le P_x\le X^k\), one for each \(x\in\mathcal E\), such that
\[
 \#\left\{(x,y)\in\mathcal E^2:
 x\ne y,\quad
 \frac{g(x,y)}{\gcd(g(x,y),P_x)}>X^\eta\right\}
                \ll_{\eta,L} n^{2/3-\sigma}.
 \tag{BC3}
\]
Each nontrivial \(P_x\) can be taken as the least common multiple
of at most \(k\) gcds from \(x\) to other family members. In
particular, apart from the stated exceptional pairs, the full
factor of \(g(x,y)\) outside the primes of \(P_x\) is at most
\(X^\eta\). For each \(x\) the part supported on \(P_x\) has
subpower many possible values, uniformly, even when its radical
exceeds every fixed power of \(n\).

Choose \(k\) sufficiently large that \(k\eta>2D+1\), where \(D\)
is the smallest integer with \(\binom{D+2}{2}\ge2k+1\).
One can take \(k\ll 1+\eta^{-2}\) since \(D=O(\sqrt{k})\).
Fix any \(k\) anchor
rows \(y_1,\ldots,y_k\). There is a nonzero integer polynomial
of degree at most \(D\), vanishing on all anchors, with coefficients
\(O_k(X^D)\). To see this, use coefficients in
\(\{0,\ldots,\lceil X^D\rceil\}\). If \(s=\binom{D+2}{2}\), there
are more than \(X^{Ds}\) coefficient vectors for large \(X\).
Their integer evaluation vectors on the anchors number at most
\((1+s\lceil X^D\rceil X^D)^k=O_k(X^{2Dk})\). A difference of two
vectors with equal evaluations supplies the polynomial.

If a query \(x\) differs from every anchor, the least common
multiple of \(g(x,y_j)\) divides that polynomial's value at \(x\).
A least common multiple exceeding \(X^{2D+1}\) forces vanishing,
since the value has absolute size \(O_k(X^{2D})\). By (BC2), for
each ordered anchor tuple there are
\(O_{\eta,L}(n^{1/3-\delta/2})\) such queries. Take the inner
\(\varepsilon\) in (BC2) smaller than \(\delta L/2\).
Writing \(E=|\mathcal E|\), the total number of these ordered
stars is at most
\[
                        O_{\eta,L}(E^k n^{1/3-\delta/2}).
 \tag{BC4}
\]

Set \(\tau=\min\{\delta/(32k),1/24\}\) and put
\(T=\lceil n^{1/3-\tau}\rceil\).
Call \(x\) exceptional if, for every list of fewer than \(k\)
anchors, at least \(T\) possible next anchors have
\(g(x,y)/\gcd(g(x,y),P)>X^\eta\), where \(P\) is the least
common multiple from the list and the empty list has \(P=1\).
Such an \(x\) gives at least \(T^k\) ordered stars in (BC4):
adjoining each chosen gcd multiplies the least common multiple
by more than \(X^\eta\). Thus the number of exceptional queries
is \(O_{\eta,L}((E/T)^k n^{1/3-\delta/2})\).

For any other query choose a list of fewer than \(k\) anchors
with fewer than \(T\) failing neighbors and take \(P_x\) from
that list. Assign \(P_x=1\) on exceptional queries. The number
of failed ordered pairs is at most
\[
 O_{\eta,L}\bigl(ET+E(E/T)^k n^{1/3-\delta/2}\bigr).
\]
Shrink \(L_0\) so the limiting small-budget count gives
\(E\le n^{1/3+\tau/4}\) for all sufficiently large \(X\) at
the fixed budget, and keep \(\kappa\) within the cap tolerance of (BC2).
The first term has exponent at most \(2/3-3\tau/4\).
The second has exponent at most
\(2/3+(5k+1)\tau/4-\delta/2\), strictly below \(2/3\).
This proves (BC3) with \(\sigma=\tau/2\), in particular with a loss
bounded below by an absolute multiple of \(\min\{1,\eta^2\}\).
It remains a power-smaller pair class on families
of exponent \(1/3-o(1)\). The subpower count of the
\(P_x\)-supported factors uses fixed-support enumeration with
\(P_x\le X^k\), at the fixed budget.

### Peeling down to the normalized radical scale

After the Lagrange interpolation form of the residue lemma,
(BC2) holds with the same absolute gap whenever
\(D L\le c_*\), for a sufficiently small absolute \(c_*>0\)
and sufficiently small \(L\). Indeed the lock gives
\(M_0\ge X^{s p_0/(4D)}\) and each of its parent prime
powers has \(e\ge s/L\). The root-height requirements
in (BC1) therefore follow from \(D L\le c_*\).
The degree dependence of the bound for the coefficient
prime pool enters the inner height constants only. For
each \(D\) fixed before height, its cardinality is at
most \(O_D(\log X/\log\log X)\) and hence sublogarithmic,
so the prime-price count uses the same normalized
allowance and small-budget threshold.

For fixed \(L\) choose \(\eta=C_* L\) with a sufficiently
large absolute \(C_*\), and then \(k\asymp 1+\eta^{-2}\),
\(D=O(\sqrt{k})\) as in (BC3). This ensures \(D L\le c_*\).
To give a relative estimate, assume
\(E\ge n^{1/3-\delta/8}\), with \(\delta\) the absolute
curve gap, and now choose
\(T=\lceil E n^{-\delta/(8k)}\rceil\) in the greedy count.
Then (BC4) bounds the number of exceptional centres by
\(O_{\eta,L}(n^{1/3-3\delta/8})\), at most a proportion
\(O_{\eta,L}(n^{-\delta/4})\) of all centres.
At every other centre the exceptional fraction of neighbors
is at most \(T/E\). Thus one can take pools of height
\(X^{O(L^{-2})}\) for which, except on a fraction
\(O_L(n^{-c L^2})\) of ordered pairs for an absolute \(c>0\),
\[
                    g(x,y)/(g(x,y),P_x)\le n^{C_*}.
 \tag{BC5}
\]
All degrees and pool height exponents are fixed once \(L\)
is fixed. The numerical tolerance in (BC5) is smaller than
in (BC3). Its normalized power loss tends to zero with \(L\),
so (BC3) is the version used for a uniform power loss
at the small-budget endpoint.

## Use with primitive direction bounds

Suppose a selected pair relation has primitive direction
\(u=(x-y)/g(x,y)\), with a specified component \(u_W\ne0\)
and \(\operatorname{rad}(|u_W|)\le n^B\) for a fixed \(B\).
After (BC3) the actual difference has the form
\[
             |W(x)-W(y)|=A\,S\,U,\qquad
             \operatorname{rad}(S)\le n^B,\quad U\le X^\eta,
\]
where \(A\) is supported on the pool \(P_x\), \(S=|u_W|\),
and \(A\) has subpower many choices conditional on \(x\).
Here \(U\) can have unrestricted radical. A scalar normalization
divides both pure letters and their difference by
\(d=\gcd(W(x),W(y))\). With \(h=\gcd(d,A)\), one has
\(d/h\mid S U\). Put \(j=\gcd(d/h,U)\). Then
\((d/h)/j\mid S\), so the normalized difference has coefficient
\(A/h\), tracked factor \(S/((d/h)/j)\), and cofactor \(U/j\).
One can fix \(d\) at divisor cost conditional on \(W(x)\).

The short-cofactor two-letter estimate (UC1) applies to this
presentation. Its threshold is chosen for the fixed radical
cap and count loss, before selecting \(k,D\) and the input
budget threshold. The fixed-support entropy loss for \(A\)
has constants depending on \(k\), with zero inner exponent
for every fixed \(k\). This gives a complete-dictionary
reduction on bounded-direction relations without a bound
on the raw radical of the scaling gcd.
[Short-cofactor transfer (UC1–UC2)](ABC_UNRECORDED_COFACTOR_SECANT_TRANSFER.md)

The anchor argument is related to the polynomial-witness step of
Walsh's inverse sieve. Walsh obtains low-degree structure from
few occupied residue classes across sufficiently many primes
for a large set of points. Here the hypothesis on a star is
the numerical least common multiple of its actual pair gcds,
which directly forces exact vanishing. The saving on its zero
set uses the three original radical caps. These bounds do not
force a near-full relation with bounded primitive direction
radicals.
[Walsh, The inverse sieve problem in high dimensions, Theorem 1.1 and Proposition 2.2](https://arxiv.org/abs/1105.1551)
· [Bernert, Browning, Lichtman and Teräväinen, Proposition 2.1 and Lemma 2.2](https://arxiv.org/html/2410.12234v2#S2)
