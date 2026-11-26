# Empirical prices for projective divisor records

A prime can be common to projective directions without dividing
any of the three original letters. Frequencies of actual projective
coincidences supply a price for such primes. The short-line lattice
bound converts the price and the numerical weight of the modulus
into a count of actual neighbours, with all three radical caps
retained.

## PF1. Directional support prices at a fixed anchor

Fix \(h_0,H,a_0,\eta>0\). At a sufficiently small fixed budget
\(L>0\), put \(n=X^L\). Let \(\mathcal F\) be a family of
\(S\ge2\) distinct positive primitive rows \(A+B=C\), with
\(X^{h_0}\le C\le X^H\) and
\(\operatorname{rad}(T)\le n^{a_0}\) for each role \(T\).
Fix an invertible real linear map \(\mathcal T\) with
\(\|\mathcal T(A,B)\|_\infty\le1\) on the family and set
\(\mathcal B=|\det\mathcal T|^{-1}\).
The ordinary complete-height box, an anisotropic box, and a
prescribed interval of ratios all have this form as in (DL2).

Write \(F_{uv}=A_u B_v-A_v B_u\) for two distinct rows. It is
nonzero by primitivity. Fix an anchor \(u\) in the family and let
\[
 \pi_u(p)=\frac{1}{S-1}\#\{v\in\mathcal F\setminus\{u\}:p\mid F_{uv}\},
 \qquad
 \mathcal C_u(D)=-\sum_{p\mid D}\log_n\max\{\pi_u(p),p^{-1}\}.
\]
For a positive integer \(D\) define
\[
                           v(D)=\bigl(\log_n(\mathcal B/D)\bigr)_+ .
\]
For every fixed \(b\ge0\) and every \(\varepsilon>0\), the number
of neighbours \(v\ne u\) admitting a divisor
\(D\mid |F_{uv}|\) with
\[
                       \mathcal C_u(D)+\tfrac12 v(D)\le b
                                                               \tag{PF1}
\]
is \(O(n^{b+\eta}X^\varepsilon)\), uniformly in the anchor and
the original family. The total ordered pair count is at most
\(O(Sn^{b+\eta}X^\varepsilon)\). Extra eligibility conditions,
including a least-repair cutoff, can be kept on the fixed original
family.

To prove this, put \(q_u(p)=\max\{\pi_u(p),p^{-1}\}\) for primes
\(p\le X^{2H}\). Distinct rows give minors of absolute value at
most \(X^{2H}\). Uniformly in \(u\),
\[
   \sum_{p\le X^{2H}}q_u(p)
     \le \frac{1}{S-1}\sum_{v\ne u}\omega(|F_{uv}|)
                              +\sum_{p\le X^{2H}}p^{-1}
     =o_X(\log X).
\]
The first term is \(O_H(\log X/\log\log X)\) by the primorial
bound, and the second is \(O_H(\log\log X)\).
A support of charge at most \(c\) contributes at least \(n^{-c}\)
to \(\prod_p(1+q_u(p))=X^{o_X(1)}\). Given its support, the
number of integers \(D\le X^{2H}\) supported on it has subpower
cost by uniform fixed-support counting. Thus there are at most
\(n^c X^{o_X(1)}\) candidate moduli of charge at most \(c\).
The proof prices a prime once, and further powers are contained
in the same supported-integer bound.
[Supported-integer and prime-price enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

For fixed \(u,D\), the condition \(D\mid F_{uv}\) puts all
eligible \(v\) in a homogeneous rank-two sublattice of the
\((A,B)\)-plane of determinant \(D\). In fact
\((A,B)\mapsto A_u B-B_u A\) maps \(\mathbb Z^2\) onto
\(\mathbb Z\) since the anchor is primitive. The preimage of
\(D\mathbb Z\) therefore has index \(D\). Its transformed
determinant is \(D/\mathcal B\). The generalized form of
(DL2) gives
\[
                       O(n^{v(D)/2+\eta'}X^{\varepsilon'})
                                                               \tag{PF2}
\]
for every prescribed inner allowance, after choosing
sufficiently small fixed budgets. This remains valid if
\(D/\mathcal B\) is below the thin-line range of (DL2). In
that case \(v(D)>L^{-\xi_1}\) and the trivial bound from two
radical-capped letters, \(O(n^{2a_0}X^{\varepsilon'})\),
is smaller for sufficiently small \(L\).
[Divisor-lattice line count and homogeneous sublattices (DL2–DL4)](ABC_DIVISOR_LATTICE_LINE_OCCUPANCY.md#dl2-slicing-the-actual-complete-value-lattice)

Bin the positive deficit \(v(D)\) in intervals of width
\(1/\log n\), with a separate zero bin. Condition (PF1)
only uses deficits at most \(2b\). In a bin with lower
endpoint \(v_0\), its moduli have charge at most \(b-v_0/2\),
and (PF2) bounds each contribution by
\(O(n^{v_0/2+\eta'} X^{\varepsilon'})\), up to a constant.
Summing the logarithmic number of bins gives (PF1).
The charge uses actual projective coincidences, which may
involve primes absent from every original coordinate.
[Charged full-divisor records (DL5)](ABC_DIVISOR_LATTICE_LINE_OCCUPANCY.md#dl5-a-charged-divisor-record-inequality)

For example, if \(S\ge n^d\) for fixed \(d>0\), then for
every fixed \(\sigma>0\), apart from
\(O(S^2 n^{-\sigma/2})\) ordered pairs for sufficiently small
fixed budgets and large heights, every permitted divisor
satisfies
\[
                        \mathcal C_u(D)+\tfrac12 v(D)>d-\sigma .
                                                               \tag{PF3}
\]
Take \(b=d-\sigma\) when this is nonnegative, and choose the
normalized allowance and \(\varepsilon/L\) smaller than
\(\sigma/4\). If \(d-\sigma<0\), the conclusion follows
from nonnegativity.

One consequence concerns the placement of prime support in an
actual small-radical minor or its divisors. Fix \(0\le\theta<1\).
Call a prime on \(D\) rare at \(u\) if
\(\pi_u(p)<p^{-\theta}\), and let \(r_{\rm rare}\) be the
normalized radical mass of those primes. With
\(r=\log_n\operatorname{rad}(D)\), the price satisfies
\(\mathcal C_u(D)\le\theta r+(1-\theta)r_{\rm rare}\).
Therefore the pairs in (PF3) require
\[
      r_{\rm rare}>
        \frac{d-\sigma-\tfrac12 v(D)-\theta r}{1-\theta}.
\]
In particular a conditional count using a high-weight projective
divisor whose primes mostly have high incidence faces a quantitative
rare-mass requirement. Extracting such divisors at low charge in
arbitrary near-saturating families remains open.

## PF4. Unconditioned projective dictionaries at full depth

For \(p^k\le X^{2H}\) let \(\mathcal P_{p^k}\) denote the
primitive projective residue classes in
\(\mathbb P^1(\mathbb Z/p^k\mathbb Z)\).
Write \(m_{p,k,\rho}\) for the number of original rows in
a class \(\rho\in\mathcal P_{p^k}\), and write
\(f_{p,k,\rho}=m_{p,k,\rho}/S\). Two primitive pairs
have the same projective class exactly when their
determinant vanishes modulo \(p^k\). Indeed a unit coordinate
of one vector determines the scalar multiple, which must
also be a unit by primitivity of the other vector.

For an ordered pair of distinct rows and \(D\mid |F_{uv}|\),
use their common class \(\rho_p\) at each \(p^k\parallel D\).
Define the unconditioned depth charge
\[
  \mathcal C(D,\boldsymbol\rho)
    =\sum_{p^k\parallel D}
           \min\{k\log_n p,-\log_n f_{p,k,\rho_p}\}.
\]
The number of ordered distinct pairs for which some such record
satisfies
\[
                 \mathcal C(D,\boldsymbol\rho)+\tfrac12 v(D)\le b
                                                               \tag{PF4}
\]
is \(O(n^{2b+\eta}X^\varepsilon)\) for each fixed \(b,\eta>0\)
at sufficiently small fixed budgets. The same proof allows
\(b=0\). This comparison is at the full recorded prime-power
depth.

For this dictionary use the distinct-pair class weights
\[
 w_{p,k,\rho}=\frac{m_{p,k,\rho}(m_{p,k,\rho}-1)}{S(S-1)},\qquad
 q_{p,k,\rho}=\frac{\max\{w_{p,k,\rho},p^{-2k}\}}{(k+1)^2}.
\]
The sum over \(\rho\) of \(w_{p,k,\rho}\) is the probability
that \(p^k\) divides a random ordered distinct-pair minor. In
particular it is bounded above by the corresponding probability
at \(k=1\). The squared decay in \(k+1\) gives
\[
 \sum_{p,k,\rho}\frac{w_{p,k,\rho}}{(k+1)^2}
      \ll \sum_p\Pr_{u\ne v}(p\mid F_{uv})=o_X(\log X).
\]
There are \(p^k+p^{k-1}\) primitive projective classes. One
chart normalizes the first coordinate to one, and the
complementary chart normalizes the second coordinate with the
first divisible by \(p\). Thus the baseline weights
\(p^{-2k}/(k+1)^2\) also sum to \(O(\log\log X)\).
Selecting at most one pair \((k,\rho)\) per prime yields a total
record weight bounded by
\(\prod_p(1+\sum_{k,\rho}q_{p,k,\rho})=X^{o_X(1)}\).

At each class used by an ordered distinct pair, \(m_{p,k,\rho}\ge2\),
so \(w_{p,k,\rho}\ge f_{p,k,\rho}^2/2\). A record of depth charge
at most \(c\) therefore has weight at least
\(n^{-2c} 2^{-\omega(D)}\tau(D)^{-2}\), where
\(\tau(D)=\prod_{p^k\parallel D}(k+1)\).
Since \(D\le X^{2H}\), both extra factors cost only a subpower.
There are at most \(n^{2c}X^{o_X(1)}\) records of charge
at most \(c\) that support pairs. The local linear functional
for a class has a unit coefficient, so its kernel has index
\(p^k\). CRT makes the simultaneous homogeneous lattice have
determinant \(D\). By (PF2), the pair count in this lattice is
\(O(n^{v(D)+2\eta'}X^{2\varepsilon'})\). Bin \(v(D)\) as
in (PF1) and sum over the records with \(c\le b-v(D)/2\).
This gives (PF4).

The two charges have different inputs. After the anchor is fixed,
the support charge of (PF1) needs no price for the prime-power lift
of a chosen direction. The unconditioned dictionary uses the
frequency of the exact class at that depth and prices all potential
directions. For a pair-supported class one has
\(\pi_u(p)\ge (m_{p,k,\rho_p}-1)/(S-1)\ge f_{p,k,\rho_p}/2\).
Thus the anchored support charge is at most the unconditioned
charge up to an inner subpower allowance. The latter proof gives
the smaller total count \(n^{2b+\eta}X^\varepsilon\) for its
own threshold when the family exponent exceeds \(b\).
Gallagher's larger sieve also uses coincidence counts of distinct
pairs at several depths. The estimates here sum individual
direction records by empirical prices and use the separate
three-cap lattice count within each record.
[Gallagher, *A larger sieve*, Theorem 2](https://matwbn.icm.edu.pl/ksiazki/aa/aa18/aa1819.pdf)

## PF5. Two caps can saturate the linear deficit

Fix \(\rho>0\) and sufficiently small fixed \(L>0\) with
\(\rho L<1\). Pick \(R=5^k\) with
\(n^\rho/50<R\le n^\rho/10\), and choose powers \(P_2\) of
\(2^4\) and \(P_3\) of \(3^4\) to be the largest such
powers at most \(X/(10R)\). Both are comparable to \(X/R\)
with absolute constants. For \(j\equiv2\bmod30\) in
\([R,2R]\) put
\[
            A_j=P_2j,\qquad B_j=P_3(j+R),\qquad C_j=A_j+B_j .
\]
These are genuine positive primitive rows, and all three values
are comparable to \(X\) and at most \(X\). The bounds on \(j\)
and \(R\) give \(\operatorname{rad}(A_j)\le4R<n^\rho\)
and \(\operatorname{rad}(B_j)\le9R<n^\rho\).
Primitivity follows since \(j\) is a unit at \(3,5\),
\(j+R\) is odd, and \(\gcd(j,j+R)=1\).
All three letters are units at \(5\), since
\(P_2\equiv P_3\equiv1\bmod5\) and \(j\equiv2\bmod5\).

There are \(\asymp R\) rows with distinct labels in each role.
Their projective minors and generated-lattice determinant are
\[
        F_{ij}=P_2P_3R(i-j),\qquad
        G_{\mathcal F}=30P_2P_3R\asymp X^2/R .
\]
Thus the ordinary linear deficit bound \(X^2/G_{\mathcal F}\)
has the true count scale \(R\). For every anchor the modulus
\(D=P_2P_3R\) costs zero in (PF1), since all its primes
appear in every off-diagonal minor. It also has zero
unconditioned charge at these depths. In particular the
large prime power \(5^k\) in this cheap projective record
comes from a prime outside all three original supports.

Here \(C_j\) runs through a unit progression of step
\(30(P_2+P_3)\) and index width \(O(R)\).
Indeed \(\gcd(P_2+P_3,P_3R)=1\), and \(C_j\) is a unit
at \(2,3,5\) by the chosen congruences. The short-block
estimate implies that only \(O(n^\eta X^\varepsilon)\)
of the \(C_j\) have radical at most \(n^\rho\), for any
fixed \(\eta,\varepsilon>0\) and sufficiently small fixed
budgets. The third cap is therefore essential in the
half-deficit improvement.
[Unit progression counts for arbitrary steps](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#arbitrary-span-progression-occupancies)

There is also a gap between the anchored support and
unconditioned depth prices on actual primitive rows.
Keep the same construction but replace \(j+R\) by \(j+5\)
in \(B_j\), keeping the same constants and range of \(j\).
The power of five in the common minor is now determined by
\(5(i-j)\). Modulo five all rows still have the same
projective direction. At a depth \(\ell\) with
\(5^{\ell-1}\asymp R^\gamma\), for any fixed \(0<\gamma<1\),
the classes are exactly the occupied \(j\)-classes modulo
\(5^{\ell-1}\). Each has relative frequency
\(\asymp5^{-(\ell-2)}\), since the fixed progression
has step divisible by five exactly once. They all
support distinct pairs for large \(X\). The record
\(D=5^\ell\) costs zero at either anchor in (PF1) and
has depth charge \(\gamma\rho+o_X(1)\) in base \(n\)
in (PF4). The two-cap example allows a genuine
difference between the two prices throughout a
polynomial range of depths.

## PF6. The fresh part must carry the projective price

In a common balanced profile with row exponent \(d=1/3\)
up to the height-first allowances, an ordered pair is
called fresh at \(p\) if \(p\) divides neither original
complete product \(ABC\). For any fixed \(\sigma>0\),
outside a power-smaller set of ordered pairs, every
\(D\mid |F_{uv}|\) satisfies
\[
           \mathcal C_u(D_{\rm fresh})
                       +\tfrac12 v(D)>d-\sigma .
                                                               \tag{PF6}
\]
Here \(D_{\rm fresh}\) is the full part of this divisor
on fresh primes and \(v(D)\) uses the weight of the
whole divisor. In particular a divisor within \(n^v\)
of the area scale requires fresh empirical price at
least \(d-v/2\) up to the chosen allowance, even when
the large numerical weight comes from original
nonfresh prime powers.

Fix the first row, give zero price to the primes of
its original complete product, and apply (DRR4) to
each of the three pairs of roles of the second row.
For any fixed positive threshold on shared radical
mass, at least two thirds of that mass belongs to some
pair of roles. The estimates therefore give a power loss
among second rows, uniformly in the anchor. After
summing over the anchors, outside a power-smaller
set the normalized shared radical mass is below any
prescribed small tolerance. For \(p\mid F_{uv}\), if
one original row has \(p\) on a complete letter then
the other row must have it on the same letter:
the three roles correspond to the distinct projective
classes \([0:1],[1:0],[1:-1]\) modulo \(p\).
Thus the nonfresh radical of every divisor \(D\) of
this minor is contained in the shared radical of the
pair. Its empirical price is at most its raw mass.
Subtract that tolerance in (PF3), taking the margins
smaller than \(\sigma\). This gives (PF6).
[Shared support in two roles (DRR4)](ABC_DIVISOR_RECORD_RECIPROCITY.md#anisotropic-intersections-for-an-arbitrary-record)
