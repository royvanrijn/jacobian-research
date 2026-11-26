# Recovery from residues at several projection orders

Complementary projected root ratios determine primitive rows at a
single modulus. A logarithmic list of projection orders can recover
almost all the original prime support even if every residue has the
same unknown unit multiplier. The modulus then only needs to be a
fixed power of the radical scale. A separate construction uses weak
limits of the exact valuation measures to make the individual ratios
much shorter on diffuse supports.

## A single mixed residue at square-root precision

Take positive pairwise coprime \(U,V,T\le X^B\) in a primitive
signed additive equation with fixed signs. Write \(n=X^\lambda\),
where the budget \(\lambda>0\) is fixed before the height limit.
Let \(P\) be a prescribed prime pool of prime product at most
\(X^{B_P}\). Define the outside projected root, radical and remainder
by
\[
 x_j(W)=\prod_{p\notin P}p^{\lfloor v_p(W)/j\rfloor},
 \qquad z(W)=\operatorname{rad}_{\overline P}(W),\qquad
 E_j(W)=\prod_{p\notin P}p^{v_p(W)\bmod j}.
\]
The part \(W_P\) of \(W\) supported on \(P\) satisfies
\(W=W_Px_j(W)^jE_j(W)\). Impose bounds \(U\le X^{h_U}\) and
\(z(W)\le n^{r_W}\) on \(W=U,T\). The integers \(D\mid V\)
and all data of the record are prescribed during a fibre count.
Write \(U\equiv sT\bmod D\), where \(s\in\{-1,1\}\) is
specified by the signs.

**MR1.** For \(j\ge2\), a prescribed unit residue \(x_j(U)/z(T)\)
modulo \(D\) admits
\[
 \ll_{B,B_P,\varepsilon}
 X^\varepsilon
 \left(1+\frac{X^{h_U/j}n^{r_T}}{D}\right)
 \left(1+\frac{n^{(j-1)r_U}}{D}\right)
                                                               \tag{MR1}
\]
original rows. In particular, for \(B=1\), \(r_U,r_T\le1\) and
\(j=\lceil\lambda^{-1/2}\rceil\), any \(D>X^{2\sqrt{\lambda}}\)
suffices at small fixed \(\lambda\).

Any two distinct fractions in the fixed residue have determinant
of absolute value at least \(D\). Scale the two coordinates by their
height bounds \(X^{h_U/j}\) and \(n^{r_T}\). The sine of the angle
between two scaled vectors is at least \(D/(2X^{h_U/j}n^{r_T})\)
when more than one fraction occurs. This gives the first factor.
Coprimality recovers \(x_j(U)\) and \(z(T)\) from the exact fraction.
Fix \(T\) on the support \(P\cup\operatorname{Supp}(z(T))\), and
fix \(U_P\), at uniform subpower cost. The original congruence fixes
\(E_j(U)\bmod D\), and \(E_j(U)\le n^{(j-1)r_U}\), giving the
second factor. Now both original free letters and the equation
recover the row. No fullness condition is needed on \(D\).
[Projective residue spacing in TL2](ABC_TWO_TARGET_SUMMAND_LOCKS.md#the-mixed-target-defect-count)
· [Support enumeration at polynomial height](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

## A logarithmic list of twisted residues at radical precision

For an integer \(Q=2^l\ge64\), let \(\mathcal J\) contain the
powers \(2^i\) for \(0\le i\le l\) and the integers
\[
     \left\lfloor\frac{2^iQ}{2^i+1}\right\rfloor,\quad
     \left\lceil\frac{2^iQ}{2^i+1}\right\rceil,\qquad
                         1\le i\le l .
\]
Include each \(j\in\mathcal J\) and \(2j\) in the order list.
The observed unit \(K_j\bmod D\) at an order in that list is
assumed to equal \(\chi x_j(U)/x_j(T)\bmod D\) for some
common unit \(\chi\) that is allowed to vary with the row.
Only the contrasts
\[
                 C_j=\frac{K_j/K_{2j}^2}{K_1/K_2^2}
                          \pmod D,\qquad j\in\mathcal J
\]
are used as records. In particular no choice of \(\chi\) is needed.
Write \(\varrho=r_U+r_T\).

**MR2.** If \(D>4n^\varrho\), each prescribed contrast record
has at most
\[
                      O_{B,B_P,Q,\varepsilon}
                  \bigl(n^{64B/(\lambda Q^2)}X^\varepsilon\bigr)
                                                               \tag{MR2}
\]
original rows. For \(Q\) the least power of two at least
\(\lambda^{-2/3}\), there are \(O(\log(1/\lambda))\) orders, every
order is \(O(\lambda^{-2/3})\), and the count per record is
\(n^{O_B(\lambda^{1/3})}X^\varepsilon\).

For a nonnegative integer \(e\), put \(b_j(e)=\lfloor e/j\rfloor
\bmod 2\), with values zero or one. The identity
\[
       \lfloor e/j\rfloor-2\lfloor e/(2j)\rfloor=b_j(e)
\]
shows that \(C_j\) reduces an exact rational number whose
\(p\)-valuation is \(b_j(e)-b_1(e)\) on \(U\), its negative
on \(T\), and zero at primes of \(P\). The common twist
cancels. This rational number has numerator times denominator
at most \(z(U)z(T)\le n^\varrho\). For each of the finitely
many contrasts, box both the numerator and denominator into
dyadic intervals at inner subpower cost. The product of their
upper endpoints is at most \(4n^\varrho\) in every nonempty
box. Any two fractions in a box that give the recorded
residue have determinant divisible by \(D\), with absolute
value less than \(D\). Thus the rational contrasts are exact.

The supports of the exact contrasts detect every outside prime
on \(U,T\) whose exponent \(e\) lies in \([1,Q^2/32]\). Here is
the parity argument. If \(b_{2^i}(e)\) are constant for
\(0\le i<l\), then \(e=kQ\) or \(e=kQ-1\) for a positive
integer \(k\). In the first case the common bit is zero, and
in the second it is one. If \(k\) is odd, \(b_Q(e)\) is
different from that common bit. If \(k\) is even, put
\(i=v_2(k)\) and \(j_*=2^iQ/(2^i+1)\). Then \(i\le l\) and
\(m=k+k/2^i\) is an odd integer. For \(e=kQ\),
\[
           0\le\frac{kQ}{\lfloor j_*\rfloor}-m
             <\frac{(3/2)k}{2Q/3-1}<1 .
\]
Thus \(b_{\lfloor j_*\rfloor}(e)=1\). For \(e=kQ-1\),
\[
             0<m-\frac{kQ-1}{\lceil j_*\rceil}
               \le\frac{(3/2)k+1}{2Q/3}<1 .
\]
Thus \(b_{\lceil j_*\rceil}(e)=0\). Both terminal inequalities
hold since \(k\le Q/32+1\le Q/16\). If the lower dyadic bits
are already unequal, one of them differs from \(b_1(e)\).
This proves the detection claim, including \(e=1\).
Dyadic probes alone up to \(Q\) miss the exponent \(2Q-1\),
whose bits on those probes are all one.

The union of the exact contrast supports and \(P\) is now
a fixed pool of polynomial prime product. Every prime
outside it on \(U,T\) occurs at exponent \(>Q^2/32\).
Since \(UT\le X^{2B}\), the product of those remaining
primes is at most \(X^{64B/Q^2}\). Enumerate that product
and then enumerate the two complete terms on the resulting
support at subpower cost. The original equation recovers
the third term. This proves (MR2).

## Selecting a divisor of small dictionary cost

In a height-first limiting box suppose the \(V\)-measure has
positive nonatomic mass in a compact band of positive scaled
degrees. Choose a compact nonatomic carrier of positive mass
in that band and retain its shrinking continuity neighborhoods.
The maximum normalized root mass of one exact slot on those
neighborhoods tends to zero. Indeed, otherwise a subsequence of
such slots converges to a point of the carrier and gives an
atom there. All retained slots have scaled degrees between
two fixed positive numbers. To obtain a divisor \(D\mid V\)
of numerical size larger than \(X^{2\sqrt\lambda}\) in (MR1)
or \(4n^\varrho\) in (MR2), add dyadically boxed full slots
from the carrier until their product has reached the desired
lower bound, with a strict margin. The cost of the raw
radical before the last slot is bounded by the target
numerical exponent divided by the lower scaled degree.
The last slot costs \(o(1)\) normalized mass. Thus the
block has at most \(n^{o(1)}\) values, up to the inner
subpower allowance. This argument uses the complete
precisions of the exact slots. The order and width limits
are taken at fixed \(\lambda\) before \(\lambda\to0\).
[Exact valuation boxing in AP3](ABC_AFFORDABLE_DIFFUSE_PROJECTED_PLANES.md#a-single-mixed-residue-suffices)
· [Anatomic factorization, BBLT Lemma 2.2](https://arxiv.org/html/2410.12234v2#S2)

For an actual balanced endpoint with such a divisor source, (MR1)
with the indicated growing order and (MR2) each require a joint
dictionary of \(D\) and the corresponding residue record of size
at least \(n^{1/3-o(1)}\). The second statement uses fewer than
a constant times \(\log(1/\lambda)\) projection orders and
remains valid with the common unit twist. An economical cover
by these records would give a saving. The size of the covering
dictionary is an additional arithmetic obligation.


## Arbitrarily small radical precision on the nonatomic support

A larger list of probes consists of adjacent parity contrasts.
Each such contrast is supported only on the slots where the
integer quotient by the order changes parity. Coincident changes
in two compact degree ranges have an arithmetic restriction
that is effective even when one of the two ranges has atoms.

Use a height-first sequence of exact valuation boxes with
\(\lambda\to0\). Write \(\nu_\lambda=\sum_e w_e\delta_{\lambda e}\)
for the sum of the normalized root-mass measures of \(U,T\)
outside the initial prescribed pool. Assume the total masses
and first moments are bounded and that the measures converge
weakly to \(\nu\). The rates in this argument are outer
allowances. Within each box, take the upper dyadic root sizes
as the weights and absorb their inner vanishing differences.
As usual one may truncate exact slots at a cutoff depending on
\(\lambda\), recording the aggregate support of still higher
orders at vanishing normalized cost, before taking the height
limit. Alternatively their total contribution to each parity
contrast is bounded by their total radical mass. Write
\(c_{\rm na}=\|\nu_{\rm na}\|\).

**MR3.** For any fixed \(\delta,\zeta>0\) and sufficiently
small budgets along this sequence, a list of contrasts on
\(O(\lambda^{-3/4})\) projection orders suffices as follows.
If \(D\mid V\) and \(D>n^\delta\), each prescribed list of
contrast residues modulo \(D\) supplies a fixed prime pool of
polynomial prime product containing at least
\(c_{\rm na}-\zeta-o(1)\) normalized radical mass from \(U,T\).
The common unit twist of all the original projected residues
is immaterial. The choice of list uses the box and its limiting
measures. The known pool can be produced after inner
subpower dyadic refinements.

To prove this, fix any \(\alpha\) with \(2/3<\alpha<3/4\),
and put \(J=\lceil\lambda^{-\alpha}\rceil\). For \(J\le j<2J\)
consider the exact adjacent contrast of parity masks
\[
 \frac{x_j(U)/x_j(T)}{\{x_{2j}(U)/x_{2j}(T)\}^{2}}\,
 \frac{\{x_{2(j+1)}(U)/x_{2(j+1)}(T)\}^{2}}
      {x_{j+1}(U)/x_{j+1}(T)} .
\]
It is read from the corresponding four residues as in (MR2).
The normalized radical mass of its support is bounded by
\[
        F_j=\sum_e w_e\,|b_j(e)-b_{j+1}(e)|.
\]

Choose a compact carrier \(K\) in a positive compact band
\([r,R]\), avoiding the atoms of \(\nu\), with
\(\nu(K)>c_{\rm na}-\zeta\). This is possible by removing
the countable atom set and using inner regularity of the
nonatomic part on bounded positive scaled degrees. If the
stated lower bound is negative the result is automatic. Use
shrinking continuity neighborhoods of \(K\), with the width
decreasing sufficiently slowly, to select subweights
\(u_e\le w_e\) of that compact part so their measures converge
weakly to \(\nu|_K\), a nonatomic measure. One way to make the
choice is to take successively smaller continuity
neighborhoods and a diagonal subsequence slow enough for
weak convergence of each restriction. Retain the exact slots
in the neighborhoods. A slightly enlarged positive compact
band still contains every retained degree for small \(\lambda\).
Put
\[
        G_j=\sum_e u_e\,|b_j(e)-b_{j+1}(e)|.
\]
Then
\[
                   \sum_{J\le j<2J} F_jG_j
                            =o\bigl((\lambda J)^{-1}\bigr).       \tag{MR3a}
\]

First restrict the \(F_j\) partner to degrees in
\([r',R']/\lambda\), with \(0<r'<R'\) fixed. For these
degrees and the retained \(G_j\) degrees, the consecutive
floors differ by at most one because \(J^2\lambda\to\infty\).
A hit from degree \(e\) therefore means there is a positive
integer \(k\) with
\[
                         j\le e/k<j+1 .
\]
For partner degree \(f\), a simultaneous hit gives a positive
integer \(l\) with \(j\le f/l<j+1\). Thus
\[
   \left|e/f-k/l\right|=O_{r,R,r',R'}(1/J),\qquad
           1\le l\ll_{R'}(\lambda J)^{-1}.
\]
For fixed \(e,f\) there is at most one reduced fraction \(k/l\)
possible for any joint hit, because two distinct fractions
with those denominator bounds are separated by
\(\gg_{R'}(\lambda J)^2\) and \(\lambda^2J^3\to\infty\).
If that fraction has reduced denominator \(d\), the number
of joint hits is \(O_{R'}((\lambda Jd)^{-1})\). In fact
\(l\) must be a multiple of \(d\) and each possible \(l\)
can give only \(j=\lfloor f/l\rfloor\).

For any fixed denominator cutoff \(M\), fractions with
\(d\le M\) lie in a finite rational set, since also
\(e/f\) stays in a fixed compact interval. The weighted
\((e,f)\) mass with those fractions tends to zero. Indeed
the two scaled degree measures converge weakly on these
bounded ranges and their product has zero mass on each of
the finitely many rational lines, as the first marginal is
nonatomic. The shrinking \(O(1/J)\) neighborhoods therefore
have vanishing mass. One can take \(r',R'\) at continuity
points for the second marginal. The contribution of larger
denominators is \(O_{R'}((\lambda JM)^{-1})\) times a bound
for the total product mass. Send \(M\to\infty\) after the
small-budget limit to bound the compact contribution.

For \(f<r'/\lambda\) the number of parity changes is at most
\(\lfloor f/J\rfloor-\lfloor f/(2J)\rfloor\le r'/(\lambda J)\).
Multiplying the sum of those weights over \(j\) by the bound
for \(\max_j G_j\) contributes \(O(r'/(\lambda J))\).
For \(f>R'/\lambda\) the first-moment bound makes the total
partner mass \(O(1/R')\). Since \(\sum_jG_j=O((\lambda J)^{-1})\)
on the fixed carrier, the contribution is
\(O((R'\lambda J)^{-1})\). Send \(r'\downarrow0\) and
\(R'\uparrow\infty\) after the compact limit. This proves (MR3a).

Take only the contrasts indexed by \(j\) with \(F_j\le\delta/2\).
The sum of \(G_j\) on the rejected \(j\)'s is
\(o((\lambda J)^{-1})\) by (MR3a). Every retained degree
\(e\asymp 1/\lambda\) has
\(\lfloor e/J\rfloor-\lfloor e/(2J)\rfloor
\gg_r(\lambda J)^{-1}\) distinct hit positions. It follows
that the mass of carrier slots hit only by rejected indices
tends to zero. In any accepted contrast, numerator times
denominator is at most \(n^{\delta/2+o_X(1)}\). Recover
the exact rational number by dyadic boxing and the
determinant argument of (MR2), using \(D>n^\delta\).
At a hit the contribution of every prime in its slot has
absolute valuation one in some recovered contrast.
The union of the recovered rational supports therefore
contains the claimed carrier mass on every row. Its prime
product divides \(UT\) on every row in the resulting
nonempty record, so its height is polynomial.

In a balanced endpoint with the usual total radical cap,
the zero price on the recovered pool gives a per-record
exponent at most \((1-c_{\rm na})/3+\zeta+o(1)\)
by (EB5). In particular \(c_{\rm na}\ge1/3\) for any pair
of free roles, so the coefficient is at most \(2/9+\zeta+o(1)\).
If both free measures are wholly nonatomic, one can instead
fix the residual support on \(U,T\) at arbitrarily small
cost, then count original powers on the enlarged support
at subpower cost. The counted residue record does not include
the powers of the recorded primes.
[Prices on a prescribed polynomial-height prime pool](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)
· [Balanced endpoint laws and their pairwise nonatomic mass](ABC_AFFORDABLE_DIFFUSE_PROJECTED_PLANES.md#the-projected-plane-criterion)


## A bounded number of small-modulus probes

**MR4.** The support conclusion of (MR3) can be obtained with
a number of contrast indices bounded independently of
\(\lambda\) along each limiting profile. The bound may depend
on \(\delta,\zeta\) and the limiting measures. The largest
actual projection order is \(O(\lambda^{-\alpha})\) for any
fixed \(1/2<\alpha<1\). As in (MR3), the contrast list is
chosen from the box before counting any of its residue fibres.

It suffices to take \(0<\zeta<1/4\). Use the weights
\(u_e,w_e\) of (MR3), this time choosing
\(\sum u_e>c_{\rm na}-\zeta/3-o(1)\) and a fixed compact
band \([r,R]\) for the support of the \(u_e\). Take
\(Q=\lceil\lambda^{-\alpha}\rceil\), and fix a small
\(\xi>0\) before the budget limit. Set
\(h=\lfloor \xi\lambda Q^2\rfloor\) and compare the parity
masks at \(q\) and \(q+h\), choosing \(q\) uniformly among
\(Q,\ldots,2Q-1\). In this proof expectations are over that
choice. Let
\[
 I_e(q)=|b_q(e)-b_{q+h}(e)|,\qquad
 G_q=\sum_e u_e I_e(q),\qquad F_q=\sum_e w_e I_e(q).
\]
The contrast again uses \(K_q,K_{2q},K_{q+h},K_{2(q+h)}\)
and cancels the common twist. For sufficiently small
\(\lambda\) we have \(h\ge \xi\lambda Q^2/2\) and \(h=o(Q)\).
For \(e\) in the carrier range,
\[
     a:=\xi r/12\ \le\ e/q-e/(q+h)
                      =\frac{eh}{q(q+h)}\ \le\ A:=\xi R .
\]
Take \(A<1\). For \(\Delta=e/q-e/(q+h)<1\) the quotient
floors differ by one exactly if \(\{e/q\}<\Delta\), where
the strict inequality also handles the case of equality
with the lower endpoint \(e/(q+h)\) integral. Thus the
hit lies between the phase events \(\{e/q\}<a\) and
\(\{e/q\}<A\).

The phase estimates are uniform in the exact slots. For
real \(c\ne0\) put \(e(t)=\exp(2\pi i t)\) in the next
display only. A Riemann sum and integration by parts give
\[
 \left|\frac1Q\sum_{Q\le q<2Q}e(c/q)\right|
          \ll \frac{|c|}{Q^2}+\frac{Q}{|c|}.                \tag{MR4a}
\]
For the first term compare the sum with the integral of
\(e(c/(Qt))\) on \(1\le t\le2\), using the bound
\(O(|c|/Q)\) on its derivative. For the second substitute
\(y=1/t\) in that integral and integrate against the
bounded smooth weight \(1/y^2\) on \(1/2\le y\le1\).
Since \(\lambda Q^2\to\infty\) and \(\lambda Q\to0\), (MR4a)
gives uniform phase equidistribution for the single
phases \(e/q\) on the carrier. More explicitly one can
minorize and majorize the interval indicators by continuous
periodic functions with integrals within any fixed
tolerance of the interval lengths, then approximate the
functions uniformly by trigonometric polynomials. Each
nonzero frequency is a fixed nonzero multiple of \(e\),
so (MR4a) applies. In particular \(\mathbb E I_e(q)\ge a/2\)
on the carrier for late budgets, and
\(\mathbb E G_q\ll \xi R\).

For the mixed expectation first restrict the partner
\(f\) to a compact degree range \([r',R']/\lambda\), with
\(r'<R'\) fixed. Its hit is contained in the phase interval
\(\{f/q\}<\xi R'\), provided \(\xi R'<1\). For \(e,f\)
in these compact ranges, every fixed joint Fourier
frequency \(c=a_1e+a_2f\) satisfies the upper bound needed
in (MR4a). The lower bound holds uniformly after excluding
\(|\lambda(a_1e+a_2f)|\le\theta\) for a fixed \(\theta>0\)
and each nonzero frequency pair \((a_1,a_2)\) in the
chosen polynomial majorants. These are finitely many
thin neighborhoods of rational lines. Their weighted
mass under \(u_e w_f\) can be made arbitrarily small,
since the limiting carrier marginal is nonatomic and
the other marginal is supported on a positive compact
band. This applies also when a coefficient \(a_1\) or
\(a_2\) is zero. First take the budget limit at fixed
\(\theta\), then shrink \(\theta\) and the majorant
integral errors. Consequently the compact partner
contributes \(O(\xi^2 R R')\) to the upper limiting
mixed expectation.

For \(f<r'/\lambda\), a parity change implies
\(q\le f/l<q+h\) for a positive integer \(l\le f/Q\).
Each \(l\) allows at most \(h+1\) positions, even when
some intervals overlap. Thus this hit has frequency
at most \((h+1)f/Q^2\le \xi r'+o(1)\) uniformly.
This range contributes \(O(\xi r')+o(1)\), using bounded
\(\max G_q\). For \(f>R'/\lambda\) the total partner
root mass is \(O(1/R')\), giving
\(O(\xi R/R')\) by the bound on \(\mathbb E G_q\).
Therefore
\[
       \limsup_{\lambda\to0}\mathbb E_q(F_qG_q)
          \ll \xi^2RR'+\xi r'+\xi R/R'.                     \tag{MR4b}
\]
The implicit constants can depend on the total mass and
first-moment bounds. The compact product estimate uses
upper interval lengths in the limit, so its constant
does not increase when \(r'\) decreases.

Reject \(q\) with \(F_q>\delta/2\). The carrier mass of
slots whose hit probability on the rejected indices is
at least \(a/4\) is bounded by
\[
                 \frac{8}{\delta a}\,\mathbb E_q(F_qG_q).
\]
Choose \(r'\) sufficiently small and \(R'\) sufficiently
large first, in terms of \(\delta,\zeta,r,R\), and then
take \(\xi\) sufficiently small. By (MR4b) the upper
limiting bound on this lost mass is less than \(\zeta/3\),
with a strict margin. The remaining carrier slots have
hit probability at least \(a/4\) in the nonrejected
portion of the original sample space. Choose
\(m=O(a^{-1}(1+\log(1/\zeta)))\) independent draws in the
original sample space and drop rejected draws. The expected
remaining missed carrier mass can be made less than
\(\zeta/3\), by bounding each such slot's miss probability
by \((1-a/4)^m\). Fix a choice of at most \(m\) nonrejected
probe positions with that bound. Together with the two
earlier carrier losses this covers
\(c_{\rm na}-\zeta-o(1)\) mass. Every selected contrast has
numerator times denominator at most
\(n^{\delta/2+o_X(1)}\), so exact recovery at \(D>n^\delta\)
uses the same dyadic determinant bound as (MR3). The list
is fixed by the box weights at the current budget, and all
sample sizes and order lists are fixed before the height
limit. The argument works equally if the finitely many
choices of probe lists are included as part of the record
at fixed \(\lambda\).
