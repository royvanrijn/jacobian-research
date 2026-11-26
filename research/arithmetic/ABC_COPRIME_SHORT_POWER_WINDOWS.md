# Coprime powerful integers in primitive short-power windows

For each fixed integer \(k\ge2\) and each fixed \(m\ge1\), a
positive proportion of integer centres \(q\) have \(m\) distinct
pairwise coprime \(k\)-full integers strictly between
\((q-1)^k\) and \(q^k\), all coprime to \(q\). This gives actual
primitive additive triples with two capped powerful roles and
a short positive displacement. One can also prescribe a growing
prime-power full divisor on every displacement, at a specified
small cost in the count of centres.

Narumi and Tachiya, Theorem 1 and Lemmas 1–2, give lower bounds
for occupancy between successive powers by simultaneous
approximation to fixed \(k\)-free classes. Their theorem allows
occupancy growing with the centre. It does not impose
coprimality of the values or coprimality to the centre. We use
the corresponding one-sided approximation and impose a
coprimality sieve for a fixed number of values.
[Narumi–Tachiya, arXiv:2611.20656v1](https://arxiv.org/html/2611.20656v1#S2)

## A positive-density primitive construction

Fix distinct primes \(a_1,\ldots,a_m\) not dividing \(k\). Put
\[
 \alpha_i=a_i^{-1-1/k},\qquad
 \ell_i(q)=\lfloor q\alpha_i\rfloor,\qquad
 T_i(q)=a_i^{k+1}\ell_i(q)^k .
\]
Fix \(0<\kappa<1\). For all sufficiently large \(Q\), at least
\(cQ\) integers \(q\in[\kappa Q,Q]\), with \(c>0\) depending
on the fixed data, obey
\[
 \begin{split}
 &0<\{q\alpha_i\}<\alpha_i/2\quad(1\le i\le m),\\
 &\gcd(q,T_i(q))=1,\qquad
 \gcd(T_i(q),T_j(q))=1\quad(i\ne j).
 \end{split}                                                   \tag{CPW1}
\]
The first inequality gives \(q-1/2<T_i(q)^{1/k}<q\).
Each \(T_i(q)\) is \(k\)-full and its valuation at \(a_i\) is
\(1\bmod k\), so it is not an exact \(k\)-th power.

To see the required independence, work over
\(\mathbb Q(\zeta_k)\). A prime above \(a_i\) has valuation
one on \(a_i\) and zero on all the other \(a_j\), since
\(a_i\nmid k\). Their classes therefore generate a copy of
\((\mathbb Z/k\mathbb Z)^m\) in the group of \(k\)-th power
classes. Kummer theory supplies an automorphism multiplying
\(a_i^{1/k}\) by \(\zeta_k\) and fixing the other roots for
each \(i\). Applying these automorphisms to a putative linear
relation proves that \(1,\alpha_1,\ldots,\alpha_m\) are
linearly independent over \(\mathbb Q\). In particular, all
\(\alpha_i\) and all ratios \(\alpha_j/\alpha_i\) for \(i\ne j\)
are algebraic irrational numbers.

Fix a prime cutoff \(P\) exceeding all \(a_i\). For primes
\(r\le P\) not among the \(a_i\), allow at most one of
\(q,\ell_1(q),\ldots,\ell_m(q)\) to vanish modulo \(r\).
For \(r=a_i\) require \(q\) and all \(\ell_j(q)\), \(j\ne i\),
to be nonzero modulo \(r\). Write \(D\) for the product of
these primes. In each progression \(q\bmod D\) the vector
\((\{q\alpha_i/D\})_{1\le i\le m}\) is equidistributed by
Weyl's criterion. For each \(i\), the integer part of
\(D\{q\alpha_i/D\}\) gives \(\ell_i(q)\bmod D\), and its
fractional part gives \(\{q\alpha_i\}\). The residue is thus
uniform and independent of the fractional-part test in the
limit. The local factor at \(r\notin\{a_1,\ldots,a_m\}\) is
\[
 (1-1/r)^{m+1}+(m+1)r^{-1}(1-1/r)^m=1-O_m(r^{-2}),
\]
and at \(a_i\) it is \((1-1/a_i)^m\). All factors are
positive. Their product has a positive lower bound
independent of \(P\). Therefore after the small-prime tests,
at least \(c_0Q+o_{P,k,m,\mathbf a,\kappa}(Q)\) centres remain
with \(c_0>0\) independent of \(P\).

For the large-prime tests one can sum over every integer
\(d>P\) dividing two of \(q,\ell_1(q),\ldots,\ell_m(q)\).
For two floors write them as \(du,dv\). Then
\[
        u\ll Q/d,\qquad
        |(\alpha_j/\alpha_i)u-v|\ll1/d .
\]
For \(q=du,\ell_i(q)=dv\), the same statement has
\(\alpha_i\) as the irrational coefficient. Each \(u,d\)
allows \(O(1)\) centres because it fixes one floor or the
centre itself. Roth's theorem bounds \(d\) by
\(O_{\varepsilon}(Q^{(1+\varepsilon)/(2+\varepsilon)})\),
using \(\|\beta u\|\gg_\varepsilon u^{-1-\varepsilon}\) for
the relevant fixed algebraic irrational \(\beta\). Also
uniformly in real intervals \(I\subset[0,1)\),
\[
 \#\{1\le u\le T:\{\beta u\}\in I\}
             =T|I|+O_{\beta,\varepsilon}(T^\varepsilon).
\]
Indeed the next partial quotient of a convergent of
denominator \(q_j\) is \(O_{\beta,\varepsilon}(q_j^\varepsilon)\)
by Roth. Each full block of length \(q_j\) has bounded
interval discrepancy in unnormalized counts. Comparing its
rotation orbit with the rational orbit of the convergent
moves each point by at most \(1/q_j\) and leaves only a bounded
number near the boundaries. The greedy decomposition of an
arbitrary orbit length into these block lengths uses at most
the sum of the relevant partial quotients. There are
\(O(\log T)\) denominators, which proves the bound after
resizing \(\varepsilon\).
[Roth, Rational approximations to algebraic numbers](https://doi.org/10.1112/S0025579300000644)

It follows that the number excluded for each pair is
\[
 \ll Q\sum_{d>P}d^{-2}
       +\sum_{d\ll Q^{(1+\varepsilon)/(2+\varepsilon)}}
                                    (Q/d)^\varepsilon
 \ll Q/P+O(Q^{7/10}),
\]
taking \(\varepsilon>0\) sufficiently small. The constants are
independent of \(P\). First take \(P\) large relative to \(c_0\)
and the fixed number of pairs, then take \(Q\) large with that
\(P\) fixed. The primes of each \(a_i\) have already been
excluded from the other values by the small-prime conditions.
This proves (CPW1).

## A growing full shift divisor

For the following variant, also fix a prime \(p\nmid k\) and
choose the \(a_i\) above congruent to \(1\bmod p\). Fix a
sufficiently small constant \(\gamma>0\) depending on the fixed
data and set
\[
 j=\lfloor\gamma\log_p Q\rfloor,\qquad B=p^j,\qquad R=pB.
\]
For each \(i\), Hensel lifting gives \(b_i\equiv1\bmod p\)
with \(a_i^{k+1}b_i^k\equiv1\bmod R\). Instead of (CPW1)
on an unrestricted centre, simultaneously impose
\[
                 q\equiv1\bmod R,\qquad
                 \ell_i(q)\equiv b_i+B\bmod R.
\]
For \(j\ge1\), the binomial theorem and \(p\nmid ka_i b_i\)
give \(v_p(q^k-T_i(q))=j\). There are
\[
                         \gg Q^{1-(m+1)\gamma}                \tag{CPW2}
\]
centres in \([\kappa Q,Q]\) satisfying these congruences and
(CPW1). In particular \(B\) is the full part of each short
displacement on the support \(\{p\}\).

Here is the growing-modulus count. In a fixed progression of
step \(S\), the vector \((\{q\alpha_i/S\})_i\) follows a
translate of the rotation with frequency vector \(\boldsymbol
\alpha\). This rotation has uniform box discrepancy
\(O(T^{-\delta})\) on \(T\) successive iterates for some
\(\delta>0\) depending on the fixed \(\alpha_i\). To see it,
linear independence and the algebraic norm give
\(\|\mathbf h\cdot\boldsymbol\alpha\|\gg
(\max_i |h_i|)^{-(D_0-1)}\) for nonzero integer vectors,
where \(D_0=[\mathbb Q(\alpha_1,\ldots,\alpha_m):\mathbb Q]\).
Clear fixed denominators before taking the norm. The
other conjugates of the relevant linear form are
\(O(\max_i |h_i|)\). Bound each exponential sum by the
geometric-series estimate and apply the multidimensional
Erdős–Turán–Koksma discrepancy inequality with frequency
cutoff \(\lfloor T^{1/(2D_0)}\rfloor\).
The estimate is uniform in the progression translate and in
the box positions. [Niederreiter, Lemma 2](https://www.impan.pl/shop/en/publication/transaction/download/product/82627)
We may decrease \(\delta\) below \(1/16\).

First take the fixed small-prime sieve above, excluding \(p\)
from its product \(D\). At \(p\) all \(q,\ell_i(q),a_i\)
are already units. Use progressions modulo \(RD\) for \(q\)
and boxes for \(\{q\alpha_i/(RD)\}\) to prescribe the combined
residue of each floor and the one-sided closeness condition.
For fixed \(P\) the number of boxes is bounded independently
of \(Q\). The main term is a positive constant times
\(Q/R^{m+1}\), bounded below independent of \(P\) after
dividing by \(Q/R^{m+1}\), by the same local product.
The error is \(O_P((Q/R)^{1-\delta})\).

Next let \(P<d\le Q^\epsilon\) be a prime. For two floors
required divisible by \(d\), split \(q\equiv1\bmod R\)
into its \(d\) progressions modulo \(Rd\). In each one
there are at most \(d^{m-2}\) boxes for the floor residues.
The \(m-2\) untested floors each have \(d\) allowable
residues modulo \(Rd\), and the two tested residues are
already fixed. The boxes have side lengths \(O((Rd)^{-1})\).
For a pair involving \(q\), there is one progression modulo
\(Rd\) and \(d^{m-1}\) boxes. Both counts, per \(d\), have
main term \(O(Q/(R^{m+1}d^2))\) and total discrepancy at
most
\(O((Q/R)^{1-\delta}d^{m-2+\delta})\).
Summing gives the error
\[
               O\!\left((Q/R)^{1-\delta}
                         Q^{\epsilon(m-1+\delta)}\right).
\]
For \(d>Q^\epsilon\), use the unrestricted Roth bound already
proved, giving \(O(Q^{1-\epsilon}+Q^{7/10})\).
Choose \(0<\epsilon<\delta/(8(m+1))\) and
\(0<\gamma<\epsilon/(4(m+1))\). Both discrepancies and the
unrestricted large-divisor tail are
\(o(Q/R^{m+1})\). The middle-divisor main terms cost
\(O(Q/(PR^{m+1}))\). Take \(P\) large relative to the fixed
local density before sending \(Q\) to infinity. This proves
(CPW2).

## Scope for the small-budget problem

Set \(X=Q^k\) and take \(\kappa=2^{-1/k}\). The constructions
give actual positive primitive triples
\[
                   T_i+(q^k-T_i)=q^k,\qquad X/2\le q^k\le X .
\]
They have displacement at most \(k X^{1-1/k}\), and
\(\operatorname{rad}(T_i)\le a_i\ell_i(q)<q\) while
\(\operatorname{rad}(q^k)\le q\). Hence for each fixed \(L>0\)
and a fixed \(k\ge 3/L\) both tall roles satisfy the absolute
cap \(X^{L/3}\). If the cap is relative to the actual target
the same inequality holds with \(X\) replaced by \(q^k\).
With \(k=\lceil 3/L\rceil\), CPW1 gives \(X^{1/k}\) distinct
centres up to a constant factor, reaching \(1/3\) in the
limiting normalized two-cap coefficient. CPW2 gives exponent
\((1-(m+1)\gamma)/k\) on its specified charts. For fixed \(k\)
and \(m\), that approaches \(1/k\) as \(\gamma\downarrow0\).
Neither theorem imposes a small radical cap on the whole
displacement. The full divisor in CPW2 has numerical weight
\(\gamma/k+o(1)\) at fixed \(k\), which tends to zero in the
small-budget limit \(k\to\infty\). The chosen positive
\(\gamma\) may also depend on \(k,m\). The construction
therefore supplies no modulus of fixed positive numerical
weight in that outer limit. It gives no lower bound for the
original three-cap exceptional set.

## Uniform upper bounds once the shift part is larger

A fixed positive numerical modulus is effective in the weighted
single-power setting, even when its support varies. Let \(k\ge2\)
and \(0<\eta\le1\) be fixed. For nonzero positive integers \(A,B\),
a sign \(\epsilon\in\{1,-1\}\), and positive coprime \(x,y\),
put
\[
 h=|Ax^k+\epsilon By^k|>0,\quad
 H=\max(x,y),\quad M=\max(Ax^k,By^k),\quad K=\max(A,B).
\]
For any prescribed finite set \(S\) of primes not dividing
\(kAB\), put \(D_S(h)=\prod_{p\in S}p^{v_p(h)}\).

**CPW3.** Uniformly in \(A,B,S\), the number of these pairs with
\[
 \frac{h}{D_S(h)}\le M H^{-2-\eta},
 \qquad
 H>\max\!\left(C_\eta,\ K^{1/[3(k+2)k^2]}\right)
                                                               \tag{CPW3}
\]
is at most \(\exp(O_{k,\eta}(1+|S|))\), where \(C_\eta\)
is fixed. In particular this statement permits varying
polynomial-height coefficients, with the displayed lower
threshold on the roots. It uses the actual numerical outside
part \(h/D_S(h)\).

Orient the roots so \(x=H\). At every finite prime in \(S\)
dividing \(h\), both \(x,y\) are units since the pair is
primitive and the prime is good. The congruence exhibits a
simple root of \(A Z^k+\epsilon B\) modulo \(p\). For its lift
\(\beta_p\in\mathbb Q_p\), Hensel's lemma gives
\[
              |x-\beta_p y|_p=|h|_p .
\]
Use the two forms \(x-\beta_p y,y\) at that prime, and use
the coordinate forms at primes in \(S\) not dividing \(h\).
At the real place use coordinate forms when \(h\ge M/2\).
Otherwise the sign is negative, \(Ax^k\ge M/2\), and one
can use \(x-\beta_\infty y,y\) with
\(\beta_\infty=(B/A)^{1/k}>0\). Indeed
\[
 \frac{|x-\beta_\infty y|}{H}
        =|1-\beta_\infty y/x|
        \le |1-(\beta_\infty y/x)^k|
        =\frac{h}{Ax^k}\le \frac{2h}{M}.
\]
At each place divide each evaluated form by the local maximum
norm of the pair. Every resulting absolute value is at most
one, and the product over both forms and all places is at
most \(2h/(M D_S(h))\). The determinant of each form system has
absolute value one. The forms are obtained from a single list
of at most \(k+2\) algebraic forms, the coordinate forms and
those using a root of \(A Z^k+\epsilon B\). For each choice of
root orientations the degrees of their coefficient fields are
at most \(k\). The inhomogeneous heights are at most
\(K^{1/k}\), since every root has \(k\)-th power of rational
height at most \(K\).

For completeness, to meet the fixed weight hypotheses of
Evertse–Ferretti, write \(J=2(1+|S|)\). If
\(H>C_\eta\), the product is at most \(H^{-2-\eta/2}\).
Take the \(J\) nonnegative numbers minus the base-\(H\)
logarithms of the individual normalized values. Cap each at
\(2+\eta/2\), which still leaves sum at least \(2+\eta/2\).
Rescale them to sum exactly \(2+\eta/2\) and round down at
mesh \(\eta/(8J)\). The resulting sum is at least
\(2+3\eta/8\). There are at most
\(\binom{O_\eta(J)+J}{J}\) possible discretized lists, by
the bound on the total sum. For each list decrease weights
proportionally to sum \(2+\eta/4\). Their negatives are the
fixed exponents \(d_{iv}\) in the theorem, all nonpositive,
and all the individual normalized bounds still hold.
The determinant factor \(A_v\) in (3.6) is one here,
and the Galois maximum in (3.7) is redundant for the rational
pair. The height threshold of Theorem 3.1 is at most the
one stated in (CPW3), using dimension two, \(k+2\) as the
number of distinct forms, and \(k\) as the degree bound.
The number of proper linear subspaces per such system
is bounded in terms of \(k,\eta\). Every line contains
at most one rational projective point arising from a
positive primitive pair. There are at most
\((k+1)^{|S|}\) assignments of the finite-place form systems.
The symmetric orientation is another constant factor.
This proves (CPW3) and keeps the bound exponential in
\(|S|\) even when the local products are very small.
[Evertse–Ferretti, Theorem 3.1 and equations (3.1)–(3.7)](https://annals.math.princeton.edu/wp-content/uploads/annals-v177-n2-p04-p.pdf#page=11)

Here is one counting consequence in the original
height-first setting. Let \(n=X^L\) and consider an original
primitive additive triple with two letters of the form
\(Ax^k,By^k\) and third letter \(h\). Suppose
\(X/4\le M\le X\) and the ordered coefficient pair ranges
over at most \(n^{f+o_X(1)}\) values of height at most \(X\).
Allow \(S\) to be chosen from the row, with the above good
prime condition and \(\operatorname{rad}(D_S(h))\le n^\mu\).
Keep only primes of \(S\) actually dividing \(h\).
There are at most \(n^\mu\) choices of the resulting
squarefree support, and every choice has
\(O(\log X/\log\log X)\) primes. For the pairs below the
root-height threshold of (CPW3) use trivial enumeration on
the whole root pair, without enumerating \(S\). For rows
satisfying the numerical inequality in (CPW3) the resulting
count is
\[
 O_{k,L,\eta,\varepsilon}\!\left(
     n^{f+o_X(1)}
       \left(X^{2/[3(k+2)k^2]}+n^\mu X^\varepsilon\right)
                            \right).                          \tag{CPW4}
\]
The diagonal \(Ax^k=By^k\) is outside the positive-third-letter
setting. No coefficient equality is needed. A physical full
divisor \(D_S(h)\ge X^w\) with \(kw>2+\eta\) by a fixed
margin is sufficient for the numerical inequality in (CPW3),
since \(H\le X^{1/k}\) and \(h\le X\). In a short window
\(h\ll X^{1-1/k}\), it suffices that \(kw>1+\eta\) by a fixed
margin. When \(k\) is fixed at each \(L\) and \(k\asymp 1/L\),
the first term in parentheses costs \(O(L^2)\) in base \(n\).
More generally that term has normalized exponent
\(2/[3(k+2)k^2L]\), tending to zero whenever
\(k L^{1/3}	o\infty\) after the fixed-budget height limit.
Thus in this two-weighted-power sector a fixed positive
numerical lock can give a power saving wherever the actual
coefficient and support prices \(f+\mu\) stay below \(1/3\)
by a fixed margin. One must supply these records in the
primitive quotient before applying the estimate. Broad exact
valuation profiles need not have two such weighted powers
with a small joint coefficient dictionary. The theorem does
not give a uniform threshold in \(L\) independently of the
chosen \(k\) and the stated margins.


## An adaptive common-power sector with a third-letter lock

**CPW5.** Fix \(0<\delta<1/3\), \(w>0\), and an integer
lower bound \(K_0(L)\) with \(L^{1/3}K_0(L)\to\infty\)
as \(L\downarrow0\). In an original primitive height shell
\(X/2<C\le X\), suppose two letters have exact representations
\(Ax^k,By^k\) with a shared order \(k\ge K_0(L)\). Allow
\(k\) and the orientation to depend on the row. Let the
joint coefficient dictionary have exponent \(f\) in
base \(n=X^L\). Require a full divisor \(D\) of the third
letter with \(D\ge X^w\), \(\operatorname{rad}(D)\le n^\mu\)
and primes coprime to \(k\). Assume \(f+\mu\le1/3-\delta\).
Then for all sufficiently small fixed \(L\) the count in
this sector is
\[
 O_{L,\delta,w,\varepsilon}
                       \bigl(X^{(1/3-\delta/2)L+\varepsilon}\bigr),
\]
with the implied constant also allowed to depend on
the selected function \(K_0\) and the dictionary bound.
The original three radical caps and any repair
eligibility condition can be retained.

Indeed two of the letters in the original shell have
maximum at least \(X/4\); this holds for any
chosen pair, since at most the smaller summand can
lie below that height. The good-prime condition
of CPW3 is automatic on \(D\) apart from the stated
restriction on \(k\), by primitivity. For
\(K_0(L)\le k\le \lceil 8/(\delta L)\rceil\), use CPW4
with \(\eta=1/2\). At sufficiently small \(L\),
\(wK_0(L)>2+\eta\) and
\[
              \frac{2}{3(k+2)k^2 L}\le \delta/4.
\]
The set of such \(k\)'s is finite at each fixed budget,
before the height limit. The coefficient charge plus
the maximum of this small-root exponent and \(\mu\)
fits inside \(1/3-\delta/2\).
If \(k>\lceil 8/(\delta L)\rceil\), and at least
one of \(x,y\) exceeds one, then \(k\le\log_2 X\).
Enumerating both roots and the order costs at most
\(X^{\delta L/4+o_X(1)}\) after fixing the coefficients.
When both roots are one, the coefficients already
determine the entire row on that orientation.
This proves the bound with no uniform dependence
of the Subspace Theorem constants on varying \(k\)
at the height limit.

For instance, require the joint radical of the
two off-grid coefficients to have exponent less
than \(1/3\) by a fixed margin. Ordinary radical
enumeration gives the required coefficient
dictionary. If the remaining third letter has
positive nonatomic mass along a putative balanced
saturating sequence, a compact positive-order
block gives the affordable full divisor at an
arbitrarily small cost \(\mu>0\) and some fixed
\(w>0\). At each fixed budget finitely many relevant
orders \(k\) in the bounded-order branch have
only a fixed set of prime divisors. Remove them
from the selected block: their total absolute
logarithmic weight within a fixed scaled band
is height-subpower. The lower numerical margin
of the block persists. The large-order branch
already has the independent enumeration bound.
Consequently this sector is excluded even when
the weak measures of the two powered letters
are broad and the remaining letter has arbitrary
valuations outside the selected divisor.
[Cheap full divisors from diffuse valuation mass](ABC_MODULAR_DIGIT_SUPPORT_AND_ATOMIC_COMPLETION.md#completing-all-balanced-profiles-after-small-modulus-support-recovery)
· [Earlier two-order bounds with unequal orders](ABC_TWO_VALUATION_GRIDS.md#count-and-range)

## Paying support on all three letters

**CPW6.** The conclusions of CPW3 and CPW4 hold for disjoint
prime sets \(S_x,S_y,S_h\) carried by \(x,y,h\), respectively,
with \(S_h\) coprime to \(kAB\), after replacing \(D_S(h)\) by
\[
 D_xD_yD_h,\qquad
 D_x=\prod_{p\in S_x}p^{v_p(x)},\quad
 D_y=\prod_{p\in S_y}p^{v_p(y)},\quad
 D_h=\prod_{p\in S_h}p^{v_p(h)} .
 \tag{CPW6}
\]
At a prime in a root set use the coordinate forms \(x,y\).
Their normalized product pays the local absolute value of
that root, since the other root is a unit. Use the same
\(k+2\) possible forms and the combined set of
places in the rounding and quantitative count. For (CPW4)
the radical cap \(n^\mu\) is on \(\operatorname{rad}(D_xD_yD_h)\). Each
squarefree support has at most \(3^{\omega}\) allocations,
which is height-subpower at these caps.

If the full parts selected on \(Ax^k,By^k,h\) have numerical
weights at least \(w_U,w_V,w_h\), respectively, then a
sufficient condition is
\[
                     w_U+w_V+k w_h>2+\eta
\]
by a fixed margin. For \(h\ll X^{1-1/k}\) the right side is
\(1+\eta\). Here \(p\nmid AB\) on the selected root pools
makes their contributions to root height exactly \(1/k\)
times their weights on the complete values. If the larger
root is already known to exceed the height threshold
in (CPW3), the small-root term in (CPW4) can be omitted.
In particular, if one weighted-power letter is at least
\(X/2\) and its coefficient is at most \(X^{1-\chi}\)
for fixed \(\chi>1/[3(k+2)k]\), the larger root exceeds
that threshold at late heights even for fixed small \(k\).
