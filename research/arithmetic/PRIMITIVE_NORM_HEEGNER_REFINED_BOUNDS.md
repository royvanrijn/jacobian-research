# Refined bounds and finite reductions for the remaining class-number-one norms

For \(D\in\{43,67,163\}\), put
\[
 Q_D(x,y)=x^2+xy+\frac{D+1}{4}y^2,\qquad
 P_D=\{Q_D(x,y)>0:\gcd(x,y)=1\}.
\]
Let \(C_D\) consist of 1 and the positive integers supported on primes
\(q\) with \((q/D)=1\). All sums below have **positive numerical values**
as summands. Define the explicitly finite set
\[
 E_D=\{2\le e\le1024:2\mid e,\ D\nmid e,\ e\notin C_D+C_D\}.
                                                               \tag{H.1}
\]
The checker evaluates this definition exactly; the complete sets appear
in Section 8.

**Theorem H1 (HeegnerCoreAndQuotientOnsets).** Set
\[
 X_{43}=2^{32},\qquad X_{67}=2^{35},\qquad X_{163}=2^{45}.
                                                               \tag{H.2}
\]
For each of the three discriminants:

1. Every even \(D\)-unit \(n\ge X_D\) all of whose prime divisors are
   inert in \(\mathbb Q(\sqrt{-D})\) belongs to \(C_D+C_D\).
2. If \(e\in E_D\), \(p\) is any split prime, and \(ep\ge X_D\),
   then \(ep\in C_D+C_D\). This includes \(p\mid e\).

**Theorem H2 (HeegnerFiniteClassificationReduction).** The following two
statements are equivalent:
\[
 n\in C_D+C_D\quad\text{for every even }2\le n\le X_D
       \text{ with }D\nmid n,\ n\notin E_D;                 \tag{H.3}
\]
\[
 C_D+C_D=\{n\ge2:2\mid n,\ D\nmid n\}\setminus E_D.        \tag{H.4}
\]
If they hold, define the finite set
\[
 R_D=E_D\setminus(P_D+P_D).
\]
Membership here only involves summands smaller than 1024. Then
\[
 P_D+P_D=\{N\ge2:2\mid N,\ D^2\nmid N\}
                    \setminus(R_D\cup D E_D).             \tag{H.5}
\]
Conversely (H.5) implies (H.4), so it has the same remaining finite
obligation.

| \(D\) | Exact analytic onset \(X_D\) | \(|E_D|\) | Candidate primitive exceptions \(|R_D\cup D E_D|\) | Largest candidate |
|---|---:|---:|---:|---:|
| 43 | 4,294,967,296 | 10 | 18 | 3182 |
| 67 | 34,359,738,368 | 21 | 38 | 9246 |
| 163 | 35,184,372,088,832 | 67 | 129 | 161696 |

**Scope.** The analytic onsets and the finite reduction are proved here.
The full finite conditions (H.3) have **not** been checked by this
packet. In particular, this note does not claim (H.4) or (H.5)
unconditionally for any of these three discriminants. The exact checker
covers the small range through 1024, not the ranges through \(X_D\).
The last range remains prohibitively large for the present segmented
sumset method. The sizes above are sufficient bounds, with no claim of
optimality or historical novelty.

The primitive value identity (H.6) is classical: it is also given by
Baake–Scharlau–Zeiner, *Similar Sublattices of Planar Lattices*,
[Proposition 5.2 and Corollary 5.3](https://doi.org/10.4153/CJM-2011-019-9),
printed pp.1227–1228. Those statements classify individual primitive
norm values; they do not supply the additive finite coverage required
in (H.3). They are context, while Section 1 proves the identity needed
here directly.

The improvements over the generic method in the
[discriminant −19 proof](PRIMITIVE_NORM_MINUS19_EXACT.md) are a lattice
estimate retaining its area and boundary terms, separate treatment of
inert primes dividing the target, and a size constraint on the Euler
factors in the spectral error. Classical spectral inputs remain the
same; there is no extra zero-free-region or leading-count hypothesis.

## 1. Primitive values and the spectral normalization

Write \(R=\mathbb Z[(1+\sqrt{-D})/2]\). The Minkowski bound
\((2/\pi)\sqrt D\) is less than 5, 6, and 9 respectively. Every
prime below the corresponding bound is inert. The only possible
integral ideal-class representatives have norm 1 or 4, and the ideal
of norm 4 is \((2)\). Thus the ideal class number is one. This uses
Milne, Theorem 4.3, with the elementary inequality \(\pi>3\);
the requisite residue symbols are checked exactly.

Primitive ideal factorization gives
\[
 P_D=C_D\ \dot\cup\ D C_D,\qquad
 Dm\in P_D+P_D\ \Longleftrightarrow\ m\in C_D+C_D.          \tag{H.6}
\]
An inert prime dividing a field norm divides both coordinates; a
primitive element has ramified norm valuation at most one. Conversely,
choose all of each split prime power on a single oriented prime ideal,
principalize the product, and optionally multiply by \(\sqrt{-D}\).
This produces a primitive generator. Also, \(C_D\) consists of odd
nonzero quadratic residues modulo \(D\). Since \(-1\) is nonsquare,
no sum of two such residues is zero modulo \(D\). These facts prove
both assertions in (H.6), including the converse at a ramified target.

Use the maximal order
\[
 B=(-D,-1)_{\mathbb Q},\quad O=R\oplus Rj,\quad
 j^2=-1,\quad ja=\bar a j,
\]
whose norm is \(F_D=Q_D\oplus Q_D\). Voight, Theorem 15.5.5 and
Example 15.5.7, apply to each prime \(D\equiv3\pmod4\): the algebra
has discriminant \(D\) and this order is maximal. Its units are
\(\{\pm1,\pm j\}\). Retain all right-ideal classes and their
inverse-projective-unit measure. The mass formula (Voight, Theorem
25.1.1) gives
\[
 \mathfrak m=\frac{D-1}{12},\qquad
 w=\frac{1/2}{\mathfrak m}=\frac6{D-1}.                  \tag{H.7}
\]
The principal class need not be the only class.

At an odd inert prime \(p\) the fine projective quotient is
\(\operatorname{PGL}_2(\mathbb F_p)\). At \(D\) it is
\(\mathbb F_{D^2}^{\times}/\mathbb F_D^{\times}\), of order \(D+1\).
All projective unit stabilizers inject. Indeed a noncentral unit has
trace \(-1,0,1\) and discriminant \(-3\) or \(-4\). Scalar reduction
at an odd unramified prime would force divisibility by \(p^2\), and
scalar reduction at the radical would force divisibility by \(D\).
Both are impossible, including at \(p=3\). The ramified fine quotient
therefore has \((D^2-1)/12\) points, each of probability
\[
                         \nu=\frac{12}{D^2-1}.           \tag{H.8}
\]
Exactly one global norm character survives at maximal level, and
exactly two at any one of these fine levels. Their eigenvalues for the
all-integral-ideal Hecke operator \(T_m\) are \(\sigma(m)\) and,
when present, \((m/p)\sigma(m)\). On their orthogonal complement,
\(\|T_m\|\le\tau(m)\sqrt m\), when \(m\) is coprime to the algebra
discriminant and the added level.

These facts use the complete
[classical Hecke transfer lemma](../../results/eisenstein-restricted-prime-sums/sections/hecke-transfer.tex),
with the normalizations proved in Sections 2–3 of the −19 note.
The actual primary inputs are Voight's all-integral-ideal degree and
recurrence, Hasse–Schilling and strong approximation; Jacquet–Langlands,
Theorem 14.4 and Section 5, Proposition 5.20; and Deligne, *Weil I*,
Theorem 8.2 (arbitrary level, outside the level). Their precise retained
versions and theorem locations are in the
[source ledger](primitive_norm_heegner_refined_bounds_source_ledger.json).
Changing 19 to \(D\) changes the mass and point probability in (H.7)–(H.8),
not this transfer argument. In particular no level-specific value of
\(a_2\), theta identity, or cancellation on an auxiliary cusp space is
assumed.

For even \(n=2^t m\), \(m\) odd, \(D\nmid n\), put
\[
 S(n)=3\cdot2^{t-1}\sigma(m)\ge\frac32n,\quad
 \Delta(n)=\tau(m)\sqrt m\,2^{t/2}\frac{t+3}{2},\quad
 \delta(n)=\frac{\Delta(n)}{S(n)}.                         \tag{H.9}
\]
Since 2 is inert, a field norm is even precisely when both its
coordinates are even. The shell with both field norms odd has size
\(G(n)=r_{F_D}(n)-r_{F_D}(n/4)\); nonintegral indices mean zero.
Both components are nonzero. The filtered operator
\(T_n-T_{n/4}\) has character eigenvalues \(S(n)\) and
\((n/p)S(n)\), and complementary norm at most \(\Delta(n)\).
The local factor at 2 is bounded by
\(2^{t/2}(2+(t-1)/2)\), valid for every \(t\ge1\).
The coarse principal return correlation then gives both inequalities
\[
             wS-(1-w)\Delta\ \le\ G(n)/4
                      \ \le\ wS+(1-w)\Delta.             \tag{H.10}
\]
Four generators correspond to each principal ideal. Fine marked-ideal
flag reversal and the inverse-unit coarse measure justify the same
self-adjoint normalization as in the −19 proof.

## 2. The local valuation branches

Let \(\mathcal I\) be the odd inert primes and \(\mu_p=2/(p(p-1))\).
Let \(b_p(n)\) count odd-shell points for which at least one component
norm is divisible by \(p\). The following branches are kept separate:
\[
 \begin{array}{ll}
 p\nmid n:& b_p(n)/4\le w\mu_p S(n)+\Delta(n),\\
 v_p(n)=1:& b_p(n)=0,\\
 v_p(n)\ge2:& b_p(n)=G(n/p^2),\quad
 b_p(n)/4\le wS(n)/p^2+(1-w)\Delta(n)/p.
 \end{array}                                             \tag{H.11}
\]
On the first branch, the bad invertible reductions form the normalizer
of a nonsplit torus, of order \(2(p+1)\). Its mean is \(w\mu_p\)
and its quadratic-norm-character coefficient is zero, because the field
norm onto \(\mathbb F_p^{\times}\) is surjective. The same return
correlation as in the −19 proof gives the first bound.

If \(p\mid n\), anisotropy of \(Q_D\) modulo \(p\) forces both
vectors to be divisible by \(p\) as soon as one component norm is.
This is impossible when \(v_p(n)=1\). Otherwise dividing both vectors
by \(p\) is a bijection onto the odd shell at \(n/p^2\). For
\(a=v_p(n)\ge2\), the exact ratios are
\[
 \frac{S(n/p^2)}{S(n)}=\frac{\sigma(p^{a-2})}{\sigma(p^a)}\le p^{-2},
 \qquad
 \frac{\Delta(n/p^2)}{\Delta(n)}=
                    \frac{a-1}{(a+1)p}\le p^{-1}.
\]
This proves the last line of (H.11). It never invokes a fine-level
operator whose index is divisible by its own level prime.

At the ramified prime the two field-component axes form the principal
fine point. Its two character coefficients are both \(\nu\), whence
\[
 b_D(n)/4\le(1+(n/D))\nu S(n)+(1-2\nu)\Delta(n).
                                                               \tag{H.12}
\]
In fact \(b_D(n)=0\) when \((n/D)=-1\): a nonzero field norm modulo
\(D\) is a square, so a point with one norm zero would force the other
to be \(n\), a square. The certificate conservatively uses the upper
bound \(2\nu S+(1-2\nu)\Delta\) in every target class.

## 3. An area and boundary estimate for odd norms

**Lemma H3 (OddNormLatticeBound).** For every real \(X\ge0\),
\[
 S_{\rm odd,D}(X):=\#\{(x,y)\in\mathbb Z^2:
                0<Q_D(x,y)\le X,\ Q_D(x,y)\text{ odd}\}
 \le \frac{3\pi}{2\sqrt D}X+
                \left(3+\frac6{\sqrt D}\right)\sqrt X+2.
                                                               \tag{H.13}
\]

**Proof.** Complete the square:
\(Q_D=(x+y/2)^2+Dy^2/4\). Set
\(h(y)=\sqrt{\max(X-Dy^2/4,0)}\) and
\(Y=2\sqrt{X/D}\). In each integer row \(|y|\le Y\), the number of
allowed \(x\)'s lies between \(2h(y)-1\) and \(2h(y)+1\).
Monotonicity on the positive half-axis gives
\[
 2\int_0^Yh(y)\,dy-h(0)
 \le\sum_{y\in\mathbb Z}h(y)
 \le2\int_0^Yh(y)\,dy+h(0).
\]
Here the sum can be restricted to \(|y|\le Y\); zero boundary values
cause no problem. The integral is \(\pi X/(2\sqrt D)\), and there
are at most \(4\sqrt{X/D}+1\) rows. Subtracting the origin, the
positive-norm count satisfies
\[
 \frac{2\pi}{\sqrt D}X-\left(2+\frac4{\sqrt D}\right)\sqrt X-2
 \le S_{Q_D}(X)
 \le\frac{2\pi}{\sqrt D}X+\left(2+\frac4{\sqrt D}\right)\sqrt X.
\]
Because \(D\equiv3\pmod8\), the even-norm vectors are exactly the
even vectors. Thus
\(S_{\rm odd,D}(X)=S_{Q_D}(X)-S_{Q_D}(X/4)\). Use the upper bound
at \(X\) and the lower bound at \(X/4\) to obtain (H.13).
The same inequalities include \(X=0\). \(\square\)

Define
\[
 M_D(X)=\max_{1\le a\le X}
       \prod_{\substack{q\mid a\\(q/D)=1}}(v_q(a)+1).
\]
Ideal factorization and the two field units give
\(r_{Q_D}(a)\le2M_D(X)\) for \(1\le a\le X\).
For an inert prime, the bad component is \(p\) times a nonzero
vector of odd norm. Count both component positions and use (H.13):
\[
 b_p(n)\le4M_D(n)\left(\alpha_D\frac n{p^2}
                         +\beta_D\frac{\sqrt n}{p}+2\right),\quad
 \alpha_D=\frac{3\pi}{2\sqrt D},\quad \beta_D=3+\frac6{\sqrt D}.
                                                               \tag{H.14}
\]
Only primes \(p\le\sqrt n\) can contribute.

## 4. A reusable valuation-resolved criterion

Choose a nonempty finite subset \(J\subset\mathcal I\), and a real cutoff
\(z\ge\max J\). For \(p\in J\), write
\(s_p=\min(2,v_p(n))\), and put
\[
 V(s)=\sum_{s_p>0}\mu_p-\sum_{s_p=2}p^{-2},\qquad
 W(s)=\#\{p\in J:s_p>0\}-(1-w)\sum_{s_p=2}p^{-1}.
\]
Also set
\[
 B_0=\sum_{p\in\mathcal I}\mu_p,\quad
 \Pi(z)=\#\{p\in\mathcal I:p\le z\},\quad
 T(z)=\sum_{\substack{p\in\mathcal I\\p>z}}p^{-2},
\]
\[
 U(z,X)=\sum_{\substack{p\in\mathcal I\\z<p\le X}}p^{-1},\qquad
 V_0(z,X)=\#\{p\in\mathcal I:z<p\le X\}.
\]

**Proposition H4 (ValuationResolvedQuaternionCriterion).** Let \(n\ge2\)
be even and \(D\nmid n\). If
\[
 \begin{split}
 \mathcal B_D(n;z,s):={}&\frac4{D+1}+B_0-V(s)\\
 &+\frac{2M_D(n)}{3w}\left(
       \alpha_DT(z)+\frac{\beta_D}{\sqrt n}U(z,\sqrt n)
                    +\frac2nV_0(z,\sqrt n)\right)\\
 &+\frac{\delta(n)}w
                  \bigl(2-w-2\nu+\Pi(z)-W(s)\bigr)<1,
 \end{split}                                             \tag{H.15}
\]
then \(n\in C_D+C_D\).

Subtract (H.11) for small inert primes and (H.12) from (H.10).
Keeping the exact branches on \(J\) saves \(wSV(s)\) in the main
cost and \(\Delta W(s)\) in the error cost. Replace
\(\sum_{p\le z}\mu_p\) by the larger \(B_0\). For the remaining
primes use (H.14), and divide by \(4wS\), with \(S\ge3n/2\).
The ramified main term becomes \(2\nu/w=4/(D+1)\); the two return
errors sum to \(2-w-2\nu\). The surviving count is at least
\(4wS(1-\mathcal B_D)\). Any surviving shell point has two positive
norms supported exclusively on split primes. This proves the criterion.
No probabilistic independence or simultaneous-level spectral estimate
has been assumed.

## 5. Complete Euler and interval bounds

The certificate uses
\[
 J_{43}=\{3,5,7,19,29\},\qquad
 J_{67}=J_{163}=\{3,5,7,11,13\}.
\]
These are inert primes; all \(3^5=243\) valuation patterns are retained.
The following global Euler bounds hold:
\[
\begin{array}{c|ccc}
 D&43&67&163\\\hline
 M_D(X)\le C_M X^{1/6}: C_M&4&3&2\\
 \delta(n)\le C_\delta n^{-1/3}
       \text{ for inert cores}: C_\delta&7&9&14.
\end{array}                                             \tag{H.16}
\]
To prove the first row, maximize \((a+1)^6/p^a\) at each split prime
\(p<64\), for \(0\le a\le10\). Exponents beyond 10 cannot improve
a factor because \((12/11)^6<2\) and the successive ratios decrease.
For \(p\ge64\), \(a+1\le2^a\) bounds the factor by 1. Multiplying
the finite exact maxima proves the stated constants.

For the second row, direct substitution in (H.9) gives
\[
 [\delta(n)n^{1/3}]^6
   =\frac{(t+3)^6}{3^6 2^t}
       \prod_{p^a\parallel m}
           \frac{(a+1)^6p^{5a}}{\sigma(p^a)^6}.             \tag{H.17}
\]
The dyadic factor is maximized at \(t=6\), where it is \((3/2)^6\).
For each odd inert \(p<64\), exponents through 4 suffice: after
\(a=4\), the ratio is less than \((6/5)^6/3<1\).
For \(p\ge64\), every factor is at most 1. The checker multiplies
these exact rational maxima and compares their product with
\(C_\delta^6\).

The interval argument uses more information than (H.16). For an inert
core \(n\le U\), truncate \(t\) to 6. For each odd inert \(p<64\),
let \(c_p\) be the first nonnegative integer with
\[
                         (c_p+2)^6<p(c_p+1)^6.
\]
The factor in (H.17) decreases for every exponent \(a\ge c_p\), since
\(\sigma(p^{a+1})>p\sigma(p^a)\). Truncate to \(c_p\), except that
at a prime in \(J\) retain at least exponent 2, to preserve its state.
Delete the primes at least 64, and also any unbranched prime whose
factors are all at most 1. These operations reduce the integer and
never reduce the product (H.17). They preserve every \(s_p\).

Consequently the following finite enumeration is an upper bound,
including its size constraint: list the truncated integers \(u\le U\),
group by state \(s\), and maximize (H.17) in each group. Let the exact
maximum be \(F_s(U)\). Empty groups cannot arise from a core in the
interval. Put
\[
 c_s(U)=\frac{1+\lfloor10^4 F_s(U)^{1/6}\rfloor}{10^4},\qquad
 d(U)=M_D(U).
\]
For \(L\le n\le U\),
\[
 \delta(n)<\frac{c_s(U)}{\lfloor L^{1/3}\rfloor},\qquad
                         M_D(n)\le d(U).                 \tag{H.18}
\]
All roots are bounded with integer root arithmetic. The coefficient
of \(\delta\) in (H.15) is positive: \(\Pi(z)\ge|J|\),
\(W(s)\le|J|\), and \(2-w-2\nu>0\). Thus replacing it by
an upper bound preserves the inequality. The truncated
enumerations have respectively 6480, 51840, and 207360 entries before
removing dominated entries; no large interval of target integers is
being searched.

The divisor maximum \(d(U)\) is also exact. Remove nonsplit factors
and place decreasing positive exponents on an initial prefix of the
smallest split primes. Exchanging \(a<b\) at \(p<q\) decreases the
integer by \((p/q)^{b-a}\) without changing its divisor product.
This proves that the recursive enumeration over such prefixes includes
a maximizer. It stops at \(2^{96}\); the available prime prefix has
product greater than this bound, checked exactly.

### Prime sums and numerical rounding

The odd inert primes lie in \(r=(D-1)/2\) classes modulo \(2D\).
For every \(z>0\), first-term-plus-integral bounds give
\[
 T(z)\le\frac r{2Dz}+\frac r{z^2},\qquad
 \Pi(z)\le\frac r{2D}z+r,\qquad
 B_0\le\sum_{p\le1024}\frac2{p(p-1)}
           +3\left(\frac r{2D\cdot1024}+\frac r{1024^2}\right).
                                                               \tag{H.19}
\]
All prime sums in this section are over odd inert primes. For
\(X>z\), the same progression bounds are
\[
 U(z,X)\le\frac r{2D}\log(X/z)+\frac rz,\qquad
 V_0(z,X)\le\frac r{2D}(X-z)+r.                            \tag{H.20}
\]
A progression has at most one endpoint contribution; integrating the
decreasing functions \(u^{-2}\) and \(u^{-1}\) over the rest proves
these formulas. For \(z\) above all residue representatives, the
checker also uses the sharper bound
\(\Pi(z)\le rz/(2D)+r-\sum a/(2D)\).

The primes through \(R=2^{20}\) are sieved exactly. Reciprocal and
square-reciprocal sums through \(R\) are bounded term by term by
rounding upward to multiples of \(2^{-80}\). Beyond \(R\), use
(H.19)–(H.20). A logarithm of a ratio in \([1,2^k]\) is bounded by
\(7k/10\). The inequality \(\log2<7/10\) follows already from
\(1+7/10+(7/10)^2/2+(7/10)^3/6>2\).
Use \(\pi<22/7\) and the lower bound
\(\lfloor10^6\sqrt D\rfloor/10^6\) for \(\sqrt D\) in (H.14).
Square-root denominators are rounded downward. Thus every acceptance
comparison is between exact rational numbers.

For each dyadic interval \([2^k,2^{k+1}-1]\), from the exponent in
(H.2) through \(k=95\), the checker substitutes (H.18), the exact
\(d(U)\), and these prime-sum bounds into (H.15), for every possible
state. The cutoff can be any \(z\ge\max J\). The deterministic choice
in the checker balances its large-prime and error terms; the chosen
cutoff and the exact worst cost for every interval are retained in the
receipt. Optimality of the cutoff is not needed. There are 42768
core-state inequalities in total. Each has cost less than \(49/50\).
The largest costs, for information only, are approximately
0.967701, 0.961743, and 0.973195. These decimals are not acceptance
criteria.

## 6. Every exceptional-quotient family and the infinite tail

Let \(e\in E_D\) and write \(e=2^t m_e\), \(m_e\) odd. For any
split prime \(p\), including one dividing \(e\), (H.9) gives
\[
 \delta(ep)\le A_e(ep)^{-1/2},\qquad
 A_e=\frac{2(t+3)}3\frac{\tau(m_e)m_e}{\sigma(m_e)}.        \tag{H.21}
\]
Indeed \(\delta(n)\sqrt n=(t+3)\tau(m)m/(3\sigma(m))\).
If \(a=v_p(m_e)\), adjoining \(p\) changes the odd factor by
\[
 \frac{a+2}{a+1}\frac{p\sigma(p^a)}{\sigma(p^{a+1})}\le2.
\]
No hypothesis \(p\nmid e\) is needed. Since every prime in \(J\) is
inert, the state of \(ep\) on \(J\) is exactly the state of \(e\).
The maxima of \(A_e\) are respectively
\[
                          7,\qquad70/9,\qquad375/31.
\]
Substituting \(A_e/\lfloor\sqrt L\rfloor\) for \(\delta(n)\) in
(H.15) proves all the thin families on the same intervals as above.
There are 5338 such inequalities, each of cost less than \(7/10\).

For completeness, neither interval enumeration is extrapolated beyond
\(2^{96}\). Put \(\gamma=r/(2D)\),
\(c_0=4/(D+1)+B_0\), and let \(a,b\) be the rational upper bounds
for the two coefficients in (H.14). Ignore the savings \(V,W\).
For cores choose \(z=n^{1/4}\), use (H.16), and bound
\(U(z,\sqrt n)\le1+\tfrac12\log n\) and
\(V_0(z,\sqrt n)\le\sqrt n\). Formula (H.15) is at most
\[
 \begin{split}
 c_0+\frac{2C_M}{3w}\bigl[&a\gamma n^{-1/12}
  +(ar+b(1+\tfrac12\log n)+2)n^{-1/3}\bigr]\\
 &+\frac{C_\delta}{w}\bigl[
    \gamma n^{-1/12}+(2-w-2\nu+r)n^{-1/3}\bigr].
 \end{split}                                             \tag{H.22}
\]
For the thin families choose \(z=n^{1/3}\). The corresponding bound is
\[
 \begin{split}
 c_0+\frac{2C_M}{3w}\bigl[&a\gamma n^{-1/6}+ar n^{-1/2}
  +(b(1+\tfrac12\log n)+2)n^{-1/3}\bigr]\\
 &+\frac{\max_e A_e}{w}\bigl[
    \gamma n^{-1/6}+(2-w-2\nu+r)n^{-1/2}\bigr].
 \end{split}                                             \tag{H.23}
\]
All displayed terms decrease for \(n\ge2^{96}\); for the logarithmic
ones use \((\log n)n^{-\beta}\) decreasing once \(\log n>1/\beta\).
At \(2^{96}\), the relevant powers are exactly \(2^{-8},2^{-32}\)
and \(2^{-16},2^{-48},2^{-32}\). Bounding \(\log(2^{96})<96\cdot7/10\),
the checker verifies that (H.22) is less than \(19/20\) and (H.23)
less than \(3/5\) in all three cases. This proves the full infinite
tails, and finishes Theorem H1.

## 7. Why these are complete finite reductions

The finite enumeration through 1024 proves exactly that all members of
\(E_D\) are absent from \(C_D+C_D\). Assume (H.3), and suppose
\(n\) is the least additional even \(D\)-unit exception. Necessarily
\(n>X_D\). If \(n\) has a split prime divisor \(p\), then
\(n/p\) is a smaller even \(D\)-unit. Unless \(n/p\in E_D\), its
representation in \(C_D+C_D\), multiplied by \(p\), represents
\(n\). If \(n/p=e\in E_D\), the second part of Theorem H1 applies
to \(n=ep\). If there is no split divisor, all prime divisors of
\(n\) are inert, and the first part of H1 applies. These cases
exhaust every possible least exception. This proves (H.4); its converse
implication to (H.3) is immediate.

For \(D\)-unit targets, (H.4) supplies primitive representations
outside \(E_D\); the finite calculation of \(R_D\) decides precisely
which members of \(E_D\) are filled by mixed \(C_D+D C_D\) sums.
For targets divisible by \(D\), (H.6) transports the exceptions to
\(D E_D\). Odd targets and multiples of \(D^2\) are impossible.
Thus (H.5) follows. Conversely apply (H.5) to \(Dm\) and use (H.6).
This establishes Theorem H2 without presuming that an exceptional
quotient is coprime to the removed split prime.

The decisive remaining task is the complete finite coverage in (H.3),
not another small list of plausible exceptions. The present checker
has no mode that reports a classification after merely passing its
analytic tests.

## 8. Finite inputs and evidence

The finite sets in (H.1), computed independently by split-factor
membership and checked against primitive coordinate enumeration, are
\[
 E_{43}=\{4,6,8,10,16,20,38,50,56,74\},
\]
\[
 E_{67}=\{4,6,8,10,12,14,16,22,26,28,32,44,50,62,68,80,
                  86,98,114,116,138\},
\]
\[
\begin{split}
 E_{163}=\{&4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,\\
 &46,50,52,56,58,60,64,66,68,70,74,76,78,80,92,110,116,120,128,\\
 &134,146,148,162,164,170,176,182,186,188,190,206,218,230,236,254,\\
 &272,278,290,300,332,338,344,362,386,392,484,568,992\}.
\end{split}
\]
The primitive candidate lists, all \(A_e\), every dyadic interval,
Euler maxima, finite-domain counts and exact rational costs are retained
in the [receipt](primitive_norm_heegner_refined_bounds_checks.json).
The deterministic checker reconstructs every input; its trace digest
is supplementary, not a substitute for the inequalities.

Run from the repository root:

```sh
python3 research/arithmetic/verify_primitive_norm_heegner_refined_bounds.py \
  --output research/arithmetic/primitive_norm_heegner_refined_bounds_checks.json
```

It uses only Python's standard library and rejects optimized execution.
The written argument proves the infinite reductions and imports the
classical spectral estimates. The checker proves the stated finite
rational certificates and small exact controls. It does not prove the
deep imported theorems, perform the outstanding large coverage, or
establish formal verification, independent external review, historical
priority, or publication.
