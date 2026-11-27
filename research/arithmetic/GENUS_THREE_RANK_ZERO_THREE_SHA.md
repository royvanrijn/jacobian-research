# Rank-zero genus-three Jacobians with unbounded elementary three-Sha

This note gives an explicit HLP gluing family, an unbounded three-isogeny
Selmer ratio, and a rank-zero quadratic-twist selection. The rank selection
uses Smith's 2025 preprint, Theorem 1.1 and Corollary 1.3, as an external
theorem. The Jacobians are split. Their entire three-primary Sha and its
finiteness are not determined. Historical priority and independent
mathematical review remain unestablished.

The general rank-controlled large-Sha construction objective is achieved.
Controlling the genus-two companion rank in the original Mordell cover
[HG8](GENUS_THREE_MARKED_TORSION_DEFORMATIONS.md) remains a separate open
task; the present elliptic factors are different.

## Precise primary inputs

* [Howe–Leprévost–Poonen](sources/howe-leprevost-poonen-split-jacobians.pdf),
  Lemma 13, Proposition 14 and Section 4.2 equation (10): non-split
  (2,2,2) gluing and its hyperelliptic double cover of a conic.
* [Bhargava–Klagsbrun–Lemke Oliver–Shnidman](sources/bhargava-klagsbrun-lemke-oliver-shnidman-three-isogeny.pdf),
  Theorem 6.2, Section 9.1 equation (9.2), Propositions 10.1–10.4 and
  Theorem 10.5: Cassels' isogeny-Selmer ratio, local factors, and
  invariance at good primes under quadratic twisting. No Selmer
  average theorem is required.
* [Williams](sources/williams-mertens-ap-1974.pdf), Theorem 1, printed
  p. 356: Mertens' product for a fixed reduced residue class.
* [Delaunay–Roblot](sources/delaunay-roblot-quadratic-twist-regulators.pdf),
  Section 1, printed pp. 603–604: for a fundamental discriminant $d$
  coprime to the conductor $N$, $w(E^{(d)})=w(E)\chi_d(-N)$.
* [Smith, arXiv:2503.17619v1](sources/smith-bsd-implies-goldfeld-2503.17619v1.pdf),
  Theorem 1.1, inequalities (1.1)–(1.2), and Corollary 1.3, PDF p. 2:
  root-number +1 twists have rank zero with density one. This is an
  unconditional Selmer-corank consequence; the BSD-conditional
  analytic-rank corollary is not used.

These actual statements were inspected. The retained Smith file is the
March 25, 2025 PDF with arXiv v1 identifier, rather than the August 2026
date shown by a later HTML rendering. No complete independent replay of
that external proof is claimed.

## Explicit gluing family and growing Selmer ratio

For an even integer $s\ge0$, put
\[
 h=-\frac{3(s^2+31)}{2(s^2+49)},\quad b=2h+3=\frac{54}{s^2+49},\quad
 a=h^2+12h+9,\quad \Delta=a^2-4b^3=(h+1)(h-3)^3.
\]
Define elliptic curves
\[
 \begin{split}
 E_h&:y^2=x^3+(hx+h+1)^2,\\
 E'_h&:y^2=U(U^2+aU+b^3),\\
 F_h&:y^2=x(x-1)(x-\lambda),\qquad\lambda=\frac{a^2}{4b^3}.
 \end{split}
\]
The normalized three-isogeny $\phi:E_h\to E'_h$ has kernel generated
by $(0,h+1)$ and explicit formulas
\[
 U=x+\frac{4h(h+1)}x+\frac{4(h+1)^2}{x^2}-(4h+3),\qquad
 y'=y\left(1-\frac{4h(h+1)}{x^2}-\frac{8(h+1)^2}{x^3}\right).
\]
Let $\psi:E'_h\to E_h$ be its dual, whose kernel is $\mu_3$.

For the hyperelliptic equation put
\[
 u=-\frac{7h+5}{2},\qquad v=\frac{(h+1)s}{2},
\]
\[
 L_s(t)=u(1-t^2)-2vt,\qquad M_s(t)=v(1-t^2)+2ut,
\]
\[
 G_s(t)=b\left(16b^2L_s(t)^2M_s(t)^2-\Delta(1+t^2)^4\right),
 \qquad C_s:y^2=G_s(t). \tag{14}
\]
Write $\omega_{\ge5}(m)$ for the number of distinct prime divisors at
least five of a nonzero integer $m$, using its absolute value, and set
\[
 W(s)=\omega_{\ge5}(s^2-5)+\omega_{\ge5}(s^2+49)
                            -\omega_{\ge5}(s^2+43).
\]

**Proposition HG13.** The curves (14) are smooth of genus three and their
Jacobians have degree-eight dual isogenies with
\[
                         E'_h\times E'_h\times F_h,
\]
whose composites are [2]. The global Selmer ratio of $\psi$ satisfies
\[
                          c(\psi)\ge3^{W(s)-2}. \tag{15}
\]
Furthermore $W(s)/\log\log T\to1$ in probability as even integers
$s\in[0,T]$ are ordered by size. In particular $c(\psi)$ is unbounded
in this rational gluing family.

**Geometry.** Substitution gives
\[
 u^2+v^2=-\frac a{2b},\quad
 h+1=\frac{5-s^2}{2(s^2+49)},\quad
 h-3=-\frac{9(s^2+43)}{2(s^2+49)},\quad
 a=-\frac{27(s^4+62s^2+529)}{4(s^2+49)^2}.
\]
Thus $a,b,\Delta$ are nonzero and all elliptic factors are smooth for
integer $s$. In HLP Proposition 14 take $A_1=A_2=a$, $B_1=B_2=b^3$,
$A_3=-1-\lambda$, $B_3=\lambda$. The discriminant product root is
$R=\Delta(1-\lambda)$ and the twisting factor is zero. Their equations are
\[
 2b^3(X^2+Y^2)+aZ^2=0,
\]
\[
 W_0^2Z^2=-\frac{a^2b^3}{2}(X^4+Y^4)
                         +\frac{a^2(a^2+4b^3)}{16b^3}Z^4.
\]
The conic has point $(u:v:b)$ and parametrization
$(X,Y,Z)=(L_s,M_s,b(1+t^2))$. Substitution into the quartic gives
$a^2G_s(t)/16$; the ordinate change $y=4W_0Z/a$ gives (14).
HLP supplies the polarized gluing and smoothness. The leading coefficient
can vanish only if $b=(u^2-v^2)^2$, whereas $v_3(b)=3$. Thus (14) has
degree eight for every integer $s$.

**Local ratios.** For $p\ge5$, the three integers
$s^2-5,s^2+49,s^2+43$ are pairwise coprime outside 2 and 3, and the table is

| Condition on $p\ge5$ | $c_p(\psi)$ |
|---|---:|
| $p\mid s^2-5$ | 3 |
| $p\mid s^2+49$ | 3 |
| $p\mid s^2+43$, with $(-3/p)=1$ | 1/3 |
| $p\mid s^2+43$, with $(-3/p)=-1$ | 1 |
| None of these | 1 |

For a direct reduction check, write $h=m/n$ and use marked integral models
\[
 \begin{array}{c|cc}
 &A&B\\ \hline
 E_h&m^2-3n^2&n^3(2m+3n)\\
 E'_h&m^2+12mn+9n^2&n(2m+3n)^3.
 \end{array}
\]
Their quadratic discriminants are $(m+n)^3(m-3n)$ and
$(m+n)(m-3n)^3$. Our parametrization has $2m+3n=108$, up to factors
supported at 2 and 3. At $m+n=0$ the multiplicative discriminant orders
have ratio 3:1, with square tangent slope. At $n=0$ they have ratio 6:2,
again with square slope. At $m-3n=0$ they have ratio 1:3, with tangent
squareclass $-3$. The models have unit $c_4$ at each such prime.
Split multiplicative Tamagawa ratios give the table; nonsplit
multiplicative reduction gives ratio one. In particular this computation
fixes the isogeny direction and places the large Sha on the glued factor.

BKLOS Proposition 10.4 and Theorem 10.5 bound each local factor at 2 and
3 below by $1/3$. The real local factor is one because $\mu_3(\mathbb R)$
is trivial. Multiplying proves (15). Counting every prime divisor of
$s^2+43$ negatively only weakens the bound; it is not an equality.

**Unboundedness.** Each of $s^2-5,s^2+49,s^2+43$ is irreducible over
$\mathbb Q$. Away from finitely many primes, its number of roots modulo
$p$ is $\rho(p)=1+\chi(p)$ for a fixed nonprincipal quadratic character.
Williams' Theorem 1, summed over fixed residue classes and logarithmically
expanded, gives
\[
                    \sum_{p\le z}\frac{\rho(p)}p=\log\log z+O(1).
\]
Restricting $s$ to even integers does not change the estimate. Set
$z=T^{1/8}$. The truncated sum of indicators $p\mid f(s)$ has this mean
and variance $O(\log\log T)+O(z^2/T)$. Indeed, for distinct primes the
Chinese remainder theorem gives joint density $\rho(p)\rho(q)/(pq)$,
with error $O(\rho(p)\rho(q)/T)$. There are only a bounded number of
distinct prime divisors greater than $z$, since $|f(s)|\ll T^2$.
Chebyshev's inequality gives
$\omega(f(s))/\log\log T\to1$ in probability. Intersecting the three
density-one sets proves the asserted normal order of $W$, without an
independence assumption between the three polynomials. This proves HG13.

## Rank zero and arbitrarily many elementary Sha directions

**Theorem HG14.** For every integer $r\ge1$, the explicit family
\[
                         C_{s,d}:d y^2=G_s(t)
\]
contains infinitely many rational isomorphism classes with Jacobians
satisfying
\[
 \operatorname{rank}J_{s,d}(\mathbb Q)=0,\qquad
                    (\mathbb Z/3)^{2r}\subseteq\Sha(J_{s,d})[3].
\]
For a suitable fixed even $s$, eligible signed squarefree twist parameters
$d$ have positive density. The theorem uses the Smith input above and
requires neither BSD nor finiteness of Sha.

**Proof.** Choose even $s$ with $W(s)\ge r+3$, using HG13; then
$c(\psi)\ge3^{r+1}$. Fix this $s$ and abbreviate $E'_h,F_h$ by $E',F$.

To arrange compatible even root numbers, set $P(s)=s^4+62s^2+529$.
For even $s$ it is odd, greater than one, and prime to three. Choose a
prime $q\ge5$ dividing it. Exact congruences give
\[
 P(s)\equiv864\pmod{s^2-5},\quad
 P(s)\equiv-108\pmod{s^2+49},\quad
 P(s)\equiv-288\pmod{s^2+43}.
\]
Thus $E_h,E'$ are good at $q$. At that prime $a$ vanishes and $b$ is a
unit, so $\lambda$ has positive valuation. The Legendre curve $F$ is
multiplicative and has conductor exponent one at $q$.

Let $S$ contain 2,3 and all bad primes of $E_h,E',F$. Put
$\tau=w(E')$, $\eta=w(F)$, and impose these conditions on squarefree $d$:

* $d\equiv1\pmod8$, and $d$ is a local square at every odd prime of
  $S\setminus\{q\}$;
* $d$ is a unit at $q$ with $(d/q)=\tau\eta$;
* $\operatorname{sgn}(d)=\tau$.

These conditions define a nonempty reduced arithmetic progression with a
sign, so squarefree members have positive density. They are fundamental
discriminants prime to both conductors. The coprime-twist formula gives
\[
 w(E'^{(d)})=\tau\operatorname{sgn}(d)=1,\qquad
 w(F^{(d)})=\eta\operatorname{sgn}(d)(d/q)=1.
\]
The prime where only the companion is bad permits independent control
of the two signs.

At every finite bad prime of $E_h,E'$, the isogeny local squareclass is
unchanged. At all other finite primes the ratio is one under every twist
by BKLOS Theorem 6.2, with 2 and 3 already fixed. Only the real factor can
change: it is one for positive $d$ and $1/3$ for negative $d$. Hence
\[
                           c(\psi_d)\ge c(\psi)/3\ge3^r.
\]
Smith's Corollary 1.3 makes the rank-positive exceptions among root-number
+1 twists of either fixed curve a set of size $o(H)$. The union is still
$o(H)$. Our local conditions have positive density among signed squarefree
integers, so a density-one subset of them has
\[
 \operatorname{rank}E'^{(d)}(\mathbb Q)
          =\operatorname{rank}F^{(d)}(\mathbb Q)=0.
\]
The glued Jacobian has rank zero. Here $s$ is fixed before the twist
density is taken; no uniformity in Smith's theorem as $s$ varies is needed.

Discard the finitely many squareclasses $d=1,-3$. The two isogeny kernel
characters are those of $d$ and $-3d$, so their rational kernels vanish.
Cassels' formula gives
\[
 |\operatorname{Sel}_{\psi_d}(E'^{(d)})|
   =c(\psi_d)|\operatorname{Sel}_{\phi_d}(E_h^{(d)})|\ge3^r.
\]
The same two characters are the composition factors of $E_h^{(d)}[3]$,
so it has no rational fixed vector. The isogenous curve $E_h^{(d)}$ has
rank zero and rational torsion prime to three. Therefore
$E_h^{(d)}(\mathbb Q)/\psi_d E'^{(d)}(\mathbb Q)$, killed by three, is zero.
The Kummer sequence identifies the displayed Selmer group with
$\Sha(E'^{(d)})[\psi_d]\subseteq\Sha(E'^{(d)})[3]$.

The degree-eight gluing and its dual have composites [2]. Their induced
maps on Sha have the same composites, and multiplication by two is an
automorphism of three-primary torsion. Thus they identify full
three-primary Sha with the direct sum for the three elliptic factors,
without finiteness. They commute with [-1], so the argument applies to
each quadratic twist, as in
[HG5](GENUS_THREE_THREE_MARKED_ELLIPTIC_QUOTIENTS.md).
The two copies of $E'^{(d)}$ give at least $2r$ independent directions.
There are infinitely many rational forms: a finite automorphism group
identifies at most finitely many quadratic twist classes. This proves HG14.

## Exact fixture and verification scope

The pair $s=8,d=57$ gives rank zero and
$(\mathbb Z/3)^4\subseteq\Sha(J_{8,57})[3]$ without invoking Smith for
this individual example. The distinct elliptic factors have minimal
$a$-invariants
\[
 E'^{(57)}:[1,1,0,-546524405,-4917455709459],
\]
\[
 F^{(57)}:[0,-1,0,-32816655601327594443617,
                         -2284648019026739769616608032943135].
\]
Exact mwrank descent, with databases disabled, gives rank and upper bound
zero for both. The rational three-isogeny kernels vanish and
$c(\psi_{57})=9$: the finite local factors are 3 at 59 and 113 and one
everywhere else, including 2,3,19,107. The real factor is one.
Cassels' formula and the zero Mordell–Weil quotient supply two directions
on each copy of $E'^{(57)}$, hence four on the Jacobian. No complete Sha
order or equality of primary groups is inferred.

The retained command is

    sage -python research/arithmetic/verify_genus_three_rank_zero_three_sha.py

It passed with SageMath 10.9. It verifies the generic isogeny and gluing
identities, rational conic and parametrization, finite local controls,
companion-prime residue controls, the exact rank-zero fixture and Selmer
ratio, and finite-field Jacobian counts against the elliptic product.
The normal-order proof, density selection and imported theorems are
written inputs; they do not follow from finite tests. Independent replay,
external review and formal verification are absent. The geometrically
simple torsion, maximal three-rank family and three virtual-direction
goals remain open.
