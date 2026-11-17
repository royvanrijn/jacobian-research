# Distinct Hasse-failing fibre algebras on the fixed target line

The distinct Hasse-failing fibre algebras on the fixed target line satisfy
`D_Phi(B) ~ C_* B^2/log B` in the actual target height. The positive
constant is a convergent limit of explicit finite rectangle-area sums.
The proof determines the minimum height in each algebra class and controls
the large-prime allocation tail. A squarefree restricted family also has
an explicit Euler-product constant.

Retain `Phi` and `y_a` from the [fixed-map manuscript](../../unpublished/fixed-map-hasse-failures/sections/02-fixed-map.tex).
For a noncube `a`, the fibre algebra is

\[
 \mathbb Q(\sqrt[3]{a})\times\mathbb Q(\zeta_3).
\]

The pure-cubic isomorphism classification in
[the local-density note](FIXED_MAP_LOCAL_DENSITIES.md) identifies the
algebras at noncube parameters `a,b` precisely when their nonzero rational
cube classes are equal or inverse. The sign is itself a rational cube.

## 1. A subfamily with exact multiplicity two

Let `S` be the positive squarefree integers whose prime divisors are all
`1 (mod 3)`; include `1`. Consider pairs

\[
 m,n\in S,\qquad (m,n)=1,\qquad m\equiv n\pmod9,
 \qquad (m,n)\ne(1,1).                               \tag{1.1}
\]

Both signs of `a=m/n` satisfy the rational local criterion. Neither is a
cube: the nonzero valuations of `m/n` are `1` or `-1`, and at least one
is nonzero. Thus these pairs give Hasse failures.

**Lemma 1.1 (unique representatives up to inversion).** Within (1.1), two
pairs give isomorphic fibre algebras if and only if they are identical or
obtained from one another by interchanging numerator and denominator.

**Proof.** At every prime, the valuation of a squarefree reduced ratio is
`0,1,-1`. These three values have distinct residues modulo three. Equal
cube classes therefore give identical valuations and identical positive
reduced ratios; inverse classes interchange the two coordinates. Conversely,
interchanging coordinates inverts the cube class and preserves the cubic
field. The pure-cubic isomorphism classification supplies necessity for
field isomorphism, and the quadratic field component has a different degree
and cannot be exchanged with the cubic field. \(\square\)

Every integer in `S` is coprime to six. Hence the content of the target's
projective coordinates at `a=-m/n` is one. Moreover,

\[
 H(y_{-m/n})=\max(32m,9n),\qquad
 H(y_{m/n})\ge H(y_{-m/n}).                           \tag{1.2}
\]

Indeed `|3n-24m|<=max(9n,32m)`. Thus allowing either sign never improves
on the negative target. Let

\[
 R=(0,1/32]\times(0,1/9],\qquad U=R\cup R^{\mathsf T}.
                                                               \tag{1.3}
\]

A class represented by (1.1) has a target from this restricted subfamily
of height at most `B` precisely
when one, equivalently both, of its two ordered representatives belongs
to `BU`. Their area is

\[
 \operatorname{area}(U)
      =\frac2{32\cdot9}-\frac1{32^2}=\frac{55}{9216}. \tag{1.4}
\]

## 2. Coprime squarefree split pairs

Let `f` be the indicator of `S`. Write `G3` for the analytic factor defined
in [the integer-count note, equation (2.1)](FIXED_MAP_INTEGER_COUNT.md), so

\[
 \prod_{p\equiv1(3)}(1-p^{-s})^{-1}=\zeta(s)^{1/2}G_3(s).
\]

For positive fixed `alpha,beta`, let `P_{alpha,beta}(B)` count the pairs
`m<=alpha B,n<=beta B` in `S` with `(m,n)=1` and `m==n (mod 9)`, including
`(1,1)` when it lies in the rectangle. Define the positive convergent product

\[
 \mathfrak P=\prod_{p\equiv1(3)}(1-3p^{-2}+2p^{-3}). \tag{2.1}
\]

**Theorem 2.1 (rectangle pair density).**

\[
 P_{\alpha,\beta}(B)
   \sim\frac{G_3(1)^2\mathfrak P}{3\pi}
                 \alpha\beta\frac{B^2}{\log B}.       \tag{2.2}
\]

### Proof of Theorem 2.1: squarefree sums and character means

Let `chi` have order three modulo nine. Its kernel on units is `{1,-1}`.
Since all elements of `S` lie in `{1,4,7}` modulo nine, the projector

\[
 \mathbf1_{m\equiv n(9)}
       =\frac13\sum_{j=0}^2\chi(m)^j\chi(n)^{-j}      \tag{2.3}
\]

is valid on `S x S`. Set `f_j(n)=f(n)chi(n)^j` and
`A_j(X)=sum_{n<=X}f_j(n)`. Their prime-power values have absolute value
at most one. Their absolute prime density is `1/2`, and their prime means
are `1/2,0,0`. Therefore the same
[de la Breteche--Tenenbaum Theorem 1.2, p. 3](https://tenenb.perso.math.cnrs.fr/PPP/On-SD.pdf)
used in the manuscript applies with `r=1/2,sigma=3/4` and fixed-modulus
prime-distribution errors of every logarithmic order. In particular,

\[
 A_0(X)\sim K_{\rm sf}\frac{X}{\sqrt{\log X}},\quad
 A_j(X)=o\!\left(\frac{X}{\sqrt{\log X}}\right)\ (j=1,2),\quad
 A_0(X)\ll\frac{X}{\sqrt{\log(2X)}},                 \tag{2.4}
\]

where

\[
 K_{\rm sf}=\frac{G_3(1)}{\sqrt\pi}
                     \prod_{p\equiv1(3)}(1-p^{-2}).   \tag{2.5}
\]

The leading factor follows from
`prod_split(1+p^{-s})=zeta(s)^{1/2}G3(s)prod_split(1-p^{-2s})`.
Both nontrivial twists have all main coefficients zero by the reciprocal
gamma factor in the cited theorem.

### Proof of Theorem 2.1: the bivariate coprimality correction

For a split prime put `z=p^{-s}`, `w=p^{-t}`, `eta=chi(p)^j`. The local
factor for coprime squarefree pairs is `1+eta z+eta^{-1}w`. Its correction
to the product of the two univariate factors is

\[
 \frac{1+\eta z+\eta^{-1}w}{(1+\eta z)(1+\eta^{-1}w)}
   =1-\frac{zw}{(1+\eta z)(1+\eta^{-1}w)}.           \tag{2.6}
\]

Let `h_j(d,e)` be the bivariate multiplicative coefficients of the
product of (2.6) over split primes. At a single prime the coefficient at
`(p^u,p^v)`, for `u,v>=1`, is
`(-1)^{u+v-1}eta^{u-v}`; the constant coefficient is one, and all other
nonconstant coefficients with either exponent zero vanish. For every
`sigma>1/2`,

\[
 \sum_{d,e\ge1}\frac{|h_j(d,e)|}{(de)^\sigma}
 =\prod_{p\equiv1(3)}
       \left(1+\frac{p^{-2\sigma}}{(1-p^{-\sigma})^2}\right)
 <\infty.                                           \tag{2.7}
\]

Coefficient comparison gives the exact twisted rectangle sum

\[
 \sum_{d\le\alpha B,\ e\le\beta B}h_j(d,e)
           A_j(\alpha B/d)A_{-j}(\beta B/e).          \tag{2.8}
\]

For `d,e<=sqrt(B)`, (2.4) and `|A_j|<=A_0` bound its summands by
`O_{alpha,beta}(B^2|h_j(d,e)|/(de log B))`. This is a summable dominating
sequence after normalization by `B^2/log B`. For the remaining terms,
use `|A_j(X)|<=X` and (2.7) with `sigma=3/4` to obtain

\[
 B^2\sum_{\max(d,e)>\sqrt B}\frac{|h_j(d,e)|}{de}
     \ll B^{15/8}=o(B^2/\log B).                     \tag{2.9}
\]

Indeed one of `d^{-1/4},e^{-1/4}` is at most `B^{-1/8}`. Dominated
convergence now gives the principal density `K_sf^2 E_sf alpha beta`,
while both nontrivial twisted densities vanish, where

\[
 E_{\rm sf}=\sum_{d,e}\frac{h_0(d,e)}{de}
   =\prod_{p\equiv1(3)}\frac{1+2/p}{(1+1/p)^2}.
\]

Finally,

\[
 K_{\rm sf}^2E_{\rm sf}
  =\frac{G_3(1)^2}{\pi}
      \prod_{p\equiv1(3)}(1-p^{-1})^2(1+2p^{-1})
  =\frac{G_3(1)^2}{\pi}\mathfrak P.                 \tag{2.10}
\]

The projector (2.3) supplies the factor `1/3`, proving (2.2).

## 3. Distinct restricted algebras and the sharp total order

Let `D_sf(B)` count isomorphism classes of fibre algebras represented by
(1.1), with either sign and either ordering giving an actual target of
height at most `B`. Let `D_Phi(B)` count all distinct Hasse-failing fibre
algebras on the original target line at that height.

**Theorem 3.1 (exact restricted-family asymptotic).**

\[
\boxed{
 D_{\rm sf}(B)\sim C_{\rm sf}\frac{B^2}{\log B},\qquad
 C_{\rm sf}=\frac{55G_3(1)^2}{55296\pi}
                    \prod_{p\equiv1(3)}(1-3p^{-2}+2p^{-3})>0.}
                                                               \tag{3.1}
\]

**Proof.** Count the union `BU` by two applications of Theorem 2.1 to
`R,R^T` and subtract their square intersection. Its area is (1.4).
Lemma 1.1 gives exactly two ordered representatives for each nontrivial
class in that union. The sole fixed pair `(1,1)` is excluded. For `B>=32`
the exact identity is therefore

\[
 D_{\rm sf}(B)
  =P_{1/32,1/9}(B)-\frac{P_{1/32,1/32}(B)+1}{2}.     \tag{3.2}
\]

Insert (2.2). The factor `1/2` times the area `55/9216`, together with
the residue factor `1/3`, gives (3.1). \(\square\)

**Corollary 3.2 (sharp order for all distinct algebras on the line).**

\[
 C_{\rm sf}\le
 \liminf_{B\to\infty}\frac{D_\Phi(B)\log B}{B^2}
 \le\limsup_{B\to\infty}\frac{D_\Phi(B)\log B}{B^2}
 \le C_\Phi,
 \qquad D_\Phi(B)\asymp B^2/\log B.                 \tag{3.3}
\]

Here `C_Phi` is the manuscript's complete rational target-line constant.
The lower bound follows from `D_sf<=D_Phi`; the upper bound follows from
counting targets, since each algebra counted has a target representative.
This proves the order of growth for distinct algebras on this one line.
Every admissible class does have a squarefree split representative. The
restricted-height condition above is stronger than admitting any target
representative at that height. Section 4 counts all representatives.

## 4. Minimum heights and the full distinct-algebra asymptotic

Write

\[
 \kappa_{\rm sf}=\frac{G_3(1)^2}{\pi}
       \prod_{p\equiv1(3)}(1-3p^{-2}+2p^{-3}).                    \tag{4.1}
\]

This is the positive principal pair constant before the modulo-nine projector.

### Theorem 4.1: canonical classes and exact minimum target height

An oriented locally admissible nontrivial cube class is uniquely represented
by a positive ratio `m/n`, where `m,n` are coprime squarefree split integers
and `m==n (mod 9)`. At a split prime, class exponent one places the prime
in `m`, and exponent two places it in `n`. An algebra is the pair of inverse
classes, so it identifies `(m,n)` with `(n,m)`.

For a class `v`, let `A(v)` be the finite set of coprime positive cubefree
split pairs `(R,S)` representing it. At a class exponent one, the choices
are `(p,1)` or `(1,p^2)`; at exponent two they are `(p^2,1)` or `(1,p)`.
Absent primes occur in neither entry. Define

\[
 h(R,S)=\min\{\max(4R,9S),\ \max(36R,3S)\}.                     \tag{4.2}
\]

The least actual target height at which the algebra of `v` occurs is exactly

\[
 H_{\min}(v)=\min_{(R,S)\in A(v)}\min\{h(R,S),h(S,R)\}.           \tag{4.3}
\]

**Proof.** The class conditions follow by reducing each valuation modulo three.
All inert-prime and 3-adic valuations vanish in the cube-class group. The
remaining exponents one and two give the indicated squarefree ratio. Both
numerator and denominator lie in `{1,4,7}` modulo nine, so the local cube
condition at three is precisely `m==n (mod 9)`. The field-isomorphism theorem
identifies exactly the two inverse nonzero classes.

A sign is a rational cube. For positive `M,N`,

\[
 H(y_{-M/N})=\frac{\max(32M,9N)}{\gcd(N,8)\gcd(M,3)}
             \le H(y_{M/N})                                  \tag{4.4}
\]

when the fraction is reduced. The absolute fourth coordinate for the negative
parameter is bounded by the displayed maximum. Thus a minimum may be sought
among negative parameters.

At primes other than two and three, removing an available numerator or
denominator cube can only decrease (4.4). Every inert-prime exponent then
vanishes, and split exponents are one or two on one side. At two, a numerator
cube has no beneficial content factor and is removed; a denominator exponent
larger than three can be reduced while retaining the content eight. Between
exponents zero and three, adding a denominator eight keeps the normalized
`9N` term unchanged and divides the normalized `32M` term by eight. Hence
there is always a minimum with denominator exponent three. At three,
denominator cubes can be removed, and numerator exponents larger than three
can be reduced while retaining the content three. Only numerator exponents
zero and three remain.

For each resulting split pair `(R,S)`, the two parameters are therefore
`-R/(8S)` and `-27R/(8S)`. Their actual heights are `max(4R,9S)` and
`max(36R,3S)`. Every reduction used above preserves the cube class and does
not increase height; each displayed candidate is an actual target on the
line. Taking the inverse class gives (4.3). \(\square\)

For example the squarefree representatives of the class of 19 have negative
target heights 608 and 171. The representative `-27/152` has height 57.
Squarefree existence therefore does not identify the least-height target.

Put

\[
 D=((0,1/4]\times(0,1/9])\ \cup\ ((0,1/36]\times(0,1/3]),
 \qquad W=D\cup D^{\mathsf T}.                                  \tag{4.5}
\]

Thus (4.3) is equivalent to some allocation `(R,S)` lying in `BW`.
The four rectangles form a staircase with area

\[
 \operatorname{area}(W)
 =\frac1{108}+\frac1{48}+\frac5{324}+\frac1{432}
 =\frac{31}{648}.                                               \tag{4.6}
\]

### Proposition 4.2: finite-prime allocation constants

Let `P` be a finite set of split primes. For a pattern
`e in {0,1,2}^P`, let `C(e)` consist of all coefficient pairs `(r,s)`
obtained by multiplying the local choices

| `e_p` | Local coefficient pairs |
|---|---|
| 0 | `(1,1)` |
| 1 | `(p,1)`, `(1,p^2)` |
| 2 | `(p^2,1)`, `(1,p)` |

Define the finite union of rectangles

\[
 V_e=\bigcup_{(r,s)\in C(e)}
          \{(x,y):(rx,sy)\in W\}.                              \tag{4.7}
\]

Let `D_P(B)` count algebra classes with a target of height at most `B`
using an allocation with no square prime factor outside `P`, allowing either
inverse class. Then

\[
 D_P(B)\sim C(P)\frac{B^2}{\log B},\qquad
 C(P)=\frac{\kappa_{\rm sf}}
              {6\prod_{p\in P}(1+2/p)}
             \sum_{e\in\{0,1,2\}^P}\operatorname{area}(V_e).     \tag{4.8}
\]

These count distinct algebras, not parameter representatives. If `P` is empty,

\[
 C(\varnothing)=\frac{31\kappa_{\rm sf}}{3888}.                  \tag{4.9}
\]

**Proof.** For a fixed pattern, remove its primes from the canonical squarefree
numerator and denominator. The remaining pair `(x,y)` is coprime, squarefree,
split-supported, and prime to `P`. Write

\[
 m_e=\prod_{e_p=1}p,\qquad n_e=\prod_{e_p=2}p.
\]

The modulo-nine condition is `m_e x==n_e y (mod 9)`. Its order-three character
projector still has principal coefficient `1/3`; both nontrivial twisted
means vanish as in Theorem 2.1. Excluding a finite prime `p` changes the
principal pair constant by `1/(1+2/p)`: the coprime squarefree local pair
factor at `s=t=1` is `1+2/p`. Equivalently remove its univariate factors
and its coprimality correction in (2.5)--(2.10). The same absolute convergence
and tail argument applies after excluding finitely many primes.

Consequently each fixed rectangle union `BV_e` has density
`kappa_sf area(V_e)/(3 product_P(1+2/p))`, by finite inclusion-exclusion
of the rectangle theorem. The canonical oriented class is unique, so taking
a union in (4.7) avoids all allocation multiplicities. Pattern negation and
interchanging `x,y` give exactly the two oriented classes of every nontrivial
algebra. Division by two supplies the remaining factor. The zero class
contributes at most one oriented pair and is discarded. This proves (4.8).
Equation (4.9) follows from (4.6). \(\square\)

### Lemma 4.3: a uniform tail for square allocations

Let `P_z={p<=z:p==1 (mod 3)}`. Uniformly for fixed `z>=2`,

\[
 0\le D_\Phi(B)-D_{P_z}(B)
 \ll \frac{B^2}{\log B}\sum_{\substack{p>z\\p\equiv1(3)}}p^{-2}
           +\frac{B^{7/4}}{\sqrt{\log B}}.                     \tag{4.10}
\]

The implied constant is independent of `z` and `B` for sufficiently large `B`.

**Proof.** Let `F(X)` count positive cubefree split integers at most `X`.
Each such integer is uniquely `a b^2`, with coprime squarefree split `a,b`.
Using (2.4), and then dropping the restrictions on `b`, gives

\[
 F(X)\le\sum_{b\le\sqrt X}A_0(X/b^2)
       \ll\frac{X}{\sqrt{\log(2X)}}.                           \tag{4.11}
\]

Indeed split the sum at `b=X^(1/4)`. The first part uses (2.4) and the
convergence of `sum b^{-2}`. The second uses `A_0(t)<=t` and is `O(X^(3/4))`,
which is bounded by the right side. This imports no new analytic theorem.

By Theorem 4.1 every class counted by `D_Phi(B)` has a cubefree split allocation
`(R,S) in BW`, hence `R,S<=B/3`. If it is not counted by `D_{P_z}(B)`,
every successful allocation has some `p>z` with `p^2|R` or `p^2|S`.
Count these pairs, dropping coprimality and the modulo-nine condition.
For `p<=B^(1/4)`, (4.11) bounds their number by

\[
 F(B/(3p^2))F(B/3)\ll\frac{B^2}{p^2\log B}.
\]

For `p>B^(1/4)` use `F(B/(3p^2))<=B/(3p^2)` and (4.11) for the other
coordinate. The sum is `O(B^(7/4)/sqrt(log B))`, since
`sum_{p>B^(1/4)}p^{-2}<<B^(-1/4)`. Interchange the two coordinates and apply
the union bound. A finite collection of pairs surjects onto the omitted
algebra classes, so counting pairs is a valid upper bound. \(\square\)

### Theorem 4.4: the full leading asymptotic

The full distinct Hasse-algebra count on the original target line satisfies

\[
 \boxed{D_\Phi(B)\sim C_*\frac{B^2}{\log B}},\qquad
 C_*=\lim_{z\to\infty}C(P_z),\qquad
 0<\frac{31\kappa_{\rm sf}}{3888}\le C_*\le C_\Phi.              \tag{4.12}
\]

The limit is monotone, and its tail satisfies

\[
 0\le C_*-C(P_z)
       \ll\sum_{\substack{p>z\\p\equiv1(3)}}p^{-2}.             \tag{4.13}
\]

This is an exact positive constant given by convergent finite rational area
sums times (4.1); no multiplicative Euler factorization of the area sum is
asserted. Some exact finite-prime coefficients are

| `P` | `C(P)/kappa_sf` |
|---|---|
| empty | `31/3888` |
| `{7}` | `1747/214326` |
| `{7,13}` | `26155/3184272` |
| `{7,13,19}` | `18113057423/2196736908912` |

**Proof.** The sets counted by `D_P` increase with `P`, and are contained in
those counted by `D_Phi`. Their constants therefore increase, and are bounded
above by `C_Phi` using the manuscript's rational target count. They have a
finite positive limit `C_*`. For every fixed `z`, (4.8) and (4.10) bound the
lower and upper normalized limits of `D_Phi` between `C(P_z)` and
`C(P_z)+O(sum_{p>z}p^{-2})`. Let `z` tend to infinity. The limits coincide,
proving (4.12); the same comparison gives (4.13). \(\square\)

## 5. Verification and dependencies

The proof uses the exact rational local criterion, the pure-cubic
isomorphism theorem, and the indicated Selberg--Delange theorem. Corollary
3.2 additionally uses the manuscript's rational target-count theorem.
The coprimality correction, exact height normalization, and both uniform tail
bounds are proved above.

Run `python3 research/arithmetic/verify_fixed_map_distinct_algebra_count.py --write-receipt`.
The [checker](verify_fixed_map_distinct_algebra_count.py) verifies finite
Euler coefficients and twisted correction factors, bounded exact
coprimality convolutions, cube-class representatives, target heights,
the two-orbit count identity (3.2), finite-prime area sums, optimized heights
against bounded cube moves, and bounded complete class counts. The
[receipt](fixed_map_distinct_algebra_count_verification.json) records its
bounds and source hashes. Full class enumeration reaches height 1000;
cube multiplier controls reach 42 in each direction, and exact area sums use
the first four split primes through 31. Finite checks do not prove the unbounded
asymptotic or replace the imported field-isomorphism and analytic proofs.
No formal verification, independent whole-proof replay, external review
or publication is claimed.
