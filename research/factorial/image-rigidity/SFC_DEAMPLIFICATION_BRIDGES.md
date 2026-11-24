# Exact de-amplification bridges and their limits

Internal research proofs and countercontrols, 21 November 2026. Repository
base: `19eb9912fb60142f8b620dac24896931d4ff3e46`.

## Status

There is now an exact arithmetic approximation to removing equal private
factorial powers in any fixed finite moment range. There is also a complete
description of why several proposed transports of regularity do not follow
from this approximation or from invertible coefficient operators. These
results do not remove private factors from general ordinary factorial
supports. Full SFC remains open.

| Result | What is proved | What is not implied |
|---|---|---|
| Digit-window approximation | Actual normalized factorial moment forms approach the undecorated forms modulo every prescribed prime power in a fixed finite range. | Nonvanishing of the approximating systems does not prove nonvanishing of the limiting system. |
| Pure-coefficient invariant | Balanced moment systems retain an explicit polynomial identity under coefficient rescaling, degree normalization, and coefficientwise limits. | Such limits cannot recover arbitrary unbalanced supports. |
| Operator algebra | Private factorial maps are exact isomorphisms to a specified twisted product; ordinary ideal transport fails even for ordinary regular factorial frames. | Invertible linear or Koszul-operator transport does not identify the ordinary moment ideals. |
| Reciprocal instability | An unbounded family satisfying the proved amplification theorem has an unstable reciprocal quadratic and two positive Hessian directions. | Reciprocal stability and a one-positive-direction Hessian cannot be assumed as universal amplification invariants. |

The positive shared-row merger theorem is in
[Strong Factorial from separating shared rows](SFC_SEPARATING_ROW_AMPLIFICATION.md).
The arguments below explain which further steps remain unsupported.
All countercontrols distinguish failure of an auxiliary mechanism from
failure of SFC. No external review, formal verification, or novelty is claimed.

Notation is local to each part: Part A uses actual moments; Part B divides
degree-`m` moments by `m!`; Part C uses actual moments and explicitly defined
coefficient reciprocals. These nonzero degree scalars preserve common zero sets.

## A. Arithmetic approximation and its limiting rank

### A. Scope

Fix a finite maximum moment order `M`. There are arbitrarily large
positive private exponents `D` for which the factorial decoration is
arbitrarily close, coefficient by coefficient, to no decoration at all.
The construction below is exact and uses integer digit arithmetic.

This does **not** transfer regularity from those large exponents to
`D=0`. An exact factorial example shows that the required determinant
valuations can grow without bound. A second example has distinct base
supports and exhibits arbitrarily poor conditioning of a prefix shorter
than its full support-size window, without violating SFC.

The large-private-exponent theorem being investigated is
`research/factorial/image-rigidity/SFC_BALANCED_PRIVATE_AMPLIFICATION.md`.
No new unamplified SFC closure is asserted here.

### A.1. Normalize by D! once per moment power

Put

\[
 K_D(k)=\frac{(Dk)!}{(D!)^k}\in\mathbb Z_{>0},\qquad K_D(0)=1.
 \tag{A1}
\]

For a phase with fresh private variables,

\[
 F_D=\sum_{j=1}^s c_jX_j^D M_j(Y),
\]

its actual moment, divided by the nonzero scalar `(D!)^m`, is

\[
 \widehat P_{m,D}(c)=\frac{\Lambda(F_D^m)}{(D!)^m}
 =\sum_{|k|=m}\binom{m}{k_1,\ldots,k_s}
          \Lambda\!\left(\prod_jM_j^{k_j}\right)
          \prod_jK_D(k_j)c^k.
 \tag{A2}
\]

This is an integral coefficient polynomial when the base monomials
have integer exponents. At `D=0`, all kernels in (A1) equal one and (A2)
is exactly the undecorated moment. There is no degree-dependent
rescaling of the coefficient variables.

### A.2. A simultaneous digit-window congruence

**Theorem.** Let `p` be a prime, `M<p`, and `r>=1`. Set

\[
 L=\varphi(p^r)=(p-1)p^{r-1},\qquad
 D_r=\sum_{i=1}^{L}p^{ri}.
 \tag{A3}
\]

Then

\[
                         K_{D_r}(k)\equiv1\pmod{p^r}
                              \quad(0\le k\le M).
 \tag{A4}
\]

The numbers `D_r` tend to infinity. In particular they eventually
exceed every fixed private-amplification threshold.

**Proof.** Define the product of the units through `N` by

\[
 W_p(N)=\prod_{1\le j\le N,\ p\nmid j}j.
\]

Splitting off multiples of `p` repeatedly gives the exact identity

\[
 N!=p^{v_p(N!)}\prod_{e\ge0}W_p(\lfloor N/p^e\rfloor).
 \tag{A5}
\]

All products indexed by `e` are finite. The digits of `D_r` are one
at positions `r,2r,...,Lr`, and zero elsewhere. Because `k<p`, adding
`D_r` to itself `k` times produces no carries. In particular

\[
 \lfloor kD_r/p^e\rfloor=k\lfloor D_r/p^e\rfloor
                      \quad\text{for every }e\ge0,
\]

and `v_p(K_(D_r)(k))=0`.

For fixed `k`, put

\[
 R_k(N)=\frac{W_p(kN)}{W_p(N)^k}\in\mathbb Z_{(p)}^\times.
\]

Modulo `p^r`, this function has period `p^r`. Indeed every interval of
`p^r` consecutive integers contains each unit residue exactly once.
If `B` denotes their product, then
`W_p(N+p^r) congruent B W_p(N)`, and the multiplier in the numerator
of `R_k(N+p^r)` is `B^k`, cancelling the denominator multiplier.
This argument does not need `B=-1` and is valid also at the prime two.

Equation (A5) now yields

\[
 K_{D_r}(k)\equiv
       \prod_{e\ge0}R_k(\lfloor D_r/p^e\rfloor)\pmod{p^r}.
 \tag{A6}
\]

Each factor depends only on the `r` consecutive base-`p` digits in
positions `e,...,e+r-1` of `D_r`. Such a window contains at most one
nonzero digit. Each digit contributes exactly once at each relative
position `0,...,r-1`: its absolute position is at least `r`, so no
window is truncated by `e<0`. The windows attached to consecutive
nonzero digits are disjoint. Therefore

\[
 K_{D_r}(k)\equiv
       \left(\prod_{j=0}^{r-1}R_k(p^j)\right)^L\pmod{p^r}.
\]

The parenthesized expression is a unit modulo `p^r`, and Euler's
theorem with `L=phi(p^r)` gives (A4). The cases `k=0,1` also follow
directly from (A1). QED.

For `r=1`, the windows have length one. Consecutive nonzero digits
are allowed; the construction has `p-1` consecutive digits one above
the units position. For `r>=2`, the spacing is part of the proof;
it cannot simply be dropped.

If `p=2` and `M<p`, only `k=0,1` occur, and (A4) is trivial. Every
nontrivial moment range `M>=2` therefore uses an odd prime `p>M`.
There is an exact reason for this restriction, not only a proof
convenience: for every `D>0`,

\[
                         v_p(K_D(p))=s_p(D)>0,
 \tag{A7}
\]

where `s_p(D)` is the base-`p` digit sum. Thus no positive `D` can
satisfy `K_D(p) congruent1 mod p`.

### A.3. Computing valuations without enormous factorials

For any `D,k>=1`, Legendre's elementary floor formula gives

\[
 v_p(K_D(k))
    =\sum_{e\ge1}\bigl(\lfloor kD/p^e\rfloor
                              -k\lfloor D/p^e\rfloor\bigr)
    =\frac{k s_p(D)-s_p(kD)}{p-1}.
 \tag{A8}
\]

For a modulus `P=p^t`, precompute the prefix products

\[
 A(j)=\prod_{1\le i\le j,\ p\nmid i}i\pmod P,
                      \quad0\le j<P,
\]

and the complete block product `B=A(P-1)`. Then

\[
 W_p(N)\equiv B^{\lfloor N/P\rfloor}A(N\bmod P)\pmod P.
 \tag{A9}
\]

Equations (A5), (A8), and (A9) compute the valuation and unit part of
`K_D(k)` using at most a digit-length number of modular products.
They never form `(Dk)!`. If `K_D(k)` is a unit and its residue differs
from one modulo `p^t`, the valuation of that residue minus one is the
exact value of `v_p(K_D(k)-1)`. Otherwise the computation certifies
the lower bound `v_p(K_D(k)-1)>=t`.

Two exact negative controls for (A3) are

\[
 \begin{array}{c|c|c}
 p^r&D&K_D(2)\bmod p^r\\\hline
 9&3+3^2+\cdots+3^6=1092&7\\
 9&3^2+3^4=90&4
 \end{array}
\]

The first keeps `L=phi(9)=6` but uses insufficient digit spacing.
The second keeps the correct spacing but uses only `L=p-1=2`
instead of `phi(9)`. Neither kernel is congruent to one.

The natural simpler sequences also have different limits. For
`D=p^h`, equation (A8) gives `v_p(K_D(k))=v_p(k!)`; if `k<p`, its
reduction is `k!`, not one in general. To see the reduction directly,
use (A6) with modulus `p`: every factor `R_k(p^j)` with `j>=1`
equals `R_k(0)=1`, while the last factor is `R_k(1)=k!`.
For `D=p^h-1` and
`1<=k<=p^h`, digit complementation gives

\[
                         v_p(K_D(k))=h(k-1).
 \tag{A10}
\]

Thus every fixed `k>=2` tends toward zero along this second sequence,
not toward the undeformed kernel one. Removing a single Morita unit
factor would miss the lower-scale factorial contributions in (A5).

### A.4. Actual moment forms do converge to the undecorated forms

For any fixed support and all `m<=M`, (A2) and (A4) imply

\[
                \widehat P_{m,D_r}(c)\equiv P_{m,0}(c)
                                      \pmod{p^r}
 \tag{A11}
\]

coefficient by coefficient. For fixed rational polynomial templates
the same statement holds after choosing `p` away from their finitely
many coefficient denominators. The congruence concerns the actual
factorial moments, not a surrogate coefficient sequence.

It removes the newly appended equal private powers in the balanced
amplification theorem. It does not automatically remove the separate
fixed private factors used to balance an arbitrary original support.
Scaling those missing-profile factors changes the balance conditions.
Even for an already balanced support, the rank issue below remains.

This distinction is forced by an exact algebraic invariant, not only
by the choice of digit construction. Suppose a family of base columns
has common pure moment `H_m` at every order. Allow arbitrary nonzero
coefficient rescalings `c_j -> lambda_j c_j` and an arbitrary nonzero
scalar `b_m` multiplying the entire degree-`m` moment. Its pure
coefficient is then

\[
 A_{m,j}=b_m H_m\lambda_j^m.
\]

For every pair of columns `j,h` it satisfies the polynomial identity

\[
 A_{m,j}(A_{1,h})^m-A_{m,h}(A_{1,j})^m=0.
 \tag{A11a}
\]

Consequently every coefficientwise limit of such systems, whether
complex or `p`-adic, satisfies the same identity. This remains true
when the factors used to balance the profiles, the variable count,
or the rescaling factors vary with the approximating system: (A11a)
is an identity in the finitely many displayed coefficients, and uses
no division.

The distinct unbalanced base columns `X` and `Y^2` violate (A11a)
already at `m=2`. Their pure coefficients at orders one and two are
`(1,2)` and `(2,24)`, and

\[
                         2\cdot2^2-24\cdot1^2=-16\ne0.
 \tag{A11b}
\]

Thus no family of exactly balanced lifts can converge to this
undecorated moment system using only independent coefficient
rescalings and degreewise scalar normalizations. Formula (A11) is
valid for any fixed base, but the known regularity theorem applies
to its large-`D` members only when that fixed base meets the balance
hypothesis. Starting with the standard balancing lift of an arbitrary
base, (A11) retains those fixed balancing factors. It cannot remove
them by an unmentioned second limiting step. This obstruction does
not address more general transformations that mix coefficient
variables or a separate de-amplification theorem.

### A.5. An exact actual-factorial determinant obstruction

Take two identical base columns `M_1=M_2=1`. The amplified phase is
`c_1 X_1^D+c_2 X_2^D`, whose terms are distinct for every `D>0`.
The normalized first two moments are exactly

\[
 P_{1,D}=c_1+c_2,\qquad
 P_{2,D}=K_D(2)(c_1^2+c_2^2)+2c_1c_2.
\]

In the degree-two column basis `(c_1^2,c_1c_2,c_2^2)`, the rows
`c_1P_1,c_2P_1,P_2` have determinant

\[
 \det\begin{pmatrix}1&1&0\\0&1&1\\K_D(2)&2&K_D(2)\end{pmatrix}
                       =2\left[\binom{2D}{D}-1\right].
 \tag{A12}
\]

It is nonzero for every `D>0`. Along (A3), for any odd `p>2`, its
valuation is at least `r`, without a uniform bound. At `D=0` it is
zero. Thus actual large-private-exponent regularity does not provide
a uniform integral conditioning bound, even when the kernels in (A11)
converge exactly as desired.

The unamplified columns in this example coincide: `c_1+c_2=0` gives
the zero phase. Therefore (A12) is **not** an SFC counterexample and
does **not** refute a removal theorem restricted to distinct original
columns. It refutes an automatic nonvanishing-to-limit argument and
a universal conditioning bound based only on the amplification theorem.

### A.6. Distinct supports with a badly conditioned shorter prefix

Let

\[
 M_1=YZ,\quad M_2=YW,\quad M_3=UV,
 \qquad(c_1,c_2,c_3)=(1,x,-1-x),\quad3x^2+4x+3=0.
\]

These are distinct balanced columns, each with profile `(1,1)`, and
all three coefficients are nonzero. Direct expansion gives the
undecorated moments

\[
              P_{1,0}=0,\qquad P_{2,0}=0,
              \qquad P_{3,0}=60+20x\ne0.
 \tag{A13}
\]

Indeed `H_2=4`, the pairwise mixed moments are `H_12=2` and
`H_13=H_23=1`, and `sum c_j^2=-2x/3`. At order three the pure
moment is `36`, the mixed moments for the special pair `(1,2)`
are `12`, those for the other pairs are `4`, and the fully mixed
moment is `2`; these give `60 c_1c_2c_3=60+20x`.

For positive private exponent `D`, normalized as in (A2),

\[
 \widehat P_{1,D}=0,\qquad
 \widehat P_{2,D}=-\frac{8x}{3}\,[K_D(2)-1].
 \tag{A14}
\]

This second value is nonzero for every `D>0`. Choose `p=7` and the
construction (A3). The two reductions of `x` are `3,5`, both units;
`60+20x` reduces respectively to `1,6`. Consequently (A14) has
valuation at least `r`, while the third moment is a unit at both
places above seven by (A11).

This supplies an actual distinct-support example with arbitrarily
poor conditioning in the first two orders. Its full support size is
three, and the third order is nonzero. It therefore supports the
sharp distinction between a shorter prefix and the SFC window,
without contradicting distinct-support de-amplification.

### A.7. The exact missing full-rank condition

Fix one support, one starting order `n`, one active support size `q`,
and one Macaulay degree. For the full window choose `p>n+q-1`;
the prime need not exceed the Macaulay degree, because its extra
coefficient degree comes only from monomial multiples. Matrix
dimensions must stay fixed as the precision varies. All orders
and all active subsupports require separate rank assertions.

The remaining issue can be stated without any factorial estimates.
Let integral matrices `A_r` of a fixed size satisfy
`A_r congruent A_0 mod p^r`. Let `mu(A)` be the least valuation among
its full-column-size minors, with `mu=+infinity` if they all vanish.

* If `A_0` has full column rank and `mu(A_0)=b<infinity`, then
  `mu(A_r)=b` for every `r>b`.
* If `A_0` lacks full column rank, then `mu(A_r)>=r` for every `r`.

Each assertion follows because determinants are integer polynomials
in the entries. It applies in particular to any fixed-degree Macaulay
matrices formed from (A11). Thus an eventual uniform bound on the
minimum valuations of full minors along the approximating sequence
is equivalent to the
missing full rank of the limiting matrix; it cannot be assumed as a
weaker consequence of large-`D` regularity. A concrete certificate with
valuation `<r` would prove full rank at that precision, but by (A11)
it is also already a direct certificate for the undecorated matrix.

The digit-window theorem solves the approximation problem. It does
not solve this determinant-conditioning problem, nor does it supply
private-factor removal or unrestricted SFC.

## B. Factorial operators and ordinary ideals

Internal exploratory proof note, 21 November 2026. This note proves the
operator statements and countercontrols below. It does **not** prove
ordinary factorial de-amplification or unrestricted Strong Factorial.

### B.1. The operator and the ordinary ideal are different objects

Let `K` be a characteristic-zero field, `R=K[c_1,...,c_s]`, and let
`A=(a_ij)` have distinct nonnegative integer columns. Write

\[
 P_m(c)=\frac{\Lambda((\sum_j c_j X^{a_j})^m)}{m!}
 =\sum_{|k|=m}\frac{\prod_i((Ak)_i)!}{\prod_jk_j!}c^k.
 \tag{B1}
\]

Projective regularity throughout means that the forms have no nonzero
common zero over an algebraic closure of `K`, not merely over `K`.

Adjoin any finite collection of private rows `d e_j`, with positive integers
`d`. Set

\[
 w(k)=\prod_{(j,d)\text{ private}}(d k_j)!,\qquad
 T(c^k)=w(k)c^k,\qquad C=T^{-1}.
 \tag{B2}
\]

The decorated moment forms are exactly `Q_m=T P_m`. Every weight is a
nonzero rational integer. Thus `T` is an invertible, degree-preserving
linear map. It is generally not a ring map, and

\[
 C\bigl((Q_n,\ldots,Q_{n+s-1})\bigr)
 \ne (P_n,\ldots,P_{n+s-1})
 \tag{B3}
\]

even for ordinary factorial frames for which both ideals are regular.
Section B.3 supplies exact witnesses to both failed inclusions in (B3).

For one private row `D e_j`, put `theta_j=c_j partial_(c_j)`. On polynomials,

\[
 T c_j T^{-1}=c_j\prod_{\ell=1}^D(D\theta_j+\ell),\qquad
 T^{-1}c_jT=c_j\left[\prod_{\ell=1}^D(D\theta_j+\ell)\right]^{-1}.
 \tag{B4}
\]

The bracketed inverse is a well-defined diagonal operator on `R`; all its
eigenvalues there are positive integers. More generally, for a monomial
multiplier `c^beta`,

\[
 C(c^\beta T f)
 =c^\beta\frac{w(\theta)}{w(\theta+\beta)}f.
 \tag{B5}
\]

For one private row this ratio is
`prod_(ell=1)^(D beta_j)(D theta_j+ell)^(-1)`. Equation (B5) is an exact
finite product, not an asymptotic Gamma formula. It explains why an
amplified ideal certificate becomes a certificate using Euler resolvents
of the undecorated forms, rather than the original forms themselves.

#### The valid twisted-algebra statement

Define `f star g = T((Cf)(Cg))`. This is an associative commutative graded
product with unit one, and

\[
 c^\alpha\star c^\beta
 =\frac{w(\alpha+\beta)}{w(\alpha)w(\beta)}c^{\alpha+\beta}.
 \tag{B6}
\]

Then `T` is an algebra isomorphism from the ordinary ring to the star
algebra, and `T((P_n,...,P_(n+s-1)))` is precisely the **star ideal**
generated by the `Q_m`. Koszul complexes transport under this isomorphism.
The amplification theorem proves regularity of the `Q_m` with **ordinary**
multiplication, so this identity does not supply the desired implication.

Equivalently, conjugating multiplication by `P_m` gives `P_m(B_1,...,B_s)`,
where `B_j=T c_j T^(-1)`; applying this operator to `1` gives `Q_m`.
The operator itself is not multiplication by `Q_m`. Evaluation on this
one vector is insufficient to identify the two Koszul complexes.

### B.2. A general classification and a localization obstruction

**Proposition 1.** Let `s>=2` and let `T(c^k)=w(k)c^k`, with all weights
nonzero and `w(0)=1`. The image `T(I)` is an ordinary ideal for every
homogeneous ideal `I` if and only if

\[
                   w(k)=a_{|k|}\prod_j b_j^{k_j}
 \tag{B7}
\]

for nonzero scalars `a_m,b_j`. One may normalize `a_0=a_1=1`.

**Proof.** Put `l=sum_j c_j` and `b_j=w(e_j)`. If `T((l))` is an ideal,
it contains `(T l)`. Both have the same dimension in every homogeneous
degree, so they are equal. For every `beta`, the polynomial
`c^beta sum_j w(beta+e_j)c_j` is therefore divisible by `sum_j b_jc_j`.
That linear form is relatively prime to the monomial `c^beta`, whence
`w(beta+e_j)/b_j` is independent of `j`. Moving one unit between
coordinates connects all compositions of any fixed degree. Consequently
`w(k)/prod_j b_j^(k_j)` depends only on `|k|`, giving (B7). Conversely,
on each homogeneous piece (B7) is a scalar times the same diagonal change
of variables, which sends homogeneous ideals to ideals. QED.

Every nonempty private factorial decoration violates (B7): for a decorated
column `j` and any different column `h`,

\[
 \frac{w(2e_j)}{w(e_j)^2}
      =\prod_{d\text{ on }j}\binom{2d}{d}>1,
 \qquad
 \frac{w(e_j+e_h)}{w(e_j)w(e_h)}=1.
 \tag{B8}
\]

This excludes universal ideal transport, not preservation of the much
more specific assertion that a consecutive moment ideal is primary.

**Proposition 2.** The conjugation in (B4) has no extension by an invertible
linear map to the Laurent ring `K[c_1^(+/-1),...,c_s^(+/-1)]` that preserves
the displayed Euler-operator formula.

**Proof.** Multiplication by `c_j` is invertible on that ring. Its conjugate
by an invertible linear map would therefore be invertible. But

\[
 \left[c_j\prod_{\ell=1}^D(D\theta_j+\ell)\right]c_j^{-1}=0,
 \tag{B9}
\]

since the factor with `ell=D` vanishes. The vector `c_j^(-1)` is nonzero.
The same argument applies to any nonempty product of private factors on
that column. QED.

Thus moving first to the coefficient torus does not remove the operator
gap: it creates a genuine resonance at negative exponent `-1`. This does
not exclude other arguments using localization; it excludes this purported
invertible extension of the factorial conjugation.

There is also no positive independent Mellin variable realizing `C`.
Already for one private row its proposed moments would satisfy
`E Z=1/D!`, `E Z^2=1/(2D)!`, with strictly negative variance. A positive
averaging interpretation of inverse markers is therefore unavailable.

### B.3. An ordinary factorial countercontrol, with exact nonmembership

Take the actual base features `X,Y`, independent unit exponentials, and
the marked features `XZ,Y`, where `Z` is another independent unit
exponential. Use coefficients `x,y`, and divide each moment by `m!`.
For starting order two the consecutive forms are

\[
 \begin{array}{ll}
 p_2=x^2+xy+y^2,&p_3=x^3+x^2y+xy^2+y^3,\\
 q_2=2x^2+xy+y^2,&q_3=6x^3+2x^2y+xy^2+y^3.
 \end{array}
 \tag{B10}
\]

Here `T(x^k y^l)=k! x^k y^l`, so `q_i=T p_i` exactly. The dehomogenized
Sylvester resultants are

\[
      \operatorname{Res}_x(p_2(x,1),p_3(x,1))=1,
 \qquad
      \operatorname{Res}_x(q_2(x,1),q_3(x,1))=36.
 \tag{B11}
\]

There is no common zero at `y=0`, since the leading pure coefficients
are nonzero. Thus both homogeneous ideals are primary to `(x,y)`.

Nevertheless

\[
 C(xq_2)=\tfrac13x^3+\tfrac12x^2y+xy^2
                  \notin I=(p_2,p_3).
 \tag{B12}
\]

An explicit linear nonmembership certificate is
`L(f)=[x^2y]f-[xy^2]f`: it vanishes on `xp_2,yp_2,p_3`, which span
`I_3`, but gives `-1/2` on (B12). The reverse inclusion also fails:

\[
 T(xp_2)=6x^3+2x^2y+xy^2\notin J=(q_2,q_3).
 \tag{B13}
\]

Here `L_tilde(f)=[x^2y]f-[xy^2]f-[y^3]f` vanishes on
`xq_2,yq_2,q_3`, while its value on (B13) is `1`.

At the critical Macaulay degree four, even a row/column diagonal
identification of the two coefficient matrices fails. The entrywise
marker multipliers for rows `x^2 p_2,xy p_2` and columns
`x^3y,x^2y^2` form

\[
                         \begin{pmatrix}1&1\\2&1\end{pmatrix}.
 \tag{B14}
\]

Its multiplicative cross-ratio is `1/2`, whereas a separable row/column
scaling would have cross-ratio one. No entries used in this comparison
vanish. This only rejects that proposed direct matrix equivalence.

**Saturation boundary.** Both `I` and `J` in this example have unit
saturation by `(x,y)`; both also become the unit ideal after inverting
`xy`. We do not claim that these saturated ideals differ. Equations
(B12)--(B14) identify the failed transport step before saturation, while
Proposition 2 shows why the proposed conjugation cannot simply be
continued across that localization.

### B.4. The precise missing finite-dimensional implication

Let `d_i=n+i`, `0<=i<s`, and put

\[
 N=1+\sum_{i=0}^{s-1}(d_i-1)
   =1+s(n-1)+\frac{s(s-1)}2.
 \tag{B15}
\]

For `s` homogeneous forms of these degrees in `s` variables, having only
the common zero zero **over an algebraic closure of `K`** is equivalent
to surjectivity of the multiplication map

\[
 \bigoplus_{i=0}^{s-1}R_{N-d_i}\longrightarrow R_N,
 \qquad (h_i)\longmapsto\sum_i h_iP_{d_i}.
 \tag{B16}
\]

For clarity, the forward direction uses the Cohen--Macaulay property of a
polynomial ring: a homogeneous system of parameters is a regular sequence.
Its successive multiplication exact sequences give Hilbert series
`prod_i(1+...+t^(d_i-1))`, so the quotient vanishes in degree `N`.
Conversely surjectivity puts every `c_j^N` in the ideal and excludes a
nonzero common zero. The standard algebraic input is the dimension
criterion in [Stacks, Tag 02JN](https://stacks.math.columbia.edu/tag/02JN);
the Koszul implications and their local hypotheses are explicitly stated
in [Stacks, Tag 062D](https://stacks.math.columbia.edu/tag/062D).

Amplification proves (B16) for the decorated forms. Applying `C` to its
rows gives (B5), not the undecorated multiplication map. Proving full row
rank for the latter remains exactly the missing implication. Neither
the invertibility of `C` on `R_N` nor the scalar factorial shift identities
identify the row spaces. The ordinary example above already prevents the
stronger row-space/ideal identities often proposed at this step.

In particular, taking the ordinary ideal generated by **all** transformed
elements `C(J)` conceals the issue: if `J` contains powers of every
coordinate, then this new ideal does too, because `C` scales monomials.
It includes extra Euler-resolvent equations, and is not known to be the
ideal generated by the original consecutive forms.

### B.5. What a uniform removal theorem would mean

Suppose one proves, for all distinct ordinary supports and all starting
orders, that deleting any private row `d e_j` preserves projective
regularity whenever the remaining columns are distinct and nonzero.
This would imply full SFC. Indeed first adjoin a common row of ones to
any target support; this multiplies its degree-`m` moments by `m!`, so
preserves every zero set and makes all columns nonzero. The retained
balanced-amplification construction then adds only private rows. Remove
those rows one at a time. The original distinct columns remain as a
submatrix throughout, and the common row keeps every column nonzero.
At the end discard that common scalar factor.

Conversely full SFC immediately gives this preservation assertion.
Thus this particular universal removal rule is an **equivalent reformulation**,
not an independently easier theorem already established by amplification.
The same warning applies to requiring the implication only for the one
explicit amplified enlargement of each support: its antecedent is now
proved for every support.

The retained positive-measure counterexample is stronger than a failure of
one isolated window: its source is regular at every starting order after
private degree `D>=12`, but the undecorated target has the two forbidden
moments. Its remaining feature law is not the ordinary factorial law.
Consequently it does not refute the integer-toric removal rule above;
it shows that a proof must use more than positive coupling, strict Gram
positivity, equal pure moments and mixed-moment bounds.

### B.6. A valid marker-profile extension

One actual extension of the existing amplification argument is immediate.
Let the balanced base monomials have common pure moments `H_m`. In every
column use an independent copy of the **same arbitrary private monomial**
with positive exponent profile `(b_1,...,b_r)` and total degree
`D=sum_h b_h`. Its variables are private to that column. If `D>=6s`,
the same all-order SFC conclusion and quantitative bound hold.

For a composition `k` of `m`, the normalized coefficient is

\[
 C_k\frac{H(k)}{H_m}
 \prod_{h=1}^r\frac{\prod_j(b_h k_j)!}{(b_hm)!}
                  \le C_k^{1-D}.
 \tag{B17}
\]

Apply the retained coloring injection separately for each `b_h`; the
rest of the Turan proof is unchanged. This includes `D` independent
unit-row markers per column, in which case the product ratio is exactly
`C_k^(-D)`. Zero coefficients are deleted as before, and characteristic-zero
descent uses the same finitely generated coefficient field argument.

This extension keeps all markers and still requires balanced base pure
moments. In particular it does not dispose of the arbitrary private
completion powers used to balance an originally unbalanced support.
It therefore does not prove that unit-row removal by itself suffices
without another bridge for those completion powers.

## C. Reciprocal instability despite arbitrarily strong amplification

Internal mathematical note, 21 November 2026. Repository base
`19eb9912fb60142f8b620dac24896931d4ff3e46`. This note proves an obstruction to
one proposed auxiliary invariant, not a counterexample to Strong Factorial.

### C.1. The auxiliary invariant and the exact conclusion

For distinct monomials `N_1,...,N_s` in independent factorial variables, write

\[
 P_m(c)=\Lambda\left((\sum_j c_jN_j)^m\right)
       =\sum_{|k|=m}p_kc^k,
 \qquad Q_m(c)=\sum_{|k|=m}\frac{c^k}{p_k}.
 \tag{C1}
\]

All `p_k` are positive integers, so the coefficientwise reciprocal `Q_m` is
well-defined. A real polynomial is real stable if it has no zero with all
variables in the open upper half-plane. This note gives an unbounded family
for which `Q_2` is not real stable, even though the original family is known
to satisfy SFC at every starting order. Its quadratic Hessian also has two
positive eigenvalues, violating the necessary quadratic Hessian condition
for a Lorentzian polynomial.

**Theorem.** For every integer `D>=18`, put `a=4D` and take

\[
\begin{aligned}
 N_1&=X_1^D Y_1^{2a}Y_2^a,\\
 N_2&=X_2^D Y_1^aY_2^{2a},\\
 N_3&=X_3^D Y_2^{2a}Y_3^a.
\end{aligned}
\tag{C2}
\]

Every nonzero coefficient combination of these monomials satisfies SFC over
every characteristic-zero field, with its active support count as the
window length. Nevertheless, its reciprocal quadratic `Q_2` has inertia
`(2 positive, 1 negative)` and is not real stable. The inertia conclusion
already holds for every `D>=2` in (C2); the threshold `18` is used solely to
import the all-order SFC theorem.

This rules out reciprocal real stability or a one-positive-eigenvalue
quadratic invariant as a universal intermediate consequence of balanced
private amplification. It makes no claim that either invariant is necessary
for SFC, and it does not rule out a more selective transport mechanism.

### C.2. Direct factorial computation

The three base exponent columns are

\[
 (2a,a,0),\qquad(a,2a,0),\qquad(0,2a,a).
\]

Their positive exponent profiles agree. Define

\[
\begin{aligned}
 H&=(4a)!(2a)!,\\
 H_{12}&=((3a)!)^2,\\
 H_{13}&=(2a)!(3a)!a!,\\
 H_{23}&=(a!)^2(4a)!.
\end{aligned}
\tag{C3}
\]

The pure coefficient of `P_2` is `(2D)!H`; its coefficient of `c_ic_j` is
`2(D!)^2H_ij` for `i!=j`. Set

\[
 b_D=\frac14\binom{2D}{D},\qquad R_{ij}=\frac{H}{H_{ij}}.
\]

Then the positively rescaled reciprocal quadratic has the exact matrix

\[
 (2D)!H\,Q_2(c)=c^{\mathsf T}Mc,\qquad
 M=\begin{pmatrix}
 1&A&B\\ A&1&C\\ B&C&1
 \end{pmatrix},
\quad
 A=b_DR_{12},\ B=b_DR_{13},\ C=b_DR_{23}.
 \tag{C4}
\]

The factor `1/4` in `b_D` accounts both for the multinomial coefficient
`2` in `P_2` and for the two symmetric off-diagonal matrix entries.

Factorial cancellation gives

\[
\begin{aligned}
 R_{12}&=\frac{(4a)!(2a)!}{((3a)!)^2},&
 R_{13}&=\binom{4a}{a},&
 R_{23}&=\binom{2a}{a},\\
 \frac{R_{13}}{R_{12}R_{23}}
 &=L_a:=\frac{a!(3a)!}{((2a)!)^2}
 =\prod_{i=1}^a\frac{2a+i}{a+i}.
\end{aligned}
\tag{C5}
\]

Each factor in the last product is at least `3/2`. With `a=4D`,

\[
 L_a\ge(3/2)^{4D}=(81/16)^D>4^D>
 \tfrac12\binom{2D}{D}=2b_D.
 \tag{C6}
\]

The second inequality is `81>64`; the last follows from the binomial
sum `sum_k binom(2D,k)=4^D`. Consequently `B>2AC`. Also `R_12>=1`
by factorial log-convexity, `R_23>=1`, and `b_D>=3/2` for `D>=2`.
Thus `A,C>1`, and

\[
 \det M=1+2ABC-A^2-B^2-C^2
       <1-A^2-C^2<0.
 \tag{C7}
\]

The leading principal minors are `1`, `1-A^2<0`, and `det M<0`.
The exact diagonal pivots of an LDL-transformation therefore have signs
`(+,-,+)`. In particular the inertia is `(2,1)`.

Here is a self-contained reason this gives a zero in the required
half-planes. Let `y=(1,1,1)`. All entries of `M` are positive, so
`q(y)=y^TMy>0`. On the two-dimensional positive subspace of the real
quadratic form, the linear condition `x^TMy=0` has a nonzero solution.
For that solution `q(x)>0`; rescale it so `q(x)=q(y)`. Then

\[
             q(x+iy)=q(x)-q(y)+2i\,x^TMy=0,
\]

and all three imaginary parts equal one. Hence `Q_2` is not real stable.

Finally, the base monomials in (C2) have the same positive exponent
profile `(a,2a)`, and the three `X_j` are private. The
[balanced private amplification theorem](SFC_BALANCED_PRIVATE_AMPLIFICATION.md)
applies with `s=3,D>=6s=18`. It proves the stated SFC conclusion for
arbitrary characteristic-zero coefficients. That imported conclusion
is logically separate from the direct reciprocal calculation above.

### C.3. A small weighted-row control

The obstruction already occurs in degree two after adding one integer row.
Start from the monomials `X_1Y_1,X_2Y_2,X_3Y_3`. Their degree-`m` moment
coefficient at `c^k` is `m! prod_j k_j!`, so

\[
                  Q_m(c)=\frac{(c_1+c_2+c_3)^m}{(m!)^2}.
 \tag{C8}
\]

Add one factorial variable `Z` with exponent row `(0,1,2)`. This replaces
the reciprocal coefficient at `c^k` by that coefficient divided by
`(k_2+2k_3)!`. At order two the new reciprocal is

\[
 Q_2(c)=\frac14c^{\mathsf T}
 \begin{pmatrix}
 1&1&1/2\\
 1&1/2&1/6\\
 1/2&1/6&1/24
 \end{pmatrix}c.
 \tag{C9}
\]

The matrix has determinant `-1/144` and exact diagonal pivots
`1,-1/2,1/72`. An explicit upper-half-plane zero is

\[
      c=\frac{(23,-94,140)}{\sqrt{56}}+i(1,1,1).
 \tag{C10}
\]

Indeed, for `x=(23,-94,140)` and `y=(1,1,1)`, direct rational arithmetic
gives `x^TMy=0`, `x^TMx=273`, and `y^TMy=39/8`, with quotient `56`.
This is a failure of an auxiliary real-stability preservation claim for
an actual nonnegative integer factorial row. It is not an SFC counterexample.

## D. Exact replay, source scope, and remaining obligation

The three standard-library checkers verify distinct parts of this note.
They pin this complete note, retain exact inputs and certificates in their
receipts, and reject optimized Python and altered proof bytes. The finite
receipts are controls for the written proofs, not substitutes for their
unrestricted arguments.

| Part | Checker | Receipt |
|---|---|---|
| A | [verify_sfc_deamplification_arithmetic.py](verify_sfc_deamplification_arithmetic.py) | [sfc_deamplification_arithmetic_verification.json](sfc_deamplification_arithmetic_verification.json) |
| B | [verify_sfc_deamplification_operators.py](verify_sfc_deamplification_operators.py) | [sfc_deamplification_operators_verification.json](sfc_deamplification_operators_verification.json) |
| C | [verify_sfc_deamplification_reciprocal.py](verify_sfc_deamplification_reciprocal.py) | [sfc_deamplification_reciprocal_verification.json](sfc_deamplification_reciprocal_verification.json) |

Run each checker from the repository root, for example:

```sh
python research/factorial/image-rigidity/verify_sfc_deamplification_arithmetic.py
python research/factorial/image-rigidity/verify_sfc_deamplification_operators.py
python research/factorial/image-rigidity/verify_sfc_deamplification_reciprocal.py
```

The arithmetic checker compares its modular-unit algorithm with 4,400
direct small-factorial computations before testing 14 large digit
constructions, the exact faulty-digit controls, actual moment forms,
the two determinant/prefix examples, and the limiting invariants.
It never forms factorials of the large constructed exponents.

The operator checker reconstructs both ordinary exponent-matrix examples,
retains exact resultant and linear-nonmembership certificates, tests
Euler conjugation and the Laurent kernel, and rejects false row/column
scalings and positive inverse-moment interpretations.

The reciprocal checker reconstructs the matrices directly from integer
exponent columns for `D=2,3,18,19,24`, retains exact determinants and
LDL pivots, and checks the explicit upper-half-plane zero using rational
bilinear identities. The universal statement for all `D>=2` follows
from the displayed factorial inequalities in Part C.

The ordinary homogeneous-ideal criterion in Part B uses the standard
Cohen--Macaulay dimension criterion. The primary sources inspected are
[Stacks, Tag 02JN](https://stacks.math.columbia.edu/tag/02JN), its proof
through [Tag 00N6](https://stacks.math.columbia.edu/tag/00N6), and the
Koszul statements and local hypotheses in
[Tag 062D](https://stacks.math.columbia.edu/tag/062D).
The local application is at the irrelevant maximal ideal, with common
zeros taken over an algebraic closure; this does not assert a converse
from having no nonzero rational points over a smaller field. The
successive graded multiplication sequences then give the stated Hilbert
series and degree bound. This source use is standard algebra, not a
claim to have formally verified those references.

Retained proof inputs are the
[balanced amplification theorem](SFC_BALANCED_PRIVATE_AMPLIFICATION.md),
[row-merging obstructions](SFC_ROW_MERGING_OBSTRUCTIONS.md), and
[unrestricted mechanism audit](UNRESTRICTED_SFC_MECHANISMS.tex).
Their exact source bytes are checked where they are used. The classical
Turan source and its proof audit remain recorded in the balanced theorem's
source ledger. No claim from an RH proof is used.

The remaining theorem is preservation of the ordinary consecutive
moment ideal's projective regularity when the private factors are
removed from actual distinct integer supports. Part B.5 explains why
the universal removal rule would already imply full SFC. Part A.7
explains why a bound on full-minor valuations cannot be supplied by
renaming that same missing rank assertion. These are open obligations.
