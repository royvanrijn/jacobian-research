# Exact factorial-moment closure under ramified polynomial pullback

Internal research theorem and exact certificates, 21 November 2026.
Repository base and latest main inspected:
`19eb9912fb60142f8b620dac24896931d4ff3e46`.

A polynomial change of the residue coordinate preserves the **entire
normalized factorial-moment sequence**, provided its critical values lie
in specified fibers. The transformed residue is again the ordinary
factorial observable of an explicitly constructed monomial times a linear
form. The construction permits unequal common exponents, repeated
coefficients and zero exponents, and multiplies the support size while
retaining every proved seed window. It is a structural closure with a
prescribed coefficient map; it does not prove arbitrary-coefficient SFC.

Write `Lambda(X^alpha)=product_i alpha_i!` for the ordinary factorial
functional. A phase has a nonzero-moment window of length `w` when, for
every integer `n>=1`, at least one of its moments at orders
`n,...,n+w-1` is nonzero. SFC asks for such a window of length equal to
the number of actual nonzero monomials. All algebraic arguments below
hold over characteristic-zero fields. Splitting-field qualifications are
explicit, and no order or positivity on those fields is assumed.

The transfer proof is self-contained. The cyclic examples prove their
seed windows directly. Separately cited existing seed theorems and the
plane-tree interpretation are not premises of the residue identity.

## 1. Statement, including repeated coefficients

Let `K` be a characteristic-zero field. Fix `t>=1` distinct nonzero scalars
`c_1,...,c_t in K`, positive integers `p_j`, and nonnegative integers
`r_(j,l)`, where `1<=l<=p_j`. Put

\[
 R_j=\sum_{l=1}^{p_j}r_{j,l},\qquad S=\sum_{j=1}^t p_j,
\]

\[
 F_0(Y)=\left(\prod_{j,l}Y_{j,l}^{r_{j,l}}\right)
                   \sum_{j=1}^t c_j\sum_{l=1}^{p_j}Y_{j,l}.
 \tag{1}
\]

Its actual monomial support is `S`, including when some common exponents
vanish or some linear coefficients repeat. Define

\[
\begin{aligned}
 P_0(z)&=\prod_j(1-c_jz)^{p_j},&Q_0(z)&=zP_0(z),\\
 A_0(z)&=\prod_j(1-c_jz)^{R_j},&a_0(z)&=zA_0(z).
\end{aligned}
 \tag{2}
\]

Choose a polynomial `h in K[z]` of degree `d>=1` satisfying

\[
 h(0)=0,\qquad h'(0)=1,\qquad h'(z)\mid Q_0(h(z)).
 \tag{3}
\]

Equivalently, over an algebraic closure every critical value of `h`
belongs to the root set of `Q_0`. This equivalence is proved below.
Set

\[
             Q(z)=\frac{Q_0(h(z))}{h'(z)},\qquad a(z)=a_0(h(z)).
 \tag{4}
\]

Take a finite splitting field `L/K` of `Q`. For every nonzero root
`beta` of `Q`, let `p_beta=ord_beta Q` and `R_beta=ord_beta a`.
Choose any nonnegative integers `r_(beta,l)`, `1<=l<=p_beta`, whose
sum is `R_beta`. Define

\[
 F_h(Z)=\left(\prod_{\beta\ne0}\prod_{l=1}^{p_\beta}
                      Z_{\beta,l}^{r_{\beta,l}}\right)
                  \sum_{\beta\ne0}\frac1\beta
                            \sum_{l=1}^{p_\beta}Z_{\beta,l}.
 \tag{5}
\]

**Theorem.** The data in (5) are well-defined and have the following
properties.

1. The actual support of `F_h` is exactly `dS`; the number of distinct
   nonzero coefficient values in its linear factor is exactly `dt`.
2. For every integer `m>=0`,

   \[
   \frac{\Lambda(F_h^m)}{m!\prod_{\beta,l}(mr_{\beta,l})!}
     =\frac{\Lambda(F_0^m)}{m!\prod_{j,l}(mr_{j,l})!}.
   \tag{6}
   \]

3. Every moment zero and every nonzero moment is preserved between the
   two phases. Any nonzero-moment window of length `w` for the seed
   remains valid with the same length `w` for the pullback. In particular,
   a proved SFC seed with support `S` produces an SFC phase of support
   `dS`, with the stronger original `S`-window.
4. The construction is closed under successive admissible pullbacks.
   Their degrees multiply, while the normalized moment sequence in (6)
   remains unchanged.

The new coefficients are prescribed reciprocals of polynomial preimages.
The theorem does not assert arbitrary freedom in all `dS` target
coefficients. If those preimages belong to `K`, the target is a polynomial
over `K`; otherwise it is over the indicated finite algebraic extension.

## 2. The exact normalized moment observable

Put

\[
 C_m=\prod_{j,l}(mr_{j,l})!,\qquad
 g_m(F_0)=\frac{\Lambda(F_0^m)}{m!C_m}.
\]

The multinomial expansion and factorial cancellation give

\[
\begin{aligned}
 g_m(F_0)
 &=\sum_{\sum k_{j,l}=m}
       \prod_{j,l}\frac{(mr_{j,l}+1)_{k_{j,l}}}{k_{j,l}!}
                              c_j^{k_{j,l}}\\
 &=[z^m]\prod_{j,l}(1-c_jz)^{-mr_{j,l}-1}\\
 &=[z^m]A_0(z)^{-m}P_0(z)^{-1}\\
 &=\operatorname{Res}_{z=0}\frac{a_0(z)^{-m}}{Q_0(z)}\,dz.
\end{aligned}
 \tag{7}
\]

All these are formal identities in characteristic zero. Each displayed
coefficient is a finite sum. This also proves that the normalized
sequence depends only on the group totals `(c_j,p_j,R_j)`, not on the
allocation of `R_j` among equal-coefficient variables. The factorial
normalization itself does depend on that allocation.

For comparison with the retained monomial-linear observable, define
the formal primitive `H_0(0)=0`, `H_0'=a_0/Q_0=A_0/P_0`, and
`B_0=H_0 composed with a_0^{-1}`. Formal Lagrange inversion gives

\[
                  g_m(F_0)=(m+1)[u^{m+1}]B_0(u).
 \tag{8}
\]

The primitive is polynomial if every `R_j>=p_j`; the transfer does not
require this restriction. The rational derivative is regular at zero,
so its normalized formal primitive always exists.

## 3. Divisibility, local multiplicities, and exact support count

For a finite point `beta`, write `e_beta` for the local degree of `h`
there, so `ord_beta(h-h(beta))=e_beta` and
`ord_beta h'=e_beta-1`. These identities use characteristic zero.

If `h(beta)` is not a root of `Q_0` but `beta` is critical, the
quotient `Q_0(h)/h'` has a pole there. If `h(beta)` is a root of
`Q_0` of multiplicity `p>=1`, the quotient has order

\[
                         e_\beta p-(e_\beta-1)
                              =e_\beta(p-1)+1\ge1.
 \tag{9}
\]

Thus (3) is equivalent to the critical-value condition, and the roots
of `Q` are exactly the preimages of the roots of `Q_0`. At zero both
`Q` and `a` have derivative one and a simple zero. At any nonzero
preimage `beta` of `alpha_j=1/c_j`, the exact multiplicities are

\[
             p_\beta=e_\beta(p_j-1)+1,
             \qquad R_\beta=e_\beta R_j.
 \tag{10}
\]

At a nonzero preimage of zero they are

\[
                       p_\beta=1,\qquad R_\beta=e_\beta.
 \tag{11}
\]

The zero fiber can be included by the convention `(p_0,R_0)=(1,1)`;
only its distinguished preimage `beta=0` is subsequently removed.
In particular every `R_beta` is a nonnegative integer, even when the
seed has zero common exponents. Every root of `a` belongs to the root
set of `Q`. Normalization at zero now gives

\[
 Q=z\prod_{\beta\ne0}(1-z/\beta)^{p_\beta},\qquad
 a=z\prod_{\beta\ne0}(1-z/\beta)^{R_\beta}.
 \tag{12}
\]

The degree computation is exact:

\[
                  \deg Q=d(S+1)-(d-1)=dS+1.
 \tag{13}
\]

Hence the number of variables, and the number of distinct monomials in
(5), is `sum_(beta!=0)p_beta=dS`. Repeated coefficients do not merge
terms indexed by different independent variables.

For the number of distinct coefficient values, replace `Q_0` by its
normalized squarefree part

\[
                       \widehat Q_0=z\prod_j(1-c_jz).
\]

The same critical-value condition makes
`widehat Q_0(h)/h'` a polynomial with exactly simple zeros, by (9).
Its root set is the root set of `Q`; its degree is `dt+1`. Removing
zero gives exactly `dt` different reciprocal coefficients.

If the seed coefficients are all distinct, then every `p_j=1`, and
(9) shows that `Q` is squarefree as well. Each target coefficient
occurs once, and its common exponent is exactly the local degree times
the corresponding seed exponent, with (11) on the extra zero fiber.
If all seed common exponents are positive, all these target exponents
are positive. Thus the construction crosses into supports whose
zero/nonzero incidence columns are all identical.

More generally, when `R_j>=p_j` for every repeated-coefficient group,
the new totals satisfy `R_beta>=p_beta`, so an allocation with every
target common exponent positive exists.

## 4. Preservation of the distinguished scalar observable

For every `m>=0`, apply (7) and (4):

\[
\begin{aligned}
 g_m(F_h)
 &=\operatorname{Res}_{z=0}
       \frac{a_0(h(z))^{-m}h'(z)}{Q_0(h(z))}\,dz\\
 &=\operatorname{Res}_{u=0}
       \frac{a_0(u)^{-m}}{Q_0(u)}\,du
 =g_m(F_0).
\end{aligned}
 \tag{14}
\]

Here the change of variable is formal and preserves the distinguished
point zero because `h=z+O(z^2)`. To check it without an analytic
assumption, expand a Laurent series in monomials `u^k`. For `k!=-1`,
`h^k h'` is the derivative of `h^(k+1)/(k+1)`, whose residue is zero.
For `k=-1`, `h'/h=1/z+O(1)` has residue one. Only finitely many negative
Laurent powers occur in (14), so the argument applies coefficientwise.

Equivalently, the normalized primitive for the target is exactly

\[
 H=H_0\circ h,\qquad
 H'=\frac{a_0(h)h'}{Q_0(h)}=\frac aQ,
 \qquad H\circ a^{-1}=H_0\circ a_0^{-1}=B_0.
 \tag{15}
\]

This does not infer scalar observability from a rank or a Gram matrix.
It transports the actual scalar residue at zero by an explicit formal
coordinate substitution. The factors removed in (6) are nonzero
integers, so (6) proves equality of the zero sets over `L`.

Although the target coefficients may lie in `L`, its normalized moments
belong to `K` by (14). This fact does not imply that the target polynomial
itself has coefficients in `K`.

## 5. Iteration

Suppose `h_1` is admissible for `Q_0`, giving `Q_1`, and `h_2` is
admissible for `Q_1`. Then

\[
 \frac{Q_1(h_2)}{h_2'}
 =\frac{Q_0(h_1(h_2))}{h_1'(h_2)h_2'}
 =\frac{Q_0(h_1\circ h_2)}{(h_1\circ h_2)'}.
 \tag{16}
\]

The left side is polynomial, so the composition is admissible for the
original frame. Also `a_2=a_0 composed with h_1 composed with h_2`.
Thus the entire construction, including moment equality, is the single
pullback by `h_1 composed with h_2`. Its degree and the support multiplier
are `deg(h_1)deg(h_2)`. The condition on the second map is relative to
the new root set; arbitrary two maps admissible for the old root set
are not silently declared composable.

## 6. An explicit family with every moment nonzero

Start from the one-term seed `F_0=cY^(r+1)`, where `c!=0` and `r>=0`.
For an integer `d>=1` take

\[
 h(z)=\frac1c\left(1-(1-cz/d)^d\right).
 \tag{17}
\]

Its sole critical value for `d>1` is `1/c`, which is an allowed root
of `Q_0=z(1-cz)`. Let `zeta_d` be a primitive `d`th root of unity in
an extension of `K`. The target phase can be written explicitly as

\[
 F_h=Z_*^{dr}\prod_{j=1}^{d-1}Z_j
       \left(\frac cd Z_*+
          \sum_{j=1}^{d-1}\frac{c}{d(1-\zeta_d^j)}Z_j\right).
 \tag{18}
\]

The preimage `d/c` of `1/c` has local degree `d`; the other nonzero
preimages of zero are `(d/c)(1-zeta_d^j)`, all simple. These prove
the coefficients and exponents in (18), including their distinctness.
The seed has

\[
                  g_m(F_0)=\binom{(r+1)m}{m}c^m.
\]

Therefore the exact ordinary moment formula is

\[
 \boxed{\Lambda(F_h^m)
   =(m!)^d(drm)!\binom{(r+1)m}{m}c^m\ne0
       \quad(m\ge0).}
 \tag{19}
\]

For `d>1,r>=1`, every variable belongs to the common monomial and its
exponents are unequal: one is `dr`, the remaining ones are one. The
coefficient pattern is prescribed by (18), not arbitrary. This
unbounded-support result follows directly from a one-term seed and the
formal transfer, without any analytic rigidity input. The case
`r=1,c=d` is the proposed example
`H(z)=(1-(1-z)^d)/d`, `J(u)=u(1-du)`, and gives
`g_m=binom(2m,m)d^m`.

For a nontrivial preserved zero set, start instead with
`F_0=XY(X-Y)`. Its normalized sequence is zero at odd orders and equals
`binom(3k,k)` at order `2k`, by (7). Every admissible pullback preserves
that exact parity pattern, including unilateral ramification above one
of the coefficient points. This observation uses the same identity,
not an additional symmetry or positivity argument.

## 7. The squarefree outer-composition class

Let `J(u)` be a polynomial with `J(0)=0`, `J'(0)=1`, degree `k>=2`,
and simple roots. Its nonzero reciprocal roots define a seed with
`p_j=R_j=1`, hence `a_0=Q_0=J` and `H_0(u)=u`. For any normalized
polynomial `H` whose critical values belong to the roots of `J`, the
pullback is

\[
       a=J\circ H,\qquad Q=\frac{J\circ H}{H'},
       B=J^{-1},\qquad S_{\rm new}=(\deg H)(k-1).
 \tag{20}
\]

The target common exponents are the local multiplicities of `H` at
the preimages of roots of `J`, hence are all positive. Every proved
moment-window bound for this seed transfers unchanged. In particular
the retained common-exponent monomial-linear theorem supplies its
`k-1` window. The refined bound by the number of nonzero critical
values of `J` uses the separately retained inverse-contact theorem;
it is not needed for the transfer identity proved here.

## 8. An unbounded family with identical zero patterns

Fix integers `q>=2`, `b>=0`, and `d>=2`. Use the cyclic seed

\[
                         P_0=1-z^q,\qquad A_0=P_0^b,
\]

and the normalized polynomial

\[
                         h_d(z)=1-(1-z/d)^d.
\]

Its unique finite critical value is one, which is a root of `P_0`.
The seed coefficient calculation gives, for all `k>=0`,

\[
              g_m=0\quad(q\nmid m),\qquad
              g_{qk}=\binom{(bq+1)k}{k}.
\]

Indeed `(1-z^q)^(-bm-1)` has no other degree classes, and at `m=qk`
its degree-`qk` coefficient is the displayed binomial coefficient.
The residue identity (14) proves the same complete zero set for the lifted phase.

The lifted coefficient polynomial `P` has `dq` distinct roots. At
the point `z=d`, its common exponent is `db`; at the `d(q-1)` points
above the other cyclic roots it is `b`; and at the `d-1` nonzero roots
of `h_d` it is one. Consequently its actual moments are

\[
 \Lambda(F^m)=
 \begin{cases}
 0,&q\nmid m,\\
 (m!)^d\,((bm)!)^{d(q-1)}\,(dbm)!
       \displaystyle\binom{(bq+1)k}{k},&m=qk.
 \end{cases}
\]

In particular every `q`-window is nonzero, while the actual support
has `dq` terms. When `b>=1`, every common exponent is positive:
all monomial columns have the same full zero/nonzero pattern. Thus
this family crosses the zero-pattern limitation of the separating-row
estimate. Its coefficients are constrained by the covering polynomial;
it is not an arbitrary-coefficient theorem for that exponent matrix.

## 9. A four-term example with exact cancellations

Take `q=d=2,b=1`, so `h=z-z^2/4`. Then

\[
 \begin{aligned}
 P&=(1-z/4)(1-z/2)(1+z-z^2/4),\\
 A&=(1-z/4)(1-z/2)^2(1+z-z^2/4).
 \end{aligned}
\]

The actual phase over `Q(sqrt(2))` is

\[
 F=XY^2UV\left(\frac X4+\frac Y2
              +\frac{-1+\sqrt2}{2}U
              +\frac{-1-\sqrt2}{2}V\right).
\]

Every odd moment is zero, whereas

\[
 \Lambda(F^{2k})=((2k)!)^4(4k)!\binom{3k}{k}\ne0.
\]

For example the first two moments are zero and `1152`, respectively.
The checker directly expands the actual four-term phase in its
quadratic coefficient field, independently of the residue calculation.

For an example that necessarily has nonreal coefficients, take
`q=3,d=2,b=1` and the same `h`. Now

\[
 P=(1-z/4)(1-z/2)(1+h+h^2),\qquad A=P(1-z/2).
\]

The four roots `beta_1,...,beta_4` of `1+h+h^2` are all nonreal,
because `4(1+h+h^2)=(2h+1)^2+3` is strictly positive on the real
axis. They are distinct by the pullback multiplicity calculation.
The six-term phase is

\[
 XY^2\prod_{j=1}^4Z_j
       \left(\frac X4+\frac Y2+\sum_{j=1}^4\frac{Z_j}{\beta_j}\right).
\]

Its complete normalized zero set is `3` not dividing the order;
in particular `g_1=g_2=0` and `g_3=4`. Its first three ordinary
moments are `0,0,134369280`. This example has identical full zero
patterns and nonreal coefficient cancellation, with a proved window
of length three despite its six monomials.

## 10. An integer-polynomial certificate with two critical values

Let

\[
 h=z-3z^2+2z^3,
 \quad D=1-6z+6z^2=h',
 \quad C=1+12z-12z^2.
\]

The complete polynomial identity

\[
                           1-108h^2=D^2C
\]

certifies that both finite critical values of `h` are roots of the
seed polynomial `P_0=1-108z^2`. Use `A_0=P_0^b`, `b>=1`. The pullback formulas (4) and (12)
then become

\[
 P=(1-z)(1-2z)DC,\qquad
 A=(1-z)(1-2z)D^{2b}C^b.
\]

The six distinct coefficients and common exponents are

\[
 \begin{array}{c|rrrrrr}
 c&1&2&3+\sqrt3&3-\sqrt3&-6+4\sqrt3&-6-4\sqrt3\\
 r&1&1&2b&2b&b&b.
 \end{array}
\]

This is another full-zero-pattern family with genuine cancellations.
Its normalized moments vanish exactly at odd orders, and

\[
                        g_{2k}=108^k\binom{(2b+1)k}{k}.
\]

The ordinary moment is obtained by multiplying by
`(m!)^3 ((2bm)!)^2 ((bm)!)^2`. The checker verifies the displayed cubic identity as a complete
integer-polynomial identity and directly expands actual cases in
`Q(sqrt(3))`.

## 11. Plane trees as ramification data

For the one-term seed `F_0=cY^(r+1)`, the condition on `h` is precisely
that all finite critical values belong to `{0,1/c}`. Over the complex
numbers, such a polynomial is a Shabat polynomial. Its preimage of the
segment joining these two values is a bicolored plane tree; vertices are
preimages of the two values and their valencies are the corresponding
local degrees. Adrianov and Zvonkin, *Composition of plane trees*,
Definition 1.3 and Theorem 1.5, give this correspondence in the
[author manuscript, actual PDF page 2](https://www.labri.fr/perso/zvonkin/Research/adrzvon.pdf).
That page and the accompanying explanation were inspected directly.

Here is the exact conditional construction, without requiring an
existence theorem for an arbitrary prescribed tree. Given such a
polynomial and a simple preimage of zero, an affine change of input
moves that point to zero and makes `h'(0)=1`. Equivalently, choose a
leaf over zero. If a leaf initially belongs to the other color, swapping
the two output values puts it in the zero fiber. For degree `d`, the
tree has `d+1` vertices. Discard the distinguished leaf, and attach one
variable to each of the `d` remaining vertices. For its normalized
position `beta_v` and valency `e_v`, put

\[
 c_v=1/\beta_v,\qquad
 r_v=\begin{cases}
 e_v,&h(\beta_v)=0,\\
 r e_v,&h(\beta_v)=1/c.
 \end{cases}
\]

Then the actual `d`-term phase

\[
 F=\left(\prod_v Z_v^{r_v}\right)\sum_v c_vZ_v
\]

has the complete moment formula

\[
 \Lambda(F^m)=m!\prod_v(mr_v)!
                \binom{(r+1)m}{m}c^m\ne0\qquad(m\ge0).
\]

This follows from the local multiplicities and residue substitution,
not from a positivity argument on the complex vertex positions.
The valencies on each color sum to `d`, so
`sum_v r_v=d(r+1)-1`, consistent with the degree invariant below.
The field must contain the chosen map and its vertex positions; no
claim that every base field already contains a prescribed tree's
coefficients is made. The explicit maps of Sections 6 and 8 require no
imported tree-existence theorem.

This is a tree encoding of admissible polynomial covers. It does not
supply an expansion or a cancellation theorem for an arbitrary given
factorial phase.

## 12. Every admissible pullback has a degree obstruction

In this section write `S` for the target's actual support size and
`R=sum_i r_i` for its total common exponent over active variables.
Write `S_0,R_0` for the analogous seed quantities. Repeated coefficient
values are counted with their variable multiplicities.

**Proposition 1.** Suppose a normalized polynomial `G` relates two
canonical frames by `a=a_0(G)` and `Q=Q_0(G)/G'`, with `d=deg G`.
Then

\[
              S=dS_0,\qquad R+1=d(R_0+1).
\]

Consequently

\[
                         d\mid\gcd(S,R+1).
\]

**Proof.** Degree multiplication in the identity `a=a_0(G)` gives the
identity `R+1=d(R_0+1)`. Since `Q_0` has degree `S_0+1` and `G'`
has degree `d-1` in characteristic zero, the identity `Q=Q_0(G)/G'` gives
`deg Q=d(S_0+1)-(d-1)=dS_0+1`, proving the first identity. QED.

This proof does not assume that the coefficients are distinct or that
the common exponents are positive. Repeated coefficients count with their
actual variable multiplicities in `P`; dropping that multiplicity would
change the ordinary factorial observable. Thus neither repeated
coefficients nor zero exponents evade the degree divisibility within
this construction.
The ratio `(R+1)/S` is an exact invariant of the pullback. Moreover, the
direction of a window transfer matters: equality of normalized moments
transports a known seed window `S_0` upwards, but knowing only the target's
ordinary SFC window `dS_0` does not produce the smaller seed window.

## 13. Positive exponents: common right factors descend within the class

Now suppose all active `r_i>=1` and all `c_i` are distinct. Then `H` is
a polynomial and

\[
                 \deg a=R+1,\qquad \deg H=R-S+1,
                 \qquad Q=zP\text{ is squarefree}.
\]

**Proposition 2.** Suppose `a=a_0(G)` and `H=H_0(G)` are common right
polynomial decompositions. Normalize their factors at zero so that
`G(0)=a_0(0)=H_0(0)=0` and all three derivatives at zero are one. Then
`Q_0=a_0/H_0'` is a normalized squarefree polynomial, every critical value
of `G` is a root of `Q_0`, and `(a_0,H_0,Q_0)` is again positive-exponent
monomial-linear data. In particular this is exactly an admissible
pullback of the type (4).

**Proof.** The chain rule gives `Q=Q_0(G)/G'` as an identity of rational
functions. At a finite value `u`, let `v=ord_u Q_0`, and take any
preimage `z` of local degree `e`. Then

\[
                         \operatorname{ord}_z Q=e v-(e-1).
\]

Work over an algebraic closure for this local calculation. Since `Q` is
polynomial, `v<0` is impossible. Thus `Q_0` has no finite poles and is
polynomial. If `v=0`, the local-order formula requires `e=1`; this places every
critical value of `G` among its roots. If `v>=1`, squarefreeness of `Q`
forces `v=1`, proving squarefreeness of `Q_0`. Normalization gives its
simple zero and derivative one at zero.

Every root of `a` is a root of `Q`, so every root of `a_0` is a root of
`Q_0`. Conversely a root of `Q_0` not shared by `a_0` would be a pole of
`H_0'=a_0/Q_0`, impossible for a polynomial `H_0`. Thus their root sets
agree, and the multiplicities of the roots of `a_0` give positive integer
common exponents. The identity `H_0'=a_0/Q_0` then gives the canonical
primitive, including its normalization. QED.

This descent statement does not assert the existence of a nontrivial
factor `G`. It classifies a common right factor if one exists.

The same conclusion, allowing a grouped seed, holds whenever every
coefficient group satisfies `R_j>=p_j`, which is exactly the condition
that `H` is polynomial. Distinctness of the individual coefficients is
unnecessary for this extension. The local-order formula still rules
out finite poles of `Q_0` and critical values of `G` outside its roots.
Every root of `a_0` is a root of `Q_0` because this is true after
composition for `a,Q`; conversely polynomiality of `H_0'=a_0/Q_0`
forces every root of `Q_0` into `a_0`. If their multiplicities at such
a root are `p_0,R_0`, then `R_0>=p_0`. Allocating `R_0` among `p_0`
positive exponents constructs a canonical seed. At the distinguished
zero the normalization gives a simple root in both. Thus every existing
polynomial common right factor descends within the grouped polynomial-
observable class as well, even if the original target allocation has
some zero individual exponents. It still does not assert that a
nontrivial factor exists.

### Rational right factors give no escape for a polynomial pair

Suppose instead `a=U(G)` and `H=V(G)` with rational functions. Since `a,H`
are polynomials, the preimages under `G` of all poles of `U,V` must lie at
infinity. There can be only one such pole value: a rational map takes
infinity to one point, and every value has a preimage. Send this pole value
to infinity by a Möbius transformation. The transformed `G` then has no
finite poles and is polynomial; both transformed outer functions are
polynomials as well. Their degrees are unchanged. Thus any nontrivial
rational common right factor still has degree dividing
`gcd(deg a,deg H)=gcd(R+1,S)`.

## 14. A coefficient-uniform primitive obstruction

For `r=(1,1,2)` and three arbitrary nonzero coefficients,

\[
 a=z(1-c_1z)(1-c_2z)(1-c_3z)^2,
 \qquad H=z-\frac{c_3}{2}z^2.
\]

Their degrees are five and two. Hence they have no common right factor
of degree greater than one. This remains true when some coefficient
values coincide: all three variables remain active, and the leading
coefficients in the displayed ternary pair remain nonzero. Equivalently,
over the coefficient field,
`K(a,H)=K(z)`, and the algebraic function `H(a^{-1})` has degree five.

Adding inactive variables leaves `a,H,P` unchanged after the nonzero
moment factors are removed. Deleting an active term changes the actual
support and is a different problem. Splitting or identifying coefficient
groups while dropping their multiplicities changes `P` and therefore
changes the actual residue observable. None of these is an unmentioned
reduction of the displayed ternary phase.

**Limit of this obstruction.** It does not refute SFC or exclude another
inverse identity. In fact every equal-exponent family `r_i=q` also has
`gcd(qS+1,S)=1`, yet the retained sparse inverse bridge already proves
its SFC statement. The degree test rejects this common-right-factor
reduction only; calling every degree-primitive case unproved would be
incorrect.

## 15. The shared derivative alone is insufficient

The zero-fiber condition cannot be discarded from a sharper contact
argument. Set

\[
 H=z-z^2/2,\quad
 D=1+3z/5+6z^2/25,\quad
 a=z-z^2/5-3z^3/25-3z^4/50,\quad
 b(u)=u-3u^2/10.
\]

These are normalized polynomials, `a'=H'D`, `deg D=s=2`, and
`deg b=n=2`. Direct exact expansion gives

\[
               H-b(a)=-\frac{27}{1250}z^5+O(z^6).
\]

Thus `ord_0(H-b(a))=5>n+s=4`. The shared critical point of `H` is
`z=1`, but `a(1)=31/50`, not zero. This is outside the special
integer-root-profile family and is not an SFC counterexample. It
refutes an argument retaining only the low degree of `D=a'/H'`.

## 16. Exact boundary controls and a recovered descent

The coefficient map in this construction is essential. A target has
coefficients `1/beta` determined by the preimages, rather than independent
free coefficients at all its terms. On a family with fixed multiplicities,
after substituting this coefficient map, the ideals of any chosen moment
window agree with those of the seed up to invertible degreewise factorial
scalars. This is not an assertion about full ambient Macaulay matrices
with freely chosen target coefficients.

Three exact controls guard the endpoint hypotheses.

- **Zero common exponent.** Seed values `(2,-3)`, multiplicities `(1,1)`,
  aggregate exponents `(0,2)`, and `h=z-z^2/2` give
  `P=(1-z/2)(1-z)(1+3z-3z^2/2)` and
  `A=(1-z/2)(1+3z-3z^2/2)^2`. The point `z=1` is still active,
  with common exponent zero. Removing it changes the first moment.
- **Repeated coefficient.** With the same values and cover, but
  `p=(2,1),R=(1,3)`, the point `z=1` carries three variables of
  coefficient one and total common exponent two. Both allocations
  `(0,0,2)` and `(1,1,0)` preserve the normalized sequence. Its correct
  multiplicity is three, not one or four.
- **Invalid cover.** For `Q_0=z(1-z)` and `h=z+z^2`, division of
  `Q_0(h)` by `h'=1+2z` leaves remainder `-5/16`. Its critical value
  `-1/4` is outside `{0,1}`. Thus `Q_0(h)/h'` has a finite pole and
  does not define the asserted ordinary frame. Equivalently the
  transported differential has an extra finite zero. Formal
  invertibility at zero alone does not establish admissibility.

One previously observed residue-rank degeneration has a genuine descent.
At `t=2` in `a=z(1-z)(1-tz)^2`, put

\[
 h=z-z^2,\qquad a_0=Q_0=u(1-4u),\qquad H_0=u.
\]

Then

\[
 a=a_0(h)=z(1-z)(1-2z)^2,\quad
 Q=Q_0(h)/h'=z(1-z)(1-2z),\quad H=h.
\]

The seed has one term and the target has two, with common exponents
`(1,2)` and coefficients `(1,2)`. Its normalized moments are
`4^m binom(2m,m)`, all nonzero. This is compatible with the persistent
`t-2` factor in the retained residue-observability determinants: that
higher-dimensional vector model is redundant at this specialization.
No scalar vanishing follows, and isolated rank degenerations at other
values do not establish a descent.

More explicitly, the source differential has residues `g_m,-g_m` at
`u=0,1/4`. Under the cover, the residues at `z=0,1,1/2` are
`g_m,g_m,-2g_m`, respectively, because residue pullback multiplies by
the local degree. Thus the two-coordinate residue vector in the retained
example is `(g_m,g_m)` at every order. This proves dependence of all
its two-by-two determinants without implying a zero factorial moment.

The new theorem compares a seed with its constructed phase. It does not
compare an arbitrary linear form with an arbitrary common-monomial
multiple. The primitive obstruction and the contact counterexample are
limits of proposed mechanisms, not counterexamples to SFC.

## 17. Primary-source comparison and remaining obligation

The factorial coefficient observable agrees with the existing
[monomial-linear contact note](SFC_MONOMIAL_LINEAR_CONTACT_BOUND.md).
The proof of the general transfer uses only the multinomial expansion,
local polynomial multiplicities and formal residues. Existing all-order
seeds can additionally come from the
[common-exponent theorem](POWER_DIAGONAL_LINEAR_SFC_TRANSFER.md) and
[weighted radial theorem](WEIGHTED_RADIAL_LINEAR_SFC_TRANSFER.md), with
all their restrictions retained.

Edo and van den Essen, *The Strong Factorial Conjecture*,
[arXiv:1304.3956v2](https://arxiv.org/pdf/1304.3956v2), Theorem 2.25
and Lemma 2.26, supply an earlier inverse-series bridge for a specified
monomial-linear class. The statements and displayed proof were inspected.
The current general residue identity and its ramification hypotheses
are proved explicitly above; the source is not cited as an
arbitrary-common-exponent SFC theorem.

Pakovich, *On the Polynomial Moment Problem*, Mathematical Research
Letters **10** (2003), 401--410,
[author PDF](https://pakovich.github.io/opmp.pdf), Section 2.1, Lemma 1,
proves the algebraic branch-degree identity

\[
 \deg_{\rm alg}\bigl(H(a^{-1})\bigr)
   =[\mathbb C(a,H):\mathbb C(a)]
   =\frac{\deg a}{[\mathbb C(z):\mathbb C(a,H)]}.
\]

The actual proof is the tower of field degrees. Corollary 1 relates a
common right factor to equality between two inverse-branch observables.
These statements fit the descent classification, but neither they nor
that paper's Theorem 1 imply composition from a finite local Taylor gap.
The latter theorem concerns every integral moment with endpoint
conditions. A finite-gap-to-composition implication has not been proved
here. No literature-priority claim follows from the search.

After removing any existing common right factors in the positive,
distinct-coefficient case, the remaining sharp contact problem is still

\[
 \operatorname{ord}_0(H-b\circ a)\le\deg b+S.
\]

It must use the prescribed zero-fiber and integer multiplicity data.
The primitive `(1,1,2)` profile is a concrete boundary for the present
reduction, while the counterexample in Section 15 shows why keeping
only `a'=H'D` cannot prove the desired inequality. A window known for
an enlarged phase also cannot automatically be shortened for its seed.
Unrestricted SFC, unrestricted common-factor multiplication and the
sharp general contact bound remain open in this investigation.

## 18. Exact replay and assurance boundary

Run from the repository root:

```sh
python research/factorial/image-rigidity/verify_sfc_ramified_pullback.py
```

The checker uses only standard-library rational polynomial arithmetic
and exact quadratic-field pairs. It authenticates this note and the
retained proof dependencies before writing
[sfc_ramified_pullback_verification.json](sfc_ramified_pullback_verification.json).
The full source is
[verify_sfc_ramified_pullback.py](verify_sfc_ramified_pullback.py).

The finite certificate retains normalized moments through order 12 for
45 polynomial pullbacks; 51 independently expanded ordinary factorial
moments, including order-zero controls, over `Q(sqrt(2))`, `Q(sqrt(3))` and `Q(sqrt(15))`; the integer
cubic identity; exact complex-cancellation data; composition consistency;
zero-exponent and repeated-coefficient controls; and rejection of invalid
covers. It also records the complete weakened-contact counterexample,
the primitive-degree controls, and the recovered `t=2` descent.
Optimized execution and changed proof bytes must be rejected.

The universal statements rest on the written algebraic proofs. Finite
examples check formulas and boundary cases; they do not prove the
unbounded quantifiers. The arguments and arithmetic received independent
internal reconstructions by other agents with workspace context, which
are not external referee review or formal verification. No RH or Alpoge
claim, or positivity on the coefficient field, is used. This is a local
research addition; no publication or remote repository update is made.
