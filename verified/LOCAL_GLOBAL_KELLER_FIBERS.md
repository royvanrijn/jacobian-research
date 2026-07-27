# Local-to-global synthesis for Keller fibers

The fixed-seed [adelic fiber theorem](ADELIC_FIBER_ENGINEERING.md) moves a
target in a two-parameter pencil and prescribes real sheets together with
unramified factorization types.  The present theorem allows the Keller map
itself to depend on the answer.  This extra freedom permits arbitrary finite
étale algebras over finitely many completions, including ramified algebras.

## 1. Compatibility data

Fix an integer `N>=3`.  Prescribe:

1. for every prime `p` in a finite set `S`, a rank-`N` finite étale
   `Q_p`-algebra `A_p`;
2. a signature `(r_1,r_2)` with `r_1+2r_2=N`; and
3. for every prime `ell` in a finite set `R`, a partition
   `lambda_ell` of `N`.

The compatibility conditions are exactly the following.

* All local ranks and the archimedean rank equal `N`.
* For each `ell`, there is a monic squarefree polynomial over `F_ell` with
  factor-degree partition `lambda_ell`.
* If `ell` lies in both `R` and `S`, then `A_ell` is unramified and its
  residue degrees are `lambda_ell`.

The second condition is automatic for most primes but is necessary at small
primes: for example, `F_2` does not contain arbitrarily many distinct linear
factors.  There is no Grunwald--Wang condition because no common Galois group
or global Galois structure is prescribed.

## 2. Local-to-global polynomial theorem

> **Theorem (finite local specifications).**  Subject to the compatibility
> conditions above, there is a monic squarefree polynomial
>
> \[
>                         P(T)\in\mathbb Q[T],\qquad \deg P=N,       \tag{1}
> \]
>
> such that
>
> \[
> \mathbb Q[T]/(P)\otimes_{\mathbb Q}\mathbb Q_p\simeq A_p
> \quad(p\in S),                                                   \tag{2}
> \]
>
> `P` has exactly `r_1` real roots, and its squarefree reduction at every
> `ell in R` has factor-degree partition `lambda_ell`.
>
> The polynomial may additionally be required to be irreducible over `Q`.
> Thus it may be chosen to define a degree-`N` number field having all the
> prescribed completions, signature, and unramified splitting types.

### Proof

Every finite étale algebra over the infinite field `Q_p` is monogenic.
Choose a monic squarefree presentation

\[
                         A_p\simeq\mathbb Q_p[T]/(P_p).           \tag{3}
\]

Factor stability and Krasner's lemma give a coefficient neighborhood
`Omega_p` of `P_p` with the property that every monic polynomial in
`Omega_p` defines the same `Q_p`-algebra.  More explicitly, first separate
the pairwise coprime factors of `P_p` by Hensel factorization.  Around each
irreducible factor choose a Krasner radius for one primitive root.  A
sufficiently small coefficient ball preserves the factor degrees and places
a root inside every selected Krasner ball.  Products of repeated isomorphic
field factors cause no problem: monogenicity lets their primitive elements
be chosen with distinct minimal polynomials.

At the real place, choose `r_1` distinct real numbers and `r_2` distinct
conjugate pairs of nonreal numbers.  The monic polynomial with these roots
has a coefficient neighborhood `Omega_infty` contained in the discriminant
complement and on which the real-root count is `r_1`.

For `ell in R-S`, choose a monic squarefree witness
`bar P_ell in F_ell[T]` of type `lambda_ell`.  The coefficient congruence

\[
                         P\equiv\bar P_\ell\pmod\ell              \tag{4}
\]

is an `ell`-adic open condition.  At a prime in `R intersect S`, choose the
presentation (3) integral with this reduction; the compatibility hypothesis
makes the two conditions the same after shrinking.

Weak approximation on the affine space of the `N` lower coefficients gives
a rational monic `P` in

\[
              \Omega_\infty\times\prod_{p\in S}\Omega_p
                    \times\prod_{\ell\in R-S}\Omega_\ell.        \tag{5}
\]

It is squarefree because `Omega_infty` was chosen inside the discriminant
complement.  Conditions (2), the signature, and all reductions follow from
the definitions of the neighborhoods.

For connectedness, choose a new prime `q` outside `S union R` and a monic
irreducible polynomial of degree `N` over `F_q`, and add its residue class
to (5).  The resulting rational polynomial is irreducible modulo `q`, hence
irreducible over `Q`.  This proves the last assertion without Hilbert
irreducibility.  Alternatively, Hilbert irreducibility with weak
approximation selects an irreducible point inside (5).  The auxiliary-prime
argument is more elementary and constructive.  \(\square\)

### A universal automatic stability radius

The Krasner neighborhoods in the proof can be made explicit without
factoring over a splitting field.  Let `f in Z_p[T]` be monic and squarefree,
and put

\[
                         D=v_p(\operatorname{Disc}(f)).          \tag{5a}
\]

Then every monic same-degree polynomial `g` satisfying

\[
                         g\equiv f\pmod {p^{\,2D+1}}             \tag{5b}
\]

defines the same finite étale `Q_p`-algebra:

\[
                         \mathbb Q_p[T]/(g)\simeq
                         \mathbb Q_p[T]/(f).                     \tag{5c}
\]

Indeed, work in a finite splitting field and let `alpha_i` be the roots of
`f`.  They are integral.  Write

\[
 d_i=v_p(f'(\alpha_i))
     =\sum_{j\ne i}v_p(\alpha_i-\alpha_j).
\]

All summands are nonnegative and

\[
                         \sum_i d_i=D.                           \tag{5d}
\]

For `m=2D+1`, coefficient congruence gives
`v_p(g(alpha_i))>=m`.  It also gives
`v_p(g'(alpha_i)-f'(alpha_i))>=m>d_i`, hence
`v_p(g'(alpha_i))=d_i`.  Strong Hensel produces a unique root `beta_i` of
`g` satisfying

\[
 v_p(\beta_i-\alpha_i)\ge m-d_i>d_i
       \ge\max_{j\ne i}v_p(\alpha_i-\alpha_j).                   \tag{5e}
\]

The root balls are disjoint, so the `beta_i` are all distinct.  Uniqueness
makes the matching `alpha_i -> beta_i` Galois-equivariant.  The two root
Galois sets, and therefore the corresponding finite étale `Q_p`-algebras,
are isomorphic.  This proves (5c).

The exponent `2D+1` is uniform and often far from minimal.  A factorwise
Hensel--Krasner calculation can certify a smaller ball, as in the optimized
quintic below.  The automatic bound is valuable because it requires only an
exact discriminant valuation.  Every finite étale `Q_p`-algebra admits a
monic integral presentation after choosing an integral primitive generator,
so the bound applies without an unramifiedness hypothesis.

## 3. Explicit coefficient CRT

The proof becomes an algorithm once each `Omega_p` is supplied as a certified
coefficient ball.  After rational truncation, write its coordinate
conditions as

\[
                c_i\equiv c_{i,p}\pmod {p^{m_p}}\qquad(0\le i<N).
                                                                    \tag{6}
\]

Choose a positive integer `D_0` clearing every finite center and put

\[
 e_p=m_p+v_p(D_0),\qquad M=\prod_p p^{e_p}.
\]

For any `u=1+kM`, put `D=D_0u`.  Coefficientwise CRT gives residues

\[
                 a_i\equiv D c_{i,p}\pmod {p^{e_p}},\qquad
                 c_i={a_i\over D}.                               \tag{7}
\]

For fixed `D`, the choices for `a_i` form grids of mesh `M/D`.  Increasing
`k` makes the grids meet any prescribed rational box inside
`Omega_infty`, while `u=1 mod M` preserves every finite condition.  This is
a terminating coefficient construction, not merely an existence appeal.
Krasner or Hensel inequalities certify how large the exponents in (6) must
be.

The construction is implemented in
`jcsearch.local_global.synthesize_monic_polynomial`.  Its input is a
dictionary

```python
prime: (certified_precision, monic_local_model)
```

together with rational open intervals for the lower coefficients.  It
returns the polynomial and a certificate recording `D_0`, `M`, `k`, and the
common denominator `D`.  Local centers may have denominators divisible by
their primes: the implementation includes `v_p(D_0)` in the CRT modulus, as
in (7).  A precision value of `None` requests the automatic radius
`2v_p(Disc)+1`; an explicit integer uses a sharper independently certified
radius.

## 4. Compilation into a complete Keller fiber

Choose `a in Q` outside the finite zero set of `P'P'''` and put

\[
 G(S)=P(a+S)-P(a)=g_1S+\cdots+g_NS^N.                \tag{8}
\]

Then `g_1g_3g_N` is nonzero.  The root-engineered quadratic-gauge compiler
gives a polynomial map

\[
 \widetilde F_G:\mathbb A^3_{\mathbb Q}\longrightarrow
                 \mathbb A^3_{\mathbb Q}
\]

with

\[
 \det D\widetilde F_G=1,\qquad
 \operatorname{gdeg}(\widetilde F_G)=N,\qquad
 \max_i\deg(\widetilde F_G)_i\le6N+2.               \tag{9}
\]

At the rational target

\[
                    y_{P,a}=\left(1,0,-{2P(a)\over P'(a)}\right) \tag{10}
\]

there is a scheme-theoretic identity

\[
 \boxed{\;
\widetilde F_G^{-1}(y_{P,a})
    \simeq\operatorname{Spec}\mathbb Q[T]/(P).
\;}                                                               \tag{11}
\]

The compiler is exposed as
`jcsearch.keller_fiber.compile_polynomial_to_keller_fiber`.  It checks
squarefreeness and degree, chooses an admissible rational translation when
one is not supplied, and returns the seed, determinant-minus-two map,
determinant-one normalization, target, inverse polynomial, geometric degree,
and exact coordinate degrees.  Together with
`synthesize_monic_polynomial`, this gives an end-to-end exact API from
certified local coefficient balls to the displayed Keller fiber.

The fiber has rank equal to the geometric degree, so it is complete.  Base
change of (11) to `R` or `Q_p` transports all the prescribed local data.
Combining the polynomial theorem with (11) proves:

> **Local-to-global Keller-fiber theorem.**  Every compatible finite
> collection of rank-`N` local étale specifications over `Q`, for `N>=3`,
> occurs in a complete fiber of an explicit Jacobian-one polynomial map
> `A^3_Q -> A^3_Q`.  The fiber may be required to be connected.

This is stronger arithmetically than realizing an already constructed
abstract number field.  It first synthesizes one global field or étale
algebra from independent local specifications and then compiles that exact
algebra into a full Keller fiber.  The price is that the map is tailored to
the synthesized polynomial; the theorem does not put every such algebra
inside one fixed map.

## 5. Adelic data in infinitely many stable classes

Combine the local-to-global theorem with the
[universal Keller-fiber multiplicity theorem](UNIVERSAL_KELLER_FIBER_MULTIPLICITY.md).
For `N>=4`, first choose a connected algebra

\[
                         A=\mathbb Q[T]/(P)                         \tag{12}
\]

realizing the prescribed local package.  Universal multiplicity gives
infinitely many stable polynomial left--right classes of determinant-one
Keller maps, every one having `Spec A` as a complete fiber.  Base change of
that same fiber preserves all its completions, signature, and Frobenius
types.  Therefore:

> **Adelic stable-multiplicity corollary.**  For every `N>=4`, every
> compatible finite collection of local étale specifications occurs in one
> connected complete fiber shared, fiberwise, by infinitely many stable
> polynomial left--right classes of determinant-one Keller maps over `Q`.

The maps and marked targets vary, but the finite étale fiber algebra—and
hence the entire prescribed arithmetic package—is fixed.

If the maps themselves must remain fixed, the compatible local packages in
the one-parameter common-fiber pencil are treated by
[locally prescribed common fibers](LOCALLY_PRESCRIBED_COMMON_FIBERS.md).

## 6. A ramified quintic certificate

Let `U_(p,d)` denote the unique unramified extension of `Q_p` of degree `d`.
Take the local algebras

\[
\begin{aligned}
 A_2&=\mathbb Q_2(\sqrt[3]{2})\times U_{2,2},\\
 A_3&=\mathbb Q_3(\sqrt3)\times U_{3,3},
 \end{aligned}                                                     \tag{13}
\]

prescribe signature `(1,2)`, an inert prime at `5`, and splitting type
`(2,2,1)` at `7`.

### Fully automatic discriminant-radius lift

Use the local models

\[
\begin{aligned}
 f_2&=(T^3-2)(T^2+T+1),&
 v_2(\operatorname{Disc}(f_2))&=2,\\
 f_3&=(T^2-3)(T^3-T+1),&
 v_3(\operatorname{Disc}(f_3))&=1.
\end{aligned}
\]

The automatic theorem selects precisions `5` at `2` and `3` at `3`.
The unramified models at `5` and `7` have discriminant valuation zero, so
their automatic precision is one.  Thus the coefficient CRT modulus is

\[
                         2^5\,3^3\,5\,7=30240.
\]

Using the same real coefficient box as the optimized lift below, the generic
constructor returns

\[
\boxed{
 Q(T)=T^5-{9855\over30241}T^4+{163265\over30241}T^3
       +{190\over30241}T^2+{113214\over30241}T
       -{7266\over30241}.
}                                                               \tag{13a}
\]

No factorwise local calculation is needed: the four discriminant
certificates prove all local algebras and Frobenius types directly.  Exact
factorization modulo `5` proves global irreducibility, and a Sturm count
gives signature `(1,2)`.  The public compiler takes `a=0` and produces a
determinant-one degree-five map with complete target

\[
                         \left(1,0,{2422\over18869}\right).       \tag{13b}
\]

### Sharper factorwise lift

The factorwise Hensel--Krasner estimates reduce the needed precisions at both
ramified primes to two.  The resulting smaller coefficient CRT lift is

\[
\boxed{
 P(T)=T^5+{225\over1261}T^4+{5765\over1261}T^3
       +{190\over1261}T^2+{4854\over1261}T+{294\over1261}.
}                                                                  \tag{14}
\]

Its primitive integral form is

\[
1261T^5+225T^4+5765T^3+190T^2+4854T+294.            \tag{15}
\]

The exact congruences are

\[
\begin{aligned}
 P&\equiv(T^3-2)(T^2+T+1) &&\pmod4,\\
 P&\equiv(T^2-3)(T^3-T+1) &&\pmod9,\\
 P&\equiv T^5-T-1 &&\pmod5,\\
 P&\equiv T(T^2+1)(T^2+T+3) &&\pmod7.               \tag{16}
\end{aligned}
\]

For the first line, the two displayed factors have unit resultant over
`Z_2`, so Hensel factorization separates degrees three and two.  The
quadratic factor is the unique unramified quadratic extension.  The cubic
factor is coefficientwise `4`-adically close to `T^3-2`.  In
`Q_2(alpha)`, `alpha^3=2`, the perturbation valuation is at least `2`, while
`v_2(3alpha^2)=2/3`; strong Hensel and Krasner give the same ramified cubic
field.

The second line is identical with degrees two and three exchanged.  For
`alpha^2=3`, the perturbation valuation is at least `2` and
`v_3(2alpha)=1/2`, so the ramified factor is `Q_3(sqrt3)`; the other factor
is the unique unramified cubic extension.

The polynomial `T^5-T-1` is an irreducible Artin--Schreier polynomial over
`F_5`, proving that (14) is globally irreducible.  The two quadratics at `7`
are distinct and irreducible.  An exact Sturm count gives one real root, so
the field has signature `(1,2)`.

Here `a=0` is admissible.  Put `t=1+xy` and

\[
 q=t^2z+{4854\over5765}y^2(1+3t),\qquad \Pi=tq.
\]

The determinant-minus-two outputs are

\[
\begin{aligned}
 B={}&y+{5765\over1618}xq+{190\over2427}tq
       +{150\over809}t^2x^2q^4
       +{6305\over4854}t^2x^3q^5,\\
 C={}&x(5-3t)-{5765\over4854}x^3z
       -{75\over809}(xq)^4-{1261\over1618}(xq)^5 .
\end{aligned}                                                     \tag{17}
\]

Thus

\[
                 \widetilde F=(\Pi,-B/2,C)                       \tag{18}
\]

has determinant one and geometric degree five, and

\[
 \boxed{\;
 \widetilde F^{-1}\left(1,0,-{98\over809}\right)
       \simeq\operatorname{Spec}\mathbb Q[T]/(P).
 \;}                                                               \tag{19}
\]

## Verification

Run

```bash
.venv/bin/python scripts/verify_local_global_keller_fibers.py
```

The checker derives the four automatic discriminant radii, reconstructs
(13a), and verifies its local balls, signature, irreducibility, and Keller
target.  It separately reproduces the sharper polynomial (14), feeds it
through the public Keller compiler, and verifies the factorwise
Hensel--Krasner inequalities, determinant-one compilation, inverse
polynomial, and scheme-theoretic quotient reconstruction.  The older
finite-étale verifier also imports the same shared compiler and replays its
degree-three, degree-four, degree-five, and Hasse-fiber regressions.
