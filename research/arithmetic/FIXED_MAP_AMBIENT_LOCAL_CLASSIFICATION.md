# Ambient arithmetic classification for the fixed quintic Keller map

17 November 2026. These are written consequences of the retained quadratic-gauge reconstruction, elementary finite permutation groups, and the classical decomposition-group criterion. The new statements have exact symbolic and finite-group controls below, but no formal verification or external review is claimed. They classify local solubility; they do not supply an ambient counting asymptotic.

The fixed map is the determinant-one map `Phi=L o F0 o A` of [the fixed-map manuscript, Section 2](../../unpublished/fixed-map-hasse-failures/sections/02-fixed-map.tex). Its generic degree is five. We use a normalized target `(u,b,c)` for `Phi` and write `pi=-u`, so `(pi,b,c)` is the target for `F0`. Its inverse polynomial is

\[
 E_{\pi,b,c}(S)=\pi^5S^5-\frac32\pi^4S^4+\frac32\pi S^3
 +\left(-\frac54\pi-\frac9{32}b\right)S^2
 +\frac9{16}S-\frac9{32}c.                 \tag{1}
\]

A finite scheme over `Q` is **Hasse-failing** if it has a point over every completion of `Q` and no `Q`-point. Empty schemes have no local points. All fibres under consideration are entire scheme-theoretic fibres, not selected subcollections of points.

## 0. Lemma FMA0: the low-rank local-global obstruction

Let `X` be a finite étale `Q`-scheme of rank at most four. If it has a point over every completion, it has a rational point. In particular every Hasse-failing finite étale `Q`-scheme has rank at least five.

**Proof.** Rank zero is empty and cannot meet the local hypothesis. A rank-one component is a rational point. Without such a component, the only possible partitions in ranks at most four are a single connected component of degree two, three or four, or two quadratic field components. For a connected component, its transitive finite Galois permutation action contains a derangement. Indeed Burnside's orbit formula gives average fixed-point count one, while the identity fixes more than one point; if every element fixed a point the average would be strictly larger than one. Chebotarev supplies an unramified prime with that Frobenius conjugacy class, so the scheme has no point over the corresponding `Q_p`.

For two quadratic fields which coincide, use the nontrivial element in their common `C_2` quotient. If they are distinct, their compositum is biquadratic and its group contains an element nontrivial in both quotients. In either case Chebotarev supplies a prime where both components lack a point. This exhausts the partitions and proves the assertion. The relation between local roots and decomposition-group fixed points is precisely Sonn's Proposition 2.1 cited below; only the elementary subgroup possibilities are needed here. ∎

For the polynomial formulation, let `f in Z[T]` have degree at most four
and a root modulo every positive integer. For each prime `p`, the nonempty
closed sets `{x in Z_p: f(x) == 0 mod p^j}` are nested and compact, so their
intersection supplies a `Q_p`-root. A polynomial with a root modulo every
positive integer cannot be a nonzero constant; the zero polynomial already
has a rational root. Replace a nonconstant `f` by its squarefree part over
`Q`, which has the same roots over every characteristic-zero field.
The proof just given uses only the finite primes, and therefore supplies a
rational root. This proves the low-degree intersective obstruction directly
from the same Chebotarev interface, without importing an unchecked statement
from Berend--Bilu.

## 1. Lemma FMA1: every possible Hasse failure is a full squarefree quintic fibre

For `pi!=0`, the complete fibre is

\[
 \Phi^{-1}(u,b,c)\simeq
 \operatorname{Spec}\bigl((\mathbb Q[S]/E_{\pi,b,c})[1/E'_{\pi,b,c}]\bigr).
                                                        \tag{2}
\]

Its geometric rank is the number of simple roots of `E`. Thus a nonsquarefree quintic has a fibre of rank at most three. For `pi=0`, the fibre has rank three if `b!=0` and `bc!=1`, and rank one otherwise; every such rational target has a rational source point. Consequently every Hasse-failing fibre of `Phi` has `pi!=0`, squarefree `E`, rank five, and rank equal to the generic map degree.

**Proof.** The source variables are

\[
 t=1+xy,\qquad q=t^2z+\frac38y^2(1+3t),\qquad \Pi=tq.
\]

The derivative-unit reconstruction in the fixed-map manuscript, Proposition `prop:gauge`, proves (2) over every commutative test algebra: when `pi` is a unit, `t` and `q` are units, and `D=(16/9)E'=t^{-1}`. Localization removes exactly the primary factors belonging to repeated roots. If a quintic has a repeated root, at least two units of its degree belong to repeated-root primary factors, leaving at most three simple roots.

For completeness the `pi=0` boundary can be computed directly, without the nonproperness theorem or its analytic hypotheses. The divisors `t=0` and `q=0` do not meet: at `t=0`, `xy=-1` and `q=(3/8)y^2` is invertible. On `q=0`, `B=y`, `S=x/t` is defined, and

\[
 C=2S-BS^2,\qquad 1-BS=t^{-1}.
\]

This component therefore contributes two points if `b!=0,bc!=1`, zero if `b!=0,bc=1`, and one if `b=0`. On `t=0`, `B=-2y=2/x`; this component contributes one rational point if `b!=0` and none if `b=0`. For `b!=0`, take `x=2/b`, `y=-b/2`; the equation for `C` determines `z` over `Q` because its coefficient is `-(8/3)x^3!=0`. For `b=0`, the normalized source point `(c,0,0)` maps under `Phi` to `(0,0,c)`.

Finally Lemma FMA0 excludes every short fibre above. This also recovers the low-rank argument used in [the minimality section](../../unpublished/fixed-map-hasse-failures/sections/05-minimality.tex), independently of any unverified application of the Berend–Bilu criterion. ∎

## 2. Theorem FMA2: classification of rank-five Hasse-failing algebras

Let `A` be a finite étale `Q`-algebra of rank five. The following are equivalent:

1. `Spec A` is Hasse-failing.
2. There are a non-Galois cubic field `L/Q` and its quadratic resolvent `K`, with

   \[
    A\simeq K\times L,
    \qquad K=\mathbb Q(\sqrt{\operatorname{disc}(L)}),       \tag{3}
   \]

   such that every finite decomposition group in the `S_3` normal closure `M/Q` of `L` is proper.
3. There are `L,K` as in (3) such that, at every rational prime `p`, either `K` splits over `Q_p` or `L` has a degree-one factor over `Q_p`.

Every proper subgroup of `S_3` is cyclic, so “proper” can equivalently be replaced by “cyclic” in condition 2. Unramified primes automatically meet the condition. The real completion is automatic because a cubic has a real root.

**Proof.** A rational-point-free rank-five étale algebra has factor-degree partition either `5` or `2+3`. A connected degree-five algebra has a transitive Galois action. Such an action contains a derangement: Burnside's orbit formula gives average number of fixed points one, whereas the identity fixes five points, so some element fixes none. Chebotarev produces an unramified prime without a local point. Thus `A=K x L`, with `K` quadratic and `L` cubic fields.

If `L` is cyclic, its compositum with `K` has Galois group `C_3 x C_2`. An element nontrivial in both factors fixes no point in either component, again contradicting Chebotarev. Suppose therefore that `L` has `S_3` normal closure `M`. Its only quadratic subfield is the discriminant field. If `K` is a different quadratic field, `K cap M=Q` and the compositum has Galois group `S_3 x C_2`; a three-cycle in the first factor and the nonidentity element in the second again give a common derangement. Hence (3) is necessary.

Under (3), the five geometric points have `S_3`-set

\[
 S_3/A_3\ \sqcup\ S_3/C_2.                              \tag{4}
\]

A point is defined over `Q_p` exactly when a decomposition group fixes one of these points, equivalently when that group is contained in `A_3` or a conjugate of `C_2`. These are exactly the proper subgroups of `S_3`. This proves both implications and the equivalence with condition 3. Neither component has a rational point, and the real assertion follows from the odd degree of the cubic. ∎

The local decomposition-group criterion is classical. Its precise primary source is J. Sonn, *Polynomials with roots in Q_p for all p*, Proposition 2.1, pp. 2–3 of [arXiv:math/0612528v4](https://arxiv.org/pdf/math/0612528v4), 28 January 2007; published in *Proc. Amer. Math. Soc.* **136** (2008), 1955–1960, [DOI](https://doi.org/10.1090/S0002-9939-08-09155-7). The source was read directly. The rank-five statement above is its elementary specialization with the subgroup possibilities proved explicitly; no novelty is claimed for the general group criterion.

## 3. Theorem FMA3: exact criterion on the entire rational target of Phi

For any `(u,b,c) in Q^3`, the fibre `Phi^{-1}(u,b,c)` is Hasse-failing if and only if all of the following hold:

- `u!=0` and `E_{-u,b,c}` is squarefree;
- over `Q`, its irreducible factor degrees are exactly `2+3`, so `E=lambda q h`, with monic irreducible quadratic `q` and cubic `h`;
- their discriminants have the same squareclass:

  \[
   \frac{\operatorname{disc}(h)}{\operatorname{disc}(q)}
        \in\mathbb Q^{\times2};                         \tag{5}
  \]

- at every prime ramified in the splitting field of `h`, either `disc(q)` is a square in `Q_p` or `h` has a root in `Q_p`.

The final condition is a finite effective test. In particular, choose an integer `d>0` such that

\[
 H(T)=d^3h(T/d)\in\mathbb Z[T]
\]

is monic. It is enough to test primes dividing the nonzero integer `Disc(H)`. At each such prime `p`, put `v=v_p(Disc(H))`. The assertion that `H` has a `Q_p`-root is equivalent to the existence of a root modulo `p^(2v+1)`. Thus the criterion contains only rational polynomial factorization, rational-square testing, and finitely many finite congruence tests, although the displayed bound is not optimized.

More generally the fibre is everywhere locally soluble if and only if either it has a rational point or it satisfies the Hasse-failure criterion above. For `u=0` it always has a rational point. For `u!=0` it has a rational point exactly when `E_{-u,b,c}` has a **simple rational root**; a repeated rational root removed by (2) does not suffice.

**Proof.** Lemma FMA1 reduces any Hasse failure to a full rank-five algebra. Theorem FMA2 gives precisely the factorization, common quadratic-resolvent condition (5), and local tests. Conversely these conditions give the algebra (3) and local points everywhere, with no rational point.

For the finite prime set, `H` defines the same cubic field as `h`. At primes outside `Disc(H)` its three roots give an unramified étale algebra over `Z_p`; its splitting field is consequently unramified. Such decomposition groups are cyclic, and Theorem FMA2 makes the local condition automatic. Condition (5) also places the quadratic field inside that same splitting field. Extra primes dividing the index of `Z[T]/H` cause harmless extra tests.

Here is the finite root bound. The integral resultant identity gives `U H+V H'=Disc(H)` up to sign, with `U,V in Z[T]`, because `H` is monic. If `H(x)=0 mod p^(2v+1)`, evaluation of this identity shows `v_p(H'(x))<=v`. The strict Hensel inequality `v_p(H(x))>2v_p(H'(x))` gives a `p`-adic root. Conversely a root of a monic integral polynomial is `p`-adically integral and reduces to a root modulo every power of `p`. The rational-point assertion is immediate from (2). ∎

This is an arithmetic decision criterion for every rational target, not a parametrization of all targets and not a counting theorem. The condition that the fibre have a specified factorization and resolvent is essential; the generic quintic covering has symmetric monodromy.

## 4. Proposition FMA4: a second completely classified target line

Put

\[
 \gamma=\frac3{32},\qquad k=\gamma^3=\frac{27}{32768},\qquad
 z_r=\left(2,\frac{56}3-\frac{1024}9r,
                    -\frac76-\frac{64}3r\right),\quad r\in\mathbb Q.
                                                               \tag{6}
\]

Every fibre over this line is full, finite étale of rank five. Its exact inverse factorization, in the unnormalized coordinate `pi=-2`, is

\[
 E_{-2,b_r,c_r}(S)=-32\left(S^2+\frac3{16}\right)
 \left[\left(S+\frac14\right)^3
       -\frac9{32}\left(S+\frac14\right)-r\right].       \tag{7}
\]

Define the already classified parameter set

\[
 \mathcal H=\left\{a\in\mathbb Q^\times\setminus\mathbb Q^{\times3}:
 \begin{array}{l}
 3\mid v_p(a)\text{ for every }p\equiv2\pmod3,\\
 a\in\mathbb Q_3^{\times3}
 \end{array}\right\}.                                  \tag{8}
\]

Then the **entire** line (6) has the exact Hasse-failure classification

\[
 \Phi^{-1}(z_r)\text{ is Hasse-failing}
 \quad\Longleftrightarrow\quad
 r=a+\frac{k}{a}\text{ for some }a\in\mathcal H.          \tag{9}
\]

For every rational `a!=0`, with `r=a+k/a`, its fibre algebra is

\[
 \mathbb Q(\sqrt{-3})\times\mathbb Q[Z]/(Z^3-a).          \tag{10}
\]

Every target in the parametrized image has exactly two rational parameters, `a` and `k/a`. There are no rational fixed points of this involution. In particular these fibres lie outside the original line whose first coordinate is `-1`.

**Proof.** Expansion proves (7). Write `U=S+1/4` and

\[
 h_r(U)=U^3-3\gamma U-r,\qquad
 \operatorname{disc}(h_r)=-27(r^2-4k).                    \tag{11}
\]

Since `k` is nonsquare in `Q`, `r^2=4k` has no rational solution, so every cubic is separable. Reduction of the cubic in (7) modulo its quadratic factor gives

\[
 -\frac9{32}S-r-\frac{25}{128}.
\]

The resultant is

\[
 \left(r+\frac{25}{128}\right)^2+\frac{243}{16384}>0.     \tag{12}
\]

Thus the two factors are always coprime. The quadratic field is `Q(sqrt(-3))`.

Suppose `r=a+k/a`. In the cubic étale algebra `B_a=Q[Z]/(Z^3-a)`, put

\[
 U=Z+\frac\gamma Z=Z+\frac\gamma a Z^2.
\]

Then `U^3-3gamma U=a+k/a`. The determinant of the columns `1,U,U^2` in the basis `1,Z,Z^2` is

\[
 1-\frac{k}{a^2}\ne0.                                  \tag{13}
\]

Thus these elements form a basis and define an isomorphism of the cubic algebras, including when `a` is a cube. This proves (10) scheme-theoretically, not only at the level of fields. The existing [rational-line classification, Theorem `thm:rational-classification`](../../unpublished/fixed-map-hasse-failures/sections/03-local-global.tex) now proves sufficiency in (9), because that theorem concerns exactly the algebra in (10).

Conversely a Hasse failure on (6) has irreducible cubic factor: otherwise that cubic supplies a rational point. Theorem FMA2 forces its quadratic resolvent to be `Q(sqrt(-3))`. By (11) this is equivalent to `r^2-4k` being a rational square. The two rational roots of `a^2-ra+k=0` are then nonzero and give `r=a+k/a`. Formula (10) and the existing rational classification force `a in H`. They are distinct because `k` is nonsquare.

Finally `k=gamma^3` is a rational cube, so `a -> k/a` preserves the inverse cubeclass and the conditions (8). The equation `a+k/a=b+k/b` factors as `(a-b)(ab-k)=0` after clearing nonzero denominators. These are exactly the two parameters claimed. ∎

For example `a=19` gives a Hasse-failing target of (6). The admissible parameters do not cover all rational `r`; rational solvability of `a^2-ra+k=0` is an essential additional condition. The complete criterion (9) includes this condition.

## 5. Corollary FMA5: arbitrary entire fibre schemes in characteristic zero

Let `K` have characteristic zero. A finite `K`-scheme occurs as an entire fibre of a Jacobian-one polynomial map `A^3_K -> A^3_K` if and only if it is finite étale, allowing the empty scheme. For **full** nonempty fibres the existing classification instead requires rank different from two.

**Proof.** Necessity follows from étaleness of Keller maps. The [common-arithmetic-fibres theorem and rank classification](../../unpublished/common-arithmetic-fibers/sections/01-keller-fibers.tex), `cor:finite-etale-realization` and `cor:rank-classification`, already realize rank one and all ranks at least three, with full fibres. It remains to realize rank two and zero without fullness.

For a rank-two algebra choose a separable presentation `K[T]/P`, with `deg P=2`, and choose `b in K` with `P(b)!=0`. Set

\[
 R(T)=P(T)(T-b)^2.
\]

Choose `a in K` such that `R'(a)R'''(a)!=0`; both derivatives are nonzero polynomials and the field is infinite. The arbitrary-seed gauge applies to `G(S)=R(a+S)-R(a)`, whose linear, cubic and top coefficients are nonzero. The same distinguished target has inverse equation `R(a+S)`. The derivative-unit fibre formula applies without squarefreeness, as stated in [the gauge section, Remark `rem:supplied-hypotheses`](../../unpublished/common-arithmetic-fibers/sections/02-quadratic-gauge.tex). By the Chinese remainder theorem,

\[
 K[T]/R\simeq K[T]/P\times K[T]/(T-b)^2.
\]

On the first factor `R'=P'(T)(T-b)^2` is a unit; on the second it is nilpotent. Localization at `R'` therefore leaves exactly `K[T]/P`. The map has geometric degree four, so this is a short rank-two fibre. The determinant-one output normalization fixes the distinguished target, exactly as in the existing construction.

For the empty scheme use `R(T)=(T-b)^3` and choose `a!=b`. Its linear and cubic coefficients after translation are nonzero. Its derivative is nilpotent modulo `R`, so the same localization is the zero ring and the entire fibre is empty; the map has geometric degree three. These arguments establish all finite étale possibilities. ∎

## 6. Theorem FMA6: sharp order of the second-line target count

Let `H` denote projective multiplicative height: a primitive integer representative `[x_0:...:x_n]` has height `max |x_i|`. Let

\[
 N_2(B)=\#\{r\in\mathbb Q:H([1:z_r])\le B,
                  \ \Phi^{-1}(z_r)\text{ is Hasse-failing}\}.
\]

There are positive constants `c_1,c_2,B_0`, independent of `B`, such that

\[
 c_1\frac{B}{\log B}\le N_2(B)\le c_2\frac{B}{\log B}
       \qquad(B\ge B_0).                               \tag{14}
\]

Thus `N_2(B)=Theta(B/log B)`. The [second-line count](FIXED_MAP_SECOND_LINE_COUNT.md) refines this to `N_2(B)~C_2 B/log B`, with the exact target-height constant, a global-solubility asymptotic, and the relative frequency of Hasse failures among locally soluble targets.

**Proof.** Write `a=m/n` in lowest terms, `m!=0,n>0`, and put `M=max(|m|,n)`. Set `D=32768` and `K=27`. Then

\[
 r=\frac{Dm^2+Kn^2}{Dmn},\qquad
 g=\gcd(Dm^2+Kn^2,Dmn).                                \tag{15}
\]

The common divisor `g` divides `DK`. Indeed a prime other than two or three cannot divide both entries, since `m,n` are coprime. At two, divisibility is possible only when `n` is even and `m` is odd. The valuations of the two summands in the numerator are then `15` and `2v_2(n)`, which cannot tie; their minimum is at most `15`. At three, divisibility is possible only when `m` is divisible by three and `n` is not. The valuations are `2v_3(m)` and `3`, which likewise cannot tie and have minimum at most `3`. These are exactly the exponent bounds in `DK`.

The numerator in (15) is positive and at least `min(D,K)M^2`; both numerator and denominator have absolute value at most `(D+K)M^2`. Dividing by `1<=g<=DK` proves the uniform comparison

\[
 H(r)\asymp H(a)^2.                                    \tag{16}
\]

The target embedding `r -> [1:z_r]` is linear over `Q`, with a rational linear inverse on its image, so clearing its fixed denominators in both directions gives `H([1:z_r])` comparable to `H(r)`. The same elementary argument gives `H([1:y_a])` comparable to `H(a)` for the original line. The notation `asymp` denotes two-sided comparison by fixed positive constants, not asymptotic ratio one.

Let `N_H(T)=#{a in H:H(a)<=T}`. The established positive-constant asymptotic `C_Phi B^2/log B` on the original line, [Theorem (All rational parameters on the line)](../../unpublished/fixed-map-hasse-failures/sections/01-introduction.tex), together with its linear height comparison, gives

\[
 N_{\mathcal H}(T)\asymp T^2/\log T.                    \tag{17}
\]

By Proposition FMA4, a target in (9) has exactly two rational parameters, both in `H`. Combining the two-sided bounds (16)–(17) bounds `N_2(B)` above and below by fixed positive multiples of

\[
 (\sqrt B)^2/\log(\sqrt B)=2B/\log B.
\]

This proves (14). In particular no new Selberg–Delange limit theorem or uncomputed leading constant has been assumed. ∎
## 7. Counting and assurance boundaries

The original-line target count, the [second-line count](FIXED_MAP_SECOND_LINE_COUNT.md), and the [full distinct-algebra count](FIXED_MAP_DISTINCT_ALGEBRA_COUNT.md) now all have leading asymptotics in their stated target heights. The last constant is given by a convergent finite-prime area limit rather than an asserted multiplicative Euler product. The [ambient bounds](FIXED_MAP_AMBIENT_COUNT_BOUNDS.md) prove `A(B)<<B^3` for Hasse failures and sharp order `Theta(B^3)` for both globally and everywhere locally soluble rational targets. The exact ambient Hasse order and the full ambient leading constants remain open, together with compatible integral-source conditions.

The [local-density note](FIXED_MAP_LOCAL_DENSITIES.md) classifies the original-line algebras and compares all, locally soluble, globally soluble, and Hasse-failing targets in the same height. The [integer-count note](FIXED_MAP_INTEGER_COUNT.md) counts every signed integer Hasse parameter, including a second asymptotic coefficient. These companion counts retain their explicit line and parameter domains.

For the second line, `r=a+k/a` is degree two, with two rational preimages throughout its rational image. Its target height cannot be replaced by the original linear target height in `a`; the second-line count explicitly treats its quadratic height region and finite content corrections. The second-line count has a different domain from the original-line count. The [fixed-quadratic rigidity theorem](FIXED_QUADRATIC_LINE_RIGIDITY.md) proves that exactly these two affine target lines can contain Hasse failures while retaining one fixed irreducible quadratic factor. Its third fixed-quadratic line has local points exactly at its rational-point parameters. Moving quadratic factors remain a separate architecture.

The script [verify_fixed_map_ambient.py](verify_fixed_map_ambient.py) checks exact polynomial identities for (7), (11)–(13), the algebra generator relation and involution; boundary reconstruction identities; and exhaustive subgroup controls for the `S_3` local criterion and the wrong-resolvent common derangements. It also checks the rank bound over all partitions of five and one rank-two localization example. These controls do not verify Chebotarev, the general Hensel lemma, all arithmetic local conditions, or the underlying all-degree Keller-map construction. No numerical sampling is used to justify the unrestricted statements.

Run from the repository root with Python 3 and SymPy:

```sh
.venv/bin/python research/arithmetic/verify_fixed_map_ambient.py --output research/arithmetic/fixed_map_ambient_checks.json
```

The proof sources and exact external application are cited above. The written reductions, the finite symbolic checks, existing formalizations and external review remain distinct assurance levels.
