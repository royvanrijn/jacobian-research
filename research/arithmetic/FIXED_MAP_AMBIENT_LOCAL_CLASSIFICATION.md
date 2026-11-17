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

is monic. It is enough to test primes dividing the nonzero integer `Disc(H)`. At each such prime `p`, put `v=v_p(Disc(H))`. The assertion that `H` has a `Q_p`-root is equivalent to the existence of a root modulo `p^(floor(v/2)+1)`, by Lemma FMA3a below. Thus the criterion contains only rational polynomial factorization, rational-square testing, and finitely many finite congruence tests. The cubic discriminant exponent is optimal in general.

More generally the fibre is everywhere locally soluble if and only if either it has a rational point or it satisfies the Hasse-failure criterion above. For `u=0` it always has a rational point. For `u!=0` it has a rational point exactly when `E_{-u,b,c}` has a **simple rational root**; a repeated rational root removed by (2) does not suffice.

**Proof.** Lemma FMA1 reduces any Hasse failure to a full rank-five algebra. Theorem FMA2 gives precisely the factorization, common quadratic-resolvent condition (5), and local tests. Conversely these conditions give the algebra (3) and local points everywhere, with no rational point.

For the finite prime set, `H` defines the same cubic field as `h`. At primes outside `Disc(H)` its three roots give an unramified étale algebra over `Z_p`; its splitting field is consequently unramified. Such decomposition groups are cyclic, and Theorem FMA2 makes the local condition automatic. Condition (5) also places the quadratic field inside that same splitting field. Extra primes dividing the index of `Z[T]/H` cause harmless extra tests.

Lemma FMA3a supplies the sharper finite cubic test. A root of a monic integral polynomial is `p`-adically integral, so its reductions give the converse direction. The rational-point assertion is immediate from (2). ∎

This is an arithmetic decision criterion for every rational target, not a parametrization of all targets and not a counting theorem. The condition that the fibre have a specified factorization and resolvent is essential; the generic quintic covering has symmetric monodromy.

## 3a. Lemma FMA3a: optimal cubic discriminant root test

Let `H in Z_p[T]` be monic, separable and cubic, and put
`d=v_p(Disc(H))`. Then

\[
 H\text{ has a root in }\mathbb Q_p
 \quad\Longleftrightarrow\quad
 H\text{ has a root modulo }p^{\lfloor d/2\rfloor+1}.
                                                        \tag{5a}
\]

The statement applies to every prime, including two and three. For
`p!=3` its exponent cannot be decreased in general: the cubics
`T^3-p^k`, with `k>0` not divisible by three, have `d=2k`, a root modulo
`p^k`, and no `Q_p`-root.

**Proof.** If there is a local root, monicity makes it integral and it
reduces modulo every power of `p`. Suppose instead that there is no local
root. A cubic with no root over a field is irreducible. In a finite Galois
splitting field, extend `v_p` with its normalization `v_p(p)=1`; its value
group need not be `Z`. Uniqueness of the extension makes it invariant
under the Galois action. This exact valuation input is
[Milne, *Algebraic Number Theory*, version 3.00, Theorem 7.38 and its proof,
printed pp.116–117](https://www.jmilne.org/math/CourseNotes/ANT300.pdf),
retained in [the source copy](sources/milne-ant300.pdf).

For any `x in Z_p`, transitivity on the three roots `alpha_i` therefore
makes all `v_p(x-alpha_i)` equal to a value `r`. Since `H` is monic,

\[
 v_p(H(x))=3r,
 \qquad
 v_p(\alpha_i-\alpha_j)\ge r.
\]

The second inequality is the nonarchimedean triangle inequality. Hence

\[
 d=2\sum_{i<j}v_p(\alpha_i-\alpha_j)
       \ge6r=2v_p(H(x)).                              \tag{5b}
\]

An integer representative of a root modulo
`p^(floor(d/2)+1)` contradicts (5b). Thus such a residue root forces a
local root. For the sharpness family the discriminant is `-27p^(2k)`;
a local root would have nonintegral valuation `k/3`, whereas the residue
class zero is a root modulo `p^k`. This proves every assertion. ∎

The older resultant/Hensel bound `p^(2d+1)` remains sufficient but is
unnecessary here. The improvement uses the special fact that a cubic
without a local root is irreducible. The same exponent is not asserted
for arbitrary-degree polynomials.

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

## 6a. Theorem FMA7: unified fibre existence and Hasse rank spectrum

Let `K` be a characteristic-zero field, and let `X` be a finite `K`-scheme.

1. `X` occurs as an entire fibre over a `K`-rational target of a Jacobian-one polynomial map `A^3_K -> A^3_K` if and only if `X` is finite étale. This includes the empty scheme.
2. If `X` has positive rank `N`, it occurs as a **full** fibre of such a map of geometric degree `N` if and only if it is finite étale and `N!=2`.
3. Over `Q`, a full Hasse-failing fibre of rank `N` exists for such a map if and only if `N>=5`. The map can always be chosen on `A^3_Q` with geometric degree exactly `N`.

The first two assertions also hold on every `A^m_K` with `m>=3`, by adjoining identity coordinates. They classify possible fibre schemes up to `K`-isomorphism while allowing the map to vary. The third assertion is an existence spectrum, not an assertion that every rank-`N` algebra is Hasse-failing.

For a given finite étale algebra over a number field there is a simultaneous arithmetic classification. Let `M/K` be a finite Galois splitting field, `G=Gal(M/K)`, and

\[
 \Omega=\operatorname{Hom}_{K\text{-alg}}(A,M).
\]

For each place `v` choose a decomposition group `D_v` in `G`, well defined up to conjugacy. Then

\[
 \operatorname{Spec}A\text{ is Hasse-failing}
 \quad\Longleftrightarrow\quad
 \Omega^G=\varnothing\quad\text{and}\quad
 \Omega^{D_v}\ne\varnothing\text{ for every place }v.    \tag{18}
\]

Equivalently, every element of `G` must fix some point of `Omega`, the finitely many ramified finite decomposition groups must each have a common fixed point, and the archimedean decomposition groups must each have a common fixed point, while `G` itself has none. This gives a finite group and local test once the splitting field and decomposition data are supplied. The element-cover condition alone does not replace the ramified-group tests.

**Proof.** Assertion 1 is Corollary FMA5. For assertion 2, étaleness is necessary. The identity map gives the unique rank-one étale scheme. The characteristic-zero Common Arithmetic Fibers theorem realizes every étale algebra of rank at least three with geometric degree equal to its rank. A rank-two full fibre would require geometric degree two; its characteristic-zero function-field extension is Galois, which is excluded by the classical Campbell–Razar–Wright Galois case of the Jacobian theorem. This is exactly [the existing rank-classification corollary](../../unpublished/common-arithmetic-fibers/sections/01-keller-fibers.tex), `cor:rank-classification`; the classical exclusion is an imported written input separate from the end-to-end formalized construction. Empty schemes are entire fibres but cannot be full, since a dominant generically finite map has positive geometric degree. Stabilization by identity coordinates retains the fibre and generic degree.

There is also a narrower written route for the degree-two exclusion, without the general Galois theorem. Over an algebraically closed characteristic-zero field `k`, [the divisor power-saturation identity (3.1)](../jacobian/monodromy-theory/THEORY.md) gives `b(F) in k(X)^{*2} => b in k(Y)^{*2}` for every whole-affine Keller map. A quadratic generic extension would have the form `k(Y)(sqrt(b))` for a nonsquare `b`, contradicting that identity. Generic degree is preserved by extension of the constant field, so this also excludes degree two over the original field `K`. The argument uses only the squarefree pullback of target prime divisors; it does not require the broader monodromy classification.

For assertion 3, Lemma FMA0 gives the lower bound five. Here are rank-five and rank-six base algebras:

\[
 A_5=\mathbb Q(\sqrt{-3})\times\mathbb Q[Z]/(Z^3-19),
 \qquad
 A_6=\mathbb Q(\sqrt{13})\times\mathbb Q(\sqrt{17})
                         \times\mathbb Q(\sqrt{221}).   \tag{19}
\]

The existing local criterion proves that `A_5` is Hasse-failing: `19` is a noncube, its sole prime divisor is `1 mod3`, and `19=1 mod9` is a three-adic cube. It has no rational factor.

None of the three factors of `A_6` is `Q`, so it has no rational point. At `13`, the field `Q(sqrt(17))` supplies a point because `17=2^2 mod13` with nonzero derivative. At `17`, `Q(sqrt(13))` supplies a point because `13=8^2 mod17`. At `2`, `Q(sqrt(17))` supplies a point: for `f(T)=T^2-17`, `v_2(f(1))=4>2v_2(f'(1))=2`, so Hensel's lemma applies. At every other odd prime all three radicands are units. Their quadratic residue characters have product one, since `13*17*221=221^2`; they cannot all be minus one, and a quadratic residue gives a local point by Hensel's lemma. At the real place all three radicands are positive. Thus `A_6` is Hasse-failing.

For odd `N>=5`, take `A_5` and `(N-5)/2` additional factors `Q(sqrt(2))`; for even `N>=6`, take `A_6` and `(N-6)/2` such factors. Products retain a local point at every place through their base algebra and gain no rational point, because every new field factor has degree two. Repeated isomorphic field factors are distinct étale components and cause no problem: the infinite-field monogenicity theorem chooses a separable presentation for the resulting product. The Common Arithmetic Fibers theorem now realizes the product as a full fibre on `A^3_Q` of geometric degree `N`. Conversely no rank below five can be Hasse-failing. More generally the fibre-rank bound in [the existing minimality proof](../../unpublished/fixed-map-hasse-failures/sections/05-minimality.tex) shows that any Keller map with a Hasse-failing fibre has generic degree at least five.

Finally the finite étale/Galois-set correspondence identifies `K`-points with `Omega^G`, and `K_v`-points with `Omega^(D_v)`, proving (18). At an unramified finite place the decomposition group is cyclic, generated by Frobenius, so it fixes a point exactly when that element does. Chebotarev realizes every conjugacy class at unramified places; therefore the element-cover condition is both necessary and sufficient for all unramified places. The remaining ramified and archimedean tests are exactly those in (18). This is the fixed-point formulation of Sonn's Proposition 2.1 cited in section 2, applied to the supplied algebra.

To see why testing individual elements is weaker at ramified places, let `G=V_4` and let `Omega` be the disjoint union of its two cosets for each of its three order-two subgroups. Every element fixes a point, because it lies in one of those subgroups. The whole group fixes no point. A decomposition group equal to `V_4` would therefore obstruct local solubility despite the element-cover condition. The explicit algebra `A_6` above avoids that obstruction, as the direct local proof established. ∎

This theorem supplies the complete characteristic-zero classification of possible fibre schemes with a varying map, together with the exact Hasse rank spectrum over `Q` and the given-algebra local-global criterion. It does not classify Keller maps up to equivalence, prescribe their generic monodromy, solve arbitrary inverse-Galois/decomposition-group realization problems, or give an ambient counting asymptotic for a fixed map. The characteristic hypotheses and the distinction between entire and full fibres are essential. No extension to all characteristics or new formal verification or external review is claimed.

## 7. Counting and assurance boundaries

The [fixed-plane local counting theorem](FIXED_MAP_PLANE_LOCAL_COUNTS.md)
shows that `A_P(B)=O_P(B^2)` on every fixed rational affine plane.
Outside the two section planes `u=0,c=0`, the same bound holds for
locally soluble targets. These two planes are therefore exactly the
planes with order `B^3` local/global counts, and both have no Hasse
failures. The argument uses degree-five-or-ten components of the
quadratic-factor cover; the constant may depend on the plane.

The [five-prime elliptic count](ELLIPTIC_EXCEPTION_JOINT_2_3_5_11_17_SIEVE.md)
improves the Hasse upper coefficient to `116302241/531619620`, with
`O_E(sqrt(log B))` error. This is the exact coefficient for tests at
`2,3,5,11,17`. It replaces the four-prime coefficient as the current
best upper bound while preserving that earlier exact prefix count.
The certified quadratic field still has upper coefficient `1/8`.

The elliptic sieve has an additional arithmetic support restriction.
[Kummer branch exclusion](ELLIPTIC_EXCEPTION_KUMMER_BRANCH_EXCLUSION.md)
gives a positive-density set of primes `p=2 mod3` which pass every
quadratic eligibility test but never divide the reduced numerator
of `128x-25865` at any rational elliptic point. Only the pole can
contribute to the local hazard there. The proof uses the independent
branch point over `Q(sqrt(51730))`, a full Kummer translation kernel,
and rational-prime Chebotarev; it asserts no recurrence of pole
returns or all-places count.

[Finite-primary separation](ELLIPTIC_EXCEPTION_FINITE_PRIMARY_SEPARATION.md)
removes any fixed finite collection of primary correlations on a
refined affine cell. Its abstract divisibility countermodel shows
why this and divergent fixed-cell return means still do not imply
recurrence: large common-quotient clustering remains uncontrolled.
The countermodel is not an elliptic curve or an arithmetic
counterexample.

[Finite obstruction certificates](ELLIPTIC_EXCEPTION_FINITE_OBSTRUCTION_CERTIFICATES.md)
describe simultaneous infeasibility as a zero-weight cover of a
certified finite rational congruence image. At a sufficiently large
largest prime, every previously surviving reduction would have to
lie in the two fatal finite zero branches, with a bad offset and
formal image depth at least two. The Kummer primes cannot terminate
such an obstruction. This isolates the remaining feasibility question
for the two classes without an all-finite-place positivity theorem.
[Exact common quotients and moments](ELLIPTIC_EXCEPTION_COMMON_QUOTIENT_MOMENTS.md)
give the actual return-correlation formulas. Divergence of the
conditional mean with variance negligible relative to its square
would force divergent hazard and a sublogarithmic Hasse count;
the required increasing-prime estimates remain open.

The elementary [small-image height bound](ELLIPTIC_EXCEPTION_SMALL_REDUCTION_IMAGES.md)
shows that only `O_E(z^2/log z)` good primes have rational reduction
image smaller than `z`. For each fixed `A>0`, this gives absolute
new-prime coefficient loss at most `3(log p)^A/sqrt(p)` on a set of
prime density one, uniformly over all finite prefixes and classes.
This improves the marginal estimate without providing the conditional
correlation or all-primes tail control.

The [explicit four-prime core](ELLIPTIC_EXCEPTION_DEPENDENCY_CLOSED_CORE.md)
gives a fixed safe-prime pool of exact density `13/16+delta_Kum`:
every finite subset leaves logarithmically many full field targets
in each of the four quadratic classes. The two unresolved classes
have separate pools of density `7/8+delta_Kum`. Fatal branches
remain the only possible first failure outside the core.
[Rank-three branch bounds](ELLIPTIC_EXCEPTION_UNIFORM_OBSTRUCTION_PRIMES.md)
give a finite computable candidate set, uniform over prefixes of a
given positive mass, and square-summable numerator-return masses
with a uniform tail for two independent Haar vectors. The
[quadratic variance obstruction](ELLIPTIC_EXCEPTION_POLE_VARIANCE_OBSTRUCTION.md)
rules out negligible relative variance for unconditional pole sums.
These results do not settle the conditioned one-vector tail or
Hasse infinitude.

[Three-torsion disjointness](ELLIPTIC_EXCEPTION_THREE_TORSION_DISJOINTNESS.md)
gives exact compatible Frobenius densities and enables
[pole opportunities in every affine coset](ELLIPTIC_EXCEPTION_AFFINE_COSET_POLE_MEANS.md).
Every fixed feasible prefix has unweighted expected pole count
`>>_(S,a) log log Y`; weighted hazard and covariance control remain
separate. An [exact bounded certificate](ELLIPTIC_EXCEPTION_BOUNDED_FATAL_PRIME_CERTIFICATE.md)
makes every finite set of tests below `100000` feasible in all
four classes. The [class-group bridge](ELLIPTIC_EXCEPTION_WIEFERICH_CLASS_GROUP_BRIDGE.md)
identifies actual equivariant ideal-class-group quotients forced
by the rational formal depth, providing a global arithmetic
constraint on possible fatal primes.

The [rank-five discriminant classification](RANK_FIVE_DISCRIMINANT_CLASSIFICATION.md) replaces the finite local tests by a condition on the maximal-order cubic field discriminant. Write `Disc(L)=d_K f^2`, where `d_K` is the fundamental discriminant of the quadratic resolvent. A matching field product `K x L` is Hasse-failing exactly when every prime divisor of `f` other than three is `1 mod3`, and either `3` does not divide `f` or `d_K=1 mod3`. This applies to every fixed-map full `2+3` field fibre and, with FMA2, classifies all rank-five finite étale Hasse failures over `Q`. A generator polynomial discriminant cannot replace the field discriminant: even `T^3-19` introduces an extraneous factor at three. For integral fixed-map Hasse targets, the earlier ramification theorem forces the stronger condition `3` does not divide `f`.

The [discriminant count](RANK_FIVE_DISCRIMINANT_COUNT.md) applies the acceptable-local-specification Davenport–Heilbronn theorem to these exact conditions. Rank-five Hasse-failing algebra classes, counted by the absolute discriminant of their cubic field factor, have an explicit positive-density leading asymptotic. Common Arithmetic Fibers realizes every such class as a full fibre of a degree-five map while allowing that map to vary. This discriminant count does not give an asymptotic for the rational targets of the one fixed map.

The original-line target count, the [second-line count](FIXED_MAP_SECOND_LINE_COUNT.md), and the [full distinct-algebra count](FIXED_MAP_DISTINCT_ALGEBRA_COUNT.md) now all have leading asymptotics in their stated target heights. The last constant is given by a convergent finite-prime area limit rather than an asserted multiplicative Euler product. The [ambient bounds](FIXED_MAP_AMBIENT_COUNT_BOUNDS.md) prove `A(B)<<B^3` for Hasse failures and sharp order `Theta(B^3)` for both globally and everywhere locally soluble rational targets. The exact ambient Hasse order and the full ambient leading constants remain open, together with compatible integral-source conditions.

The [general resolvent fibration](FIXED_MAP_RESOLVENT_FIBRATION.md) covers every noncentered Hasse target, retaining its singular specializations; centered targets give exactly the three known lines. Its generic fibres are genus two, with rational two-torsion on their Jacobians. For each fixed smooth base pair, Faltings' theorem makes the Hasse-target count finite, even though all factorized targets on the associated rational cubic curve have order `B^(2/3)`. This is an individual finiteness theorem, with no uniform or effective ambient bound. The [integral local analysis](FIXED_MAP_INTEGRAL_LOCAL_OBSTRUCTIONS.md) supplies exact three-adic conditions, two-adic exclusions and automatic local points, and proves that every integral Hasse target has a quadratic resolvent ramified at three. These restrictions leave integral existence unresolved.

On the explicit elliptic exceptional slice, the [arithmetic classification](ELLIPTIC_EXCEPTION_ARITHMETIC.md) gives exactly four quadratic fields, rank two, torsion of order two, and a genus-four cubic-root cover. Consequently only finitely many matching targets have reducible cubic factor, and each quadratic field occurs with infinitely many distinct cubic-field and fibre-algebra classes. The [elliptic target count](ELLIPTIC_EXCEPTION_COUNT.md) proves a leading `C_E log B` asymptotic in actual target height, with each quadratic field contributing one quarter. This counts matching field products; the Hasse subset still requires its local tests. An [explicit certified Hasse target](ELLIPTIC_EXCEPTION_HASSE_EXAMPLE.md) has quadratic resolvent `Q(sqrt(-81275190))`, proving that the pure-cubic/cyclotomic-quadratic families do not exhaust the fixed map's Hasse failures. The [two-adic refinement](FIXED_MAP_TWO_ADIC_REFINEMENT.md) further gives exact conditional local criteria for even first coordinate and for odd first coordinate with even `c`, and narrows the remaining odd/odd integral branch to two patterns.

The [elliptic local sieve](ELLIPTIC_EXCEPTION_LOCAL_SIEVE.md) shows that any fixed finite collection of local conditions retains a full field-fibre subfamily with `c_S log B+O_S(sqrt(log B))` targets and `c_S>0`. The [full finite-place asymptotic](ELLIPTIC_EXCEPTION_FINITE_LOCAL_ASYMPTOTICS.md), EFA1--EFA2, counts the entire prescribed-local set in each quadratic class: `N_(S,a)(B)=mu_(S,a) C_E log B+o(log B)`. Its constant is Haar measure in the actual joint closure of the rational elliptic group, and is positive for every finite `S` in the certified field `K0=Q(sqrt(-81275190))`. EFA3 gives an exact tail criterion for passage to all primes. The limiting measure can be zero despite positivity at every finite stage; if it is positive, a vanishing normalized large-prime tail would give a positive Hasse asymptotic. Neither assertion has been proved. The exact all-places condition on this slice is that the maximal-order factor `f` be squarefree and supported on primes `1 mod3`, since all four quadratic fields ramify at three. Controlling the varying additional bad primes on this existing slice remains the arithmetic task.

The [explicit elliptic criterion](ELLIPTIC_EXCEPTION_EXPLICIT_LOCAL_CRITERION.md) makes those varying tests concrete: outside `2,3,5`, they are exactly divisibility by three of `v_p(128x-25865)` for primes `p=2 mod3`, where `x=X+3552`. Any one cubeclass of this rational function gives a genus-four curve and only finitely many rational elliptic points. The [complete three-adic test](ELLIPTIC_EXCEPTION_COMPLETE_THREE_ADIC_CLASSIFICATION.md) and [two/five-adic test](ELLIPTIC_EXCEPTION_TWO_FIVE_ADIC_CLASSIFICATION.md) make the remaining small-prime conditions explicit and determine the rational local images. Exact finite-local logarithmic proportions are `55/156` at three, `4/7` at two, and `289/1092` jointly at two and three. The [joint sieve at two, three and seventeen](ELLIPTIC_EXCEPTION_JOINT_2_3_17_SIEVE.md) gives Hasse upper coefficient `1687/6936<1/4` with `O_E(sqrt(log B))` error, and exact local-only coefficient `5321/23946` with `o(log B)` error. J2174 separates this last coefficient among the four quadratic fields; the certified field has coefficient `1/8` relative to the whole-slice `C_E log B`. These finite-place counts give Hasse upper bounds; all-places infinitude remains open.

The [quantitative finite-place theorem](ELLIPTIC_EXCEPTION_QUANTITATIVE_FINITE_LOCAL_ASYMPTOTICS.md) proves that every finite-place elliptic coefficient is rational. The [sharp error theorem](ELLIPTIC_EXCEPTION_SHARP_FINITE_LOCAL_ERROR.md) uses Hirata-Kohno--Takada's bound for two p-adic elliptic logarithms to give `O_(E,S)(sqrt(log B))` error for every fixed finite `S`. The [complete joint prefix at two, three, five and seventeen](ELLIPTIC_EXCEPTION_JOINT_2_3_5_17_SIEVE.md) has exact coefficient `2300089/10392564`. The [five-prime refinement](ELLIPTIC_EXCEPTION_JOINT_2_3_5_11_17_SIEVE.md) improves the Hasse upper coefficient to `116302241/531619620` with that error; the certified field retains coefficient `1/8`. Five changes only the class-one coefficient, multiplying it by `16/31`. The earlier coefficient `1687/6936` retains its direct finite-quotient proof. A nonempty finite-local field-target set has a positive rational coefficient. [Class-one positivity](ELLIPTIC_EXCEPTION_FINITE_LOCAL_CLASS_ONE.md) establishes this for every finite `S` in `Q(sqrt(-3))` without an all-places seed; two quadratic classes have that property. The [one-prime theorem](ELLIPTIC_EXCEPTION_ONE_PRIME_POSITIVITY.md) proves positivity at every individual place in all four classes, and cofinite one-prime extensions of any feasible finite set, without asserting simultaneous positivity for every finite set in the remaining classes. The [finite-quotient saturation theorem](ELLIPTIC_EXCEPTION_FINITE_QUOTIENT_SATURATION.md) computes the actual rational image and coefficients once a certified finite ambient quotient is supplied, without a full Mordell--Weil basis.

The [class-specific eligibility sieve](ELLIPTIC_EXCEPTION_PRIME_ELIGIBILITY.md) erases every inert-prime test where the global class `a` is not locally square; numerator obstructions also require `51730` locally square. For each nontrivial class these eligible prime sets have densities `1/4` and `1/8`. The [ambient Haar diagnostic](ELLIPTIC_EXCEPTION_AMBIENT_HAAR_DIAGNOSTIC.md) proves that both the full-product all-soluble set and the rational closure have ambient measure zero. Its finite-prime support index grows at least as `2^(s-1)` for `s` good odd primes. The [triangular sieve](ELLIPTIC_EXCEPTION_TRIANGULAR_PRIME_SIEVE.md) proves independence of the fresh formal `p`-coordinate from all previous tests. Its all-depth conditional hazard kernel retains the correlated finite reduction and branch offset. The all-primes soluble set has finite accumulated predictable hazard almost surely, and a finite-cutoff hazard inequality gives explicit conditional count bounds. Proving weighted divergence for typical vectors remains open; no separate depth-one distribution assumption is required by that reduction. Katz's simultaneous elliptic Wieferich conjecture supplies context for a conditional depth-one-density statement, not an input to these unconditional results. The [branch-twist theorem](ELLIPTIC_EXCEPTION_BRANCH_TWIST.md) supplies a nontorsion branch point independent of the rational rank-two group over `Q(sqrt(51730))`, giving a rank-at-least-three input for future division-field work. The all-primes infinitude problem remains unresolved.

The one-prime theorem also proves `|im(E(Q)->E(F_p))|>>_E log p` at large good primes. Consequently the absolute loss of a finite-local coefficient on adding `p` is `O_E(1/log p)`, uniformly over the existing finite conditions and class. A fixed positive prefix therefore has retention tending to one as the new prime grows. The estimate is nonsummable over primes and does not control the relative loss of a shrinking prefix or the growing-prime counting error.

The [local-density note](FIXED_MAP_LOCAL_DENSITIES.md) classifies the original-line algebras and compares all, locally soluble, globally soluble, and Hasse-failing targets in the same height. The [integer-count note](FIXED_MAP_INTEGER_COUNT.md) counts every signed integer Hasse parameter, including a second asymptotic coefficient. These companion counts retain their explicit line and parameter domains.

For the second line, `r=a+k/a` is degree two, with two rational preimages throughout its rational image. Its target height cannot be replaced by the original linear target height in `a`; the second-line count explicitly treats its quadratic height region and finite content corrections. The second-line count has a different domain from the original-line count. The [fixed-quadratic rigidity theorem](FIXED_QUADRATIC_LINE_RIGIDITY.md) proves that exactly these two affine target lines can contain Hasse failures while retaining one fixed irreducible quadratic factor. Its third fixed-quadratic line has local points exactly at its rational-point parameters.

The [moving-factor quadric](FIXED_MAP_MOVING_QUADRATIC_FAMILY.md) supplies a rational surface of translated-pure-cubic factorizations, its exact local criterion, and an explicit Hasse target outside both lines. No Hasse target in that particular branch is three-adically integral; this does not settle integral existence for general cubic generators or quadratic fields. The [plane-section theorem](FIXED_MAP_PLANE_SECTIONS.md) proves that `u=0` and `c=0` are the only rational affine planes with generic rational sections. Their union doubles the known ambient lower coefficient to `8/zeta(3)`; every other fixed plane contributes `O_P(B^2)` globally soluble targets. Neither result determines the full ambient leading constants.

The script [verify_fixed_map_ambient.py](verify_fixed_map_ambient.py) checks exact polynomial identities for (7), (11)–(13), the algebra generator relation and involution; boundary reconstruction identities; and exhaustive subgroup controls for the `S_3` local criterion and the wrong-resolvent common derangements. It also checks the rank bound over all partitions of five, one rank-two localization example, the rank-six `V_4` fixed-point action, and the displayed residue controls. These controls do not verify Chebotarev, the general Hensel lemma, all arithmetic local conditions, or the underlying all-degree Keller-map construction. No numerical sampling is used to justify the unrestricted statements.

Run from the repository root with Python 3 and SymPy:

```sh
.venv/bin/python research/arithmetic/verify_fixed_map_ambient.py --output research/arithmetic/fixed_map_ambient_checks.json
```

The proof sources and exact external application are cited above. The written reductions, the finite symbolic checks, existing formalizations and external review remain distinct assurance levels.

The [strict Selmer converse](ELLIPTIC_EXCEPTION_CLASS_GROUP_SELMER_CONVERSE.md)
identifies the level-one class-group maps exactly. At good primes
`p>=7` with full `GL_2(F_p)` image, excluding the finite bad-local
orders specified there,
`Hom_G(Cl(Q(E[p]))/p,E[p])` is the Selmer kernel of localization
at `p`. In particular the class group has an actual equivariant
quotient `E[p]^2` exactly when `e_p>1` or `Sha(E/Q)[p]` is
nonzero. No finiteness of `Sha` is assumed. This explains the
possible extra source in the density problem without proving a
prime-density estimate or a return-correlation bound.
