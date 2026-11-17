# Whole-plane stable multiplicity

The common power shifts and fiber-invisible cubic lifts were introduced to
keep one selected finite etale fiber fixed while changing the stable
polynomial left--right class of the ambient map.  Their fiber invisibility is
not point-specific.  On the source hypersurface where the first map
coordinate is one, the deformed maps themselves agree.  Consequently they
share the entire squarefree inverse cover over the two-dimensional target
plane `P=1`.

This strengthens common-fiber multiplicity in every noninvertible rank and
gives one countable family of stably inequivalent quintic maps carrying the
same dense family of degree-optimal Hasse-principle failures.

## 1. The common inverse plane

Let `K` have characteristic zero and let

\[
 G(S)=\sum_{j=1}^N g_jS^j,\qquad g_1g_3g_N\ne0.
 \tag{1.1}
\]

On the target plane `P=1`, put

\[
 \mathcal E_G(B,C;S)
 =G(S)-\frac{g_1}{2}(BS^2+C)
 \tag{1.2}
\]

and let

\[
 \Delta_G(B,C)=\operatorname{disc}_S\mathcal E_G(B,C;S).
 \tag{1.3}
\]

Since the leading coefficient of (1.2) is the nonzero constant `g_N`, the
open set

\[
 U_G=D(\Delta_G)\subset\mathbb A^2_{B,C}
 \tag{1.4}
\]

is exactly its degree-`N` squarefree locus.  Define the finite etale root
cover

\[
 \mathcal X_G
 =\operatorname{Spec}
 K[B,C,\Delta_G^{-1},S]/(\mathcal E_G)
 \longrightarrow U_G.
 \tag{1.5}
\]

The cover has rank `N`.  Its fibers retain the complete quotient algebra,
not merely the geometric root count.

## 2. Whole-plane theorem in ranks at least four

Assume `N>=4` and

\[
 g_1g_3g_4\cdots g_N\ne0.
 \tag{2.1}
\]

For `m>=0`, the common power-shifted inverse seed is

\[
 G_{P,m}(S)
 =g_1S+P(g_2S^2+g_3S^3)
  +\sum_{j=4}^Ng_jP^{j+m}S^j.
 \tag{2.2}
\]

Let

\[
 F_m=(P,B_m,C_m):\mathbb A^3_K\longrightarrow\mathbb A^3_K
 \tag{2.3}
\]

be the associated determinant-minus-two map.  Write the common first source
coordinate as

\[
 t=1+xy,\qquad
 q=t^2z+\frac{g_1}{g_3}y^2(1+3t),\qquad P=tq.
 \tag{2.4}
\]

Every shifted higher term in `B_m` and `C_m` satisfies

\[
\begin{aligned}
 t^{m+2}x^{j-2}q^{j+m}
 &=P^m t^2x^{j-2}q^j,\\
 t^mx^jq^{j+m}
 &=P^m x^jq^j.
\end{aligned}
\tag{2.5}
\]

Thus, for polynomials `H_B,H_C` independent of `m`,

\[
 F_m=(P,B_{\mathrm{low}}+P^mH_B,
         C_{\mathrm{low}}+P^mH_C).
 \tag{2.6}
\]

In particular,

\[
 \boxed{F_m\equiv F_{m'}\pmod{P-1}}
 \tag{2.7}
\]

coordinatewise for all `m,m'`.  This is equality of the restricted
polynomial maps on the common source divisor `P=1`.

At a target `(1,B,C)`, the inverse equation of every `F_m` is (1.2).
If `(B,C)` lies in `U_G`, its derivative is a unit in the quotient algebra.
The scheme-theoretic reconstruction formulas therefore identify the literal
source fiber with that quotient.  They globalize after inverting
`\Delta_G`, and give

\[
\boxed{
 F_m^{-1}(\{1\}\times U_G)
 \simeq\mathcal X_G
 \quad\text{over }U_G
}
\tag{2.8}
\]

for every `m`.  The isomorphisms use the same root coordinate and the same
source reconstruction; equivalently, (2.8) is one literal common cover, not
only a collection of abstractly isomorphic fibers.

The ambient maps nevertheless remain pairwise stably inequivalent.  The
normalized ramified-stratum Fitting Newton polygon has area

\[
 2N-3+(N-2)m,
 \tag{2.9}
\]

which is strictly increasing in `m` and is preserved by stable polynomial
left--right equivalence.

### Theorem 2.1

For every characteristic-zero field `K`, every `N>=4`, and every seed
satisfying (2.1), infinitely many pairwise stably inequivalent
geometric-degree-`N` Keller maps restrict to the same finite etale cover
`\mathcal X_G->U_G` over the whole squarefree `P=1` target plane.

A fixed output scaling makes every map determinant one without changing the
conclusion; it only applies the same linear coordinate change to the common
target plane.

## 3. Whole-plane theorem in rank three

Let

\[
 G(S)=g_1S+g_2S^2+g_3S^3,\qquad g_1g_3\ne0.
 \tag{3.1}
\]

For `n>=4`, the cubic lift has inverse seed

\[
 G_{P,n}(S)
 =g_1S+g_2PS^2+
   g_3P(1+P^{n-1}-P^2)S^3.
 \tag{3.2}
\]

The paired source corrections can be factored as

\[
\begin{aligned}
 \Delta B_n
 &=3\frac{g_3}{g_1}t^2xq^3(P^{n-3}-1),\\
 \Delta C_n
 &=-\frac{g_3}{g_1}x^3q^3(P^{n-3}-1).
\end{aligned}
\tag{3.3}
\]

They vanish identically on `P=1`, and (3.2) becomes `G(S)`.  Hence all
cubic lifts restrict to the same map on `P=1` and share the rank-three
version of (1.5).  Their complete canonical boundary has `n` geometric
target components, so distinct `n` give distinct stable classes.

### Theorem 3.1

The conclusion of Theorem 2.1 holds in rank three for the fiber-invisible
cubic lifts.

Combining Theorems 2.1 and 3.1 gives whole-plane stable multiplicity in every
possible noninvertible geometric degree:

\[
\boxed{\text{every }N\ge3.}
\tag{3.4}
\]

## 4. Cyclotomic strengthening

The equality is not limited to `P=1` after passing to an infinite
subfamily.  Fix `d>=1`.

For power shifts, choose one residue class

\[
 m\equiv r\pmod d.
 \tag{4.1}
\]

For two such shifts,

\[
 P^m-P^{m'}\quad\text{is divisible by}\quad P^d-1.
 \tag{4.2}
\]

Equations (2.2) and (2.6) show that the corresponding maps agree on the
entire source divisor `P^d=1`, and their inverse equations agree over the
target divisor `P^d=1`.  They remain pairwise stably inequivalent by (2.9).

For cubic lifts, choose `n` in one residue class modulo `d`.  Both (3.2) and
(3.3) depend on `n` through a power of `P`, so the same divisibility proves
equality on `P^d=1`.  On components where the common cubic leading
coefficient is nonzero, delete its repeated-root discriminant.  The
remaining cover is finite etale of rank three and is common to the infinite
subfamily.

### Theorem 4.1

For every `d>=1` and every rank `N>=3`, there is an infinite family of
pairwise stably inequivalent Keller maps sharing the full squarefree inverse
cover over the appropriate open subset of the cyclotomic target divisor

\[
 P^d=1.
 \tag{4.3}
\]

Over a field containing the `d`-th roots of unity, this is a simultaneous
common-cover statement on a union of parallel target planes `P=\zeta`.
Over the ground field it is naturally a statement over the finite etale
scheme `\mu_d`.

## 5. Fiberwise transfer principle

Let `T` be any `K`-scheme with a morphism `T->U_G`.  Base change of (1.5)
gives one cover

\[
 \mathcal X_G\times_{U_G}T\longrightarrow T
 \tag{5.1}
\]

shared by every stable class in the relevant family.

Consequently every property determined by this finite etale cover is
simultaneously shared:

- connectedness and field decomposition;
- arithmetic and geometric monodromy;
- splitting fields and deck groups;
- real signatures;
- local decomposition algebras and Frobenius cycle types;
- fiber zeta functions at good places;
- existence or failure of rational and local points.

This supplies a general use pattern.  Construct an arithmetic curve,
surface, or thin set inside one inverse plane once; it is then carried
unchanged by infinitely many stable ambient maps.  Conversely, no invariant
of the restricted finite etale cover can recover the stable class.  A
separating invariant must see how the map changes away from the common
target divisor, exactly as the boundary Fitting polygon and cubic
degree-drop components do.

## 6. The fixed quintic Hasse line

Take

\[
 G(S)=S^5-\frac32S^4+\frac32S^3-\frac54S^2+\frac9{16}S.
 \tag{6.1}
\]

Its coefficients satisfy the power-shift hypotheses.  On the line

\[
 (P,B,C)
 =\left(1,\frac{32a}{9},\frac{8a+1}{3}\right),
 \tag{6.2}
\]

write `X=S-1/2`.  The common inverse polynomial is

\[
 \mathcal E_G
 =(X^3-a)(X^2+X+1).
 \tag{6.3}
\]

Its discriminant is

\[
 \boxed{81a^2(a-1)^4.}
 \tag{6.4}
\]

Thus the line lies in the common squarefree plane away from `a=0,1`.

Let

\[
 \mathcal A=
 \{a>1:a\equiv1\pmod9,\ a\notin\mathbb Q^3,\
       p\mid a\Rightarrow p\equiv1\pmod3\}.
 \tag{6.5}
\]

For every `a` in `\mathcal A`, the fiber (6.3) has a point over every
completion of `\mathbb Q` and has no rational point.  The maps have geometric
degree five, which is the least possible degree for such a finite regular
fiber.

For `N=5`, (2.9) becomes

\[
 \boxed{7+3m.}
 \tag{6.6}
\]

Hence the maps are pairwise stably inequivalent while every one carries the
same Hasse-failing fibers at the same targets.  After applying the fixed
paper-facing determinant-one source and target normalizations, the common
target line is

\[
 \left(-1,\frac{32a}{9},\frac{8a+1}{3}\right).
 \tag{6.7}
\]

The set of targets with `a` in `\mathcal A` is infinite and therefore
Zariski dense in this line.  It is not claimed to be dense in the whole
target plane.  Its established height count is

\[
 \#\{a\in\mathcal A:H(y_a)\le B\}
 \sim
 \frac{G_3(1)}{96\sqrt\pi}\frac{B}{\sqrt{\log B}}.
 \tag{6.8}
\]

### Corollary 6.1

There are infinitely many pairwise stably inequivalent determinant-one
degree-five polynomial maps over `\mathbb Q` sharing the same Zariski-dense
family of degree-optimal Hasse-principle failures on one rational target
line.  More strongly, they share the entire squarefree inverse cover over
the ambient two-dimensional target plane containing that line.

## 7. Strength and limitations

The quantifier improvement is

\[
\begin{array}{c|c}
\text{common-fiber multiplicity}
 &\text{one prescribed target fiber is common}\\
\text{whole-plane multiplicity}
 &\text{every squarefree fiber in one target plane is common}.
\end{array}
\tag{7.1}
\]

This is strictly stronger for a fixed seed.  It does not say that one fixed
two-parameter plane contains every rank-`N` finite etale algebra.  Universal
realization still chooses or varies a presentation.

The result also does not identify the full three-dimensional generic inverse
covers of the ambient maps.  Moving `P` away from the common divisor exposes
the shifted exponents and is precisely what lets the stable boundary
invariants separate the maps.

The coordinate equality on `P=1`, the determinant identities, and selected
represented fibers are formalized in Lean.  The arbitrary-`(B,C)` literal
fiber equivalence is included in `StableGaugeFiber.lean`.

`WholePlaneStableMultiplicity.lean` now packages all target points at once.
It defines the common source-divisor functor, its morphism to the `(B,C)`
plane, the universal inverse polynomial and its discriminant, and the
discriminant-principal-open subfunctor.  For both stable families, the
pairwise equivalence is literally the identity on source-divisor points; Lean
proves that it lies over the common target morphism and is natural under
every test-algebra homomorphism.  It also proves that specializing the
universal plane equation gives the scalar inverse equation used by the
fiberwise reconstruction theorem.

`StableSeparationCertificates.lean` formalizes the numerical part of stable
separation.  It proves translation invariance and `GL_2(Z)` invariance of
normalized quadrilateral area, evaluates the four Fitting vertices as
`2N-3+(N-2)m`, proves strict injectivity in `m`, and proves injectivity of the
complete cubic boundary count `1+(n-1)=n`.  The final implications accept
only the corresponding geometric preservation statements as hypotheses.

Two geometric interfaces therefore remain outside Lean:

1. identify the written normalized ramified-stratum Fitting support with the
   four certified vertices and prove its functorial transport under stable
   polynomial left--right equivalence;
2. formalize the cubic finite-normalization boundary exhaustion and its
   stability, which supplies the certified count.

The explicit `Spec`-level isomorphism between the discriminant-localized root
quotient and the common source restriction is also not yet constructed as an
`AlgEquiv` over the localized two-variable coordinate ring.  The new
functor-of-points theorem proves that every deformation has the same
restricted morphism; the existing scalar reconstruction proves all of its
separable fibers.

## 8. Verification

Run

```bash
.venv/bin/python scripts/verify_whole_plane_stable_multiplicity.py
.venv/bin/python scripts/verify_universal_cubic_gauge_multiplicity.py
.venv/bin/python scripts/verify_universal_power_shifted_gauge_multiplicity.py
.venv/bin/python scripts/verify_infinite_hasse_keller_fibers.py
.venv/bin/python scripts/verify_multiplicative_hasse_artifact.py
cd formal/finite-etale-keller
lake build FiniteEtaleKeller.WholePlaneStableMultiplicity
lake build FiniteEtaleKeller.StableSeparationCertificates
```

The first checker verifies the whole-plane factorization for both stable
families, the cyclotomic divisibility, the exact Hasse-line factorization and
discriminant, and the quintic stable-area specialization.  The family
checkers independently replay the determinant, inverse, reconstruction, and
stable-separation data.  The two Lean builds certify the relative
discriminant-open functor and the exact separation arithmetic described
above.
