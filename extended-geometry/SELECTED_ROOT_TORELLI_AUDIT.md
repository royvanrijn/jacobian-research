# Selected-root Torelli: scope audit and low-degree tests

## Verdict

There are three different statements hidden in the phrase “the decorated
normalization determines the selected-root cover.”

1. **Normalized polynomial pencils.**  For
   \[
   E_H(W;s,t)=H(W)-sW+t
   \]
   with \(H\) in the characteristic-zero Hessian-clean normalized seed
   locus, the statement is already a theorem.  The full Fitting divisor and
   the boundary marks \(0,\infty\) leave exactly the \(N-2\) normalized
   rerootings, and the distinguished affine root sheet selects the original
   seed.  This is `D1/F2`, with the coefficient-space proof in
   [D1_F2_COEFFICIENT_SPACE_PROOF.md](D1_F2_COEFFICIENT_SPACE_PROOF.md).

2. **Weighted polynomial maps.**  The same statement is a stable-equivalence
   theorem on the exact-double, boundary-clean weighted locus.  The
   nonformal input is `WB1`: the two target boundary images, the
   discriminant-normalization stratum, the center \(0\), and the unique
   affine root sheet are selected from the map-intrinsic Zariski--Main
   package.  Stable functoriality then transports them.  This is the theorem
   proved in the
   [decorated-normalization paper](../papers/decorated-discriminant-normalization/main.tex).
   The
   [intrinsic-selector attack](INTRINSIC_SELECTOR_ATTACK.md) verifies on the
   symmetric quintic that the coarse \(\mu_3\) symmetry cannot be promoted to
   an automorphism of the full package: affine-versus-boundary membership is
   stable under arbitrary polynomial stabilization.

3. **Arbitrary finite polynomial covers with monodromy \(S_N\).**  This is
   not a consequence of the repository results and is too broad as stated.
   Full symmetric monodromy gives connectedness, simple generic inertia, and
   trivial target-fixed deck group, but it does not identify a tangent
   pencil, a twice-primitive operator, a boundary center, or a
   regular-reconstruction open.  The genuinely new problem is therefore an
   **intrinsic pencil-recognition theorem**.  Once a pencil has been
   recognized, Torelli reconstruction is formal and already proved.

The proposed “first failure”—two non-affinely-equivalent normalized
characteristic-zero seeds with the same marked Hessian/Fitting data and the
same affine sheet—cannot occur on the clean pencil locus.  The coefficient
calculation below excludes it in every degree, not only in degrees four and
five.  A first failure of the broader claim must instead be a failure of
intrinsic pencil or mark recognition.

## 1. The precise pencil Torelli theorem

Let \(k\) have characteristic zero, let \(N\geq4\), and put \(n=N-2\).
Write
\[
 H(W)=W^2(1-W)Q(W),\qquad \deg Q\leq N-3,\qquad Q(1)=1.
\]
On the clean open require:

- exact degree \(N\);
- an exact double root at \(0\) and \(n\) distinct nonzero roots;
- \(H''\) squarefree of degree \(n\), disjoint from \(0,\infty\);
- the weighted reconstruction formulas to be taken on their maximal regular
  open.

The minimal decorated object is
\[
 \mathcal T(H)=
 \bigl(\mathbb P^1;\operatorname{div}(H''),0,\infty,
       \mathcal A_{1,H}\bigr),
\]
where \(\mathcal A_{1,H}\) is the unique unramified degree-one component over
the second target boundary that lies in the affine reconstruction open
outside the degree-two zero cluster.

> **Pencil selected-root Torelli theorem.**  The functor
> \(H\mapsto\mathcal T(H)\) is faithful on the clean normalized seed locus.
> Forgetting \(\mathcal A_{1,H}\) is finite etale of degree \(N-2\) after
> shrinking to the rerooting-stable open.

Here the Fitting divisor is the full effective divisor:
\[
 \operatorname{Fitt}_0
 \Omega_{\widetilde D_H/D_H}=(H''(r)).
\]
Reduced Hessian support is enough only on the squarefree open.  The theorem
should retain the full ideal because multiplicities become essential at the
boundary.

### Coefficient proof

If two marked Hessian divisors are isomorphic, the marks \(0,\infty\) force
the normalization-coordinate change to be \(r\mapsto ar\).  Equality of
effective Fitting divisors gives
\[
 G''(r)=cH''(ar).
\]
The double-zero conditions integrate this to
\[
 G(r)=\frac{c}{a^2}H(ar).
\]
The endpoint condition \(G(1)=0\) forces \(a\) to be a nonzero root of
\(H\), and \(G'(1)=-1\) fixes the scalar.  Thus the unmarked fiber consists
exactly of the normalized rerootings
\[
 G_a(w)=-\frac{H(aw)}{aH'(a)}.
\]
There are \(N-2\) of them on the clean locus.  Under the corresponding
incidence isomorphism, the distinguished root-one sheet of \(G_a\) maps to
the root-\(a\) sheet of \(H\).  Only \(a=1\) maps the affine sheet to the
affine sheet, so \(G=H\).

This also shows why the symmetric-monodromy hypothesis is not load-bearing
for reconstruction.  The repository's universal-monodromy theorem says
that every characteristic-zero pencil \(H(W)-sW+t\) already has monodromy
\(S_N\).  Monodromy is useful for deck rigidity, but the twice-primitive
identity and the affine mark do the Torelli work.

## 2. The selected-root algebra

For a degree-\(n\) Hessian polynomial \(K\), let
\[
 J_K''=K,\qquad J_K(0)=J_K'(0)=0,\qquad
 P_K(W)=J_K(W)/W^2.
\]
On the open where \(P_K\) has \(n\) distinct nonzero roots, the algebra
\[
 \mathscr R_K=
 \mathcal O[\rho,\rho^{-1}]/(P_K(\rho))
\]
is finite etale of rank \(n\).  Its spectrum is the rerooting cover, and the
tautological root reconstructs
\[
 H_{K,\rho}(W)=
 -\frac{J_K(\rho W)}{\rho J_K'(\rho)}.
\]
This is the clean algebraic formulation of the selected root.

Calling the mark a “primitive idempotent” requires care.  At the generic
point a connected \(S_n\)-cover is a field extension and has no nontrivial
idempotents.  Primitive idempotents appear only after passing to a splitting
cover:
\[
 \mathscr R_K\otimes\widetilde K\simeq\widetilde K^n.
\]
Equivalently, after base change to the marked cover, the diagonal section in
the self-product defines the relevant idempotent.  Before splitting, the
correct intrinsic descriptions are:

- a point/section of the finite rerooting cover after base change;
- a linear factor of the universal resolvent after that base change; or
- an \(S_{n-1}\)-reduction of the \(S_n\)-torsor of ordered roots.

The resolvent-factor formulation is therefore preferable to a primitive
idempotent in the unsplit generic algebra.

## 3. When the valuation profile selects the sheet

For the weighted clean locus the canonical boundary profile over the second
target image \(Z_0\) is:

| component | residue degree | ramification | position |
|---|---:|---:|---|
| double-zero cluster | \(2\) | generically unramified over \(Z_0\) | affine |
| root-one sheet | \(1\) | unramified | affine |
| each other nonzero root | \(1\) | unramified | reconstruction boundary |

The root-one component is consequently the unique degree-one component in
the distinguished affine open outside the degree-two cluster.  Over the
other target image \(Z_\Delta\), the unique ramified boundary prime has
\((e,f)=(2,1)\).  This orders \((Z_\Delta,Z_0)\), and the closure of the
discriminant stratum over their intersection has the unique finite valuation
center \(r=0\).

These statements use more than a multiset of valuations.  They use:

- the target image of each valuation;
- incidence with the distinguished Zariski--Main open;
- the residue degrees of the components; and
- the finite discriminant-normalization stratum.

If “reduced boundary valuations” omits any of this incidence data, uniqueness
of the affine sheet does not follow.

## 4. Degrees four and five

The checker
[verify_selected_root_torelli_low_degree.py](../scripts/verify_selected_root_torelli_low_degree.py)
replays the following calculations.

### Degree four

With the extra nonzero root \(r\),
\[
 H_r(W)=-\frac{W^2(W-1)(W-r)}{1-r}.
\]
The quadratic Hessian divisor always has its affine reflection.  Its midpoint
is
\[
 m=\frac{r+1}{4},
\]
so the reflection is \(W\mapsto 2m-W\).  It preserves the intrinsic boundary
center \(0\) only when \(r=-1\).  At that centered point
\[
 H_4(W)=\frac12W^2(1-W^2),\qquad H_4''=1-6W^2.
\]
The surviving involution \(W\mapsto-W\) exchanges the primitive roots
\(1\) and \(-1\).  It therefore sends the affine root-one sheet to the
extra-root boundary sheet and is killed by the affine mark.

Thus quartic Hessian symmetry is generic, but decorated symmetry is trivial.
This is why testing only the unmarked Hessian divisor would give a misleading
quartic obstruction.

### Degree five

A generic cubic Hessian divisor has trivial affine stabilizer.  The two
nontrivial clean symmetry types are:

1. A reflection about a center \(c\neq0\).  One explicit clean seed is
   \[
   H(W)=
   -\frac{W^2(W-1)(5W^2-45W+106)}{66}.
   \]
   Its Hessian is
   \[
   H''(W)=-\frac{(W-2)(50W^2-200W+53)}{33},
   \]
   whose effective divisor is invariant under \(W\mapsto4-W\).
   The reflection sends the boundary center \(0\) to \(4\), so the boundary
   mark already removes it.

2. The centered \(\mu_3\) seed
   \[
   H_5(W)=\frac13W^2(1-W^3),\qquad
   H_5''(W)=-\frac23(10W^3-1).
   \]
   Scaling by a cube root of unity preserves \(0,\infty\) and the Hessian
   divisor, but cycles the three nonzero primitive roots.  The affine root
   sheet removes this stabilizer.  The exact valuation calculation, including
   pole order two on the two extra-root sheets, is in
   [INTRINSIC_SELECTOR_ATTACK.md](INTRINSIC_SELECTOR_ATTACK.md).

The general affine-symmetry classification shows there are no further
finite symmetry types on the clean quintic locus.  A centered involution
would require the exact-double seed exponents to occupy one residue class
modulo \(2\), which is incompatible with exact degree five; the reflection
center must therefore move the boundary mark.

## 5. Why the unrestricted conjecture needs another axiom

The phrase “canonical finite normalization” must be fixed before a Torelli
claim can be nontrivial.

- If it means the full package
  \[
  \mathcal B(F)=
  (\overline X_F\to Y,\ U\hookrightarrow\overline X_F),
  \]
  then it already determines \(F\) by restriction.  Calling this a Torelli
  theorem is tautological.
- If it means only the abstract normal source together with reduced
  divisorial data, it does not remember the finite morphism or the affine
  reconstruction open and is too weak for arbitrary covers.
- The useful intermediate object is the intrinsically selected
  discriminant-normalization stratum, its map to the target boundary image,
  the full Fitting divisor, the boundary centers, and the distinguished
  affine component.

Even this intermediate object must be recognized as a polynomial tangent
pencil.  In pencil coordinates the discriminant normalization satisfies
\[
 (S(r),T(r))=(H'(r),rH'(r)-H(r)),\qquad dT=r\,dS.
\]
Conversely, once an intrinsic affine coordinate \(r\) and target quotient
\((S,T)\) satisfy \(dT=r\,dS\), the seed is forced:
\[
 H(r)=rS(r)-T(r),\qquad H'(r)=S(r).
\]
Full symmetric monodromy alone does not imply this Legendre/tangent relation.
Nor does a Fitting divisor: for a general parametrized plane curve it records
the common vanishing of \(S'\) and \(T'\), not the identity \(T'=rS'\).

Coarse branch data are also insufficient in general.  The repository's
marked-zero-fiber Lyashko--Looijenga map has generic degree
\[
 (N-2)N^{N-3},
\]
so even in the polynomial Hurwitz category a branch-value divisor plus one
marked zero-fiber point has many preimages.  This is not a counterexample to
the full discriminant-normalization decoration, but it rules out replacing
that decoration by branch values alone.

## 6. Refined new conjecture

The reusable statement should be formulated as follows.

> **Intrinsic pencil-recognition conjecture.**  Let
> \(F:U\to Y\) be a dominant quasi-finite polynomial map of geometric degree
> \(N\ge4\) in characteristic zero.  Assume its generic monodromy is \(S_N\).
> On a clean locus suppose its Zariski--Main package intrinsically has:
>
> 1. exactly two ordered target boundary images, the first uniquely carrying
>    a ramified prime of index two;
> 2. a finite discriminant-normalization stratum with rational projective
>    normalization, a full degree-\(N-2\) Fitting divisor, and intrinsically
>    selected points \(p_0,p_\infty\);
> 3. over the second boundary, a degree-two affine cluster and a unique
>    unramified degree-one component in \(U\) outside that cluster; and
> 4. an intrinsically selected rank-two target quotient for which the
>    normalized discriminant satisfies the tangent relation \(dT=r\,dS\),
>    where \(r\) is the affine coordinate normalized by
>    \((p_0,p_1,p_\infty)=(0,1,\infty)\).
>
> Then there is a unique normalized polynomial seed \(H\), recovered by
> \(H=rS-T\), and the selected finite cover is the marked incidence cover of
> \(H(W)-sW+t\).  Stable left--right equivalence preserves this seed.

Items 1--3 are proved for the weighted family.  Item 4 is supplied there by
the tangent-map core, not reconstructed for an arbitrary cover.  Therefore
the sharp new research problem is:

\[
 \boxed{\text{derive the tangent-pencil quotient and }dT=r\,dS
 \text{ intrinsically from the decorated package.}}
\]

If that derivation is impossible, the first counterexample should consist of
two nonisomorphic \(S_N\)-covers satisfying the same boundary and Fitting
profile but inducing different tangent/contact structures on the same
normalized discriminant.  Searching for two seeds with the same marked
Hessian data will not find it: `F2` excludes that mechanism.

## 7. Concrete next tests

1. **Recognition, not reconstruction.**  Enumerate degree-four and
   degree-five finite covers having the weighted boundary degree table and
   test whether the target quotient satisfying \(dT=r\,dS\) is unique.
2. **Symmetry intersection.**  On the quartic reflection locus and the
   quintic \(\mu_2,\mu_3\) loci, compute automorphisms of the *entire finite
   stratum map*, not only of \(\operatorname{div}(H'')\), and intersect them
   with preservation of \(U\).
3. **First deformation obstruction.**  Deform the normalized plane
   parametrization by
   \((S,T)\mapsto(S+\varepsilon A,T+\varepsilon B)\) while preserving its
   Fitting divisor and boundary centers.  The quotient
   \[
   B'-rA'
   \]
   is the first-order obstruction to remaining a tangent pencil.  Determine
   whether the full finite-cover algebra and affine-sheet condition force it
   to vanish.
4. **Rerooting comparison.**  If recognition holds, the unmarked decorated
   moduli map must pull back to the universal nonzero-root algebra
   \(\mathscr R_K\), hence have stack degree \(N-2\).  On a Hessian stabilizer
   \(\mu_g\), the coarse fiber has \(n/g\) points while the stack degree stays
   \(n\); the affine mark must kill the stabilizer.
5. **Boundary of the clean locus.**  Treat collisions on the marked
   admissible-cover stack, where the selected-root cover already extends
   finitely.  Do not demand a reduced primitive idempotent at a collided
   fiber: its coarse fiber can be \(k[T]/(T^\mu)\).
