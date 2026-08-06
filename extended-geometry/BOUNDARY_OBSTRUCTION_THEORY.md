# Boundary-supported obstruction, conductor descent, and finite effectivity

> **Status and scope.**  This note proves an affine version-one package for
> bounded perfect complexes and states the additional filtered hypotheses
> needed for nonlinear Kuranishi towers.  The derived boundary theorem and
> the perfect-complex conductor theorem are unconditional under their stated
> hypotheses.  The finite-effectivity theorem is conditional on explicit
> strictness and tail-elimination data.  No application below is promoted in
> mathematical status: each one must still construct its declared complex,
> class, filtration, and comparison maps.

The recurring global problem has four logically separate layers:

\[
 \boxed{
 \text{solvability on an easy open}
 \longrightarrow
 \text{boundary obstruction}
 \longrightarrow
 \text{conductor descent}
 \longrightarrow
 \text{finite effectivity}.}
\tag{1}
\]

The
[support-saturation principle](../verified/SUPPORT_SATURATION_PRINCIPLE.md)
settles the module shadow of the second layer.  The
[unified correction complex](UNIFIED_DEFORMATION_COMPLEX.md) supplies the
common linear interface, and the
[Kuranishi cutoff theorem](DEFECT_SYMBOL_APOLARITY.md#2-kuranishi-nilpotence-cutoff-theorem)
supplies one source of finite termination.  The purpose of this note is to
state precisely how these inputs fit together and, just as importantly,
where one input does not imply the next.

### Calibration from backward cubic reduction

The [backward cubic theorem](BACKWARD_CUBIC_REDUCTION.md#3-the-whole-homogenizing-line)
supplies a clean exact calibration of the boundary distinction.  Its easy
open is \(t\ne0\), where the off-diagonal collision scheme is stably

\[
 \operatorname{Coll}_{\ne\Delta}(F)\times\mathbb G_m\times\mathbb A^r.
\]

The ordinary special fiber \(t=0\) is triangular and has no off-diagonal
collision.  Nevertheless a collision arc has pole ledger

\[
 \operatorname{wt}(t,x,y)=(1,-1,-3),
\]

and becomes finite in the weighted coordinates \(X=tx\), \(Y=t^3y\).
Thus direct specialization loses the component, saturation by \(t\)
recovers its easy-open support, and the weighted Rees boundary remembers the
dehomogenized collision.  This is an explicit example of why boundary
specialization, support saturation, and Rees closure are three different
operations.

## 1. Boundary obstruction data

A **boundary obstruction datum** is a tuple

\[
 \mathfrak D=
 (A,\widetilde A,\mathfrak b,\mathfrak c,
 C^\bullet,F_\bullet,\mathcal K)
\tag{2}
\]

with the following components.

1. \(A\) is a Noetherian coefficient algebra containing all source,
   geometric, and genuine parameter variables.  Parameter specialization is
   a later base change, not part of the definition.
2. \(\mathfrak b\subset A\) is the boundary ideal and
   \[
     j:U=\operatorname{Spec}A\setminus V(\mathfrak b)
       \longrightarrow \operatorname{Spec}A
   \]
   is the easy open.
3. \(A\subseteq\widetilde A\) is a finite normalization, or more generally
   a declared finite extension, and
   \[
     \mathfrak c=\operatorname{Ann}_A(\widetilde A/A)
   \]
   is its conductor.  It is an ideal of both \(A\) and \(\widetilde A\).
4. \(C^\bullet\) is a perfect \(A\)-complex.  In the basic correction
   problem it is represented by
   \[
     C^0\xrightarrow{d_0}C^1\xrightarrow{d_1}C^2,
     \qquad d_1d_0=0,
   \tag{3}
   \]
   with gauges, corrections, and relation defects in degrees \(0,1,2\).
5. \(F_\bullet C^\bullet\) is an exhaustive separated filtration preserved
   by the differentials.  It may be polynomial degree, PBW order, pole
   order, Rees weight, a boundary valuation, or a Newton filtration.
6. \(\mathcal K\) is a base-change-compatible Kuranishi system.  If
   \(\mathcal L_{<m}\) is the scheme of lower lifts, its next obstruction is
   a section
   \[
     \overline{\mathcal O}_m
       \in\Gamma(\mathcal L_{<m},p^*H^2_m).
   \tag{4}
   \]

The datum is **linear** when only (3) is used.  It is **effective through
\(N\)** when the filtration, conductor maps, and Kuranishi operations are
defined through degree \(N\), together with a certified termination or
tail-elimination statement after \(N\).

The cohomology of (3) is

\[
 H^0=\ker d_0,\qquad
 H^1=\ker d_1/\operatorname{im}d_0,\qquad
 H^2=\operatorname{coker}d_1.
\tag{5}
\]

The complex is the common object.  It is not a universal matrix: its terms,
relations, filtration, and parameter ring depend on the problem.

## 2. Two meanings of localized vanishing

Let \(e\in H^2(C^\bullet)\).  There are two different assertions which are
often both called “\(e\) vanishes away from the boundary.”

### 2.1 Stalkwise vanishing

The first assertion is

\[
 e_{\mathfrak p}=0
 \quad\text{in }H^2(C^\bullet_{\mathfrak p})
 \quad\text{for every }\mathfrak p\notin V(\mathfrak b).
\tag{6}
\]

Since \(H^2(C^\bullet)\) is finite, (6) is equivalent to

\[
 e\in H^0_{\mathfrak b}(H^2(C^\bullet)).
\tag{7}
\]

For a three-term complex this is the familiar saturation quotient

\[
 H^0_{\mathfrak b}(\operatorname{coker}d_1)
 \simeq
 \frac{\operatorname{im}d_1:\mathfrak b^\infty}
      {\operatorname{im}d_1}.
\tag{8}
\]

This says that local primitives exist stalkwise.  It does not say that the
primitives glue on \(U\).

### 2.2 Derived vanishing on the open

The stronger assertion is that the image of \(e\) is zero in

\[
 \mathbb H^2(U,j^*\widetilde C^\bullet).
\tag{9}
\]

This includes a coherent choice of local primitives and their lower-degree
compatibilities.  The discrepancy between (6) and (9) can be carried by
Čech classes of \(H^1(C^\bullet)\), by gauge torsors one degree lower, or
by still lower cohomology.  If \(U\) is a principal affine open, or if the
relevant lower cohomology vanishes, the distinction disappears.  It must
not be suppressed for a general nonaffine \(U\).

## 3. Theorem A: derived boundary localization

Write \(R\Gamma_{\mathfrak b}\) for derived sections with support in
\(V(\mathfrak b)\).

> **Theorem A (derived boundary localization).**  Let \(A\) be Noetherian,
> let \(\mathfrak b\subset A\), and let \(C^\bullet\) be a perfect
> \(A\)-complex.  There is a functorial triangle
> \[
>  R\Gamma_{\mathfrak b}(C^\bullet)
>  \longrightarrow C^\bullet
>  \longrightarrow R\Gamma(U,j^*\widetilde C^\bullet)
>  \longrightarrow .
> \tag{10}
> \]
> If \(e\in H^2(C^\bullet)\) has zero image in (9), then
> \[
>  e\in\operatorname{im}\left(
>  H^2(R\Gamma_{\mathfrak b}(C^\bullet))
>  \longrightarrow H^2(C^\bullet)\right).
> \tag{11}
> \]
> In particular,
> \[
>  H^2(R\Gamma_{\mathfrak b}(C^\bullet))=0
>  \quad\Longrightarrow\quad e=0.
> \tag{12}
> \]

**Proof.**  Triangle (10) is the localization triangle on the affine scheme
\(\operatorname{Spec}A\).  The degree-two part of its long exact cohomology
sequence is

\[
 H^2(R\Gamma_{\mathfrak b}C^\bullet)
 \longrightarrow H^2(C^\bullet)
 \longrightarrow
 \mathbb H^2(U,j^*\widetilde C^\bullet).
\]

Exactness gives (11), and (12) follows.  \(\square\)

The theorem deliberately assumes the derived vanishing (9).  Under only
the stalkwise hypothesis (6), the module theorem gives (7), but a separate
open-set descent calculation is still required.

### 3.1 The local-cohomology spectral sequence

Because \(C^\bullet\) is bounded with finite cohomology, there is a
convergent spectral sequence

\[
 E_2^{p,q}
 =
 H^p_{\mathfrak b}(H^q(C^\bullet))
 \Longrightarrow
 H^{p+q}(R\Gamma_{\mathfrak b}(C^\bullet)).
\tag{13}
\]

For a complex concentrated in degrees \(0,1,2\), the degree-two diagonal is

\[
\begin{array}{c|c}
(p,q)&E_2^{p,q}\\ \hline
(0,2)&H^0_{\mathfrak b}(H^2)\\
(1,1)&H^1_{\mathfrak b}(H^1)\\
(2,0)&H^2_{\mathfrak b}(H^0).
\end{array}
\tag{14}
\]

Thus the following is a practical sufficient condition for (12):

\[
 H^0_{\mathfrak b}(H^2)=0,\qquad
 H^1_{\mathfrak b}(H^1)=0,\qquad
 H^2_{\mathfrak b}(H^0)=0.
\tag{15}
\]

Equivalently, it is enough to prove

\[
 \operatorname{grade}(\mathfrak b,H^q(C^\bullet))
 \ge 3-q
 \qquad(q=0,1,2),
\tag{16}
\]

with the convention that the grade of the zero module is infinite.
Regular sequences of lengths \(3-q\) in \(\mathfrak b\) give (16).

If \(H^0=H^1=0\), or more generally if the lower rows of (13) vanish in
the relevant range, then

\[
 H^2(R\Gamma_{\mathfrak b}(C^\bullet))
 \simeq H^0_{\mathfrak b}(H^2(C^\bullet)).
\tag{17}
\]

This is the precise situation in which Theorem A reduces to the
single-module saturation calculation (8).  Without a lower-row vanishing
statement, (8) is only one piece of the derived obstruction group.

### 3.2 Checkable depth criteria

For a finite module \(M\),

\[
 \operatorname{grade}(\mathfrak b,M)
 =
 \inf_{\mathfrak p\in V(\mathfrak b)\cap\operatorname{Supp}M}
 \operatorname{depth}_{A_{\mathfrak p}}M_{\mathfrak p}.
\tag{18}
\]

Consequently (16) can be certified in any of the following ways.

* Produce the required \(H^q\)-regular sequence inside \(\mathfrak b\).
* For the \(q=2\) row, exclude every associated prime of \(H^2\) which
  contains \(\mathfrak b\).  Associated-prime exclusion alone does not
  settle the \(q=1,0\) rows.
* Prove the necessary local projective-dimension bounds and apply
  Auslander--Buchsbaum over regular localizations of \(A\).
* Prove an \(S_k\) or Cohen--Macaulay condition together with the required
  dimension lower bound at every prime in
  \(V(\mathfrak b)\cap\operatorname{Supp}H^q\).
* Use an explicit exact complex, for example a certified
  Buchsbaum--Eisenbud complex, to obtain the depth bound.
* Prove
  \(V(\mathfrak b)\cap\operatorname{Supp}H^q=\varnothing\), which kills all
  local cohomology of that module.  For a finite presentation,
  \(\operatorname{Supp}H^q\) is read from its zeroth Fitting ideal.

The qualifications in these bullets matter.  Cohen--Macaulayness of an
ambient algebra does not by itself give the required depth for a particular
cohomology module.  A projective-dimension statement must be checked at the
boundary primes at which it is used.

### 3.3 Finite jets and uniform boundary torsion

Finite-jet saturation does not commute with passage to a formal limit
without a uniform torsion bound.  The following criterion isolates the
missing hypothesis.

> **Theorem A2 (uniform-exponent finite-jet convergence).**  Let \(M\) be a
> finite \(A\)-module, let \(\mathfrak m,\mathfrak b\subset A\), and put
> \[
>  M_n=M/\mathfrak m^nM,\qquad
>  T_n=H^0_{\mathfrak b}(M_n),\qquad
>  \widehat M=\varprojlim_n M_n.
> \tag{18a}
> \]
> Suppose there is one integer \(r\) such that
> \[
>  \mathfrak b^rT_n=0
>  \qquad\text{for every }n.
> \tag{18b}
> \]
> Then the natural inclusion of inverse limits identifies
> \[
>  \varprojlim_nT_n
>  =
>  H^0_{\mathfrak b}(\widehat M).
> \tag{18c}
> \]

**Proof.**  Hypothesis (18b) gives
\[
 T_n=0:_{M_n}\mathfrak b^r,
\]
because every element killed by \(\mathfrak b^r\) is
\(\mathfrak b\)-torsion.  Inverse limits commute with the kernel of the
fixed map
\[
 M_n\longrightarrow\operatorname{Hom}_A(\mathfrak b^r,M_n),
\]
or equivalently with the simultaneous kernels of a finite generating set
of \(\mathfrak b^r\).  Hence
\[
 \varprojlim_nT_n=0:_{\widehat M}\mathfrak b^r.
\]
This is contained in \(H^0_{\mathfrak b}(\widehat M)\).  Conversely,
\(H^0_{\mathfrak b}(\widehat M)\) maps into every \(T_n\); its image in the
inverse limit is killed by \(\mathfrak b^r\) because the inverse-limit
element belongs to the displayed kernel.  This proves (18c).  \(\square\)

The uniform exponent is essential.  For

\[
 A=k[x],\qquad M=A,\qquad
 \mathfrak m=\mathfrak b=(x),
\tag{18d}
\]

one has \(T_n=A/(x^n)\), the transition maps are surjective, and the least
boundary exponent is \(n\).  Thus

\[
 \varprojlim_nT_n=k[[x]],
 \qquad
 H^0_{(x)}(k[[x]])=0.
\tag{18e}
\]

Consequently a bounded collection of successful jet calculations is not a
formal saturation certificate.  A reusable finite-jet compiler must record
the transition maps and boundary-annihilation exponents, and an all-order
argument must prove a uniform bound or a genuine termination theorem.

## 4. Theorem B: perfect conductor descent

Put

\[
 A_0=A/\mathfrak c,\qquad
 \widetilde A_0=\widetilde A/\mathfrak c.
\]

The conductor square is Cartesian:

\[
\begin{array}{ccc}
A&\longrightarrow&\widetilde A\\
\downarrow&&\downarrow\\
A_0&\longrightarrow&\widetilde A_0.
\end{array}
\tag{19}
\]

Equivalently, there is an exact sequence of \(A\)-modules

\[
 0\longrightarrow A
 \longrightarrow\widetilde A\oplus A_0
 \longrightarrow\widetilde A_0
 \longrightarrow0,
\tag{20}
\]

where the last map is the difference of the two restrictions.

> **Theorem B (perfect-complex conductor descent).**  Let
> \(C^\bullet\) be a bounded complex of finite projective \(A\)-modules.
> Then the natural map
> \[
> C^\bullet
> \longrightarrow
> \operatorname{Cone}\left(
> C^\bullet_{\widetilde A}\oplus C^\bullet_{A_0}
> \longrightarrow C^\bullet_{\widetilde A_0}
> \right)[-1]
> \tag{21}
> \]
> is a quasi-isomorphism.  Equivalently,
> \[
> C^\bullet(A)
> \simeq
> C^\bullet(\widetilde A)
> \times^h_{C^\bullet(\widetilde A_0)}
> C^\bullet(A_0).
> \tag{22}
> \]

**Proof.**  Tensor (20) with each finite projective term \(C^i\).
Projectivity preserves exactness and gives

\[
 0\longrightarrow C^i
 \longrightarrow C^i_{\widetilde A}\oplus C^i_{A_0}
 \longrightarrow C^i_{\widetilde A_0}
 \longrightarrow0.
\tag{23}
\]

The last map is surjective in every degree.  Its ordinary kernel therefore
computes its homotopy fiber, and (23) identifies that kernel with
\(C^\bullet\).  \(\square\)

No additional Tor-independence hypothesis is needed for Theorem B because
the complex is represented by projective \(A\)-modules.  If an application
uses underived base change of nonflat modules, it must instead prove the
needed Tor vanishing or use derived tensor products.  Likewise, a complex
defined separately on each chart is covered by the theorem only after
proving that its terms, relation maps, and gauge maps are obtained by the
declared base changes.

### 4.1 What the mapping cone records

Let

\[
 f:C^\bullet_{\widetilde A}\oplus C^\bullet_{A_0}
 \longrightarrow C^\bullet_{\widetilde A_0}.
\]

A degree-\(n\) cocycle in \(\operatorname{Cone}(f)[-1]\) consists of

1. degree-\(n\) cocycles on \(\widetilde A\) and \(A_0\); and
2. a degree-\((n-1)\) homotopy identifying their restrictions.

For correction complexes, the homotopy is exactly the gauge-compatible
identification on the overlap.  If only the first two cocycles have been
chosen, their difference is a closed element over \(\widetilde A_0\).  Its
cohomology class is the first mismatch obstruction to choosing the
homotopy.  Different choices of homotopy form a torsor under one-lower
cohomology, modulo restrictions from the two charts.

This separates three statements which should not be conflated:

* the obstruction class vanishes on each chart;
* chosen chartwise solutions have zero mismatch on the conductor overlap;
* the mismatch admits a compatible gauge homotopy.

### 4.2 Nonlinear solutions

Suppose the correction equations, Kuranishi operations, and gauge action are
given by polynomial maps of finite \(A\)-modules and commute with all four
base changes in (19).  Then their affine solution functor preserves the
pullback (19).  With gauges retained, the solution groupoid is the
two-fiber product of the two chartwise solution groupoids over the overlap
groupoid.

This is the affine content of nonlinear conductor descent.  It requires
the operations themselves to commute with base change; a set-theoretic
agreement of solution loci or normalizations is not sufficient.

### 4.3 Rees-projective conductor descent

Filtered conductor descent has a useful sufficient condition which removes
the need for a separate overlap calculation.

> **Proposition 4.3 (Rees-projective descent).**  Give
> \(C^\bullet\) an exhaustive filtered \(A\)-module structure.  Suppose
> every \(\mathcal R_F(C^i)\) is a finite projective \(A[t]\)-module and the
> filtrations on the three conductor base changes are defined by
> \[
>  \mathcal R_F(C^\bullet_R)
>  \simeq
>  \mathcal R_F(C^\bullet)\otimes_{A[t]}R[t]
>  \quad
>  (R=\widetilde A,A_0,\widetilde A_0).
> \tag{23a}
> \]
> Then the Rees conductor square is homotopy-Cartesian.  In particular,
> conductor descent is strict in every filtration degree.

**Proof.**  Polynomial extension of (20) is exact:
\[
 0\longrightarrow A[t]\longrightarrow
 \widetilde A[t]\oplus A_0[t]\longrightarrow
 \widetilde A_0[t]\longrightarrow0.
\]
Tensor degreewise with the finite projective
\(\mathcal R_F(C^i)\).  The proof of Theorem B over \(A[t]\) gives the Rees
homotopy pullback.  Taking the degree-\(n\) part gives strict filtered
descent.  \(\square\)

The
[conductor-first node/cusp theorem](CONDUCTOR_FIRST_ONE_CHART_OBSTRUCTION.md)
shows why this is still only descent.  Its finite marked-root algebras and
discriminants descend through the conductor, but the reconstruction
coordinate cannot simultaneously have a conductor pole and be polynomial.
Thus perfect gluing data can coexist with failure of effectivity.

## 5. Theorem C: filtered finite effectivity

The Rees complex of an increasing filtration is

\[
 \mathcal R_F(C^\bullet)
 =
 \bigoplus_n F_nC^\bullet\,t^n.
\tag{24}
\]

A filtered map \(d\) is **strict** when

\[
 d(F_nC^i)=d(C^i)\cap F_nC^{i+1}
\quad\text{for every }n.
\tag{25}
\]

Under the usual good-filtration hypotheses, strictness is equivalent to
the expected Rees image/cokernel identities.  Unwanted \(t\)-torsion in a
Rees cokernel is a diagnostic for failure of strictness.  The exact
equivalence used in a certificate must state the filtration convention and
the Rees presentation; it should not be inferred from a specialization
alone.

### 5.1 A finite strictness obstruction

Entrywise comparison of two filtered matrices is not a strictness test.
Let
\[
 d=d^{(0)}+d^{(>0)}:P\longrightarrow Q
\]
be a filtered deformation whose perturbation raises filtration.  If
\(s\in\ker d^{(0)}\) is homogeneous, apply \(d\) to the unchanged lift of
\(s\), take its first nonzero filtered component, and reduce it modulo
\(\operatorname{im}d^{(0)}\).  Denote the resulting class by
\[
 \kappa(s)\in\operatorname{coker}d^{(0)}.
\tag{25a}
\]
The precise degree shift is the filtration jump of \(d^{(>0)}s\).

> **Proposition 5.1 (leading-syzygy obstruction).**  If \(s\) is the
> initial form of an actual filtered syzygy \(\widetilde s\in\ker d\), then
> \(\kappa(s)=0\).  Consequently, a nonzero class (25a) proves that the
> chosen central syzygy does not lift with that initial form.  In
> particular, entrywise positive-order perturbation of a presentation does
> not by itself identify its Rees image with the scalar extension of the
> central image.

**Proof.**  Write
\(\widetilde s=s+s_{>0}\).  In the first filtration degree above
\(d^{(0)}s=0\), the equation \(d\widetilde s=0\) reads
\[
 \operatorname{in}(d^{(>0)}s)
 +d^{(0)}\operatorname{in}(s_{>0})=0.
\]
Thus the first term is zero in
\(\operatorname{coker}d^{(0)}\).  \(\square\)

This is a necessary test, not a general converse.  Vanishing on a chosen
syzygy set still requires compatible higher corrections, and changing a
presentation or its generator lifts introduces the corresponding gauge
problem.  A complete strictness certificate may be obtained by a
terminating filtered standard-basis computation, by a lifted finite free
resolution, or directly by proving that the relevant Rees cokernel has no
\(t\)-torsion.

If \(f\in\mathfrak b\) is regular on
\(\operatorname{coker}d^{(0)}\), it is regular on the submodule generated
by the classes \(\kappa(s)\).  Hence a nonzero leading-syzygy class in that
submodule is not boundary torsion and cannot be removed by Theorem A.
It must instead be removed by changing filtered lifts or gauges.  This is
the first formal separation between a support obstruction and a
presentation-strictness obstruction inside the version-one package.

### 5.2 The bounded linear lemma

> **Lemma 5.1 (strict bounded lifting).**  Let
> \(d_1:C^1\to C^2\) be strict.  If
> \(y\in F_NC^2\) and \([y]=0\) in
> \(\operatorname{coker}d_1\), then there is
> \(x\in F_NC^1\) with \(d_1x=y\).

**Proof.**  The hypotheses say
\[
 y\in\operatorname{im}d_1\cap F_NC^2.
\]
Strictness (25) identifies this intersection with \(d_1(F_NC^1)\).
\(\square\)

Thus boundary vanishing proves existence, while Rees strictness converts
existence into a bound.  Neither assertion implies the other.

For conductor descent, the corresponding filtered hypothesis is that
(23), or equivalently (21), is strict in every degree.  Applying the Rees
functor should again give a homotopy-Cartesian square.  This guarantees that
a bounded mismatch can be killed by a bounded homotopy.  It is stronger
than unfiltered conductor descent and must be checked.

### 5.3 Why obstruction cutoff is not enough

Suppose a deformation problem has no equations and no gauges.  Every formal
series is a solution, every obstruction group vanishes in every degree, and
all boundary and conductor tests are vacuous.  A formal series with
infinitely many nonzero terms nevertheless has no bounded representative.

Therefore

\[
 \text{“no associated-graded obstructions above \(N\)”}
 \quad\not\Longrightarrow\quad
 \text{“every formal solution is bounded by \(N\)”}.
\tag{26}
\]

One also needs **tail rigidity**: high-order homogeneous solutions must be
removable by compatible gauges, or the deformation algebra itself must
terminate.

### 5.4 A version-one effectivity theorem

Let \(\operatorname{Sol}_{\le n}/\mathcal G_{\le n}\) denote the groupoid of
solutions through filtration order \(n\), and let
\(\operatorname{Sol}^{\wedge}/\mathcal G^{\wedge}\) be the inverse-limit
formal groupoid.

> **Theorem C (filtered finite effectivity).**  Let \(\mathfrak D\) be a
> boundary obstruction datum and fix \(N\).  Assume:
>
> 1. **derived boundary vanishing:** at every Kuranishi stage through \(N\),
>    the chartwise obstruction section has a coherently zero restriction to
>    \(U\), and the relevant degree of
>    \(R\Gamma_{\mathfrak b}(C^\bullet_m)\) vanishes;
> 2. **Rees strictness:** the correction, gauge, and mismatch maps used
>    through \(N\) are strict;
> 3. **filtered conductor descent:** the Rees correction complex satisfies
>    the homotopy-pullback analogue of (22), and the Kuranishi operations
>    commute with the four filtered base changes;
> 4. **finite tail elimination:** for every \(m>N\), there is a
>    filtration-preserving operation which removes an order-\(m\)
>    homogeneous solution by gauge without changing lower orders, and these
>    operations are compatible with the Kuranishi equations and conductor
>    restrictions.
>
> Here “compatible formal solution” includes coherent null-homotopies on
> the easy open and conductor-compatible chart data at every finite
> truncation.  Then every such solution is formally gauge-equivalent to a
> global solution of filtration degree at most \(N\).  If the
> tail-elimination operations are zero after finitely many steps, the gauge
> can also be chosen finite.

**Proof.**  First apply the tail-elimination operations, chartwise, to the
given compatible formal solution.  Compatibility in hypothesis 4 preserves
the Kuranishi equations, the open-set null-homotopies, and the conductor
identifications.  Iteration gives a compatible formal gauge whose
transformed chartwise solution has no terms above \(N\).

Now induct through the remaining Kuranishi stages \(m\le N\).  Theorem A
kills the global obstruction at each stage.  Lemma 5.1 chooses the
correction in the required filtration piece.  The filtered form of
Theorem B glues the two normalization/conductor corrections and their gauge
homotopy without increasing the bound.  Since the transformed chartwise
solution already satisfies every higher equation with zero tail, polynomial
base-change compatibility shows that the glued bounded solution satisfies
the full system.  If only finitely many tail operations are nonzero, their
product is a finite gauge.  \(\square\)

The fourth hypothesis is the genuinely new effectivity input.  It can be
certified without a derived-stack framework in several concrete ways.

* **Nilpotent Kuranishi ideal.**  If the relative normal ideal
  \(\mathfrak n\) satisfies \(\mathfrak n^{N+1}=0\), the completed
  Kuranishi algebra is already represented by polynomials of normal degree
  at most \(N\).  The relative Kuranishi nilpotence cutoff theorem gives
  this from a finite homogeneous initial envelope.  This controls the named
  normal filtration only; a separate comparison is needed if PBW degree,
  pole order, or another effectivity filtration is the actual target.
* **Contractible tail.**  Give filtration-preserving contracting
  homotopies for the correction/gauge complexes above \(N\), and verify
  that the nonlinear Kuranishi operations respect the resulting recursive
  elimination.
* **Artin quotient.**  Show that the effective solution torsor factors
  through a finite-dimensional filtered quotient whose associated graded
  vanishes above \(N\).
* **Geometric degree bound.**  Prove that the relevant regularity, Newton
  support, normal ideal, or boundary valuation permits no independent
  correction above \(N\), and turn that statement into the explicit
  elimination operation in hypothesis 4.

Finite generation of a Rees module, by itself, is not a finite cutoff.
Neither is vanishing of \(H^2\) above \(N\) without control of \(H^1\) and
the gauge action.

### 5.5 Formal-gauge transfer of boundary saturation

The tail-elimination hypothesis has a particularly useful module
consequence.

> **Proposition 5.5 (formal-gauge saturation transfer).**  Let \(R\) be
> Noetherian, \(\mathfrak m\subset R\), and let \(M,M_0\) be finite
> \(R\)-modules obtained functorially from two deformation objects.  Suppose
> an \(\mathfrak m\)-preserving formal gauge and base automorphism give
> \[
>  \widehat M\simeq\sigma^*\widehat M_0
> \tag{26a}
> \]
> over the \(\mathfrak m\)-adic completion.  If
> \(H^0_{\mathfrak m}(M_0)=0\), then
> \(H^0_{\mathfrak m}(M)=0\).
>
> A sufficient triangular criterion for (26a) is that, above some order
> \(N\), every compatible homogeneous correction lies in the image of the
> actual linearized gauge action, and that a gauge of order \(d-N\) removes
> its order-\(d\) image while changing only higher orders.

**Proof.**  Under the triangular criterion, choose a gauge killing the
first nonzero correction.  Its nonlinear terms have strictly higher order,
so iteration converges in the complete separated topology and gives
(26a).  Functoriality transports the completed module.

Now let \(q\in H^0_{\mathfrak m}(M)\), say
\(\mathfrak m^nq=0\).  Its image in \(\widehat M\) is zero because the
right side of (26a) has no boundary torsion.  But
\[
 (R/\mathfrak m^n)\otimes_R\widehat R
 \simeq \widehat R/\mathfrak m^n\widehat R
 \simeq R/\mathfrak m^n,
\]
so completion is faithful on every module killed by
\(\mathfrak m^n\).  Hence \(q=0\).  \(\square\)

The condition “actual linearized gauge action” is essential.  A matrix
whose columns merely have the expected degree or compatibility need not be
the differential of an integrable finite action.  A robust certificate
should derive the differential over the dual numbers and exhibit the
triangular order estimate.  The smooth universal cubic calculation does
both, and also gives an explicit first gauge lift for all 24 quartic
directions.

## 6. Base-change and functoriality ledger

Every certificate should record which of the following maps are ordinary
and which are derived.

| Operation | Safe statement | Additional check |
|---|---|---|
| localization | perfect complexes base-change; local cohomology is tested on the declared complement | whether local primitives glue on a nonaffine \(U\) |
| flat completion | finite-module \(H^0_{\mathfrak b}\) commutes with flat completion | higher local cohomology and the chosen Kuranishi operations in the required range |
| parameter base change | use \(C^\bullet\otimes_A^{\mathbf L}A'\) | Tor vanishing before replacing derived by ordinary tensor |
| normalization | use Theorem B for the base-changed perfect complex | finiteness and compatibility of the separately defined relation/gauge maps |
| good-prime reduction | exact integral matrices reduce functorially | vertical components, rank drops, and primes dividing denominators |
| associated graded | use the Rees complex | strictness and \(t\)-torsion |
| gauge quotient | retain the solution groupoid or explicit gauge homotopies | coarse orbit sets need not satisfy descent |
| change of chart | compare the complexes and Kuranishi maps, not only their zero sets | filtered and conductor compatibility |
| finite étale extension | flat base change preserves the perfect-complex calculation | normalization and the selected boundary datum must also commute with the extension |

For modular discovery, components should be labelled as:

* **horizontal** if reconstructed over the characteristic-zero parameter
  space;
* **vertical** if created only in special characteristic;
* **boundary** if supported on a declared boundary divisor; and
* **embedded** if visible only through nilpotents, finite jets, or associated
  primes not detected on reduced support.

Matching dimensions or ranks across primes is not enough.  Freeze bases,
filtrations, term orders, and boundary ideals; compare component
fingerprints; reconstruct over \(\mathbb Q\); then verify exact containments,
saturations, and identities.

## 7. Five repository prototypes

### 7.1 Cubic normalization

Take

\[
 M=\Omega_{B/A},\qquad \mathfrak b=I=(x,y,z),
\]

with a finite free presentation \(F_1\to F_0\to M\to0\) and
\(N=\operatorname{im}(F_1\to F_0)\).  The obstruction is exactly

\[
 H^0_I(M)\simeq(N:I^\infty)/N.
\tag{27}
\]

This is the module prototype at the final support step.  A regular element,
associated-prime exclusion, positive depth, or a saturated presentation
closes this layer.  The canonical-different calculation shows why this test
is consequential: cotangent saturation implies the required
annihilator/Fitting equality.  The smooth-symbol 24-parameter family is now
closed by a formal-gauge argument; singular leading symbols remain
separate.

For the smooth central symbol, the globally unit-pruned presentation has a
minimal cokernel resolution of ranks
\[
 0\longrightarrow A^7\longrightarrow A^{13}
 \longrightarrow A^6\longrightarrow\Omega_0\longrightarrow0,
\tag{27a}
\]
and \(x+y+z\) is \(\Omega_0\)-regular.  A literal lift of this resolution
would prove the desired depth statement, but it is not how the universal
theorem closes.

The filtered-syzygy frontier checker also shows why the existing two-jet
agreement is insufficient.  In the reduced \(6\)-by-\(25\) matrix, every
parameter perturbation of a nonzero central entry raises collision degree,
but applying the 24 unchanged central input syzygies to the universal
matrix leaves 12 nonzero exact remainders modulo the central image.
Proposition 5.1 therefore exposes a genuine next correction problem:
one must change the syzygies or the presentation gauges and prove that the
corrected Rees complex is strict.  Their span is a submodule of
\(\Omega_0\), so \(x+y+z\) remains regular on it.  The remainders are
therefore horizontal presentation-gauge mismatches, not boundary torsion.
This computation neither produces boundary torsion nor disproves universal
saturation.

The corrected operation acts one level earlier, on the generalized
triple-cover tensor.  For the determinant-twisted gauge differential \(G\),
the exact identity
\[
 \ker C/\operatorname{im}G\simeq\mathbb Q(-3)
\tag{27b}
\]
shows that every compatible tensor correction of collision degree at least
four is gauge.  A direct dual-number expansion verifies that \(G\) is the
actual differential of the finite gauge action, and an explicit
\(9\)-by-\(24\) matrix lifts all universal quartic directions.  Successive
homogeneous gauges therefore identify the completed universal algebra with
the completed central algebra.  Central saturation and faithful detection
of boundary-power torsion under completion then prove (27) universally for
the smooth symbol.  This is exactly the distinction required by
Proposition 5.1: the fixed central syzygies do not lift, but the nonlinear
presentation gauge removes the underlying tensor deformation before the
support theorem is applied.

The same mechanism also has a sharp, exactly computed limit.  For each
ternary-cubic symbol \(h\), let \(G_h\) be its actual determinant-twisted
gauge differential.  The graded gauge cokernel
\[
 Q_h=\ker C/\operatorname{im}G_h
\tag{27c}
\]
is \(\mathbb Q(-3)\) for the smooth symbol, but has support dimension two
for every singular squarefree symbol and support dimension three for the
double-line, triple-line, and zero symbols.  The exact quartic nongauge
dimensions are
\[
 0;\quad 2,4,4,6,6,8;\quad 11,16,24,
\tag{27d}
\]
in the smooth, six singular-squarefree, and three non-squarefree rows,
respectively.  More precisely, the singular-squarefree annihilators are
\[
 (x),\ (x^2),\ (yz),\ (y^3),\ (xyz),\ (x^3),
\tag{27e}
\]
while the last three annihilators are zero and the corresponding generic
ranks are one, two, and four.  Hence smooth is the unique orbit for which
Proposition 5.5 trivializes every higher tensor correction.  This does not
create a boundary obstruction: known singular planes are saturated.  It says that
their next theorem must apply Theorem A to a deformation-dependent
cotangent complex rather than first identifying the family formally with
its central fiber.

For the nodal row this next layer now begins with an exact first-stage
slice.  If \(\eta\) is the tensor attached to \(Z^3\), then
\[
 Q_{\mathrm{nod}}\simeq A/(x)(-3),
\tag{27f}
\]
and its quartic piece splits as a 22-dimensional gauge image plus
\(\langle y\eta,z\eta\rangle\).  Both the first-two-coordinate plane and
the full-support sum/alternating-sum plane map isomorphically to this
quotient and have saturated cotangent presentations.  This identifies the
essential first obstruction coordinates, but it does not yet satisfy the
tail-elimination hypothesis of Proposition 5.5: removing the 22 gauge
coordinates creates higher terms, while
\(\dim(Q_{\mathrm{nod}})_d=d-2\).  For the stored row-reduced gauge lift,
the complete degree-five curvature is now an exact three-component
quadratic in the basis
\(\langle y^2\eta,yz\eta,z^2\eta\rangle\).  It vanishes on the coordinate
slice but is nonzero on the dense transverse slice, even though that slice
is cotangent-saturated.  This separates normal-form curvature from boundary
torsion.  The five-dimensional ambiguity in the quartic gauge lift now has
an exact rank-four action on the six slice--gauge curvature coordinates.
Its two-dimensional cokernel gives the intrinsic linear forms
\[
 \frac34(u_3+2u_9+2u_{11}),\qquad
 \frac34(3u_6+2u_{10}+2u_{12}).
\]
It acts trivially on the three pure-gauge quadrics.  Their reduced common
zero locus is two rational planes, while the unreduced ideal has one
degree-two socle class at the origin.  Thus the lift-independence layer is
complete through degree five.

The next calculation has now reached degree six on the two reduced
pure-curvature planes.  With
\[
 (u_3,u_4,u_5,u_6)=(3\epsilon p,p,\epsilon(p+q),q),
 \qquad \epsilon=\pm1,
\]
the stored row-reduced quartic lift has an exact quadratic correction
removing its degree-five term.  The correction space has a
15-dimensional ambiguity, but its action on the four-dimensional
degree-six quotient is zero.  The resulting relative class is
\[
 \frac{27}{8}(qy+\epsilon pz)^3\eta,
\]
so it vanishes only at the origin of either reduced plane.  This is a
clean example of a higher Kuranishi stage depending on two different
lower-lift layers: independence from the degree-five correction is proved,
whereas independence from the earlier five-dimensional quartic lift is
still open.  The embedded degree-two socle and the full slice--gauge locus
also remain to be continued before testing cotangent boundary depth.

### 7.2 Relative canonical symplectic complex

The relative polynomial de Rham complex in canonical coordinates is the
positive control.  Its bounded pieces give exact
gauge-correction-defect complexes.  Here the ordinary correction and
obstruction groups \(H^1,H^2\) vanish.  Stabilizers in \(H^0\) can still
contribute to derived local cohomology, but they do not create an ordinary
degree-two obstruction in this exact control problem.  The remaining work
in a quantization application is not abstract solvability but preservation
of the Ore/boundary filtration.

### 7.3 Ritt intersections

The actual derived Hessian intersection now has a proved full bar
cotangent-descent presentation.  A perfect coefficient functor on the Ritt
face category also has a proved face-bar-to-cellular quasi-isomorphism.
What remains is the geometric coefficient-effectivity map: coherent
transport of the actual completed coefficients through the finite Ritt face
diagram and homotopy cofinality of the bar-to-face functor.

The degree-forty-two prototype explains why this hypothesis cannot be
replaced by its associated graded.  Its split cellular complex has
\(H^2=0\) in every filtration degree, but the sector--spectator extension
first becomes non-split at order three and the completed cotangent
transitivity morphism is nonzero.  Thus the obstruction lies in the
filtered coefficient/Postnikov direction, not in the topological
cohomology of the filled braid.  The extension-retaining tower is the
correct input to Theorem A; the direct-sum associated graded is not.

The degree-thirty conormal fibers and one completed degree-forty-two flag
are effective through \(H_1\).  Transport of that flag to the other five
factor-order charts, and the higher internal cotangent rows needed for a
global \(H^2\), remain open.

### 7.4 Rank-two quantization

At order \(m\), take

\[
 C^0_m\longrightarrow C^1_m\longrightarrow C^2_m,
 \qquad E_m=\operatorname{coker}d_{1,m}.
\]

Ore localization solves the easy chartwise problem.  The required order is:

1. retain the relative classical-symbol and parameter algebra;
2. construct a coherent solution on the easy open, not only stalkwise
   primitives;
3. kill the derived boundary class;
4. minimize boundary valuations modulo \(H^1\) and gauges;
5. glue the marked-root charts through the filtered conductor square;
6. prove finite tail elimination; and
7. separately certify the image subalgebra.

A fixed-symbol obstruction can terminate one branch but cannot replace the
relative Kuranishi section over the classical-symbol base.

### 7.5 Plane Jacobian boundary

The abstract finite-determinacy step is now available in
[`CONDUCTOR_JET_TRUNCATION.md`](../plane-jc/CONDUCTOR_JET_TRUNCATION.md).
On a completed branch with conductor exponent \(c_i\), an expression with
certified contact loss \(\lambda_i\) is determined in the conductor quotient
by inputs modulo \(t_i^{c_i+\lambda_i}\), and the bound is sharp in general.
One must still construct a functorial residue/conductor matching map, compile
its expression losses and band-to-normal valuation orders, form its coherent
cokernel \(M\), and define the actual class

\[
 \rho\in H^0_Z(M).
\]

Only then can boundary saturation or conductor descent be invoked.  The
truncation theorem can replace a full missing Newton tail by a finite jet
ledger.  Its dependency-sensitive form isolates the exact input-to-output
deficit, and its normal-valuation adapter derives available orders from a
certified frontier of omitted Newton support.  It still does not manufacture
the ambient module or residue class.

<!-- status-consumer: CJT1 afb70f90ff10f3d7 -->

## 8. Certificate record

The module layer is implemented in
[`support_saturation.py`](../jcsearch/support_saturation.py).  It compiles
finite presentations through Singular and records the saturation quotient,
regular-element tests, optional associated-prime decompositions,
distinguished boundary classes, finite normal jets, transition images, and
the boundary-annihilation exponent at every requested jet.  It deliberately
does not infer a uniform all-order exponent from a finite list.

The remaining derived, conductor, and Kuranishi layers should extend that
record rather than introduce a universal matrix.  A complete minimal record
is:

```text
BoundaryObstructionProblem
  coefficient_ring
  parameter_variables
  boundary_ideal
  normalization_ring
  conductor_ideal
  complex
    immutable_bases
    monomial_orders
    d0
    d1
    composition_certificate
  cohomology
    H0
    H1
    H2
    fitting_ideals
  boundary
    saturation_quotients
    regular_element_candidates
    grade_or_associated_prime_tests
    derived_local_cohomology_tests
  conductor
    four_base_changes
    mismatch_complex
    tor_or_projectivity_certificate
  filtration
    rees_matrices
    strictness_tests
    t_torsion_diagnostics
    filtered_conductor_test
  kuranishi
    lower_lift_base
    obstruction_sections
    tail_elimination_or_nilpotence_cutoff
  arithmetic
    good_primes
    component_fingerprints
    reconstructed_rational_identities
  provenance
    commands
    software_assumptions
    hashes
    verified_conditions
    unverified_conditions
```

The exact matrices, identities, and immutable bases belong in a generated
JSON artifact only when a checker exists and the reproducing command is
documented.  A modular fingerprint is discovery metadata, not a proof
certificate.

## 9. Reproduction

Run

```bash
.venv/bin/python scripts/verify_boundary_obstruction_theory.py
```

The checker uses the reusable support-saturation compiler for a regular
module, a genuine boundary-torsion module, and the tower
\(\mathbb Q[x]/(x^n)\), whose exact boundary exponents are
\(1,2,3,4,5,6\) while all transition maps are surjective.  Exact rational
matrices independently verify the node and cusp conductor kernels and their
tensor products with a rank-three free coefficient block.  A strict
filtered matrix and a minimally non-strict control verify the bounded
lifting distinction.  The result is
[`boundary_obstruction_theory.json`](../artifacts/generated-results/boundary_obstruction_theory.json).

These are exact regressions of the general proofs.  They do not certify the
open universal cubic saturation, global Ritt coefficient effectivity,
rank-two quantization algebraization, or the missing plane-JC residue
module.

## 10. Failure ledger

The package distinguishes the following outcomes.

\[
\begin{array}{ll}
\textbf{Obstruction:}&
\text{the class is nonzero already on the easy open or in global }H^2;\\[1mm]
\textbf{Support:}&
\text{the class is nonzero torsion supported on }V(\mathfrak b);\\[1mm]
\textbf{Jet convergence:}&
\text{finite-jet torsion has no uniform boundary exponent};\\[1mm]
\textbf{Open descent:}&
\text{stalkwise primitives exist on }U\text{ but do not glue there};\\[1mm]
\textbf{Conductor descent:}&
\text{normalization-chart solutions have incompatible conductor data};\\[1mm]
\textbf{Effectivity:}&
\text{a formal glued solution exists but has no bounded representative}.
\end{array}
\tag{28}
\]

Success in one row is not evidence that the next row has been settled.
Theorem A controls derived boundary support, Theorem B controls finite
normalization gluing, and Theorem C controls polynomial effectivity only
after strictness and tail rigidity have been proved.
