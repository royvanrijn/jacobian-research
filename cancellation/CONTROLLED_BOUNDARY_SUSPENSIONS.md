# Controlled-boundary suspensions

This note formalizes the common determinant mechanism behind the weighted,
cancellation, and quadratic-gauge constructions, proves the elementary
plane-core normal forms that can already be classified, and isolates the
first obstruction to a genuinely independent two-boundary cancellation
chart. It is a scoped starting point, not a classification of all polynomial
Keller maps.

Work over a characteristic-zero field `k`.

## 1. Boundary-cancelled incidence lemma

The reusable object is more than a Jacobian identity.

### Definition 1.1 -- boundary-cancelled incidence datum

Let `X` and `Y` be smooth affine `d`-folds, with
`\mathcal O(X)^*=k^*`, and let `F_rat:X dashrightarrow Y` be dominant. A
**boundary-cancelled incidence datum** for `F_rat` consists of the following
seven pieces.

1. **Finite root cover.** A separable degree-`N` equation
   \[
    E(\tau;y)=0
   \]
   over `k(Y)`, the finite normalization
   `pi:Xbar->Y` in `k(Y)[tau]/(E)`, and rational reconstruction functions
   with two-sided inverse identities identifying `k(X)` with that
   marked-root extension.
2. **Source chart valuation.** Dominant rational charts in a commutative
   square

   ```text
   X  -------- F_rat ------>  Y
   |                          |
   | alpha                    | beta
   v                          v
   Z  -------- Phi -------->  T
   ```

   together with the height-one valuations of `J_alpha` and
   `J_beta o F_rat`.
3. **Controlled divisor.** An irreducible equation `Delta` on `Z`, an
   integer `e>=1`, and a constant `u in k^*` such that
   \[
    J_\Phi=u\Delta^e.                                  \tag{1.1}
   \]
   In marked-root coordinates this is the derivative identity
   `partial_tau E = unit times Delta^e`.
4. **Polynomiality condition.** An explicit regularity or divisibility
   condition `P` under which every coordinate of `F_rat` extends to `X` and
   the reconstruction functions identify `X` with their common regularity
   locus in `Xbar`. The determinant ledger alone does not imply `P`.
5. **Determinant ledger.** A constant `kappa in k^*` for which
   \[
    (\Delta\circ\alpha)^eJ_\alpha
      =\kappa(J_\beta\circ F_{\rm rat})                \tag{1.2}
   \]
   in `k(X)`. Equivalently, at every height-one prime `A` of `X`,
   \[
    \operatorname{ord}_A(J_\alpha)
    +e\operatorname{ord}_A(\Delta\circ\alpha)
    =\operatorname{ord}_A(J_\beta\circ F_{\rm rat}).   \tag{1.3}
   \]
6. **Boundary normalization.** A proof that the normalization of
   `V(Delta)` maps finitely and birationally to the reduced repeated-root
   discriminant. A prime on which a reconstruction function has negative
   valuation is then a prime of `Xbar-X`, rather than an affine critical
   divisor of the Keller map.
7. **Collision fibre.** A target `y_0` for which `E(tau;y_0)` is squarefree,
   all of its roots lie in the regular-reconstruction open, and the
   reconstruction formulas give one source point for each root.

The word “condition” in item 4 is deliberate. For cancellation it is a
finite jet congruence; for the weighted and quadratic-gauge families it is a
support or denominator-clearing identity. A rational zero--pole
cancellation without this step is not a polynomial Keller construction.

### Lemma 1.2 -- boundary-cancelled incidence

For a datum as above, assume the polynomiality condition `P`. Then `F_rat`
extends to a polynomial morphism `F:X->Y` and

\[
 \boxed{J_F=u\kappa\in k^*.}                            \tag{1.4}
\]

The map `F` is the restriction of the finite root cover to its
regular-reconstruction open. The normalized controlled divisor is the
normalization of the reduced repeated-root discriminant. At the declared
collision target,

\[
 \boxed{F^{-1}(y_0)\longleftrightarrow
        \{\tau:E(\tau;y_0)=0\},}                        \tag{1.5}
\]

so the collision fibre has exactly `N` reduced points.

#### Proof

On the common domain of the square, the chain rule gives

\[
 (J_\beta\circ F)J_F=(J_\Phi\circ\alpha)J_\alpha.
\]

Substituting (1.1) and (1.2) gives `J_F=u kappa`. Polynomiality extends this
identity from the dense chart to all of `X`. The root equation and
reconstruction identify the generic function-field extension. Item 6 is
exactly the finite birational criterion for normalization. At `y_0`,
squarefreeness makes the `N` roots reduced, regular reconstruction supplies
one source point for each, and the inverse identities exclude any additional
source point. QED

### Divisor form

Taking divisors in (1.2) gives the presentation-independent ledger

\[
 \boxed{
 \operatorname{div}(\det D\alpha)
 +e\operatorname{div}(\Delta\circ\alpha)
 =F^*\operatorname{div}(\det D\beta).}                 \tag{1.6}
\]

Pullback through a rational chart is understood in the common function
field, so negative coefficients record chart poles. Conversely, divisor
equality makes the quotient in (1.2) a global unit. It is a constant when
`\mathcal O(X)^*=k^*`; on a smooth affine variety with nonconstant units,
the conclusion is only that `J_F` is a unit.

## 2. Coordinate-preserving plane cores

Let

\[
 \phi:\mathbb A^2_{w,q}\longrightarrow\mathbb A^2_{q,t},
 \qquad \phi(w,q)=(q,T(w,q)).                           \tag{5}
\]

Then

\[
 \det D\phi=-T_w.                                      \tag{6}
\]

This gives the following elementary but useful normal-form lemma.

### Proposition 2.1 -- simple section normal form

Suppose the critical divisor is the reduced section

\[
 E=(q-h(w)=0)
\]

and

\[
 \det D\phi=u(h(w)-q),\qquad u\in k^*.
\]

Choose `H` with `H'=h`.  Up to a target shear preserving `q` and a nonzero
scaling of `t`, the map is

\[
 \boxed{(w,q)\longmapsto(q,wq-H(w)).}                  \tag{7}
\]

Indeed (6) gives `T_w=u(q-h(w))`, hence

\[
 T=u(wq-H(w))+g(q).
\]

The target automorphism
`(q,t)\mapsto(q,u^{-1}(t-g(q)))` gives (7).  Thus the
universal marked-root incidence is forced by four hypotheses: a preserved
coordinate, a critical divisor which is a section, reduced ramification,
and a constant Jacobian unit along that section.

### Proposition 2.2 -- primitive form for an arbitrary boundary power

If instead

\[
 \det D\phi=uD(w,q)^r,
\]

then, up to a target shear preserving `q`,

\[
 T(w,q)=-u\int_0^wD(v,q)^r\,dv.                        \tag{8}
\]

Formula (8) classifies coordinate-preserving cores once `D` is fixed, but
does not classify the divisors `D` that admit a polynomial Keller suspension.
That suspension condition is the restrictive part of the problem.

## 3. The three established families are instances

The lemma separates the common proof from the family-specific input. The
complete dictionaries are as follows.

### 3.1 Weighted tangent family

- The finite root cover is
  \[
   E_H(W;A,B,C)=H(W)-BCW+cAC^2=0.
  \]
- The controlled divisor is `Delta=gamma=-E_H'/c`, with exponent `e=1`,
  and the tangent core has `J_Phi=-c^2 gamma`.
- For
  \[
   \alpha=(W,\gamma,C),\qquad
   \beta=(BC,cAC^2,C),
  \]
  one has
  \[
   J_\alpha=b_0x^3\gamma^2,\qquad
   J_\beta=-cC^3,\qquad C=x\gamma.
  \]
  Hence
  \[
   \gamma J_\alpha
    =-\frac{b_0}{c}(J_\beta\circ F),
  \]
  so `u=-c^2`, `kappa=-b_0/c`, and Lemma 1.2 gives
  `J_F=b_0c`.
- Admissibility of the weighted seed is the polynomiality condition.
  After suppressing the free suspension coordinate `C`, the critical graph
  `s=H'(W)` is an affine line, and
  \[
   W\longmapsto(H'(W),WH'(W)-H(W))
  \]
  is its finite birational discriminant parameterization. Thus the
  plane critical normalization is `A^1`; the full controlled divisor is its
  product with the free `C`-line.
- Whenever the specialized pencil is squarefree and its roots satisfy the
  reconstruction condition `C/E_H'(W)` regular, the fibre consists exactly
  of those roots. The explicit rational collision in every degree is the
  [all-degree rational-fibre theorem](../verified/ALL_DEGREE_RATIONAL_FIBERS.md).

### 3.2 Cancellation family

- The finite root cover is
  \[
   E_{m,r}(s;P,Q,R)
    =C\int_0^s(1-t(Q-Pt)^m)^r\,dt-R=0.
  \]
- Put
  \[
   \Delta=D=1-s(Q-Ps)^m.
  \]
  Then `partial_s E=C D^r`, so `e=r` and `u=C`.
- The rational source chart
  `alpha=(s,P,Q)` has
  \[
   J_\alpha=-A^r=-D^{-r},
  \]
  while `beta=id`. Thus `D^rJ_alpha=-1`, `kappa=-1`, and
  Lemma 1.2 gives `J_F=-C`.
- The finite congruence `L_(m,r)(h)=0` is exactly the polynomiality
  condition. At fixed `P`, the plane controlled divisor has normalization
  \[
   s=Y^{-m},\qquad Q=Y+PY^{-m},
  \]
  hence `G_m`; the full divisor carries the free `P`-parameter, and the
  exponent `r` gives ramification index `r+1`.
- At `(P,Q,R)=(1,0,0)`, the specialized root polynomial is squarefree and
  every root has `D!=0`. The reconstruction formulas therefore give the
  complete `N=r(m+1)+1` point collision fibre.

### 3.3 Root-engineered quadratic gauge

- The finite root cover is
  \[
   E_G(S;P,B,C)
    =G_P(S)-\frac{g_1}{2}(BS^2+C)=0.
  \]
- In the reciprocal chart
  \[
   D=1-SQ+PS^2,\qquad J_\alpha=D^{-1},
  \]
  and on the incidence `partial_S E_G=g_1D`. After the normalized target
  scaling, the marked-line core has `J_Phi=-2D`; hence `e=1`, `u=-2`,
  `beta=id`, and `kappa=1`. Lemma 1.2 gives `J_F=-2`.
- The coefficient-weight identity defining `G_P` is the polynomiality
  condition: it removes every apparent reciprocal-chart denominator.
  For fixed `P!=0`, the tangent-line parameter has missing points
  `S=0,infinity`, so the reduced plane discriminant normalization is `G_m`;
  over the punctured `P`-base the full normalization is `G_m x G_m`.
- At `(P,B,C)=(1,0,0)`, the root equation is `G(S)=0`. If `G` is
  squarefree, all roots have `D=G'(S)/g_1!=0`, and reconstruction gives the
  complete prescribed root fibre.

The table makes the shared mechanism and the genuine classification
variables visible:

| family | controlled exponent | chart cancellation | polynomiality gate | critical normalization |
|---|---:|---|---|---|
| weighted | `1` | distributed source/target ledger | weighted admissibility | `A^1` |
| cancellation | `r` | reciprocal source valuation `D^{-r}` | finite jet congruence | `G_m` |
| quadratic gauge | `1` | reciprocal source valuation `D^{-1}` | coefficient-weight identity | `G_m` |

Thus determinant cancellation itself is no longer the classification
problem. The remaining problem is to classify controlled divisors, normalized
critical curves with their marked valuations, and polynomial algebraizations
of their source/target ledgers.

## 4. A first independent two-boundary ansatz

The current cancellation chart has one source valuation
`A=1+xf(y)`.  To test the smallest genuine enlargement, take two distinct
functions `f_1,f_2` for which

\[
 A_1=1+xf_1(y),\qquad A_2=1+xf_2(y)
\]

are algebraically independent.  Consider the multiplicative triangular
chart

\[
 B=A_1^{b_1}A_2^{b_2}z+g(x,y),
 \quad P=A_1^{a_1}A_2^{a_2}B,
\]

\[
 s=xA_1^{d_1}A_2^{d_2},
 \qquad Q=y+sP.                                       \tag{10}
\]

This is the direct two-divisor analogue of the one-inverse-variable
triangular reconstruction skeleton.

### Proposition 4.1 -- third-divisor obstruction

The chart Jacobian is

\[
 \det\frac{\partial(s,P,Q)}{\partial(x,y,z)}
 =-A_1^{a_1+b_1+d_1-1}A_2^{a_2+b_2+d_2-1}N_{d_1,d_2}, \tag{11}
\]

where

\[
 N_{d_1,d_2}
 =(1+d_1+d_2)A_1A_2-d_1A_2-d_2A_1.                   \tag{12}
\]

The polynomial `N_(d_1,d_2)` is a monomial in `A_1,A_2` exactly for

\[
 (d_1,d_2)=(0,0),\quad(-1,0),\quad(0,-1).             \tag{13}
\]

To prove (11), replace the final coordinate `Q=y+sP` by `y`, which is a
target-triangular determinant-one operation.  Then the determinant is
`-P_z s_x`, and differentiating `s` gives (11)-(12).  Since (12) is supported
on the three distinct monomials `A_1A_2,A_1,A_2`, it is a monomial only when
two of their coefficients vanish, giving exactly (13).

The first case in (13) has no inverse factor in `s`; each of the other two
inverts only one of `A_1,A_2`.  If both exponents are genuinely active, the
Jacobian acquires the additional divisor `N_(d_1,d_2)=0`.  Hence this
smallest shared-reconstruction-variable ansatz cannot produce a clean
two-boundary determinant ledger.  A new construction must do at least one of
the following:

1. accept and control a third boundary divisor;
2. use two reconstruction variables;
3. use a nonmultiplicative source chart;
4. use a nontrivial target ledger to absorb `N_(d_1,d_2)`.

This proposition is an exact obstruction for (10), not evidence that all
two-boundary suspensions are impossible.

## 5. A realistic classification target

The first attainable theorem is an **elementary one-boundary suspension
dichotomy** under the following explicit restrictions:

1. the plane core preserves a coordinate;
2. the critical divisor is smooth, rational, and has at most two places at
   infinity;
3. the suspension square is divisor-minimal;
4. reconstruction uses one primitive variable;
5. the vertical charts are either polynomial weighted charts or rank-one
   triangular affine modifications.

Under these hypotheses, Proposition 2.1 supplies the weighted core when the
critical curve is an affine line, while the rigidity theorem for the current
cancellation ansatz supplies the punctured-line branch.  The missing global
step is a classification of divisor-minimal one-valuation rational charts of
`A^3` as rank-one affine modifications.  Without that step, exhaustiveness
must not be claimed.

## 6. Where the present method stops

The next concrete search should allow a target chart in (10) or add a second
primitive reconstruction variable and solve the resulting mixed finite-jet
conditions.  In contrast, immediately attempting to classify arbitrary
birational charts of `A^3` with prescribed Jacobian divisor is an open-ended
affine-modification problem.  That is the point at which this pathway reaches
diminishing returns unless a new structural theorem or a surviving low-degree
example supplies additional leverage.

The first target-ledger example now exists after leaving affine space:
the [three-linear-factor Cox ledger](../extended-geometry/COX_LEDGER_THREE_FACTOR.md)
keeps its rank-one boundary unit and cancels it with one primitive coordinate,
giving a constant-Jacobian finite etale morphism between affine fourfold
boundary complements.  Its reciprocal coordinate does not extend across the
collision divisor, so it does not solve the polynomial `A^4` chart problem.
The [all-arity extension](../extended-geometry/COX_LEDGER_LINEAR_FACTORS.md)
constructs the corresponding separated ledger for every `s>=3` and shows
that the unit-rank dimension bound is false without the separate-character
reconstruction requirement.
In the first arity, the
[oriented cubic target chart](../extended-geometry/ORIENTED_CUBIC_COX_CHART.md)
absorbs the primitive discriminant character polynomially, extends across
one collision divisor, and leaves the other two ordered collision branches
as distinct dicritical divisors.

The exact determinant identities and bounded exponent audit are checked by
[`verify_controlled_boundary_suspensions.py`](../scripts/verify_controlled_boundary_suspensions.py).
The later
[minimal-boundary gateway and classification note](MINIMAL_BOUNDARY_CLASSIFICATION.md)
incorporates the saturated-link, boundary-monotonicity, and reciprocal-branch
advances and is now the canonical statement of the open program.  Its
gateway conditions are proposed intrinsic hypotheses, not yet an
independently checkable classification criterion.
