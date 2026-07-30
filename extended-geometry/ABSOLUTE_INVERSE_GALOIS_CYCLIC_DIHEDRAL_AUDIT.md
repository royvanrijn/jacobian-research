# Absolute inverse-Galois Keller maps: cyclic and dihedral audit

## 1. Status and scope

This is the first group-by-group audit for
[Programme 3](../KELLER_INVERSE_GALOIS_PROGRAM.md).  It computes the
invariant, discriminant, orientation, unit, class-group, derivative-unit, and
first affine-modification ledgers for cyclic groups and for the natural
degree-\(n\) action of the dihedral group
\[
 D_n=\langle r,s:r^n=s^2=1,\ srs=r^{-1}\rangle
\]
of order \(2n\), with \(n\ge3\).

There are two conclusions.

1. A nontrivial cyclic group has no faithful nonregular transitive action.
   Therefore an absolute polynomial Keller map with cyclic generic inverse
   monodromy is impossible in every degree.  This is an unconditional
   consequence of the Galois case of the Jacobian theorem, not a failure of a
   particular chart.
2. The natural degree-\(n\) dihedral action survives this obstruction.  Its
   root-incidence cover is an explicit map of affine planes, and a
   derivative-unit suspension gives a determinant-minus-one morphism of
   smooth affine threefold charts.  The remaining absolute problem is exactly
   the polynomial removal of its factored boundary denominator.  A
   one-coordinate triangular removal is impossible.  In degree five, the
   first product and separated Cox fills and every affine-linear auxiliary
   coupling which retains \(u\) are also excluded; nonlinear coupled affine
   modifications which feed back into \(u\) remain open.

All invariant and factorization statements below are geometric, or hold over
a characteristic-zero field containing the required \(n\)-th cosines.  Over a
smaller arithmetic field the geometric and arithmetic monodromy must be
separated.  This is already visible for \(n=5\): the rational Dickson
polynomial has nonsquare discriminant, so its arithmetic group cannot be the
natural \(D_5\subset A_5\).

The bounded symbolic replay
[`verify_cyclic_dihedral_keller_audit.py`](../scripts/verify_cyclic_dihedral_keller_audit.py)
checks the displayed identities through degree twelve.  The uniform formulas
are proved below from \(P_n(x+y,xy)=x^n+y^n\); the bounded replay is a
regression certificate, not the proof of the all-\(n\) statement.

## 2. The common test card

For a transitive action \(G\curvearrowright G/H\), the audit card is:

| item | question |
|---|---|
| invariant ring | Are \(k[V]^G\) and \(k[V]^H\) polynomial or at least factorial? |
| discriminant | What is the reduced branch divisor and its pullback? |
| orientation | What is the sign double cover, its normalization, and its class group? |
| boundary units | What are \(\mathcal O(B^\circ)^\times/k^\times\), \(\mathcal O(X^\circ)^\times/k^\times\), and the pullback lattice? |
| class group | Are the source, target, and oriented charts factorial? |
| derivative units | Is the primitive derivative a boundary unit, and what is its lattice vector? |
| essential dimension | Is the chosen parameter base dimension-minimal? |
| affine modifications | Which pole clearings are possible, excluded, or genuinely open? |

The action, not only the abstract group, is part of the card.  A regular
action and a core-free point action have different Keller outcomes.

## 3. Cyclic groups: a complete absolute no-go

Assume first that \(k\) contains a primitive \(n\)-th root of unity and let
\(C_n\) act on \(\mathbb A^1_t\) by \(t\mapsto\zeta t\).  Then
\[
 k[t]^{C_n}=k[q],\qquad q=t^n.                       \tag{3.1}
\]
The discriminant and derivative are
\[
 \operatorname{Disc}_T(T^n-q)
 =(-1)^{n(n-1)/2}n^nq^{n-1},\qquad
 \frac{dq}{dt}=nt^{n-1}.                            \tag{3.2}
\]
On \(q\ne0\), the source and target are both \(\mathbb G_m\).  Hence
\[
\begin{aligned}
 \mathcal O(\mathbb G_{m,q})^\times/k^\times&=\mathbb Z[q],\\
 \mathcal O(\mathbb G_{m,t})^\times/k^\times&=\mathbb Z[t],\\
 q^*[q]&=n[t],
\end{aligned}                                       \tag{3.3}
\]
and both class groups vanish.

The geometric orientation cover is split when \(n\) is odd.  When \(n\) is
even its normalization is the Kummer cover \(d^2=q\), again a copy of
\(\mathbb G_m\) after the relevant constant extension.  Its class group is
zero.

The primitive derivative is a unit on the source open.  It gives the chart
\[
 (t,z)\longmapsto
 \left(t^n,\frac{z}{nt^{n-1}}\right),                \tag{3.4}
\]
whose Jacobian is one.  This is a determinant-one map
\(\mathbb G_m\times\mathbb A^1\to
\mathbb G_m\times\mathbb A^1\), not a polynomial self-map of affine space.

### Proposition 3.1 -- cyclic absolute obstruction

No nontrivial absolute polynomial Keller map in characteristic zero has
faithful cyclic generic inverse monodromy.

### Proof

Every subgroup of a cyclic group is normal.  If a transitive cyclic action
has stabilizer \(H\), its kernel is the core of \(H\), namely \(H\) itself.
Faithfulness therefore forces \(H=1\), so the action is regular and its
degree equals \(|C_n|\).  The source function field is consequently Galois
over the target function field.  The Campbell--Razar--Wright Galois-case
theorem makes a Keller map with such a normal extension invertible, contrary
to degree \(n>1\).  \(\square\)

Thus no Cox filling or affine modification can solve the absolute cyclic
case while retaining its generic cover.  Over a splitting field,
\(\operatorname{ed}(C_n)=1\); the obstruction is not parameter dimension.

## 4. The dihedral invariant and root-incidence rings

Over a splitting field, diagonalize the reflection representation as
\[
 r(x,y)=(\zeta x,\zeta^{-1}y),\qquad s(x,y)=(y,x).
                                                               \tag{4.1}
\]
Put
\[
 a=x+y,\qquad u=xy,\qquad v=x^n+y^n.                \tag{4.2}
\]
If \(H=\langle s\rangle\), then
\[
 k[x,y]^H=k[a,u],\qquad k[x,y]^{D_n}=k[u,v].         \tag{4.3}
\]
Define the Dickson power sum by
\[
 P_0=2,\quad P_1=a,\quad
 P_m=aP_{m-1}-uP_{m-2}.                              \tag{4.4}
\]
Equation (4.2) says \(P_n(a,u)=v\).  The quotient by the point stabilizer
therefore gives the affine incidence map
\[
 \pi_n:\mathbb A^2_{a,u}\longrightarrow\mathbb A^2_{u,v},
 \qquad (a,u)\longmapsto(u,P_n(a,u)).                \tag{4.5}
\]
It has degree \(n\), and its Galois closure has geometric group \(D_n\) in
the natural action on \(D_n/H\).

The incidence polynomial is
\[
 f_n(A)=P_n(A,u)-v.                                  \tag{4.6}
\]
Write
\[
 J_n(a,u)=\partial_aP_n(a,u).
\]
Differentiating \(P_n(x+y,xy)=x^n+y^n\) at fixed \(u=xy\) gives
\[
 J_n
 =n\,\frac{x^n-y^n}{x-y}.                            \tag{4.7}
\]
With source coordinates \((a,u)\) and target coordinates \((u,v)\),
\[
 \det D\pi_n=-J_n.                                   \tag{4.8}
\]

## 5. Discriminant and orientation cover

The exact discriminant is
\[
 \operatorname{Disc}_A(f_n)=
 \begin{cases}
 n^n(4u^n-v^2)^{(n-1)/2},&n\ \text{odd},\\[2mm]
 n^n(2u^m-v)^{m-1}(2u^m+v)^m,&n=2m.
 \end{cases}                                         \tag{5.1}
\]
In either parity the reduced geometric branch support is
\[
 \mathcal D_n:\quad v^2-4u^n=0.                     \tag{5.2}
\]
Its pullback has the decisive factorization
\[
 P_n(a,u)^2-4u^n
 =(a^2-4u)\left(\frac{J_n(a,u)}n\right)^2.           \tag{5.3}
\]
The factor \(a^2-4u\) is unramified at its generic point: it records a source
point whose stabilizer is already the chosen reflection \(H\).  The factors
of \(J_n\) are the actual ramification divisors.

Geometrically, the sign orientation cover has the following normalization.

| parity | square class of (5.1) | normalized orientation chart | geometric class group |
|---|---|---|---:|
| \(n\equiv1\pmod4\) | constant | split/trivial | \(0\) on each component |
| \(n\equiv3\pmod4\) | \(4u^n-v^2\) | \(XY=u^n\), type \(A_{n-1}\) | \(\mathbb Z/n\) |
| \(n=2m\) | one of \(2u^m-v\), \(2u^m+v\) | \(\mathbb A^2\) | \(0\) |

For odd \(n\equiv3\pmod4\), the change
\(X=d+v,\ Y=d-v\), followed by harmless nonzero scalar rescaling, gives
\(XY=u^n\).  The stated class group is the standard geometric divisor class
group of the normal \(A_{n-1}\) surface.

The orientation cover is diagnostic, not a free Keller simplification.
Base change to it replaces \(D_n\) by the kernel of the sign character.  In
particular, one must not claim preservation of exact dihedral monodromy after
silently adjoining the orientation.

## 6. Boundary unit lattice and class groups

Let
\[
 B_n^\circ
 =\operatorname{Spec}k[u,v,(v^2-4u^n)^{-1}]
                                                               \tag{6.1}
\]
and let \(X_n^\circ=\pi_n^{-1}(B_n^\circ)\).  By (5.3),
\[
 X_n^\circ
 =\operatorname{Spec}
 k[a,u,\{(a^2-4u)J_n\}^{-1}].                        \tag{6.2}
\]
Both rings are localizations of polynomial UFDs, so
\[
 \operatorname{Cl}(B_n^\circ)
 =\operatorname{Cl}(X_n^\circ)=0.                   \tag{6.3}
\]

Over a fully split reflection field, (4.7) becomes
\[
\frac{J_n}{n}=
\begin{cases}
\displaystyle
\prod_{j=1}^{(n-1)/2}
\left(a^2-(2+\zeta^j+\zeta^{-j})u\right),
&n\ \text{odd},\\[4mm]
\displaystyle
a\prod_{j=1}^{n/2-1}
\left(a^2-(2+\zeta^j+\zeta^{-j})u\right),
&n\ \text{even}.
\end{cases}                                         \tag{6.4}
\]
Consequently
\[
 \operatorname{rank}
 \frac{\mathcal O(X_n^\circ)^\times}{k^\times}
 =1+\left\lfloor\frac n2\right\rfloor,               \tag{6.5}
\]
where the extra generator is \(a^2-4u\), while
\[
 \operatorname{rank}
 \frac{\mathcal O(B_n^\circ)^\times}{k^\times}
 =
 \begin{cases}
 1,&n\ \text{odd},\\
 2,&n\ \text{even}.
 \end{cases}                                         \tag{6.6}
\]
Over a nonsplitting arithmetic field, Galois-conjugate factors in (6.4)
coalesce and the rank is the number of irreducible boundary primes over that
field.  Formula (5.3), rather than the geometric rank alone, is the
base-independent pullback ledger.

For odd \(n\), if the source unit basis is
\[
 [C],[R_1],\ldots,[R_{(n-1)/2}],
 \quad C=a^2-4u,\quad J_n/n=\prod R_i,
\]
then
\[
 \pi_n^*[\mathcal D_n]=(1,2,\ldots,2),\qquad
 [J_n]=(0,1,\ldots,1).                               \tag{6.7}
\]
For even \(n\), the two target branch components partition the ramification
factors in (6.4); one pullback also contains \(C\) with multiplicity one,
and every ramification factor occurs with multiplicity two in exactly one of
the two pullbacks.  Thus the derivative is again the half-ramification
vector, but the unramified \(C\)-color prevents it from being half of a
target pullback inside the integral unit lattice.

## 7. Derivative unit and Keller chart

On \(X_n^\circ\), \(J_n\) is a unit.  The universal derivative suspension is
\[
 \widehat\pi_n:
 X_n^\circ\times\mathbb A^1_z
 \longrightarrow B_n^\circ\times\mathbb A^1_Z,
 \qquad
 (a,u,z)\longmapsto
 \left(u,P_n(a,u),\frac z{J_n(a,u)}\right).           \tag{7.1}
\]
Equations (4.8) and (7.1) give
\[
 \det D\widehat\pi_n=-1.                             \tag{7.2}
\]
Replacing \(Z\) by \(-Z\) gives determinant one.  The function-field
extension and its dihedral Galois closure are unchanged, and every regular
root fiber is complete.

The derivative exponent vector is
\[
 [J_n]=(0,1,\ldots,1)                                \tag{7.3}
\]
in the split source boundary basis.  Thus the chart correction exists
precisely because the full ramification product is a unit after deleting the
target discriminant.

## 8. Essential dimension

Over an algebraically closed characteristic-zero field,
\[
 \operatorname{ed}(C_n)=1,\qquad
 \operatorname{ed}(D_n)=
 \begin{cases}
 1,&n\ \text{odd},\\
 2,&n\ \text{even}.
 \end{cases}                                         \tag{8.1}
\]
This is the cyclic/odd-dihedral classification in
[Buhler--Reichstein](https://doi.org/10.1023/A:1000144403695).
Over an arithmetic base the value depends on the available roots of unity
or cosines.

The two-parameter invariant chart (4.5) is therefore dimension-minimal for
even \(n\), but not for odd \(n\).  Over a field containing the \(n\)-th
cosines, Hashimoto--Miyake give one-parameter generic polynomials for odd
dihedral groups; Ledet's
[dihedral generic-polynomial construction](https://www.math.ttu.edu/~aledet/papers/dihedral.pdf)
also treats the divisible-by-four cases over appropriate cosine fields.
Auditing those compressed torsors is the next parameter-minimal branch.
The present Dickson chart remains preferable for the first boundary
calculation because both invariant rings and the complete pullback ledger
are polynomial and uniform in \(n\).

## 9. Affine-modification candidates and exclusions

### 9.1 One triangular coordinate is closed

Consider a polynomial triangular suspension which retains the two cover
coordinates:
\[
 (a,u,z)\longmapsto
 \bigl(u,P_n(a,u),A(a,u)z+B(a,u)\bigr).               \tag{9.1}
\]
Its Jacobian is
\[
 -J_n(a,u)A(a,u).                                    \tag{9.2}
\]
For \(n\ge2\), \(J_n\) is nonconstant.  Hence (9.2) cannot be a nonzero
constant with polynomial \(A\).  The rational choice \(A=J_n^{-1}\) is
exactly (7.1).  This closes direct one-coordinate pole clearing, but not
coupled modifications of two or more outputs.

### 9.2 Orientation/Cox candidates

- For odd \(n\equiv3\pmod4\), the oriented target is the toric
  \(A_{n-1}\) surface \(XY=u^n\), with class group \(\mathbb Z/n\).
  Its Cox cover is explicit, but it changes the base and the monodromy must
  be descended back through the sign involution.
- For even \(n\), the normalized orientation surface is affine space.
  Nevertheless the base change replaces the exact \(D_n\)-cover by the
  sign-kernel cover, so the affine orientation alone does not solve the
  stated problem.
- The source lattice has more colors than the target lattice by
  \[
  \begin{cases}
  (n-1)/2,&n\ \text{odd},\\
  n/2-1,&n\ \text{even}.
  \end{cases}
  \]
  Any Cox or affine-modification completion must account for those colors
  without merging inverse sheets or enlarging arithmetic monodromy.

### 9.3 Surviving absolute search

The first viable unrestricted ansatz is a coupled two-output modification:
new source-dependent masks must enter at least two target coordinates so
that cross terms can cancel \(J_n\).  It must pass four independent tests:

1. polynomiality on all of affine space;
2. constant full Jacobian, not only residue Jacobian;
3. reconstruction of \(P_n(A,u)-v\) on a dense target open;
4. unchanged geometric and arithmetic \(D_n\)-closure.

No such absolute map is constructed here.  For \(D_3=S_3\), \(D_4\), and
\(D_5\), these equations are the next bounded affine-modification search.

For \(D_5\), the
[dedicated modification frontier](ABSOLUTE_DIHEDRAL_D5_MODIFICATION_FRONTIER.md)
now performs the first such precomputation.  It factors
\(J_5/5=R_+R_-\), proves that the product Cox fill is singular with class
group \(\mathbb Z\), and proves that the separated two-color fill is a
factorial but singular quadric cone.  More generally, retaining \(u\) and
making all remaining outputs affine-linear in any number of auxiliary
coordinates forces generic degree one.  An arbitrary nonlinear thickening
which retains the incidence as a zero section also keeps \(J_5\) as a
Jacobian divisor.  Thus the first surviving \(D_5\) ansatz must modify \(u\),
depend nonlinearly on the new coordinates, and preserve the degree-five
field only through a nontrivial birational elimination.

The
[nonlinear obstruction classification](D5_NONLINEAR_MODIFICATION_OBSTRUCTION_CLASSIFICATION.md)
separates that residual stratum into eight necessary gates.  In particular,
the branch-supported log ledger has the unique order pattern
\((m,2m-1,2m-1)\), and the primitive diagonal realization still leaves
the pole \(1/(R_+R_-)\).
The subsequent
[canonical two-mask blowdown](D5_TWO_MASK_BLOWDOWN_OBSTRUCTIONS.md)
has determinant \(\Delta\), but its two inverse adjugate divisibilities
exclude all constant-linear assemblies.  A generic genus-two-fibre argument
upgrades the base-mixing obstruction to polynomial automorphic recharts of
every degree, and the first nonautomorphic normalized-cusp chart fails by a
contraction-divisor mismatch.

## 10. Low-degree cards

| group/action | incidence polynomial | discriminant | source unit rank | orientation class group | absolute status |
|---|---|---|---:|---:|---|
| \(C_n\), regular degree \(n\) | \(T^n-q\) | \(\sim q^{n-1}\) | \(1\) | \(0\) | impossible |
| \(D_3=S_3\), degree \(3\) | \(A^3-3uA-v\) | \(27(4u^3-v^2)\) | \(2\) | \(\mathbb Z/3\) | open |
| \(D_4\), degree \(4\) | \(A^4-4uA^2+2u^2-v\) | \(256(2u^2-v)(2u^2+v)^2\) | \(3\) | \(0\) | open |
| geometric \(D_5\), degree \(5\) | \(A^5-5uA^3+5u^2A-v\) | \(5^5(4u^5-v^2)^2\) | \(3\) | \(0\) | open over a split base |

The \(D_5\) row is deliberately labeled geometric.  Over \(\mathbb Q(u,v)\)
its displayed discriminant is not a square, whereas the natural \(D_5\)
action lies in \(A_5\).  Exact arithmetic monodromy therefore requires a
cosine-field model or a different rational generic polynomial.

## 11. Next computations

1. For \(D_5\), search only nonlinear two-mask ansatzes which modify \(u\);
   product/separated Cox fills and affine-linear masks are closed.
2. Run the coupled quadratic two-mask search first for \(D_3\) and \(D_4\).
3. Audit a one-parameter odd-dihedral generic polynomial and compare its
   boundary lattice with the two-parameter Dickson chart.
4. Compute the arithmetic normalizer of the rational Dickson family for
   each \(n\), rather than reporting only geometric \(D_n\).
5. For the odd \(A_{n-1}\) orientation surface, compute the Cox pullback of
   every ramification color and test whether the class
   \([J_n]\) becomes a primitive character.
6. Treat \(V_4\) next as an imprimitive degree-four action, not as a regular
   cyclic-style cover; the action choice is decisive.
