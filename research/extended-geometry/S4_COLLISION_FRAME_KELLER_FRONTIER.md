# The \(S_4\) collision-frame Keller frontier

## 1. Outcome and status

Work over a characteristic-zero field.  There are two distinct outcomes.

First, the group-only part of the absolute proper-monodromy problem has a
literal solution.  If \(F:\mathbb A^3\to\mathbb A^3\) is the
[foundational map](../verified/FOUNDATIONAL_GEOMETRY.md), then

\[
 H=F\circ F
\]

is an absolute polynomial Keller map with

\[
 \det DH=4,\qquad \operatorname {gdeg}(H)=9,
\]

and nonabelian proper generic inverse monodromy.  This example is
decomposable and makes no new atomicity or priority claim.

Second, factoring a depressed quartic into two quadratics gives a new
repository-level collision-frame core.  It has generic monodromy \(S_4\) in
the six-edge action, an exact discriminant ledger, a two-normal rank
obstruction, a determinant-one rational cotangent lift, and a polynomial
two-mask logarithmic lift.  It is not yet an ordinary polynomial Keller map.

This note starts from the existing
[generic-polynomial chart theorem](../KELLER_BECKMANN_BLACK_SPECIALIZATION.md),
which already realizes the \(A_4,D_5,F_{20},A_5\) point actions on smooth
affine boundary complements.  The stricter absolute search specification is:

- use a collision frame rather than another marked-one-root compiler;
- use two genuinely coupled normal directions;
- distribute conductor cancellation between source and target outputs;
- make the affine modification asymmetric and nonautomorphic; and
- finish with an actual isomorphism of the modified source with affine
  space.  Units, class group/UFD, ML, Derksen, and flexibility are only early
  gates; motivic and topological tests are also required and are still not
  sufficient without a recognition argument.

Three wording restrictions are important.

1. The six-edge action is imprimitive, so the collision cover is not
   group-theoretically atomic.  "New" below means outside the repository's
   marked-one-root compiler classes, not new in the literature and not free
   of intermediate fields.
2. The displayed \(2\times2\) target matrix is a relative Saito matrix only
   after \(p\) is inverted; globally its determinant has the extra factor
   \(p\).
3. The remaining coefficient search is finite only after a degree and
   support ansatz is fixed.  No universal degree bound is proved here.

## 2. The decomposable absolute checkpoint

The foundational map has determinant \(-2\), geometric degree three, and,
in target coordinates \((A,B,C)\), generic inverse polynomial

\[
 CT^3-2T^2+BT-2A.
\]

Its generic group is \(S_3\).  For example, the specialization
\(T^3-2T^2-1\) is irreducible with nonsquare discriminant \(-59\), while the
generic discriminant is

\[
 -4(27A^2C^2-18ABC+16A+B^3C-B^2),
\]

which is not a square.

The chain rule gives

\[
 \det D(F\circ F)=(-2)^2=4.
\]

The function-field tower

\[
 k(H)\subset k(F)\subset k(x,y,z)
\]

has successive degrees three and three, hence total degree nine.  The
intermediate field gives three blocks of three sheets.  Therefore

\[
 \operatorname {Mon}(H)\le S_3\wr S_3<S_9.
\]

Restriction to the action on blocks is the generic monodromy of the outer
copy of \(F\), hence has quotient \(S_3\).  The total group is consequently
nonabelian and proper.  The exact subgroup of the wreath product is not
needed here and is not determined by this note.

## 3. The quartic collision frame

Put

\[
\begin{aligned}
 A(T)&=T^2+aT+\frac{p+a^2-m}{2},\\
 B(T)&=T^2-aT+\frac{p+a^2+m}{2}.
\end{aligned}
\]

Then

\[
 A(T)B(T)=T^4+pT^2+qT+r
\]

for

\[
 q=am,\qquad
 r=\frac{(p+a^2)^2-m^2}{4}.
\]

Thus the coefficient map is

\[
 \Phi(p,a,m)=(p,q,r).
\]

Let

\[
 \boxed{J=m^2+2a^2(p+a^2).}
\]

Direct differentiation gives

\[
 \boxed{\det D\Phi=-\frac J2.}                     \tag{3.1}
\]

Moreover \(J=\operatorname {Res}(A,B)\).  It is the genuine ramification
factor of the factorization-to-coefficient map.

Eliminating \(m\) from \(q=am\) gives the primitive equation

\[
 \boxed{
 a^6+2pa^4+(p^2-4r)a^2-q^2=0.
 }                                                   \tag{3.2}
\]

For roots \(\alpha_1,\ldots,\alpha_4\) of the depressed quartic, choosing
\(A\) amounts to choosing one two-subset of the four roots, and \(a\) is the
negative sum of that pair.  Hence the six values of \(a\) are indexed by the
six edges of a tetrahedron.

The generic group is exactly \(S_4\).  One exact specialization suffices:
\(T^4-T-1\) is irreducible modulo \(2\), giving a \(4\)-cycle, and has
factorization type \((1,3)\) modulo \(7\), giving a \(3\)-cycle.  Its
discriminant is \(-283\), so both primes are unramified.  A transitive subgroup
of \(S_4\) containing both cycle types is \(S_4\), so specialization forces
the generic depressed-quartic group to be \(S_4\).  Its faithful action on
two-subsets has stabilizer

\[
 S_2\times S_2
\]

and gives

\[
 \boxed{S_4<S_6.}                                   \tag{3.3}
\]

Complementary edges form three blocks.  Equivalently, \(u=a^2\) satisfies

\[
 u^3+2pu^2+(p^2-4r)u-q^2=0,                         \tag{3.4}
\]

the classical cubic resolvent.  Thus this cover has a cubic intermediate
field.  The complementary-edge involution is
\((a,m)\mapsto(-a,-m)\), whose invariant ring is

\[
 k[p,a^2,am,m^2]
 \simeq k[p,u,q,v]/(q^2-uv).
\]

Thus the corresponding quotient chart is a singular quadric cone times an
affine line, not affine three-space.  This may still make the construction a
useful test of polynomial atomicity, but no atomicity theorem is claimed.

## 4. Discriminant and primitive conductor

Let \(\Delta(p,q,r)\) be the discriminant of
\(T^4+pT^2+qT+r\).  The two quadratic discriminants are

\[
 \operatorname {disc}(A)=-(a^2+2p-2m),\qquad
 \operatorname {disc}(B)=-(a^2+2p+2m).
\]

The product formula for the discriminant, or direct expansion, gives

\[
 \boxed{
 \Delta(\Phi)
 =(a^2+2p-2m)(a^2+2p+2m)J^2.
 }                                                   \tag{4.1}
\]

Differentiating (3.2) at its selected root gives

\[
 \partial_T E(T)|_{T=a}=2aJ.                        \tag{4.2}
\]

The factor \(a\) in (4.2) is not a ramification divisor of \(\Phi\).  At
\(a=0\) and \(m\ne0\), equation (3.1) is nonzero.  Instead the primitive
\(a\) assigns the same value \(0\) to the two source points with parameters
\(m\) and \(-m\).  This is a primitive-element/conductor defect.  The genuine
ramification divisor of \(\Phi\) is \(J=0\).

## 5. Tiny normal form

To avoid reusing the later normal coordinate \(\sigma\), write

\[
 c=4r-p^2.
\]

Make the source substitution

\[
 a=x,\qquad m=xy,\qquad
 p=\frac{z-2x^2-y^2}{2}.
\]

Then the collision map becomes

\[
 \boxed{
 \Theta(z,x,y)=
 \left(
 \frac{z-2x^2-y^2}{2},
 x^2y,
 x^2(z-x^2-2y^2)
 \right),
 }                                                   \tag{5.1}
\]

where the target coordinates are \((p,q,c)\), not \((p,q,r)\).  Its
Jacobian is

\[
 \boxed{\det D\Theta=-x^3z.}                        \tag{5.2}
\]

The primitive equation is

\[
 \boxed{T^6+2pT^4-cT^2-q^2=0.}                     \tag{5.3}
\]

Here

\[
 J=x^2z.
\]

The strict transform \(z=0\) is the genuine resultant ramification.  The
factor \(x=0\) is exceptional for the birational rechart
\((z,x,y)\mapsto(p,a,m)\), whose determinant is \(x/2\); it records the
primitive/conductor discrepancy and contributes the third power of \(x\) in
(5.2).  Put

\[
 \Delta_c(p,q,c)=
 \operatorname {disc}_T\left(T^4+pT^2+qT+\frac{c+p^2}{4}\right).
\]

In target coordinates \((p,q,c)\), the sextic discriminant is

\[
 \operatorname {disc}_T(E)=64q^2\Delta_c(p,q,c)^2,
\]

which independently exposes the primitive defect \(q=0\).

## 6. Why one zero-section normal cannot work

At a deepest collision point \((p,0,0)\), the derivative of \(\Phi\) has rank
one.  Consider any polynomial stabilization by one variable \(w\) with four
outputs such that the first three outputs agree with \(\Phi\) modulo \(w\).
At \(a=m=w=0\):

- the three old source columns contribute rank at most one to the first
  three rows;
- the \(w\)-column contributes at most one further rank;
- the last output row contributes at most one further rank.

Thus the \(4\times4\) derivative has rank at most three.  Consequently

\[
 \boxed{
 \text{no one-normal zero-section-preserving stabilization is Keller.}
 }                                                   \tag{6.1}
\]

This is an all-degree rank obstruction for the stated architecture.  It does
not exclude one-normal maps which move the incidence zero section.

For the two-normal search, introduce a tangential parameter \(\tau\) and
normal variables \(\sigma,n\) by

\[
 m=-2a\tau+\sigma,\qquad
 p=-a^2-2\tau^2+n.
\]

Then

\[
 \boxed{J=\sigma^2-4a\tau\sigma+2a^2n.}             \tag{6.2}
\]

If a two-mask map retains all three outputs \(p,q,r\), its derivative is
block triangular and

\[
 \det D\widetilde\Phi
 =-\frac J2\det M
\]

for the mask block \(M\).  Hence at least one coefficient output must receive
mask feedback.  The stronger assertion that two coefficient outputs must be
modified is plausible by analogy with the \(A_4\) cone, but is not proved for
this collision frame.

## 7. Rational determinant-one cotangent lift

The relative coefficient derivative is

\[
 A_0=\frac{\partial(q,r)}{\partial(a,m)}
 =\begin{pmatrix}
 m&a\\
 a(p+a^2)&-m/2
 \end{pmatrix},
 \qquad \det A_0=-J/2.
\]

The inverse-transpose mask gives

\[
 \boxed{
\begin{aligned}
 K_{\rm rat}(p,a,m,z_1,z_2)=\bigg(&p,q,r,\\
 &\frac{mz_1+2a(p+a^2)z_2}{J},\\
 &\frac{2az_1-2mz_2}{J}\bigg).
\end{aligned}
 }                                                   \tag{7.1}
\]

The derivative is block triangular, so

\[
 \boxed{\det DK_{\rm rat}=\det(A_0)\det(A_0^{-T})=1.} \tag{7.2}
\]

The masks reconstruct as \(z=A_0^T\eta\).  Adding them therefore preserves
geometric degree six and the exact \(S_4\) edge-action normal closure.  The
only pole of (7.1) is \(J\).

## 8. Relative Saito matrix and the polynomial log lift

Define

\[
 \mathsf T(p,q,r)=
 \begin{pmatrix}
 2p(p^2-4r)+9q^2&q(p^2+12r)\\
 q(p^2+12r)&16r^2-4p^2r+\frac32pq^2
 \end{pmatrix}.                                    \tag{8.1}
\]

It satisfies

\[
 \boxed{\det\mathsf T=-\frac{p\Delta}{2}.}          \tag{8.2}
\]

If the two columns act as derivations in the relative variables \((q,r)\),
then

\[
 \delta_1(\Delta)=36q\Delta,\qquad
 \delta_2(\Delta)=-4(p^2-12r)\Delta.                \tag{8.3}
\]

Thus \(\mathsf T\) is a relative logarithmic matrix globally and a relative
Saito basis over \(k[p,p^{-1}]\).  It is not a global Saito basis for
\(\Delta=0\), because \(p\) is not a unit in \(k[p,q,r]\).

There is an exact polynomial source factorization

\[
 \boxed{
 \mathsf T(\Phi)=A_0\mathsf H A_0^{\mathsf T},
 \qquad
 \mathsf H=
 \begin{pmatrix}
 a^2+2p&2am\\
 2am&-2(a^2p-2m^2+2p^2)
 \end{pmatrix}.
 }                                                   \tag{8.4}
\]

Let

\[
 \beta_{\mathsf T}(p,q,r,\eta)
 =(p,q,r,\mathsf T(p,q,r)\eta).
\]

Then

\[
 \beta_{\mathsf T}\circ K_{\rm rat}
 =(p,q,r,A_0\mathsf H(z_1,z_2)^{\mathsf T})          \tag{8.5}
\]

is polynomial, and its Jacobian is the target pullback

\[
 -\frac{p\Delta(\Phi)}2.
\]

This is the canonical two-mask polynomial log-Keller lift for the augmented
divisor \(p\Delta=0\).  Ordinary Keller factorization is exactly the problem
of replacing this determinant blowdown by a polynomial affine-space
modification without losing the six-sheet field.

## 9. Affine-source recognition gates

The following small models explain why class group, units, locally nilpotent
derivations, and flexibility are rejection tests rather than recognition
theorems.

### 9.1 The symmetric unimodular completion

The relation

\[
 mv-2au=1                                             \tag{9.1}
\]

is exactly \(SL_2\), using the matrix

\[
 \begin{pmatrix}m&2a\\u&v\end{pmatrix}.
\]

Its coordinate ring is smooth, has only constant units, and is factorial.
For example, localizing \(m\) gives the Laurent polynomial ring
\(k[m,m^{-1},a,u]\); the sole divisor over \(m=0\) is prime and principal,
so Nagata's theorem gives class group zero.

The upper and lower unipotent left actions have kernels containing
respectively the two entries of the bottom and top rows.  Hence the Derksen
invariant is the whole coordinate ring.  Their conjugates generate the left
\(SL_2\)-action, so the common invariant ring of all locally nilpotent
derivations is \(k\); the Makar--Limanov invariant is trivial.  The variety is
also flexible after base change to an algebraically closed field.

Nevertheless

\[
 [SL_2]=\mathbb L^3-\mathbb L,                       \tag{9.2}
\]

so its compactly supported Hodge--Deligne polynomial is
\((uv)^3-uv\), not \((uv)^3\).  Topologically, \(SL_2(\mathbb C)\) retracts
onto \(SU_2\simeq S^3\), so it is not contractible.  Both obstructions
survive multiplication by affine space.  Thus smoothness, UFD, constant
units, trivial ML, full Derksen invariant, and flexibility do not recognize
affine space.

### 9.2 Resultant levels

Write two general monic quadratics in difference coordinates as

\[
 A=T^2+aT+b,\qquad B=A+hT+e.
\]

Their resultant is

\[
 \rho=e^2-aeh+bh^2.                                 \tag{9.3}
\]

The level \(X_1=\{\rho=1\}\) is smooth.  On \(h\ne0\), solve for \(b\),
giving class \(\mathbb L^2(\mathbb L-1)\).  On \(h=0\), one has
\(e=\pm1\) and two affine planes.  Therefore

\[
 \boxed{[X_1]=\mathbb L^3+\mathbb L^2.}             \tag{9.4}
\]

Its Hodge--Deligne polynomial excludes \(\mathbb A^3\), before or after
ordinary stabilization.

The asymmetric level \(X_h=\{\rho=h\}\) is also smooth and has

\[
 [X_h]=\mathbb L^3.                                  \tag{9.5}
\]

Thus its motivic and Hodge tests agree with affine three-space.  They do not
recognize it.  The prime divisor \(E=(e,h)\) satisfies

\[
 \operatorname {div}(h)=2E.
\]

After inverting \(h\), the ring is
\(k[a,e,h,h^{-1}]\), whose units are \(k^*h^{\mathbb Z}\).  Nagata's theorem
shows that the class group is generated by \(E\) with \(2[E]=0\).  If \(E\)
were principal, its generator would become \(ch^n\) after localization and
would have even valuation along \(E\), a contradiction.  Hence

\[
 \boxed{\operatorname {Cl}(X_h)\simeq\mathbb Z/2,}   \tag{9.6}
\]

so \(X_h\) is not affine space despite passing the motivic test.  Polynomial
extension preserves the class group here, so ordinary stabilization does not
repair it.

### 9.3 Repeated-root Bézout slices

Let \(A,B\) be monic quadratics and \(U,V\) linears.  Put

\[
 D=B-A,\qquad W=U+V.
\]

For

\[
 UA+VB=T^3,
\]

the leading coefficient forces \(W=T-\tau\).  Write

\[
 V=\alpha(T-\tau)+x,\qquad
 D=\beta(T-\tau)+y.
\]

Then

\[
 A=\frac{T^3-VD}{T-\tau}
\]

is polynomial exactly when

\[
 xy=\tau^3.
\]

Consequently the full slice is

\[
 \boxed{\{xy=\tau^3\}\times\mathbb A^2_{\alpha,\beta},} \tag{9.7}
\]

an \(A_2\) rational double point times an affine plane.  It is singular and
therefore not affine four-space, even though its Grothendieck class is
\(\mathbb L^4\).

For \(UA+VB=T^2\), write \(V=vT+u\), \(D=hT+e\), and \(W=w\).  Eliminating
\(w=1-vh\) gives

\[
\begin{aligned}
 (1-vh)a+ve+uh&=0,\\
 (1-vh)b+ue&=0.
\end{aligned}                                      \tag{9.8}
\]

The locus

\[
 vh=1,\qquad u=e=b=0
\]

with \(a\) free is contained in the singular locus.  Hence the equal-
quadratic version of the cubic repeated-root slice does not reproduce the
smooth \((1,2)\) Bézout miracle.

## 10. Relation to the existing boundary exclusions

The collision frame does not reopen the routes already closed elsewhere.

- The [controlled two-boundary calculation](../cancellation/CONTROLLED_BOUNDARY_SUSPENSIONS.md#4-a-first-independent-two-boundary-ansatz)
  proves that a separated one-reconstruction-variable suspension introduces
  a third divisor.
- Tensor-separated cusp boundaries have a non-finitely-generated invariant
  ring by the
  [multiboundary Hilbert--14 theorem](MULTIBOUNDARY_HILBERT14_CONTROL.md#6-arbitrarily-many-independent-cusp-boundaries).
- The minimal symmetric three-boundary Cox relation retains a conductor pole
  and fails affine-space recognition by the
  [three-boundary Cox-fill obstruction](CONDUCTOR_THREE_BOUNDARY_COX_FILL_OBSTRUCTION.md).
- For the advanced \(A_4\) route, zero-section-preserving and automorphic
  incidence recharts are excluded by the
  [two-mask frontier](A4_TWO_MASK_FACTORIZATION_FRONTIER.md) and
  [normalized-boundary audit](A4_NORMALIZED_BOUNDARY_ASSEMBLY_AUDIT.md).

These results leave the same qualitative architecture: asymmetric mixed
Rees relations and a conductor ledger distributed between source and target.
They do not prove a universal no-go theorem for that architecture.

## 11. Exact remaining target

The rational map (7.1) isolates the two raw adjugate numerators

\[
 N_1=mz_1+2a(p+a^2)z_2,\qquad
 N_2=2az_1-2mz_2.                                   \tag{11.1}
\]

With the coefficient block fixed, polynomiality would require

\[
 J\mid N_1,\qquad J\mid N_2,                        \tag{11.2}
\]

which fails for independent masks.  A live construction must therefore:

1. feed the masks into at least one of (p,q,r), and test genuinely coupled
   feedback rather than an unchanged zero section;
2. replace (11.2) by the corresponding adjugate divisibilities for the
   modified coefficient block;
3. make the full Jacobian a nonzero constant;
4. prove by elimination that the function-field extension remains the
   degree-six \(S_4\) edge cover and that both masks reconstruct;
5. prove that every affine modification used as source is actually affine
   space.

For item 5, smoothness, units, class group/UFD, ML, Derksen, and flexibility
are early gates.  Hodge--Deligne, point-count, homology, and homotopy tests
must follow.  Even agreement of all these invariants is not by itself a
recognition theorem: the final step must give explicit affine coordinates or
invoke a theorem whose hypotheses prove an isomorphism with affine space.

For every fixed monomial support and degree bound, items 1--4 form a finite
mixed-incidence Gröbner problem in the conductor-adapted coordinates (6.2).
Without such a bound, the ordinary polynomial Keller factorization remains
an open, unbounded search.

## 12. Literature and novelty discipline

The sextic (3.2), or equivalently the cubic (3.4) in \(a^2\), is the
classical midpoint/resolvent construction for a depressed quartic; see
[Javier Sánchez-Reyes, *The Midpoints Between Roots Reveal the Quartic
Equation*](https://doi.org/10.1080/00029890.2020.1697589).  The relative
logarithmic terminology in Section 8 follows
[Kyoji Saito, *Theory of logarithmic differential forms and logarithmic
vector fields*](https://doi.org/10.15083/00039637).  Flexibility of
semisimple algebraic groups is included in
[Arzhantsev--Flenner--Kaliman--Kutzschebauch--Zaidenberg, *Flexible varieties
and automorphism groups*](https://arxiv.org/abs/1011.5375).

Accordingly, no historical novelty is claimed for the quartic resolvent or
the \(S_4\) edge action.  The new repository checkpoint is their exact
integration with the two-normal Keller ledger, rational cotangent lift, and
affine-source rejection tests.

## 13. Reproduction

Run

```bash
.venv/bin/python scripts/verify_s4_collision_frame_keller_frontier.py
```

The checker expands \(F\circ F\) and verifies its determinant, checks the
quartic product, primitive, Jacobian, resultant and discriminant identities,
enumerates the six-edge action, verifies the tiny normal form and rank
obstructions, checks the determinant-one rational lift and relative Saito
factorization, and replays the resultant and repeated-root obstruction
models.  Its finite-field point counts are regressions for the displayed
motivic classes; the class-group and all-degree rank arguments are the
written proofs above.
