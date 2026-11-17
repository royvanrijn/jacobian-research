# The Huq--Kuruvilla map, Mondello's plane theorem, and repository extensions

This note begins with Irit Huq-Kuruvilla's characteristic-two threefold
counterexample to the separable Jacobian conjecture.  The threefold map, its
determinant-one property, its three-point collision, its generic inverse
cubic, and the rational reconstruction on the generic sheet are due to
Huq-Kuruvilla:

> I. Huq-Kuruvilla, *An Explicit Characteristic-\(2\) Counterexample to the
> Separable Jacobian Conjecture*,
> [arXiv:2607.20968](https://arxiv.org/abs/2607.20968) (2026).

Romy Mondello subsequently exposed a preserved polynomial coordinate after
permuting the source and target coordinates.  Its zero fiber is the explicit
plane map

\[
 \begin{aligned}
 P_2&=x+x^2y+x^4+x^6y^2,\\
 Q_2&=y+x^5+x^6y+x^7y^2+x^8y^3.
 \end{aligned}                                      \tag{0.1}
\]

The primary source for the plane theorem is:

> R. Mondello, *A Dimension-Two Counterexample to the Separable Jacobian
> Conjecture in Characteristic Two*,
> [arXiv:2608.02634v1](https://arxiv.org/abs/2608.02634) (2026).

Mondello proves the theorem over \(k=\overline{\mathbb F}_2\): the displayed
map has determinant one and is etale; its target coordinates are
algebraically independent; its function-field extension is separable of
exact degree three; its geometric generic fibre has three points; and the
three rational points \((0,1),(1,0),(1,1)\) collide at \((0,1)\).  The paper
also derives the map as the invariant fibre of a skew-product presentation of
Huq--Kuruvilla's threefold construction.

The paper's proof architecture includes the hidden variables, recovery
formulas, actual-target irreducibility argument, and geometric-generic-fibre
bridge reproduced in Section 8 below.  Its appendix records a Lean kernel
check and a computationally separate Harmonic Aristotle replay.  Those are
formal and computational assurance signals; neither constitutes independent
human peer review, and neither settles novelty, exposition, or the later
boundary and lifting claims in this repository.

The provenance boundary is therefore sharp:

- **Contained in Mondello:** the plane reduction, determinant one, collision,
  etaleness, algebraic independence, hidden cubic, actual-target
  irreducibility, exact separable degree three, geometric generic-fibre
  interpretation, skew-product derivation, and the stated Lean/replay scope.
- **Added by this repository:** the arbitrary-characteristic-two base-field
  corollary in Proposition 8.1; the explicit finite normalization and the
  distinction between reconstruction boundary and primitive-order conductor;
  conductor, different, nodal intersections, Frobenius boundary map,
  \(S_3\) monodromy and wild inertia; and the modulo-four obstruction,
  stabilized finite-Witt tower, and nonpolynomial inverse-limit phenomenon.

Accordingly `HKM2` is an external theorem with an internal exact replay, not
a repository theorem.  The threefold normalization in Sections 2--7 and the
plane extensions following Proposition 8.1 are repository results.  No
priority claim is made here for either external counterexample.

Throughout, \(k\) is a field of characteristic \(2\).

## 1. Huq--Kuruvilla's map and inverse cubic

Put

\[
\begin{aligned}
 P&=x+x^2y,\\
 Q&=y+xz+x^2yz,\\
 R&=z+x^2z^2.
\end{aligned}                                                   \tag{1.1}
\]

Huq-Kuruvilla proves

\[
 \det\frac{\partial(P,Q,R)}{\partial(x,y,z)}=1
 \quad\text{and}\quad
 F(0,1,0)=F(1,1,0)=F(1,1,1)=(0,1,0).              \tag{1.2}
\]

Use the triangular target coordinates

\[
 U=P,\qquad V=Q+PR,\qquad W=R,                       \tag{1.3}
\]

so that \(Q=V+UW\).  The generic inverse polynomial found in the cited
paper is

\[
 \boxed{D(T)=UT^3+T^2+VT+W.}                         \tag{1.4}
\]

If

\[
 t=\frac{1+x^2z}{x},\qquad a=D'(t)=Ut^2+V,
\]

then Huq-Kuruvilla's reconstruction is

\[
 \boxed{x=a^{-1},\qquad z=a(t+a),\qquad y=Q+Uz.}      \tag{1.5}
\]

The polynomial \(D\) is irreducible over \(k(U,V,W)\), so the function-field
degree is three.  Formula (1.5) also shows directly that it is separable.

## 2. Discriminant and the full projective collision

For a cubic \(AT^3+BT^2+CT+D\), the discriminant in characteristic two is

\[
 B^2C^2+A^2D^2.
\]

Consequently

\[
 \boxed{\operatorname{disc}_T(D)=V^2+U^2W^2
        =(V+UW)^2=Q^2.}                              \tag{2.1}
\]

Thus the discriminant scheme is the doubled hyperplane \(V(Q^2)\), while
its reduced support is the smooth hyperplane

\[
 \mathcal D=V(Q).
\]

The inverse polynomial must be homogenized to retain roots at infinity:

\[
 D_h(T,S)=UT^3+T^2S+VTS^2+WS^3.                     \tag{2.2}
\]

At the collision target \((P,Q,R)=(0,1,0)\), equivalently
\((U,V,W)=(0,1,0)\),

\[
 D_h(T,S)=TS(T+S).
\]

Its three projective roots \(0,1,\infty\) reconstruct respectively to

\[
 (1,1,1),\qquad(1,1,0),\qquad(0,1,0).
\]

The root at infinity is why the affine specialization of (1.4), which drops
to degree two at \(U=0\), does not by itself display the complete collision.

## 3. The finite normalization

Define

\[
 \overline X=
 V(D_h)\subset\mathbb P^1_{[T:S]}\times\mathbb A^3_{U,V,W},
 \qquad
 \pi:\overline X\longrightarrow\mathbb A^3_{U,V,W}. \tag{3.1}
\]

Every fiber is a degree-three divisor on \(\mathbb P^1\): the coefficient
of \(T^2S\) is \(1\), so the binary cubic never vanishes identically.
Therefore \(\pi\) is finite flat of degree three.

The two standard charts are explicit:

\[
\begin{array}{ll}
S\ne0,\ t=T/S:
 &X_t=\operatorname{Spec}
   k[U,V,W,t]/(Ut^3+t^2+Vt+W),\\[2mm]
T\ne0,\ s=S/T:
 &X_s=\operatorname{Spec}k[V,W,s],\quad
   U=s+Vs^2+Ws^3.
\end{array}                                                   \tag{3.2}
\]

They are glued by \(s=t^{-1}\).  The projective hypersurface is smooth:
the partial derivatives with respect to \(U\) and \(W\) are \(T^3\) and
\(S^3\), which cannot vanish simultaneously on \(\mathbb P^1\).
It is integral as well, for its two integral affine charts have a common
dense overlap.  Hence \(\overline X\) is normal.  It has the same function
field as the source by (1.4)--(1.5), so

\[
 \boxed{\overline X=
 \operatorname{Norm}_{\mathbb A^3_{U,V,W}}k(x,y,z).} \tag{3.3}
\]

This projective incidence is therefore the finite normalization, rather
than merely a convenient compactification of an affine inverse polynomial.

## 4. Reconstruction open and its pole

On \(X_t\), put

\[
 a=Ut^2+V.
\]

The projection has Jacobian

\[
 \det\frac{\partial(U,V,W)}{\partial(U,V,t)}
 =\frac{\partial(Ut^3+t^2+Vt)}{\partial t}
 =a.                                                   \tag{4.1}
\]

The reconstruction formulas are

\[
 x=a^{-1},\qquad z=a(t+a),\qquad
 y=(V+UW)+Ua(t+a).                                    \tag{4.2}
\]

Only \(x\) has a pole along \(a=0\).  The apparent restriction to the
\(t\)-chart causes no loss.  On \(X_s\), set

\[
 \delta=1+Ws^2.
\]

Then

\[
 x=\frac{s}{\delta},\qquad z=W\delta,\qquad
 y=(V+UW)+UW\delta,                                   \tag{4.3}
\]

and on the overlap \(a=\delta/s\).  The divisor \(a=0\) is wholly contained
in \(X_t\), whereas (4.3) regularly includes the simple root at infinity.
It follows that

\[
 \boxed{\mathbb A^3_{x,y,z}\simeq\overline X\setminus E,\qquad
 E=V(a).}                                             \tag{4.4}
\]

In particular, \(E\) is the entire normalization boundary and the exact
reconstruction-pole divisor.

On the source chart \(\alpha=(U,V,t)\), direct differentiation gives

\[
 J_\alpha=x=a^{-1}.
\]

Together with (4.1), this is the determinant ledger

\[
 \boxed{a\,J_\alpha=1.}                               \tag{4.5}
\]

Thus the constant Jacobian does have a marked-root zero--pole explanation.
It is a direct reciprocal cancellation, even though the geometry of the
controlled divisor differs from the established reciprocal family.

## 5. The wild boundary and the two components over \(Q=0\)

The derivative divisor has the parameterization

\[
 E:\quad V=Ut^2,\qquad W=t^2.
\]

Hence

\[
 E\simeq\mathbb A^2_{U,t},\qquad
 \pi|_E:(U,t)\longmapsto(U,W=t^2)\in\mathcal D.        \tag{5.1}
\]

After suppressing the free \(U\)-coordinate, its critical normalization is
\(\mathbb A^1_t\), but its map to the reduced discriminant is Frobenius:

\[
 k(U,W)\ \subset\ k(U,t),\qquad t^2=W.                \tag{5.2}
\]

It is generically purely inseparable of degree two, not birational.

The complete factorization over the discriminant is most transparent in
the two charts:

\[
\begin{aligned}
 Q=V+UW&=(Ut^2+V)(1+Ut)=a(1+Ut) &&\text{on }X_t,\\
 Q&=(1+Ws^2)(V+Ws)=\delta(V+Ws) &&\text{on }X_s.
\end{aligned}                                         \tag{5.3}
\]

Thus, generically,

\[
 \pi^*\mathcal D=E+A.                                 \tag{5.4}
\]

The companion component \(A\) is \(V+Ws=0\) on \(X_s\).  There
\(U=s\) and \(V=UW\), so \(A\to\mathcal D\) is an isomorphism.
At \(E\), by contrast,

\[
 e(E/\mathcal D)=v_E(Q)=1,\qquad
 [k(E):k(\mathcal D)]_{\mathrm{insep}}=2.              \tag{5.5}
\]

The Jacobian in (4.1) vanishes simply, so the different exponent is one.
This is wild defectless ramification of the ``fierce'' or
residually-inseparable kind: the ramification index remains one and the
missing factor of two occurs in the residue extension.  The exponent two
in (2.1) is the norm-level shadow of this degree-two radicial different.

## 6. Monodromy

Irreducibility of (1.4) makes the geometric monodromy a transitive subgroup
of \(S_3\), hence either \(C_3\) or \(S_3\).

The extension is not Galois.  Over the prime \(Q=0\), (5.3)--(5.5) give one
prime with separable residue degree one and another with inseparable residue
degree two.  A cyclic cubic Galois action would act transitively on primes
over a base prime and preserve their ramification and residue data, which
is impossible here.  Therefore

\[
 \boxed{G_{\mathrm{geom}}=S_3.}                       \tag{6.1}
\]

The same argument applies after extending the constant field to an
algebraic closure.  The arithmetic monodromy contains the geometric group,
so it too is

\[
 \boxed{G_{\mathrm{arith}}=S_3.}                      \tag{6.2}
\]

In the Galois closure, the wild inertia responsible for (5.2) acts as a
transposition.  In characteristic two this order-two inertia is wild,
not the tame transposition occurring in the characteristic-zero weighted
pencil.

## 7. Comparison with the established marked-root mechanisms

The answer has two levels.

First, the map belongs to the broad marked-root mechanism:

1. its inverse is selection of a root of the finite projective cubic
   (2.2);
2. affine source space is exactly the regular-reconstruction open
   \(\overline X\setminus E\);
3. the marked root \(t\) stays finite on \(E\), while \(x=a^{-1}\) has a
   pole; and
4. the root derivative \(a\) cancels the reciprocal chart Jacobian
   \(a^{-1}\).

Second, it is **not** a reduction of either established branch.

| feature | weighted | cancellation | Huq--Kuruvilla |
|---|---|---|---|
| critical normalization after the free coordinate | \(\mathbb A^1\) | \(\mathbb G_m\) | \(\mathbb A^1\) |
| determinant ledger | distributed | reciprocal | reciprocal |
| generic critical-to-discriminant map | separable birational | separable birational | radicial degree \(2\) |
| generic boundary data | tame \(e=2,f=1\) in the cover | prescribed tame data in characteristic zero | \(e=1,f_{\rm insep}=2\), different \(1\) |

The first two columns require the normalized controlled divisor to map
finitely and birationally to the reduced discriminant.  Equation (5.2)
violates that condition in an invariant way.  Nor is the displayed map a
naive good reduction of a characteristic-zero Keller map: over
\(\mathbb Z\), the same formulas have Jacobian

\[
 1+2xy+2x^2z+4x^3yz+2x^4z^2+2x^5yz^2,
\]

whose nonconstant terms disappear only in characteristic two.

Therefore the structural verdict is

\[
\boxed{\text{marked-root algebraization, but a new wild--radicial
boundary type.}}                                      \tag{7.1}
\]

The existing puncture-rank dichotomy \(\mathbb A^1/\mathbb G_m\) is not
false, but it is incomplete in positive characteristic: it must be refined
by the separable/inseparable residue degree and different of the selected
boundary map.  The minimal new entry is an \(\mathbb A^1\) critical
normalization whose discriminant parameter is its Frobenius square.

## 8. The preserved coordinate and the dimension-two theorem

### 8.1 Mondello's theorem (external)

Mondello's Theorem 1.2 proves, over \(k=\overline{\mathbb F}_2\), that
the map (0.1) is etale with determinant one, has algebraically independent
target coordinates and exact separable function-field degree three, has
three points in its geometric generic fibre, and sends the three distinct
points \((0,1),(1,0),(1,1)\) to \((0,1)\).  Section 10 of the paper gives
the preserved-coordinate skew-product derivation from Huq--Kuruvilla's
threefold map.  The calculations below are an internal exact replay of that
external theorem and make the base-field dependence visible.

### 8.2 Preserved-coordinate derivation and hidden cubic

Permute the source and target coordinates of (1.1) by

\[
 \sigma(x,y,z)=(x,z,y),\qquad
 \tau(X,Y,Z)=(X,Z,Y).
\]

Then

\[
 \Phi=\tau\circ F\circ\sigma
 =\bigl(x+x^2z,\ y+x^2y^2,\ z+xy+x^2yz\bigr).       \tag{8.1}
\]

Introduce polynomial source coordinates

\[
 a=x,\qquad b=z+y^2,\qquad c=x+y+x^2b.              \tag{8.2}
\]

They have polynomial inverse

\[
 x=a,\qquad y=c+a+a^2b,\qquad
 z=b+(c+a+a^2b)^2.                                  \tag{8.3}
\]

If the target coordinates of (8.1) are denoted by `(X,Y,Z)`, make the
linear target change

\[
 A=X,\qquad B=Z,\qquad C=X+Y.
\]

Substitution gives `C=c`.  Hence the re-coordinatized threefold map is a
skew product

\[
 (a,b,c)\longmapsto (A(a,b,c),B(a,b,c),c).           \tag{8.4}
\]

The fiber `c=0` is preserved, and its first two coordinates are

\[
 \begin{aligned}
 A(a,b,0)&=a+a^2b+a^4+a^6b^2,\\
 B(a,b,0)&=b+a^5+a^6b+a^7b^2+a^8b^3.
 \end{aligned}                                      \tag{8.5}
\]

After renaming `(a,b)` as `(x,y)`, this is exactly (0.1).  Because the third
coordinate in (8.4) is unchanged, the threefold determinant is the plane
determinant on every fiber.

> **Preserved-coordinate reduction rule.** For every new counterexample,
> search for a polynomial source coordinate \(c(x)\) whose image is a target
> coordinate or a linear combination \(L(F)\) of target coordinates.  Run
> this test alongside grading detection, locally nilpotent derivation
> searches, tangent-quotient calculations, and finite-normalization analysis.

Here \(L=X+Y\) and \(L(\Phi)=c\), so the rule exposes a genuine
two-dimensional invariant fibre rather than merely a convenient
specialization.  Directly,

\[
 \frac{\partial P_2}{\partial x}=1,\quad
 \frac{\partial P_2}{\partial y}=x^2,\quad
 \frac{\partial Q_2}{\partial x}=x^4+x^6y^2,\quad
 \frac{\partial Q_2}{\partial y}=1+x^6+x^8y^2,
\]

so

\[
 \det\frac{\partial(P_2,Q_2)}{\partial(x,y)}=1.      \tag{8.6}
\]

Moreover

\[
 (0,1),\ (1,0),\ (1,1)\longmapsto(0,1),             \tag{8.7}
\]

which proves noninjectivity over every characteristic-two field.

It remains to check the prime-to-characteristic condition rather than infer
it from the ambient threefold.  Put

\[
 r=1+xy,\qquad u=1+x^3r,\qquad w=ru^2.
\]

In characteristic two,

\[
 P_2=xru,\qquad Q_2=y+x^5r^3,
\]

and the exact identities

\[
 xQ_2=1+w,\qquad P_2^2=x^2rw,
 \qquad P_2Q_2+P_2^3=ru+w^2                       \tag{8.8}
\]

hold.  Therefore `w` is a root of

\[
 H(T)=T^3+T^2+(P_2Q_2+P_2^3)T+P_2^3               \tag{8.9}
\]

and the reconstruction formulas

\[
 x=\frac{w+1}{Q_2},\qquad
 r=\frac{P_2^2}{x^2w},\qquad
 y=\frac{r+1}{x}                                    \tag{8.10}
\]

show that `k(x,y)=k(P_2,Q_2)(w)`.  In particular the target field has
transcendence degree two, so `P_2,Q_2` are algebraically independent.

For completeness, set `K=k(P_2)` and regard `Q_2` as transcendental over
`K`.  The polynomial (8.9) has degree one in `Q_2`.  If it factored in
`K(Q_2)[T]`, Gauss's lemma over `K[Q_2]` would make one positive-degree
monic factor independent of `Q_2`.  That factor would divide both

\[
 P_2T\quad\hbox{and}\quad T^3+T^2+P_2^3T+P_2^3.
\]

Since `P_2` is a unit in `K`, its only possible nonconstant divisor from the
first polynomial is `T`, but the second has nonzero constant term `P_2^3`.
This contradiction proves that `H` is irreducible.  Hence

\[
 [k(x,y):k(P_2,Q_2)]=3.                              \tag{8.11}
\]

Finally,

\[
 H'(w)=w^2+P_2Q_2+P_2^3=ru\ne0,                    \tag{8.12}
\]

so the extension is separable.  Equations (8.6), (8.7), and (8.11)--(8.12)
also prove the following base-field extension of Mondello's theorem.

> **Proposition 8.1 (`HKM2-ALLFIELDS`, repository corollary).** For every
> field \(k\) of characteristic two, (0.1) is an etale, noninvertible
> polynomial endomorphism of \(\mathbb A^2_k\) with separable generic degree
> three.  Its geometric generic fibre has three points.  After any extension
> of the constant field, including passage from a nonperfect field to an
> algebraic closure, the same conclusions hold.

**Proof.** All coefficients lie in \(\mathbb F_2\), and every calculation
above uses only the characteristic-two identities and the fact that
\(k[x,y]\) is a domain.  The recovery formulas make \(k(x,y)\) algebraic
over \(k(P_2,Q_2)\), hence prove algebraic independence of \(P_2,Q_2\).
The degree-one-in-\(Q_2\) Gauss-lemma argument is valid over
\(k(P_2)[Q_2]\) for every field \(k\), so it proves irreducibility over the
actual target field without any perfectness hypothesis.  The nonzero value
\(H'(w)=ru\) proves separability directly.  Repeating the same
actual-target argument after an arbitrary constant-field extension keeps
the degree equal to three; over an algebraic closure of the target function
field, separability gives exactly three embeddings and therefore three
geometric-generic-fibre points.  Finally, determinant one gives etaleness,
and the three rational collision points remain distinct over every such
field. \(\square\)

The new content of Proposition 8.1 is the base-field scope.  Its formulas
and proof mechanism are Mondello's; the repository does not claim a second
discovery of the plane theorem.

This does not address the characteristic-zero plane Jacobian conjecture.  For
the same formulas over the integers, the Jacobian is

\[
\begin{aligned}
1&+2xy+4x^3-4x^6-2x^7y+4x^9-2x^9y^3-2x^{10}y\\
 &+6x^5y^2+6x^{11}y^2-2x^{12}y^3+2x^{13}y^4,
\end{aligned}                                       \tag{8.13}
\]

whose nonconstant terms vanish only after reduction modulo two.

The
[plane normalization and wild-boundary continuation](HUQ_KURUVILLA_PLANE_BOUNDARY_NORMALIZATION.md)
computes the integral closure of the hidden cubic order, proves that the
source plane is its complement of one radicial boundary line, gives the
complete generic ledger `(e,f_sep,f_insep,different,sheet loss)=(1,1,2,1,2)`,
proves that the three retained/missing intersections are ordinary completed
nodes, and proves geometric and arithmetic monodromy `S_3`.

<!-- status-consumer: HKM2B1 c197aa4165670dc5 -->

The
[modulo-four obstruction theorem](HUQ_KURUVILLA_PLANE_W2_OBSTRUCTION.md)
strengthens the final lift check below.  It proves that arbitrary polynomial
corrections divisible by two cannot repair the Jacobian.  More precisely, the
full cokernel of the first variation is
`H^2_dR(F_2[x,y]) = xy F_2[x^2,y^2]`, and the integral error has class
`xy(1+x^3(1+xy))^2`.  Its Cartier image is `1+x^3(1+xy)`, so the class is
nonzero and remains nonzero under every polynomial plane left--right
equivalence.  Thus no equivalent representative has a Keller lift through
`W_2(F_2)`, without any correction-degree bound.  Every first obstruction
becomes exact after one identity stabilization; for this map the theorem
gives an explicit compatible determinant-one polynomial lift through every
finite Witt level.  Its inverse limit is restricted two-adic, not a
polynomial map over `Z_2`.  At the first stable level the lift has sharp
degree 18: the two-coefficient functional
`[x^13 y^4]+[x^14 y^5]` excludes every degree-at-most-17 lift, with arbitrary
corrections in all three variables.
At the next level an explicit plane correction gives a degree-25 lift over
`W_3(F_2)`, improving the elementary geometric-series degree 35.  A complete
1,083-variable, 1,639-constraint Boolean coefficient system excludes every
degree-18 lift.  A second exact system produces a degree-19 lift and a direct
polynomial replay gives an odd constant Jacobian modulo 8, normalizable to
determinant one, so `d_3=19`.  These Z3
certificates have no independent implementation or external human review.
Degree 25 is sharp only among extensions whose reduction modulo four is the
canonical degree-18 lift.
The preferred frozen degree-19 witness already has determinant one and 440
correction coefficients.  Fixing its first digit leaves an exact affine
second-completion space of dimension 804 (1,485 variables and rank 681).
The independent `z0`, `z1`, and `z2` blocks split further into incidence
components and have exact support minima 92, 54, and 14.  Thus the fixed-first
minimum is 160, also under the determinant-one constraint.  Together with
first support 280 this gives total 440; no global support minimum over all
first digits is claimed.

At `W_4`, the preferred degree-19 representative is a dead end: its next
determinant digit has a nonzero 48-term class in
`H^3_dR(F_2[x,y,z])`, so no correction of any degree extends that fixed map.
The block-triangular degree-25 representative has vanishing top class.  Its
35-term next error `L` has exact fixed-representative extension degree 52:
the functional `[x^39 y^12]+[x^40 y^13]` excludes degrees at most 51, and
the correction `8zL` attains degree 52.  Hence the unrestricted result is
only `19 <= d_4 <= 52`; exact `d_4` remains open.
The hoped-for cohomological UNSAT strengthening is false within the complete
degree-19 existence ansatz: 241 next-class equations have a SAT model with
support `327+491=818` and zero `H^3_dR` class.  Its full error has support
1,027 and degree 52; exact affine elimination proves that this fixed map also
has extension degree 52, attained by a `314+98+674` correction.  This does
not decide whether another degree-19 representative extends in degree 19.
An exact joint search now eliminates the full affine correction system for
each nonlinear master model and feeds dual certificates back as cuts.  All
sampled duals through 64 cuts are singleton codomain holes.  A master seeded
with all 599 such holes active on the pinned zero-class model returns
`unknown` at 600 seconds, not UNSAT.  A shared-minor quotient compiler also
reaches all 4,340 structural holes; the full master remains `unknown`, while
the `z1` and `z3` layers are individually SAT and `z0`,`z2` remain undecided.
The unrestricted problem remains open.

<!-- status-consumer: HKM2W1 904c57385ac0b0dd -->

<!-- status-consumer: HKM2W2 474e0d677133ee23 -->

<!-- status-consumer: HKM2W3 55d99efeee0298af -->

<!-- status-consumer: HKM2W4 6075dd4fbb9cb89c -->

## 9. Follow-up experiments

Items 1, 4, and 5, together with the first-Witt-level part of item 6, are
closed.  The others are research tasks, not claims supported by the current
certificates.

1. **Arbitrary base fields (closed by Proposition 8.1).** Keep the
   algebraic-independence, actual-target irreducibility, nonperfect-base, and
   geometric-fibre arguments separate and explicit.
2. **Independent normalization audit.** Recompute the integral closure of
   the primitive order, the singular locus, conductor, reconstruction open,
   two components above \(Q=0\), completed nodes, different, and residue
   degrees in a second CAS implementation that does not reuse the current
   elimination order or intermediate formulas.
3. **Independent \(S_3\) monodromy routes.** Add a discriminant/resolvent
   calculation and a target-line specialization with explicit Galois group;
   compare arithmetic and geometric monodromy and identify the
   characteristic-two wild transposition inertia precisely.
4. **Cohomological first lifting obstruction (closed by `HKM2W1`).** For
   \(\mathcal D_F(A,B)=A_xQ_y+P_xB_y+A_yQ_x+P_yB_x\), compute the full
   cokernel and the class of the integral Jacobian error \(K\).  The answer is
   `H^2_dR(F_2[x,y]) = xy F_2[x^2,y^2]`, with class `xyu^2` and Cartier image
   `u`; one identity coordinate kills every first obstruction.
5. **Left--right equivalence before lifting (closed by `HKM2W1`).** The
   obstruction class is functorial under source automorphisms and unchanged,
   up to an invertible pullback, by arbitrary polynomial source/target
   equivalence.  Jung--van der Kulk tameness supplies lifts of all plane
   automorphisms, so no bounded search is needed.
6. **Degree growth in the stabilized Witt tower.** The first two levels are closed:
   `HKM2W2` proves `d_2=18` and
   `18 <= d_n <= 17(n-1)+1`; its degree-25 construction also gives
   an upper bound at the next level.  `HKM2W3` excludes degree 18 without a
   sparse ansatz and supplies a directly replayed degree-19 witness, proving
   `d_3=19`.  `HKM2W4` shows that this preferred degree-19 representative
   cannot extend, proves exact extension degree 52 for a fixed degree-25
   representative, and gives `19 <= d_4 <= 52`.  Determine exact `d_4` by
   imposing the bounded degree-19 `W_4` correction equations jointly—the
   weaker next-class-zero system is SAT—then compute later
   levels, prove growth bounds, and decide algebraicity of the inverse limit
   over the restricted power-series ring.
7. **Odd-characteristic analogues.** Use prescribed hidden-root covers and
   boundary residue data in characteristics \(3\) and \(5\), initially in
   generic degrees \(2,3,4\), instead of a random coefficient search.
8. **Preserved-coordinate detector.** Search degree by degree for a
   polynomial coordinate \(c(x)\) and a target linear combination \(L\)
   with \(L(F)=c\), using linear algebra and Groebner reduction; apply it to
   the foundational, weighted, cancellation, quadratic-gauge, cubic BCW,
   and positive-characteristic examples.

## 10. Exact reproduction

The accompanying dependency-light SymPy audit checks the threefold map,
cubic, discriminant, projective collision, normalization-chart identities,
reconstruction, pole divisor, discriminant pullback, and determinant ledger.
It also checks the coordinate permutations, polynomial source change, skew
product, plane fiber, plane determinant and collision, hidden cubic and
recovery identities, separability witness, and failure of the naive integral
plane lift:

```bash
.venv/bin/python scripts/verify_huq_kuruvilla_characteristic_two.py
```
