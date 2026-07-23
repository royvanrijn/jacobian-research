# Post-coordinate attacks for the absolute Sunada pair

The plane-coordinate route is closed.  The
[node-separation audit](DAVENPORT_NODE_SEPARATING_AFFINE_MODIFICATION.md)
proves that no polynomial coordinate line can separate the Davenport
derivative node while preserving its full affine normalized boundary.

The surviving problem is higher-dimensional.  This note ranks the remaining
attacks and gives an explicit success or failure gate for each one.

Work over

\[
K=\mathbb Q(a),\qquad a^2+a+2=0.
\]

## 1. Certified starting point

There are two useful affine-plane charts.

The secant chart is

\[
(2Y+1)w=4T+1.
\]

It separates the node but omits one additional torus point.  The parabola
chart is

\[
(T+Y^2)z=2Y+1.
\]

It has affine-plane ambient, Jacobian \((T+Y^2)/2\), and

\[
\pi_{\rm par}^*J=(T+Y^2)^2S_{\rm par}
\]

with smooth strict transform.  It omits only the two branches over the
coordinate origin.

The two charts cover the full normalized derivative boundary.  Their
exceptional transition is

\[
z=\frac4{w-2},
\]

so their natural ambient gluing is the blowup

\[
B=\operatorname {Bl}_p\mathbb A^2
\]

with exceptional curve \(E\simeq\mathbb P^1\).  Thus

\[
[B]=\mathbb L^2+\mathbb L.                            \tag{1}
\]

Any surviving construction must remove the extra \(\mathbb L\) without
identifying the two distinguished points of \(E\).

## 2. Attack A: globalize the determinant parameter

The fixed-parameter determinant threefold from the weighted-gluing audit is

\[
D^2U-CR=\beta.
\]

For fixed \(\beta\ne0\), its reduced modification center is
\(\mathbb G_m\), producing the unit-rank mismatch with the normalized
Davenport derivative boundary.

Do not fix \(\beta\).  Instead use the polynomial map

\[
\Psi:\mathbb A^4_{D,C,R,U}\longrightarrow
\mathbb A^4_{D,C,R,\beta},
\qquad
\beta=D^2U-CR.                                       \tag{2}
\]

Both source and target are affine four-space, and

\[
\operatorname {Jac}(\Psi)=D^2.                       \tag{3}
\]

The reduced center over \(D=0\) is now

\[
D=0,\qquad \beta=-CR,
\]

which is \(\mathbb A^2_{C,R}\), not \(\mathbb G_m\).  The stable unit-rank
obstruction of AS6 therefore disappears at the total-family level.

### The direct splice is not Keller

There is an immediate but unsuccessful way to combine (2) with the
parabola chart.  Put

\[
Y=\frac{Dz-1}{2},\qquad T=D-Y^2
\]

and form

\[
\begin{aligned}
F:\mathbb A^5_{D,z,C,R,U}&\longrightarrow
  \mathbb A^5_{T,Y,C,R,\beta},\\
(D,z,C,R,U)&\longmapsto
\left(D-Y^2,Y,C,R,D^2U-CR\right).                   \tag{4}
\end{aligned}
\]

The two relative blocks have determinants

\[
\det\frac{\partial(T,Y)}{\partial(D,z)}=\frac D2,
\qquad
\det\frac{\partial(D,C,R,\beta)}
         {\partial(D,C,R,U)}=D^2.
\]

Because they share only the preserved coordinate \(D\), the full
determinant is

\[
\boxed{\operatorname {Jac}(F)=\frac{D^3}{2}.}        \tag{5}
\]

Thus the direct fiber-product splice is not Keller.  Polynomial
automorphisms on either side have constant nonzero Jacobian, and adjoining
identity variables does not change (5), so neither recoordination nor
stabilization repairs it.

This is a useful new restriction: promoting \(\beta\) removes the center
unit-rank obstruction, but the two modifications must be **coupled**.
Their raw Jacobian factors cannot merely be multiplied.

### The exceptional equation cannot be exposed

There is also a generic-degree obstruction which precedes the Jacobian
calculation.  Put

\[
k=K(T,u),\qquad L_g=K(T,Y),\qquad u=g_T(Y).
\]

The polynomial \(g_T(Y)-u\) is irreducible over \(k\), so

\[
[L_g:k]=7.                                           \tag{6}
\]

If a stabilized target retaining \(k\) exposed

\[
D=\delta=T+Y^2
\]

as a target function, then \(Y\) would satisfy

\[
Y^2=D-T.
\]

The modified generic degree would therefore be at most two, contradicting
(6).  Equivalently, adjoining purely transcendental target variables does
not make the algebraic element \(T+Y^2\) descend from \(L_g\) to \(k\).

Thus the phrase “identify \(D\) with \(T+Y^2\)” is not compatible with
retaining the old Davenport base as a subfield.  Either the exceptional
equation must be masked by source variables, or the common target field
must be replaced and its degree-seven Gassmann closure proved again.

The simplest degree-preserving mask is

\[
d=s+\delta.
\]

Indeed

\[
L_g(s)=L_g(d),\qquad [L_g(d):k(d)]=7.                \tag{7}
\]

But the map

\[
(T,Y,s)\longmapsto(T,g_T(Y),s+\delta)
\]

still has Jacobian \(J=g_T'(Y)\).  Masking solves the field problem, not
the Keller problem.

### Affine-linear masks always trivialize when \(T\) is fixed

More generally, let \(f,h,a,b\in K[T,Y]\) and consider

\[
\Theta(T,Y,s)=
\bigl(T,\ f(T,Y)+a(T,Y)s,\ h(T,Y)+b(T,Y)s\bigr).     \tag{8}
\]

Its Jacobian is

\[
\operatorname {Jac}(\Theta)
=bf_Y-ah_Y+s(ba_Y-ab_Y).                             \tag{9}
\]

Suppose (9) is a nonzero constant.  Its \(s\)-coefficient gives

\[
ba_Y-ab_Y=0.
\]

Over \(K(T)\), this says that \(a/b\) is independent of \(Y\).  Removing
their common \(Y\)-factor in the constant term of (9) shows that in fact

\[
a=A(T),\qquad b=B(T),
\]

up to one common nonzero scalar.  Moreover \(A\) and \(B\) are coprime in
\(K[T]\), and

\[
\frac{\partial}{\partial Y}(Bf-Ah)=c\in K^*.
\]

Hence

\[
Bf-Ah=cY+q(T).                                      \tag{10}
\]

The same mask-cancelling combination of target coordinates recovers
\(Y\):

\[
B(f+As)-A(h+Bs)=cY+q(T).
\]

Bézout for \(A,B\) then recovers \(s\) polynomially.  Consequently every
constant-Jacobian map of the form (8) is a polynomial automorphism.

The same argument works for arbitrarily many affine-linear masks.  Let
\(\mathbf s=(s_1,\ldots,s_m)\), let
\[
A(T)\in\operatorname {Mat}_{m+1,m}(K[T]),
\]
and consider
\[
\Theta(T,Y,\mathbf s)
=\left(T,\ \mathbf h(T,Y)+A(T)\mathbf s\right).      \tag{11}
\]
The relative Jacobian is
\[
\det[\mathbf h_Y\mid A].
\]
Let \(\boldsymbol\lambda(T)\) be the signed vector of maximal minors of
\(A\).  Then
\[
\boldsymbol\lambda^tA=0,\qquad
\det[\mathbf h_Y\mid A]
=\boldsymbol\lambda^t\mathbf h_Y
\]
up to the fixed orientation sign.  If this determinant is
\(c\in K^*\), then
\[
\boldsymbol\lambda^t\mathbf h=cY+q(T),
\qquad
\boldsymbol\lambda^t
  (\mathbf h+A\mathbf s)=cY+q(T).                   \tag{12}
\]
Thus the targets recover \(Y\).  The maximal minors generate the unit ideal
because their polynomial combination equals \(c\); consequently \(A\) is
a split injection over \(K[T]\), and its polynomial left inverse recovers
\(\mathbf s\).  The map (11) is a polynomial automorphism.

This closes every fixed-\(T\), affine-linear multi-mask attack, not only
the one-auxiliary case.  A degree-seven construction must use nonlinear
auxiliary incidence or a target replacement which does not retain \(T\).

### The two natural nonlinear placements also fail

First suppose the shifted exceptional mask itself is a target coordinate:

\[
\Theta(T,Y,s)=
\left(T,\ P(T,Y,s),\ d=s+T+Y^2\right).              \tag{13}
\]

The Keller equation is

\[
P_Y-2Y P_s=c\in K^*.                                \tag{14}
\]

The polynomial invariants of
\(\partial_Y-2Y\partial_s\) are \(T\) and \(s+Y^2\).
Consequently every polynomial solution of (14) is

\[
P=cY+F(T,s+Y^2).                                    \tag{15}
\]

But \(s+Y^2=d-T\), so (15) recovers \(Y\) polynomially from the targets,
and then

\[
s=d-T-Y^2.
\]

Thus every Keller map of the form (13) is a polynomial automorphism.  The
shifted mask (7) cannot itself be one of the final target coordinates.

There is a complementary obstruction if the mask is used to perturb the
Davenport value.  Let

\[
P=g_T(Y)+(T+Y^2)H(T,Y,s)
\]

and suppose that \((T,P,Q)\) had constant Jacobian for some polynomial
\(Q\).  Reducing its relative Jacobian modulo
\(\delta=T+Y^2\) gives

\[
\left(
J(-Y^2,Y)+2Y\,\overline H
\right)\overline{Q_s}=c.                            \tag{16}
\]

The exact parabola restriction is

\[
J(-Y^2,Y)=(a-2)Y^4(2Y+1)^2.
\]

The first factor in (16) is therefore divisible by \(Y\), so it cannot
divide a nonzero constant in \(K[Y,s]\).  No such \(Q\) exists.

Equations (13)--(16) close both canonical one-variable nonlinear uses of
the exceptional equation.  The next ansatz must mix the \(T\)-coordinate
with the determinant variables \(C,R,U\); merely placing \(\delta\) in a
shifted output or multiplying a correction by \(\delta\) cannot work.

### The unshifted determinant parameter cannot be a Keller output

There is a final elementary obstruction to Attack A as originally stated.
If a polynomial affine-space map has

\[
\beta=D^2U-CR                                       \tag{17}
\]

as one of its target coordinates, then

\[
d\beta=(2DU)\,dD-R\,dC-C\,dR+D^2\,dU.
\]

On the nonempty linear subspace

\[
D=C=R=0
\]

this entire row of the Jacobian matrix vanishes.  No choice of the other
target coordinates can make the full Jacobian nonzero there.  Thus:

\[
\boxed{\text{No affine-space Keller map can retain (17) literally as a
target coordinate.}}                                \tag{18}
\]

The total \(\beta\)-family removes the fixed-fiber unit-rank mismatch, but
it reintroduces a critical axis.  It is a useful residue model, not an
absolute Keller map.

The minimal repair is a translated determinant incidence

\[
\widehat\beta=S+D^2U-CR.                             \tag{19}
\]

Now \(\partial\widehat\beta/\partial S=1\), so its differential is
unimodular.  However, if \(S\) is also separately recoverable from the
other target coordinates, (19) is merely a triangular suspension and
does not cancel the Davenport derivative.  The variable \(S\) must be
used to replace, rather than accompany, one old target coordinate.

### Exact gate for the translated determinant attack

Couple (19) to the parabola chart using a masked version of the exceptional
equation

\[
\delta=T+Y^2.
\]

The required polynomial diagram must satisfy:

1. the two Davenport covers map to one common affine target containing
   \(\widehat\beta\), not the singular output \(\beta\);
2. on the slice \(S=0\), the determinant incidence recovers the \(D^2\)
   residue contribution matching the doubled derivative column;
3. the construction is not equivalent to the failed product (4), and after
   including the parabola-chart Jacobian \(D/2\), every remaining
   \(D\)-valuation cancels in the residue form;
4. the common open retains the point/line
   \(\mathrm {GL}_3(\mathbb F_2)\) Gassmann closure;
5. \(S\) is not separately recoverable from the remaining target
   coordinates;
6. neither \(T+Y^2\) nor its line-cover analogue is exposed in the target;
7. the auxiliary coupling lies outside the affine-linear multi-mask class
   (11);
8. it is not one of the nonlinear placements (13) or (16), so at least one
   target coordinate mixes \(T\) with the determinant variables.

Failure of item 5 makes (19) a triangular suspension and leaves the
Davenport derivative uncancelled.  Specializing \(S\) instead returns to
the non-affine fixed-\(\beta\) threefold.  This remains the highest-priority
attack because it keeps the exact \(J^2\) residue model while removing both
the fixed-center unit mismatch and the critical axis of (17).

### Normal form of the translated incidence

The translation (19) is itself a triangular source coordinate:

\[
\Sigma=S+D^2U-CR,\qquad
(D,C,R,U,S)\longleftrightarrow(D,C,R,U,\Sigma),      \tag{20}
\]

and this change has Jacobian one.  This gives two further reductions.

First suppose a proposed map retains \(T,D,C,R,\Sigma\) among its target
coordinates.  After (20), its remaining two outputs are polynomials
\[
P(T,Y,U,\Sigma,D,C,R),\qquad
Q(T,Y,U,\Sigma,D,C,R).
\]
The full Jacobian, up to a harmless coordinate-ordering sign, is exactly

\[
\det\frac{\partial(P,Q)}{\partial(Y,U)}
\quad\text{with }T,D,C,R,\Sigma\text{ fixed}.        \tag{21}
\]

Thus the translated determinant block supplies no cancellation: the
candidate is simply a family of plane Keller maps in \((Y,U)\).  A
nontrivial degree-seven solution of this form would already be a
two-variable Jacobian-conjecture counterexample over
\(K(T,D,C,R,\Sigma)\), independently of the determinant geometry.

There is a similar reduction after replacing \(T\) if the remaining
dependence on \(U\) is affine-linear with parameter-only coefficients.
Over any coefficient ring, a unimodular vector
\((a_0,a_1,a_2)\) can be completed to an invertible matrix.  Hence a map

\[
(T,Y,U)\longmapsto
\bigl(H_0(T,Y)+a_0U,\,
      H_1(T,Y)+a_1U,\,
      H_2(T,Y)+a_2U\bigr)
\]

is polynomially equivalent to

\[
(T,Y,U)\longmapsto
\bigl(P_0(T,Y),P_1(T,Y),U+P_2(T,Y)\bigr),            \tag{22}
\]

and its Jacobian and generic degree are those of the plane pair
\((P_0,P_1)\).

Consequently the translated incidence becomes genuinely new only if:

1. \(T\) is not retained as a target coordinate;
2. the post-translation \(U\)-dependence is nonlinear or has essential
   \(T,Y\)-dependent coefficients; and
3. the resulting three-variable block is not stably equivalent to a plane
   Keller pair.

### Every \(Y\)-dependent coefficient curve fails

The first case outside (22) lets the \(U\)-coefficient depend on the
Davenport sheet but not on \(T\):

\[
\Phi(T,Y,U)=
\bigl(T+a(Y)U,\ g_T(Y)+b(Y)U,\ h(T,Y)+U\bigr).       \tag{23}
\]

Its \(U^2\)-Jacobian coefficient vanishes automatically.  Put

\[
W=ba'-ab'.
\]

The coefficient of \(U\) is

\[
Wh_T-a'\partial_Tg+b'.                              \tag{24}
\]

If \(W=0\), the nontrivial \(T\)-dependence of \(\partial_Tg\) forces
\(a'=b'=0\), returning to the parameter-linear plane reduction (22).
Suppose \(W\ne0\).  For (24) to vanish, \(W\) must divide
\(a'\partial_Tg-b'\) in \(K[T,Y]\).  The positive-\(T\) coefficients of
\(\partial_Tg\) are coprime in \(K[Y]\), so

\[
W\mid a',\qquad W\mid b'.
\]

Write

\[
a'=Wr,\qquad b'=Ws.
\]

The definition of \(W\) then gives the Bézout identity

\[
br-as=1,                                             \tag{25}
\]

and (24) forces

\[
h(T,Y)=r(Y)g_T(Y)-s(Y)T+k(Y).                       \tag{26}
\]

Exact substitution into the constant term of the Jacobian yields the
complete factorization

\[
\operatorname {Jac}(\Phi)
=\bigl(a\,\partial_Tg-b\bigr)
 \bigl(r'g_T-s'T+k'\bigr).                           \tag{27}
\]

Here the right side is a product, not a sum.  The first factor depends
nontrivially on \(T\), because \(a\ne0\) and \(\partial_Tg\) has degree two
in \(T\).  It is not a unit of \(K[T,Y]\), so (27) cannot be a nonzero
constant.

Thus every projective \(U\)-coefficient curve depending only on \(Y\) is
eliminated, in every polynomial degree.  The next coefficient search must
have genuine simultaneous \(T,Y\)-dependence.

The complementary \(T\)-only case also closes.  For

\[
\Phi_T(T,Y,U)=
\bigl(T+a(T)U,\ g_T(Y)+b(T)U,\ h(T,Y)+U\bigr),       \tag{28}
\]

put \(W=ab'-ba'\), now differentiating in \(T\).  Its \(U\)-coefficient is

\[
Wh_Y+a'J.                                           \tag{29}
\]

If \(W\ne0\), monicity of \(J\) in \(Y\) forces \(W\mid a'\).  Write
\(a'=Wr\).  Equation (29) gives

\[
h=-r(T)g_T(Y)+k(T).
\]

The constant term of the Jacobian becomes

\[
J\bigl(ar'g-ak'+br+1\bigr),
\]

which is divisible by the nonunit \(J\) and cannot be a nonzero constant.
If \(W=0\), the same equation either returns to constant parameter
coefficients or to the fixed-\(T\) plane reduction.

Hence dependence through either base coordinate separately is impossible.
Only genuinely simultaneous \(T,Y\)-dependence, or nonlinear dependence on
\(U\), remains.

The geometrically distinguished mixed coordinate also fails.  Put

\[
x=T+Y^2,\qquad y=Y,\qquad
G(x,y)=g_{x-y^2}(y).
\]

For coefficient curves depending only on \(x\), consider

\[
\Phi_x(x,y,U)=
\bigl(x-y^2+a(x)U,\ G(x,y)+b(x)U,\ h(x,y)+U\bigr).
                                                               \tag{30}
\]

With \(W=ab'-ba'\), the \(U\)-coefficient is

\[
Wh_y+a'G_y+2yb'.                                    \tag{31}
\]

The coefficient of \(y^6\) in \(G_y\) is a nonzero constant.  Therefore
(31) forces \(W\mid a'\), and then \(W\mid b'\).  Write

\[
a'=Wr,\qquad b'=Ws,\qquad as-br=1.
\]

Equation (31) gives

\[
h=-r(x)G(x,y)-s(x)y^2+k(x).
\]

The constant term of the Jacobian factors exactly as

\[
\operatorname {Jac}(\Phi_x)
=\bigl(r'G+s'y^2-k'+s\bigr)
 \bigl(aG_y+2yb\bigr).                               \tag{32}
\]

Again the right side is a product.  Its second factor has the same nonzero
\(y^6\)-coefficient multiplied by \(a\ne0\), so it is not a unit.
Consequently (32) cannot be a nonzero constant.

Thus even the parabola coordinate \(T+Y^2\) singled out by the unique clean
node chart cannot support the projective \(U\)-coefficient curve.  A live
affine-in-\(U\) ansatz must use a different genuinely mixed pencil, not
\(T\), \(Y\), or \(T+Y^2\).

There is a coordinate-invariant compiler for the remaining pencils.  Let
\((x,y)\) be polynomial coordinates on the Davenport source plane, write

\[
T=T(x,y),\qquad G(x,y)=g_{T(x,y)}(Y(x,y)),
\]

and consider

\[
\Phi_{x}(x,y,U)=
\bigl(T+a(x)U,\ G+b(x)U,\ h(x,y)+U\bigr).            \tag{33}
\]

With \(W=ab'-ba'\), the \(U\)-coefficient is

\[
Wh_y+a'G_y-b'T_y.
\]

Whenever the coefficient content forces
\(a'=Wr,\ b'=Ws\), one has \(as-br=1\) and

\[
h=-rG+sT+k(x).
\]

The constant Jacobian then factors invariantly as

\[
\boxed{
\operatorname {Jac}(\Phi_x)
=\bigl(aG_y-bT_y\bigr)
 \bigl(r'G-s'T-k'\bigr).
}                                                     \tag{34}
\]

Again the displayed right side is a product.  This gives the fast gate

\[
\boxed{aG_y-bT_y\in K^*.}                            \tag{35}
\]

Every affine-linear pencil fails this gate.  If \(Y_y\ne0\), then along an
affine \(y\)-line the term \(Y^7/7\) makes \(G_y\) have degree six in \(y\),
whereas \(T_y\) is constant.  If \(Y_y=0\), then \(T_y\ne0\) and
\(\partial_Tg\) has nonzero quadratic dependence on \(T\), so \(G_y\) is
again nonconstant.  Since \(W\ne0\) implies \(a\ne0\), no cancellation can
make (35) a unit.

Thus all affine-linear mixed pencils are eliminated at once.  The next
search begins with nonlinear polynomial coordinates beyond the already
excluded parabola \(T+Y^2\), and applies the unit gate (35) before solving
for \(h\).

## 3. Attack B: coupled affine replacement of the exceptional line

A plain stabilization of \(B\) cannot help: \(E\times\mathbb A^m\) still
contains the complete curve \(E\).

A plain affine-line or Jouanolou torsor is also insufficient.  Any
Zariski-locally trivial \(\mathbb A^1\)-bundle \(Q\to B\) has

\[
[Q]=\mathbb L[B]
=\mathbb L^3+\mathbb L^2\ne\mathbb L^3.              \tag{36}
\]

Thus even when \(Q\) is affine, it is not affine three-space.

The torsor must be coupled to a second affine modification which removes
the extra \(\mathbb L^2\).  The determinant family (2) supplies a natural
candidate: replace the complete fiber \(E\) by a determinant incidence
whose fixed fibers are \(SL_2\)-type but whose total \(\beta\)-family is
affine space.

### Exact gate

Construct an affine modification

\[
\widetilde B\longrightarrow B\times\mathbb A^m
\]

such that:

1. the inverse image of \(E\) is affine and still contains two disjoint
   labelled sections representing the node branches;
2. \([\widetilde B]=\mathbb L^{m+2}\);
3. its residue Jacobian contributes exactly the missing power of
   \(\delta\);
4. contracting either labelled section does not identify it with the other;
5. the construction is compatible with the transition
   \(z=4/(w-2)\).

Equation (36) is a quick rejection test: any proposal which is merely an
affine bundle over \(B\) fails before Jacobians are considered.

## 4. Attack C: descend a finite correspondence, not an ambient involution

The AS7 involution is an automorphism of the normalized three-puncture
boundary and has quotient \(\mathbb G_m\), but it does not lift over either
node contraction.  Requiring an ambient involution may be unnecessarily
strong.

Instead, retain the finite degree-two algebra of the normalized boundary
over its quotient and extend that algebra to an affine neighborhood.  The
ambient object would be a finite correspondence rather than a group action.

### Exact gate

Let \(C^\nu\) be the normalized derivative boundary and
\(q:C^\nu\to\mathbb G_m\) its canonical quotient.  Construct a finite flat
rank-two algebra \(\mathcal A\) on an affine neighborhood \(V\) of the
proposed determinant center such that

\[
\mathcal A|_{\mathbb G_m}=q_*\mathcal O_{C^\nu}.
\]

It must additionally satisfy:

1. trivial relative dualizing determinant, so it introduces no new
   Jacobian character;
2. compatibility with both Davenport point/line covers on the common open;
3. no finite re-gluing of the old node conductor;
4. affine-space total spectrum after the determinant parameter is included.

This route avoids the proved impossibility of lifting the involution over a
single exceptional fiber.

## 5. Attack D: change the non-symmetric monodromy core

If A--C fail, the obstruction is specific to the Davenport derivative
boundary rather than to Sunada theory.

The replacement core should be searched under three simultaneous filters:

1. a common Galois closure with an almost-conjugate nonconjugate subgroup
   pair;
2. a branch pullback whose doubled center is already affine or has a
   polynomial-coordinate normalization;
3. a Cox ledger admitting an affine-space realization without fixed
   determinant parameters.

The Cox-ledger multi-boundary construction remains the natural source.
Generic tangent pencils are excluded because their full monodromy is
symmetric.

## 6. Priority order

The concrete order is:

1. **A — translated determinant incidence.**  Use (19), with \(S\)
   replacing an old target coordinate, replace \(T\), and use essential
   nonlinear or sheet-dependent \(U\)-coupling.  Do not retain the singular
   output (17), expose \(T+Y^2\), take the failed product (4), return to the
   affine-linear mask classes (8) and (11), use the nonlinear placements
   (13) and (16), fall into the plane reductions (21)--(22), or return to
   any one-coordinate coefficient curve (23)--(29), or use the exceptional
   parabola pencil (30)--(32) or any affine-linear pencil (33)--(35).
2. **B — coupled exceptional-fiber modification.**  Use (1) and (36) as
   class filters before any coefficient search.
3. **C — finite algebra descent.**  Attempt this only if ambient
   equivariance remains the sole obstruction after A.
4. **D — new monodromy core.**  Preserve the completed Sunada group theory
   and change only the boundary geometry.

The project is no longer waiting for an unspecified affine normalization.
Its first unresolved object is a translated determinant splice based on
(19), with \(S\) and \(T\) both genuinely absorbed into a nonlinear
three-variable block, while avoiding the field obstruction (6), the
affine-mask theorems (8)--(12), the nonlinear obstructions (13)--(16), the
critical-axis obstruction (18), and the plane reductions (21)--(22).  The
first live affine-in-\(U\) coefficient map must also have simultaneous
\(T,Y\)-dependence beyond the three natural pencils \(T\), \(Y\), and
\(T+Y^2\), and beyond every affine-linear pencil, by (23)--(35).
