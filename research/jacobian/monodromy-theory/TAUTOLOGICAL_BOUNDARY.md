# Tautological sweep rank, local inertia and missing boundary branches

Work over algebraically closed characteristic-zero constants unless a
coefficient field is specified. Put

\[
 S(\gamma,w,z)=X(w,z)+\gamma\Delta,\qquad
 \Delta=(1,w,z),\quad \det JS=\gamma^2. \tag{1}
\]

These statements concern normalized Gao sweep constructions.
They do not classify every polynomial Keller action.

**KMFI38.** If X has rank at most one, the full geometric monodromy
of S is a natural symmetric group (including the trivial degree-one
case). For rank one, X factors through an affine linear polynomial
in w,z. After an affine source change and matching linear target
change, X=X(z), X3'=zX1'. Its exact field pair is a spectator
stabilization of the generic linear pencil
Y3=Y1*z+X3(z)-zX1(z). If deg X1=d>=1, the degree is d+1 and the
group S_(d+1); if X1 is constant, the map is birational.

**KMFI39.** Finiteness is unnecessary for the following local cycle
gate. For a generically finite normalized sweep S with
det JS=c*gamma^k, if gamma=0 maps onto a divisor with generic
degree delta, its valuation supplies delta cycles of length k+1
in the full inertia permutation. Missing normalization sheets may
supply additional cycles. In an actual polynomial Gao target twist,
the full inertia must still have a fixed sheet.

Consequently the natural degree-twelve PSL2(11) action still needs
k=4 or10 and affine architecture dimension at least6, respectively12,
for every such divisor-preserving sweep, even when nonfinite.
The exceptional degree-eleven action still needs k=1,2 or4.
The finite equalities for delta in KMFI27 become upper bounds
when boundary sheets are omitted; the missing ramified branches
must complete the same allowed passport.

In particular **no tautological four-dimensional sweep (1), of any
rank and with or without finiteness or a multiplier chart, can have
a polynomial Gao lift with the natural degree-twelve PSL2(11) action**.
Rank at most one is symmetric by KMFI38. At rank two, a present
three-cycle in inertia is incompatible with a fixed sheet in that
natural action. This is an obstruction to the declared construction,
not an absolute action-level exclusion under arbitrary birational charts.

For the exceptional degree-eleven action at rank two, the full
inertia must be 3^3*1^2. KMFI32 gives delta<=2, so at least one
ramified index-three normalization branch must be missing from S.

**KMFI40.** Suppose the nonzero normalization matrix of (1) is

\[
 N JX=\lambda
 \begin{pmatrix}-\rho&1\\-\rho^2&\rho\end{pmatrix},\quad
 N=\begin{pmatrix}-w&1&0\\-z&0&1\end{pmatrix},\quad
 \lambda\in k[w,z]\setminus0,\quad \rho\in k(w)[z].
 \tag{2}
\]

Any genuine pole in w of rho collapses the corresponding source
line under X, so S is nonfinite. In particular every finite sweep
in this chart has rho a nonconstant polynomial in w and falls
under KMFI28.

If rho depends on z, then X is birational onto its image,
for every allowed lambda, including positive z degree.
Thus any attempted exceptional degree-eleven PSL2(11) Gao lift
in this whole rational chart would require **exactly two additional
index-three ramification branches above the same branch divisor,
both outside the polynomial sweep source**. They cannot be inferred
from the Jacobian or from the source ramification surface.
This is a necessary boundary requirement, not a proof that the
remaining exceptional action is realizable or impossible.

[Exact checker](verify_tautological_boundary.py) ·
[model](tautological_boundary_model.json) ·
[receipt](tautological_boundary_receipt.json).

## 1. A rank-one polynomial surface has a linear parameter

Here is the needed polynomial factorization argument.
Normalize the affine image curve of X. Its coordinate ring embeds
in k[w,z]: an element integral over the image ring is integral over
the normal polynomial ring and lies in its fraction field.
A general affine line in the source maps nonconstantly to this curve.
Its smooth projective completion is therefore rational, and has
at most one omitted point: every omitted point must pull back to
the single infinity of that affine line. The normalized image is A1.
Thus X_i=f_i(T) for a polynomial T.

Take the relative algebraic closure of k(T) in k(w,z).
It is finite over k(T). The integral closure of k[T] in that field
again embeds in k[w,z], and the same affine-line argument makes it
k[U] for a polynomial U. Therefore X_i=g_i(U), with k(U)
relatively algebraically closed in k(w,z).
The generic fiber U=t is geometrically irreducible in characteristic
zero. This step avoids assuming that an arbitrary common parameter
already has irreducible generic fibers.

The gamma coefficient in the determinant of (1) is zero:

\[
 (g_2'(U)-w g_1'(U))U_w+
 (g_3'(U)-z g_1'(U))U_z=0. \tag{3}
\]

Pass to the algebraic closure of k(t) on the generic fiber U=t.
If g1'(t)!=0, (3) says that the Euler vector field centered at

(g2'(t)/g1'(t),g3'(t)/g1'(t))

is tangent to that fiber.
The irreducible polynomial U-t divides its Euler derivative.
The latter has degree at most d=deg U, so the quotient is the
constant d, by comparing highest homogeneous parts.
After shifting the center, U-t is a homogeneous binary form of
degree d. Over algebraically closed constants it is irreducible
only when d=1.

If g1'=0, a nonzero constant vector field is tangent instead.
Its derivative of U-t has degree at most d-1, so divisibility
forces the derivative to vanish identically.
Then U-t is a polynomial in one linear form; irreducibility
again gives d=1. Rank one guarantees that the constant vector
is not zero. Hence U is affine linear in both cases.

A matching affine change in w,z and linear change of the three
target coordinates preserves Delta=(1,w,z) and det JS=gamma^2.
We may therefore take U=z.

## 2. The exact rank-one group

Write X=(f(z),g(z),h(z)). The gamma coefficient of (1) gives
h'=z f'. Eliminating gamma and w gives

\[
 \gamma=Y_1-f(z),\quad
 w={Y_2-g(z)\over Y_1-f(z)},\quad
 Y_3=Y_1z+P(z),\quad P=h-zf,\quad P'=-f.
 \tag{4}
\]

If f is nonconstant of degree d, P has degree d+1 with nonzero
constant leading coefficient. Equation(4) is the exact source
and target field pair with independent spectator Y2.

A generic linear pencil P(z)+a*z has full symmetric group.
For clarity, its critical values are generically distinct:
their curve is parametrized by
a=-P'(z), b=P(z)-zP'(z), and db/da=z wherever P''!=0.
The derivative ratio belongs to the function field of that curve,
so this parametrization is birational.
The generic critical points are simple and their inertia consists
of single transpositions. These generate the transitive group
of the rational source, hence S_(d+1).

If f is constant, h is constant, and (4) recovers z rationally.
If X is constant from the start, S is also birational.
This proves KMFI38.

## 3. A local ramification cycle survives a nonfinite source

Now let S:A^m->A^m be generically finite with
det JS=c*gamma^k, gamma a source coordinate.
Assume D0={gamma=0} dominates a target divisor D.
At its generic point the source is quasi-finite over D.
Zariski's Main Theorem identifies this neighborhood with an open
part of the finite normalization; it does not add the omitted
normalization branches to the polynomial source.

In characteristic zero, the valuation is tame and its different
exponent is e-1. The Jacobian order is k, so its ramification
index is e=k+1. After passing to the geometric generic point of D,
the residue degree delta gives delta distinct e-cycles in full
inertia. Other normalization points over D can contribute further
cycles; this is exactly why the uniform finite formula cannot
simply be reused without finiteness.

For a Gao target twist, the field pair is the sweep with independent
C and the target substitution is Y_i=C^(e_i)*F_(i+1), C=F1.
On C!=0 this is an isomorphism with A^m_Y times G_m.
Thus D times G_m is still a divisor, with the same inertia, and
its closure is an affine target prime divisor.
If the lifted F is Keller, KMFI1 forces a fixed sheet in this
full inertia. Further independent stabilization and polynomial
coordinate equivalence preserve this argument.

The complete PSL2(11) fixed-sheet profiles give:

| action | present index e | k=e-1 | total ramified cycles | allowed present delta | missing ramified cycles |
|---|---:|---:|---:|---:|---:|
| natural degree12 |5|4|2|1 or2|2-delta|
| natural degree12 |11|10|1|1|0|
| exceptional degree11 |2|1|4|1..4|4-delta|
| exceptional degree11 |3|2|3|1..3|3-delta|
| exceptional degree11 |5|4|2|1 or2|2-delta|

These follow from full cycle dictionaries, not from assuming the
missing points are unramified. Together with k<=n-2 they give the
dimension claims in KMFI39.

For (1), rank two makes D0 dominate a surface and e=3.
No fixed-sheet permutation in the natural degree12 action has
a cycle of length3. The exceptional degree11 action has exactly
one such type, 3^3*1^2. KMFI32 bounds the present delta by2.
Rank at most one has the symmetric groups of KMFI38.
This proves all conclusions in KMFI39.

## 4. A rational pole forces a collapsed line

The derivative identities implied by (2) are

\[
 (X_1)_z=-\lambda_w-(\rho\lambda)_z,\quad
 (X_2)_z=w(X_1)_z+\lambda,\quad
 (X_3)_z=z(X_1)_z+\rho\lambda. \tag{5}
\]

At a pole w=a of order t>=1 in the reduced rational expression
for rho, polynomiality of rho^2*lambda forces
(w-a)^(2t) to divide lambda in k[w,z].
In the first identity of (5), the two summands are divisible
respectively by (w-a)^(2t-1) and (w-a)^t.
The last two identities then show all (Xi)_z vanish at w=a.
Thus X(a,z) is constant and S(0,a,z) collapses an entire line.
A finite morphism cannot have that fiber.

Consequently finite S in (2) forces rho polynomial.
KMFI33 makes it independent of z.
A constant rho makes im(N JX) a constant source line; more
directly the identities in KMFI28 with h'=0 give all Xi as
functions of sigma=h*w-z. Thus X has rank at most one, so
gamma=0 cannot map finitely onto a surface.
Hence all finite sweeps in this chart lie in the
nonconstant-polynomial-one-parameter setting of KMFI28.

## 5. The full Gauss cover of a non-affine Laurent multiplier

Suppose rho depends on z. KMFI33, after translation and scaling,
gives

\[
 \rho={\kappa z+R(w)\over w},\quad
 \lambda_{\rm top}=Cw^r z^s,\quad
 r(r-1)+2\kappa(r-1)(s+1)+
             \kappa^2(s+1)(s+2)=0. \tag{6}
\]

Here r>=2,s>=0 and R is Laurent polynomial.
The residue discriminant is negative rational, so
kappa!=0,1,2 and kappa is nonreal in its complex embedding.

X has rank two. Indeed KMFI38 makes any rank-one X a function
of a linear parameter U; the rank-one matrix N JX then has
image ker(dU), a constant line. In (2) that image is (1,rho),
which would force rho constant.

The normalized Gauss normal is

\[
 (\eta,\beta,1)=((\kappa-1)z+R,\;-(\kappa z+R)/w,\;1).
\]
\[
 k(\eta,\beta)\subset k(X)\subset k(w,z)=k(\beta,w),\qquad
 \eta={R(w)-(\kappa-1)\beta w\over\kappa}. \tag{7}
\]

An affine part r0+r1*w of R is removed by
z'=z+r0/kappa+r1*w/(kappa-1), with the matching linear
target change. This is a polynomial coordinate isomorphism;
the argument may use the transformed polynomial X.
If no R remains, (7) already recovers w,z rationally, giving delta1.

Otherwise R''!=0 and consider the rational Gauss cover in (7).
Its full group over k(eta,beta) is S_d, where d is its generic
degree. At generic fixed beta the finite critical points satisfy
R'(w)=(kappa-1)beta and are simple.
Their critical-value curve is

\[
 \beta={R'(w)\over\kappa-1},\qquad
 \eta={R(w)-wR'(w)\over\kappa},\qquad
 {d\eta\over d\beta}=-{\kappa-1\over\kappa}w.
\]

The derivative ratio recovers w, so the critical values are
generically distinct. Infinity and, if present, zero are poles
of the rational w-curve; there is no additional finite critical
value there. The finite inertia is single transpositions and
generates the transitive full S_d. This also prevents an
arithmetic enlargement of this already maximal group.

The natural S_d action has maximal point stabilizer, hence (7)
has no proper intermediate field.
KMFI32 bounds delta=[k(w,z):k(X)] by2.
If d>=3 this forces delta1. The only d2 cases, after removing
affine terms, are R=a*w^2 and R=a/w, a!=0.
They are settled next without a bounded group search.

## 6. Polynomial invariants rule out both quadratic exceptions

For R=a*w^2, source coordinates w,eta are polynomial coordinates:
eta=(kappa-1)z+a*w^2. The quadratic Gauss involution is

\[
 \tau(w)=-{\kappa\eta\over aw},\qquad \tau(\eta)=\eta .
\]

Any tau-invariant polynomial in k[w,eta] is in k[eta]:
a positive w power becomes a negative power, and their distinct
Laurent exponents cannot cancel.
If delta2, (7) gives k(X)=k(eta,beta), so each polynomial Xi
would be tau-invariant and X would have rank at most one.
This contradicts rank two.

For R=a/w, the Gauss involution is

\[
 \tau(w)={aw\over(\kappa-1)(\kappa wz+a)},\qquad
 \tau(z)=-(\kappa-1)z+
                {a(2-\kappa)\over(\kappa-1)w}. \tag{8}
\]

Along E={kappa*w*z+a=0}, tau(w) has a simple pole, whereas
tau(z) is finite with nonconstant residue
a/[kappa(kappa-1)w].
Thus a tau-invariant polynomial cannot have positive degree
in w: its highest w coefficient, a polynomial in z, cannot
vanish identically at that nonconstant residue.
It must be a polynomial in z only.
But then the pole of tau(z) at w=0 forces it constant.
So the invariant polynomial ring is k, again incompatible
with a rank-two polynomial X.

Therefore delta1 in every case. This proves the birationality
assertion of KMFI40 for every allowed positive z degree in lambda.
It does not assert that other normalization ramification branches
are absent. In the exceptional PSL2(11) case, the full
3^3*1^2 inertia has one present ramified branch and exactly two
ramified branches outside this source.

Universal rank, Gauss and valuation arguments are written here.
Exact controls are not formal verification or independent whole-proof review.
The abstract PSL2(11) action problems and full classification remain open.

[KMFI41--43](RATIONAL_CHART_CLASSIFICATION.md) now settle the
conditional boundary requirement above for this whole z-dependent
rational chart. The resultant has constant leading coefficient and,
after its exact extraneous w power is removed, constant term a
nonzero multiple of Y2^(s+2). Thus the sweep is finite on Y2!=0;
all missing branches lie over Y2=0, a different divisor from X(A2).
The two missing ramified branches above the ramification surface
cannot occur. The complete groups are natural alternating/symmetric.
The conditional necessity in KMFI40 remains valid; it is no longer
an open construction gate in this chart.
