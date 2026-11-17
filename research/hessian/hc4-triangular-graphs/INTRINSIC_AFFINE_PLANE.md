# An intrinsic plane and a three-direction HC4 criterion

14 November 2026. The previous parallel-trace theorem extends across
different hyperplane orientations. A hypothetical nonautomorphic
four-variable unit-Hessian polynomial has one **unique** constant
jointly affine two-plane. Its annihilator is exactly the two-dimensional
space of covectors admitting even one pure Hodge trace.

Consequently three linearly independent pure-trace directions force a
polynomial gradient inverse, without any sign condition on their Hodge
coefficients, graph shape or degree bound. Two such directions already
determine the only possible plane, which can be tested by three exact
directional second-derivative identities. The remaining cotangent case
retains precisely its plane Keller problem.

The uniqueness argument below is elementary and does not require Tang.
Existence of a plane under nonautomorphy uses the retained reproduced
Tang theorem. Hodge recognition uses the retained single-trace and
parallel-trace arguments; no new external classification is imported.

## 1. Two distinct jointly affine planes force invertibility (HC4MYHF0)

Let f be a polynomial on a four-dimensional vector space over any
characteristic-zero field, with nonzero constant Hessian determinant.
A constant two-dimensional vector plane W is **jointly affine** if

\[
              D_w D_{w'}f=0\quad\hbox{for all }w,w'\in W.    \tag{1}
\]

This means f is affine on every translate of W. An individual direction
with D_w^2 f=0, or one translate on which f happens to be affine, does
not by itself establish (1).

If two distinct such planes W,W' exist, **grad f is a polynomial
automorphism**. Their intersection has dimension zero or one.

If it is zero, take their direct sum as source coordinates. Joint
affineness on both blocks makes f of degree at most one in each block,
hence of total degree at most two. Its invertible constant Hessian
makes its gradient an affine-linear automorphism.

If the intersection is a line, choose constant coordinates (s,U,P,Q)
so that W=span(partial_P,partial_Q) and
W'=span(partial_P,partial_U). The first plane gives

\[
                 f=D(s,U)+P F(s,U)+Q G(s,U).                \tag{2}
\]

The second imposes F_U=0 and D_UU=G_UU=0, so

\[
 F=F(s),\quad G=g_1(s)U+g_0(s),\quad
 D=d_1(s)U+d_0(s).                                        \tag{3}
\]

For any form (2), the exact block identity is

\[
               \det\operatorname{Hess}f=J(F,G)^2.           \tag{4}
\]

Here J(F,G)=F'(s)g_1(s). Since (4) is a nonzero constant, both factors
are nonzero field constants. Thus F=alpha*s+beta and
G=gamma*U+g_0(s), alpha*gamma nonzero: the plane pair is a triangular
polynomial automorphism.

More generally, whenever the pair in (2) has a polynomial inverse,
the gradient outputs F,G recover s,U. Its other two outputs, after
subtracting grad D, recover P,Q using the inverse transpose of the
unit plane Jacobian. This proves the assertion.

In particular, every nonautomorphic unit-Hessian f has **at most one**
jointly affine two-plane, without any external structural theorem.

## 2. Pure traces recognize exactly that plane (HC4MYHF1)

Now work over C and assume the gradient of f is not an automorphism.
Retained [HCTANG2](../tang-hc4-audit/STRUCTURAL_REVIEW.md) supplies a
jointly affine two-plane W. HC4MYHF0 makes it unique.

For a nonzero constant source covector ell, call its direction **pure**
if for at least one scalar a the entire generic value fibre of the trace
on the actual affine hyperplane ell=a has

\[
                  E(f|_{\ell=a}=\lambda)=t^2+\kappa t,
                           \qquad t=uv.                    \tag{5}
\]

Neither the integer kappa nor the hyperplane offset is prescribed.
Write P_f for the set of such nonzero covectors. Then

\[
                         P_f\cup\{0\}=W^\perp.              \tag{6}
\]

To prove this, adapt the same plane W to the free affine coordinate
s=ell. If ell annihilates W, the rank-zero form is (2), after a constant
linear change of quotient coordinates. At each scalar s=a its trace
is polynomially equivalent to D_a(U)+g_a(U)V with a free W-coordinate,
where g_a is the gcd of the two specialized coefficient rows. Thus

\[
                 E(f|_{\ell=a}=\lambda)=t^2-m_{\ell,a}t,
                         \qquad m_{\ell,a}\ge0.              \tag{7}
\]

This holds for every a, so ell belongs to P_f.

If ell does not annihilate W, its rank on W is one. The same
[HC4MYHC0 rank-one argument](HODGE_TRACE_CLOSURE.md#2-a-single-trace-criterion-with-no-critical-point-hypothesis-hc4myhc0)
used in [HC4MYHE0](TWO_TRACE_COMPATIBILITY.md#1-a-pure-trace-fixes-the-whole-parallel-profile-hc4myhe0)
writes the trace polynomial as

\[
 t^2+(d_a-c)t+G(u+v)+(R-c),                                \tag{8}
\]

with a fixed smooth plane Keller curve having c components, total
completion genus G and total punctures R. If any such trace were pure,
G=0 and R=c would force actual closed A1 components. HC4MYGR0 then
makes the plane Keller pair and the full gradient automorphisms,
contradicting the hypothesis. Thus ell is not in P_f. This proves (6).

The projectivized pure directions are therefore exactly one projective
line. This assertion is conditional on nonautomorphy. For an automorphic
gradient the pure directions need not form a projective line.

### Three independent directions close arbitrary graphs and mixing

Let ell_1,ell_2,ell_3 be three **linearly independent** source covectors.
Suppose for each i there is one actual affine hyperplane ell_i=a_i
whose entire generic value fibre has form (5). Then **unit Hessian
forces a polynomial gradient inverse**. Under hypothetical
nonautomorphy, all three covectors would lie in the two-dimensional
space W^perp, an impossibility.

All three kappa_i may be zero or negative. Pairwise nonparallel is
insufficient: ell_1,ell_2,ell_1+ell_2 still span only two dimensions.
The offsets and generic value parameters may differ. The statement
applies to any polynomial graph pullback, actual embedding or broader
polynomial substitution for which these traces are established in the
actual Hessian variables; no global substitution inverse is needed.

### Two directions determine a concrete remaining algebraic test

For two independent pure directions, put

\[
                  W_0=\ker\ell_1\cap\ker\ell_2.             \tag{9}
\]

Choose a basis w_1,w_2 of W_0. If at least one of the three polynomial
identities

\[
       D_{w_1}^2f=0,\qquad D_{w_1}D_{w_2}f=0,
                         \qquad D_{w_2}^2f=0                \tag{10}
\]

fails, unit Hessian already forces gradient invertibility: (6) would
make W_0 the required jointly affine plane under nonautomorphy.
If all identities hold, constant coordinates produce (2), and (4)
gives a plane Keller pair. This is an exact cotangent reduction, not a
properness theorem for that pair.

The universal assertion that **two independent pure directions alone**
always imply invertibility would be equivalent to JC2: every plane
Keller pair has a cotangent lift (2) with the two quotient directions
pure. Invertibility of its full gradient implies a plane inverse. To
see the latter directly, set the two complementary gradient targets
to zero and read the quotient coordinates of the full inverse; these
give both compositions of a polynomial inverse to (F,G). Conversely
JC2 implies HC4 by HCTANG2. The third-direction criterion therefore
does not remove the known unresolved plane endpoint.

### Stability under affine-linear additions

Adding any affine-linear polynomial to f leaves its Hessian and
jointly affine planes unchanged, and translates its gradient. Thus,
under nonautomorphy, (6) proves that the entire pure-direction set is
unchanged by such additions. No invariance under nonlinear source
changes is inferred. Constant affine source changes transport the
annihilator relation covariantly.

## 3. Parallel profiles count the exact residual fibre (HC4MYHF2)

In (2), let Z=V(F,G) in the quotient affine (s,U) plane. The Keller
condition makes Z finite and reduced, possibly empty. For every pure
covector ell on the quotient and every scalar a,

\[
                    m_{\ell,a}=\#(Z\cap\{\ell=a\}).         \tag{11}
\]

Indeed, the roots of the specialized gcd g_a are exactly those common
zeros. Each root is simple: if both U derivatives vanished there,
J(F,G) would vanish. Thus (7) counts precisely these reduced points,
with no hidden multiplicities. In particular, for any one nonzero ell,

\[
                    \sum_{a\in C}m_{\ell,a}=\#Z,            \tag{12}
\]

and the sum has finite support and is independent of ell. The same
integer is the number of vertical A1 components in the entire relative
critical scheme for that base.

The full profiles over all ell also determine Z itself. A point z
belongs to Z exactly when every affine line through z has positive
count in (11). If z is not in the finite set, choose a line through z
avoiding it. Because Z is reduced, recovering this set recovers its
scheme. The quotient and the ideal are intrinsic up to constant affine
quotient coordinates and an invertible constant basis change of W:
Z is the vanishing locus of df restricted to W.

In fact, after n=#Z is obtained from one complete profile, any n+1
pairwise nonparallel quotient covectors suffice to reconstruct Z.
Intersect their positive-count line unions. A point outside Z belongs
to at most n of them: the corresponding n+1 distinct lines through
that point would otherwise require n+1 distinct points of Z, since each
point of Z lies on at most one of those lines. Points of Z belong to
every union. This is a finite-set reconstruction statement, not a claim
that an unknown plane Keller map realizes a prescribed finite set.

Consequently zero total profile means (F,G)=(1), while total profile
one means (F,G) is the maximal ideal of one reduced point. These are
the existing distinct unit-ideal and singleton-zero endpoints. A single
coefficient m_{ell,a}=0 does not mean the whole profile is zero; a single
coefficient one does not count zeros at other offsets. Neither total
one nor an ideal equality supplies properness.

The affine-linear perturbations f-alpha*P-beta*Q replace the pair by
(F-alpha,G-beta). Their planes and pure directions are unchanged.
Applying (11),(12) to each perturbation counts every fibre of the same
plane Keller map. This identifies exactly what these profiles measure;
it does not assert any particular fibre-cardinality pattern for an
unknown nonautomorphic plane pair.

## Scope, descent and checks

HC4MYHF0 is elementary over any characteristic-zero field. The
three-direction and two-direction-with-failed-identity consequences
over such a field follow by embedding the finitely generated coefficient
field into C and descending a polynomial inverse by uniqueness. Hodge
and profile descriptions above are stated over C. No source hyperplane
is straightened nonlinearly, and no finite-cover Hodge multiplicativity
is used.

[verify_intrinsic_plane.py](verify_intrinsic_plane.py) and its
[receipt](receipt-intrinsic-plane.json) check both relative plane positions,
their exact Hessian/gradient inverses, covector ranks and annihilators,
nonpositive pure-trace controls, affine additions and finite-set profiles.
The arbitrary-degree proofs and imported theorems are separate from
these finite controls. Multiple-point finite-set controls do not assert
the existence of a plane Keller map with that origin fibre.

Full HC4 and JC2 remain open, as do the unit-ideal surjectivity and the
distinct singleton-zero properness endpoints. No formal verification,
independent replay, external review or novelty priority is claimed.
