# Confluent degeneration and collision persistence protocol

## Status and purpose

This is a reusable research protocol, not a new Keller-map theorem or a
construction.  It records the part of the degeneration method in
Arathoon--Ball--Kvalheim, *The Maxwell Conjecture is False*,
arXiv:2607.27197, that can transfer to polynomial Keller-map searches:

\[
 \text{tune a leading cancellation}
 \longrightarrow
 \text{rescale to a finite limit}
 \longrightarrow
 \text{certify a separated local model}
 \longrightarrow
 \text{persist it in the original family}.
\]

The protocol must not be cited as showing that a proposed ansatz has a
collision.  The determinant identity, the convergence of the rescaled family,
and the separated limit collision are independent obligations.

## 1. Local collision-persistence lemma

Let \(S\) be a neighbourhood of \(0\) in \(\mathbb R^a\) or
\(\mathbb C^a\).  Let

\[
 G_s:U\longrightarrow \mathbb K^n \qquad (s\in S)
\]

be a \(C^1\) family on an open set \(U\subset\mathbb K^n\), where
\(\mathbb K=\mathbb R\) or \(\mathbb C\).  Suppose there are distinct
points \(X_0,Y_0\in U\) and a value \(b_0\) such that

\[
 G_0(X_0)=G_0(Y_0)=b_0,
 \qquad
 \det DG_0(X_0)\det DG_0(Y_0)\ne0.
 \tag{1}
\]

Then, after shrinking parameter and target neighbourhoods, there are two
disjoint local inverse branches

\[
 X(s,b),\quad Y(s,b),
 \qquad
 G_s(X(s,b))=G_s(Y(s,b))=b.
 \tag{2}
\]

In particular, every sufficiently small \(s\) has a collision near
\((X_0,Y_0)\).  This is the inverse-function theorem applied separately at
\(X_0\) and \(Y_0\), followed by choosing disjoint source neighbourhoods.
The statement is elementary, but it is the precise replacement for an
unjustified assertion that a limiting collision "obviously persists."

For a Keller family, \(\det DG_s\) is automatically nonzero wherever the
rescaled map is defined.  Nevertheless, (1) must still be checked for the
limit map: a singular limiting collision has no such persistence conclusion.

## 2. Rescaled-family form

The useful setup is a polynomial or analytic family \(F_\varepsilon\) for
\(\varepsilon\ne0\), together with invertible source and target scalings
\(S_\varepsilon,T_\varepsilon\), and a target translation \(c_\varepsilon\).
Define

\[
 G_\varepsilon(X)=
 T_\varepsilon^{-1}
 \bigl(F_\varepsilon(S_\varepsilon X)-c_\varepsilon\bigr).
 \tag{3}
\]

The required analytic certificate is that (3) extends to a \(C^1\) family
at \(\varepsilon=0\), preferably with an explicit expansion

\[
 G_\varepsilon=G_0+O(\varepsilon^\delta)
 \quad\text{in }C^1(K)
 \tag{4}
\]

on compact source sets containing the two limit points.  If \(G_0\)
satisfies (1), the lemma gives distinct points \(X_\varepsilon,Y_\varepsilon\)
with equal \(G_\varepsilon\)-value.  Since both scalings are invertible for
\(\varepsilon\ne0\),

\[
 F_\varepsilon(S_\varepsilon X_\varepsilon)
 =F_\varepsilon(S_\varepsilon Y_\varepsilon),
 \qquad
 S_\varepsilon X_\varepsilon\ne S_\varepsilon Y_\varepsilon.
 \tag{5}
\]

Thus the rescaled collision is an actual collision of the original map.

If the aim is a determinant-one family, verify that identity before taking a
limit.  The chain rule only relates the Jacobian determinant of
\(G_\varepsilon\) to that of \(F_\varepsilon\); it does not manufacture a
Keller identity from a promising limit model.

## 3. Construction checklist

1. Choose a symmetry or boundary stratum that forces a low-order coefficient
   which can be cancelled by one or more controlled parameters.
2. Solve the Jacobian/determinant ledger identically in those parameters.
   Do not infer it from finitely many parameter values.
3. Select source and target weights, and calculate (3) through the first
   nonzero order.  Record the exact parameter tuning that makes the desired
   finite limit appear.
4. Solve the limit collision exactly and certify both Jacobians are
   invertible.  Distinctness of the two source points is part of the
   certificate.
5. Prove the uniform \(C^1\) estimate (4), or an equivalent analytic
   extension argument, on two disjoint source neighbourhoods.
6. Apply the persistence lemma and translate the two branches back using
   (5).
7. Separately audit polynomiality, boundary denominators, and every claimed
   generic-fibre or nonproperness property.  Persistence supplies only the
   local collision.

For an exact algebraic presentation, Steps 3--5 may alternatively be done
over a Henselian local ring: simple roots of the two local inverse equations
lift uniquely.  That is a replacement for, not an automatic consequence of,
the complex analytic argument.

## 4. Repository targets

The highest-value applications are:

| Target | What the protocol could contribute | Extra obligation |
|---|---|---|
| [local coefficient components and decorated stable moduli](JELONEK_COEFFICIENT_COMPONENTS.md) | A collision open set inside an explicitly constructed non-automorphism component | Build a determinant-one family that is not merely a source reparametrization |
| [cancellation contact resultants](../cancellation/CONTACT_RESULTANT.md) | A tuned Rees degeneration whose simple reconstructed roots lift from a boundary model | The root and valuation ledgers must remain polynomial after lifting |
| [controlled boundary suspensions](../cancellation/CONTROLLED_BOUNDARY_SUSPENSIONS.md) | A source/target scaling that converts a boundary-cancelled formula into a finite collision model | Verify the polynomiality gate independently of the local model |
| [tagged GVC(3) lift](THREE_VARIABLE_GVC_TAGGED_LIFT.md) | Leading-face coordinates and finite-prefix screening | No generic perturbation: it destroys the infinite pure-moment premise |

The first three are genuine collision-persistence applications.  The GVC row
uses only the degeneration/screening half of the protocol.

## 5. Boundaries of the method

This method does not overcome the existing source-triviality issue for naive
polynomial Keller deformations.  If a family is only a source-coordinate
reparametrization, its collision behaviour is transported rather than newly
constructed.  A useful application therefore needs deformation of decorated
incidence or boundary data, with a separate proof that the relevant family is
not exhausted by that trivial action.

Nor does the final generic-perturbation step from the electrostatic paper
transfer wholesale.  It can remove accidental degeneracies of a finite
critical-point problem, but it usually destroys exact identities such as
\(\Lambda^m(P^m)=0\) for all \(m\).  Use it only where the target property is
open under perturbation, such as the separated local collision in (1).
