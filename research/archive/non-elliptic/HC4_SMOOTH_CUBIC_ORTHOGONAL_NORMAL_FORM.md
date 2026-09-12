# Smooth-cubic orthogonal normal form

## Status

This note proves `HC4NHM19`.  It is the first invariant reduction of the
clean quartic-denominator packet when the denominator has a smooth cubic
component.  It gives four intrinsic degree packets and excludes one of them.
It does **not** yet exclude the remaining three packets and therefore does
not prove that smooth cubic components are impossible.

Replay the universal matrix identities and degree ledger with

```bash
.venv/bin/python scripts/verify_hc4_smooth_cubic_orthogonal_normal_form.py
```

The line-bundle splitting and the final Hessian-integrability exclusion are
the written geometric argument below.

## 1. Clean restriction to the elliptic curve

Let

\[
 C=\operatorname{Hess}(h_5),\qquad
 \det C=P^2\ell,
 \tag{1.1}
\]

be a clean generic-corank-one packet from `HC4NHM1`.  Suppose

\[
 P=qm,
 \tag{1.2}
\]

where \(Q=V(q)\) is a smooth plane cubic and \(m\) is linear.  After base
change to an algebraic closure, the clean kernel condition is

\[
 \mathcal K_Q(3)\simeq\mathcal O_Q,
 \qquad
 \mathcal K_Q\simeq\mathcal O_Q(-3).
 \tag{1.3}
\]

Put

\[
 0\longrightarrow\mathcal O_Q(-3)
 \longrightarrow\mathcal O_Q^3
 \xrightarrow{\pi}\mathcal F\longrightarrow0.
 \tag{1.4}
\]

The defect-free hypothesis says that the form induced by \(C|_Q\) is a
symmetric isomorphism

\[
 b:\mathcal F\xrightarrow{\sim}\mathcal F^\vee(3).
 \tag{1.5}
\]

Taking determinants in (1.4) gives

\[
 \det\mathcal F\simeq\mathcal O_Q(3).
 \tag{1.6}
\]

This is the relevant line-bundle datum.  It is not directly a classical
linear symmetric determinantal representation of the cubic: the bundle in
(1.5) has rank two and the entries of the ambient Hessian have degree three.
Classical symmetric linear representations are instead controlled by
non-effective theta characteristics; see
[Beauville, *Determinantal hypersurfaces*, Section 4](https://math.univ-cotedazur.fr/u/beauvill/pubs/det.pdf).

## 2. Orthogonal splitting

For any rank-two bundle, wedge product gives the canonical alternating
isomorphism

\[
 j:\mathcal F\xrightarrow{\sim}
 \mathcal F^\vee\otimes\det\mathcal F.
 \tag{2.1}
\]

Use (1.6) and compare \(b\) with \(j\).  There is an automorphism
\(T\in\operatorname{Aut}(\mathcal F)\) such that

\[
 b(u,v)=j(Tu,v).
 \tag{2.2}
\]

In rank two, symmetry of \(b\) is equivalent to

\[
 \operatorname{tr}T=0.
 \tag{2.3}
\]

The determinant of \(T\) is a nonzero global function on the projective
integral curve, hence a nonzero scalar.  Cayley--Hamilton gives

\[
 T^2=-\det(T)I.
 \tag{2.4}
\]

Over the algebraic closure, \(T\) has two distinct constant eigenvalues.
The two projectors split \(\mathcal F\) into line bundles:

\[
 \boxed{
 \mathcal F\simeq L\oplus M,
 \qquad LM\simeq\mathcal O_Q(3).
 }
 \tag{2.5}
\]

Both summands are isotropic for \(b\), and \(b\) pairs them perfectly.

Write the two components of \(\pi\) as basepoint-free triples

\[
 \alpha=(\alpha_1,\alpha_2,\alpha_3)\in H^0(Q,L)^3,
 \qquad
 \beta=(\beta_1,\beta_2,\beta_3)\in H^0(Q,M)^3.
 \tag{2.6}
\]

After scaling the perfect pairing, the complete restriction normal form is

\[
 \boxed{
 C_{ij}|_Q=\alpha_i\beta_j+\beta_i\alpha_j.
 }
 \tag{2.7}
\]

The clean kernel generator is the cross product

\[
 e|_Q=\alpha\times\beta,
 \tag{2.8}
\]

whose entries are sections of \(LM=\mathcal O_Q(3)\).  In particular,

\[
 \operatorname{adj}(C)|_Q=-ee^{\mathsf T}
 \tag{2.9}
\]

after the same normalization.  This is an invariant normal form, not a
coefficient slice.

## 3. Finite degree packets

The quotient map in (1.4) is surjective at every point, so each triple in
(2.6) generates its line bundle.  A globally generated line bundle on an
elliptic curve has degree zero only when it is trivial, cannot have degree
one, and is globally generated in every degree at least two.  Since

\[
 \deg L+\deg M=\deg\mathcal O_Q(3)=9,
 \tag{3.1}
\]

the unordered degree pair is one of

\[
 \boxed{(0,9),\ (2,7),\ (3,6),\ (4,5).}
 \tag{3.2}
\]

Thus the elliptic Picard parameter has not disappeared, but it is confined
to four discrete degree packets.  The last three retain respectively a
degree-two pencil, a degree-three plane model, and a degree-four linear
series on \(Q\).

## 4. Exclusion of the trivial summand

Suppose the degree pair is \((0,9)\).  Then \(L\simeq\mathcal O_Q\), and a
constant change of the three ambient coordinates makes

\[
 \alpha=(1,0,0).
 \tag{4.1}
\]

Equation (2.7) gives

\[
 C_{22}|_Q=C_{23}|_Q=C_{33}|_Q=0.
 \tag{4.2}
\]

These are cubic forms and \(q\) is an irreducible cubic, so there are
constants \(a,b,c\) with

\[
 C_{22}=aq,\qquad C_{23}=bq,\qquad C_{33}=cq.
 \tag{4.3}
\]

Hessian integrability gives

\[
 a q_z=b q_y,
 \qquad
 b q_z=c q_y.
 \tag{4.4}
\]

For a smooth plane cubic, \(q_y\) and \(q_z\) are linearly independent
after every constant coordinate change.  Otherwise a nonzero constant
direction annihilates \(q\), so \(q\) is a binary cubic and is reducible
over the algebraic closure.  Equation (4.4) therefore forces

\[
 a=b=c=0.
 \tag{4.5}
\]

The lower \(2\)-by-\(2\) block of \(C\) then vanishes identically, making
\(\det C=0\), a contradiction.

## 5. Result and exact frontier

> **Theorem `HC4NHM19` -- Smooth-cubic orthogonal normal form.**  A clean
> smooth cubic component in the quartic denominator forces the restricted
> rank-two quotient bundle to split orthogonally as in (2.5), with Hessian
> restriction (2.7).  Its unordered isotropic degrees are exactly the four
> pairs in (3.2).  Hessian integrability excludes the pair \((0,9)\).

The target theorem “no smooth irreducible cubic component occurs” is
therefore reduced to the three packets

\[
 (2,7),\qquad(3,6),\qquad(4,5).
 \tag{5.1}
\]

The next invariant calculation should impose the ambient Hessian Codazzi
identities on (2.7), packet by packet.  In degree two, the triple \(\alpha\)
has one constant relation \(v\), and (2.7) gives the concrete necessary
condition

\[
 D_v^2h_5=cq
 \tag{5.2}
\]

for a scalar \(c\).  The datum is the pair \((Q,[v])\), not just the
elliptic \(j\)-invariant; normalizing \(Q\) to a Hesse cubic does not remove
the remaining direction moduli.  In degree three \(L\) is a plane cubic
model; in degree four its first multiplication relations give the smallest
remaining syzygy presentation.  None of these three cases is claimed empty
here.
