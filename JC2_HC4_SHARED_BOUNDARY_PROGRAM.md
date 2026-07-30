# A shared boundary programme for \(\mathrm{JC}_2\) and \(\mathrm{HC}_4\)

## Status and first outcome

This is a research programme, not a proof of either conjecture.  It isolates
one exact common normal form and one candidate boundary obstruction.

The exact part is:

> A two-step isotropic cotangent chart in four Hessian variables is
> equivalent to a plane Keller map, including collisions.  Along an
> index-two finite-normalization boundary, the square root of the leading
> Hessian term is the reduced conormal residue.  In the completed quartic
> cusp and connector charts already computed in the plane programme, half
> of that residue is exactly the odd-square multiplier \(\ell\) whose zero
> scheme carries the conductor endpoints.

The open part is to attach such a polynomial isotropic flag canonically to
an arbitrary hypothetical \(\mathrm{HC}_4\) counterexample, or to construct
the corresponding residue-decorated flag on a compactification, and then
prove that its conductor-paired initial-residue class is nonzero.  None of
those recognition or nonvanishing statements is claimed here.

The exact symbolic identities are checked by
[`scripts/verify_jc2_hc4_isotropic_boundary_bridge.py`](scripts/verify_jc2_hc4_isotropic_boundary_bridge.py).

## 1. The exact two-step isotropic chart

Work over a characteristic-zero field.  Let

\[
 F=(P,Q):\mathbb A^2_{x,y}\longrightarrow\mathbb A^2
\]

and introduce dual variables \(t,m\).  For an arbitrary
\(H\in k[x,y]\), put

\[
 \Psi_F(x,y,t,m)=tP(x,y)+mQ(x,y)+H(x,y).             \tag{1.1}
\]

The plane

\[
 L=\langle\partial_t,\partial_m\rangle
\]

is a constant isotropic two-plane for the Hessian form.  In the ordered
coordinates \((x,y,t,m)\),

\[
 \operatorname{Hess}\Psi_F=
 \begin{pmatrix}
 A & (DF)^{\mathsf T}\\
 DF&0
 \end{pmatrix},
 \qquad
 A=t\,\operatorname{Hess}P
   +m\,\operatorname{Hess}Q+\operatorname{Hess}H.
\tag{1.2}
\]

The zero lower-right block has the same size as \(A\), so direct block
expansion gives

\[
 \boxed{\det\operatorname{Hess}\Psi_F=J(P,Q)^2.}     \tag{1.3}
\]

The equality is independent of \(A\).  Thus \(F\) is Keller exactly when
\(\Psi_F\) has constant nonzero Hessian determinant.

There is also an exact collision comparison.  The gradient is

\[
 \nabla\Psi_F=
 \bigl((DF)^{\mathsf T}(t,m)^{\mathsf T}+\nabla H,\,
       P,Q\bigr).                                    \tag{1.4}
\]

If \(F(a)=F(b)\) with \(a\ne b\), the invertibility of \(DF(a)\) and
\(DF(b)\) lets one choose unique dual vectors over \(a\) and \(b\) giving
the same first two entries of (1.4).  Hence \(\nabla\Psi_F\) has a
collision.  Conversely, a collision of (1.4) gives \(F(a)=F(b)\), and
if \(a=b\), invertibility of \(DF(a)\) forces the dual vectors to agree.
Therefore

\[
 \boxed{\nabla\Psi_F\text{ is injective }
 \Longleftrightarrow F\text{ is injective}.}         \tag{1.5}
\]

This is stronger than using only the implication
\(\mathrm{HC}_4\Rightarrow\mathrm{JC}_2\): inside the isotropic cotangent
chart the two collision problems are the same problem.

## 2. The Schur remainder that recognizes the cotangent branch

Start one step more generally with

\[
 \Psi=tP(x,y)+\Phi(x,y,m).                           \tag{2.1}
\]

Let

\[
 R(P)=
 (\nabla P)^{\mathsf T}
 \operatorname{adj}(\operatorname{Hess}P)\nabla P.
\tag{2.2}
\]

Expanding \(\det\operatorname{Hess}\Psi\) in \(t\) gives

\[
 [t^2]\det\operatorname{Hess}\Psi=0,\qquad
 [t]\det\operatorname{Hess}\Psi=-\Phi_{mm}R(P).
\tag{2.3}
\]

The first equality is the conormal shadow of \(P\) being independent of
the second isotropic coordinate.  If the full Hessian determinant is a
nonzero constant, (2.3) gives the exact dichotomy

\[
 \Phi_{mm}=0
 \quad\text{or}\quad
 R(P)=0.                                             \tag{2.4}
\]

The first branch integrates to

\[
 \Phi=mQ(x,y)+H(x,y)
\]

and is precisely (1.1).  The second is the binary bordered-degeneracy
branch which already appears in the support-free \(\mathrm{HC}_4\)
reductions.  Thus the positive \(t\)-coefficients of the determinant form
a small **isotropic Schur class**: its vanishing is forced by constant
Hessian determinant, and one irreducible factor of that vanishing produces
the plane Keller map.

This gives a concrete recognition target:

> **Isotropic-flag recognition problem.**  Show that every minimal
> four-variable constant-Hessian collision admits, after an allowed
> polynomial or boundary-adapted symplectic rechart, a two-step isotropic
> flag for which (2.4) reaches the cotangent branch; or classify the
> \(R(P)=0\) boundary so completely that it cannot carry the collision.

Existing coordinate reductions establish this pattern on large homogeneous
strata, but not for an arbitrary potential or a non-coordinate coisotropic
embedding.

## 3. Finite normalization appears without further choices

On the cotangent branch, put

\[
 A=k[P,Q],\qquad K=k(x,y),\qquad
 B=\operatorname{Norm}_{A}(K).
\tag{3.1}
\]

The plane finite-normalization theorem makes \(B\) finite free over \(A\).
The canonical factorization is

\[
 \mathbb A^2\hookrightarrow
 \bar X=\operatorname{Spec}B\longrightarrow\mathbb A^2. \tag{3.2}
\]

Consequently an isotropic reduction of a hypothetical \(\mathrm{HC}_4\)
counterexample does not merely output two polynomials.  It outputs the full
finite-normalization packet

\[
 \mathcal N(F)=
 \bigl(B/A,\ \mathbb A^2\subset\bar X,\
       \{E_i,\pi(E_i),e_i,f_i\}\bigr)                \tag{3.3}
\]

used by the plane programme.  Conversely, every hypothetical plane Keller
counterexample supplies (1.1), its constant isotropic flag, and the same
packet.

The four-dimensional package should therefore retain both languages:

\[
 \boxed{
 \mathcal H(\Psi,L)=
 \bigl(\Psi,L,\mathcal S_L;\,
       B/A,E_i,\rho_{E_i},\mathfrak c,
       \text{endpoint pairing}\bigr),}
\tag{3.4}
\]

where \(\mathcal S_L\) is the isotropic Schur class, \(\rho_E\) is the
reduced conormal residue defined below, and \(\mathfrak c\) is the relevant
normalization conductor.

## 4. The conormal residue is the boundary Hessian square root

Let \(E\) be a tame ramification prime of index \(e\) in (3.2).  At its
generic smooth point choose adapted source parameters \((r,z)\), a target
branch equation \(g\), and a target tangential parameter \(q\), normalized
so that \(q\circ\pi=z\bmod r\).  Write

\[
 g\circ\pi=u\,r^e,\qquad u\notin(r).                 \tag{4.1}
\]

After fixing the two local volume frames, define the reduced conormal
residue by

\[
 \rho_E=
 \left.
 r^{1-e}\frac{\pi^*(dg\wedge dq)}{dr\wedge dz}
 \right|_E.                                         \tag{4.2}
\]

Changing adapted frames multiplies \(\rho_E\) by a unit, so intrinsically
it is a nonzero rational section of the corresponding conormal/different
line on \(E\), not a preferred scalar.

Equation (1.3) shows that the leading boundary coefficient of the
cotangent Hessian is

\[
 \rho_E^2.                                          \tag{4.3}
\]

This is the desired common carrier: the Hessian side naturally sees its
square, while the finite-normalization side sees the unsquared different
residue and its endpoint behavior.

## 5. The quartic conductor chart identifies the carrier

The completed one-boundary quartic cusp and connector packets have the
adapted factorization

\[
 g=r^2\ell.                                         \tag{5.1}
\]

Hence

\[
 \left.r^{-1}\frac{\partial g}{\partial r}\right|_{r=0}
 =2\ell|_{r=0}.                                     \tag{5.2}
\]

Thus

\[
\boxed{\rho_E=2\ell}
\]

in the adapted frames.  In the clean cusp coordinates already used by the
quartic packet,

\[
 \ell=4r-9T^2,
\]

so

\[
 \rho_E|_E=-18T^2.                                  \tag{5.3}
\]

The cusp endpoint therefore has conormal contact two and initial coefficient
\(-18\) in this normalization.  At a connector branch with
\(\ell|_E=c\tau^m+\cdots\), the corresponding initial residue is
\(2c\tau^m\).  The target-side quadratic base change is

\[
 as^2=r^2\ell,
\]

and its normalization adjoins \(z=r\ell/s\), with

\[
 z^2=a\ell,\qquad
 \mathfrak c=(r,s),\qquad
 N/O\simeq O/(r,s).                                 \tag{5.4}
\]

The same \(\ell\) therefore has three meanings:

1. half of the reduced conormal residue;
2. the square root, up to the factor \(2\), of the leading cotangent-Hessian
   coefficient; and
3. the odd-square multiplication coefficient whose zero divisor is the
   affine-companion/conductor-endpoint carrier.

This identification is exact.  It explains why a useful common obstruction
must retain a square root and endpoint residues: the scalar Hessian
determinant sees only \(4\ell^2\), while the conductor pairing distinguishes
the individual initial forms of \(\ell\).

The known transition laws make the relevant character explicit.  If

\[
 r_i=u_{ij}r_j,\qquad \ell_i=u_{ij}^{-2}\ell_j,
\]

then

\[
 \rho_i=u_{ij}^{-2}\rho_j.                           \tag{5.5}
\]

Thus the conormal residue carries the inverse-square boundary character.
The normalized generator \(z_i=u_{ij}^{-1}z_j\) is its primitive
half-character, while the scalar Hessian leading term carries the
inverse-fourth character.

## 6. Candidate obstruction: the paired initial-conormal class

Ordinary values of \(\rho_E\) vanish at the cusp and connector endpoints,
so the unfiltered conductor quotient forgets the required information.
Retain instead the first nonzero normal/conductor initial form

\[
 \operatorname{in}_p(\rho_E)
\]

at every endpoint \(p\).  If the conductor identifies endpoint branches
\(p\sim p'\), use its induced identification of the conormal/different
lines and form the mismatch

\[
 \delta_{p,p'}(\rho_E)=
 \operatorname{in}_p(\rho_E)
 -\tau_{p,p'}\operatorname{in}_{p'}(\rho_E).         \tag{6.1}
\]

Collecting these mismatches gives a class

\[
 \operatorname{Obs}_{\mathrm{pair}}(\Psi,L)
 \in
 \operatorname{coker}\!\left(
 \operatorname{gr}_{\mathfrak c}\mathcal L
 \longrightarrow
 \nu_*\operatorname{gr}_{\widetilde{\mathfrak c}}
             \nu^*\mathcal L\right),                \tag{6.2}
\]

where \(\mathcal L\) is the conormal/different line carrying \(\rho_E\).
Formula (6.2) is schematic until the global packet fixes the filtration and
the transition character of \(\mathcal L\).

The intended vanishing mechanism is the same on both sides:

- a plane Keller counterexample supplies a global finite-normalization
  different section, so its initial residues must descend through the
  actual conductor pairing;
- an \(\mathrm{HC}_4\) counterexample in the recognized isotropic chart
  supplies the same section as the square root of the boundary-leading
  Hessian term, and polynomial Hessian gluing forces the same descent.

Thus either counterexample would force

\[
 \operatorname{Obs}_{\mathrm{pair}}(\Psi,L)=0.       \tag{6.3}
\]

What is not proved is that every \(\mathrm{HC}_4\) counterexample reaches
this chart, that the square root has a canonical sign/character globally,
or that every surviving quartic packet makes (6.2) nonzero.  These are the
three genuine research gates.

## 7. Why coarser candidates cannot work

The existing calculations already exclude several tempting shortcuts.

1. The generic residue \(\rho_E\) is nonzero by tame ramification, so the
   desired statement cannot be generic vanishing along \(E\).
2. The completed cusp and connector rings have the same normal
   determinantal overring, conductor \((r,s)\), and canonical module.
   The conductor ideal alone does not distinguish them.
3. The primitive boundary/different character glues compatibly across the
   packet.  Its divisor class alone gives no contradiction.
4. The scalar Hessian determinant records \(\rho_E^2\), losing the sign,
   transition character, and conductor pairing of the initial residue.

The first viable obstruction is therefore filtered, square-root-sensitive,
and pairing-sensitive, as in (6.2).

## 8. Immediate calculations

The programme now has four narrow next steps.

1. **Completed comparison.**  Compute (6.2) in the standard \(3+1\) cusp
   chart and in both branches of a \(2+2\) connector, using
   \(\rho_E=2\ell\) and the transition laws
   \(r_i=u_{ij}r_j\), \(\ell_i=u_{ij}^{-2}\ell_j\).
2. **Sign/character descent.**  The cotangent flag selects the
   inverse-square different residue \(\rho_i\), while the scalar Hessian
   remembers only its inverse-fourth-character square.  Decide whether an
   arbitrary Hessian boundary package canonically selects that root, and
   whether the primitive inverse character carried by \(z_i\) is also
   needed to compare paired endpoints.
3. **Schur-to-boundary comparison.**  Express the first nonzero filtered
   isotropic Schur remainder in the same associated-graded conormal line
   as \(\rho_E\).  This is the missing map from (2.3) to (6.2).
4. **Packet nonvanishing.**  Test the actual cusp/connector endpoint
   pairing, not an arbitrary matching, for a forced residue mismatch.  A
   nonzero class would simultaneously exclude the plane packet and its
   four-dimensional cotangent Hessian lift.

The first and third steps are local symbolic calculations.  The fourth is
global and must use the degree-four monodromy and the two-generator
degree-zero algebra; the existing no-finiteness examples show that an
unbounded abstract conductor matching is insufficient.

## 9. Reproduction

Run

```bash
.venv/bin/python scripts/verify_jc2_hc4_isotropic_boundary_bridge.py
```

The checker verifies (1.3), (2.3), (5.2), and the cusp specialization
(5.3).  The finite-normalization, completed conductor, and odd-square
statements used in Sections 3 and 5 remain canonically sourced in
[`plane-jc/FINITE_NORMALIZATION_PROGRAM.md`](plane-jc/FINITE_NORMALIZATION_PROGRAM.md)
and
[`plane-jc/JC2_GLOBAL_COX_PACKET_ATTACK.md`](plane-jc/JC2_GLOBAL_COX_PACKET_ATTACK.md);
they are not re-proved by this small checker.
