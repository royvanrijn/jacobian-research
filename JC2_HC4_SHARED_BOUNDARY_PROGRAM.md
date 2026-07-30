# A shared boundary programme for \(\mathrm{JC}_2\) and \(\mathrm{HC}_4\)

## Status and outcome

This is a research programme, not a proof of either conjecture.  It isolates
one exact common normal form and computes the first candidate boundary
obstruction.

The exact part is:

> A two-step isotropic cotangent chart in four Hessian variables is
> equivalent to a plane Keller map, including collisions.  Along an
> index-two finite-normalization boundary, the square root of the leading
> Hessian term is the reduced conormal residue.  In the completed quartic
> cusp and connector charts already computed in the plane programme, half
> of that residue is exactly the odd-square multiplier \(\ell\) whose zero
> scheme carries the conductor endpoints.  However, the proposed paired
> initial-conormal class in the associated-graded conductor cokernel is
> identically zero: at a cusp and in positive degree at a node, the
> conductor map to the normalization is an isomorphism.  At the node only
> degree-zero values are paired, and both residue values vanish.

The open part is to attach such a polynomial isotropic flag canonically to
an arbitrary hypothetical \(\mathrm{HC}_4\) counterexample, or to construct
the corresponding residue-decorated flag on a compactification.  Any
nonzero endpoint obstruction must add global jet transport supplied by the
two-generated degree-zero algebra, a dualizing-residue structure, or an
equivalent datum.  The conductor pairing and sheet monodromy alone do not
provide that transport.

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

More precisely, take \(F=(g,q)\) in the adapted coordinates.  If \(e=2\),
then (4.2) says

\[
 J(g,q)=r\rho_E\pmod {r^2}.
\]

Consequently

\[
 \operatorname{in}_{(r)}
 \bigl(\det\operatorname{Hess}(tg+mq+H)\bigr)
 =
 [r]^2\rho_E^2.                                    \tag{4.4}
\]

The cotangent volume frames select the oriented square root
\([r]\rho_E\); division by the normal generator \([r]\) gives precisely
\(\rho_E\).  Thus the filtered cotangent Hessian and the reduced conormal
residue occupy the same associated-graded line, with the former seeing the
square of the latter.  This is the exact Schur-to-boundary comparison on
the cotangent branch selected by \(\Phi_{mm}=0\) in (2.4).  It does not
assert that the general coefficient \(-\Phi_{mm}R(P)\) is itself a boundary
residue before that branch has been selected.

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

## 6. The paired initial-conormal class

Let \(R_p\subset\widetilde R_p\) be the completed boundary-curve ring and
its normalization at a packet endpoint, and let \(\mathfrak c_p\) be their
conductor.  After trivializing the inverse-square different line
\(\mathcal L\), the proposed class is

\[
 \operatorname{Obs}_{\mathrm{pair}}(\Psi,L)
 =
 [\operatorname{in}(\rho_E)]
 \in
 \operatorname{coker}\!\left(
 \operatorname{gr}_{\mathfrak c_p}(R_p)\otimes\mathcal L
 \longrightarrow
 \operatorname{gr}_{\mathfrak c_p\widetilde R_p}
       (\widetilde R_p)\otimes\mathcal L\right).     \tag{6.1}
\]

This makes the schematic formula precise, but the completed calculation
shows that its hoped-for positive-degree mismatch does not exist.

### 6.1 The \(3+1\) cusp

For an ordinary cusp,

\[
 R_{\rm cusp}=k[[\tau^2,\tau^3]]
 \subset
 \widetilde R_{\rm cusp}=k[[\tau]],
 \qquad
 \mathfrak c=(\tau^2,\tau^3)=\tau^2\widetilde R_{\rm cusp}.
\tag{6.2}
\]

For every \(n\ge1\),

\[
 \mathfrak c^n/\mathfrak c^{n+1}
 \longrightarrow
 (\tau^{2n})/(\tau^{2n+2})
\tag{6.3}
\]

is an isomorphism, with common basis
\(\tau^{2n},\tau^{2n+1}\).  The only cokernel occurs in degree zero and
detects the missing linear term \(\tau\).  Since
\(\rho_E=-18\tau^2+\cdots\) in the standard cubic frame, its first nonzero
initial lies in conductor degree one and has class zero in (6.1).

### 6.2 Both branches of the \(2+2\) connector

For a node with normalized branch parameters \(\tau_+,\tau_-\),

\[
\begin{aligned}
 R_{\rm node}
 &=
 \{(f_+,f_-)\in k[[\tau_+]]\oplus k[[\tau_-]]:
       f_+(0)=f_-(0)\},\\
 \widetilde R_{\rm node}
 &=k[[\tau_+]]\oplus k[[\tau_-]],\\
 \mathfrak c
 &=(\tau_+)\oplus(\tau_-).
\end{aligned}                                      \tag{6.4}
\]

In every positive conductor degree the map to the normalization is the
identity on the two branch coefficients.  In degree zero its cokernel is
the value-difference map

\[
 (a_+,a_-)\longmapsto a_+-a_- .                    \tag{6.5}
\]

The residue vanishes at both connector endpoints, so (6.5) is zero.  Its
two first nonzero initials are independent elements of

\[
 (\tau_+)/(\tau_+^2)\ \oplus\
 (\tau_-)/(\tau_-^2),                              \tag{6.6}
\]

and every such pair already descends from
\(\mathfrak c/\mathfrak c^2\).  Therefore

\[
 \boxed{\operatorname{Obs}_{\mathrm{pair}}(\Psi,L)=0}
 \tag{6.7}
\]

for the cusp and for the actual paired connector, before any choice of
sheet labels.

### 6.3 Exact quartic coefficients

The quartic spectator model makes all three initials explicit.  Put

\[
 f(X)=X^4-X^3+uX-v,\qquad
 v=T^4-T^3+uT,
\]

and let \(\Delta(u,v)=\operatorname{disc}_X(f)\).  With

\[
 r=4T^3-3T^2+u
\]

one has

\[
 \Delta(T,u)=r^2\ell(T,u),
\qquad
 \ell|_{r=0}
 =
 -9T^2(2T-1)^2(8T^2-4T-1).                        \tag{6.8}
\]

Thus \(\rho=2\ell\).  The full quartic discriminant frame gives

\[
 \operatorname{in}_{T=0}(\rho)=18T^2.              \tag{6.9}
\]

This differs from the standard cubic-frame coefficient \(-18\) by the
unit relating the two target branch equations.  At

\[
 T_\pm=\frac{1\pm\sqrt3}{4}
\]

the two connector branches give

\[
\begin{aligned}
 \operatorname{in}_{T_+}(\rho)
   &=-\frac{9\sqrt3}{2}(T-T_+),\\
 \operatorname{in}_{T_-}(\rho)
   &= \frac{9\sqrt3}{2}(T-T_-).
\end{aligned}                                      \tag{6.10}
\]

The opposite displayed scalars do not define a mismatch.  Replacing the
two branch parameters independently by
\(\tau_+=\alpha\tau'_+\) and
\(\tau_-=\beta\tau'_-\) rescales them independently.  The conductor
identifies endpoint values, not tangent parameters.  Hence an expression
such as \(c_+-c_-\) is not intrinsic without an additional jet-transport
isomorphism.

### 6.4 Sign, character, and monodromy

The boundary transition

\[
 r_i=u_{ij}r_j
\]

gives

\[
 \rho_i=u_{ij}^{-2}\rho_j,\qquad
 \rho_i^2=u_{ij}^{-4}\rho_j^2.                     \tag{6.11}
\]

Thus the cotangent volume frames select the inverse-square root seen only
up to sign by the scalar Hessian determinant.  This character descent is
compatible; it does not identify the two independent tangent lines in
(6.6).  The primitive inverse character of \(z_i\) likewise concerns
boundary-chart overlap, not node-branch jet transport.

There are 24 ordered nondegenerate cusp braid pairs and three perfect
matchings at the connector.  All 72 labelled packets generate \(S_4\).
Changing these sheet labels does not change (6.4), and therefore every
monodromy-compatible packet has the same zero class (6.7).  Monodromy
selects which sheets collide; it supplies no missing identification of
\((\tau_+)/(\tau_+^2)\) with
\((\tau_-)/(\tau_-^2)\).

## 7. Consequence

The four requested local steps now have exact answers.

1. The completed cusp and connector comparison makes the candidate
   cokernel precise.
2. The inverse-square residue character descends, while the Hessian square
   alone retains a sign ambiguity.
3. On the cotangent branch, the oriented initial Hessian square root divided
   by the normal generator is exactly \(\rho_E\).
4. For the actual connector pairing and every compatible \(S_4\) labelling,
   the conductor class is zero.

This is a negative theorem about the proposed obstruction, not a quartic
exclusion.  A forced nonzero mismatch cannot come from the ordinary
associated-graded conductor.  A viable refinement must add a canonical
transport between endpoint jets.  The two live sources already isolated by
the plane programme are:

- the two-generated degree-zero algebra \(k[x,y]\), which may select a
  global parameter or a constrained jet frame; and
- a dualizing/residue pairing whose sign rule is proved to agree with the
  cotangent square root and the global meridian relation.

Neither structure is contained in (6.1).  Merely expanding the quartic
coefficient search cannot repair this missing datum.

## 8. Remaining flagship gate

The next sharply posed statement is a **global jet-transport theorem**:
construct, from the factorial degree-zero algebra and the marked
degree-four cover, an intrinsic map

\[
 (\tau_+)/(\tau_+^2)\otimes\mathcal L_{p_+}
 \longrightarrow
 (\tau_-)/(\tau_-^2)\otimes\mathcal L_{p_-},        \tag{8.1}
\]

and prove its compatibility with both the \(S_4\) meridian packet and the
oriented cotangent Hessian root.  Only after (8.1) exists does a scalar
initial mismatch become well-defined.  The exact coefficients (6.9)--(6.10)
then provide the first regression case: depending on the proven transport,
they either force nonvanishing or select the anti-diagonal structured
packet.

## 9. Reproduction

Run

```bash
.venv/bin/python scripts/verify_jc2_hc4_isotropic_boundary_bridge.py
```

The checker verifies (1.3), (2.3), (4.4), (5.2), the cusp specialization
(5.3), the quartic discriminant factorization (6.8), both connector
initials (6.10), the cusp and node associated-graded conductor maps, and
all 72 monodromy-compatible labellings.  The finite-normalization,
threefold conductor, and odd-square statements used in Sections 3 and 5
remain canonically sourced in
[`plane-jc/FINITE_NORMALIZATION_PROGRAM.md`](plane-jc/FINITE_NORMALIZATION_PROGRAM.md)
and
[`plane-jc/JC2_GLOBAL_COX_PACKET_ATTACK.md`](plane-jc/JC2_GLOBAL_COX_PACKET_ATTACK.md);
they are not re-proved by this small checker.
