# The nonreduced Hessian--Schur module

## Status

This note proves `HC4NHM1`, a structural reduction for the nonaligned rank-three
degree-five packet in the direct `HC4` filtration.  It does **not** construct
or exclude an `HC(4)` counterexample.

The main new point is a normalization degree gate.  For a ternary quintic it
splits every nonzero Schur section into one of two types:

1. a quartic-denominator packet
   \(\det C=P^2\ell\), with \(\deg P=4\) and \(\ell\) linear; or
2. a packet supported on an explicitly large corank-two defect divisor of a
   repeated Hessian component.

Thus repeated smooth cubics cannot occur in the defect-free packet, a smooth
conic can occur there only with multiplicity four, and a smooth quartic can
occur there only with multiplicity two and a trivial degree-zero kernel
class.  The gradient condition is imposed only after this module-theoretic
split.

Replay the exact algebraic calibrations with

```bash
.venv/bin/python scripts/verify_hc4_nonreduced_hessian_schur_module.py
```

The squarefree starting point is
[`HC4_DIRECT_SQUAREFREE_HESSIAN_OBSTRUCTION.md`](HC4_DIRECT_SQUAREFREE_HESSIAN_OBSTRUCTION.md),
and the existing repeated-linear analysis is
[`HC4_DIRECT_DOUBLE_LINEAR_HESSIAN_GATE.md`](HC4_DIRECT_DOUBLE_LINEAR_HESSIAN_GATE.md).

Throughout, `K` is a field of characteristic zero,
\(S=K[x,y,z]\), and

\[
C=C^{\mathsf T}\in\operatorname{Mat}_3(S_r),\qquad
\Delta=\det C\ne0,\qquad D=V(\Delta)\subset\mathbf P^2.
\tag{0.1}
\]

For the minimal quintic packet,

\[
r=3,\qquad C=\operatorname{Hess}(h_5),\qquad
d=\nabla s_3\in S_2^3,
\tag{0.2}
\]

and the first open face is

\[
\Delta\mid d^{\mathsf T}\operatorname{adj}(C)d.
\tag{0.3}
\]

## 1. The symmetric matrix factorization

Let \(i:D\hookrightarrow\mathbf P^2\) and define \(\mathcal L\) by

\[
0\longrightarrow \mathcal O_{\mathbf P^2}(-r)^3
 \xrightarrow{C}\mathcal O_{\mathbf P^2}^3
 \longrightarrow i_*\mathcal L\longrightarrow0.
\tag{1.1}
\]

The identities

\[
C\operatorname{adj}(C)=\operatorname{adj}(C)C=\Delta I_3
\tag{1.2}
\]

are a symmetric graded matrix factorization of \(\Delta\).  Dualizing
(1.1) against \(\mathcal O_{\mathbf P^2}(-r)\) gives the same presentation.
Grothendieck duality for the degree-\(3r\) Cartier divisor gives

\[
\mathcal E xt^1_{\mathbf P^2}
 (i_*\mathcal L,\mathcal O_{\mathbf P^2}(-r))
 \simeq i_*\mathcal H om_D
 (\mathcal L,\mathcal O_D(2r)).
\tag{1.3}
\]

Consequently the symmetric presentation supplies a self-duality

\[
\mathcal L\simeq\mathcal L^\vee(2r).
\tag{1.4}
\]

On the locus where \(\mathcal L\) is generically a rank-one locally free
module, this is equivalently

\[
\mathcal L^{[2]}\simeq\mathcal O_D(2r).
\tag{1.5}
\]

There is an important qualification in (1.5): it applies to the
generic-corank-one Smith stratum.  If \(C\bmod q\) has rank at most one along
an entire component, the cokernel remains a self-dual maximal
Cohen--Macaulay module, but it need not be a generically invertible module.
The radial sextic in Section 6 is of this lower-Smith type.

The pairing in coordinates is exactly

\[
B_C(\bar u,\bar v)
=u^{\mathsf T}\operatorname{adj}(C)v\pmod\Delta.
\tag{1.6}
\]

It descends to the cokernel because replacing \(u\) by \(u+Ca\) changes
the numerator by \(\Delta a^{\mathsf T}v\).

A vector \(d\in S_{r-1}^3\) defines

\[
\sigma_d\in H^0(D,\mathcal L(r-1)).
\tag{1.7}
\]

Indeed, after twisting (1.1) by \(r-1\), both
\(H^0(\mathbf P^2,\mathcal O(-1))\) and
\(H^1(\mathbf P^2,\mathcal O(-1))\) vanish, so

\[
H^0(D,\mathcal L(r-1))\simeq S_{r-1}^3.
\tag{1.8}
\]

Under (1.6), condition (0.3) is precisely

\[
B_C(\sigma_d,\sigma_d)=0.
\tag{1.9}
\]

On the generically invertible locus this may be called
\(\sigma_d^{[2]}=0\).  The pairing formulation (1.9) is the one that also
makes sense on lower Smith strata.

If \(D\) is reduced and \(C\) has generic corank one on every component,
then a rank-one torsion-free module embeds in the product of the component
function fields.  The square of a nonzero element cannot vanish there.
Thus (1.9) gives \(\sigma_d=0\), and (1.8) gives \(d=0\).  This recovers
the squarefree theorem without a coefficient calculation.

## 2. Local nilpotence and the minimal denominator

Let \(q\) be an irreducible factor of \(\Delta\), with
\(m=v_q(\Delta)\), and assume that \(C\bmod q\) has generic rank two.
Over the DVR at the generic point of \(q=0\), symmetric congruence gives

\[
C\sim\operatorname{diag}(u_1,u_2,u_3q^m),
\qquad u_i\in R^\times.
\tag{2.1}
\]

Writing \(d=(d_1,d_2,d_3)\) in this frame, (0.3) says

\[
q^m\mid
u_2u_3q^m d_1^2+u_1u_3q^m d_2^2+u_1u_2d_3^2.
\tag{2.2}
\]

Hence

\[
v_q(d_3)\ge\left\lceil\frac m2\right\rceil,
\qquad
v_q(C^{-1}d)\ge-\left\lfloor\frac m2\right\rfloor.
\tag{2.3}
\]

This is the local nilradical filtration: the kernel coordinate of the
section lies at least halfway into the thickening.

Assume generic corank one along every component of \(D\).  Let \(P\) be
the componentwise-minimal homogeneous denominator for \(C^{-1}d\):

\[
e:=P C^{-1}d\in S^3,
\qquad
q\mid P\Longrightarrow q\nmid\gcd(e_1,e_2,e_3).
\tag{2.4}
\]

The squarefree components contribute no pole.  Equation (2.3) gives

\[
P^2\mid\Delta.
\tag{2.5}
\]

Homogeneity is especially restrictive:

\[
\deg(C^{-1}d)=-1,qquad
\rho:=\deg P,qquad
e\in S_{\rho-1}^3.
\tag{2.6}
\]

If \(P=1\), then \(C^{-1}d\) would be a polynomial vector of degree
\(-1\), so \(d=0\).  Every nonzero Schur section therefore has
\(\rho\ge1\).

## 3. Normalization and the corank-two defect divisor

Fix an irreducible component \(Q=V(q)\) of degree \(k\) with \(q\mid P\),
and let

\[
\nu:\widetilde Q\longrightarrow Q
\tag{3.1}
\]

be its normalization.  Pull back \(C|_Q\) and let \(\mathcal K_Q\) be the
saturated kernel line in \(\mathcal O_{\widetilde Q}^3\):

\[
\mathcal K_Q
=\ker\left(\mathcal O_{\widetilde Q}^3
 \xrightarrow{\nu^*C}
 \mathcal O_{\widetilde Q}(r)^3\right)^{\rm sat}.
\tag{3.2}
\]

Put \(\mathcal F_Q=\mathcal O_{\widetilde Q}^3/\mathcal K_Q\).
Symmetry induces a generically invertible rank-two map

\[
\bar C:\mathcal F_Q\longrightarrow\mathcal F_Q^\vee(r).
\tag{3.3}
\]

Let \(B_Q\) be the zero divisor of \(\det\bar C\).  It is the normalized,
scheme-theoretic corank-two defect: away from singularities of \(Q\), its
support is exactly the set where all \(2\)-by-\(2\) minors of \(C\) vanish
on \(Q\).

Since \(\det\mathcal F_Q=\mathcal K_Q^{-1}\), taking determinants in
(3.3) gives

\[
\boxed{
\mathcal K_Q^{\otimes2}\simeq
\mathcal O_{\widetilde Q}(B_Q-2r).
}
\tag{3.4}
\]

Here \(\deg\mathcal O_{\widetilde Q}(1)=k\), so

\[
\deg\mathcal K_Q=\frac{\deg B_Q-2rk}{2}.
\tag{3.5}
\]

The residue of the vector \(e\) from (2.4) is nonzero on \(Q\) and
satisfies \(Ce=0\bmod q\).  It is therefore a nonzero section

\[
e|_{\widetilde Q}\in
H^0\left(\widetilde Q,
\mathcal K_Q(\rho-1)\right).
\tag{3.6}
\]

Every line bundle with a nonzero section has nonnegative degree.  Combining
(3.5) and (3.6) proves the reusable normalization gate

\[
\boxed{
\deg B_Q\ge 2k(r+1-\rho).
}
\tag{3.7}
\]

When the right side is zero and \(B_Q=0\), the section is nowhere vanishing
and

\[
\mathcal K_Q(\rho-1)\simeq\mathcal O_{\widetilde Q}.
\tag{3.8}
\]

Formula (3.4) packages both the conductor correction from normalization and
the finite corank-two failures of the generic Smith form.  It is the global
counterpart of (2.3).

## 4. The ternary-quintic dichotomy

Now take \(r=3\).  Since \(\deg\Delta=9\), (2.5) gives

\[
1\le\rho\le4.
\tag{4.1}
\]

For every component \(Q\mid P\), (3.7) becomes

\[
\boxed{
\deg B_Q\ge2\deg(Q)(4-\rho).
}
\tag{4.2}
\]

Consequently, if even one essential denominator component has
\(B_Q=0\), then \(\rho=4\).  Equation (2.5) and the degree-nine determinant
then force

\[
\boxed{
\Delta=P^2\ell,qquad \deg P=4,qquad\deg\ell=1.
}
\tag{4.3}
\]

Thus every nonzero generic-corank-one quintic Schur section lies in one of
the following two packets.

> **Quintic normalization dichotomy.** Either
>
> 1. \(\rho=4\) and \(\Delta=P^2\ell\); or
> 2. \(1\le\rho\le3\), and every irreducible component on which the
>    denominator is essential carries the positive defect required by
>    (4.2).

For one smooth repeated factor and a squarefree coprime residual factor,
the weakest consequences (using the largest locally permitted pole) are:

| repeated component | possible multiplicity in \(\deg\Delta=9\) | defect-free consequence |
|---|---:|---|
| line | \(2\le m\le9\) | requires \(m\ge8\); for \(m=2,3\), \(\deg B\ge6\); for \(m=4,5\), \(\deg B\ge4\); for \(m=6,7\), \(\deg B\ge2\) |
| smooth conic | \(m=2,3,4\) | \(m=2,3\) require \(\deg B\ge8\); the clean case is only \(m=4\) |
| smooth cubic | \(m=2,3\) | no clean case; \(\deg B\ge6\) |
| smooth quartic | \(m=2\) | possible only in the quartic-denominator packet |

In the clean \(\rho=4\) case, (3.8) says

\[
\mathcal K_Q(3)\simeq\mathcal O_{\widetilde Q}
\tag{4.4}
\]

on every essential component.  For a smooth conic this degree-zero class is
automatically trivial because its normalization is \(\mathbf P^1\).  For a
smooth quartic it is a genuine restriction: (3.4) makes
\(\mathcal K_Q(3)\) a two-torsion line bundle, and (4.4) selects the trivial
two-torsion class.

The table does not classify the positive-defect rows.  It replaces the
phrase “corank-two locus” by an exact minimum length that can be computed as
the common divisor of the restricted \(2\)-by-\(2\) minors.

## 5. Exact equations on the quartic-denominator packet

Assume (4.3) and retain the minimal denominator definition (2.4).  Then

\[
e=P C^{-1}d\in S_3^3
\tag{5.1}
\]

and the Schur problem is equivalent to the existence of a linear form
\(a\in S_1\) satisfying

\[
\boxed{
\det C=P^2\ell,qquad
Ce=Pd,qquad
d^{\mathsf T}e=Pa.
}
\tag{5.2}
\]

Indeed, multiplying \(Ce=Pd\) by \(\operatorname{adj}(C)\) gives

\[
\operatorname{adj}(C)d=P\ell e,
\tag{5.3}
\]

and then

\[
d^{\mathsf T}\operatorname{adj}(C)d
=P\ell\,d^{\mathsf T}e
=P^2\ell a
=\Delta a.
\tag{5.4}
\]

Conversely, (0.3) and (5.1) give \(P\mid d^{\mathsf T}e\), so (5.2) loses
no point of this packet.

Only now impose Hessian and gradient integrability:

\[
C=\operatorname{Hess}(h_5),qquad d=\nabla s_3.
\tag{5.5}
\]

This order matters.  For fixed \(P\), the equation \(Ce\in P S^3\) is a
kernel/conductor calculation, while the curl equations for \(Ce/P\) are
linear in \(e\).  In characteristic zero a closed homogeneous quadratic
vector is automatically the gradient of the cubic

\[
s_3=\frac13(x,y,z)\cdot d.
\tag{5.6}
\]

Euler's identity also turns (5.2) into useful first-order vector-field
relations.  Writing \(D_e=e\cdot\nabla\), one obtains

\[
4D_eh_5=3Ps_3,qquad D_es_3=Pa,
\tag{5.7}
\]

The curl-free condition for \(Ce/P\) has the matrix form

\[
C\operatorname{Jac}(e)-\operatorname{Jac}(e)^{\mathsf T}C
=\nabla s_3(\nabla P)^{\mathsf T}
 -\nabla P(\nabla s_3)^{\mathsf T}.
\tag{5.8}
\]

Multiplying (5.8) by the Euler vector, or directly differentiating the
first identity in (5.7), gives

\[
(\operatorname{Jac}e)^{\mathsf T}\nabla h_5
=\frac34s_3\nabla P-\frac14P\nabla s_3.
\tag{5.9}
\]

Equations (5.2), (5.7), and (5.8) are a smaller invariant target than the
original universal coefficient divisibility.  The factor partitions of the
quartic \(P\), rather than all ternary quintics at once, give the natural
finite chart decomposition.

### 5.1 The clean smooth-quartic subpacket

Suppose \(P=q\) is a smooth irreducible quartic and \(B_Q=0\).  Then
\(e|_Q\) is a nowhere-vanishing generator of \(\mathcal K_Q(3)\).  Both
\(\operatorname{adj}(C)|_Q\) and \(e e^{\mathsf T}|_Q\) are nonzero
symmetric rank-one matrices with the same image.  After base change to an
algebraic closure, their ratio is a global regular function on the projective
integral curve \(Q\), hence a nonzero scalar \(\lambda\).  Therefore

\[
\operatorname{adj}(C)=\lambda e e^{\mathsf T}+qA,
\qquad A=A^{\mathsf T}\in\operatorname{Mat}_3(S_2).
\tag{5.10}
\]

Using \(Ce=q d\), \(d^{\mathsf T}e=qa\), and
\(C\operatorname{adj}(C)=q^2\ell I\) gives the lower-degree paired
factorization

\[
\boxed{
\begin{aligned}
CA+\lambda d e^{\mathsf T}&=q\ell I,\\
AC+\lambda e d^{\mathsf T}&=q\ell I,\\
Ad&=(\ell-\lambda a)e.
\end{aligned}}
\tag{5.11}
\]

Thus the smooth-quartic row is not a general determinant-factorization
problem: it carries an auxiliary symmetric quadratic matrix \(A\) and a
degree-\((3,2)\) paired factorization.  The subsequent theorem `HC4NHM14`
sets \(\mu=\ell-\lambda a\), proves
\(\det A=q\ell\mu\) and
\(\operatorname{adj}(A)=\mu C+\lambda dd^{\mathsf T}\), and splits the
residual-line boundary into nine simple-line and ten doubled-line gradient
types.  The smooth-conic multiplicity-four row has
\(P=q^2\); it has the same first residue
\(\operatorname{adj}(C)\equiv\lambda ee^{\mathsf T}\bmod q\), but requires
one more normal layer before a quadratic matrix such as \(A\) appears.

## 6. Calibration against the exact Schur atlas

### 6.1 Fermat

For the normalized Fermat form of degree \(r+2\),

\[
C=\operatorname{diag}(x^r,y^r,z^r),qquad
\Delta=x^ry^rz^r.
\tag{6.1}
\]

Along \(x=0\), the kernel line is constant, while the induced rank-two
determinant is \(y^rz^r\).  Hence

\[
\mathcal K_x\simeq\mathcal O_{\mathbf P^1},qquad
\deg B_x=2r.
\tag{6.2}
\]

For the quintic \(r=3\), a single diagonal cubic channel has \(P=x\),
\(\rho=1\), and (4.2) is the equality \(\deg B_x=6\).  For the sextic
\(r=4\), the same calculation gives \(\deg B_x=8\).  Thus the Fermat
Schur spaces survive because of their finite corank-two defect divisors,
not because the normalization degree gate misses them.

### 6.2 Radial

For \(R=x^2+y^2+z^2\) and \(h_6=R^3/30\),

\[
C=\frac R5(RI+4xx^{\mathsf T}),qquad
\det C=\frac{R^6}{25}.
\tag{6.3}
\]

At the generic point of \(R=0\), the Smith valuations are \((1,2,3)\), so
\(C\bmod R\) has rank zero.  The unique radial Schur section therefore
belongs to the lower-Smith module stratum excluded from (2.1), not to the
generic-corank-one normalization theorem.  This is why a separate lower
Smith analysis is essential.

## 7. Relation with the ternary-quintic Hessian map

Ciliberto and Ottaviani prove birationality of the ternary Hessian map for
all degrees \(d\ge4\), \(d\ne5\), and explicitly state that they expect the
exception to be removable.  In their odd-degree argument, \(d=2k+1\), the
missing graph-boundary exclusion is

\[
(\ell^{2k+1},q^{3(k-1)}\ell^3).
\tag{7.1}
\]

Their degeneration theorem forces only \(\ell^{d-3}\) in a limiting
Hessian.  For \(d\ge7\), this has order strictly larger than three and
excludes (7.1); for \(d=5\), it gives only \(\ell^2\), so the boundary pair

\[
(\ell^5,q^3\ell^3)
\tag{7.2}
\]

is not excluded.  This identifies the exact reason degree five falls out of
their proof, rather than merely recording a numerical exception.  See
[Ciliberto--Ottaviani, Theorem 2.23 and Corollary 3.6](https://ems.press/content/serial-article-files/53074?nt=1).

Beorchia's graph-of-the-polar-map approach identifies the Hessian as the
ramification divisor and gives a conditional recovery result when the ruled
surface joining the two lifted ramification curves is a product.  The
revised paper explicitly leaves the non-product ruled-surface case open; it
does not independently settle generic injectivity in degree five.  See
[Beorchia, Sections 3--4](https://arxiv.org/pdf/2406.05423).

The two boundary problems are related but not identical.  The Hessian-map
problem studies fibers of \(h_5\mapsto\det\operatorname{Hess}(h_5)\); the
`HC4` Schur problem adds the self-dual section \(\sigma_d\).  A useful common
target is therefore the pair

\[
(D,\mathcal L,\sigma_d)
\tag{7.3}
\]

over the nonreduced Hessian locus, not just the divisor \(D\).

## 8. Next exact calculations

The module reduction suggests the following order.

1. Compute \(P\) and the normalized defect divisors \(B_Q\) from the gcd of
   the restricted \(2\)-by-\(2\) minors.  The bounds (4.2) eliminate every
   insufficient defect before any Schur quotient is introduced.
2. On the clean packet, split the quartic \(P\) by partitions
   \(4\), \(3+1\), \(2+2\), \(2+1+1\), and \(1+1+1+1\), and solve the
   kernel equations \(Ce=Pd\).
3. Impose the curl-free condition on \(d=Ce/P\), then the single scalar
   equation \(d^{\mathsf T}e=Pa\).
4. Treat positive-defect rows and generic corank-two/three Smith partitions
   separately.  Fermat and radial show that these rows are necessary, but
   their defect lengths and Smith partitions make them much smaller than
   the full nonsquarefree discriminant.

For the first positive-defect continuation, the septuple-line pole/defect
ladder and the exclusion of its extremal quadratic-pencil kernel are proved
in
[`HC4_DIRECT_SEPTUPLE_LINEAR_HESSIAN_GATE.md`](HC4_DIRECT_SEPTUPLE_LINEAR_HESSIAN_GATE.md).
The subsequent defect-free cubic-kernel calculation closes the octuple and
nonuple generic-corank-one line packets, equivalently the quartic-denominator
partition \(P=x^4\), in
[`HC4_OCTUPLE_NONUPLE_LINEAR_HESSIAN_GATE.md`](HC4_OCTUPLE_NONUPLE_LINEAR_HESSIAN_GATE.md).
The first two-component partition \(P=x^3y\), with residual line equal to
\(x\), is classified in
[`HC4_TWO_LINE_QUARTIC_DENOMINATOR_PACKET.md`](HC4_TWO_LINE_QUARTIC_DENOMINATOR_PACKET.md).
Its nonzero module section is then excluded from four-variable prolongation
at the next determinant face in
[`HC4_TWO_LINE_QUARTIC_DENOMINATOR_PROLONGATION.md`](HC4_TWO_LINE_QUARTIC_DENOMINATOR_PROLONGATION.md).
The other two residual-line incidences in the same \(3+1\) partition are
empty already as ternary Hessian boundaries by
[`HC4_REMAINING_THREE_ONE_QUARTIC_DENOMINATOR_GATE.md`](HC4_REMAINING_THREE_ONE_QUARTIC_DENOMINATOR_GATE.md).
The \(2+2\) and \(2+1+1\) partitions are empty by
[`HC4_TWO_TWO_QUARTIC_DENOMINATOR_GATE.md`](HC4_TWO_TWO_QUARTIC_DENOMINATOR_GATE.md)
and
[`HC4_TWO_ONE_ONE_QUARTIC_DENOMINATOR_GATE.md`](HC4_TWO_ONE_ONE_QUARTIC_DENOMINATOR_GATE.md),
respectively. For the last squarefree partition, pole order one forces four
constant kernel directions and the double-line determinant equations become
linear multiple-polar conditions; see
[`HC4_SQUAREFREE_QUARTIC_DENOMINATOR_FRONTEND.md`](HC4_SQUAREFREE_QUARTIC_DENOMINATOR_FRONTEND.md).
The exact concurrence, general-position, and tangent-fourth closures
`HC4NHM10--12` eliminate all forty-eight resulting flag rows.  Therefore the
complete split-linear clean quartic-denominator packet is empty.  Clean
denominators with irreducible conic, cubic, or quartic components remain,
followed by the positive-defect and lower-Smith strata in item 4.  The first
smooth-cubic continuation `HC4NHM19` gives an orthogonal normal form on the
elliptic component.  The rank-two quotient splits into isotropic line bundles
of degrees `(0,9)`, `(2,7)`, `(3,6)`, or `(4,5)`, and Hessian integrability
excludes `(0,9)`.  Thus the cubic-plus-line packet is reduced invariantly to
the last three degree types; see
[`HC4_SMOOTH_CUBIC_ORTHOGONAL_NORMAL_FORM.md`](HC4_SMOOTH_CUBIC_ORTHOGONAL_NORMAL_FORM.md).
The first
smooth-quartic continuation `HC4NHM14` converts the irreducible-quartic row
to a reciprocal quadratic factorization and a finite residual-line atlas;
see
[`HC4_SMOOTH_QUARTIC_RECIPROCAL_FRONTEND.md`](HC4_SMOOTH_QUARTIC_RECIPROCAL_FRONTEND.md).
Its first exact elimination, `HC4NHM16`, excludes the generic point of the
basepoint-free squarefree-line type $d_0=(x^2,y^2,0)$, while leaving its
exceptional parameter locus and the other residual-line types open; see
[`HC4_SMOOTH_QUARTIC_SQUAREFREE_LINE_GENERIC_GATE.md`](HC4_SMOOTH_QUARTIC_SQUAREFREE_LINE_GENERIC_GATE.md).
The follow-up `HC4NHM17` enters the first visible exceptional divisor and
closes nine generic or algebraic slices; see
[`HC4_SMOOTH_QUARTIC_SQUAREFREE_LINE_EXCEPTIONAL_SLICES.md`](HC4_SMOOTH_QUARTIC_SQUAREFREE_LINE_EXCEPTIONAL_SLICES.md).
The visible pivot in those two calculations is identified invariantly by
`HC4NHM20` as a first polar of the binary resultant
`Res(s^3+t^3,H)`.  Its generic coefficient-space fiber is a smooth conic,
and its squarefree degree-fifteen degeneration locus has only two-line
fibers; see
[`HC4_SMOOTH_QUARTIC_PIVOT_POLAR_GEOMETRY.md`](HC4_SMOOTH_QUARTIC_PIVOT_POLAR_GEOMETRY.md).
The follow-up `HC4NHM22` uses the universal point `[1:3:1]` to parametrize
the generic polar conic.  Ten linear pivots and eight selected reciprocal
coefficients put the sixth power of every active deformation coordinate in
the exact ideal, leaving only the determinant-zero boundary matrix.  Thus
the generic polar-conic component is empty; see
[`HC4_SMOOTH_QUARTIC_POLAR_CONIC_GENERIC_GATE.md`](HC4_SMOOTH_QUARTIC_POLAR_CONIC_GENERIC_GATE.md).
The symmetry continuation `HC4NHM23` proves that the fifteen non-generic
line-fiber slopes have only three Fermat-automorphism normal forms of sizes
`3+6+6`.  It transports the exact `tau=-1` certificates across the complete
size-three orbit `tau^3=-1`, leaving only two genuine line-fiber types; see
[`HC4_SMOOTH_QUARTIC_FERMAT_SYMMETRY_ORBITS.md`](HC4_SMOOTH_QUARTIC_FERMAT_SYMMETRY_ORBITS.md).
The conic continuation `HC4NHM13` excludes the complete conic-divisible top
subrow for the double-conic denominator `P=q^2`; see
[`HC4_SMOOTH_CONIC_DIVISIBLE_TOP_GATE.md`](HC4_SMOOTH_CONIC_DIVISIBLE_TOP_GATE.md).
For the complementary nonzero-restriction row, `HC4NHM15` gives the exact
four-layer harmonic decomposition of the Hessian determinant into binary
covariants. It excludes all decics with at most three distinct roots and the
harmonic-cross-ratio four-root locus before any Schur equation is imposed.
It also closes eight complete arbitrary-cross-ratio partitions and confines
the last row \((3,3,2,2)\) to a finite exceptional cross-ratio locus; see
[`HC4_DOUBLE_CONIC_NORMAL_LAYERS.md`](HC4_DOUBLE_CONIC_NORMAL_LAYERS.md).
The follow-up `HC4NHM18` removes that locus by endpoint normal-layer chains
and a two-equation middle contradiction. Thus every restriction supported
on at most four points is empty before the Schur equations; see
[`HC4_DOUBLE_CONIC_BALANCED_FOUR_ROOT_CLOSURE.md`](HC4_DOUBLE_CONIC_BALANCED_FOUR_ROOT_CLOSURE.md).
For the remaining many-root row, `HC4NHM21` shows that invariant elimination
must first saturate the four normal-layer ideal by the three coefficients of
the nonzero residual line `Phi_2`.  Without that saturation the proposed
discriminant membership is false, even for the squarefree decic
`s^10+t^10`; after saturation, discriminant membership closes only the
squarefree open, while nullcone containment is the correct all-stable target.
See
[`HC4_DOUBLE_CONIC_INVARIANT_SATURATION_GATE.md`](HC4_DOUBLE_CONIC_INVARIANT_SATURATION_GATE.md).

This is the precise self-dual matrix-factorization version of the surviving
minimal rank-three packet.
