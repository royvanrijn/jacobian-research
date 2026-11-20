# The third syzygy derivation and two polynomial restart obstructions

All fields have characteristic zero. Work in the exact old coordinate algebra
\[
 A_0=K[u,w,p_1,p_2,p_3,A,B,E,D,C,R],\qquad L=K(u,w).
\]
Use the components and pre-\(S\) derivations \(\delta,\epsilon\) of
[the complete scalar criterion](COMMUTING_SYZYGY_SLICES_AND_SCALAR_OBSTRUCTION.md).
The results below concern this specified algebra and two restart
architectures, not arbitrary parallel constructions.

## The third derivation

Put \(T=up_2-2wp_1\), \(\omega=w^2dp_1-2uwdp_2+u^2dp_3\), and
\[
 \gamma=9uw\partial_D-36w^2\partial_C
                     -9(6wp_2-4up_3)\partial_R .
\]
It is locally nilpotent and commutes with \(\delta,\epsilon\), as direct
evaluation on generators proves; in particular
\(\delta(6wp_2-4up_3)=0\).
It is defined on the OLD algebra, not by an assumed stable table on later
coordinates.

**Lemma ThirdFibreCoefficient.** For \(\delta Q=\epsilon Q=0\), the fibre
part of the third Jacobian row of a coordinate with component \(Q\) is
\(\gamma(Q)\omega\).

Indeed the relevant second-power fibre rows are
\[
\begin{array}{c|ccc}
 &p_1&p_2&p_3\\\hline
D&9uw^3&-3u^2w^2&-u^3w\\
C&-36w^4&48uw^3&-20u^2w^2\\
R&-9w^2a&3uwa+3w^2b&u^2a-2uwb
\end{array}
\]
where \(a=6wp_2-4up_3,\ b=-12up_2+24wp_1\); all earlier rows vanish
in these columns. Their product with \((u^2,2uw,3w^2)^t\) is zero,
and with \((2uw,w^2,0)^t\) gives the \(\epsilon\)-coefficients.
Thus \(\epsilon Q=0\) puts the resulting row in the span of \(\omega\);
its first entry divided by \(w^2\) is \(\gamma Q\).

For the seven prior joint generators, direct calculation gives
\[
 \gamma\kappa_C=-108w^2/5,\qquad
 \gamma\kappa_R=108wI_2/5,
\]
and \(\gamma I_2=\gamma I_3=\gamma I_A=\gamma I_B=\gamma I_E=0\).
Here \(\kappa_C=\pi(C)+(8w/5u)\pi(D)\) and
\(\kappa_R=\pi(R)-f\pi(D)/(15u^2w^4)\).

**Theorem TripleLocalizedKernel.**
\[
 \ker\delta\cap\ker\epsilon\cap\ker\gamma
   =L[I_2,I_3,I_A,I_B,I_E,\kappa_R+(I_2/w)\kappa_C].
\]
These six generators are algebraically independent. Apply the slice
\(-5\kappa_C/(108w^2)\) to the prior seven-generator polynomial common
kernel. This is a complete LOCALIZED result, without automatic denominator
clearing or ordinary cubic degree.

## Complete ordinary cubic image

Let \(V\subset A_0\) be the homogeneous ordinary cubics killed by
\(\delta,\epsilon\).

**Theorem CubicGammaImage.** Over every characteristic-zero field,
\[
 \dim V=21,\qquad
 \gamma(V)=w^2\operatorname{span}_K\{u^2,uw,T\},\qquad
 \dim(V\cap\ker\gamma)=18.
\]
Every member of the last cubic kernel ignores \(D,C,R\) and is affine-linear
in \(A,B,E\).

For an exact finite proof, enumerate all \(\binom{13}{3}=286\) cubic
monomials. Exact column elimination of their \((\delta,\epsilon)\)-images
has rank 265 and gives a complete 21-element kernel basis. Its
\(\gamma\)-images lie in the displayed span and have rank three.
Elimination for \((\delta,\epsilon,\gamma)\) has rank 268 and gives a
complete 18-element basis with the stated support restrictions.
The checker reconstructs these rational matrices, ranks and ALL basis
vectors. Scalar extension preserves the exact ranks in characteristic zero.
This classifies one finite-degree space, not an infinite closure.

Two helpers in \(V\) are
\[
\begin{aligned}
Y_0={}&5u^2C+8uwD-60up_1B-6up_2A+6up_2E\\
 &-20wp_1A-12wp_1E-120p_1^2p_3+105p_1p_2^2,\\
Y_1={}&5uwC+46up_2B+\tfrac{74}{3}up_3A+\tfrac83up_3E\\
 &+8w^2D-152wp_1B-53wp_2A-4wp_2E+10p_1p_2p_3 .
\end{aligned}
\]
They satisfy
\[
 \gamma Y_0=-108u^2w^2,\quad
 \gamma Y_1=-108uw^3,\quad
 \gamma P_S=-108w^2T.
\]
All cubic helpers with either prescribed image are
\(Y_i+(V\cap\ker\gamma)\); the full 18-dimensional freedom is retained.

## Two restart families

Introduce a NEW coordinate \(y\) with component \(Y\), and consider
\[
 Q=aS-Ty+Z,\qquad Z\in A_0,
\]
where either \(a=u^2,\gamma Y=-108u^2w^2\), or
\(a=uw,\gamma Y=-108uw^3\).
Here \(S,y\) are coordinates, not cubic components, so the first two terms
are cubic.

Let \(\eta\) be the old derivation defined by
\[
 \eta z=\sum_{j\ {\rm moving}}(\partial_jP_z)(\partial_uP_j),
 \qquad \eta u=\eta w=0.
\]
It records the \(du\)-coefficient in second-power rows.
For every pre-\(S\) component \(P_z\), the finite component table gives
\(\delta z=(P_z)_w\), \(\delta P_z=0\), and \(\epsilon P_z=0\).
Its second-power row is
\[
 e_zJ^2=(\delta z)s+\eta z\,du+(\delta P_z)dw+F_z,
 \qquad s=dwJ ,
\]
where \(F_z\) is exactly the three-column fibre row displayed earlier;
all its other new-coordinate columns are zero. This identity follows
by expanding \(dP_zJ\): the \(dw\) coefficient is
\(\sum_j(P_z)_j(P_j)_w=\delta P_z\).
In particular that coefficient vanishes individually.

For a valid component \(H\in V\), sum these rows against \(H_z\), and
include \(H_w\,dwJ^2=H_wr\). Its \(s\) coefficient is
\(\sum_jH_j(P_j)_w=\delta H=0\), its \(dw\) coefficient is zero by
the preceding component identities, and \(\epsilon H=0\) puts its
remaining fibre row in \(\gamma(H)\omega\). Thus the COMPLETE third row is
\[
 dH\,J^2=H_wr+\eta(H)du+\gamma(H)\omega .
\]
The checker verifies all ten component invariance identities, all
\(9\cdot58=522\) pre-\(S\) second-row entries, and all
\(21\cdot58=1218\) entries for the full cubic joint-kernel basis.
The derivation proves the identity for every ordinary joint invariant;
the finite full-basis comparisons additionally cover the complete
cubic space used below.

The new \(\omega\) part cancels. The old gradients of \(-Ty\) give
\(2wy(u^2r)-uy(2uwr)=0\), while the \(u,w\) third rows vanish.
Consequently there are NO new \(S,y\) variables in the residual:
\[
 (\alpha+\delta Z)r+(\beta+\epsilon Z)du,\qquad
 \alpha=a(P_S)_w-TY_w,\quad
 \beta=a\eta(P_S)-T\eta(Y).
\]
Actual-base independence of \(r,du\) makes cancellation equivalent to
\(\delta Z=-\alpha,\epsilon Z=-\beta\).

**Theorem BothSimpleRestartsFail.** For EACH value of \(a\) above, EVERY
ordinary homogeneous cubic helper with the stated \(\gamma\)-image, and
EVERY old ordinary polynomial \(Z\), the candidate fails to preserve
fourth Jacobian power zero.

For \(a=u^2,Y_0\), expansion gives
\[
 \beta\equiv48u^4w^3B+63u^3w^3p_2^2\pmod{(w^4,p_3)}.
\]
Every \(\epsilon Z\) lies in \((w^4,p_3)\), by the coefficients of
\(\epsilon\). Also
\([u^3w^3p_2^2]\,T\eta K=0\) on ALL 18 cubic triple-kernel basis
vectors. The coefficient 63 therefore persists under every allowed helper
variation and cannot be cancelled. The complete projected image has
rank one, versus augmented rank two.

For the second family use the quartic coefficient functional
\[
\begin{aligned}
\Lambda(F)={}&4[wp_1p_2E]F+48[up_2^2E]F+33[up_1p_3E]F\\
 &-144[uwE^2]F-36[uwAE]F+42[u^2BE]F.
\end{aligned}
\]
Expansion on all 286 cubics proves \(\Lambda(\delta Z)=0\) for every
cubic \(Z\). Since \(\delta\) increases ordinary degree by one, it holds
for EVERY ordinary polynomial \(Z\).
Expansion on the complete 18-element triple kernel gives
\(\Lambda(TK_w)=0\), whereas
\(\Lambda(uw(P_S)_w-T(Y_1)_w)=1312\ne0\).
Thus every helper variation still contradicts \(\delta Z=-\alpha\).

In both families \(\delta\beta=\epsilon\alpha\).
An old correction DOES exist in the localized algebra
\(L[p_1,\ldots,R]\): the prior two slices turn these equations into
a closed polynomial one-form in two variables, which can be integrated
successively in characteristic zero.
The obstruction is ordinary polynomiality at parameter divisors, not
rational compatibility.

The theorem leaves general mixed helper blocks, arbitrary quadratic
coefficient syzygies, nonlinear new-coordinate terms and iterated
triple-invariant constructions open. It is not a universal parallel-block
impossibility theorem. No new inverse-degree computation is used.

## Exact replay

The standard-library checker pins the prior slice proof/checker and full
literal degree-207 data. It reconstructs the complete cubic kernels, image
rank, third derivation, triple normal form, both helper families, all
obstruction functionals and compatibility identities. It also checks the
all1218 joint-basis row entries, all522 pre-S second-row entries, and the
full59-coordinate helper block and60-coordinate residual for both
representatives, without computing inverse degrees.

Reproduce with
python3 research/jacobian/subtree-shuffling/verify_third_covector_restart_obstructions.py --repo-root . --output /tmp/third_covector_restart_review.json.
The source ledger pins every imported local module and literal input.
The receipt includes the complete21/18 cubic bases as exact rational vectors.

20 November 2026. Internal written proof and exact finite replay only.
