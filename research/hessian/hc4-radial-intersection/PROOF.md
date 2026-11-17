# A local intersection test for omitted pointed-Keller sheets

13 November 2026. HC4MYRB1 is the local intersection/canonical gate and
its chain exclusion. HC4MYRB2 is a calibrated radial-resolution calculation
on the already retained tame-node family. Work in characteristic zero,
over an algebraically closed field; use C for the cited surface geometry.
No Tang theorem is used. These are applications of established surface
techniques, without a novelty-priority claim.

## 1. Fix the actual finite normalization and resolve the radial ideal

Let G:A2 -> A2 be Keller with (G1,G2)=(x,y), of generic degree d.
Let N -> A2_target be the finite normalization in the source function
field. By [HC4MYOR2](../hc4-relative-row-orbit/PROOF.md), N is flat of
rank d and contains the actual source as an open subset. Its zero fibre
contains the source origin of length one and omitted points q with
lengths ell_q satisfying sum ell_q=d-1.

Choose a smooth resolution W -> N on which the pullback of the ideal
(G1,G2) is invertible, resolving further to SNC if needed. Near each
omitted q, write the exceptional curves as C1,...,Cs and put

\[
 (G_1,G_2)\mathcal O_W=\mathcal O_W(-\sum_i e_i C_i),
 \qquad Q_{ij}=C_i\cdot C_j.
\]

All e_i are positive integers. Include **all** exceptional curves over
q, including those contracted by the radial pencil. The universal
property of the blowup gives p:W -> Y=Bl_0(A2_target). Its exceptional
line E has normal bundle O_P1(-1). Let f_i be the degree of C_i -> E
when this map is nonconstant, and put f_i=0 for a contracted curve.
The resulting labels use a resolved pencil, not just a smooth surface.

The matrix Q is negative definite: it is an exceptional intersection
matrix over a normal surface point. The precise import is
[Mumford1961, SectionI, p.6](https://www.dam.brown.edu/people/mumford/alg_geom/papers/1961a--TopSingSurface-NumDam.pdf).
The inherited complete-source boundary machinery is
[A2B1](../../foundations/plane-jc/INTRINSIC_A2_BOUNDARY_GATE.md).
It does not permit replacing a complete global boundary by this block
when computing the global canonical coefficients.

## 2. The intersection equation and each point's missing length

Restrict O_W(sum e_i*C_i)=p^*O_Y(E) to C_j. If C_j maps to E with
degree f_j, that line bundle has degree -f_j; if C_j is contracted,
its degree is zero. Hence, directly from degrees of line bundles,

\[
 \boxed{Qe=-f,\qquad e=(-Q)^{-1}f.}\tag{I}
\]

This uses no guessed ramification support. Curves with f_i=0 still
affect the matrix and the integrality of e. At least one f_i is positive.

There is also a pointwise refinement of the previous total radial
degree identity:

\[
 \boxed{\ell_q=e^Tf=f^T(-Q)^{-1}f=-e^TQe.}\tag{L}
\]

Here is why its left side really is local fibre length. Complete the
target at zero, R=C[[u,v]]. The finite flat normalization algebra splits
as a product of the completed local algebras at the points above zero.
Each factor B_q is a direct summand of a finite free R-module, hence is
free; reducing modulo (u,v) shows its rank is ell_q. These completions
are normal domains, since the original normal surface is excellent.
Apply the DVR degree formula at the generic point of the exceptional
line of Bl_(u,v)(Spec R), separately to this rank-ell_q extension.
It gives ell_q=sum e_i*f_i, counting just the components over q.
Additional resolution curves contracted to a point of E have f_i=0.
Completion leaves the exceptional curves and their direction maps
unchanged. Equation (I) then gives the other forms in (L).

The algebraic inputs are the henselian finite-algebra decomposition
[Stacks tag03QH](https://stacks.math.columbia.edu/tag/03QH), normality of
completed finite-type normal surface local rings
[tag0ADX](https://stacks.math.columbia.edu/tag/0ADX), and the already
pinned [DVR degree formula, tag0EXQ](https://stacks.math.columbia.edu/tag/0EXQ).
These are written geometric/algebraic arguments, not outputs of the checker.

Consequently every proposed local block must pass an integer condition:
(-Q)^(-1)f is a strictly positive integer vector. Its quadratic value
in (L) must equal the actual local length, and those lengths must sum
to d-1. Passing these conditions does not construct a normal cover or
a global Keller map.

## 3. Canonical coefficients distinguish horizontal and vertical curves

Let k_i=ord_Ci(dx wedge dy), computed on the resolved **source** model.
The Keller equation gives p^*(du wedge dv)=c*(dx wedge dy).
Near a point of the target exceptional line use radial/slope coordinates
(r,t), so du wedge dv is, up to sign, r*dr wedge dt. At a general
point of C_i choose s=0 as its local equation.

If f_i>0, write r=s^e_i*alpha with alpha a unit. The restriction of t
to C_i is nonconstant and separable. Thus dr wedge dt has order e_i-1,
and

\[
 \boxed{f_i>0\ \Longrightarrow\ k_i=2e_i-1.}\tag{K1}
\]

If f_i=0, write t-t0=s^n*beta, n>=1, with beta a unit at a general
point. The order of dr wedge dt is at least e_i+n-1, allowing leading
cancellation. Hence

\[
 \boxed{f_i=0\ \Longrightarrow\ k_i\ge2e_i.}\tag{K2}
\]

All these omitted radial curves therefore have positive canonical order.
These equations can be compared with A2B1's canonical vector reconstructed
from the **full** source boundary. Adjunction on the small block alone
does not reconstruct that vector: other boundary curves can meet the
block and contribute to its canonical intersection equations. Omitting
those attachments would create a false contradiction.

There is a further bookkeeping distinction at the actual source origin.
Apply A2B1 first on a completion preserving the original A2 source,
resolving the omitted points and infinity there. The subsequent blowup
of the source origin creates a disjoint exceptional P1 of canonical
order one. It leaves every omitted block and every infinity canonical
coefficient unchanged. The union of that new P1 with the infinity
boundary is not an A2 boundary; do not apply A2B1's connected-tree
test to that union. Its extra canonical coefficient is tracked separately.

## 4. A complete chain congruence and an excluded boundary stratum

Assume the resolved block is a chain C1--...--Cs with
Ci^2=-b_i, b_i>=2. Set A=-Q, so A has diagonal b_i and neighbouring
entries -1. Define prefix determinants

\[
 K_0=1,\quad K_1=b_1,\quad
 K_j=b_jK_{j-1}-K_{j-2},\quad \Delta=K_s.
\]

In the lattice cokernel Z^s/A Z^s, the first basis class generates:
the relations imply [epsilon_j]=K_(j-1)*[epsilon_1], and the last
relation is Delta*[epsilon_1]=0. Since det A=Delta, there are exactly
Delta classes. Thus (I) has an integral solution precisely when

\[
 \boxed{\sum_{j=1}^s K_{j-1}f_j\equiv0\pmod\Delta.}\tag{C}
\]

The inverse matrix has positive entries, so a nonzero f>=0 satisfying
(C) automatically gives e>0. For completeness, if T_(j+1) is the
determinant of the tail on j+1,...,s, with T_(s+1)=1, its entries are

\[
 (A^{-1})_{ij}=\frac{K_{\min(i,j)-1}T_{\max(i,j)+1}}\Delta.
\]

If just one component C_j is horizontal, of degree a, its allowed
degrees are precisely multiples of

\[
 h_j=\frac\Delta{\gcd(\Delta,K_{j-1})},\qquad
 \ell_q=\frac{a^2K_{j-1}T_{j+1}}\Delta.\tag{S}
\]

Because every b_i>=2, the positive prefix determinants strictly increase:
1=K0<K1<...<Ks=Delta. Thus 0<K_(j-1)<Delta and h_j>=2.
In particular **one radial horizontal component of degree one is
impossible on any such resolved chain, in every chain length**.
This is a closed stratum; no bounded chain enumeration is used.

For example, on the three-component (-2,-2,-2) chain, f=(0,1,0)
has the integral formal length f^T*A^(-1)*f=1, but its multiplicity
vector is (1/2,1,1/2). It is impossible despite passing the scalar
length test. Doubling f gives e=(1,2,1) and length four, realized by
the local N4 model below. The degree-two threshold is sharp among
these local controls, for every odd chain length in the even-n family.

Even for a single exceptional curve of self-intersection -b, allowing
b=1, (I) says f=b*e and (L) says ell_q=e*f. If f=1, necessarily
b=e=ell_q=1. A length-one completed factor B_q is R itself by Nakayama,
so this is an unramified missing normalization point. A ramified omitted
point with degree-one radial map needs a more complicated resolved block.

The chain hypothesis concerns the model where the radial ideal is
invertible. It must not be imposed on a minimal surface resolution that
still has pencil base points. The following all-degree control shows
exactly how the distinction matters.

## 5. Retained tame nodes: resolving the surface versus the radial pencil

The normal surface and finite cover

\[
 N_n=\{uv=w^n\}\longrightarrow\mathbb A^2_{u,v},\qquad n\ge2,
\]

are already used in [PF2K1TN1](../../foundations/plane-jc/F2_AFFINE_K1_TAME_NODE_PACKET.md).
That result supplies the minimal toric resolution and its (-2) chain;
we do not claim to rediscover them. The calculation here resolves the
additional radial ideal (u,v), computes e,f and the local fibre length,
and audits the parity-dependent scope of (C).

The cover has free basis 1,w,...,w^(n-1) over C[u,v], so the zero fibre
has length n. The surface is normal: it is an integral hypersurface,
smooth away from the isolated origin, hence Cohen-Macaulay and regular
in codimension one. Its source is a surface germ, not an affine-plane
Keller map. It is ramified over the target axes.

The retained minimal resolution has charts U_i, 0<=i<=n-1:

\[
 u=a^{i+1}b^i,\quad v=a^{n-i-1}b^{n-i},\quad w=ab.
\]

Adjacent chart changes are a'=1/b, b'=a*b^2. They show the exceptional
curves are P1 with normal degree -2. On C_j the valuations of (u,v)
are (j,n-j), so their minimum is e_j=min(j,n-j). The chart determinant
of (u,v) is n*a^(n-1)*b^(n-1).

### Even n=2m

The two monomial generators are comparable on every chart, so the
radial ideal is already invertible on the minimal resolution. Its
multiplicities are

\[
 e=(1,2,...,m-1,m,m-1,...,2,1).
\]

There is just one horizontal curve, C_m. On U_(m-1) the radial ratio
v/u=b^2, giving f_m=2. Every other f_j is zero. Directly,
(-Q)e=2*epsilon_m and e^T f=2m=n. The canonical order at C_m is
n-1=2e_m-1; all other curves have order n-1>=2e_j.
For n=2 this is Q=[-2], e=[1], f=[2]: radial ramification e=1
must be retained even though the finite surface cover is branched.

### Odd n=2m+1

The middle chart U_m has

\[
 u=a^{m+1}b^m,\quad v=a^mb^{m+1},\qquad
 (u,v)=a^mb^m(a,b).
\]

The ideal is **not** invertible at a=b=0. Its blowup is necessary.
On the chart b=a*t it becomes principal, with

\[
 u=a^n t^m,\quad v=a^n t^{m+1},\quad w=a^2t.
\]

The new exceptional curve has multiplicity e_new=n, maps to E by
v/u=t with degree f_new=1, and has self-intersection -1. Its two
neighbours change from -2 to -3. All original exceptional curves are
vertical. On this complete graph Qe=-f and e^Tf=n. Its new canonical
order is 2n-1, while the original curves keep order n-1>=2e_j.
The other blowup chart handles t=infinity and gives the same labels.

For n=3 the correctly ordered chain is

\[
 Q=\begin{pmatrix}-3&1&0\\1&-1&1\\0&1&-3\end{pmatrix},
 \quad e=(1,3,1),\quad f=(0,1,0).
\]

It is negative definite and has local length three. The one degree-one
horizontal curve is fully consistent because the resolved chain has
a -1 component. Before that blowup, the two (-2) curves are both
vertical and the pencil has a base point. Substituting just their
valuation minima in Qe=-f would falsely label them as horizontals.

These controls work for every n>=2. They show that local normal covers
can meet all the new budgets with arbitrary local fibre length. They
do not realize an omitted point of an actual global plane Keller map.
Global boundary attachment, the source A2 condition and the actual
constant-Jacobian map remain essential obligations.

## 6. Verification and remaining work

The checker uses A2B1's exact inertia routine, checks representative
chain congruences against direct integer matrix inversion, and replays
the node charts for n=2,...,12. Negative definiteness, the chain
congruence and both parity formulas have written proofs for all lengths
and all n; the finite regressions are not their quantifiers.

The next usable input is a complete candidate boundary block with its
radial ideal already principalized, its direction degrees and its global
canonical labels. The local chain stratum above is excluded, but an
extra blowup, another horizontal component or a non-chain graph can
escape that stratum. The singleton-zero plane Keller problem,
OP-HC4-SCHUR and JC2 remain open. No formal verification, external
review, new plane degree bound or novelty priority is claimed.
