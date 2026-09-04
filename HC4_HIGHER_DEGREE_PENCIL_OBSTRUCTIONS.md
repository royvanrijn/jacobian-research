# Higher-degree constant-Hessian pencil obstructions for HC4

## Status

This note continues the nonzero-corner scalar reverse-Schur branch after the
complete quadratic classification. It gives a global nilpotent formulation,
closes every pencil direction of generic Hessian rank one without a degree
bound, closes the genuinely moving leading-rank-three cubic direction, and
develops the constant-kernel tangent-ruling analysis through the complete
degree-six and degree-seven closures `HC4RSD32` and `HC4RSD40`.  The
degree-at-least-six statement in `HC4RSD24` is the historical frontier at that
stage, not the current endpoint of this note.  The later master reduction
`HC4MR1` consolidates the all-degree reductions but is partial after the
[transport correction](HC4_MOTION_FRAME_TRANSPORT_AUDIT.md): its final
negative maximal-motion sign remains open.

> **Theorem HC4RSD17 (nilpotent Hessian-pencil equivalence).** Let (K) have
> characteristic zero and let
>
> \[
> S=\operatorname{Hess}\psi,
> \qquad T=\operatorname{Hess}A,
> \qquad \det S=\delta\in K^\times.
> \tag{0.1}
> \]
>
> Then
>
> \[
> \det(S+sT)=\delta\quad\text{for every }s
> \tag{0.2}
> \]
>
> if and only if
>
> \[
> N=S^{-1}T=\delta^{-1}\operatorname{adj}(S)T
> \tag{0.3}
> \]
>
> is a nilpotent polynomial matrix. It automatically satisfies
>
> \[
> N^{\mathsf T}S=SN,
> \qquad S=\operatorname{Hess}\psi,
> \qquad SN=\operatorname{Hess}A.
> \tag{0.4}
> \]
>
> Thus the nonlinear scalar branch is exactly the classification of
> polynomial nilpotent endomorphisms that are self-adjoint for a unit
> Hessian metric and whose metric product is again Hessian-integrable.

> **Theorem HC4RSD18 (all-degree rank-one direction reduction).** Under
> (0.1)--(0.2), suppose
>
> \[
> \operatorname{rank}_{K(x)}\operatorname{Hess}A=1.
> \tag{0.5}
> \]
>
> Then every member of the pencil has injective gradient unless it is in the
> exact cotangent lift of a plane Keller map. More precisely, after scalar
> extension and a constant affine change,
>
> \[
> A=h(x)+\text{affine},
> \tag{0.6}
> \]
>
> and the two possible passive ternary packets are
>
> \[
> \psi=xw+C(x,y,z),
> \qquad \det\operatorname{Hess}_{y,z}C\in K^\times,
> \tag{0.7}
> \]
>
> or
>
> \[
> \psi=zP(x,y)+wQ(x,y)+R(x,y),
> \qquad \det J(P,Q)\in K^\times.
> \tag{0.8}
> \]
>
> The first is injective by HC2. The second is precisely the JC2 cotangent
> packet. There is no degree bound on (h), and no rationally moving ternary
> chart survives.

> **Theorem HC4RSD19 (moving leading-rank-three cubic obstruction).** Let
> (deg A\le3), suppose (operatorname{Hess}A) has generic rank three, and
> suppose the Hessian of the leading cubic part (A_3) also has rank three.
> If (0.2) holds, then (operatorname{Hess}A) has a constant kernel line.
> Equivalently, the sole genuinely moving normal form
>
> \[
> A=wz+y\,b(z)+G(x,z),
> \qquad \deg b\le2,
> \qquad b''\ne0,
> \tag{0.9}
> \]
>
> cannot occur in a constant-Hessian pencil. Its kernel is
>
> \[
> v=(0,1,0,-b'(z))^{\mathsf T}.
> \tag{0.10}
> \]
>
> The only residual leading-rank-three cubic packet is therefore
>
> \[
> A=a(x,y,z),
> \qquad
> \psi=w\,c(x,y,z)+D(x,y,z),
> \tag{0.11}
> \]
>
> with
>
> \[
> (\nabla c)^{\mathsf T}
> \operatorname{adj}(\operatorname{Hess}a)\nabla c=0.
> \tag{0.12}
> \]
>
> This ternary Hessian-eikonal packet is reduced further below; no moving
> cubic kernel remains in the stratum covered by the theorem.

> **Theorem HC4RSD20 (tangent-ruling synchronization and fixed-ruling
> collapse).** In the residual packet (0.11), impose the complete pencil
> equation, not only (0.12), and put
>
> \[
> p=\nabla c,qquad {\cal T}=\ker(p^{\mathsf T})
> \tag{0.13}
> \]
>
> over the fraction field. The restrictions
>
> \[
> E=\operatorname{Hess}D|_{\cal T},\qquad
> F=\operatorname{Hess}a|_{\cal T},\qquad
> G=\operatorname{Hess}c|_{\cal T}
> \tag{0.14}
> \]
>
> satisfy
>
> \[
> \det(E+sF+wG)=\det E\ne0
> \tag{0.15}
> \]
>
> up to the square of an invertible choice of tangent basis. The two singular
> binary forms (F) and (G) are proportional and have a common kernel line;
> that line is null for (E). If the common ruling is constant, then (c) is a
> cylinder and the packet reduces either to HC2 or to the exact JC2 cotangent
> lift. Thus no fixed-ruling packet is HC4-specific.

> **Theorem HC4RSD21 (polynomial fixed-ruling criteria).** In HC4RSD20 the
> ruling is constant in either of the following degree-unbounded/bounded
> classes:
>
> 1. (c) is homogeneous, with no bound on its degree;
> 2. (\deg c\le3), without a homogeneity assumption.
>
> Consequently both classes reduce to HC2 or the exact JC2 cotangent packet.
> The constant-kernel leading-rank-three cubic branch now begins only with a
> genuinely moving, nonhomogeneous ruling having (\deg c\ge4).

> **Theorem HC4RSD22 (quartic moving-ruling obstruction).** In HC4RSD20,
> every polynomial (c) of degree at most four is a fixed cylinder. For a
> non-pure binary quartic top, the first transverse face is a finite binary
> Schur equation. Its sole nonzero solution is the double--double chart, and
> the next face kills it by a pure fourth power. For a fourth-power top, the
> passive cubic is either a cube or zero; all three resulting correction
> charts end in an exact square or a rank-one passive quadratic. Therefore
> every quartic border coefficient again reduces to HC2 or the exact JC2
> cotangent packet. The genuinely moving nonhomogeneous residual begins at
> (\deg c\ge5).

> **Theorem HC4RSD23 (non-pure quintic moving-ruling obstruction).** In
> HC4RSD20, let (\deg c=5). If the leading quintic is not a fifth power of a
> linear form, then (c) is a fixed cylinder. At a simple root of its binary
> top, the first Schur numerator becomes a square and forces the transverse
> cubic to vanish there. The four- and five-distinct-root charts therefore
> vanish immediately. The (2+2+1) and (3+1+1) charts have zero radical; the
> only first-face solutions in the (3+2) and (4+1) charts are killed by the
> next face, and their possible lower transverse quadratics are killed by a
> pure square. Thus the only degree-five border packet not reduced to
> HC2/JC2 has pure fifth-power leading part.

> **Corollary HC4RSD24 (complete quintic ruling obstruction).** Every
> polynomial (c) of degree at most five in HC4RSD20 is a fixed cylinder.
> HC4RSD23 handles every non-pure leading quintic, while the pure-fifth
> chart is the quintic bordered lemma HC4BL5. Consequently every such
> residual packet reduces to HC2 or the exact JC2 cotangent packet. The
> first unresolved moving nonhomogeneous ruling has (\deg c\ge6).

> **Theorem HC4RSD25 (all-degree squarefree-top obstruction).** In
> HC4RSD20, let (c) have arbitrary degree (d), and suppose its leading
> binary form is squarefree. Then (c) is a fixed cylinder. At each simple
> root, the first transverse Schur numerator is a nonzero scalar times the
> square of the transverse coefficient. Since that coefficient has degree
> at most (d-2), the (d) simple roots force it to vanish. The same argument
> at the first remaining transverse layer inductively removes every
> dependence on the third variable. Thus, in every degree, the dense
> squarefree leading packet reduces to HC2 or the exact JC2 cotangent lift.

> **Theorem HC4RSD26 (generic-discriminant obstruction).** HC4RSD25
> remains true when the leading binary form has exactly one double root and
> every other root is simple. At the double root, the bordered determinant
> has local order four while the transverse Schur numerator has order two.
> Its leading coefficient is a nonzero characteristic-zero scalar times
> the square of the transverse coefficient at that root. Thus the double
> root joins the (d-2) simple roots in the vanishing set, giving (d-1)
> distinct roots for a polynomial of degree at most (d-2). Induction again
> kills every transverse layer. Hence a residual leading form must have at
> least two double roots or one root of multiplicity at least three.

> **Theorem HC4RSD27 (root-valuation resonance sieve).** Let a degree-(d)
> binary leading form have a root of multiplicity (m<d), and let (g) be a
> degree-(e\le d-2) coefficient in the first nonbinary weighted face. If
> (n=\operatorname{ord}g<m/2), polynomial Schur divisibility is possible
> only when
> \[
> C_{d,m,e,n}=d^2m+d^2n^2-2demn-2dem-dm^2-dm
>               +e^2m^2+2em^2+m^2=0.
> \]
> Thus every root supplies an explicit minimum vanishing weight, and the sum
> of those weights cannot exceed (e) unless (g=0). This finite sieve closes
> the sextic root partition (3+1+1+1) in every transverse degree.

> **Theorem HC4RSD28 (non-pure sextic moving-ruling obstruction).** In
> HC4RSD20, let (\deg c=6). If the leading binary sextic is not a sixth
> power, then (c) is a fixed cylinder. HC4RSD25--HC4RSD27 remove the
> squarefree, one-double, and (3+1+1+1) strata. For each of the seven
> remaining partitions, the valuation sieve leaves a finite weighted Schur
> face. The (2+2+1+1), (4+1+1), and leading (3+2+1) ideals are units; all
> apparent solutions in the (5+1), (2+2+2), (4+2), and (3+3) charts are
> killed by a positive-(z) coefficient of the same weighted initial face.
> In every weight-two chart the complete initial face includes the available
> scalar (z^3)-term; retaining it leaves the same zero transverse radicals.
> The resonant lower (3+2+1) face also has zero transverse radical. Hence
> the only unresolved scalar degree-six leading chart is a pure sixth power.

> **Theorem HC4RSD29 (pure-sixth passive-flag stabilization).** In the
> remaining chart of HC4RSD28, normalize (c_6=x^6). After a constant change
> of the two passive coordinates, the quintic correction (c_5) is binary.
> Indeed its passive Hessian is singular, so initially
> (c_5=H_5(x,y)+kx^4z). If both (k\ne0) and ((H_5)_{yy}\ne0), the next two
> faces leave only three normalized moving-direction charts. Their following
> homogeneous faces contain the immutable coefficients (1), (6), and
> (-10/9), respectively. Thus the moving collision is impossible. This
> stabilizes one passive flag through degree five; lower components may still
> try to break that flag, so HC4RSD29 is a narrowing theorem, not yet the
> complete pure-sixth bordered lemma.

> **Theorem HC4RSD30 (curved-quintic pure-sixth closure).** Continue in the
> pure-sixth chart after HC4RSD29, and write the stabilized quintic correction
> as (c_5=H_5(x,y)). If ((H_5)_{yy}\ne0), then the complete potential (c)
> is a fixed cylinder. For a first quartic direction break
> (c_4=R_4(x,y)+zP_3(x,y)), the next face is
> ((H_5)_{yy}(c_3)_{zz}=(P_3)_y^2). Unique factorization reduces every
> nonzero solution to five repeated-root cubic charts or the (y)-constant
> escape. Four charts have immediate immutable coefficients; the fifth and
> the constant escape die one face later by (-1) and (-4). Once (c_4) is
> binary, the later transverse cubic, quadratic, and affine tails vanish
> successively by pure-square faces. Hence the only scalar degree-six chart
> still open has passive-affine quintic correction
> (c_5=ax^5+x^4L(y,z)).

> **Theorem HC4RSD31 (curved completed-quartic closure).** Continue on the
> passive-affine boundary of HC4RSD30. After removing the unary quintic term,
> write (c_5=x^4L) for a constant passive linear form and define
> (\widehat c_4=c_4-x^2L^2/4). The first face makes (\widehat c_4) passive
> singular-Hessian. If its nonlinear part has nonzero passive curvature,
> then (c) is a fixed cylinder. In its normal form
> (\widehat c_4=Q_4(x,y)+kx^3z), every misaligned ruling dies in one of two
> ratios with terminal coefficients (2) and (5/108). Aligned transverse
> charts reduce to a localized two-equation resonance whose next face makes
> its ideal ((D,s)); aligned lower breaks die by (36), (-12), or the
> zero-Schur square cascade. The same charts with (L=0) close as well. Thus
> the only remaining degree-six scalar boundary has
> \[
> c_5=ax^5+x^4L,\qquad
> c_4=bx^4+x^3M+\frac{x^2L^2}{4}
> \]
> for constant passive linear forms (L,M).

> **Theorem HC4RSD32 (complete scalar pure-sixth closure).** Every potential
> in the two-linear-form tower left by HC4RSD31 is a fixed cylinder. If
> (L,M) are independent, normalize them to (y,z). The degree-twelve cubic
> packet is an irreducible three-parameter variety covered by two rational
> charts. In the finite chart its degree-ten descendants force (au=ap=0);
> the aligned degree-nine face then forces (p=u=0). The infinity chart dies
> by (-v^2/1296) and (-4r^2). At the sole base point, degree eight requires
> (q_5=1/12), while degree six has the immutable coefficient (-1/46656).
> If (L\ne0) and (M) is dependent on it, an (x)-translation reduces to
> (L=y,M=0). Five square faces make the cubic tail binary; the nonzero
> (z^2)-branch dies by (-4q_5^2), and the zero branch would require both
> (b_6=1/72) and (b_6=1/108) unless its last nonlinear transverse term
> vanishes. The remaining form is a cylinder. The (L=0,M\ne0) and
> (L=M=0) endpoints close by the same two cubic charts and the global
> identities
> \[
> [J(h(x,w)+vD(x))]_v=D h_{ww}(DD''-2(D')^2),\qquad
> J(h(x)+yf(x)+zg(x))=-(fg'-f'g)^2.
> \]
> Thus every scalar degree-six leading direction is fixed.

> **Theorem HC4RSD33 (non-pure septic closure).** In the synchronized
> constant-kernel packet of HC4RSD20, every degree-seven border coefficient
> whose leading binary septic is not a seventh power is a fixed cylinder.
> The root-valuation sieve directly removes the squarefree, one-double, and
> (3+1+1+1+1) partitions. It leaves eleven root partitions and exactly
> seventeen weighted faces, all in transverse degrees four and five. The
> complete degree-five face is (f+zg+z^2q/2+z^3r_1), so no same-weight
> cubic tail is omitted. Exact characteristic-zero coefficient ideals,
> localized at the one- and two-cross-ratio discriminants, contain an eighth
> power of every coefficient of (g,q,r_1). Thus all seventeen faces have
> the transverse origin as radical. The only scalar degree-seven leading
> chart still open is a pure seventh power.

> **Theorem HC4RSD34 (pure-septic two-face opening).** Normalize the sole
> scalar septic residue after HC4RSD33 to (c_7=x^7). Its degree-twenty face
> is
> \[
> 49x^{12}\det\operatorname{Hess}_{y,z}(c_6).
> \]
> Hence passive singular-Hessian classification gives
> (c_6=H_6(x,y)+kx^5z) after constant passive coordinates. The complete
> degree-nineteen face, with an arbitrary quintic correction retained, is
> \[
> \frac{49}{2}x^{12}(H_6)_{yy}
> \left(2(c_5)_{zz}-\frac87k^2x^3\right).
> \]
> Thus, when ((H_6)_{yy}\ne0),
> (c_5=R_5(x,y)+zP_4(x,y)+(2/7)k^2x^3z^2). This is a narrowing theorem:
> the descendants of this moving chart and the passive-affine
> ((H_6)_{yy}=0) boundary remain open.

> **Theorem HC4RSD35 (pure-septic degree-eighteen square obstruction).**
> Continue on the curved chart of HC4RSD34 and write
> (c_5=R_5+zP_4+(2/7)k^2x^3z^2). With every quartic coefficient retained,
> the complete next face is
> \[
> [J(c)]_{18}=x^8\left\{(H_6)_{yy}
> \left(49x^4(c_4)_{zz}-42kx^2P_4+18k^2H_6-6k^3x^5z\right)
> -\left(7x^2(P_4)_y-4k(H_6)_y\right)^2\right\}.
> \]
> It kills the (z^4)-coefficient of (c_4), fixes its complete cubic-passive
> tail to (k^3xz^3/49), and reduces the remainder to the binary identity
> \[
> (H_6)_{yy}(49x^4D_2-42kx^2P_4+18k^2H_6)
> =\left(7x^2(P_4)_y-4k(H_6)_y\right)^2.
> \]
> If (k\ne0), the coefficients below the (x^4D_2)-threshold force
> (a_6=a_5=0), so ((H_6)_{yy}) is divisible by (x^2). For (a_4\ne0) the
> two possible last ratios are (p_4=a_4k/7) and (p_4=5a_4k/14); the latter
> also forces (p_3=5a_3k/14). At (a_4=0), one has (p_4=0) and
> (p_3=2a_3k/7). Solving the three (D_2)-coefficients then decomposes each
> nonzero-(a_4) ratio into a generic resonance component and a double-root
> discriminant component; the (a_4=0) boundary has the (a_3\ne0) and
> (a_3=0) endpoints. Thus the moving chart has exactly six degree-eighteen
> packets. This is again a narrowing theorem: their degree-seventeen
> descendants, the (k=0) curved chart, and the passive-affine boundary of
> HC4RSD34 remain open.

> **Theorem HC4RSD36 (nonzero-(k) pure-septic closure).** All six packets
> of HC4RSD35 are impossible. The four nonzero-(a_4) packets and the
> (a_3\ne0) endpoint have immutable degree-seventeen coefficients; the
> pure-(x^4) endpoint leaves two ratios, killed in degree sixteen by
> (-24/49) and (256/1225).

> **Theorem HC4RSD37 (zero-Wronskian recursive-square closure).** On
> (k=0,(P_4)_y=0), the next face is the exact recurrence (11.63). Its
> nonzero-(p) solutions have three global UFD types, all killed by the
> normalized coefficients (-648/49), (16), (12/7), and (48/7). The
> zero-(p) tails close by (-4), (-24), and exact nonlinear-coordinate
> identities.

> **Theorem HC4RSD38 (coupled Wronskian closure).** On
> (k=0,(P_4)_y\ne0), retaining the cubic (z^3)-tail gives the coupled
> equations (11.58)--(11.59). Their complete solution is
> 
> \[
> (H_6)_{yy}=LE^3,\quad (P_4)_y=LE^2,\quad D_2=LE,
> \quad c=L_y/6.
> \]
> 
> All five ordered projective charts close. Thus every curved pure-septic
> chart is impossible.

> **Theorem HC4RSD39 (passive-affine curved closure and tower reduction).**
> On ((H_6)_{yy}=0), the shifted quintic face is (11.69). Every curved
> shifted quintic closes: the misaligned ratios die by (-20/21) and
> (36/7), while the aligned ordered-line charts die by degree fifteen or
> the incompatible factors (490k^2-35k+1) and (1449k-125). The remaining
> two-linear-form tower has the shifted passive-Hessian face (11.74), which
> leaves exactly eight degree-fifteen quartic direction packets.

> **Theorem HC4RSD40 (complete scalar pure-septic closure).** Every one of
> the eight packets left by HC4RSD39 is impossible or a fixed cylinder. Their
> degree-fourteen equations share a rank-one quartic-polar system. Its only
> nonzero resonance is
> 
> \[
> 3A_3^2=8A_2A_4,
> \]
> 
> equivalently ((Q_4)_{NN}) is a square. Five resonant direction charts die
> by the immutable coefficients (-3q^2/(49r)), (4b_0u), (-6u^2),
> (-u^6/(4b_0^2)), and (-147t_5^4/(4r^2)). Three transverse charts first
> force (A_4=A_3=0) and then die by (12/7) or (-12/2401). On the zero
> stratum, two charts contradict passive curvature and the other three have
> the global form (c=h(x,y)+z(a(x)y+d(x))) with (\deg a\le1); the complete
> bordered identity makes this a fixed cylinder. Together with HC4RSD33,
> every scalar degree-seven leading direction is therefore fixed and reduces
> to HC2 or the exact JC2 cotangent packet.

The exact determinant identities are replayed by
[scripts/verify_hc4_higher_degree_pencil_obstructions.py](scripts/verify_hc4_higher_degree_pencil_obstructions.py),
which writes
[artifacts/generated-results/hc4_higher_degree_pencil_obstructions.json](artifacts/generated-results/hc4_higher_degree_pencil_obstructions.json).
That artifact ends at the `HC4RSD28` stage.  Its `open_frontier` field is
historical discovery provenance, not the live frontier after `HC4RSD40` and
`HC4MR1`.

## 1. The global nilpotent deformation

Since (det S=delta) is a unit, (S^{-1}=delta^{-1}\operatorname{adj}S)
has polynomial entries. Hence (N=S^{-1}T) in (0.3) is polynomial and

\[
\det(S+sT)=\det S\det(I+sN).
\tag{1.1}
\]

Equation (0.2) is equivalent to

\[
\det(I+sN)=1.
\tag{1.2}
\]

Over the fraction field this says that every coefficient of the
characteristic polynomial of (N), except its leading term, vanishes.
Equivalently (N) is nilpotent. Conversely nilpotence gives (1.2).

Symmetry of (T=SN) gives

\[
SN=N^{\mathsf T}S.
\tag{1.3}
\]

The point is that (1.3) is not the complete structure: both (S) and
(SN) obey all third-derivative mixed-partial identities. This separates the
HC4 pencil problem from the classification of arbitrary self-adjoint
nilpotent matrices over a function field.

After shifting the pencil parameter, any collision may be placed at
(s=0). Thus the remaining marked problem is

\[
\nabla\psi(p)=\nabla\psi(q),
\qquad p\ne q,
\tag{1.4}
\]

together with the polynomial nilpotent and double-integrability equations
(0.3)--(0.4).

## 2. Rank one is degree-independent

The characteristic-zero rank-one polynomial-Hessian classification gives,
after scalar extension and constant affine coordinates,

\[
A=h(x)+\ell,
\tag{2.1}
\]

with (ell) affine. This is the rank-one case of the small-rank Hessian
classification in
[de Bondt's polynomial-Hessian paper](https://arxiv.org/abs/1609.03904).
Affine terms do not affect the pencil, so discard (ell).

Write

\[
S=
\begin{pmatrix}
k&d^{\mathsf T}\\
d&E
\end{pmatrix},
\qquad E=\operatorname{Hess}_{y,z,w}\psi.
\tag{2.2}
\]

The matrix determinant lemma, without any division by (h''), gives

\[
\det(S+s\operatorname{Hess}A)
=\det S+s\,h''(x)\det E.
\tag{2.3}
\]

Since (h''\ne0) and the polynomial ring is a domain,

\[
\det E=0.
\tag{2.4}
\]

The full determinant is a unit, so (E) has generic rank two. The argument
of HC4RSD15 now applies over (K(x)), and the unit-transverse globalization
obstruction HC4RSD16 applies verbatim: neither argument used that (h'')
was constant. The constant-kernel chart gives (0.7), while the exceptional
chart gives (0.8). This proves HC4RSD18.

In particular, raising the degree of a one-variable pivot direction cannot
produce an HC4-specific counterexample.

## 3. Classification of the leading-rank-three cubic direction

Write (A=A_3+A_{\le2}). The homogeneous four-variable Hesse theorem says
that the singular cubic (A_3), whose Hessian has rank three by hypothesis,
depends on three constant linear forms. Normalize

\[
A_3\in K[x,y,z].
\tag{3.1}
\]

Write the quadratic part involving (w) as

\[
\frac{\alpha}{2}w^2+w(a_1x+a_2y+a_3z)+Q(x,y,z).
\tag{3.2}
\]

The top-degree coefficient of (det\operatorname{Hess}A=0) is

\[
\alpha\det\operatorname{Hess}_{x,y,z}A_3,
\tag{3.3}
\]

so (alpha=0). Put (a=(a_1,a_2,a_3)). The next coefficient is

\[
-a^{\mathsf T}
\operatorname{adj}(\operatorname{Hess}_{x,y,z}A_3)a=0.
\tag{3.4}
\]

If (a=0), the (w)-direction is a constant Hessian kernel. Suppose
(a\ne0) and normalize (a=e_z), including a rescaling of (w). For the
complete ternary polynomial

\[
F=A-wz
\tag{3.5}
\]

the exact determinant, not only its leading part, is

\[
\det\operatorname{Hess}A
=-\det\operatorname{Hess}_{x,y}F.
\tag{3.6}
\]

The three entries of (operatorname{Hess}_{x,y}F) are affine linear
polynomials. In the UFD (K[x,y,z]),

\[
F_{xx}F_{yy}=F_{xy}^2
\tag{3.7}
\]

forces the nonzero entries to be constant multiples of one affine linear
form. Mixed-partial integrability fixes their projective direction. A
constant change of (x,y) therefore gives

\[
F=G(x,z)+y\,b(z).
\tag{3.8}
\]

Generic rank three forces (G_{xx}\ne0) and (b'\ne0). We have obtained
(0.9). If (b''=0), then

\[
wz+y\,b(z)=z(w+b'y)+\text{affine},
\tag{3.9}
\]

so the kernel is constant. The only genuinely moving case is (b''\ne0).

## 4. The moving cubic determinant faces

Put

\[
g=G_{xx},\qquad c=G_{xz},\qquad
d=G_{zz}+y b'',\qquad q=b'.
\tag{4.1}
\]

Then

\[
T=\operatorname{Hess}A=
\begin{pmatrix}
g&0&c&0\\
0&0&q&0\\
c&q&d&1\\
0&0&1&0
\end{pmatrix},
\qquad Tv=0,
\quad v=(0,1,0,-q)^{\mathsf T}.
\tag{4.2}
\]

Let (S=operatorname{Hess}\psi). The (s^3)-coefficient of
(det(S+sT)) is

\[
-g\,v^{\mathsf T}Sv.
\tag{4.3}
\]

Thus

\[
v^{\mathsf T}Sv=0.
\tag{4.4}
\]

Define the three remaining components of (Sv) by

\[
X=\psi_{xy}-q\psi_{xw},\qquad
U=\psi_{yw}-q\psi_{ww},\qquad
Z=\psi_{yz}-q\psi_{zw}.
\tag{4.5}
\]

Modulo (4.4), exact expansion of the (s^2)-face gives

\[
\boxed{(X-cU)^2+gU(2Z-dU)=0.}
\tag{4.6}
\]

This compact identity is the main cancellation obstruction.

## 5. Integrating the null direction

Let

\[
D_v=\partial_y-q(z)\partial_w,
\qquad r=w+q(z)y.
\tag{5.1}
\]

Since (D_vq=0), equation (4.4) is exactly (D_v^2\psi=0). Therefore

\[
\psi=yC(x,z,r)+D(x,z,r).
\tag{5.2}
\]

Write (m=b''\in K^\times). Direct differentiation gives

\[
U=C_r,\qquad X=C_x,
\qquad Z=C_z+2myC_r+mD_r.
\tag{5.3}
\]

Substitute (5.3) and (d=G_{zz}+my) in (4.6). Its coefficient of (y)
in the polynomial coordinates ((x,y,z,r)) is

\[
3m g C_r^2.
\tag{5.4}
\]

Characteristic zero, (m\ne0), and (g\ne0) give

\[
C_r=0.
\tag{5.5}
\]

Equation (4.6) now gives (C_x=0), so (C=C(z)), and

\[
Z=C'(z)+mD_r.
\tag{5.6}
\]

The (s^1)- and (s^0)-faces reduce exactly to

\[
-gZ^2D_{rr}=0,
\tag{5.7}
\]

\[
\delta=Z^2(D_{xr}^2-D_{rr}D_{xx}).
\tag{5.8}
\]

Because (delta) is a unit, (Z\ne0). Hence (5.7) gives (D_{rr}=0),
and (5.8) becomes

\[
\delta=(ZD_{xr})^2.
\tag{5.9}
\]

Both (Z) and (D_{xr}) are polynomial units. Write

\[
D=rL(x,z)+M(x,z).
\tag{5.10}
\]

Then (L_x=D_{xr}\in K^\times), while (5.6) says

\[
Z=C'(z)+mL(x,z)\in K^\times.
\tag{5.11}
\]

Differentiating (5.11) in (x) gives

\[
0=mL_x,
\tag{5.12}
\]

contradicting (m,L_x\in K^\times). This proves HC4RSD19.

## 6. The residual cubic packet

Only the constant-kernel leading-rank-three form remains. Normalize

\[
A=a(x,y,z).
\tag{6.1}
\]

The (s^3)-face is

\[
\det(\operatorname{Hess}a)\,\psi_{ww}=0,
\tag{6.2}
\]

so

\[
\psi=w\,c(x,y,z)+D(x,y,z).
\tag{6.3}
\]

The (s^2)-face is the ternary Hessian-eikonal equation (0.12). The complete
four-variable Hessian is

\[
 \operatorname{Hess}(\psi+sA)=
 \begin{pmatrix}
  \operatorname{Hess}(D+w c+s a)&p\\
  p^{\mathsf T}&0
 \end{pmatrix},
 \qquad p=\nabla c,
\tag{6.4}
\]

and hence

\[
 \det\operatorname{Hess}(\psi+sA)
 =-p^{\mathsf T}\operatorname{adj}
   (\operatorname{Hess}(D+w c+s a))p.
\tag{6.5}
\]

In particular the constant unit belongs to the square of the ideal generated
by the components of (p), so (p) is a unimodular row. Equation (0.12) alone
is not asserted to make (c) affine. The next section uses all three
directions (D,a,c) in (6.5), rather than only that one face.

## 7. Tangent-pencil synchronization

Work first over the fraction field (L=K(x,y,z)). Choose any basis of

\[
 {\cal T}=\ker(p^{\mathsf T}).
\tag{7.1}
\]

The bordered determinant (6.5) is, up to a nonzero square depending only on
that basis, the determinant of the restriction of the ternary Hessian to
({\cal T}). Thus (6.5) becomes the binary identity (0.15).

Write

\[
 F=\begin{pmatrix}f_0&f_1\\f_1&f_2\end{pmatrix},
 \qquad
 G=\begin{pmatrix}g_0&g_1\\g_1&g_2\end{pmatrix}.
\tag{7.2}
\]

The (s^2,w^2,sw)-coefficients of (0.15) are

\[
 f_0f_2-f_1^2=0,qquad
 g_0g_2-g_1^2=0,qquad
 f_0g_2+f_2g_0-2f_1g_1=0.
\tag{7.3}
\]

Every two-by-two minor of

\[
 \begin{pmatrix}f_0&f_1&f_2\\g_0&g_1&g_2\end{pmatrix}
\tag{7.4}
\]

has square in the ideal generated by (7.3). Therefore (F) and (G) are
linearly dependent over (L). Since (\operatorname{Hess}a) is generically
nonsingular, its restriction to a plane cannot vanish identically: a
nondegenerate ternary quadratic form has no two-dimensional totally
isotropic subspace. Hence (F) has rank one. The restrictions have one common
kernel line (\mathcal R\subset\mathcal T).

The coefficient linear in (s) says more. If

\[
 F=\binom{r}{t}(r\ t),
 \qquad k=(-t,r)^{\mathsf T},
\tag{7.5}
\]

then

\[
 [s]\det(E+sF)=k^{\mathsf T}Ek=0.
\tag{7.6}
\]

Thus the synchronized ruling is also null for the unit tangent metric (E).
This is the rank-two nilpotent shadow of HC4RSD17.

Suppose now that the ruling is represented by a nonzero constant vector in
(K^3). Since it lies in (\mathcal T), one has (k\cdot\nabla c=0), so after a
constant change

\[
 c=c(x,y).
\tag{7.7}
\]

Let (r=(-c_y,c_x,0)) be the other tangent direction and put

\[
 \kappa=r^{\mathsf T}\operatorname{Hess}(c)r.
\tag{7.8}
\]

The coefficient linear in (w) in (6.5) is exactly

\[
 D_{zz}\,\kappa.
\tag{7.9}
\]

If (\kappa=0), the two-variable straight-level lemma used in HC4RSD4 makes
(c) a polynomial in one linear form. Unimodularity of (\nabla c) then makes
that polynomial affine. The coordinate (\psi_w=c) recovers one source
variable, the remaining binary gradient has constant Hessian, and HC2 gives
injectivity.

If (\kappa\ne0), equation (7.9) gives (D_{zz}=0), so

\[
 \psi=w c(x,y)+zL(x,y)+M(x,y).
\tag{7.10}
\]

Direct expansion gives

\[
 \det\operatorname{Hess}\psi
 =\bigl(L_xc_y-L_yc_x\bigr)^2.
\tag{7.11}
\]

The plane map ((c,L)) is Keller, and (7.10) is exactly its cotangent lift.
This proves HC4RSD20.

## 8. When a polynomial ruling must be fixed

The equation supplied by the (w^2)-face is the universal-field equation

\[
 (\nabla c)^{\mathsf T}\operatorname{adj}(\operatorname{Hess}c)\nabla c=0.
\tag{8.1}
\]

Introduce one new variable (\tau). There is an exact identity

\[
 \det\operatorname{Hess}_{\tau,x,y,z}(\tau c)
 =-\tau^2(\nabla c)^{\mathsf T}
   \operatorname{adj}(\operatorname{Hess}c)\nabla c.
\tag{8.2}
\]

If (c) is homogeneous, then (\tau c) is a homogeneous form in four
variables with singular Hessian. The four-variable homogeneous Hesse theorem
gives a nonzero constant Hessian-kernel vector ((\alpha,q)). Splitting the
equation

\[
 \begin{pmatrix}0&(\nabla c)^{\mathsf T}\\
 \nabla c&\tau\operatorname{Hess}c\end{pmatrix}
 \binom{\alpha}{q}=0
\tag{8.3}
\]

in powers of (\tau) gives (\alpha=0) and

\[
 q\cdot\nabla c=0,qquad \operatorname{Hess}(c)q=0.
\tag{8.4}
\]

Thus (c) is a fixed cylinder in every homogeneous degree.

There is also a complete nonhomogeneous calculation through degree three.
For a homogeneous leading part (c_d), Euler's identities give

\[
 (\nabla c_d)^{\mathsf T}\operatorname{adj}(\operatorname{Hess}c_d)
 \nabla c_d
 =\frac d{d-1}c_d\det\operatorname{Hess}c_d.
\tag{8.5}
\]

Hence a cubic leading part has singular Hessian and, by the ternary Hesse
theorem, is binary after a constant change. Over an algebraic closure its
three binary root charts are

\[
 xy(x-y),\qquad x^2y,\qquad x^3/3.
\tag{8.6}
\]

Write the quadratic and linear tails as

\[
 \frac12(Ax^2+2Bxy+2Cxz+Dy^2+2Eyz+Fz^2)
 +Lx+My+Nz.
\tag{8.7}
\]

For the first two charts, the degree-five face is respectively

\[
 -6Fxy(x-y)(x^2-xy+y^2),\qquad -6Fx^4y.
\tag{8.8}
\]

Thus (F=0). The degree-four faces then become

\[
 -(Cx^2-2Cxy-2Exy+Ey^2)^2,qquad
 -x^2(Cx-2Ey)^2,
\tag{8.9}
\]

so (C=E=0). The remaining bordered equation is (N^2) times the
nonzero binary Hessian determinant, and hence (N=0). Both charts are fixed
cylinders.

For the cube chart, the passive ((y,z))-quadratic block first has rank at
most one. In its rank-zero normalization the complete equation is

\[
 -(BN-CM)^2=0,
\tag{8.10}
\]

so the quadratic and linear passive directions are proportional and (c)
again depends on two linear forms. In the rank-one normalization, with
passive quadratic part (y^2/2), the coefficients of (y^2) and then of (x)
are (-C^2) and (2N^2) after (C=0). Thus (C=N=0), giving another fixed
cylinder. For completeness, if (c=q_2+\ell) is quadratic and
(Q=\operatorname{Hess}c), rank three is excluded by the quadratic top of
(8.1). In rank two, normalize (Q=\operatorname{diag}(1,1,0)); equation
(8.1) is the square of the linear coefficient in the kernel direction, so
that coefficient vanishes and (c) is a cylinder. In rank at most one, the
quadratic direction and the linear direction span at most two constant
forms. Affine (c) is immediate. The resulting constant kernel descends from
the algebraic closure because it is the nullspace of a linear system over
(K).
This proves HC4RSD21.

## 9. The quartic transverse cancellation

Let

\[
 c=c_4+c_3+c_2+c_1,qquad \deg c_i=i.
\tag{9.1}
\]

Equation (8.5) and the ternary Hesse theorem make (c_4) binary. Choose its
top kernel to be (\partial_z). First suppose that the binary quartic (f=c_4)
is not a fourth power. Its binary bordered Hessian

\[
 {cal B}_f=
 \begin{pmatrix}
 0&f_x&f_y\\
 f_x&f_{xx}&f_{xy}\\
 f_y&f_{xy}&f_{yy}
 \end{pmatrix}
\tag{9.2}
\]

is nonsingular. The first transverse face forces ((c_3)_{zz}=0), so

\[
 c_3=r_3(x,y)+z g_2(x,y).
\tag{9.3}
\]

Put (q=(c_2)_{zz}) and (b_g=(g,g_x,g_y)^{\mathsf T}). The next face is the
finite binary Schur equation

\[
 q\det({\cal B}_f)=b_g^{\mathsf T}\operatorname{adj}({\cal B}_f)b_g.
\tag{9.4}
\]

Over an algebraic closure, the four non-pure root charts give

| root partition | (f) | solutions ((g,q)) |
|---|---|---|
| (1+1+1+1) | (xy(x-y)(x-\lambda y)), (\lambda(\lambda-1)\ne0) | ((0,0)) |
| (2+1+1) | (x^2y(x-y)) | ((0,0)) |
| (3+1) | (x^3y) | ((0,0)) |
| (2+2) | (x^2y^2) | ((bxy,b^2/4)) |

For the first three rows, applying the same Schur numerator to a remaining
linear (z)-coefficient (Ux+Vy) forces (U=V=0), and a constant coefficient
(N) gives (N=0). In the last row, include all binary cubic and quadratic
tails and all lower (z)-linear terms. The next homogeneous face contains

\[
 -\frac34 b^4x^2y^2z^2,
\tag{9.5}
\]

with no other contribution to that monomial. Thus (b=0), after which the
linear and constant calculation applies. Every non-pure quartic top is a
fixed cylinder.

It remains to treat the fourth-power top. Normalize

\[
 c_4=x^4/12,qquad
 c_3=P_3(y,z)+xQ_2(y,z)+x^2\ell(y,z)+k x^3.
\tag{9.6}
\]

Write

\[
 P_3=p_0z^3+p_1yz^2+p_2y^2z+p_3y^3,qquad
 Q_2=q_0z^2+q_1yz+q_2y^2.
\tag{9.7}
\]

The first nonzero face consists exactly of

\[
\begin{aligned}
3p_0p_2-p_1^2&=0,&9p_0p_3-p_1p_2&=0,
&3p_1p_3-p_2^2&=0,\\
3p_0q_2-p_1q_1+p_2q_0&=0,
&p_1q_2-p_2q_1+3p_3q_0&=0,
&4q_0q_2-q_1^2&=0.
\end{aligned}
\tag{9.8}
\]

The first row says that (P_3) is a cube or zero. If it is nonzero, normalize
(P_3=y^3/3); the second row then gives (Q_2=a y^2). Write

\[
 \ell=uy+vz,qquad
 c_2=\frac12(Ax^2+2Bxy+2Cxz+Dy^2+2Eyz+Fz^2).
\tag{9.9}
\]

The next two faces are

\[
 F=3v^2,qquad
 -\frac{x^2}{9}
 \bigl((-E+3uv)x^2+6avxy+3vy^2\bigr)^2=0.
\tag{9.10}
\]

They force (v=E=F=0). The following face has coefficient (-C^2) at
(y^4), so (C=0); the remaining equation is (N^2) times a nonzero binary
Hessian determinant, hence (N=0).

If (P_3=0) but (Q_2\ne0), normalize (Q_2=a y^2), (a\ne0). Equations
(9.10) are replaced by

\[
 F=3v^2,qquad
 -\frac{x^4}{9}
 \bigl((-E+3uv)x+6avy\bigr)^2=0.
\tag{9.11}
\]

Again (v=E=F=0). The next coefficient is
(2aC^2x^5/3), followed by the nonzero binary Hessian multiple (N^2), so
(C=N=0).

Finally take (P_3=Q_2=0). If (\ell\ne0), normalize (\ell=y). The first
passive quadratic relation is

\[
 DF-E^2-3F=0.
\tag{9.12}
\]

Were (F\ne0), a passive shear would give (E=0,D=3), but the next face has
coefficient (9F) at (x^2y^2), a contradiction. Hence (E=F=0). The
coefficient at (x^4) is (C^2(D-3)/3). If (D\ne3), this kills (C); if
(D=3), the coefficient (-9C^2) at (y^2) does so. In either case the next
binary-Hessian coefficient kills (N).

If also (\ell=0), translate (x) to remove (kx^3). The passive quadratic
block has determinant zero. In passive rank zero the complete equation is

\[
 -(BN-CM)^2=0,
\tag{9.13}
\]

so the two passive directions are proportional. In passive rank one,
normalize the passive term to (y^2/2); the coefficients (C^2x^4/3) and,
after (C=0), (N^2x^2) force (C=N=0). These cases exhaust (9.8) and prove
HC4RSD22.

## 10. The non-pure quintic top

Let (\deg c=5). The leading part is again binary; write it as (f(x,y)). If
(f) is not a fifth power, its bordered Hessian (\mathcal B_f) from (9.2) is
nonsingular. The first transverse face gives

\[
 c_4=r_4(x,y)+z g_3(x,y),qquad
 q_1\det(\mathcal B_f)
 =b_g^{\mathsf T}\operatorname{adj}(\mathcal B_f)b_g,
\tag{10.1}
\]

where (q_1=(c_3)_{zz}) is binary linear and
(b_g=(g,g_x,g_y)^{\mathsf T}).

There is a useful all-degree local obstruction behind (10.1). At a simple
root, choose coordinates with (f=xh), (h\ne0), and evaluate at (x=0). Then

\[
 \mathcal B_f=
 \begin{pmatrix}0&h&0\\h&a&b\\0&b&0\end{pmatrix},
 \qquad
 b_g^{\mathsf T}\operatorname{adj}(\mathcal B_f)b_g
 =-(bg-hg_y)^2.
\tag{10.2}
\]

Since (g) has degree three while (f) has degree five, Euler's identity makes
the last square a nonzero scalar multiple of (g^2) at the root. Thus every
simple root of (f) is a root of (g). If (f) has four or five distinct roots,
then (g=0), and (10.1) gives (q_1=0). The same argument at the lower
transverse layers then kills every remaining (z)-coefficient.

The other non-pure root partitions give the exact table

| root partition | (f) | first-face solutions ((g,q_1)) |
|---|---|---|
| (2+2+1) | (x^2y^2(x-y)) | ((0,0)) |
| (3+1+1) | (x^3y(x-y)) | ((0,0)) |
| (3+2) | (x^3y^2) | ((b x^2y,11b^2x/30)) |
| (4+1) | (x^4y) | ((b x^2y,b^2y/5)) |

For the first two rows, the coefficient ideal has radical
((g,q_1)); the subsequent transverse quadratic has degree two and vanishes
at all three distinct roots, so all lower layers vanish as well.

In the (3+2) row, include arbitrary binary quartic corrections, the next
binary quadratic (z)-coefficient, and the available scalar (z^2)-term. The
next face has the immutable coefficient

\[
 -\frac8{15}b^3x^7y^3z.
\tag{10.3}
\]

In the (4+1) row the corresponding coefficient is

\[
 \frac65b^3x^8y^2z.
\tag{10.4}
\]

Hence (b=0) in both rows. A remaining transverse quadratic must vanish at
the two simple-root directions and is therefore a multiple (a xy). Its
next face is respectively

\[
 -a^2x^6y^4,qquad -9a^2x^8y^2,
\tag{10.5}
\]

so (a=0). Linear and constant transverse tails then vanish at the same two
roots. This exhausts every non-pure binary quintic and proves HC4RSD23.

The sole quintic residue is

\[
 c_5=\ell^5
\tag{10.6}
\]

for one constant linear form (\ell), with nonhomogeneous lower corrections.

This residue is closed by the pure-fifth calculation in
[HC4_QUINTIC_COMMON_DIRECTION.md](HC4_QUINTIC_COMMON_DIRECTION.md). Its
first face is the passive Hessian determinant of the quartic correction.
The curved correction charts have unit full lower-tail ideals; the
passive-affine and transverse charts force every lower term to use one
constant passive linear form. This is HC4BL5 and proves HC4RSD24.

The simple-root computation is not specific to degree five. If the binary
top (f) has degree (d), then at a simple root choose (f=xh). For a
homogeneous transverse coefficient (g) of degree (e\le d-2), equation
(10.2) and Euler give

\[
 bg-hg_y=(d-1-e)hg.
\tag{10.7}
\]

The scalar is nonzero in characteristic zero. If (f) is squarefree, all
(d) roots divide (g), so (g=0). At the first nonzero lower transverse
layer the same equation applies with a still smaller (e), and induction
kills that layer as well. Hence every transverse layer vanishes. This is
the degree-unbounded HC4RSD25.

At a double root, work in the affine chart (y=1), write
(f=x^2(h_0+h_1x+\cdots)), and let (g) have degree (e\le d-2). The bordered
determinant starts in order four, whereas its Schur numerator starts with

\[
 -2h_0^2g(0)^2E_{d,e}x^2,
\qquad
 E_{d,e}=d^2-2de-3d+2e^2+4e+2.
\tag{10.8}
\]

If (k=d-1-e\ge1), then

\[
 E_{d,e}=2(k-d/2)^2+d(d-2)/2>0.
\tag{10.9}
\]

Thus the double root also divides (g). With one double root and all other
roots simple, (g) has (d-1) distinct roots but degree at most (d-2), so it
vanishes. The same induction proves HC4RSD26.

The same local computation works at a root of arbitrary multiplicity. Write

\[
 f=x^m(h_0+O(x)),\qquad g=x^n(u_0+O(x)),\qquad h_0u_0\ne0.
\]

Euler reconstruction of the homogeneous derivatives gives

\[
 \det\mathcal B_f
 =d m(d-m)h_0^3x^{3m-2}+O(x^{3m-1}),
\tag{10.10}
\]

and

\[
 b_g^{\mathsf T}\operatorname{adj}(\mathcal B_f)b_g
 =-h_0^2u_0^2 C_{d,m,e,n}x^{2n+2m-2}
   +O(x^{2n+2m-1}),
\tag{10.11}
\]

where

\[
 C_{d,m,e,n}=d^2m+d^2n^2-2demn-2dem-dm^2-dm
              +e^2m^2+2em^2+m^2.
\tag{10.12}
\]

If (n<m/2) and (C_{d,m,e,n}\ne0), the right side of the Schur equation has
smaller root order than its determinant factor, which is impossible. Define

\[
 \rho_{d,e}(m)=\min\left(
 \left\lceil\frac m2\right\rceil,
 \{\,n<m/2:C_{d,m,e,n}=0\,\}
 \right),
\tag{10.13}
\]

where the minimum of the empty set is omitted. Every transverse coefficient
must vanish to order at least (\rho_{d,e}(m)) at that root. Consequently

\[
 \sum_i\rho_{d,e}(m_i)>e\quad\Longrightarrow\quad g=0.
\tag{10.14}
\]

This records the exceptional cancellations rather than hiding them in a
genericity assumption. For the first layer (e=d-2), there are no resonances
when (d\ge4) and (m<d). As a quadratic in (n), (C) is decreasing on
(n<m/2). At the last possible odd value,

\[
 4C_{d,m,d-2,(m-1)/2}=(dm-d-2m)^2,
\tag{10.15}
\]

and at the last even value, writing (d=m+a) with (a>0),

\[
 4C_{d,m,d-2,m/2-1}
 =a^2(m^2+4)+2am(m^2-2m+2)+m^2(m-2)^2>0.
\tag{10.16}
\]

For a sextic, the root-weight sums for transverse degrees (e=0,1,2,3,4)
are

| partition | root-weight sums |
|---|---|
| (3+1+1+1) | (5,5,5,4,5) |
| (2+2+1+1) | (4,4,4,4,4) |
| (4+1+1) | (4,4,4,4,4) |
| (3+2+1) | (4,4,4,3,4) |
| (5+1) | (4,4,4,4,4) |
| (2+2+2) | (3,3,3,3,3) |
| (4+2) | (3,3,3,3,3) |
| (3+3) | (4,4,4,2,4) |

The first row vanishes immediately. In every other row, (10.14) leaves only
the top degree (e=4), together with the displayed resonant (e=3) cases.
Assign weights

\[
 \operatorname{wt}(x)=\operatorname{wt}(y)=1,
 \qquad \operatorname{wt}(z)=d-e.
\]

At the first nonbinary layer the complete initial potential is

\[
 c_{\mathrm{in}}=f+zg+\frac12z^2q+
 \sum_{k\ge3}z^k r_k,
 \qquad \deg q=2e-d,\quad \deg r_k=d-k(d-e),
\tag{10.17}
\]

with any term omitted when its displayed degree is negative. In particular,
the sextic top layer (e=4) includes a scalar (z^3)-term. The constant-(z) part of
(J(c_{\mathrm{in}})=0) is the Schur equation, but every positive-(z)
coefficient belongs to the same initial face and must vanish too. Exact
coefficient ideals give

| partition/layer | locally admissible (g) | Schur result |
|---|---|---|
| (2+2+1+1), (e=4) | squarefree radical | unit ideal, uniformly in the cross-ratio |
| (4+1+1), (e=4) | (x^2y(x-y)) | unit ideal |
| (3+2+1), (e=4) | (x^2y(x-y)) | unit ideal |
| (3+2+1), (e=3) | (xy(x-y)) | zero transverse radical |
| (5+1), (e=4) | (b x^3y) | (q=11b^2xy/30), then (b=r_3=0) in the full face |
| (2+2+2), (e=4) | (xy(x-y)(ux+vy)) | radical ((u,v)) |
| (2+2+2), (e=3) | (bxy(x-y)) | (q=b^2/6), then (b=0) |
| (4+2), (e=4) | (x^2y(ux+vy)) | unique (q), then radical ((u,v,r_3)) |
| (4+2), (e=3) | (b x^2y) | (q=b^2/6), then (b=0) |
| (3+3), (e=4) | (b x^2y^2) | (q=b^2xy/2), then radical ((b,r_3)) |
| (3+3), (e=3) | (xy(ux+vy)) | (q=2uv/3), then radical ((u,v)) |

Lower layers are excluded by (10.14). These are all partitions of a
non-pure binary sextic, proving HC4RSD28.

## 11. The pure-sixth collision

It remains to take (c_6=x^6). For an arbitrary homogeneous quintic (h), the
degree-sixteen face is

\[
 [J(x^6+h)]_{16}=36x^{10}\det\operatorname{Hess}_{y,z}(h).
\tag{11.1}
\]

The binary Hesse theorem over (K(x)), together with ordinary homogeneity,
puts the singular passive Hessian in the constant-direction form

\[
 c_5=H_5(x,y)+kx^4z.
\tag{11.2}
\]

For a general quartic correction, the next face factors globally as

\[
 [J(c)]_{15}
 =18x^{10}(H_5)_{yy}\bigl(2(c_4)_{zz}-k^2x^2\bigr).
\tag{11.3}
\]

Suppose (k(H_5)_{yy}\ne0). Then

\[
 c_4=R_4(x,y)+zP_3(x,y)+\frac{k^2}{4}x^2z^2.
\tag{11.4}
\]

Write (H_5=\sum_{i=0}^5a_ix^{5-i}y^i) and let (p_3) be the coefficient of
(y^3) in (P_3). Four coefficients of the degree-fourteen face give, in
order,

\[
 -25a_5^2k^2,qquad -24a_4^2k^2,qquad
 -3(18p_3-7a_3k)(6p_3-a_3k),
\tag{11.5}
\]

and

\[
 (a_2,a_3)(108(c_3)_{z^3}-k^3)=0.
\tag{11.6}
\]

After passive translations and the scaling torus, the remaining cases are
one quadratic-passive chart and two cubic-passive charts, corresponding to
(p_3=a_3k/6) and (p_3=7a_3k/18). Keep every arbitrary quartic, cubic, and
quadratic tail capable of entering the next face. The coefficients at
(x^9y^4), (x^8y^4z), and (x^8y^4z) in those three charts are respectively

\[
 1,qquad 6,qquad -\frac{10}{9}.
\tag{11.7}
\]

Thus no chart survives. If (k=0), (c_5) was binary already. If
((H_5)_{yy}=0), its passive part is affine and a constant passive change
again makes (c_5) binary. This proves HC4RSD29.

### 11.1. The lower flag when the quintic correction is curved

Assume now that (c_5=H_5(x,y)) and ((H_5)_{yy}\ne0). Equation (11.3) with
(k=0) gives

\[
 c_4=R_4(x,y)+zP_3(x,y).
\tag{11.8}
\]

Write ((c_3)_{zz}=Q_1(x,y)). The next face factors without a remainder:

\[
 [J(c)]_{14}=36x^{10}
 \left((H_5)_{yy}Q_1-(P_3)_y^2\right).
\tag{11.9}
\]

This is a useful second root-valuation obstruction. If ((P_3)_y\ne0), the
binary cubic ((H_5)_{yy}) has at most one irreducible factor of odd
multiplicity. Under the affine group preserving the pure direction (x),
the possibilities are

\[
 y^3,\quad x^3,\quad xy^2,\quad x^2y,\quad y^2(y-x).
\tag{11.10}
\]

Keeping all terms capable of entering degree thirteen, four charts contain
the immutable coefficients

\[
 \frac5{24},\qquad 36,\qquad -36,\qquad 36.
\tag{11.11}
\]

The remaining triple-root-at-infinity chart has a resonant initial package
containing (x(xy+z)^2/2), but after solving its degree-thirteen equations
the coefficient at (x^8y^4) in degree twelve is (-1).

If ((P_3)_y=0), homogeneity gives (P_3=px^3). For (p\ne0), the next face
first forces ((H_5)_{yy}) to be proportional to (x^3). After normalization,
the following face has coefficient (-4) at (x^{10}z). Thus (P_3=0), and
(Q_1=0) follows from (11.9).

The rest is triangular. Write the next possible transverse cubic as
(zT_2(x,y)). Successive homogeneous faces are

\[
\begin{aligned}
 [J(c)]_{13}&=72x^{10}(H_5)_{yy}(c_2)_{z^2},\\
 [J(c)]_{12}&=-36x^{10}(T_2)_y^2,\\
 [J(c)]_{11}&=6x^8(H_5)_{yy}[T_2]_{x^2}^2,\\
 [J(c)]_{10}&=-36x^{10}[c_2]_{yz}^2,\\
 [J(c)]_9&=18x^6(H_5)_{yy}[c_2]_{xz}^2.
\end{aligned}
\tag{11.12}
\]

They remove every nonlinear (z)-term. Finally, for
(c=h(x,y)+\ell z),

\[
 J(c)=\ell^2\det\operatorname{Hess}_{x,y}(h),
\tag{11.13}
\]

whose leading term is nonzero when ((H_5)_{yy}\ne0). Hence (\ell=0), and
(c) is binary. This proves HC4RSD30.

### 11.2. Completing the passive-affine quintic square

It remains to take

\[
 c_6=x^6,\qquad c_5=x^4L(y,z)
\tag{11.14}
\]

after an (x)-translation removes the unary quintic term. The complete next
face is

\[
 [J(c)]_{14}=36x^{10}
 \det\operatorname{Hess}_{y,z}
 \left(c_4-\frac{x^2L^2}{4}\right).
\tag{11.15}
\]

Thus, after passive coordinates,

\[
 \widehat c_4:=c_4-\frac{x^2L^2}{4}
 =Q_4(x,y)+kx^3z,qquad L=Ay+Bz.
\tag{11.16}
\]

Suppose ((Q_4)_{yy}\ne0). If (B\ne0), shear (z) to set (A=0) and scale
(B=1). The next face removes the (y^4) and (y^3) coefficients of (Q_4),
and the two surviving curvature ratios have immutable coefficients

\[
 2,qquad \frac5{108}.
\tag{11.17}
\]

If (B=0) and (k\ne0), normalize (A=k=1). Put (D=12q-1), where (q) is the
available scalar (z^2)-coefficient. The (D=0) face loses all curvature in
the next two degrees. On (D\ne0), a finite square-root direction is killed
by (12vD). At the root at infinity, the next two faces give

\[
\begin{aligned}
 E_1&=D^2-2Ds+2s,\\
 E_2&=D^2-4Ds+D-4s^2,\\
 F&=4D^3-10D^2s+3D^2+16Ds^2-4Ds-4s^2.
\end{aligned}
\tag{11.18}
\]

Their coefficient ideal is exactly ((D,s)), contradicting the localized
chart. If (k=0), a first lower cubic break satisfies another square-factor
equation. Its finite and infinite charts die by (36) and (-12). The
zero-scalar chart has possible ratios cut out by

\[
 s(6s-1)^2(18s+1)=s(12s-1)=0,
\tag{11.19}
\]

so no nonzero (s) survives. Its later quadratic and affine tails vanish by
the same pure-square cascade as (11.12).

When (L=0), the identical split applies with (D=12q-k^2). The (D=0)
coefficient row directly removes every curved coefficient; the localized
(D\ne0) row has the square (-108s^4/D^2). With (k=0), the terminal finite,
infinite, and zero-scalar coefficients are (36), (-1), and (-648s^4).
Therefore ((Q_4)_{yy}\ne0) always gives a fixed cylinder. This proves
HC4RSD31. Its complement is precisely

\[
 c_5=ax^5+x^4L,qquad
 c_4=bx^4+x^3M+\frac{x^2L^2}{4}.
\tag{11.20}
\]

### 11.3. Closing the two-linear-form tower

Unary terms in (11.20) do not enter the argument. If (L\ne0) and
(M=kL), translation in (x) changes the coefficient of (x^3L) by (4s), so
it also lets us put (M=0). There are consequently four normalizations:

\[
(L,M)=(y,z),\quad(y,0),\quad(0,z),\quad(0,0).
\tag{11.21}
\]

First take the independent case. The relevant degree-twelve coefficient
ideal is radical-set-theoretically one irreducible three-dimensional packet.
Its finite chart is

\[
\begin{aligned}
c_3={}&b_0x^3+b_1x^2y+b_2x^2z+p x(y+az)^2
       +\frac13xyz+\frac1{108}y^3\\
&+\frac{u}{108}(y+az)^3,
\end{aligned}
\tag{11.22}
\]

and its chart at infinity replaces the last curved terms by

\[
r x(hy+z)^2+\frac13xyz+\frac1{108}y^3
 +\frac{v}{108}(hy+z)^3.
\tag{11.23}
\]

On (11.22), degree eleven is the single equation

\[
E=-2a^2b_1+12a^2q_3+2ab_2-12aq_4+12q_5-1=0
\tag{11.24}
\]

whenever ((p,u)\ne(0,0)). After using (E), put

\[
R=-2ab_1+12aq_3+b_2-6q_4.
\]

The degree-ten face contains

\[
-R^2,\qquad -\frac{a^2u^2}{1296},\qquad
\frac{auR}{18}-4a^2p^2.
\tag{11.25}
\]

Thus (a=0). Then (R=0), and degree nine contains (-4p) and (-u/9),
a contradiction. On (11.23), degree ten instead contains
(-v^2/1296), followed at (v=0) by (-4r^2). The only cubic point left is

\[
c_3=b_0x^3+b_1x^2y+b_2x^2z+\frac13xyz+\frac1{108}y^3.
\tag{11.26}
\]

Its degree-eight coefficient at (x^4y^4) is
(5(12q_5-1)/1296), whereas its degree-six coefficient at (y^6) is
((9q_5-1)/11664). The first fixes (q_5=1/12), and the second becomes the
immutable (-1/46656). Hence the independent tower is empty.

For ((L,M)=(y,0)), successive degree-ten coefficients are

\[
-9b_9^2,quad-112b_8^2,quad-109b_7^2,quad-4b_5^2,quad-36b_4^2.
\tag{11.27}
\]

The cubic tail is therefore binary apart from (b_2x^2z). If (q_5\ne0),
the next faces set (b_3=0,b_6=1/108) and then contain (-4q_5^2). If
(q_5=0), they first set (q_4=b_2/6), and two later coefficients are

\[
-\frac12b_2^3(72b_6-1),\qquad
-\frac1{18}b_2^3(108b_6-1).
\tag{11.28}
\]

Thus (b_2=0). What remains has the form (h(x,y)+z(\alpha x+\beta)). The
coefficient of (z) in its bordered invariant is

\[
-2\alpha^2(\alpha x+\beta)h_{yy}.
\tag{11.29}
\]

Here (h_{yy}) has leading term (x^2/2), so (\alpha=0) and the potential is
a cylinder.

For ((L,M)=(0,z)), the unshifted cubic packet has finite chart

\[
c_3=b_0x^3+b_1x^2y+b_2x^2z+p x(y+az)^2+u(y+az)^3.
\tag{11.30}
\]

Writing (S=a^2q_3-aq_4+q_5), a curved member first requires (S=1/12).
The next face contains (-36u^2), and after (u=0) it contains (-12p^2).
The infinity chart has the corresponding coefficients
(-36h^2v^2) and (-12h^2r^2). Its (h=0) endpoint successively loses the
(y)-coefficients (b_1,q_1,\ell_1) and is a cylinder. At the base cubic,
the two passive-quadratic equations are

\[
12q_3q_5-q_3-3q_4^2=0,qquad
20q_3q_5-2q_3-5q_4^2=0.
\tag{11.31}
\]

Their difference after multiplying the first by (5/3) is (-q_3/3), so
(q_3=q_4=0). Two pairs of square coefficients then give (b_1=q_1=0),
again leaving an affine cylinder direction.

Finally take ((L,M)=(0,0)). A nonzero packet (11.30) makes the passive
quadratic and the (x^2)-linear passive term use the same constant form
(w=y+az). Any attempted remaining break is (h(x,w)+vD(x)), and direct
adjugate expansion gives

\[
[J(h+vD)]_v=D h_{ww}\bigl(DD''-2(D')^2\bigr).
\tag{11.32}
\]

For a quadratic (D), the coefficients of the last factor begin with
(-6\gamma^2x^2) and, after (\gamma=0), (-2\delta^2). Thus (D) is constant.
If the passive quadratic is zero as well, write
(c=h(x)+yf(x)+zg(x)). Then

\[
J(c)=-(fg'-f'g)^2.
\tag{11.33}
\]

The Wronskian vanishes only when (f,g) are constant-proportional, so this
last endpoint is also a cylinder. This proves HC4RSD32 and completes every
scalar pure-sixth chart.

### 11.4. The complete non-pure septic row

For (d=7), the valuation weight is independent of the transverse degree
(0\le e\le5) on every partition that survives HC4RSD25--HC4RSD27. The
three partitions

\[
(1^7),\qquad(2+1^5),\qquad(3+1^4)
\tag{11.34}
\]

have weight sums (7,6,6), respectively, and therefore vanish directly.
The remaining partitions split into

| root partitions | root-weight sum | possible transverse degrees |
|---|---:|---|
| (6+1), (5+2), (4+3), (4+2+1), (3+2+2), (2+2+2+1) | 4 | (4,5) |
| (5+1+1), (4+1+1+1), (3+3+1), (3+2+1+1), (2+2+1+1+1) | 5 | (5) |

Thus there are exactly seventeen faces. At a root of multiplicity (m),
let (R) contain that root to order (\lceil m/2\rceil). The admissible
transverse coefficient is

\[
g=R G_{e-\deg R},
\tag{11.35}
\]

where (G) is a generic binary scalar or linear form. The complete weighted
potential is

\[
\begin{array}{ll}
e=4:&f+zg+\frac12z^2q_1,\\[2mm]
e=5:&f+zg+\frac12z^2q_3+z^3r_1.
\end{array}
\tag{11.36}
\]

In particular, the binary-linear cubic tail (r_1) is not discarded. For
up to three distinct roots, normalize them to (0,\infty,1). Four distinct
roots add a cross-ratio (a) and localize at (a(a-1)); five distinct roots
add (a,b) and localize at

\[
ab(a-1)(b-1)(a-b).
\tag{11.37}
\]

For each of the seventeen cases, the exact characteristic-zero coefficient
ideal of (J(c_{\rm in})) contains the eighth power of every coefficient of
(g,q_1), or of (g,q_3,r_1), after the indicated localization. Therefore its
radical is the complete transverse coefficient origin. Lower degrees have
root-weight sum greater than (e), so no later break remains. This proves
HC4RSD33.

### 11.5. Opening the pure-seventh chart

Set (c_7=x^7) and retain a completely generic sextic correction. Direct
homogeneous extraction gives

\[
[J(c)]_{20}=49x^{12}\det\operatorname{Hess}_{y,z}(c_6).
\tag{11.38}
\]

The passive singular-Hessian normal form is

\[
c_6=H_6(x,y)+kx^5z.
\tag{11.39}
\]

With every coefficient of (c_5) independent, the next complete face is

\[
[J(c)]_{19}=\frac{49}{2}x^{12}(H_6)_{yy}
\left(2(c_5)_{zz}-\frac87k^2x^3\right).
\tag{11.40}
\]

Consequently, on the curved chart ((H_6)_{yy}\ne0),

\[
c_5=R_5(x,y)+zP_4(x,y)+\frac27k^2x^3z^2.
\tag{11.41}
\]

Unlike HC4RSD33, this does not yet close the chart. The next calculation is
the degree-eighteen collision face of (11.39)--(11.41), including its
quartic and lower same-weight tails. The passive-affine boundary
((H_6)_{yy}=0) must be split separately. This is HC4RSD34.

Retain a completely generic homogeneous quartic (c_4). Exact extraction of
the degree-eighteen face gives

\[
\boxed{
[J(c)]_{18}=x^8\left[
(H_6)_{yy}\left(49x^4(c_4)_{zz}-42kx^2P_4+18k^2H_6-6k^3x^5z\right)
-\left(7x^2(P_4)_y-4k(H_6)_y\right)^2
\right].}
\tag{11.42}
\]

This is the next global Schur-square obstruction. Write

\[
c_4=U_4(x,y)+zV_3(x,y)+\frac12z^2D_2(x,y)
+z^3(\ell_xx+\ell_yy)+qz^4.
\tag{11.43}
\]

The (z^2)-coefficient of (11.42) is
(588q x^4(H_6)_{yy}), so (q=0) on the curved chart. Its (z)-coefficient is

\[
6x^4(H_6)_{yy}\bigl(49\ell_xx+49\ell_yy-k^3x\bigr),
\tag{11.44}
\]

and therefore (\ell_x=k^3/49) and (\ell_y=0). The remaining equation is
the binary divisibility identity

\[
(H_6)_{yy}\left(49x^4D_2-42kx^2P_4+18k^2H_6\right)
=\left(7x^2(P_4)_y-4k(H_6)_y\right)^2.
\tag{11.45}
\]

Thus every irreducible factor of ((H_6)_{yy}) occurring to odd order must
divide the resonant form (7x^2(P_4)_y-4k(H_6)_y). This is HC4RSD35. It
turns the degree-eighteen descendant into finitely many UFD root charts.
There is a sharper conclusion on the genuinely moving chart. Use
(z\mapsto z+\alpha x+\beta y) to remove (a_0x^6+a_1x^5y) from (H_6).
Write its remaining coefficients as

\[
H_6=a_2x^4y^2+a_3x^3y^3+a_4x^2y^4+a_5xy^5+a_6y^6.
\tag{11.46}
\]

For (k\ne0), the coefficients of (y^{10}) and (x^2y^8) in (11.45) are

\[
-36a_6^2k^2,\qquad -40a_5^2k^2,
\tag{11.47}
\]

after the first has vanished. Hence (a_6=a_5=0), so
((H_6)_{yy}) is divisible by (x^2). The next immutable coefficient factors
as

\[
-8(-5a_4k+14p_4)(-a_4k+7p_4).
\tag{11.48}
\]

If (a_4\ne0), this leaves two ratios. The second,
(p_4=5a_4k/14), has next coefficient
(-18a_4k(-5a_3k+14p_3)), so it also fixes
(p_3=5a_3k/14). If (a_4=0), then (p_4=0), and the following coefficient is

\[
-9(-2a_3k+7p_3)^2,
\tag{11.49}
\]

which fixes (p_3=2a_3k/7). Thus only three aligned moving packets remain:
the first (a_4\ne0) ratio with free (p_3), the second fully resonant ratio,
and the (a_4=0) rank-drop endpoint. Their degree-seventeen faces have not
yet been classified.

The rest of (11.45) admits a uniform exact decomposition. On the
(a_4\ne0) chart the three coefficients involving (D_2,D_1,D_0) determine
them uniquely. For (p_4=a_4k/7), the penultimate coefficient is a product
(A_LA_R). On (A_L=0), the last coefficient is a nonzero scalar times
(-A_G^2). On (A_R=0), it is a nonzero scalar times

\[
(8a_2a_4-3a_3^2)A_L^2.
\tag{11.50}
\]

Thus this ratio has one generic resonance component and one exceptional
double-root component. For (p_4=5a_4k/14) the same calculation gives
(B_LB_R), followed by (-B_G^2) on (B_L=0) and

\[
(8a_2a_4-3a_3^2)B_L^2
\tag{11.51}
\]

on (B_R=0). This gives two more components. When (a_4=0,a_3\ne0), the
three (D)-coefficients are again unique and the last equation is

\[
-\frac{(6a_2^2k-14a_2p_2+21a_3p_1)^2}{9a_3^2}=0.
\tag{11.52}
\]

Finally, at (a_4=a_3=0) and (a_2\ne0), all three (D)-coefficients are
uniquely determined with no residual equation. Consequently the nonzero-(k)
degree-eighteen locus has exactly six packets: two generic resonance
packets, two repeated-root packets, one (x^3L) packet, and one pure-(x^4)
packet. The (k=0) curved chart and the passive-affine boundary remain
separate at this face.

### 11.6. Closing all six moving packets

The six packets do not survive the next two faces. On the two
(a_4\ne0) ratios, independently of whether the last degree-eighteen
component is generic or lies on (8a_2a_4-3a_3^2=0), the coefficient of
(x^7y^{10}) in degree seventeen is respectively

\[
 \frac{72}{7}a_4^3k^2,
 \qquad
 \frac{18}{7}a_4^3k^2.
\tag{11.53}
\]

Thus all four nonzero-(a_4) packets vanish. On the
(a_4=0,a_3\ne0) packet, the coefficient of (x^{10}y^7) is

\[
 \frac{24}{7}a_3^3k^2,
\tag{11.54}
\]

so that packet also vanishes.

It remains to take (a_4=a_3=0,a_2\ne0). The complete degree-seventeen
face has two decisive (z)-coefficients

\[
 -\frac{12p_1(-3a_2k+7p_2)^2}{a_2},
 \qquad
 -\frac{12(-4a_2k+7p_2)(-3a_2k+7p_2)^2}{7a_2}.
\tag{11.55}
\]

Together with its last binary coefficient, they leave exactly

\[
 p_2=\frac{3a_2k}{7},
 \qquad
 p_2=\frac{4a_2k}{7},\quad p_1=0.
\tag{11.56}
\]

After normalizing (a_2=k=1) and solving every remaining
degree-seventeen equation, the two degree-sixteen descendants contain

\[
 -\frac{24}{49}x^{11}y^4z,
 \qquad
 \frac{256}{1225}x^{10}y^6.
\tag{11.57}
\]

Both are immutable. This proves HC4RSD36: every nonzero-(k) curved packet
in HC4RSD35 is impossible.

### 11.7. The zero-(k) coupled Wronskian and recursive square

Put (k=0). The degree-eighteen binary identity becomes

\[
 (H_6)_{yy}D_2=(P_4)_y^2.
\tag{11.58}
\]

Write (F=(H_6)_{yy}), (G=(P_4)_y), and let (c) be the coefficient of
(z^3) in (c_3). The complete (z)-coefficient of degree seventeen is

\[
 49x^{12}\bigl(D_2G_y-2(D_2)_yG+6cF\bigr).
\tag{11.59}
\]

This extra (6cF) term is essential; omitting it loses projective charts.
When (G\ne0), equations (11.58)--(11.59) have a global classification.
In (K(x)(y)),

\[
 \left(\frac{D_2^2}{G}\right)_y=6c.
\tag{11.60}
\]

Homogeneity makes (D_2^2/G) a linear form (L). Polynomiality of both
(G=D_2^2/L) and (F=D_2^3/L^2) forces (L\mid D_2). Hence, up to nonzero
scalars,

\[
 F=LE^3,\qquad G=LE^2,\qquad D_2=LE,\qquad c=L_y/6
\tag{11.61}
\]

for a second linear form (E). Relative to the marked pure form (x), the
ordered pair ((L,E)) has five projective configurations. Three have
degree-seventeen immutable coefficients (7/36), (21/200), or a nonzero
four-term constant packet. The reversed chart (L=y,E=x) has coefficient
(-1/16) in degree sixteen. The coincident chart (L=E=x) first forces its
last quartic coefficient to zero and then has coefficient (-1). Thus every
(G\ne0) chart vanishes.

It remains to take (G=0). Then (P_4=px^4), (D_2=0), and degree seventeen
fixes

\[
 c_3=T_3(x,y)+zW_2(x,y)+\frac{p^2}{7}xz^2.
\tag{11.62}
\]

For (q=(c_2)_{zz}), the complete next face is the shifted recurrence

\[
 [J]_{16}=x^6\left[
 (H_6)_{yy}(49x^6q-14px^3V_3+6p^2H_6)
 -(7x^3(V_3)_y-3p(H_6)_y)^2
 \right].
\tag{11.63}
\]

On (p\ne0), its first coefficients force (H_6=x^3K_3). Setting
(G_3=7V_3-3pK_3) reduces (11.63) to

\[
 (K_3)_{yy}(49x^3q-2pG_3)=(G_3)_y^2.
\tag{11.64}
\]

If ((G_3)_y=(K_3)_{yy}L_1), differentiation gives

\[
 L_1\bigl(2(K_3)_{yy}(p+(L_1)_y)+(K_3)_{yyy}L_1\bigr)=0.
\tag{11.65}
\]

This globally leaves the zero-(L_1) packet, the proportional repeated-root
packet, and the ((K_3)_{yyy}=0) endpoint. Their normalized next-face
coefficients are

\[
 -\frac{648}{49},\qquad 16,\qquad \frac{12}{7},\qquad \frac{48}{7},
\tag{11.66}
\]

where the last two are the two endpoint subcharts. Hence (p=0).

For (p=0,q\ne0), equation (11.63) is (q(H_6)_{yy}=(V_3)_y^2).
Degree fifteen aligns ((V_3)_y) with (x^2), and degree fourteen has the
immutable coefficient (-4). The global identity

\[
 J\bigl(h(x,z+D(x)y)\bigr)=-D'(x)^2h_z^4
\tag{11.67}
\]

explains this obstruction. If (q=0) and the (x^3z) quartic tail is nonzero,
degree fourteen leaves two algebraic ratios
(49\rho^2-35\rho+7=0); after the degree-thirteen triangular equations,
degree twelve contains (-24x^{11}z). If that quartic tail vanishes, a
nonzero cubic tail collects with all lower linear tails as
(h(x,y)+zD(x,y)), where the leading part of (D) is (x^2). Its exact
(z^2)-coefficient first removes the (y)-term of (D), and its
(z)-coefficient then forces (h_{yy}=0). Linear and affine terminal tails
similarly force a constant cylinder direction. This proves HC4RSD37 and
HC4RSD38: every curved zero-(k) packet is closed.

### 11.8. The passive-affine septic boundary

The sole residue of the curved analysis is ((H_6)_{yy}=0). After removing
unary terms and normalizing its nonzero passive linear form, write

\[
 c_7=x^7,\qquad c_6=x^5L.
\tag{11.68}
\]

For an arbitrary quintic correction, the complete degree-eighteen face is

\[
 [J]_{18}=49x^{12}\det\operatorname{Hess}_{y,z}
 \left(c_5-\frac{2}{7}x^3L^2\right).
\tag{11.69}
\]

Suppose the shifted quintic is curved. In coordinates where it is
(H_5(x,y)+kx^4z), there are two charts. If (L=y), degree seventeen forces
((c_4)_{zz}=0), and degree sixteen is

\[
 x^{12}\left[
 14(H_5)_{yy}(7E_1-k^2x)
 -(7(V_3)_y-3kx^2)^2
 \right],
\tag{11.70}
\]

apart from the coefficient which first removes (z^3) from (c_3). Unique
factorization writes its nonzero solutions as

\[
 (H_5)_{yy}=LE^2,\qquad
 7(V_3)_y-3kx^2=LE,\qquad
 14(7E_1-k^2x)=L.
\tag{11.71}
\]

Four ordered-line charts die in degree fifteen. The coincident pure-root
chart reaches (490k^2-35k+1=0), but its next face contains the incompatible
factor (1449k-125).

If (L=z), degree seventeen instead fixes

\[
 (c_4)_{zz}=\frac67kx^2+\frac6{49}xz.
\tag{11.72}
\]

Degree sixteen removes the last two coefficients of (H_5) and leaves the
ratios (v_3=h_3/3) and (v_3=h_3/7). Their degree-fifteen coefficients are
(-20/21) and (36/7). Thus every curved shifted quintic is impossible.

The remaining two-linear-form tower is

\[
 c_6=x^5L,\qquad
 c_5=\frac27x^3L^2+x^4M.
\tag{11.73}
\]

All four rank patterns of ((L,M)) have zero degree-seventeen face. The next
face is uniformly

\[
 [J]_{16}=49x^{12}\det\operatorname{Hess}_{y,z}(c_4-S),
\tag{11.74}
\]

where

\[
 S=0,\qquad \frac1{49}xy^3,\qquad
 \frac1{49}xy^3+\frac37x^2yz
\tag{11.75}
\]

according to the zero, dependent, or independent rank chart. Passive
singular-Hessian classification reduces (c_4-S) to a binary quartic in one
passive form plus (x^3) times the complementary form. Relative to ((L,M)),
this leaves exactly eight degree-fifteen direction packets. This is the
input to HC4RSD40 and the final conclusion of HC4RSD39.

### 11.9. The common quartic-polar obstruction

Write the nonlinear part of the passive binary quartic as

\[
 Q_4=A_2x^2N^2+A_3xN^3+A_4N^4,
 \qquad (A_2,A_3,A_4)\ne(0,0,0).
\tag{11.76}
\]

The degree-fifteen equations first fix the transverse part of (c_3). In
each of the five packets not killed immediately by a transverse power, the
next three equations are one rank-one system. For example, in the
independent chart with (N=L=y), put

\[
 r=s-7b_0,qquad q=2s-7t_5,qquad p=49t_2-1.
\]

The complete degree-fourteen face is equivalent to

\[
 q^2=-28A_2r,qquad qp=147A_3r,qquad
 p^2=-2058A_4r.                                  \tag{11.77}
\]

The other four charts give (11.77) after nonzero diagonal rescaling. If the
transverse scalar is zero, both polar coordinates vanish. If it is nonzero,
eliminating them gives

\[
 3A_3^2=8A_2A_4.                                  \tag{11.78}
\]

This is not an accidental discriminant: (11.78) says exactly that the
binary quadratic ((Q_4)_{NN}) is a square. Thus all eight packets split
globally into a zero-polar stratum and one square-Hessian resonance.

In the three charts where (N) is transverse to the marked tower direction,
degree fourteen already contains (-16A_4^2), followed after (A_4=0) by
(-12A_3^2). Normalize the remaining (A_2) to one. Two scalar roots remain.
For (N=z), they are (t_1=1/7,2/7); the first has coefficient (12/7) in
degree thirteen, while the second has (-12/2401) in degree twelve. The same
two constants close the dependent chart. For (N=y+z), the two roots are
(q=0,-7), and they die by the same constants. Hence all three transverse
direction packets are empty.

The zero-polar strata of (N=y,L=M=0), (N=y,L=0,M=y), and
(N=y,L=y,M=0) have

\[
 c=h(x,y)+z\bigl(a(x)y+d(x)\bigr),qquad \deg a\le1. \tag{11.79}
\]

For (D=a(x)y+d(x)), direct expansion gives

\[
 [J(c)]_{z^2}=-D\left[-3ya(a')^2+2a^2d''-4aa'd'+d(a')^2\right]. \tag{11.80}
\]

If (a') is nonzero, the coefficient of (y) in the bracket is impossible;
thus (a) is constant. For (a\ne0), equation (11.80) gives (d''=0). An
affine normalization makes (D=y), after which

\[
 [J(h+yz)]_z=-2yh_{xx}.
\]

Consequently (h=xq(y)+r(y)) and the remaining equation is
(-(yq'-q)^2=0), so (q) is linear homogeneous and
(c=r(y)+y(\alpha x+z)) is already a fixed cylinder. If (a=0), the coefficient
linear in (z) is

\[
 d,h_{yy}\bigl(dd''-2(d')^2\bigr).               \tag{11.81}
\]

Passive curvature makes (h_{yy}\ne0), and the leading coefficient of the
parenthesis for (\deg d=n>0) is (-n(n+1)); hence (d) is constant. The
two-variable straight-level lemma then again makes (c) a fixed cylinder.
The other two zero strata instead force (A_4=A_3=A_2=0) in degree twelve,
contradicting (11.76).

It remains only to inspect the square-Hessian resonance. Its degree-thirteen
face always removes the non-(x) root of ((Q_4)_{NN}). The five aligned
charts then have the following immutable descendants:

| packet | degree | coefficient |
|---|---:|---:|
| independent ((N=y)) | 12 | (-3q^2/(49r)) |
| dependent ((N=y)) | 12 | (4b_0u) |
| (L=0,M=y,N=y) | 10 | (-6u^2) |
| (L=M=0,N=y) | 8 | (-u^6/(4b_0^2)) |
| (L=0,M=y,N=z) | 12 | (-147t_5^4/(4r^2)) |

Every displayed denominator and numerator parameter is nonzero on its
localized resonance. This closes all eight packets and proves HC4RSD40.
Combined with the non-pure result HC4RSD33, every scalar degree-seven
leading direction is fixed.

## 12. Revised scalar frontier

The scalar reverse-Schur branch now has the following exact boundary:

- zero-corner quadratic pivots are completely excluded by HC4RSD12;
- nonzero-corner quadratic directions reduce to HC2 or JC2 by
  HC4RSD13--HC4RSD16;
- arbitrary-degree rank-one directions reduce to the same endpoints by
  HC4RSD18;
- genuinely moving leading-rank-three cubic directions are impossible by
  HC4RSD19;
- the residual tangent pencil synchronizes to one ruling by HC4RSD20;
- fixed rulings, all homogeneous border coefficients, and every border
  coefficient through degree four reduce to HC2 or JC2 by
  HC4RSD20--HC4RSD22;
- a degree-five border coefficient with non-pure leading quintic reduces to
  the same endpoints by HC4RSD23;
- the remaining pure-fifth chart is a fixed cylinder by HC4BL5, so every
  degree-five border coefficient closes by HC4RSD24;
- in every degree, a squarefree leading binary border form closes by
  HC4RSD25;
- the generic discriminant stratum with one double root and all remaining
  roots simple closes by HC4RSD26.
- the arbitrary-multiplicity valuation sieve HC4RSD27 exposes every local
  resonance, and all non-pure sextic leading forms close by HC4RSD28.
- in the remaining pure-sixth chart, HC4RSD29 stabilizes a fixed passive
  direction through the quintic correction.
- if that binary quintic correction is curved, HC4RSD30 closes the entire
  potential; only its passive-affine boundary remains.
- if the completed quartic on that boundary is curved, HC4RSD31 closes it;
  only the two-linear-form tower (11.20) remains.
- HC4RSD32 closes every rank pattern of that tower. Consequently every
  scalar degree-six leading direction is fixed.
- HC4RSD33 closes every non-pure scalar degree-seven leading direction.
- HC4RSD34 puts the remaining pure-seventh chart into the exact two-face
  normal form (11.39)--(11.41).
- HC4RSD35 computes its complete degree-eighteen face, forces the quartic
  (z^4,z^3)-tail, leaves the binary square identity (11.45), and aligns the
  nonzero-(k) curvature with (x^2), leaving exactly six packets after the
  final two coefficient factorizations.
- HC4RSD36 closes all six nonzero-(k) packets by the degree-seventeen and
  degree-sixteen coefficients (11.53)--(11.57).
- HC4RSD37 closes the complete ((P_4)_y=0) part of the zero-(k) chart by
  the recursive square (11.63), its global cubic classification, and the
  identities for nonlinear and linear transverse coordinates.
- HC4RSD38 corrects and closes the ((P_4)_y\ne0) part: the cubic (z^3)
  tail upgrades the Wronskian to (11.59), whose complete solution is the
  ordered two-linear-form packet (11.61).
- HC4RSD39 closes every curved shifted-quintic chart on the passive-affine
  boundary and reduces its last two-linear-form tower to eight explicit
  degree-fifteen quartic direction packets.
- HC4RSD40 gives those eight packets one common quartic-polar rank-one
  equation. Its square-Hessian resonance dies in degrees twelve through
  eight, and its only nonempty zero strata are fixed cylinders. Consequently
  every scalar degree-seven leading direction is fixed.

At the endpoint `HC4RSD40`, the next scalar targets were:

1. repeated-root leading forms in degree at least eight;
2. cubic directions whose leading Hessian has rank at most two;
3. higher-degree directions of generic Hessian rank two or three.

The continuation consolidated as
[HC4MR1](HC4_RELATIVE_NILPOTENT_MASTER_REDUCTION.md) consolidates the scalar
and lower-rank continuations. The final regular rank-three negative
maximal-motion sign remains open after HC4MRA1; the complete auxiliary
branch is not closed. The old degree-based list is superseded by that more
precise frontier. This does not
imply unrestricted HC4 or JC2: direct degree-five resonance and nonlinear or
multiple repeated factors remain on the unrestricted direct route, while
matrix pivots and mixed/coisotropic transformations remain separate routes.

## 13. Reproduction

Run:

~~~bash
.venv/bin/python scripts/verify_hc4_higher_degree_pencil_obstructions.py
# committed `HC4RSD17--28` stage artifact only, without symbolic replay:
.venv/bin/python scripts/verify_hc4_higher_degree_pencil_obstructions.py --audit-existing-only
.venv/bin/python scripts/verify_hc4_pure_sextic_collision.py
.venv/bin/python scripts/verify_hc4_pure_sextic_lower_flag.py
.venv/bin/python scripts/verify_hc4_pure_sextic_affine_quartic.py
.venv/bin/python scripts/verify_hc4_pure_sextic_two_linear_tower.py
.venv/bin/python scripts/verify_hc4_nonpure_septic.py
.venv/bin/python scripts/verify_hc4_pure_septic_opening.py
.venv/bin/python scripts/verify_hc4_pure_septic_degree18.py
.venv/bin/python scripts/verify_hc4_pure_septic_moving_closure.py
.venv/bin/python scripts/verify_hc4_pure_septic_kzero.py
.venv/bin/python scripts/verify_hc4_pure_septic_kzero_wronskian.py
.venv/bin/python scripts/verify_hc4_pure_septic_passive_affine.py
.venv/bin/python scripts/verify_hc4_pure_septic_quartic_packets.py
.venv/bin/python scripts/verify_hc4_quintic_bordered_lemma.py
~~~

The maintenance-only mode preserves the exact stage artifact and explicitly
reports that its degree-seven `open_frontier` field is superseded; it neither
recomputes nor rewrites the identities.

The checker verifies the nilpotent characteristic coefficients, the
all-degree rank-one pencil face, the exceptional cubic Hessian and its
moving kernel, all three reduced determinant faces, the invariant
((X,U,Z)) factorization, and the decisive (3b''G_{xx}C_r^2) coefficient. It
also verifies the binary null-pair synchronization, the fixed-cylinder and
cotangent determinants, identity (8.2), the homogeneous Euler face, and all
three degree-three root charts. The quartic continuation checks the four
non-pure binary Schur charts, the double--double fourth-power obstruction,
the six fourth-power-top equations, and every normalized lower correction
face used in Section 9. The quintic continuation checks the simple-root
square, all four repeated-root Schur ideals, the two immutable exceptional
coefficients (10.3)--(10.4), and the lower squares (10.5).
The same simple-root square and Euler deficit verify the all-degree
squarefree induction HC4RSD25.
The double-root valuation coefficients (10.8)--(10.9) verify HC4RSD26.
The companion pure-fifth checker verifies the passive quartic-Hessian face,
both curved unit ideals, every passive-affine projective chart, the complete
transverse alignment radical, and the terminal binary Schur square.
