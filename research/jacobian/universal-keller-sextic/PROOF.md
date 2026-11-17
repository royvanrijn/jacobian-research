# Universal sextics and stable moduli across coefficient boundaries

14 November 2026. A different coprime-weight section extends the N-2 universal construction to rank six. The resulting determinant-one map on A4 realizes every rank-six finite etale algebra over every infinite field of characteristic other than two.

A second result extends the stable-orbit classification of the three-dimensional quadratic-gauge slices across all intermediate coefficient-zero boundaries in characteristic zero. It uses the full two-dimensional support of the intrinsic Fitting divisor instead of its one-dimensional gap pattern. As a consequence, the trace-normalized universal families have N-5 independent stable-class parameters for every N>=6, and each fixed algebra occurs densely among those classes. The universality and stable-inequivalence assertions have different characteristic scopes.

## 1. Polynomial sections from arbitrary coprime weights

Retain the normalized F_a of [UACOMP1](../universal-keller-compression/PROOF.md). Its exact covariance is

\[
 a'_j=\tau^{j+1}a_j,\quad
 (x',y',z')=(\tau x,\tau^{-1}y,\tau^{-2}z),\quad
 F_{a'}(x',y',z')=(\tau^{-2}\Pi,\tau^{-1}b,\tau c).
 \tag{1}
\]

The weights available for a4,a5,a6 are 5,6,7.

Let two retained nonzero coefficients A,B have coprime positive weights r,s. Choose nonnegative integers u,v with

\[
 rv-su=\epsilon,\qquad \epsilon\in\{1,-1\}.
\]

Such a choice exists by Bezout and adding multiples of (r,s). Then

\[
 \tau=(B^u/A^v)^\epsilon,\qquad
 h=(B^r/A^s)^\epsilon
 \tag{2}
\]
give
\[
 \tau^r A=h^u,\qquad \tau^s B=h^v.
 \tag{3}
\]

Conversely A=tau^(-r)*h^u and B=tau^(-s)*h^v. These are Laurent inverse isomorphisms of two tori over Z. Other coefficients transform by their prescribed weights. Setting A=h^u,B=h^v is therefore a root-free quotient section, and nonnegative u,v make its coefficient substitution polynomial on the whole affine h-line. The covariance proves that the generic cover is unchanged after adjoining the independent variable tau.

The previous section a4=a5=h uses (r,s,u,v)=(5,6,1,1), epsilon=-1. After trace zero at rank six, a5=0, but a4 and a6 survive. Use

\[
 (r,s,u,v)=(5,7,2,3),\quad
 \tau=a_6^2/a_4^3,\quad h=a_6^5/a_4^7.
 \tag{4}
\]

This replaces the inappropriate requirement that the surviving weights be consecutive by the actual requirement that they be coprime.

At rank five, trace zero leaves only a5 of weight six. A rational scaling sending every remaining a5 to1 does not exist: in k(a5), the equation tau^6*a5=1 would have valuation 6*v(tau)+1=0 at a5=0. This blocks this normalization of every presentation; it does not rule out a different universal quintic map on A3.

## 2. The universal four-coordinate sextic

Set

\[
 a_4=h^2,\qquad a_5=0,\qquad a_6=h^3
\]
and define
\[
 \boxed{\mathcal S_6(h,x,y,z)=(h,F_{h^2,0,h^3}(x,y,z)).} \tag{5}
\]

It is a polynomial map on A4 over Z[1/2]. The source formulas are

\[
 t=1+xy,\qquad q=t^2z+y^2(1+3t),
\]
\[
 \begin{aligned}
 \Pi&=tq,\\
 B&=y+3xq+4h^2t^2x^2q^4+6h^3t^2x^4q^6,\\
 C&=x(5-3t)-x^3z-2h^2(xq)^4-4h^3(xq)^6,
 \end{aligned}
\]
with F=(Pi,-B/2,C). Its ordinary coordinate degrees are (1,7,41,39). The determinant-one identity follows by substitution into the universal identity, then by the identity diagonal block for h.

For pi!=0 its complete fibre is the derivative-localized algebra of

\[
 E(S)=S+bS^2+\pi S^3+h^2\pi^4S^4+h^3\pi^6S^6-c/2.
 \tag{6}
\]

**UASEXT1.** Over every infinite field of characteristic other than two, every rank-six finite etale algebra is a complete fibre of S6. Over every field of these characteristics, S6 has degree six, geometric and arithmetic monodromy S6 and absolute stable atomicity.

To prove universality, choose a primitive trace-zero theta in the algebra with monic characteristic polynomial
\[
 P_\theta(T)=T^6+d_4T^4+d_3T^3+d_2T^2+d_1T+d_0
\]
and d1*d3*d4!=0. Such generators form a nonempty open in ker Tr. After splitting the algebra, the relevant coefficient functions are elementary symmetric functions e5,e3,e2 on sum(ri)=0; the explicit nonzero-monomial argument in [UATRACE1](../universal-keller-trace/PROOF.md) applies. The trace kernel has dimension five even if the characteristic divides six. An infinite field supplies a rational point of the open.

The original coefficients are a4=d4*d1^3/d3^4, a6=d1^5/d3^6. Equations (4) become

\[
 \tau=d_1/d_4^3,\qquad h=d_1^4/(d_3^2d_4^7).
\]

A realizing target is
\[
 \boxed{
 h=\frac{d_1^4}{d_3^2d_4^7},\quad
 \pi=\frac{d_3d_4^6}{d_1^3},\quad
 b=\frac{d_2d_4^3}{d_1^2},\quad
 c=-\frac{2d_0}{d_4^3}.}                                \tag{7}
\]

Then
\[
 E(S)=d_4^{-3}P_\theta(d_4^3S/d_1).
 \tag{8}
\]
The complete fibre is the desired algebra, with S=(d1/d4^3)*theta. No root extraction occurs.

The generic trace-zero monic sextic has group S6 in every characteristic, by the faithful permutation action on the root hyperplane, as proved in UATRACE1. Dividing by the linear coefficient and using the torus isomorphism (4) identifies that generic cover with (6) after adjoining tau. Hence the sextic map has the same group and degree, geometrically and arithmetically. The characteristic-free primitive-monodromy argument gives absolute stable atomicity.

Over odd F_q, universality is guaranteed by q>25. Indeed the power-basis determinant on ker Tr has degree15, and d1*d3*d4 has degree5+3+2=10. Their product is nonzero by the good-generator argument. The polynomial zero bound yields a good F_q point when q>25. This bound is sufficient, not optimal.

Combining this theorem with UATRACE1 gives, in characteristic zero,

\[
 \boxed{\operatorname{ktdim}_k(S_N)\le N-2
                   \quad\text{for every }N>=6.}       \tag{9}
\]

The rank-five bound remains four, and rank four is unchanged.

### The complete rank-six realization space

Let T6 be the open h*pi*Disc(E)!=0. This is exactly the locus of complete rank-six finite etale fibres over fields: at pi=0 the rank is at most three by UABOUND1; at h=0 with pi!=0 the polynomial degree is three; a repeated sextic root is removed by derivative localization and lowers the rank.

For a rank-six finite etale algebra A, let H_A be the good-generator open in P(ker Tr), now using d1*d3*d4!=0. Since tau(theta)=d1/d4^3 has scaling weight -1, beta=tau(theta)*theta is a canonical generator of the line and satisfies tau(beta)=1. Its characteristic coefficients obey d1=d4^3.

Consequently, exactly as in UALIFT1,
\[
 H_A\simeq\operatorname{Isom}_{T6\text{-alg}}
            (B_6,A\otimes O_{T6}),                     \tag{10}
\]
where B6=O_(T6)[S]/E. This is a finite etale rank720 torsor under the group scheme Aut(A), not a claim that A has720 rational automorphisms. The lift space is an open in projective four-space.

This description is functorial over Z[1/2]-schemes and commutes with base change. For a family of rank-six finite etale algebras over a base S, its good-generator-line space gives a canonical rational realization correspondence. Good lines exist etale-locally on S, and Zariski-locally when the residue fields of S are infinite. These choices refer to S; the finite Isom torsor over T6 is in general only etale-locally split. Over an infinite field, realizing targets for every fixed A are Zariski dense in T6. Compatible weak approximation over number fields follows from the same projective-space argument as UALIFT1.

## 3. Stable classification with arbitrary intermediate zeros

Work in characteristic zero. Normalize a3=1 and write F_a for the degree-N three-dimensional compiler slice, with N>=4, aN!=0 and arbitrary a4,...,a_(N-1), including zeros.

**UASPARSE1.** Two such maps F_a,F_b are stably polynomially left-right equivalent over k if and only if
\[
 \boxed{b_j=s^{j+1}a_j\quad(4<=j<=N)
                   \text{ for some }s\in k^\times.}  \tag{11}
\]
In particular a stable equivalence preserves every zero coefficient. Existence of a stable equivalence is already existence of an ordinary diagonal equivalence. This is an orbit classification, not a classification of all self-equivalences or the vertical stabilizer.

The geometric inputs are already proved for every admissible seed, not only the clean coefficient torus: [RQG2](../../foundations/verified/QUADRATIC_WEIGHTED_STABLE_SEPARATION.md) intrinsically selects the target prime P=0 and the normalized ramified stratum Spec k[P^+-1,r^+-1]. Its hypotheses are g1*g3*gN!=0. [Stable normalization functoriality](../../foundations/verified/STABLE_NORMALIZATION_FUNCTORIALITY.md) transports the stratum and its Fitting ideal through identity stabilization. Vanishing intermediate coefficients does not change those hypotheses, the Newton endpoints, or the normalization.

The relative Fitting ideal has generator
\[
 J_a(P,r)=-1+3Pr^2+
       \sum_{j=4}^N j(j-2)a_jP^jr^{j-1}.              \tag{12}
\]
This follows from Omega_(A/R)=A*dr/(B_r,C_r) and C_r=-r^2*B_r, exactly as in [RQG3](../../foundations/verified/QUADRATIC_GAUGE_STABLE_MODULI.md). Its coefficients retain their multiplicities; r and P are units.

A stable equivalence must induce
\[
 P'\mapsto\beta P,\quad r'\mapsto\alpha P^m r^\epsilon,
 \quad m\in Z,\ \epsilon=\pm1
\]
and therefore
\[
 J_b(\beta P,\alpha P^m r^\epsilon)
             =\lambda P^u r^vJ_a(P,r).               \tag{13}
\]
Units after polynomial stabilization are still constant times monomials in P,r. Because P=0 is an intrinsically selected prime in the affine target, P' maps to a constant times P.

Here is the replacement for the clean-torus gap argument. The extreme r-exponents of J are always 0 and N-1, with unique endpoint monomials -1 and N(N-2)*aN*P^N*r^(N-1). Suppose epsilon=-1. Matching the two extreme r-exponents in (13) forces v=-(N-1). Matching their P-exponents first forces u=-N and then
\[
 N+m(N-1)=-N,\qquad m=-2N/(N-1).
\]
This cannot be an integer for N>=4, since N-1 would have to divide2. Thus epsilon=1 for every intermediate zero pattern.

Now the r-extremes give v=0. The unique r^0 term gives u=0 and lambda=1. The nonzero 3Pr^2 term gives m=0 and beta*alpha^2=1. Comparing the remaining r-exponents gives b_j=alpha^(j+1)*a_j, with simultaneous zeros. Conversely (1) realizes precisely this equivalence. This proves (11).

The argument works over the ground field: all scalar units in the transported K-torus are in k^times, so s lies in k. Extending constants gives the geometric version. It does not extend this orbit theorem to small positive characteristic, where factors such as3 or N(N-2) in (12) may vanish.

At rank six the one-dimensional r-support {0,2,3,5} is reflection symmetric; using only that projection would miss the obstruction. The full (P,r) endpoints still force the impossible m=-12/5. This explains precisely why the earlier all-nonzero proof could not simply be quoted on the trace boundary, and supplies the missing argument.

## 4. Stable trace moduli and fixed-algebra non-descent

**UATRACEMOD1.** For every characteristic-zero field and N>=6, the trace-normalized universal family has an explicit (N-5)-dimensional parameter space of pairwise stably inequivalent, absolutely stably atomic three-dimensional maps. Every fixed rank-N finite etale algebra occurs as a complete fibre on a Zariski-dense set of that parameter space.

For N=6 use Y6=Gm with the coefficients of (5). Under the only possible equivalences (11),
\[
 s^5h^2=(h')^2,\qquad s^7h^3=(h')^3.
\]
The invariant a6^5/a4^7 is h itself, so h=h'. Equivalently raising the two equations to powers3 and2 makes s=1, after which both h^2,h^3 agree.

For N>=7 use the trace family of UATRACE1, with a4=a5=w!=0, a_(N-1)=0, and aN!=0. Its parameter space is
\[
 Y_N\simeq G_m^2\times A^{N-7},
 \qquad\dim Y_N=N-5.
\]
If two members satisfy (11), the equal nonzero a4,a5 coordinates force s=1. All their parameters then agree. Thus the displayed family is separated by stable equivalence even on its intermediate coefficient-zero strata. This is a faithful parameter family, not a component theorem for an ambient moduli stack.

For a fixed algebra A, its canonical lift space H_A maps finitely etale and surjectively onto the corresponding finite target chart T_N. The projection T_N->Y_N is smooth, surjective and has relative dimension3: for every geometric seed the independent pi,b,c coordinates have a nonempty discriminant open. The composite
\[
 H_A\longrightarrow Y_N
\]
is therefore smooth surjective of relative dimension3.

Since H_A is an open in projective (N-2)-space, its k-points are Zariski dense over an infinite field. Their image in Y_N is dense. It follows in particular that each fixed A is realized by infinitely many pairwise stably inequivalent members, and that their stable-class parameters vary in the entire displayed N-5-dimensional family. The unrestricted infinitude theorem already existed in the repository; the new assertion locates this variation inside the trace-normalized family and proves its exact parameter dimension.

Thus the trace normalization does not make this compiler descend through the abstract algebra. Distinct good trace-zero generator presentations of the same A still vary in stable class. The old clean-torus defect was N-4; the trace-family defect is N-5 and remains positive even at rank six. This is not a lower bound for arbitrary universal Keller maps and does not exclude a different compiler or a stack/correspondence receiver retaining the choices.

## Assurance and remaining gates

The statements above have written universal proofs with exact controls. Stable classification imports the existing canonical-boundary and functoriality theorems; it is not independently formalized. The checker verifies root-free sections, the sextic Jacobian, complete fibres including six explicit source points over F27, normalization identities, sparse Fitting transport and failure controls. It is not a replacement for the written all-rank arguments. External review and literature priority remain unclaimed.

The next dimension questions are rank five on A3 and N-3 coordinates in larger ranks. A single weight-six residual coefficient cannot be normalized generically without extracting a root, and the earlier real-split obstruction still rules out simply deleting the next highest coefficient. The stable-orbit gap on the trace boundary is now closed in characteristic zero; effective polynomial landing for a different canonical receiver remains open.

