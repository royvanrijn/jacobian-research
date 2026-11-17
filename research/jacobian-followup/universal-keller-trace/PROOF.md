# Trace compression and canonical realization spaces

14 November 2026. The universal Keller map admits a second coordinate reduction for N>=7. Choosing a good trace-zero generator sets the penultimate coefficient to zero and gives a universal map on A^(N-2). The generator argument also removes the earlier characteristic p>N hypothesis over infinite fields: every characteristic other than two is allowed.

Keeping all good generator lines, rather than choosing one, gives a canonical realization correspondence with rational geometric fibres. This proves dense realization targets for each fixed algebra and a weak-approximation statement. It does not identify the different three-dimensional compiler slices up to stable equivalence.

## Inputs and notation

Use the polynomial F_a=(Pi,-B/2,C) of [UACOMP1](../universal-keller-compression/PROOF.md), with t=1+xy and q=t^2*z+y^2*(1+3t). Its determinant is one over Z[1/2]. For pi!=0 its complete fibre is

\[
 \operatorname{Spec}(k[S]/E)[1/E'],\qquad
 E=S+bS^2+\pi S^3+\sum_{j=4}^N a_j\pi^jS^j-c/2.       \tag{1}
\]

The explicit inverse formulas and the rational scaling a_j -> alpha^(j+1)*a_j are proved there over arbitrary coefficient rings with the displayed units. Set a4=a5=w as before.

Define two parameterized polynomial maps, always promoting the remaining coefficients as unchanged coordinates:

- W_N on A^(N-1), N>=5: a4=a5=w, with a6,...,aN free.
- Z_N on A^(N-2), N>=7: additionally a_(N-1)=0, with a6,...,aN except a_(N-1) free.

The second definition leaves N-5 parameters and three moving coordinates. The notation Z_N refers to a map, not a closed subscheme.

Let T_N denote the target open where w, pi, the top coefficient a_N and Disc(E) are units. For Z_N this includes the penultimate-zero equation in its definition; for N=5 the top coefficient is w. The universal finite etale algebra on T_N is B_N=O_(T_N)[S]/E, with E divided by its unit leading coefficient when a monic presentation is needed. Both maps may have other finite fibres outside this chosen chart.

The trace-zero generator technique already appears in the repository's [multiplicity audit](../../archive/non-elliptic/verified/UNIVERSAL_MULTIPLICITY_ADVERSARIAL_AUDIT.md). Its use here removes a promoted coefficient and supplies an explicit realization atlas; no priority claim is made for the underlying primitive-element or trace methods.

## Good generators in every odd characteristic

For a rank-N finite etale algebra A/k and theta in A, write its monic characteristic polynomial as

\[
 P_\theta(T)=\sum_{j=0}^N d_j(\theta)T^j.
\]

Call theta good if it is primitive and d1*d3*d4*d5!=0. The degrees of these four homogeneous coefficient functions in theta are N-1,N-3,N-4,N-5.

**Generator lemma.** For every infinite field k and N>=5, good generators exist. For N>=7, good generators also exist in V_A=ker(Tr_(A/k)). This latter statement does not require N invertible in k.

Here is a direct proof, including small characteristic. After extending to an algebraic closure, identify A with kbar^N and theta with (r1,...,rN). The primitive condition is that all ri are distinct. Each required dj is, up to sign, the elementary symmetric polynomial e_(N-j); none vanishes identically on affine N-space. Thus the good locus is a nonempty open and is defined over k. Nonempty opens in an affine space have k-points when k is infinite.

For trace zero restrict to rN=-(r1+...+r_(N-1)). For every 2<=m<=N-1 the restricted e_m is nonzero: its coefficient of r1^2*r2*...*r_(m-1) is -1. The required indices m=N-1,N-3,N-4,N-5 all lie in this range precisely when N>=7. Every difference ri-rj also restricts to a nonzero linear polynomial. Their product is nonzero because the hyperplane coordinate ring is a domain. The primitive and coefficient conditions therefore define a nonempty open in this hyperplane in every characteristic.

The trace map of a finite etale algebra is surjective, even when Tr(1)=N vanishes: its trace pairing is perfect, so the nonzero vector 1 defines a nonzero linear functional. Consequently V_A is a vector space of dimension N-1 over k, and its good open has a k-point over every infinite k. This avoids dividing by N.

For a good theta put

\[
 \alpha(\theta)=\frac{d_3d_4}{d_1d_5}.
\]

The exact scaling laws are

\[
 d_j(\lambda\theta)=\lambda^{N-j}d_j(\theta),\qquad
 \alpha(\lambda\theta)=\lambda^{-1}\alpha(\theta).       \tag{2}
\]

Hence beta=alpha(theta)*theta depends only on the generator line, and alpha(beta)=1.

The target formulas from UACOMP1, used at translation zero, are

\[
 \begin{aligned}
 w&=\frac{d_3d_4^6}{d_1^2d_5^5},&
 a_j&=\frac{d_jd_3d_4^{j+1}}{d_1^2d_5^{j+1}}\quad(j>=6),\\
 \pi&=\frac{d_1d_5^2}{d_3d_4^2},&
 b&=\frac{d_2d_5}{d_3d_4},&
 c&=-\frac{2d_0d_3d_4}{d_1^2d_5}.
 \end{aligned}                                                   \tag{3}
\]

All are invariant under scaling theta. The corresponding inverse polynomial is

\[
 E(S)=\frac{\alpha}{d_1}P_\theta(S/\alpha).              \tag{4}
\]

It is squarefree of degree N, and S maps to beta in A. Thus (1) gives the whole fibre, not only its geometric points.

For theta trace zero, d_(N-1)=0, so a_(N-1)=0. The extra condition is compatible with all four required nonzero coefficients when N>=7.

## A second dimension reduction

**UATRACE1.** For every N>=7, Z_N is an explicit determinant-one polynomial map over Z[1/2] on A^(N-2). Over every infinite field of characteristic other than two it realizes every rank-N finite etale algebra as a complete fibre. Over every field of characteristic other than two it has geometric degree N, geometric and arithmetic monodromy S_N, and absolute stable atomicity.

The determinant identity is inherited by setting a_(N-1)=0 in the parameterized polynomial identity, before promotion. The generator lemma and (3)-(4) prove fibre universality.

For monodromy, take independent r1,...,r_(N-1) and put rN=-sum ri. Their function field is rational of dimension N-1. Permuting the N roots gives a faithful S_N action in every characteristic for N>=3: no difference ri-rj is identically zero on the trace hyperplane. The polynomial product(T-ri) is separable and its splitting field is this function field. Its coefficient field has splitting degree at most N!, whereas the displayed automorphisms already number N!. Therefore its group is S_N and its coefficient field is the invariant field. Its N-1 nonleading coefficients other than the penultimate coefficient are algebraically independent, since the root field is algebraic over them.

Dividing by the nonzero linear coefficient identifies this generic depressed monic polynomial with the generic coefficient chart h1=1,h_(N-1)=0. Setting pi=h3 and a_j=h_j/pi^j gives (1) with a_(N-1)=0. The same root-free scaling quotient a4/a5 used in UACOMP1 identifies that cover with the Z_N cover after adjoining one independent scaling variable. Purely transcendental extension preserves the Galois group and root action. Thus Z_N has degree N and monodromy S_N, over the field and over its algebraic closure.

The natural S_N action is primitive. The existing intermediate-field argument then proves atomicity and stable atomicity. Its degree-one-factor step is characteristic-free: an etale birational self-map of affine space is an open immersion by Zariski's Main Theorem; factoriality, constant units and normal Hartogs extension make it surjective. The same proof survives extensions of constants and identity stabilizations.

In characteristic zero this improves the programme's universal upper bound to

\[
 \boxed{\operatorname{ed}_k(S_N)\le
       \operatorname{ktdim}_k(S_N)\le N-2\quad(N>=7).}   \tag{5}
\]

The previous N-1 construction remains the proved bound for N=5,6. No optimality of these bounds, or improved rank-four bound, follows.

At N=5 trace zero forces d4=0; at N=6 it forces d5=0. Those are exactly coefficients required by (3). The threshold is structural, not a failed finite search.

### Removing the old characteristic bound

**UAODD1.** W_N, N>=5, is universal over every infinite field of characteristic other than two; p>N is unnecessary. Its degree, monodromy and absolute stable atomicity hold over every such field, finite or infinite.

The good-generator lemma applies without trace zero. Equations (3)-(4) prove universality. The old generic-polynomial and scaling-quotient argument works in every characteristic, and the preceding characteristic-free atomicity argument completes the proof.

The old small-characteristic failure P=T^p-T is a failure of its translation-only chart: changing the primitive generator can enter the good open. It is not an obstruction to universality of W_N over infinite fields.

For finite F_q of odd characteristic, either W_N or Z_N in its indicated rank range is universal whenever

\[
 q>D_N,\qquad D_N=\frac{N(N-1)}2+4N-13.                \tag{6}
\]

Choose any vector-space basis of A. The determinant of the N power-basis columns 1,theta,...,theta^(N-1) has degree N(N-1)/2 in the coordinates of theta. Multiplying by d1*d3*d4*d5 gives a nonzero polynomial of total degree D_N. For Z_N restrict to any basis of the trace kernel; the generator lemma proves this restriction still nonzero. The elementary finite-field polynomial bound gives at most D_N*q^(r-1) zeros in q^r points, with r=N or N-1. Thus q>D_N supplies a good rational generator. Formula (3) applies.

This is a sufficient bound in all odd characteristics, not a sharp threshold. The earlier p>N, q>4N-13 bound for W_N is still useful where it is smaller. Small-field nonmonogenicity remains a real limitation: F_q^N has no generator if N>q.

## A canonical correspondence retaining all choices

**UALIFT1.** For W_N set V_A=A. For Z_N set V_A=ker Tr. Define H_A to be the open subset of P(V_A) consisting of good generator lines. Then H_A is canonical under algebra isomorphisms, geometrically integral and rational over k, of dimension N-1 or N-2 respectively. Moreover

\[
 \boxed{
 H_A\simeq
 \operatorname{Isom}_{T_N\text{-alg}}
       (B_N,A\otimes_k O_{T_N}).
 }                                                               \tag{7}
\]

The right side is a scheme over T_N, not just a set of k-points.

For the forward map take a generator line, form the canonical beta=alpha(theta)*theta, and apply (3). Equation (4) gives the algebra isomorphism S->beta. To construct its inverse, an isomorphism at a target sends S to an element beta in A. Its characteristic polynomial is the monic form of E. Its coefficients satisfy d3*d4=d1*d5, so alpha(beta)=1. The line k*beta is therefore good and recovers the target. Its trace is zero for Z_N. The constructions use only characteristic polynomials, units and scaling, so hold over arbitrary base rings and are inverse as functors.

In particular H_A->T_N is finite etale of degree N!, and is a torsor under the finite etale group scheme Aut_(k-alg)(A). This group scheme is a form of S_N. It must not be replaced by the possibly much smaller group of k-rational automorphisms. The cover need not be a Galois field extension over k. The H_A projective-space description proves rationality of these specific twists directly, without assuming a general rationality theorem for twists.

Over any Z[1/2]-scheme S and rank-N finite etale algebra A, the same construction makes H_A an open in the projective bundle P(V_A), compatible with every base change. The trace pairing makes the trace map locally split surjective, hence its kernel locally free of rank N-1 even on fibres whose characteristic divides N. Each geometric fibre has a nonempty good open by the generator lemma.

Consequently

\[
 S\ \longleftarrow\ H_A\ \longrightarrow\ T_N
\]

is a canonical realization correspondence, with a smooth surjective left arrow and rational geometric fibres. Pulling back the Keller map along the right arrow gives Spec(A) over H_A. Smooth surjectivity supplies etale-local sections. If every residue field of S is infinite, the open-projective-bundle description supplies sections Zariski locally by lifting homogeneous coordinates of a good residue-field point. The same holds at finite residue fields satisfying (6).

This retains the choices in a parameter space instead of choosing them functorially. There is still no section selecting one target naturally for every algebra. Nor does this prove that the distinct three-dimensional compiler slices become stably equivalent: the earlier coarse-boundary obstruction remains separate.

The old GTND1 stable-boundary dimension theorem assumes that all seed coefficients are nonzero. The trace slice a_(N-1)=0 leaves that clean torus. Its dimension count and stable-orbit classification are therefore not imported onto Z_N without a new boundary argument. The monodromy proof of stable atomicity does not use that classification and is unaffected.

### Dense targets and weak approximation

For every fixed A over an infinite k, the targets in T_N(k) with complete fibre isomorphic to A are Zariski dense. Indeed H_A(k) is dense in its nonempty open projective-space chart, and the finite etale map in (7) is surjective. In particular any prescribed nonempty Zariski open target set contains a realization of A. This is stronger than the existence of one realizing target, but does not say that its isomorphism class holds on a Zariski open.

Over a number field k, let V be a finite set of places. Prescribe targets t_v in T_N(k_v) such that the complete local fibre is A tensor k_v. Choose a local isomorphism and hence a lift in H_A(k_v). Weak approximation on projective space gives a point in H_A(k) arbitrarily close to all those lifts; regularity makes its target arbitrarily close to every t_v. Its complete global fibre is A.

This preserves the entire chosen global algebra, not only local factorization types. It uses compatible local lifts; it does not assert strong approximation, integrality, arbitrary local prescriptions for an unrelated global algebra, or a new general theorem about rational points on twists. Earlier repository adelic results treat rational-chart local conditions; here the new ingredient is the exact fixed-algebra lift space for the compressed maps.

## The entire pi-zero fibre and a limit of further depression

**UABOUND1.** For the parameterized map F_a in every characteristic other than two, a target with pi=0 has complete fibre

\[
 \operatorname{Spec}\bigl(k[S]/(S+bS^2-c/2)\bigr)[1/(1+2bS)]
 \ \sqcup\
 \begin{cases}\operatorname{Spec}k,&b\ne0,\\
               \varnothing,&b=0.
 \end{cases}                                                     \tag{8}
\]

Thus its geometric fibre has rank one or three, independently of N and the parameters: rank one if b=0 or 1+2bc=0, rank three otherwise.

The ideals (t) and (q) are comaximal, since q is y^2 modulo t and 1=(1-xy)t+x^2*y^2. Therefore tq=0 splits into the two components t=0 and q=0. On t=0, xy=-1 and q=y^2. The second output gives y=b, so this component is empty when b=0 and otherwise one rational point, with

\[
 x=-1/b,\quad y=b,\quad
 z=b^3c+5b^2+b^3\sum_{j=4}^N(j-2)a_j(-b)^j.
\]

On q=0, t is a unit. The same S=x/t reconstruction as in (1) gives the quadratic algebra in (8); all higher-coefficient terms vanish. If b=0 this is one point. If b!=0 its discriminant is 1+2bc. A nonzero discriminant gives a rank-two etale algebra; a zero discriminant gives a double root entirely deleted by derivative localization. This proves (8). In particular a rank-N fibre, N>=5, must have pi!=0.

This closes a potential boundary escape for one proposed next compression. Impose both a_(N-1)=0 and a_(N-2)=0 in the same compiler. Such a subfamily cannot realize R^N for N>=7. A rank-N fibre cannot lie at pi=0 by (8); on pi!=0 it must have a squarefree degree-N polynomial with N distinct real roots. Its last two nonleading coefficients would be zero. But Newton's identity gives

\[
 0=e_1=\sum r_i,\qquad
 e_2=-\frac12\sum r_i^2<0,
\]

contradicting e2=0. A degree drop or a repeated root also cannot produce a rank-N fibre. Thus simple double depression is globally insufficient for real-split universality, not merely inconvenient for the chosen finite chart. This is not a lower bound against different N-3-dimensional universal constructions.

## Assurance and next gates

The all-rank arguments and the correspondence identification are written proofs. Exact controls check coefficient normalization, trace-zero examples, source reconstruction in quotient algebras, small odd characteristic, the pi-zero components, and the real-split obstruction. They do not constitute formal verification, independent whole-proof replay, external review or a novelty-priority certificate.

The current characteristic-zero bound is N-2 for N>=7 and N-1 for N=5,6. Removing another coordinate requires a different construction or normalization; simply killing the next coefficient is ruled out above. A natural next conceptual gate is a projective-root or stack-valued receiver with effective polynomial landing. The rational lifting correspondence retains all generator choices and does not solve that further descent problem.

