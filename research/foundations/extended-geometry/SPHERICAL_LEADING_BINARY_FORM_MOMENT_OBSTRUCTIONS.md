# Spherical leading binary forms and nonvanishing power moments

## 1. The functional and the conclusions

Let K be a field of characteristic zero, and put
\[
 A_K=K[x,y,t]/(t^2+2xy-1),\qquad
 \Lambda(f)=\operatorname{CT}_x\frac12\int_{-1}^1
       f\left(x,\frac{1-t^2}{2x},t\right)dt.                 \tag{1}
\]
Integration is the rational polynomial rule, and (a)_n=a(a+1)...(a+n-1) denotes the rising factorial, with (a)_0=1. The degree of a class is the least total degree of a representative. A phase is pure if Lambda(f^m)=0 for every positive integer m, and eventually pure if this holds for every sufficiently large m.

For an exact degree-d phase, d>=1, let H_d be its leading class on the homogeneous cone t^2+2xy=0. Its associated nonzero binary form is
\[
                 h_f(u,v)=H_d(u^2,-v^2/2,uv),             \tag{2}
\]
of degree2d. Roots and multiplicities below are geometric points of P1, including a possible root at infinity. A geometric SL2 change in u,v induces an actual linear trace-preserving sphere automorphism, as proved below.

**Theorem 1 (elementary leading-form obstructions).** The phase f has infinitely many nonzero power moments in either of the following cases.

1. After a geometric SL2 change, its leading binary form is
\[
 \gamma u^A v^B(a u^k+b v^k)^C,\quad \gamma ab\ne0,\quad k\ge1,
 \quad A+B+kC=2d,\quad 0\le A,B,C\le d.                    \tag{3}
\]
2. The leading binary form has a root of exact multiplicity d, with no restriction on the other roots.

In particular, if h_f has at most three distinct roots and every multiplicity is at most d, f has infinitely many nonzero moments. These elementary conclusions require none of the ordinary-reduction inputs in Section6.

**Theorem 2 (additional quintic leading forms).** Let d=5. The same infinite nonvanishing conclusion holds if either

- h_f has at most four distinct roots, each of multiplicity at most five; or
- h_f has exactly five distinct roots, each of multiplicity two.

The proof uses the exact ordinary-reduction theorems stated in Section6, with the full endomorphism hypotheses verified for ALL members of the relevant cyclic-cover families.

Every conclusion allows arbitrary lower sphere-degree terms over K. Thus these phases are neither pure nor eventually pure. The statements concern specified leading configurations; the global primitive degree-five problem remains open. They do not assert that every pure leading form is geometrically unstable, or normalize the full quintic space into a smaller subalgebra.

## 2. The full sphere filtration and trace-preserving changes

Every sphere class has a unique normal form consisting of monomials x^a t^b, a,b>=0, and y^a t^b, a>=1,b>=0. Replacing xy by(1-t^2)/2 gives existence without increasing degree. After substituting y=(1-t^2)/(2x), distinct positive and negative angular weights give uniqueness. The degree-at-most-d space has dimension(d+1)^2; for d=5 it is36-dimensional.

On gr(A_K)=K[x,y,t]/(t^2+2xy), the degree-d normal basis maps under(2), up to nonzero scalars, to exactly the2d+1 binary monomials u^jv^(2d-j). Thus h_f is nonzero and depends only on the leading cone class, not the chosen representative.

Here is a local proof of orthogonal trace invariance over every characteristic-zero extension. Define the formal Gaussian functional G by
\[
             G(e^{ux+vy+wt})=e^{uv+w^2/2}.
\]
Its covariance matrix is the inverse of the matrix of t^2+2xy, which is its own inverse. Hence G is invariant under every linear isometry of that quadratic form. For a homogeneous polynomial V of degree2n, literal monomial integration gives
\[
                    G(V)=2^n(3/2)_n\Lambda(V).             \tag{4}
\]
Both functionals vanish in odd total degree. For(4), both sides vanish on x^a y^b t^c unless a=b,c=2r. In that case
\[
 G(x^ay^at^{2r})=a!(2r)!/(2^rr!),\quad
 \Lambda(x^ay^at^{2r})=
 \frac{2^aa!(2r)!(a+r)!}{r!(2a+2r+1)!},
\]
which verifies the identity. Consequently Lambda is invariant on every homogeneous component and on every polynomial. No complex contour is assumed invariant under an arbitrary complex matrix.

Explicitly, a matrix M=((a,b),(c,e)) with ae-bc=1 induces
\[
 \begin{split}
 x'&=a^2x-2b^2y+2abt,\\
 y'&=-c^2x/2+e^2y-cet,\\
 t'&=acx-2bey+(ae+bc)t.
 \end{split}                                               \tag{5}
\]
Taking the determinant of M[[x,t],[t,-2y]]M^T proves t'^2+2x'y'=t^2+2xy. On the null cone these are precisely the coordinates associated to M(u,v). These transformations preserve the filtration and Lambda.

## 3. A leading coefficient detected at good primes

For any degree-at-most-d phase, write F(x,t) for its literal Laurent polynomial after eliminating y. Every occurring monomial x^i t^j satisfies
\[
                         i+j\le d,\qquad j\ge0.           \tag{6}
\]
For x^a t^b this is immediate. A term y^a t^b contributes x^-a t^(b+2r), 0<=r<=a, for which -a+b+2r<=a+b<=d.

Therefore CT_x F^m has t-degree at most dm. Only the upper face i+j=d can contribute to [x^0t^(dm)]F^m. This face is
\[
 F_\infty(x,t)=H_d(x,-t^2/(2x),t)=t^dL(x/t),\qquad
                         L(z)=z^{-d}h_f(z,1).             \tag{7}
\]
In particular
\[
                   [x^0t^{dm}]F^m=\operatorname{CT}_zL(z)^m.\tag{8}
\]
All lower-degree terms disappear from this coefficient, even when the literal Laurent substitution has individual t-exponents greater than d before angular extraction.

Suppose for the moment that all coefficients lie in a number field. Outside finitely many rational primes they are integral at every prime ideal above p; any specified nonzero coefficient may also be kept a unit. For an odd good prime p=1 modulo d and m=(p-1)/d, (6) gives CT_x F^m of degree at most p-1. Normalized interval integration keeps even t^j with multiplier1/(j+1). The ONLY denominator divisible by p occurs at j=p-1. Thus
\[
 p\Lambda(f^m)\equiv\operatorname{CT}_zL(z)^m
                                      \pmod{\mathfrak p}. \tag{9}
\]
The right side is integral at that prime. If it is nonzero modulo the prime ideal, the original moment is nonzero.

Infinitely many primes1 modulo any fixed positive integer n are sufficient here. This special progression has an elementary cyclotomic proof. For n>1, if only finitely many such primes were available, take a large positive integer r divisible by n and those primes, together with all finite exclusions. Every prime ell dividing Phi_n(r) is new and does not divide n, because Phi_n(0)=1. The polynomial X^n-1 is separable modulo ell, so its distinct cyclotomic factors show that r has exact order n in the residue field. Thus ell=1 modulo n, a contradiction. Phi_n(r)>1 for sufficiently large r. For n=1 use ordinary infinitude of primes.

For arbitrary characteristic-zero coefficients, work after any required finite geometric normalization and form the finitely generated Q-algebra of the transformed coefficients and specified nonzero parameters, with their product inverted. ONE maximal-ideal quotient is a number field. Every rational moment specializes to the corresponding moment, and the required top identity and nonzero factors survive. Infinitely many nonzero specialized moments imply nonzero original moments at those orders. This neither embeds the entire field K in C nor changes the specialization from order to order.

## 4. Proof of the elementary obstructions

Under(3), formula(7) becomes
\[
                  L(z)=\gamma z^{A-d}(az^k+b)^C.
\]
Put s=d-A. The degree relation and B<=d give0<=s<=kC. Choose arbitrarily large good primes p=1 modulo dk and put m=(p-1)/d. Then m is divisible by k and j=sm/k is an integer with0<=j<=Cm<=p-1. Hence
\[
 \operatorname{CT}_zL(z)^m=
                   \gamma^m\binom{Cm}{j}a^jb^{Cm-j}       \tag{10}
\]
is a unit at every sufficiently good prime ideal above p. Every factorial in the binomial has index at most p-1. Formula(9) proves infinite nonvanishing over number fields, and Section3 transfers it to every K with gamma,a,b kept nonzero. Zero exponents C=0 or s=0 cause no exception.

If one leading root has exact multiplicity d, move it to u=0 using(5). Then h_f=u^dP_d(u,v), P_d(0,1)!=0, so L(z)=P_d(z,1) has only nonnegative exponents and nonzero constant c. Thus CT L^m=c^m, and the same prime argument proves Theorem1(2), regardless of the other roots.

Finally, at most three roots can be placed at0,infinity,-1 by a projective linear change; over an algebraic closure rescale its matrix to determinant one. If every multiplicity is at most d, the form becomes gamma*u^A*v^B*(u+v)^C with A+B+C=2d and A,B,C<=d. Apply(3) with k=1. This proves the three-root consequence, including omitted roots with zero exponents.

A root of multiplicity M>d itself supplies a geometric contracting cocharacter for the LEADING binary form: move it to u=0, then (u,v)->(z*u,z^-1*v) gives every monomial positive weight at least2M-2d. This is not a statement about the entire inhomogeneous phase.

## 5. Cyclic curves and the exact Frobenius line

For the remaining quintic cases all relevant root multiplicities a_j are in{1,2,3,4}, summing to ten. Consider the smooth projective normalization C of
\[
                            w^5=h_f(u,v),                \tag{11}
\]
with w of homogeneous weight two. With s distinct roots it is a connected degree-five cover of P1, fully ramified at them, of genus2s-4. The nontrivial deck action gives E=Q(zeta_5) multiplication on its Jacobian; deck invariants have zero H^1 because the quotient is P1.

The normalized cyclic algebra gives the character-i eigensheaf
\[
 O(-2i)\left(\sum_j\lfloor i a_j/5\rfloor\,P_j\right)
                         \subset\pi_*O_C.                \tag{12}
\]
Locally this follows from the integral generator w^i divided by a branch parameter to power floor(i a_j/5). In particular, since every a_j<5, the character-one eigensheaf is EXACTLY O(-2), with one-dimensional H^1.

At a good absolute degree-one prime p=1 modulo5, Frobenius preserves this eigenspace. On an affine chart its Cech generator w/z maps to
\[
                z^{-p}w^p=z^{-p}w h_f(z,1)^m,
                         \qquad m=(p-1)/5.
\]
Projection onto H^1(O(-2)) is multiplication by
\[
                         c_p=[z^{p-1}]h_f(z,1)^m.         \tag{13}
\]
This is on the normalized cover; the character-one summand has no floor correction. If the Jacobian has ordinary reduction, Frobenius is bijective on H^1(O_C), so c_p!=0. Equations(8)--(9) identify exactly this scalar with the prime obstruction to a zero moment.

The eigensheaves also give, up to simultaneous reversal of the character indices,
\[
             \dim H^1(O_C)_{\chi^i}
                       =-1+\sum_j\{ia_j/5\}.             \tag{14}
\]
All our arithmetic arguments may first be made over a number field containing the branch points, covering data and zeta_5. After omitting finitely many primes, distinct branch points retain their multiplicities and tame Kummer normalization gives a smooth proper model. Ordinary reduction is preserved by isogenies and products away from finitely many primes. An ordinary density-one set remains density one after a further finite extension; absolute degree-one primes have density one and then satisfy p=1 modulo5.

## 6. The precise ordinary-reduction inputs

The following classical statements are imported at their exact scopes.

**Surface input.** Every abelian surface has density-one ordinary primes after a finite extension. This follows from Sawin, *Ordinary Primes for Abelian Surfaces*, arXiv1506.04784v2, Theorem1, Corollary2 and Theorem3. The explicit CM cases of Theorem3 become density one after the required endomorphism and CM fields are defined. The same input covers an elliptic factor by applying it to its square.

**Simple-signature input.** Cantoral Farfan--Li--Mantovan--Pries--Tang, arXiv2508.11174v1, Theorem4.18: if the FULL geometric endomorphism algebra is a CM FIELD and its signature is simple, ordinary primes have density one after finite extension. Simple means that some CM type has Lie multiplicities zero except for a single multiplicity one. Mere multiplication by a CM subfield is insufficient. This proof does not apply their Corollary5.1 without its full-endomorphism hypothesis.

**CM input.** Fite, arXiv2005.10185v3, Remark12: a potential-CM variety is ordinary at good absolute degree-one primes whose rational primes split completely in its endomorphism field. Adjoining its Galois closure makes these density one over a finite extension. Potential CM here means that a number field of degree twice the abelian dimension embeds in geometric End^0.

For the source audit, the complete supplied four-page Sawin proof, the exact definition and supplied Theorem4.18 proof plus its stated Case(1) of Theorem4.19, and Fite's complete Remark12 argument were inspected in the actual primary texts. Their preceding deep imports are not re-proved here. The local arguments below VERIFY the endomorphism hypotheses at all special members instead of assuming the generic endomorphism algebra.

## 7. A fourfold endomorphism lemma

Let J be an abelian fourfold with E-action and signature(1,2,0,1) or its simultaneous conjugate. It has density-one ordinary primes after finite extension.

If J is geometrically simple, D=End^0(J) is a division algebra and dim_Q D divides eight, by its faithful module on H_1(J,Q). A noncommutative possibility must have dimension eight: dimension four would equal the embedded E and be commutative. Its quadratic center F must lie in E. Otherwise E and the central F generate an eight-dimensional commutative FIELD inside the division algebra, again equal to D. Thus F=Q(sqrt5), the unique quadratic subfield of E.

At each complex embedding of F, D becomes M2(C). On every such module the two E-embeddings above that real place have equal multiplicities. They are complex conjugate, contradicting the unequal pair(2,0). Therefore D is a field: either E, giving the EXACT simple-signature input, or degree eight, giving potential CM. Both yield the required ordinary primes.

If J is nonsimple, E acts faithfully on each geometric isotypic block, since central projections of a unital field action are injective. Each block's H_1 dimension is divisible by four, so its abelian dimension is even. In dimension four all simple factors then have dimension at most two, and the surface input applies simultaneously.

We also need the signature(0,2,0,2), up to conjugation. Here the E-action is of pure CM type. The complex structure is scalar on each E⊗R component, so every rational E-line in the rank-two E-module H_1 is a Hodge substructure. Weight-one Hodge substructures correspond to abelian subvarieties, giving an isogeny to two abelian surfaces of potential CM. The surface input again suffices. These elementary E-line and isotypic-decomposition facts are also the classical method used in the existing compact-field ordinary arguments.

## 8. Four-root leading decics

The exact-fivefold-root case was already handled without ordinary inputs. With four roots and every multiplicity at most four, the exhaustive partitions of ten are
\[
 (4,4,1,1),\ (4,3,2,1),\ (4,2,2,2),\
 (3,3,3,1),\ (3,3,2,2).                                \tag{15}
\]
The first, second and fifth admit a pairing of complementary exponents summing to five. Every four-point configuration has a projective involution interchanging the two chosen pairs: in coordinates0,infinity,1,lambda it is z->lambda/z. Since this negates branch exponents modulo5, the involution lifts to tau on C, inverting the deck generator g. One constructs the lift by taking a fifth root of a rational function whose divisor is divisible by five, after adjoining a scalar root.

The square tau^2 is a deck transformation. Conjugation by tau both fixes its own square and inverts the deck group; its odd order forces tau^2=1. The two fixed base points are unramified, and inversion followed by translation on each cyclic fibre of size five has exactly one fixed point. Thus tau fixes two points on the genus-four cover. Riemann--Hurwitz gives genus two for C/<tau>, and likewise for C/<g*tau>.

The invariant H^1-subspaces for tau and g*tau have dimension four and intersect trivially: the intersection is g-invariant, but the deck quotient is P1. Their pullbacks therefore give an isogeny from the product of the two genus-two Jacobians onto J(C). The surface input proves ordinary reduction at infinitely many suitable primes.

For the other two partitions, formula(14) gives signatures(1,2,0,1) and(1,0,2,1), respectively. Section7 handles EVERY such Jacobian, including all members with extra endomorphisms. Formula(13), followed by(9), gives infinitely many nonzero moments. This proves the four-root part of Theorem2.

## 9. Five double roots: the all-special sixfold gate

For five roots all of multiplicity two, the normalized cover has genus six and signature
\[
                              (1,3,0,2).                \tag{16}
\]
Let J be its Jacobian.

If J is geometrically simple, D=End^0(J) has dimension dividing twelve. A noncommutative division algebra among these dimensions could only have dimension four or twelve. Dimension four would equal its quartic subfield E. In dimension twelve the center F has degree three, and the commuting fields E and F generate a degree-twelve field inside D, because E/Q is Galois and their degrees are coprime. It would equal D. Both alternatives contradict noncommutativity. Hence D is a field, equal to E or of degree twelve. The EXACT simple-signature input applies in the first case: choose the dimensions1 and0 from the two conjugate pairs. The second case is potential CM.

Now suppose J is nonsimple. Its E-isotypic blocks have even dimensions. If there are several blocks, their dimensions are2+2+2 or4+2. Surface blocks are settled. A dimension-four block has conjugate-pair sums two, its multiplicities are bounded by(16), and the third multiplicity is zero. It therefore has signature(1,2,0,1) or(0,2,0,2), both covered by Section7.

For a sole isotypic block write J isogenous to B^m with B geometrically simple. Factors of dimension one or two are covered by the surface input; the simple dimension-six case was just proved. The only remaining possibility is J~B^2 with dimB=3. Its division endomorphism algebra F has dimension dividing six and hence is a FIELD. The E-action embeds E into M2(F).

This matrix-algebra step needs a separate field justification. Let T be the image of zeta_5. Its F-minimal polynomial divides the separable Phi_5 and has degree at most two by Cayley--Hamilton. It cannot be linear or a product of linear factors: no primitive fifth root belongs to F, since4 does not divide[F:Q], which divides six. Thus the minimal polynomial is irreducible quadratic, and F[T] is an actual quadratic FIELD containing E and the scalar F. Only now may one use the compositum identity to obtain[E:E intersect F]=2. Thus E intersect F=Q(sqrt5).

At each complex embedding of F, this quadratic field acts in the standard two-dimensional M2(C)-module with its two E-characters of equal multiplicity. They are precisely the conjugate E-embeddings above the chosen real embedding of Q(sqrt5). Summing the modules forces equal Lie multiplicities in BOTH E-conjugate pairs. This contradicts(16), whose pairs are(1,2) and(3,0). Thus the remaining B^2 case is impossible.

Every special sixfold has consequently been treated without a generic-End assertion. Its ordinary primes make(13) nonzero, proving the five-double-root part of Theorem2.

## 10. Arbitrary coefficient fields and the signature boundary

For the curve arguments over an arbitrary characteristic-zero K, split the finite root configuration over an extension and spread out the phase coefficients and branch points over a finitely generated Q-algebra. Invert the leading scalar, differences of distinct branch points and required coordinate denominators. ONE number-field specialization preserves the exact root partition and exact sphere degree. After the finite extensions in Sections6--9, omit finitely many bad primes and take ordinary absolute degree-one primes. Equations(9) and(13) give infinitely many nonzero specialized moments, hence infinitely many nonzero original moments. The lower-degree coefficients were unrestricted throughout.

There is a precise finite boundary for the simple-signature part of this method. If a leading decic has s>=4 roots with multiplicities a_j in{1,2,3,4}, sum ten, and u counts those with multiplicity at least three, formula(14) becomes
\[
                      (1,\ 3-u,\ s-5+u,\ s-3).           \tag{17}
\]
For s>=4 a CM type has exactly one multiplicity one and another zero precisely when one of3-u or s-5+u vanishes. Feasible partitions give only the two simple-signature four-root cases in(15), and the five-double-root case. For s>=6, u=3 requires total degree at least s+6>10, while u=5-s is negative. The complementary four-root cases use the separate genus-two factor argument.

This is a boundary of the PRESENT signature argument. It says neither that ordinary primes fail in the other families nor that all remaining purity questions reduce to an ordinary-reduction conjecture. Configurations with at least five roots remain outside the listed exclusions except for the elementary packet and exact-half cases and the five-double-root theorem. Unstable leading forms and the global quintic problem remain unclassified here.

## 11. Classical comparison and exact source scope

The specialization, prime-denominator and binomial arguments are classical tools already used in retained height-box proofs. The ordinary curve and E-action methods have close antecedents in the existing [field-extension manuscript](MATHIEU_CLASSIFICATION_FIELD_AND_COMPACT_EXTENSIONS.md), particularly its degree-five branch-cover and simple-signature passages. The present account proves the full-sphere leading coefficient and the dimension-four/six all-special gates locally; it does not assert worldwide priority or duplicate those older results as new assumptions.

The primary ordinary inputs are [Sawin1506.04784v2](https://arxiv.org/abs/1506.04784v2), [Cantoral Farfan--Li--Mantovan--Pries--Tang2508.11174v1](https://arxiv.org/abs/2508.11174v1), and [Fite2005.10185v3](https://arxiv.org/abs/2005.10185v3), at the exact theorem/proof scopes stated in Section6. In particular the full-End hypothesis of Theorem4.18 is always checked; the generic hypothesis of its Corollary5.1 is never silently discarded.

Pakovich, [arXiv0910.2105v1](https://arxiv.org/abs/0910.2105v1), Theorem4.2 and Section6 Theorems6.1/6.5, concerns fixed rational or Laurent phases of ONE variable and a fixed differential on a contour. Those actual statements and supplied proofs were inspected. They are useful method context, but averaging an additional height variable does not give their fibrewise vanishing hypothesis. No Laurent monodromy classification is imported as a full-sphere theorem.

The exact global sphere purity-recognition theorem remains available, without a practical numerical degree-five cutoff asserted here. The positive conclusions proved above are infinite nonvanishing results for the specified leading forms. They supply necessary restrictions on a potential pure quintic; they do not construct one or settle global degree five.
