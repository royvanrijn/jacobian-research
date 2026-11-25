# Endomorphism and field-of-definition reduction for the two residual sixfolds

This private note concerns a specific signature. It does not prove an ordinary-prime theorem for the remaining absolutely simple case.

Let E be Q(zeta7) or Q(zeta9), a cyclic CM field of degree 6, with its unique cubic subfield Eplus totally real. Let J be an abelian variety of dimension 6 over a number field, equipped geometrically with an embedding E into End^0(J). Assume its Lie multiplicities in the relevant character order are

    (1,1,2,0,1,1),                                        (1)

up to an automorphism of E and simultaneous conjugation. For E7 the order is1,2,3,4,5,6; for E9 it is1,2,4,5,7,8. In each case conjugate pairs have multiplicities (1,1),(1,1),(2,0). The signature is exactly that of the C7 Jacobian and primitive-order9 factor in the residual construction. No generic endomorphism hypothesis is imposed.

## 1. Precise dichotomy

**Theorem.** If J is NOT absolutely simple with full geometric End^0(J)=E, then its ordinary primes have density one after a finite extension of its number field. Thus any residual sixfold with ordinary-prime density zero after every finite extension must be absolutely simple and have FULL geometric endomorphism field E.

More precisely, every case other than that remaining field case reduces to potential CM or to powers/products of factors of dimension at most two. The converse is NOT asserted: an absolutely simple variety with End^0=E may perfectly well have ordinary primes.

The ordinary inputs are only the known potential-CM and dimension-at-most-two statements. For potential CM a degree-2g field embeds into geometric End^0 of a g-dimensional factor; after defining the action and adjoining its Galois closure, good absolute-degree-one primes are ordinary, giving density one. This is the precise split-prime argument in Fite2005.10185v3 Remark12. Every abelian surface, and every elliptic factor by applying the statement to its square, has density-one ordinary primes after a finite extension by Sawin1506.04784v2. Products and isogenies preserve ordinarity outside finitely many primes. No simple-signature theorem is used for the remaining E-field case.

## 2. Absolutely simple case: the quaternion alternative is impossible

Suppose J is geometrically simple and put D=End^0(J). It is a division algebra. Its faithful action on H_1(J,Q), a12-dimensional rational vector space, gives

    dim_Q D divides12.

The embedded field E makes D a left E-vector space, so6 divides dim_Q D. Hence its dimension is6 or12. In dimension6 it equals E. In dimension12, if commutative, it is a degree12 field, giving potential CM.

Consider the remaining noncommutative dimension12 possibility. The dimension of a central division algebra is a square times the degree of its center. Thus it is a quaternion division algebra over a cubic center F. The commuting E and central F generate a COMMUTATIVE DOMAIN inside D, hence a field. This field justification uses division and must not be transplanted to a matrix algebra.

Since E/Q is Galois, if F is not contained in E then their intersection is Q and their compositum has degree18. It cannot embed into D of dimension12. Consequently F is the unique cubic subfield Eplus. At each complex embedding of F the algebra D becomes M2(C). Every module over this matrix algebra is a sum of standard modules. Restricting a standard module to the quadratic extension E/F gives the two complex E-characters above that real embedding with EQUAL multiplicities. They are a complex-conjugate pair.

Applying this to the Lie representation of J forces equality in every conjugate pair, contradicting the pair(2,0) in(1). The quaternion case is impossible. The only unresolved simple case is therefore D=E; every other simple case is potential CM.

## 3. Several isotypic blocks

Over an algebraic closure, decompose J up to isogeny into its distinct isotypic blocks. The geometric endomorphism algebra is a product of their matrix algebras. Projection of the unital field E to every block is injective, since its kernel is an ideal of E and the identity maps to the identity of that block. Thus each block's rational H_1 dimension is divisible by6, so its abelian dimension is divisible by3.

If there is more than one block, their dimensions must be3+3. Each dimension-three block contains the degree-six field E in its endomorphism algebra, hence is potential CM. This conclusion does not require the block to be simple and does not require E to be central in its matrix endomorphism algebra. Both blocks consequently have ordinary density one after a common finite extension.

## 4. A sole nonsimple block: the B^2 threefold case

Suppose J is isogenous to B^m with B geometrically simple. Its dimension is1,2,3 or6. Dimensions1 and2 are covered by the small-factor ordinary input, and dimension6 is Section2. It remains to consider

    J isogenous to B^2, dim B=3.

Let F=End^0(B), a division algebra of dimension dividing6. Since any noncommutative central-division degree is a square at least4, it cannot occur among divisors of6. Thus F is a FIELD, of degree1,2,3 or6. If its degree is6, B is potential CM, and so is its square for the present ordinary-prime conclusion.

For degree f<6, the E-action is an embedding into M2(F), possibly NONCENTRAL. Let T be the image of a primitive root zeta7 or zeta9. Its F-minimal polynomial divides the separable cyclotomic polynomial Phi7 or Phi9 and has degree at most2 by Cayley--Hamilton. It cannot be linear or a reducible quadratic: a root in F would be a primitive seventh or ninth root, forcing a degree-six cyclotomic field into F, impossible for f<6. Hence the polynomial is irreducible quadratic, and F[T] is an ACTUAL FIELD of dimension2f over Q.

This field contains the image of E, so6 divides2f. For f=1 or2 that is impossible. For f=3 the field has dimension6 and equals E. Therefore the scalar F is contained in E, and must be Eplus, its unique cubic subfield.

Now M2(F) splits over each complex embedding of F, and each Lie module restricts to E/F with equal multiplicities in its two conjugate E-characters, exactly as in Section2. This again contradicts(2,0). Thus the non-potential-CM B^2 case is impossible.

In particular, no unsupported rule that commuting fields in a MATRIX algebra automatically have a field compositum has been used. Irreducibility of the explicit degree-at-most-two minimal polynomial supplies the missing field property. The degree-two possibility for F, whether real or imaginary and regardless of its intersection with E, is excluded by the dimension6 divides2f calculation. This completes the theorem.

## 5. Defining the E-action already forces the base to contain E

**Lemma.** Let J and the E-action have signature(1) at a complex embedding of a number field L. If the action E->End^0(J) is defined over L, then E is contained in L, inside the chosen algebraic closure.

Indeed choose eta=zeta7 or zeta9. Its derivative acts L-linearly on Lie(J/L), so its characteristic polynomial belongs to L[X]. Over an algebraic closure it is

    product_sigma (X-sigma(eta))^(n_sigma).

Exactly ONE multiplicity n_sigma is2. Therefore this characteristic polynomial has exactly one double root, while all its other roots are simple. Every automorphism fixing L must fix that unique double root. Thus sigma0(eta) belongs to L. As E/Q is Galois and sigma0(eta) generates E, this proves E subset L.

Equivalently the Galois action on the multiplicity vector has trivial stabilizer, because it must fix the unique multiplicity-two embedding. The claim does not assume that an abelian variety descending as a variety has all geometric endomorphisms defined over its descent field. It only constrains a field that DEFINES the action. The same argument survives geometric isogeny and a relabeling of the E embeddings.

Consequently, in the remaining simple End^0=E case, ANY number field L whose endomorphisms contain a degree-six field K must have K=E and must contain E. It cannot be Q. At every unramified rational prime below an absolute-degree-one prime of L, the prime splits completely in E, since E/Q is Galois. Such a prime is neither inert nor a product of only two equal-degree prime ideals in E.

## 6. The exact Wang--Zhang2026 transfer fails in the remaining case

The actual published primary is Tian Wang--Pengcheng Zhang, *Ordinary primes for GL2-type abelian varieties and weight2 modular forms*, Mathematika72(2026),e70081, DOI10.1112/mtk.70081. The downloaded published PDF's Theorems1.1--1.3 and adjacent Remarks1.4--1.5 were read, and the entire theorem page was visually checked. No whole32-page proof audit or theorem import is claimed.

Theorem1.1 covers A/Q of dimension q or2q, q prime, with a degree-equal-to-dimension field in End_Q^0. Dimension6 fits2q, but Section5 shows that the necessary End_Q degree-six action CANNOT occur in our remaining signature. This is an intrinsic obstruction, not merely a missing descent construction.

Theorem1.2 treats prime dimension over a quadratic field; dimension6 does not satisfy its dimension hypothesis. It supplies no direct sixfold conclusion.

Theorem1.3 requires absolute simplicity, a degree-g field in End_L^0, and EXPLICITLY excludes Albert typeIV. Our remaining full endomorphism algebra is the CM field E, hence is precisely Albert typeIV (Schur index1). This geometric obstruction cannot be removed by changing the number field. In addition, any field defining that degree-six action must contain E, by Section5, so its unramified absolute-degree-one primes have completely split E, not the inert or two-equal-prime splitting required by Theorem1.3. The ramified cyclotomic conductor prime is totally ramified and is not an inert prime or a splitting into two distinct primes. Neither transfer is available.

The theorem's End_L condition is therefore not interchangeable with geometric End over an algebraic closure, and its NOT-typeIV clause cannot be dropped. Applying the theorem to Eplus instead also fails its degree-g requirement: Eplus has degree3, not6.

## 7. Exact remaining arithmetic question

All special extra-endomorphism and nonsimple cases are now accounted for by the theorem above. The remaining actual candidate sixfolds are ABSOLUTELY SIMPLE, have FULL End^0=Q(zeta7) or Q(zeta9), relative dimension2, and signature(1,1,2,0,1,1).

The three conjugate-pair minima are1,1,0. Thus their minimum CM-type multiplicity sum is2, not1, and the precise simple-signature theorem of Cantoral Farfan et al.4.18 does not apply. The Wang--Zhang hypotheses fail as just proved. No claim about their Mumford--Tate group, generic endomorphism algebra, or ordinary-prime density follows solely from these failures.

To rule out the corresponding hypothetical moment phases by the ordinary-density reverse, one still needs positive-density ordinary primes for these specific simple typeIV sixfolds (or an independent moment argument). This note proves the algebraic reduction and the transfer obstruction; it does not prove that remaining theorem or produce an exceptional phase.
