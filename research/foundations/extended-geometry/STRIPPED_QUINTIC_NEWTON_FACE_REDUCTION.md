# Newton-face reduction for mixed quintics in the stripped sphere surface

Work over a characteristic-zero field K, with

    A=K[x,y,t]/(t^2+2xy-1),
    Lambda(F)=CT_x (1/2) integral[-1,1] F(x,(1-t^2)/(2x),t) dt.

Integration is the rational polynomial rule. Pure means Lambda(f^m)=0 for every positive integer m; eventual purity means the same for all sufficiently large m. Put q=xt^2. The stripped degree-at-most-five space is

    f=P(y,t)+xt^2 Q(t),   deg P<=5, deg Q<=2.                 (1)

Its dimension is 24. Its even-t subspace has dimension 14. The normal monomials y^j t^r, x^i t^r (j>=0,i>=1) are independent, and their filtration degrees are j+r,i+r: their highest parts remain independent modulo the homogeneous cone relation 2xy+t^2. Thus degree in (1) is the maximum of the displayed monomial degrees, without an assumption about weights of B,C.

## 1. Proven conclusions and the remaining face

**Theorem 1 (mixed-winding necessary conditions).** Suppose the positive part xt^2 Q is nonzero and f has least sphere degree five.

(a) If deg(t^2 Q)=3, f is not eventually pure. In fact it has infinitely many nonzero moments.

(b) If deg(t^2 Q)=2 and f is eventually pure, then Q=B is a nonzero constant, and P has the form

    P=sum_{j+r<=5, j+2r<=7} p_{jr} y^j t^r,                 (2)

with

    A=p_{1,3}!=0, C=p_{3,2}!=0, A^2+2BC!=0.               (3)

The maximal crossing face of its literal Laurent expansion is exactly

    B x t^2 -(A/2)x^{-1}t^5 -(C/8)x^{-3}t^8.              (4)

No exclusion of all phases satisfying (2)--(3) is asserted.

(c) Consequently, in the EVEN 14-dimensional space, every eventually pure genuine quintic with a nonzero positive part must have a nonzero xt^4 coefficient. The condition Q=0 is separate: it includes one-sided pure safe quintics such as y^5.

All excluded cases have infinitely many nonzero moments, not merely a failed finite test. The residual vector space in (2), together with the coefficient B, has dimension 18. Its face contains the odd-t term yt^3, so it has no point in the even-t subspace with A!=0.

Theorem 1 uses the precisely stated surface and E5-fourfold inputs, together with the potential-CM input for a repeated-root degeneration. Section 9 additionally uses the exact threefold theorem stated there. The distinct-root genus-six degree-seven curve arising from (4) is recorded as an unclosed gate; no generic ordinary-prime assertion is used.

## 2. Two-term faces and good-prime isolation

For a general phase P(y,t)+xS(t), let v=deg S and B=[t^v]S!=0. A monomial p_{jr} y^j t^r contributes the upper Laurent monomial

    p_{jr}(-1/2)^j x^{-j}t^{2j+r}.                         (5)

All its other Laurent monomials have strictly smaller t height. For j>=1 set

    sigma(j,r)=(2j+r-v)/(j+1).

If the upper line through (1,v) and (-j,2j+r) is supported by exactly these two monomials and crosses angular zero, its constant-term height per factor is

    H=(j(v+2)+r)/(j+1)=L/(j+1), L=j(v+2)+r.               (6)

At order m=(j+1)k the top angular-constant coefficient is exactly

    binom((j+1)k,k) [p_{jr}(-1/2)^j]^k B^{jk}.             (7)

Take a good prime p=1 mod L, k=(p-1)/L. Since L>=2j>=j+1, m<=p-1. All factorials in (7) are therefore prime to p. Every angular-constant height is at most p-1. The only denominator divisible by p in the interval rule is that of t^{p-1}, so p*Lambda(f^m) reduces to the nonzero coefficient (7). This proves infinitely many nonzero moments over a number field. The sign of the angular coefficient in an exposing weight is irrelevant.

If the maximal height at angular zero is instead a uniquely exposed core monomial a t^h, h>=1, use p=1 mod h and m=(p-1)/h. Its coefficient a^m gives the same conclusion. A face containing only that core and positive angular terms is handled identically: its angular-zero coefficient must use only the core.

Infinitely many primes 1 modulo a fixed positive integer are available by the elementary cyclotomic Euclidean argument: for n divisible by all excluded primes and the modulus, any prime dividing Phi_l(n) has exact order l, away from primes dividing l. The case l=1 is immediate. No finite moment prefix is used.

## 3. Exhausting the possible high faces in (1)

Assume v=2 or 3 and f has degree five. Then P contains a monomial of total degree five. The upper line through (1,v) with largest sigma has intercept H=v+sigma. At least one such monomial has

    sigma=1+(5-v-1)/(j+1)>1.                              (8)

Thus the maximal sigma is >1. Lower terms of each y expansion cannot lie on this face.

For v=3, sigma>1 forces j+r=5, and

    sigma=1+1/(j+1), 0<=j<=5.

These values are all distinct. There is only one P monomial on the maximal line. It is either a core endpoint or the negative endpoint of a two-term crossing face. Section 2 proves (a).

For v=2, sigma>1 forces j+r=4 or 5. The corresponding values are

    1+1/(j+1) for degree 4,  1+2/(j+1) for degree 5.

Values within either row are distinct. A collision between the two rows requires j_5=2j_4+1. The exhaustive possibilities are

    sigma=2:   (j,r)=(0,4),(1,4), H=4;
    sigma=3/2: (j,r)=(1,3),(3,2), H=7/2;
    sigma=4/3: (j,r)=(2,2),(5,0), H=10/3.                  (9)

If only one of the two possible P coefficients is nonzero, Section 2 applies. It remains to treat these three genuine trinomial faces. The treatment depends only on the face; lower terms are arbitrary below its supporting line.

## 4. The exact normalized-cover Frobenius line

The following elementary cyclic-cover calculation fixes the arithmetic input. Let ell>=2 and let h(z) be a polynomial such that EVERY zero of its raw homogeneous degree-2ell completion has multiplicity strictly between 0 and ell; in particular, there is no hidden ellth-power factor. Normalize the connected cover

    w^ell=h(z).

The character-one eigensheaf of its normalized direct image is O(-2): locally the integral generator is w divided by a branch parameter to floor(a/ell), which is zero here. Consequently H^1 in this character is one-dimensional, generated in Cech notation by w/z. At an absolute degree-one good prime p=1 mod ell, put k=(p-1)/ell. Frobenius maps this generator with scalar

    [z^{p-1}] h(z)^k.                                     (10)

This follows by reducing z^{-p} w h(z)^k in H^1(O(-2)). It is on the normalized cover, not on a possibly singular plane equation. If the Jacobian is ordinary, the scalar is nonzero.

All curves and endomorphism factors may be defined after a finite extension. Good ordinary primes of density one then intersect the density-one set of absolute degree-one primes; adjoining zeta_ell makes their rational primes 1 mod ell. Distinct nonzero branch points and tame normalization are preserved away from finitely many primes. Repeated quadratic roots below are treated separately by their normalized genus.

The exact imported ordinary results are: Sawin, arXiv1506.04784v2, Theorem 1/Corollary 2/Theorem 3, gives density-one ordinary primes after finite extension for every abelian surface, and covers an elliptic factor by its square. Cantoral Farfan--Li--Mantovan--Pries--Tang, arXiv2508.11174v1, Theorem 4.18, applies when the FULL geometric End^0 is a CM FIELD and some CM type has multiplicity sum one. A CM subfield alone is insufficient. Fite, arXiv2005.10185v3, Remark 12, supplies ordinary primes after finite extension for potential-CM factors. These are the same precisely scoped inputs used in the reviewed spherical-leading proof; no Corollary 5.1 generic-End substitution occurs.

## 5. The E5 fourfold gate, including extra endomorphisms

We recall its proof to make the new use explicit. If an abelian fourfold J has E=Q(zeta5)-multiplication with signature (1,2,0,1), up to simultaneous conjugation, then J has density-one ordinary primes after a finite extension.

If J is simple, D=End^0(J) is a division algebra of dimension dividing 8. Since E embeds, its dimension is 4 or 8. A noncommutative possibility has dimension 8 and quadratic center F. The commuting E and F generate a FIELD inside the division algebra. If F is not contained in E, that field has dimension 8 and equals D, a contradiction. Hence F=Q(sqrt5). Over each complex embedding of F the algebra D is M2(C); its standard modules give equal multiplicities for the two E embeddings over that real place. This contradicts the conjugate pair (2,0). Thus D is a field: D=E meets the EXACT simple-signature condition (choose multiplicities 1 and 0), while a field of degree 8 gives potential CM.

If J is nonsimple, the E-action on each isotypic block is faithful. Its H_1 dimension is divisible by 4, so each block has even abelian dimension. In dimension 4 the simple factors then have dimension at most 2. The surface input applies simultaneously. This proves the assertion at every special member, rather than assuming the generic endomorphism ring.

## 6. Excluding the sigma=2 and sigma=4/3 faces

For sigma=2 the face is

    B x t^2+D t^4 -(A/2)x^{-1}t^6,  AB!=0.

Set R(z)=Bz^2+Dz-A/2. At order m=k, its top constant term is

    [z^k]R(z)^k,   top t height 4k.

Use the normalized degree-four cover w^4=z^3 R(z), with homogeneous branch exponents (3,1,1,3) when R has distinct roots. It has genus 3. The involution w->-w has a genus-one quotient (normalize v^2=z^3R to (v/z)^2=zR), and the complementary abelian factor has dimension 2. Rational projector endomorphisms give their product isogenous to the Jacobian. Both factors have density-one ordinary primes after finite extension. If R has a double root, the normalization has genus 1. Formula (10) is exactly [z^k]R^k at p=4k+1. Prime isolation then excludes eventual purity.

For sigma=4/3 write its face as

    B x t^2+(A/4)x^{-2}t^6-(C/32)x^{-5}t^10,  BC!=0.

Set R(z)=Bz^2+(A/4)z-C/32, and let Rstar(z)=z^2R(1/z). At order m=3k the top constant term is

    [z^{5k}]R(z)^{3k}=[z^k]Rstar(z)^{3k},

at height 10k. Use the normalized degree-ten cover

    w^10=z^9 Rstar(z)^3.                                  (11)

For distinct roots the homogeneous branch exponents are (9,3,3,5). Riemann--Hurwitz gives genus 7. The eigenspace dimensions for characters i=1,...,9 are

    (1,1,2,0,1,1,0,0,1).

The rational idempotents for deck-character orders 2,5,10 give abelian factors of dimensions 1,2,4. Their images sum isogenously to the Jacobian; the invariant order-one factor is zero. On the order-ten factor, Q(zeta10)=Q(zeta5) has signature (1,2,0,1) in character order 1,3,7,9. Section 5 handles this factor at ALL special members. The other factors are covered by the surface input. If Rstar has a double root, the normalized cover has branch exponents (9,6,5) and genus 2, again covered.

At p=10k+1, formula (10) equals [z^k]Rstar^{3k}. The original moment order is 3k<=p-1 and its constant-term height is p-1. Thus prime isolation again proves infinitely many nonzero moments. This completes both unconditional exceptional-face exclusions.

## 7. The sole remaining sigma=3/2 face and its exact boundary

The middle face in (9) is (4), with R(z)=Bz^2-(A/2)z-C/8. At order m=2k its top constant term is

    [z^{3k}]R(z)^{2k}=[z^k]Rstar(z)^{2k},

at height 7k. The normalized degree-seven cover

    w^7=z^6 Rstar(z)^2                                    (12)

has homogeneous branch exponents (6,2,2,4) and genus 6 when R has distinct roots. Formula (10) is exactly this top coefficient at p=7k+1. If R has a double root, the branch exponents become (6,4,4), the genus is 3, and its Jacobian admits Q(zeta7) of degree 6; it is potential CM. This excludes A^2+2BC=0. Zero A or C already reduces to a two-term face, so (3) is necessary.

For distinct roots, the Q(zeta7) signature from the normalized eigensheaves is

    (1,1,2,0,1,1).                                        (13)

The conjugate pairs have minimum multiplicities 1,1,0. EVERY CM type therefore has multiplicity sum at least 2. The exact simple-signature hypothesis of Cantoral et al. Theorem 4.18 FAILS, even if the full endomorphism field were Q(zeta7). A generic ordinary locus or an assertion about only a CM subfield would not repair this. No all-member ordinary result for (12) is imported here.

The supporting inequality for this face is j+2r<=7 on the original P monomials. Its equality cases inside total degree five are exactly (1,3),(3,2). This proves (2)--(3), including the dimension count 6+5+4+2+1=18 (P rows r=0,1,2,3, plus B). In the even-t space the coefficient A must vanish, so no such residual is possible. This proves Theorem 1(b)--(c).

## 8. Coefficient fields, known cancellation, and scope

For arbitrary characteristic-zero K, take the finite-type Q algebra generated by the phase coefficients and inverses of the finitely many required nonzero face coefficients and discriminants. ONE maximal-ideal specialization has a number-field residue field and preserves the face, degree, and branch pattern. Each of the infinite nonzero specialized moments is then a nonzero original moment. No embedding of all K in C or order-dependent specialization is used.

The case Q=0 is not part of the mixed-winding exclusions. Here moments are those of the height core P(0,t). The one-variable interval-moment theorem (Mueger--Tuset, arXiv1910.12148, Theorem 2.1) implies pure/eventual purity exactly when P(0,t)=0. Then P has strictly negative angular support and is safe against every fixed full-sphere marker, with the usual finite angular-support cutoff.

Known pure phases explain why an arbitrary collinear face cannot be declared nonzero. The orthogonal transport

    q_c=(x-ct-c^2y/2)(t+cy)^2

is pure and safe. Its maximal face is t^3 L(x/t), with

    L(z)=(z-c/2)^4/z^3.

At p=3m+1 the top coefficient is binom(4m,3m)(-c/2)^m, which vanishes modulo p by Lucas's digit criterion. Thus the unit-binomial argument really needs its stated two-term face hypotheses. Fibrewise Laurent purity likewise does not follow from height-averaged purity.

The earlier proper-Laurent fixed-contour results of Pakovich (arXiv0910.2105v1, Theorems 3.4, 6.1, 6.6 and their pole/density hypotheses) do not classify this two-variable height average. The cyclic-cover arguments here are applied only to the explicit selected face coefficients above. No whole S24 purity theorem, unsafe phase construction, general ordinary-prime assertion, or global degree-five conclusion is claimed. In particular, the v=4 leading face and the residual (12) remain outside the proved exclusions.

## 9. A threefold input excludes the lower-degree even perturbation sector

The following further reduction uses one additional, precisely scoped classical input.

**Theorem 2.** In (1), suppose v=4, deg P<=4, and P contains a nonzero negative-angular monomial. If f is eventually pure, then

    P=sum_{j+r<=4, 3j+2r<=9} p_{jr} y^j t^r,
    A=p_{1,3}!=0, C=p_{3,0}!=0, B=[t^2]Q!=0,
    A^2+2BC!=0.                                           (14)

The maximal face is Bxt^4-(A/2)x^{-1}t^5-(C/8)x^{-3}t^6. This residual vector space, including all three coefficients of Q, has dimension 15. It is not claimed to contain no pure phases.

In particular, an EVEN pure genuine quintic in (1) that has BOTH positive and negative angular terms must have B=[t^2]Q!=0, and at least TWO of its three degree-five negative coefficients

    [yt^4]P, [y^3t^2]P, [y^5]P                           (15)

must be nonzero. A phase with only positive or only negative angular support remains a safe possibility and is not excluded.

Here is a proof of every new ordinary and face gate. The complete high-face calculation is still elementary. For v=4 and j+r<=4, the only sigma values attained by two or more P points, at or above the maximal sigma of any negative-angular term, are

    sigma=-1, H=3: (j,r)=(0,3),(1,0);
    sigma=0,  H=4: (j,r)=(0,4),(1,2),(2,0);
    sigma=1/2,H=9/2: (j,r)=(1,3),(3,0).                  (16)

Indeed solve (2j+r-4)/(j+1)=sigma over 0<=j<=4, 0<=r<=4-j. A negative-angular point has sigma>=-1. Values attained at only one P point are excluded by Section 2, including a uniquely exposed core. This finite list is an exhaustion of the displayed integer range, not a numerical purity test.

At sigma=-1, the face is Bxt^4+Ut^3-(A/2)x^{-1}t^2. Write R(z)=Bz^2+Uz-A/2. At order m=k its top constant coefficient is [z^k]R^k, at height 3k. The normalized cover w^3=z^2R(z) has branch exponents (2,1,1,2), genus 2; if R has a double root it has exponents (2,2,2), genus 1. Sawin's surface input therefore supplies ordinary primes p=3k+1 after finite extension. Its character-one scalar is exactly the selected coefficient. This excludes the first row of (16).

At sigma=0, the face is

    Bxt^4+Ut^4-(A/2)x^{-1}t^4+(C/4)x^{-2}t^4,

so its Laurent angular polynomial is

    L(z)=Bz+U-A/(2z)+C/(4z^2).

If C=0, this is the quadratic Laurent case treated by the degree-four cover in Section 6 (or a core-only coefficient when A=0). Suppose C!=0. Put R(z)=Bz^3+Uz^2-(A/2)z+C/4. At order m=k the selected coefficient is [z^{2k}]R^k. It equals the character-one Frobenius scalar at p=4k+1 for the normalized cover

    w^4=z^2R(z).                                         (17)

For three distinct roots of R, all nonzero, the branch exponents are (2,1,1,1,3). The genus is 4. The quotient by w->-w is the elliptic curve obtained by normalizing v^2=z^2R(z), namely (v/z)^2=R(z). The complementary rational projector has abelian dimension 3; the order-four deck action on that factor gives a faithful Q(i) action. This statement requires no simplicity or full-endomorphism assertion.

The exact additional primary input is Fite, *Ordinary primes for some varieties with extra endomorphisms*, arXiv2005.10185v3, Theorem 1: an abelian THREEFOLD whose geometric End^0 contains an imaginary quadratic field has a positive-density set of ordinary primes. It applies to this complementary factor at EVERY special member. Enlarge the number field first to define the factors, include i, and make the elliptic factor ordinary at density one by the surface input. Applying Fite's theorem over that enlarged base gives a positive-density intersection with the elliptic ordinary set and the absolute degree-one primes. Their rational primes are 1 mod 4. Thus the whole cover has infinitely many ordinary primes in the required residue class.

If R has a double and a distinct simple root, the normalized exponents are (2,2,1,3), with genus 2. If it has a triple root, the exponents are (2,3,3), with genus 1. These cases are again covered by the surface input. The endpoints B,C remain nonzero, so no root collides with 0 or infinity. Formula (10) now reads

    [z^{4k}]z^{2k}R^k=[z^{2k}]R^k,

and the original top height is 4k=p-1. Prime isolation excludes eventual purity in the entire second row of (16), not just in the generic threefold case.

The last row of (16) gives R(z)=Bz^2-(A/2)z-C/8, with both negative coefficients nonzero. At order m=2k the top constant term is [z^{3k}]R^{2k}=[z^k]Rstar^{2k}, at height 9k. The corresponding cover is

    w^9=z^8Rstar(z)^2.                                   (18)

For distinct roots its branch exponents are (8,2,2,6), with genus 7. The order-three character factor is an elliptic curve, while the primitive order-nine factor is a sixfold with Q(zeta9)-signature (1,1,2,0,1,1). The latter has the SAME sum-two obstruction to the exact simple-signature condition as (13). No all-member ordinary assertion for that sixfold is used. If R has a double root, its branch exponents become (8,4,6), the normalized genus is 3 and the degree-six field Q(zeta9) acts faithfully; potential CM excludes this degeneration. This accounts for A^2+2BC!=0 in (14).

The supporting inequality in this last row is 2h+i<=9. On a monomial y^jt^r it becomes 3j+2r<=9. Its equality cases inside degree four are exactly (1,3),(3,0), proving (14). There are 12 admissible P monomials and three Q coefficients, so the stated dimension is 15.

For the EVEN consequence, Theorem 1 first forces v=4 whenever the positive part is nonzero and both winding signs occur in a pure genuine quintic. If degP<=4, the last residual would require the odd-t coefficient [yt^3]P!=0; this is impossible in the even space. Hence degP=5. Its only even-t degree-five monomials are the three in (15). If exactly one were present, together with Bxt^4 it would be a two-term supported crossing face at height five per factor; all lower-degree terms lie strictly below. Section 2 excludes that possibility. This proves the assertion without using ordinarity of any generic sixfold.

The primary source scope for the new import is explicit: the actual Theorem 1 statement and its ENTIRE supplied Section 4 proof in arXiv2005.10185v3 were freshly read. That proof invokes its earlier representation-theoretic Propositions 8 and 10, which remain declared imports rather than newly re-proved originals. It does not assume that End^0 is the displayed quadratic field. No conclusion about the remaining degree-five leading face or either distinct-root sixfold (12),(18) is made.
