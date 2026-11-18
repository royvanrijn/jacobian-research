# The exact moment ideal of zero-boundary Gaussian invariant quintics

Use the original Gaussian Wick functional with E[ZW]=E[Y^2]=1 and
all other basic pairings zero. Put S=ZW and consider the full
zero-boundary invariant quintic slice

    P0=S(aY^3+bY^2+cY+d)+S^2(alpha Y+beta).          (M1)

All six coefficients are independent. In the polynomial ring
Q[a,b,c,d,alpha,beta], write mu_n=E[P0^n],
I_N=(mu_1,...,mu_N), I_infinity=(mu_n:n>=1), and let m be the ideal
of all six coefficient parameters. Then, globally,

    I_infinity=I_10, I_9!=I_10,
    m^11 subset I_10, m^10 not subset I_10.          (M2)

The quotient has length 546 and graded Hilbert vector

    (1,5,14,29,49,71,90,100,95,70,22),              (M3)

starting in degree zero and followed by zero. Thus ten is the least
generating moment prefix for the ORDINARY ideal and eleven is the
exact coefficient-nilpotence index. No least set-theoretic cutoff is
asserted. These identities require neither localization nor completion.

Over every characteristic-zero field, P0 is pure if and only if all
six coefficients vanish. If an arbitrary strict one-sided angular
remainder R is added, purity is still equivalent to zero core P0=0.
These pure-locus conclusions are algebraic consequences of the finite
certificate below, independent of every analytic period comparison.

The certificate concerns the pure moment ideal. It is not used to
infer a full finite-image classification or a criterion for a unit
that is merely polynomial. Such conclusions require additional
arguments. Internal exact verification, formal verification and
external review are distinct; the latter two are not claimed here.

## 1. Exact first-moment elimination and original moments

The first moment is mu_1=b+d+2beta. Eliminate it by

    d=-b-2beta, aa=alpha, bb=beta.

The reduced coefficient ring is A=Q[a,b,c,aa,bb], with maximal
coefficient ideal n=(a,b,c,aa,bb), and its original phase is

    P=a S Y^3+b S(Y^2-1)+c S Y+aa S^2Y
                          +bb(S^2-2S).             (M4)

The quotient by mu_1 is a graded polynomial-ring isomorphism with A,
and the image of m is n. Thus no local slice or formal substitution
is involved. Every reduced mu_n is homogeneous of coefficient degree n.
Let g_n be its primitive integer polynomial, dividing by its positive
integer content; this does not change ideals over Q. For example

    g_2=15a^2+2b^2+6ac+c^2+18a aa+6c aa
                                  +12aa^2+4bb^2,  (M5)

and mu_2=2g_2. Put J=(g_2,...,g_10) in A.

Here is an all-order formula binding the generators to the ORIGINAL
Wick functional. For i+j+k+l+q=n, the coefficient of a^i b^j c^k
aa^l bb^q in mu_n is

    n!/(i!j!k!l!q!)
     *sum_(p=0)^j sum_(v=0)^q
       binomial(j,p)binomial(q,v)(-1)^(j-p)(-2)^(q-v)
       *(n+l+v)! G_(3i+k+l+2p),                    (M6)

where G_h=0 for odd h and G_(2h)=(2h-1)!!. Expand
S^n S^l(S-2)^q Y^(3i+k+l)(Y^2-1)^j and apply
E_S[S^h]=h!, E_Y[Y^h]=G_h to prove (M6). This is an identity at every
order, not a conclusion from a finite moment list. The verifier
regenerates g_2,...,g_10 from these finite sums and separately checks
finite fixtures by multiplying the literal original bivariate S,Y
polynomial. The canonical replay has no runtime scratch-file dependency.

## 2. A literal rational Groebner superset

The compressed
[rational certificate](gaussian_zero_boundary_quintic_moment_ideal_certificate.json.gz)
retains 131 explicit homogeneous rational polynomials G, in graded
reverse lexicographic order a,b,c,aa,bb. These literal polynomials,
rather than a CAS success flag, are the finite witness.

After monic normalization the verifier covers all 8,515 S-pairs:
8,126 reduce to zero by exact rational polynomial division, and 389
have coprime leading monomials. For the latter pairs, write
G_i=m_i+r_i and G_j=m_j+r_j. The identity

    S(G_i,G_j)=r_i G_j-r_j G_i

gives the product-criterion standard representation with every product
leading monomial strictly below m_i m_j. Buchberger's criterion proves
that the literal G is a Groebner basis of its own ideal B=(G).
Verification invokes no Groebner-generation oracle.

The same exact division reduces each ORIGINAL g_2,...,g_10 to zero,
proving J subset B. Counting standard monomials of the leading ideal
of G gives (M3), followed by zero in degree eleven. The monomial
bb^10 is among the 22 degree-ten standard monomials. Therefore

    dim_Q(A/J)_r >= dim_Q(A/B)_r=h_r.               (M7)

This is only one direction. In particular J subset B does not by
itself prove that the certificate's G belongs to J. The independent
literal minors supply the opposite dimension inequality.

## 3. Nonzero integer minors and the dimension squeeze

For each coefficient degree r=2,...,11, form the literal homogeneous
Macaulay matrix whose rows are all monomial multiples u g_n with
2<=n<=min(10,r), deg u=r-n, and whose columns are all degree-r
monomials in the five variables. The
[minor certificate](gaussian_zero_boundary_quintic_macaulay_minors.json)
fixes row and column order and selects square minors of sizes

    1,6,21,55,120,230,400,645,979,1365.              (M8)

The determinants modulo the prime 32003 are respectively

    15,28534,20175,8495,27015,18918,30804,
                                  12475,2364,9437. (M9)

The verifier reconstructs every selected integer entry from (M6),
checks the prime, row/column ranges and uniqueness, then computes the
literal determinants. It does not trust a stored rank or RREF flag.
A determinant nonzero modulo a prime is a nonzero integer, so it gives
a rank LOWER bound for the rational matrix. Consequently

    dim_Q(A/J)_r<=h_r, 2<=r<=11.                    (M10)

Together with (M7), and the immediate degree-zero/one dimensions,
this proves equality through degree eleven. In degree eleven the
1365-by-1365 minor has full column rank and determinant 9437 modulo
32003. Thus J_11=A_11 and n^11 subset J. Homogeneity makes every
higher quotient degree zero. The componentwise dimension equalities
then give J=B globally. Since bb^10 is standard for B, it is nonzero
modulo J, proving the nilpotence index is exactly eleven. Summing the
Hilbert vector gives length 546.

In fact the degree-eleven full minor alone already proves the upper
nilpotence bound n^11 subset J and the infinite-tail equality below.
The rational superset is required for the exact Hilbert vector,
length and SHARPNESS of the nilpotence index. This distinguishes what
each part of the certificate establishes; a modular rank deficiency
alone would not prove the rational lower dimension bound.

## 4. The complete infinite tail and least ordinary prefix

For EVERY n>=11, the original reduced mu_n has coefficient degree n,
so it lies in n^11 subset J. The lower moments form J itself.
Therefore (mu_n:n>=2)=(mu_2,...,mu_10) in A, and lifting through
the exact linear first-moment quotient proves I_infinity=I_10 in
the original six-parameter polynomial ring. The first ten moments
and the automatic homogeneous tail meet consecutively: no intervening
moments remain to check.

The rational certificate also retains a linear functional lambda on
the degree-ten coefficient space, with 494 nonzero rational entries.
The verifier checks

    lambda(u g_n)=0 for EVERY monomial u of degree10-n,
                                           2<=n<=9;
    lambda(g_10)=1.                                  (M11)

There are exactly 1,286 displayed-type monomial products, and every
one is checked against the original regenerated generators. Since
the first-nine ideal is homogeneous, these products span its entire
degree-ten part. Thus g_10 does not belong to that ideal. This proves
I_9!=I_10 and least ordinary generating prefix ten, without relying on
how lambda was generated. No least set-theoretic prefix is inferred.

At the zero base phase the original mixed cutoff is C=0. The
[mixed-cutoff finite-jet and certified-tail theorem](MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md)
therefore has moment r-jet budget r. The actual coefficient-nilpotence
index eleven leaves a tail threshold of ten here. This is consistent
with (M2); neither a generic cubic cutoff nor a bespoke asymptotic
tail argument is needed for this coefficient slice.

## 5. Field purity and all coefficient algebras

Over every characteristic-zero field, vanishing of the first ten
moments makes the ideal generated by the six actual coefficients
have eleventh power zero. A field has no nonzero nilpotents, so all
six coefficients vanish. Conversely zero is pure. Thus the field
pure locus is exactly the zero core, proved entirely algebraically.

Let R have strictly positive Z-minus-W weights, or reverse the pair
for a strictly negative remainder. Every term involving R in the
original UNIT moment has nonzero weight, while P0 has weight zero.
Hence Phi_(P0+R)[1]=Phi_P0[1] coefficientwise over every coefficient
ring. Over a field, purity of P0+R is therefore equivalent to P0=0.
That surviving strict phase is safe for every original marker by
angular balance. This proves purity implies safety on this slice;
the certificate does not separately prove that safety implies purity.

All ideal identities in (M2) remain valid after specialization to ANY
commutative Q-algebra T. If its first ten moments vanish, the ideal
J_T of the six actual coefficients has J_T^11=0. It follows that

    P0^11=0 as a spatial polynomial.                 (M12)

Consequently every ORIGINAL marked moment E[Q P0^m] vanishes for
m>=11, regardless of the degree of Q. With a strict remainder R of
minimum angular weight gamma>=1, an original marker of degree D has

    E[Q(P0+R)^m]=0 for m>10+floor(D/gamma).          (M13)

Indeed expansion uses only j<=10 surviving core factors by (M12),
and at most floor(D/gamma) strict factors by angular balance. For
R=0 the cutoff is m>=11. These are nilpotent-coefficient conclusions
over T, not assertions that the coefficients are zero or that every
finite vacuum is primitive over a nonreduced coefficient algebra.

## 6. A sharp all-marker tail in the universal coefficient algebra

The exact phase nilpotence index, and the last surviving marked
power, have explicit witnesses. Work in

    T=A/J=Q[a,b,c,aa,bb]/(g_2,...,g_10).

The coefficient of S^20Y^0 in P^10 is bb^10, nonzero by the standard
monomial certificate. Together with P^11=0 this proves that P itself
has exact nilpotence index eleven, not only its coefficient ideal.

An original marker detects that last power. Define the usual Laguerre
polynomial by its finite formula

    L_20(S)=sum_(i=0)^20 (-1)^i binomial(20,i) S^i/i!,

and let He_j(Y) denote the probabilists' Hermite polynomials. Put

    B(Y)=sum_(j=0)^5 (-1)^j He_(2j)(Y)/(2^j j!),
    Q(S,Y)=L_20(S)B(Y).                              (M14)

For Gamma(1), finite differences of the polynomial (i+1)...(i+h)
give E_S[L_20(S)S^h]=0 for h<20 and 20! for h=20. For ordinary
Gaussian Y, integration by parts gives E[He_(2j)F]=E[F^(2j)]. Hence
for every F of degree at most ten,

    E_Y[B F]=E_Y[exp(-partial_Y^2/2)F]=F(0).        (M15)

The exponential here is a finite operator on that polynomial space;
the last identity follows from the inverse Gaussian heat operator.
The highest S coefficient of P^10 is (aa Y+bb)^10; all other S
coefficients have degree below twenty. Equations (M14)-(M15) yield

    E[Q P^10]=20! bb^10 !=0 in T.                   (M16)

After S=ZW, Q is an ORIGINAL polynomial marker of spatial degree fifty.
Thus the uniform all-marker threshold m>=11 under the first-ten
moment equations is sharp in this universal Q-algebra. This is not a
nonzero field counterexample: every field specialization annihilates
the coefficient ideal. No optimal marker degree is asserted.

## 7. Scope, context and replay

The [quadratic radial affine-coefficient theorem](GAUSSIAN_QUADRATIC_RADIAL_AFFINE_COEFFICIENT_FINITE_IMAGES.md)
gives a separate full finite-image classification and, at zero
boundary, a polynomial-unit criterion using complete relative period
comparison. It is context, not a premise of this algebraic proof.
In particular no analytic comparison theorem, boundary exponent
classification, or universal safety-to-purity statement is used to
deduce (M2)-(M16).

The [axial cubic moment ideal](GAUSSIAN_AXIAL_CUBIC_MOMENT_IDEAL.md)
has the same Hilbert vector, length and nilpotence index. The literal
generators and basis here are different. Equality of these numerical
invariants does not assert an isomorphism of the coefficient algebras.

Run

    python verify_gaussian_zero_boundary_quintic_moment_ideal.py \
      --output gaussian_zero_boundary_quintic_moment_ideal_verification.json

from this directory, with Python and python-flint. The
[portable verifier](verify_gaussian_zero_boundary_quintic_moment_ideal.py)
rebuilds the moment generators, checks exact rational S-pair division,
the literal integer minors, the rational separator and original
finite-source controls. Default replay calls no Singular, Groebner
generator, ideal-membership oracle, RREF or rational reconstruction.
Singular 4.3.2 was used only to generate candidate rational basis and
separator entries. The compressed literal inputs are retained because
they are indispensable witnesses; the generation transcript is not
needed for verification.

The [verification receipt](gaussian_zero_boundary_quintic_moment_ideal_verification.json)
pins the exact proof, checker, certificates and regenerated moment
digest. Corruption controls reject altered source metadata,
basis, separator or minor indices. The infinite tail follows from
homogeneity and the proved eleventh coefficient-ideal power, rather
than extrapolation from a tested prefix. The result is an internal
finite-certificate proof, with no formal-verification or external-review
claim.
