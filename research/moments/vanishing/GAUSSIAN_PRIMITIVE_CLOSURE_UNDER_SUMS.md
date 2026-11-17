# Finite Gaussian images: independent sums and Gaussian automorphisms

The actual polynomial-time image of an independent sum is finite exactly
when both block images are finite. Its original vacuum is primitive exactly
when both block vacua are primitive. These statements close every established
finite-image primitive-vacuum architecture under independent sums and then
under Gaussian-preserving polynomial automorphisms. The blocks need not be
safe: their polynomial exponents can cancel.

This is an algebraic consequence of the retained finite-image classification,
with the closure principle proposed by the user. The proofs below retain the
original unit, all fixed original markers and their spatial degrees. They do
not identify lattices only after rational localization. General Gaussian
safety-to-purity remains open; independent mathematical review, formal
verification and external review are not claimed.

## 1. Original images and primitive vacua

Work over C with a normalized nonsingular Gaussian Wick functional E; all
identities are formal in t. Let X and Y be independent Gaussian blocks and let

    P(X,Y)=P_1(X)+P_2(Y).

For a phase Q and a fixed original polynomial marker A, set

    Phi^Q_A(t)=E[A exp(tQ)]=sum_(m>=0) E[A Q^m] t^m/m!,
    f_Q=Phi^Q_1,
    V_Q={Phi^Q_A:A is an original polynomial marker},
    O=C[t],  M_Q=O V_Q subset C[[t]].

Every f_Q has f_Q(0)=1. A primitive original vacuum means

    M_Q=f_Q O.                                             (1)

In particular primitivity implies finiteness over O. Safety means that every
Phi^Q_A is polynomial in t, and purity means f_Q=1. The module M_Q and the
literal C-linear fixed-marker image V_Q are distinct definitions.

The [general finite-image theorem](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md#7-every-finite-actual-image-and-its-original-fixed-marker-image)
states that every finite actual M_Q has the unique normalization

    M_Q=exp(psi_Q) O,
    f_Q=p_Q exp(psi_Q),
    psi_Q in tC[t], p_Q in C[t], p_Q(0)=1.                 (2)

For nonconstant Q, also V_Q=M_Q literally. For a constant Q=c, instead
V_Q=C exp(ct), while M_Q=exp(ct)O and p_Q=1. None of the module or vacuum
claims below needs nonconstancy.

## 2. Independent sums preserve and reflect finiteness

For separated original markers A(X), B(Y), independence gives

    Phi^P_(AB)=Phi^(P_1)_A Phi^(P_2)_B.                    (3)

At each t order this is a finite Cauchy product, justified directly by the
binomial expansion of (P_1+P_2)^m. Every original polynomial in X,Y is a
finite C-linear combination of separated monomials. Consequently

    V_P=span_C{uv:u in V_1, v in V_2},
    M_P=span_O{uv:u in M_1, v in M_2},                    (4)

where the second expression means finite sums of products inside C[[t]].
It is an actual multiplication image; no tensor-product injectivity is
asserted.

The exact finiteness equivalence is

    M_P is finite over O
      iff M_1 and M_2 are finite over O.                  (5)

For the forward implication, multiplication by f_2 defines an O-linear
injection

    M_1 -> M_P, h |-> h f_2.                              (6)

Its image lies in M_P by (3)-(4). It is injective because f_2(0)=1, so f_2
is a unit in C[[t]]. A submodule of a finite O-module is finite. Here is
the needed Noetherian argument explicitly. The ring C[t] is a PID. For a
cyclic module O/I its submodules are J/I for ideals J containing I, hence
are finitely generated. Induction on the number of generators, using a
cyclic quotient and the intersection with the remaining generated
submodule, proves the assertion for every finite O-module. Thus f_2 M_1,
and hence M_1, is finite. Interchanging the blocks proves the same for M_2.
This is also the submodule assertion and proof of
[Stacks Project, Lemma 10.51.1](https://stacks.math.columbia.edu/tag/00IJ).

Conversely, if u_1,...,u_r generate M_1 and v_1,...,v_s generate M_2, then
the finitely many products u_i v_j generate M_P by (4). This proves (5)
without any classification or safety hypothesis.

## 3. The exact lattice and its distinguished vacuum

When the equivalent conditions (5) hold, apply (2) to each block:

    M_i=exp(psi_i) O, f_i=p_i exp(psi_i), p_i(0)=1.

Equation (4) gives both inclusions in the actual equality

    M_P=exp(psi_1+psi_2) O,
    f_P=p_1 p_2 exp(psi_1+psi_2).                         (7)

For example, the reverse inclusion uses exp(psi_i) in M_i itself, as a
finite O-linear combination of actual periods. Their product belongs to
M_P by (4). This does not replace an original marker by a rational or
time-dependent marker: it establishes the equality for the stated
O-spans.

By (1)-(2), a finite-image vacuum is primitive exactly when p_Q=1.
Since p_1 p_2=1 in C[t] forces both factors to be constants, and their
values at zero are one, (7) proves

    P has primitive original vacuum
      iff P_1 and P_2 have primitive original vacuum.      (8)

The equivalence can be read without an advance finiteness assumption:
primitivity of P first implies finiteness of M_P and then (5); primitive
blocks likewise have finite images. In contrast, rational localization
would make every nonzero p_i invertible and lose precisely this test.

For finite block images there is also the exact safety criterion

    P is safe iff psi_1+psi_2=0.                          (9)

Indeed (7) is contained in O precisely when its exponential generator is
a polynomial. If exp(eta) is a polynomial for eta in tC[t], differentiation
gives g'=eta' g; polynomial degrees force eta'=0, hence eta=0.
Safety is equivalent to M_P subset O because M_P is the polynomial-time
span of all original periods. Conversely every safe M_P is a submodule
of O and is finite, so (5) and (9) apply. Thus for arbitrary blocks:

    P is safe iff both M_i are finite and psi_1+psi_2=0.   (10)

Under (10), P is pure exactly when p_1=p_2=1. The same statements hold
for any finite number of independent blocks, by induction.

## 4. Gaussian automorphisms transport the full original marked map

Let R=C[X,Y] and let phi be a polynomial algebra automorphism preserving
the SAME Gaussian functional:

    E[phi(A)]=E[A] for every A in R.                       (11)

This hypothesis is algebraic Gaussian invariance; a general polynomial
coordinate change is not assumed to satisfy it. Write Q=phi(P). Applying
(11) coefficientwise and using phi^-1 as well gives the exact identity

    Phi^Q_A=Phi^P_(phi^-1(A)) for every original A.         (12)

For example A phi(P)^m=phi(phi^-1(A)P^m), so (12) holds at every order.
Since phi^-1 bijects original polynomials, it follows that

    V_Q=V_P, M_Q=M_P, f_Q=f_P,
    ker Phi^Q=phi(ker Phi^P).                             (13)

The kernels here are C-linear subspaces, not necessarily ideals. In
particular finiteness, primitivity, safety and purity are each preserved
and reflected. In the finite case the normalized psi and p are unchanged.

The degree information is also retained. Put r=deg(phi) and
q=deg(phi^-1), where the degree of an automorphism is the maximum total
degree of its coordinate images. Both are positive integers. For the
literal marker-degree spaces

    V_P(d)={Phi^P_A:deg A<=d},

equation (12) gives

    V_Q(d)=Phi^P(phi^-1(R_(degree<=d))),
    V_P(floor(d/r)) subset V_Q(d) subset V_P(qd).          (14)

The first inclusion follows by taking A=phi(B) for every deg B<=floor(d/r).
These statements concern fixed original markers; there are no denominators
and no changes to the Gaussian vacuum.

This general transport identity specializes to the explicit heat-transform
kernels, marker lifts and sharp affine-shear filtration in the
[existing automorphism formulas](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md#8-exact-original-kernels-for-gaussian-automorphism-orbits).
For example the established Gaussian-preserving shear

    phi(Z)=Z, phi(Y)=Y-Zg(Z),
    phi(W)=W+g(Z)Y-Zg(Z)^2/2

has E[ZW]=E[Y^2]=1 and can be extended by the identity on independent
spectator blocks. Its Gaussian invariance is proved there by a locally
nilpotent divergence-free vector field preserving 2ZW+Y^2.

## 5. Closing the established primitive-vacuum architectures

Suppose each starting class has the following proved property: whenever
a phase in that class has finite actual image, its original vacuum is
primitive. The class need not consist only of finite-image phases.

Take finitely many such phases on independent blocks and apply a
Gaussian-preserving polynomial automorphism to their sum. If the resulting
phase is safe, (13) makes the independent sum safe, hence finite. Equation
(5) makes every block image finite. The starting theorems make their
vacua primitive, and (8), followed by (13), makes the resulting vacuum
primitive. For a safe phase, (2) has psi=0, so a primitive vacuum is
f=1. Therefore safety implies purity throughout this closure.

More strongly, every finite actual image in this closure has primitive
vacuum, whether or not it is safe. Once finite block images are known,
their transformed sum is safe exactly when sum_i psi_i=0.

One substantial starting class is the
[multibase polynomial quadratic-fibre theorem](GAUSSIAN_MULTIBASE_QUADRATIC_FIBRE_PRIMITIVITY.md#5-every-finite-actual-image-has-primitive-vacuum):

    P_i=W_i^T H_i(Z_i)+Y_i^T A_i(Z_i)Y_i
          +B_i(Z_i)^T Y_i+G_i(Z_i).

Here H_i depends only on Z_i. Its exact finite-image criterion is that
the source branch r_i=tH_i(r_i) is polynomial and

    det(I-t DH_i(r_i))=det(I-2t A_i(r_i))=1.

In that case p_i=1 and its exponent is the polynomial

    psi_i=tG_i(r_i)
       +t^2 B_i(r_i)^T (I-2t A_i(r_i))^-1 B_i(r_i)/2.      (15)

The closure therefore has a concrete branch test followed by the exact
cancellation test sum_i psi_i=0. The same closure principle applies to the
proved mixed zero-root and weakly one-sided quadratic zero-weight classes
in [Sections 9-10 of that note](GAUSSIAN_MULTIBASE_QUADRATIC_FIBRE_PRIMITIVITY.md#9-a-separate-mixed-source-extension-with-zero-root),
and to the precisely covered finite-vacuum branches of the
[W-affine classification](GAUSSIAN_W_AFFINE_FINITE_IMAGE_OBSTRUCTIONS.md).
Each application retains the hypotheses and assurance limits of its
starting theorem.

The practical test for a coupled phase is consequently to exhibit a
Gaussian-preserving phi for which phi^-1(Q) separates into classified
independent blocks, verify their finite-image criteria, and add their
normalized exponents. This exact closure principle makes no assertion
that every phase has such a decomposition. The existing
[divergent-unit obstruction to universal auxiliary encoding](GAUSSIAN_MULTIBASE_QUADRATIC_FIBRE_PRIMITIVITY.md#8-spectators-classical-mechanisms-and-scope)
remains: the unit for (ZW)^2 has radius zero, whereas all original periods
of the multibase architecture converge near zero.

## 6. Exact degree-bounded cutoffs and exponent cancellation

For a safe phase define the integer

    K_P(d)=max{deg_t Phi^P_A : deg A<=d, Phi^P_A!=0},
    d in Z_(>=0).                                        (16)

This maximum exists: the space of markers of degree at most d is finite
dimensional and every image is polynomial. It is nonempty because the
unit marker has nonzero period. Thus K_P(d) is the last possible nonzero
moment index among these markers: E[A P^m]=0 for every deg A<=d and
m>K_P(d), with equality attained by some marker at m=K_P(d). The convention
uses m>K, rather than m>=K.

For independent safe blocks, the exact formula is

    K_(P_1+P_2)(d)
      =max_(d_1+d_2<=d) (K_(P_1)(d_1)+K_(P_2)(d_2)),      (17)

where d_1,d_2 are nonnegative integers. For the lower bound choose markers
attaining the two block maxima. Their product has spatial degree at most
d_1+d_2 and, by (3), period degree exactly the sum, since nonzero polynomial
leading coefficients multiply. For the upper bound expand any deg A<=d
as separated monomials X^alpha Y^beta. Each product period has degree
at most K_1(|alpha|)+K_2(|beta|), and |alpha|+|beta|<=d. A sum cannot
increase the largest degree. This proves both directions in (17).
Induction gives the corresponding formula for any finite number of blocks.

For Q=phi(P), the exact spaces (14) imply

    K_P(floor(d/r))<=K_Q(d)<=K_P(qd),
    K_P(d)<=K_Q(rd).                                      (18)

General automorphisms give these bounds; an exact transformed cutoff can
be recovered from (12) when its original-degree marker images are known.
No equality of degree filtrations is assumed.

The blocks in (7)-(10) need not be safe. For independent standard
one-variable Gaussians X,Y, take P_1=X and P_2=iY. Their actual images
and unit periods are

    M_1=exp(t^2/2)O, f_1=exp(t^2/2),
    M_2=exp(-t^2/2)O, f_2=exp(-t^2/2).

Both vacua are primitive and both blocks are unsafe. Their sum has
M_(X+iY)=O and f_(X+iY)=1, hence is safe and pure. More explicitly, put

    H_a(z)=exp(-z^2/2) (d/dz)^a exp(z^2/2).

Then Phi^(X+iY)_(X^a Y^b)=H_a(t)H_b(it), a polynomial of exact degree
a+b. Thus K_(X+iY)(d)=d. Formula (17) was stated for safe blocks; this
example instead uses their cancelling exponential factors. Its proof
also applies to finite blocks after replacing every Phi^P_A by its
normalized polynomial exp(-psi_P)Phi^P_A and defining the corresponding
normalized degree maxima.

The unrestricted results above use the written factorization, module,
classification and transport arguments. No finite experiment is used to
infer them, and no mathematical computation is required by this proof.

