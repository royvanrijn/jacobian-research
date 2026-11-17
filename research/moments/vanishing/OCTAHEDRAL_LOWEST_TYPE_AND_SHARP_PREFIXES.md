# The full lowest octahedral type: sharp pure and ideal prefixes

15 November 2026. Internal proof with exact universal coefficient
certificates. This executes a complete finite coefficient-space calculation
on SL_2/2O, following the previously theoretical
[Haar-period decision framework](ENTIRE_HAAR_PERIODS_AND_EXCEPTIONAL_DECISION.md).
The [checker](verify_octahedral_minimal_type.py),
[frozen identities](octahedral_moment_relations.json),
[receipt](octahedral_minimal_verification.json) and
[source audit](octahedral_minimal_source_audit.json) delimit its evidence.
Classical invariant-ring inputs remain explicit. No novelty, independent
specialist review or formal verification is claimed.

**Result.** On SL_2/2O every pure-null function of matrix-entry degree at
most 11 lies in the SL_2-nullcone and is eventually marked-null for every
fixed regular marker. The first seven pure moments suffice, sharply. The
full ordinary pure-moment ideal has least generating prefix ten. The
optimal uniform zero-offset marked bound on this pure-null space is

    h(QP^m)=0 whenever 2m>deg_matrix_entries(Q).            (1)

Thus a Mathieu--Zhao counterexample on the octahedral model must have
matrix-entry degree at least 12, improving the retained bound 8. The
all-degree octahedral question and the icosahedral lowest-type question
remain open. This is a result about all multipliers in the whole coordinate
algebra for P in a specified full vector space, not an MZ assertion for
a finite-dimensional subalgebra.

## 1. The full coefficient space, without a support restriction

Use coordinates g=[[x,y],[z,w]], xw-yz=1, and the compact Klein form

    p=x^8+14x^4 y^4+y^8,
    F=z partial_x+w partial_y,
    p_j=F^j p/j!       (j=0,...,8).

Right binary-octahedral invariance of p was certified in the
[Haar-shear proof](HAAR_SHEARS_AND_GENERAL_NULLCONE_CLASSIFICATION.md).
F commutes with right translations, so every p_j is invariant. They form
one full irreducible left type V_8, with substitution weights 8-2j.

The exact character calculation gives, through highest weight 12,

    dim V_n^(2O)=1 for n=0,8,12 and 0 for all other n<=12. (2)

Algebraic Peter--Weyl and the matrix-degree filtration consequently give

    A_(degree<=11)=C plus span(p_0,...,p_8).               (3)

Indeed degree<=d has only types V_n with n<=d; each such type has a
homogeneous matrix-coefficient representative of degree n. Multiplicities
are the fixed dimensions in (2). Constants are removed by the first
moment. No coefficients of V_8 are set to zero: write

    P_a=sum_(j=0)^8 a_j p_j,
    f_a(X,Y)=sum_(j=0)^8 binom(8,j) a_j X^(8-j)Y^j.       (4)

The latter is the binary octavic model of the parameter action. At the
identity the function is a_0+14a_4+a_8. Polarization verifies that
P_a(g) is evaluation of this same linear functional on the transformed
binary form; equivalently one can check both unipotent generator actions.
The checker verifies the ladder, weights, determinant-harmonicity, and
identity values of all nine matrix-coefficient basis functions.

## 2. Compute exact universal moment polynomials by Reynolds projection

Let H_m(a)=h(P_a^m). On the parameter polynomial ring put

    U=sum_(i=1)^8 i a_(i-1) partial_(a_i),
    D=sum_(i=0)^7 (8-i) a_(i+1) partial_(a_i),
    T=sum_(i=0)^8 (8-2i) a_i partial_(a_i).

They satisfy [U,D]=T, [T,U]=2U and [T,D]=-2D. An initial torus average of
(a_0+14a_4+a_8)^m leaves the weight-zero polynomial

    B_m=sum_(r=0)^floor(m/2)
          m! 14^(m-2r)/(r!^2 (m-2r)!)
          a_0^r a_4^(m-2r) a_8^r.

The remaining Reynolds projection is the finite formula

    H_m=sum_(k=0)^(4m) (-1)^k/(k!(k+1)!) U^k D^k B_m.    (5)

Here is a proof of the normalization. Every irreducible summand of
Sym^m(V_8) containing weight zero is V_(2l), with 0<=l<=4m. On its
zero-weight line U^k D^k acts by (l+k)!/(l-k)! for k<=l and by zero
otherwise. For l>0 the sum of these scalars in (5) is

    (1/l) sum_(k=0)^l (-1)^k binom(l,k) binom(l+k,l-1)=0,

the l-th finite difference of a polynomial of degree l-1. For l=0 the
sum is one. Complete reducibility proves (5). This is an exact projection
for every coefficient, not a fit to finitely many parameter samples.

The standard-library checker expands (5) for m=2,...,10 over Q and checks
U H_m=D H_m=0. It also detects the wrong denominator (k!)^2: that false
projector preserves the noninvariant zero-weight vector 14a_4 at m=1.

## 3. Match the entire small invariant spaces, including degree ten

For reviewable coordinates define primitive transvectants as follows.
For binary forms A,B use the unnormalized operation

    (A,B)_r=sum_(k=0)^r (-1)^k binom(r,k)
        (partial_X^(r-k) partial_Y^k A)
        (partial_X^k partial_Y^(r-k) B).

After every such operation divide all coefficients by their positive
rational content, leaving primitive integer coefficients and preserving
the sign. Starting with f=f_a, put

    g=(f,f)_4, k=(f,f)_6, h=(k,k)_2, m=(f,k)_4,
    n=(f,h)_4, r=(g,k)_4, s=(g,h)_4,

and, with the same primitive convention,

    J_2=(f,f)_8, J_3=(f,g)_8,
    J_4=(k,k)_4, J_5=(m,k)_4, J_6=(k,h)_4,
    J_7=(m,h)_4, J_8=(r,h)_4, J_9=(n,h)_4, J_10=(s,h)_4.

These agree with the classical octavic transvectants up to nonzero
rational scalars. Their definitions are recorded in
[Shaska, equations (9)--(10)](https://arxiv.org/abs/1209.1237v1).
We do not import that paper's syzygy corrections, moduli equations, or
omitted elimination details in its nullcone proof.

The [frozen certificate](octahedral_moment_relations.json) gives the exact
identity H_d=gamma_d J_d plus a polynomial in earlier J_i, for every
d=2,...,10. All monomial coefficients are checked. The diagonal values are:

| d | gamma_d |
|---:|---:|
| 2 | 16/15 |
| 3 | 6272/715 |
| 4 | -200704/182325 |
| 5 | -16859136/1154725 |
| 6 | -89915392/47805615 |
| 7 | -4405854208/462120945 |
| 8 | 40282095616/71628746475 |
| 9 | 31581162962944/2650263619575 |
| 10 | 288742061375488/108660808402575 |

For example H_2=(16/15)J_2, H_3=(6272/715)J_3 and

    H_4=(580864 J_2^2-200704 J_4)/182325,
    H_5=(72052736 J_2 J_3-16859136 J_5)/1154725.

The checker also establishes that all monomials in J_2,...,J_10 of each
weighted degree d<=10 are linearly independent. Their counts at d=2,...,10
are 1,1,2,2,4,4,7,8,12. These are the **full** invariant dimensions:
count monomials in Sym^d(V_8) of weights 0 and2 and subtract the latter
count from the former. Every nontrivial even irreducible contributes one
to both counts; only a trivial summand contributes to their difference.
This independently checks the dimensions rather than inferring them from
the existence of an identity.

Thus these products span every invariant space of degree at most ten.
Because the diagonal values are nonzero, the same is true for products of
the H_i, and

    (H_2,...,H_d)=(J_2,...,J_d) for 2<=d<=10.              (6)

Two classical facts are imported from binary invariant theory: the
octavic ring is generated by nine basic invariants in degrees 2 through10,
and it has a homogeneous system of parameters in degrees 2 through7.
See [Brouwer--Draisma--Popoviciu, Section 6.6 and Proposition 2.1](https://arxiv.org/abs/1404.5722v1).
The cited source states these as classical inputs from Shioda. Our checker
does not independently prove those all-degree invariant-ring facts.

Since our degree spaces through ten are complete, inducting on degree
shows that J_2,...,J_10 (and H_2,...,H_10) generate the invariant ring.
The same low-degree completeness implies that every element of a
classical degree-(2,...,7) system of parameters belongs to (H_2,...,H_7).
Their common zero set is therefore the octavic nullcone. This avoids
requiring any particular normalization of the classical system of parameters.

## 4. The sharp pure prefix and complete marked safety

By (6) and the classical parameter criterion,

    H_2=...=H_7=0
      iff f_a has a root of multiplicity at least five
      iff P_a is SL_2-unstable.                           (7)

The binary-form nullcone characterization follows from Hilbert--Mumford;
it is also stated in Section 2 of the cited Brouwer--Draisma--Popoviciu
paper. Every unstable element is pure-null. Hence (7) is equivalent to
all positive pure moments vanishing, and supplies the sufficient prefix.
For c+P_a, its first moment is c, so the same statement applies to the
whole space (3).

Seven is least. The explicit P_*=p_0+p_7 has weights 8 and -6. Weight
balance is impossible at positive orders below seven. Direct integration
and the universal certificate agree on

    h(P_*^7)=-8811708416/154040315 !=0.                    (8)

Equivalently its parameter octavic is X^8+8XY^7, which has eight distinct
roots. The value in (8) is also computed independently as
35 h(p_0^3 p_7^4) using the determinant-one Haar monomial integral.
This witness has six zero moments followed by a nonzero seventh; it is
not an all-order pure-null function.

If f_a has a root of multiplicity at least five, move that root to a
coordinate point. All surviving substitution weights of P_a then have
one strict sign, and their absolute values are at least two. A marker
of matrix-entry degree e has weights between -e and e. Thus 2m>e
forces its product with P_a^m to have no weight zero, proving (1).
Left translation preserves matrix-entry degree, so this is uniform on
the entire pure-null space, not just a fixed root chart.

The zero-offset coefficient 1/2 is optimal: take P=p_3 and Q=p_8.
Then P has weight2 and is pure-null, Q has exact degree8, and

    h(QP^4)=-67436544/230945 !=0.                         (9)

Both the full fifth-moment coefficient and direct Haar multiplication
give (9). Any smaller coefficient would wrongly force this value to zero.
This proves optimality of the uniform zero-offset bound; it does not
assert an optimal affine asymptotic slope with an arbitrary added offset.

## 5. The full ordinary moment ideal and its least prefix

Let R=C[a_0,...,a_8]. Since all H_m are positive-degree invariants and
the invariant ring is generated by H_2,...,H_10, its full moment ideal is

    I_infinity=(H_m:m>=1)=(H_2,...,H_10)
               =(J_2,...,J_10) in R.                    (10)

This is an ordinary-ideal identity, not merely an equality of radicals.
The least generating prefix is ten. To see strictness at the last step,
suppose H_10 were in (H_2,...,H_9). Homogeneous coefficients could be
chosen in the required complementary degrees. Averaging those coefficients
over the parameter SL_2 makes them invariant without changing the identity.
Each has degree at most eight, where the complete invariant bases above
express it as a polynomial in earlier H_i. This would make H_10 a linear
combination of products of H_2,...,H_9 of degree ten. The checked full
column rank and nonzero gamma_10 exclude exactly that relation.

Adding the constant coefficient c gives the full ideal

    (c,J_2,...,J_10) in C[c,a_0,...,a_8],                  (11)

again with least prefix ten. The binomial expansion of h((c+P_a)^m)
and the first-moment generator c prove this without specialization loss.
The sharp set prefix remains seven. In particular I_7<I_infinity even
though their radicals agree. No assertion that these ordinary ideals are
reduced or a complete intersection is made.

## 6. Verification and remaining model questions

The checker rebuilds the universal moment polynomials and primitive
transvectants, verifies every frozen coefficient identity and full column
rank through degree ten, and checks the invariant dimensions by weight
counts. It independently integrates the ten powers of p_4 by matrix
monomials, as well as both sharp controls (8)--(9). It checks the exact
binary-octahedral group and fixed dimensions through degree12. The
projector normalization mutation is rejected. The two classical all-degree
invariant-ring inputs remain imported, not formally certified here.

The next full octahedral degree space is

    A_(degree<=12)=C plus V_8 plus V_12,

with one copy of each nonconstant type. For the icosahedral model the
first full space remains C plus V_12. Neither space nor either all-degree
MZ model is decided in this note. Both models now require counterexample
degree at least12. The nullcone classification remains closed by Haar
shears; the current result shows that its octahedral counterexamples need
additional representation types beyond this lowest module.
