# Complete high-contact quintic cusp locus and the excess-one singleton gate

14 November 2026. PF2R7B5CUSP2 completes the finite cusp locus left by
PF2R7B5CUSP1. There are161 distinct normalized(C,B) polynomial-image
pairs over its160 first-image coefficient points. At159 points the
unique nonimmersive B has one ordinary cusp and three ordinary nodes.
At the special d0 point there are two B choices, each with one ordinary
cusp and one ordinary triple point. Every cusp is clean and lies off C.
All these pairs have infinity contact38.

PF2R7B5B1 consequently excludes singleton fibres in the immersed
excess-one row B whenever the additional image has same-coordinate
polynomial normalization degrees(3,5) and infinity contact at least38.
Together with the earlier nonimmersed-row result PF2R7B5N1, this leaves
only the immersed rows A or C as singleton explanations within that
explicit degree/contact chart.

Neither the degree chart nor high contact is forced for actual B.
No global Keller cover or source A2 is constructed, and absence of a
singleton is not itself a contradiction to a Keller map. The proof is
authored algebraic geometry with exact certificates, without formalization
or independent external review. Whole r7 and JC2 remain open.

## 1. Retain the complete finite coefficient algebra

PF2R7B5REC1 supplies the rank-seven coefficient algebra E over QQ[d]
or QQ[z], with its monic weighted Groebner basis and fixed standard
monomials. Its uniquely normalized image parametrization is

    P(s)=s^3+A*s+M,
    Q(s)=s^5+b*s^4+c*s^3+e*s^2+f*s+n, c=5*A/3.

All remaining coefficients are polynomial expressions in A,M and the
first-image parameter. Let D be the seven-by-seven multiplication
matrix of Delta=Res(P',Q'). Its rational entries and determinant are
reconstructed by the canonical norm code, not accepted from a printed
solver basis. PF2R7B5CUSP1 gives

    det D=k_L*d^2*F12(d^3),  det D=k_K*z^15*G41(z),

with nonzero rational scalars and primitive squarefree factors, nonzero
at0. Write H=F12(d^3) of degree36, or H=G41(z) of degree41. These
cover every nonspecial nonimmersion parameter; d0 is treated separately.
The three shared L/K coefficient points do not meet either set.

## 2. Recover the unique nonimmersive candidate without choosing roots

Work over the finite reduced base T=QQ[parameter]/H. Let ell be row0
of the adjugate of the rational matrix D. Its jth entry is the signed
minor obtained by deleting row j and column0. The universal identity

    ell*D=(det D,0,...,0)

places ell in the left kernel over T. The constant standard monomial
has index0. The certificate proves ell0 is a unit on the entire base T.
Consequently D has rank exactly6 at every complex root of H: det D0,
but its indicated six-by-six minor is nonzero.

For each such root, E/(Delta) is a one-dimensional unital algebra.
It therefore has exactly one closed coefficient point, and the unique
nonimmersive candidate is recovered by

    A*=ell_(A)/ell0, M*=ell_(M)/ell0.                (1)

Here ell_(A),ell_(M) are indexed by those standard monomials. These
are rational functions on T defined by minors of the displayed rational
matrix, not numerical algebraic-root approximations. No irreducibility
of H or of the seven-sheet correspondence is required.

The unit and subsequent incidence tests are certified modulo103 on L
and modulo101 on K. All rational input denominators and the leading
coefficient of H survive these reductions. The code computes the actual
minor reductions, checks the left-kernel identity and proves
gcd(ell0,H)=1 modulo the chosen prime. This gives a nonzero rational
resultant, so ell0 is a unit over T in characteristic zero.

Every further denominator inverted below is likewise certified a unit
modulo that prime. Its characteristic-zero inverse exists, with reduction
equal to the inverse used in the code. Thus the later matrices are
reductions of canonically defined rational finite-algebra matrices.
An arbitrary modular unit ideal would not justify this lifting.

The initial prime101 on L gave a degree3 common factor in the chosen
minor test. That reduction supplied no characteristic-zero conclusion.
The successful103 certificate covers every complex root of H. The
rejected reduction is retained in discovery.json; no prime-dependent
exception is silently discarded.

## 3. A single clean ordinary cusp at every nonspecial parameter

Evaluate the five Q coefficients at(1). Modulo P'=3s^2+A*, the
derivative Q' has linear remainder

    L*s+K,
    L=2e-4A*b/3, K=5A^2/9-cA+f.

The certificate proves L is a unit on T. Hence the unique critical
parameter is u=-K/L. The critical resultant zero implies
A*=-3u^2, and the critical equation gives

    f=10u^4-4b*u^3-2e*u.

There is no other common derivative root. The ordinary-cusp jet
determinant is

    P''(u)Q'''(u)-P'''(u)Q''(u)
       =12*(20u^3+6b*u^2-e).                       (2)

This factor is a unit on T. The image branch is therefore an ordinary
cusp. For its complete target fibre,

    P(s)-P(u)=(s-u)^2*(s+2u),
    Q(-2u)-Q(u)=9u^2*(e+3b*u^2-2u^3).             (3)

The last parenthesis is also a unit. If u is nonzero, the potential
additional parameter -2u is absent. If u0, the fibre of P has only
the root0, and(2) gives e nonzero. In both cases the cusp is clean:
there is exactly one normalization point over its value. In particular
it cannot supply the b(y*)=2 singleton option of the omitted-prime ledger.

## 4. The remaining three singularities are ordinary nodes

The collision quartic R(v) from PF2R7B5GEN1 has the diagonal root2u.
The code divides by v-2u exactly in T[v] and obtains a monic cubic R3.
Multiplication by R3' in T[v]/R3 has full modular rank108 on L and123
on K. Thus all three roots are simple, over every complex root of H.
The ordinary-cusp jet factor in(2) is also R'(2u), so these roots differ
from the removed diagonal root. The unique critical parameter already
excludes any other diagonal pair.

The cubic-fibre test

    A^2+A*b^2-b*e-A*c+f

is a unit on T. Thus the three remaining unordered pairs give distinct
target values rather than a triple value. They also avoid the clean
cusp target by Section3. Simple off-diagonal roots give transverse
normalization pairs. Hence these are precisely three ordinary nodes.
The collision quartic exhausts all finite identified pairs and the
derivative calculation exhausts all nonimmersions, so no further
affine singularities are missing.

The U^2 coefficient k6 of G_C-F_B is another certified unit on T.
The first-image pole basis therefore gives infinity contact exactly38
and total finite intersection multiplicity12 for each of these pairs.

## 5. The cusp value does not lie on C

Its coordinates in T are

    U*=M*-2u^3, V*=Q(u).

Use the monic finite-free algebra

    T[t]/(p_C(t)-U*).

Its dimension over QQ is216 on L and246 on K. Multiplication by
q_C(t)-V* has full rank modulo the respective good prime. By the unit
lifting in Section2 its rational determinant is nonzero. There is no
complex solution of p_C(t)=U*,q_C(t)=V* at any root of H. The polynomial
normalization of C is finite and onto its image, so the cusp target is
off C, not merely absent from a selected parametrization chart.

This tests all36 L parameters and all41 scaled K parameters at once.
It does not assert that all other C/B incidences satisfy the global
covering constraints.

## 6. Complete special fibre d0

The first two coefficient equations specialize to

    A^3-9M^2-(12/5)M-1/15=0,
    A^2*(M+2/15)=0.

If A is nonzero, then M=-2/15 and A^3=-7/75. Its critical resultant
is25A^4, nonzero, so this branch contributes no nonimmersive image.
If A0, the exact equation is

    135M^2+36M+1=0.

It has two distinct roots. Put e=5M/3+2/9, a unit in this quadratic
algebra. The two normalizations are

    P=s^3+M, Q=s^5+e*s^2.

Each has a unique clean ordinary cusp at s0, target(M,0). Its other
three identified pairs, cut out by v^3=e in the collision quartic
v*(v^3-e), form a single ordinary triple point at(M-e,0). The three
source parameters have s^3=-e. Their tangent slopes -e/s are distinct,
so the triple is ordinary. There are no other finite singularities.

Direct resultant computation on the known C at d0 gives the monic
V0 restriction

    Res_t(p_C(t)-M,q_C(t))
      =M*(6561M^3+2916M^2+378M+14)^3/282429536481.

It is coprime to the quadratic above, so both cusp targets are off C.
The k6 coefficient is also a unit; its quadratic remainder is
-14M/3645-7/10935. Contact is again exactly38.

## 7. Counts and the immersed excess-one row

There is one nonimmersive B over each of the36 nonspecial L parameters
and over each of the123 nonzero K x parameters, plus two over d0.
There are therefore161 distinct normalized(C,B) coefficient pairs over
160 first-image coefficient points. The uniqueness of the centered
monic normalization identifies distinct coefficient choices with
distinct B curves. These are counts in the retained coefficient
normalization, not inequivalent affine moduli. This proves PF2R7B5CUSP2.

In residual immersed row B, PF2OMIT6/PF2OMIT7 require exactly one
nonimmersive parameter, an ordinary cusp y*. Every high-contact
degree(3,5) possibility is in the complete list above. At y* the
normalization branch counts are c0 and b1: the cusp is clean and off C.
The retained source cardinality is consequently

    6-2c-2b-1=3.

At every other value the count is6-2c-2b and is even. Thus no singleton
fibre occurs anywhere in this entire row under the stated additional
degree/contact hypothesis. This is PF2R7B5B1.

Combined with PF2R7B5N1, a singleton in this chart must come from an
immersed first image and one of rows A or C. Those rows remain open;
generic local incidences already pass. Other degree/contact charts,
the optional H image, actual joint complement relations and the
normal-cover/source construction remain. A Keller map without a
singleton is not excluded merely by that absence. Neither unrestricted
singleton-zero properness nor JC2 is proved here.
