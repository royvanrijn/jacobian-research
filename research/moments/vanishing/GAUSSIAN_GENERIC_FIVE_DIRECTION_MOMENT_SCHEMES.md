# Generic Gaussian moment schemes in five centered directions

The theorem below gives an explicit generic coefficient open, exact
ordinary and geometric moment prefixes, and a finite flat family.
All original polynomial markers and coefficient-algebra qualifications
are retained. This is an internal written proof with exact finite
determinant verification; formal verification and external review are
not asserted.

Let K be ANY characteristic-zero field, let E be the normalized Gaussian
Wick functional in n>=3 variables with nonsingular covariance, and fix
d>=3. Let V_d=ker(E) inside the spatial polynomials of degree<=d. For
Q=(Q1,...,Q5) in V_d^5 put

    P_Q(u,X)=sum_(i=1)^5 u_i Q_i(X),
    mu_j(Q,u)=E[P_Q^j], I_N=(mu_1,...,mu_N),
    m=(u_1,...,u_5).

The first moment is identically zero. There is an explicitly specified
nonempty Zariski open U in V_d^5 on which, globally in K[u],

    I_infinity=I_10 != I_9,
    m^11 subset I_8, m^10 not subset I_10.          (G1)

The moment quotient has length546 and graded Hilbert vector

    h=(1,5,14,29,49,71,90,100,95,70,22),            (G2)

followed by zero. Thus the ORDINARY least generating prefix is ten and
the sharp coefficient-ideal nilpotence index is eleven. The least
GEOMETRIC set-theoretic prefix, equivalently the least N for
radical(I_N)=m, is six. This geometric statement concerns algebraically
closed points; least cutoffs for only K-rational points can be smaller.

Every pure phase on these five-direction slices has u=0. This does not
settle safety-to-purity: safety need not make the first unit moments zero.
All conclusions concern slices through the zero phase, not automatically
affine transverse slices at an arbitrary nonzero pure phase.

## 1. Explicit universal determinant opens

Choose any K-basis of V_d and use its coefficients q for the five Q_i.
The coefficient ring is R_0=K[q]. The moments are polynomial in q,u and
homogeneous of degree j both in u and in q. For each r form the literal
degree-r Macaulay matrix: rows are u^beta mu_j, with j<=N and
abs(beta)=r-j; columns are the degree-r monomials in u. Use descending
lexicographic order within every fixed degree, with moment degrees
ordered increasingly. No coefficient content is divided out universally.

For r=2,...,11 use the selected row/column index lists in the retained
[axial minor certificate](gaussian_axial_cubic_macaulay_minors.json). Call the
corresponding UNSCALED universal determinants Delta_r(q). Also use the
degree16, prefix6 selected index lists in
[complete-intersection minor](gaussian_generic_five_direction_ci_minor.json); call that unscaled
universal determinant epsilon(q). Define

    Omega=epsilon * product_(r=2)^11 Delta_r,
    U={Omega!=0}.                                  (G3)

This is an explicit straight-line/determinant description. Expanding
the very large polynomial Omega is not required to define or verify it.

The degree-eleven selected minor is square of size1365 and uses ONLY
mu_2,...,mu_8. Its selected row counts are

    j:       2   3   4   5   6   7   8
    count: 715 285 149  96  71  46   3.

Consequently Delta_11 has homogeneous q-degree4133. The adjugate identity
for this matrix gives the UNIVERSAL containment over R_0[u]

    Delta_11 * m^11 subset (mu_2,...,mu_8).          (G4)

Indeed, if z is the column vector of degree-eleven monomials and f is
the vector of its selected moment multiples, then f=M z and
adj(M) f=det(M) z. Every entry of f belongs to I_8. Multiplication by
arbitrary polynomials gives (G4) for the entire ideal m^11.

Likewise the degree-sixteen prefix-six full-column matrix gives

    epsilon * m^16 subset (mu_2,...,mu_6).          (G5)

These identities hold before specializing q. On U their determinants
are units, giving m^11 subset I_8 and m^16 subset I_6 over every
coefficient algebra in which Omega is invertible.

## 2. Why the open is nonempty for every covariance and d>=3

Over an algebraic closure of K, any nonsingular Gaussian in dimension
at least three admits normalized linear coordinates Z,W,Y with
E[ZW]=E[Y^2]=1 and all their other pairings zero. Use the centered tuple

    Q1=Y^3-3Y, Q2=(ZW-1)Y, Q3=Y^2-1,
    Q4=ZW-1, Q5=Y.                                 (G6)

Its members have degree<=3. In the axial canonical data, primitive
moments g_j are obtained from mu_j by integer contents c_j=j-1,
2<=j<=10. Each primitive selected determinant is nonzero modulo32003.
The degree-eleven primitive determinant is17570; restoring its row
contents gives the unscaled determinant5815 modulo32003.

The exact original moment formula at this witness is as follows. For
i+j+h+k+l=N, the coefficient of u1^i u2^j u3^h u4^k u5^l in mu_N is

    N!/(i!j!h!k!l!) D_(j+k)
       *E_Y[He_3(Y)^i He_2(Y)^h Y^(j+l)],

where D_0=1,D_1=0 and D_N=(N-1)(D_(N-1)+D_(N-2)) are the centered
Gamma moments E[(S-1)^N]. This follows by the multinomial theorem and
independence, with E[S^N]=N! and ordinary Gaussian moments. It binds
all determinant entries to the original unweighted Wick functional.
The checker regenerates these finite sums and compares21 small direct
original S,Y expansion fixtures independently.

For epsilon, the retained literal degree-sixteen minor has4845 columns,
4845 selected rows among9626, and determinant6771 modulo32003 for the
primitive mu_2,...,mu_6. Restoring row contents also stays nonzero.
Thus every factor in (G3) is nonzero at (G6), over Q and therefore over
every characteristic-zero field after this Gaussian coordinate change.

The polynomial Omega was defined over K, so its nonzero value after
extension to Kbar proves it was nonzero already over K. Since K is
infinite, its nonvanishing open has K-points as well. This argument
does not require a split Witt pair or square roots in K itself.

## 3. Complete intersection supplies a rational LOWER Hilbert bound

Fix any field-valued tuple in U. By (G5), I_6 is primary to m and its
quotient has dimension zero. In the regular local ring K[u]_m of
dimension five, the five positive-degree forms mu_2,...,mu_6 are
therefore a regular sequence. The exact primary inputs are Stacks
[Lemma10.106.3, tag00NQ](https://stacks.math.columbia.edu/tag/00NQ)
(regular local rings are Cohen--Macaulay), and
[Lemma10.104.2, tag02JN](https://stacks.math.columbia.edu/tag/02JN)
(the dimension criterion for a regular sequence).

This regularity is global in the graded polynomial ring. A nonzero
graded kernel of multiplication by a homogeneous form has a nonzero
homogeneous element; localization at m cannot kill it, since an
annihilator outside m has a nonzero scalar lowest-degree part.

Successive injective multiplication exact sequences give

    Hilb(C,z)=product_(j=2)^6 (1-z^j)/(1-z)^5,
    C=K[u]/I_6,
    c=(1,5,14,29,49,71,90,101,101,90,71,49,29,14,5,1).

In each r<=10 the additional ideal is the image of

    direct_sum_(j=7)^min(10,r) C_(r-j) -> C_r.

Its rank is at most its domain dimension. Hence

    dim(K[u]/I_10)_r >= c_r-sum_(j=7)^min(10,r)c_(r-j)=h_r. (G7)

For r<=6 this is c_r. For r=7,8,9,10 it gives100,95,70,22.
This elementary complete-intersection bound is unconditional; no
generic-Hilbert-series conjecture or modular rank deficiency is used.

## 4. Upper bounds, the entire tail, and both sharp prefixes

The Delta_r in (G3) are the retained-size nonzero Macaulay minors, so
they give the reverse dimension bounds <=h_r for r<=10. The full
Delta_11 gives quotient degree11 zero. Combined with (G7) this proves
(G2), length546 and sharp nilpotence index11. Homogeneity puts EVERY
moment of order>=11 in I_8, proving I_infinity=I_10 globally.

For I_9 the complete-intersection argument in degree ten gives

    dim(K[u]/I_9)_10 >= c_10-c_3-c_2-c_1
                         =71-29-14-5=23.

But (G2) gives dimension22 for I_10. Thus I_9!=I_10 and the least
ordinary ideal prefix is ten. A rational dual separator is unnecessary.

The first six moments already have radical m by (G5). The first five
moments generate the first four forms of the same regular sequence.
Their quotient has dimension one by the same Stacks dimension
criterion, so its radical cannot equal m.
This proves the least geometric set-theoretic prefix is six. It makes no
claim of a least prefix for points rational over a nonclosed field.

## 5. Coefficient algebras, original markers and a universal free basis

Over any commutative K-algebra T with Omega(q) a unit, (G4) is an actual
ideal identity. If the first eight moments vanish then the ideal J of
the five actual direction coefficients has J^11=0. Therefore P^11=0 as
a spatial polynomial and E_T[A(X)P^m]=0 for every ORIGINAL polynomial
marker A and every m>=11. The first six moments alone similarly give
J^16=0 and the original marked cutoff m>=16. Coefficients need not
vanish in a nonreduced T; no all-ring primitive-vacuum claim is made.

The universal quotient R_0[Omega^-1][u]/I_10 is finite FREE of rank546
over R_0[Omega^-1]. For every geometric fibre its graded dimensions are
the fixed h_r. The coefficient ring is reduced, so all larger Macaulay
minors vanish identically. The selected Delta_r is already a unit on
this entire open. Its fixed pivot block therefore reduces the degree-r
matrix globally to an identity block and zeros. The nonselected
monomial columns give a free basis of the graded quotient of size h_r.
Only degrees0..10 occur, so the union of these monomials is a universal
free basis of size546. In particular the moment family is finite flat;
no further moment computation is required.

## 6. Evidence and limits

The new CI minor is a literal selected determinant. Its verification
does not call Singular, a Groebner generator, RREF or rational
reconstruction. A modular Groebner computation was used only as a
generation probe; the stored nonzero integer minor is the proof witness.
The other ten minors are the existing exact axial inputs. No rational
Groebner superset or degree-ten dual certificate is needed for this
generic proof: the complete intersection provides the lower dimension
bound, and its degree-ten estimate provides the strict-prefix proof.

The [primary manifest](sources/generic-five-direction/manifest.json)
pins the inspected Stacks
[regular-local lemma](sources/generic-five-direction/00NQ.html) and
[dimension criterion](sources/generic-five-direction/02JN.html).
Their exact hypotheses are applied in Section3. The Hilbert series
and rank/adjugate arguments are derived above; no analytic comparison
or universal Gaussian safety-to-purity theorem is invoked.

Run

    python verify_gaussian_generic_five_direction_moment_schemes.py \
      --output gaussian_generic_five_direction_moment_schemes_verification.json

from this directory, with Python and python-flint. The
[standalone verifier](verify_gaussian_generic_five_direction_moment_schemes.py)
rebuilds all eleven literal minors, compares21 independent original
moment controls, checks the finite Hilbert arithmetic and rejects
duplicate-row and wrong-prefix metadata corruptions. Its
[receipt](gaussian_generic_five_direction_moment_schemes_verification.json)
pins this proof, code, literal inputs and primary mirrors. Default
replay calls no Groebner generator, RREF or rational reconstruction
and has no runtime scratch-file dependency. The finite checks verify
the witnesses; the unrestricted generic and coefficient-algebra
conclusions use the written proof. Independent replay, formal
verification and external review are not asserted.
