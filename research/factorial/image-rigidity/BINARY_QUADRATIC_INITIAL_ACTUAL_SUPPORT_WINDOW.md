# Every initial actual-support window for a binary quadratic

For a characteristic-zero field K let Lambda:K[X,Y]->K be the K-linear map
Lambda(X^iY^j)=i!j!. For a polynomial f write N(f) for its number of
nonzero monomials. This note proves an INITIAL-window statement, with a
finite exact polynomial certificate over all coefficient values.

**Theorem 1 (initial binary-quadratic window).** If f is a nonzero
polynomial in K[X,Y] of total degree at most2, and s=N(f), then

             (Lambda(f),Lambda(f^2),...,Lambda(f^s)) !=0.       (1)

This includes every actual support, including constants, affine phases,
all quadratic ranks, and all complex coefficient values. The proof does
not establish the analogous window at every positive starting index.
In particular it does not settle full binary-quadratic SFC or SFC(2).

## 1. The exact support domain and coefficient conventions

Order the six monomials and their exponent pairs by

 index:       0      1      2      3      4      5
 monomial:    1      X      Y      X^2    XY     Y^2
 exponent:   (0,0)  (1,0)  (0,1)  (2,0)  (1,1)  (0,2).

A mask b in{1,...,63} selects precisely the indices whose bits are1.
Let S be the selected set and s=|S|. Write m_i=X^(nu_i1)Y^(nu_i2) and

 L_i=Lambda(m_i)=nu_i1!nu_i2! in{1,2},
 f=sum_(i in S)c_i m_i.                                      (2)

Here c_i are the ACTUAL monomial coefficients. Equivalently, using the
natural divided-factorial coordinates v_i=L_i c_i, one can write
f=sum v_i m_i/L_i, with Lambda(f)=sum v_i. The two coefficient conventions
are related by invertible diagonal scalars, so the support is unchanged.

We prove a stronger assertion separately for every S: the simultaneous
equations Lambda(f^j)=0 for1<=j<=s have only the zero coefficient vector
in the entire span of S, allowing some of its c_i to vanish. Applying
this assertion to the actual support proves(1). A six-term certificate
alone would not justify every smaller support-size window; the data
therefore includes all63 nonempty masks independently.

## 2. Injective elimination of the first moment

If s=1, then Lambda(f)=L_i c_i and L_i!=0, proving the assertion directly.
These are the six singleton cases.

Suppose s>=2. Choose the least index p in S, put R=S without p and r=s-1,
and introduce parameters u_i, i in R. Set

 c_i=L_p u_i for i in R,       c_p=-sum_(i in R)L_i u_i.       (3)

This is an injective K-linear map onto the entire hyperplane Lambda(f)=0.
Indeed the first moment vanishes identically, and its inverse on that
hyperplane is u_i=c_i/L_p. The equation for c_p then follows from
L_p c_p+sum L_i c_i=0. Thus no nonzero coefficient vector is lost, and
no assumption that an individual unknown c_i is a unit or nonzero is made.
The fixed scalar L_p is1 or2, hence invertible in characteristic zero and
modulo the certificate prime101.

The resulting phase and homogeneous native moments are

 f_S(u;X,Y)=sum_(i in R)u_i(L_p m_i-L_i m_p),
 Phi_(S,j)(u)=Lambda(f_S(u)^j), 2<=j<=s.                    (4)

Every Phi_(S,j) belongs to Z[u] and is homogeneous of degree j.

For clarity the checker reconstructs its coefficients directly, without
symbolic integration or a computer algebra dependency. In the original
support multinomial expansion let alpha=(alpha_i)_(i in S), |alpha|=j,
t=alpha_p, and beta=(beta_i)_(i in R), |beta|=t. The contribution is

 [j!/product_i alpha_i!] *
 (sum_i alpha_i nu_i1)! (sum_i alpha_i nu_i2)! *
 L_p^(j-t) (-1)^t *
 [t!/product_(i in R)beta_i!] *
 product_(i in R)L_i^beta_i *
 product_(i in R)u_i^(alpha_i+beta_i).                      (5)

Summing over all these alpha and beta gives EXACTLY Phi_(S,j).
This follows by the ordinary multinomial theorem and substituting(3).
All quotients displayed as multinomial coefficients are integers.
The degree is j because sum_(i in R)(alpha_i+beta_i)=j.

## 3. The finite modular polynomial certificates

The accompanying file
[binary_quadratic_initial_support_certificates.json](binary_quadratic_initial_support_certificates.json)
contains one record for each mask1 through63. For every s>=2 it specifies,
over F101, a square (s-1)-by-(s-1) matrix of multiplier polynomials U,
and uses the exact degree

                       D_s=s(s-1)/2+1.                    (6)

With the R variables listed in increasing original-index order, every
column h of that matrix satisfies the identity

 sum_(j=2)^s Phi_(S,j)(u) U_(j-2,h)(u)=u_h^(D_s).           (7)

The multiplier in row j-2 is homogeneous of degree D_s-j.
Coefficients are canonical integers1 through100, reduced modulo101.
No parameter specialization is present: these are identities in the
whole polynomial ring F101[u_1,...,u_(s-1)].

The file records57 nonsingleton supports,129 identities and60,779 nonzero
multiplier terms. The identity count also follows from
sum_(s=2)^6 binom(6,s)(s-1)=129. All certificates use the same prime101.
The singleton records contain their nonzero scalar first moments.

The standard-library checker
[verify_binary_quadratic_initial_support_certificates.py](verify_binary_quadratic_initial_support_certificates.py)
does not import Singular or any computer algebra package. It independently
reconstructs(5) with integer arithmetic, reduces modulo101, multiplies
every polynomial on the left of(7), and compares EVERY coefficient with
the stated pure power. It also checks all masks, coefficient domains,
degrees, dimensions and variable conventions. Its separate literal native
form hash catches changes to the factorial transcription.

Singular was used only to propose these finite multipliers. A reported
Groebner-basis dimension, numerical rank, finite list of coefficient points,
or unverified computer algebra assertion is not a premise.
The exact finite polynomial identities(7), verified coefficient by
coefficient, are the certificate.

## 4. Why these identities prove the characteristic-zero assertion

Here is an elementary transfer using finite linear algebra; it needs
neither a projective properness theorem nor a p-adic approximation.

Fix a nonsingleton support, put r=s-1 and D=D_s, and set
q=r(D-1)+1. Every monomial of degree q in r variables is divisible by
some u_h^D: otherwise its degree would be at most r(D-1).
Multiply the corresponding identity(7) by the complementary monomial.
Consequently, over F101, the degree-q part of the ideal generated by
Phi_(S,2),...,Phi_(S,s) is the entire space of degree-q forms.

Consider the finite INTEGER matrix whose columns are the coefficient
vectors of all products Phi_(S,j) times degree-(q-j) monomials, expressed
in the monomial basis of degree q. The preceding identities prove that
this matrix has full row rank modulo101. Hence one of its square maximal
minors is nonzero modulo101 and is therefore a nonzero integer.
That same minor is nonzero over Q and over every characteristic-zero
field K. The matrix has full row rank over K as well.

It follows that every degree-q monomial belongs to the ideal generated
by the native Phi_(S,j) over K. In particular it contains every u_h^q.
A simultaneous K-valued zero of all the Phi_(S,j) must therefore have
u_h=0 for every h. By the injective parameterization(3), its original
coefficient vector is zero. Together with the singleton cases, this
proves Theorem1.

The matrix in this argument need not be materialized: identities(7)
already supply its full-rank certificate modulo101. For the full
six-term case s6, D6=16 and q=76. We do NOT mistake the integer lifts of
the F101 multipliers for rational ideal-membership identities; only
the full-rank implication is used for the transfer.

For comparison, the standard projective explanation gives the same
conclusion. A nonzero complex projective zero of rational homogeneous
equations would yield a Qbar projective zero, hence a number-field point.
Embed that field into Q101bar, scale its coordinates so their minimum
valuation is0, and reduce in the algebraic closure of the residue field.
Homogeneity and the integral native coefficients give a NONZERO F101bar
zero of the Phi polynomials. Identities(7) force all its coordinates to
vanish, a contradiction. This explanation uses the usual algebraic
closed-field point and number-field completion facts; the preceding
integer-matrix proof is an independent elementary transfer.

## 5. Replay and limits

From the directory containing this note, run a fresh verification with

    python3 verify_binary_quadratic_initial_support_certificates.py --output fresh_initial_quadratic_check.json

The checker finds the certificate beside itself by default. It rejects
optimized Python execution and refuses to overwrite an existing output.
The retained
[fresh exact receipt](binary_quadratic_initial_support_exact_verification.json)
covers all129 identities and seven deliberate corruptions, including
a changed multiplier coefficient, a missing face, a wrong native moment,
a nonprime modulus and an attempted all-start scope extension.

The theorem quantifies over every characteristic-zero coefficient field
and EVERY coefficient choice, via the finite polynomial identities.
It covers the initial actual-support window only. It is not an
inference from finitely many starting indices, and no later-start,
full SFC(2), SFC(3), arbitrary degree, or worldwide-priority conclusion
is asserted.

## 6. Bounded primary-source comparison

Edo--van den Essen, The Strong Factorial Conjecture, arXiv1304.3956v2,
Definition2.7 and Conjecture2.8, distinguish the initial window from the intersection
over every positive start. Their Theorem2.25 and Corollary2.28 treat the
specified monomial-times-linear family E^[m], rather than every binary
quadratic coefficient:
https://arxiv.org/pdf/1304.3956 .

Rocks, On the Incompatability of Diophantine Equations Arising from the
Strong Factorial Conjecture, author preprint dated6September2015,
Theorems3.3,3.10,3.12 and Corollary3.13, treat powers of linear forms,
separate-variable prime powers and independent-variable adjunction:
https://people.math.sc.edu/rocks/BradyRocksSFC2015.pdf .

Liu--Sun, The factorial conjecture and images of locally nilpotent
derivations, Bull.Aust.Math.Soc.101(2020),71--79,
DOI10.1017/S0004972719000546, proves FC for homogeneous binary phases,
a statement about vanishing at every positive order. That theorem is
not an initial support-size window for arbitrary inhomogeneous quadratics:
https://doi.org/10.1017/S0004972719000546 .

These actual primary statements and the repository's retained initial
monomial-times-linear resultant theorem were checked to separate their
scope from Theorem1. This bounded comparison is not an exhaustive
literature or priority claim. None of these papers is a mathematical
premise of the native modular certificate proof.
