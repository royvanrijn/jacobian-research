# Every eventually mixed-safe Gaussian polynomial has a linear degree cutoff

15 November 2026. Internal all-degree, all-dimension consequence of the
integration b-function already used in the
[finite multiplier theorem](FINITE_MULTIPLIER_MIXED_DECISION.md).
No novelty, independent review or formal verification is claimed.
The universal construction is written, not a general implemented solver.
The [checker](verify_linear_mixed_cutoff.py) verifies finite Weyl identities,
normal-order signs, exact recurrences and rational cutoff extraction on
explicit Gaussian certificates. No full coefficient campaign is run.

## 1. The qualitative property already forces a linear bound

Let P in C[X_1,...,X_n] be arbitrary. Suppose that, for every polynomial Q,

    F_Q(t)=E(Q exp(tP))

is a polynomial in t. This is exactly eventual mixed vanishing, with its
cutoff initially allowed to depend arbitrarily on Q.

**Theorem.** There are nonnegative rational numbers A(P),C(P) such that

    deg_t F_Q <= A(P)+C(P) deg Q                              (1)

for every nonzero Q, with the zero series ignored in degree comparisons.

If, in addition, all positive pure moments of P vanish, there is a
nonnegative rational C_pure(P) such that

    E(QP^m)=0 whenever m>C_pure(P) deg Q.                      (2)

This holds for arbitrary complex coefficients. Over a specified effective
finitely generated characteristic-zero field, the constants can be computed
by a terminating exact construction. The constants obtained are sufficient;
their optimality is not claimed.

For an individual pure-null P, therefore, qualitative GMC safety is
equivalent to having a linear degree-Q cutoff with a coefficient depending
on P. This does not prove that every pure-null P is safe: the retained
quartic failure remains a failure. Nor is the coefficient independent of
P or of its degree. The previous finite multiplier theorem did not claim
linear growth; the additional argument here supplies it.

## 2. Retain a b-function identity, including its lower-weight terms

Use the polynomial Weyl algebra in time t and Gaussian variables X.
Put

    Psi=tP(X)-sum_i X_i^2/2,
    theta=-sum_i partial_(X_i) X_i,
    wt(X_i)=1, wt(partial_(X_i))=-1,
    wt(t)=wt(partial_t)=0.

The exponential connection is cyclic and holonomic. The classical
integration b-function statement gives a nonzero b(s) in C[s] and an
operator R of weight at most -1 such that

    (b(theta)+R) exp(Psi)=0.                                  (3)

It is important to retain R, not just the integer roots of b. In normal
order write

    R=sum_l c_l t^(k_l) partial_t^(j_l)
                  X^(u_l) partial_X^(v_l),                    (4)
    r_l=|v_l|-|u_l|>=1.

The coefficients c_l are scalars, and all exponents are nonnegative.
Terms are collected before extracting degree bounds.

The same classical statement appears in
[Oaku 2018, Algorithms for D-modules, integration, and generalized functions
with applications to statistics, theorem 5.4](https://www.lab.twcu.ac.jp/oaku/MSJSI2015oakuLecture.pdf),
pp.311--314, with the algorithmic lifting framework in algorithm 5.6.
This is the enlarged 2018 publication of the MSJ SI 2015 lectures, not a
new 2026 paper. Theorem 5.4 and its proof were read; its graded identity
says b(theta) lies in the annihilator modulo weight -1.
Lifting gives (3). The corresponding identity was already present in the
[earlier Oaku integration source](https://www.lab.twcu.ac.jp/oaku/Oaku2013.pdf).
No analytic integration assertion is needed.

An exact certificate can be verified by applying b(theta)+R to exp(Psi)
and dividing out exp(Psi). Every derivative has polynomial amplitude, so
this is a finite polynomial identity. A general b-function can be obtained
by the classical Weyl Groebner algorithm and a lift. Alternatively enumerate
finite boxes for b and the normal-order coefficients of R, enforce one
coefficient of b to be one, and solve the finite linear identity over the
effective coefficient field. Existence of (3) guarantees eventual success.
No such general search is executed here.

## 3. A scalar recurrence that lowers multiplier degree

Let alpha be a multi-index and d=|alpha|. Multiply (3) on the left by
X^alpha and pass to the algebraic integration quotient, or apply formal
Gaussian integration by parts coefficientwise.

Because

    X^alpha theta=(theta+d)X^alpha

and every leftmost total X derivative integrates to zero,

    Lambda(X^alpha b(theta))=b(d) F_(X^alpha).

For each normal-order term of R, moving the X derivatives to the left gives

    Lambda(X^(alpha+u) partial_X^v)
      =(-1)^|v| (alpha+u)_v F_(X^(alpha+u-v)),                 (5)

where (eta)_v is the product of falling factorials. A term with a negative
coordinate in alpha+u-v is zero and is omitted. The time operators commute
with the X integration and stay in the polynomial Weyl algebra. Hence

    b(d) F_(X^alpha)
      =-sum_l c_l (-1)^|v_l| (alpha+u_l)_(v_l)
            t^(k_l) partial_t^(j_l) F_(X^(alpha+u_l-v_l)).    (6)

Every right-side amplitude has degree d-r_l<d. The sign in b(d) is
positive: replacing it by b(-d) is wrong.

In fact b(0)=0. At alpha=0 every term of R vanishes in (5), because
|u_l|<|v_l|. Thus b(0) F_1=0, whereas F_1 has constant term one.
Let B>=0 be the largest nonnegative integer root of b. For d>B,
b(d) is nonzero, so (6) determines the image from lower-degree amplitudes.
This also reproves the finite multiplier criterion using the certificate
(3), without needing to expand a full integration-module presentation.

## 4. Extracting an affine cutoff, and then a homogeneous one

For the finite nonzero terms of R define

    C_rec=max(0, max_l (k_l-j_l)/r_l).                        (7)

Take C_rec=0 when R=0. Since differentiating a polynomial j times lowers
degree by j or kills it,

    deg(t^k partial_t^j f)<=deg f+k-j

when the result is nonzero.

Test the finitely many base series F_(X^alpha), |alpha|<=B, for polynomiality
by the earlier exact endpoint-indicial test. If any is nonpolynomial, P
fails the eventual mixed property and this construction makes no finite
cutoff claim. If they all pass, put

    A=max(0, max_(|alpha|<=B,F_alpha!=0)
                        (deg F_alpha-C_rec |alpha|)).         (8)

Induction on d in (6) gives

    deg F_(X^alpha)<=A+C_rec d

for all alpha: a right-side term has degree at most

    A+C_rec(d-r_l)+k_l-j_l <= A+C_rec d.

Finite linear combination then proves (1) for every Q.

For a pure-null P, F_1=1 has degree zero. Increase the slope, rather than
introducing an offset:

    C_pure=max(C_rec,
         max_(1<=|alpha|<=B,F_alpha!=0) deg F_alpha/|alpha|),   (9)

with the second maximum omitted if there is no such alpha.
Every base satisfies deg F_alpha<=C_pure |alpha|, including alpha=0.
The same induction proves (2). If an integer coefficient is preferred,
ceil(C_pure) is sufficient.

The distinction between affine (1) and homogeneous (2) is deliberate:
the proof of a zero offset uses the all-order pure-null seed.
Eventual pure vanishing alone is not substituted for that hypothesis.

Over an effective field, the b-function/lift, its finite integer-root list,
base polynomial series and rational maxima are all computable. The
construction either returns a certified nonpolynomial finite-degree marker,
or returns A,C and, when separately pure-null, C_pure. It is not a claim
that the full general solver has been implemented.

## 5. Exact finite certificates for the zero and isotropic linear polynomials

These examples verify signs and extraction; they are not a growing
low-degree research campaign.

For P=0 in standard Gaussian coordinates, take

    b(s)=s,             R=-sum_i partial_(X_i)^2.

Indeed theta exp(-|X|^2/2)=sum_i partial_(X_i)^2 exp(-|X|^2/2).
There is only the degree-zero base, C_rec=0 and C_pure=0.
Every F_Q is constant. Recurrence (6) is the Euler Gaussian moment
recurrence.

In circular coordinates E(ZW)=1, for P=W take

    Psi=tW-ZW,
    theta=-partial_Z Z-partial_W W,
    b(s)=s,
    R_0=-2 partial_Z partial_W+t partial_Z.                   (10)

Direct application shows (theta+R_0) exp(Psi)=0. The lower-weight terms
have degree drops two and one, with time degree shifts zero and one.
Thus B=0, C_rec=1 and the only base is F_1=1. The theorem returns
C_pure=1, which is sharp because F_(Z^d)=t^d.

A nonminimal certificate exercises roots, mixed normal ordering and time
derivatives. Set b_1(s)=s(s-1) and

    R_1=(theta-1)R_0+partial_Z(partial_t+partial_Z).

The last summand annihilates exp(Psi), since
(partial_t+partial_Z)exp(Psi)=0. Expanded in normal order,

    R_1=2Z partial_Z^2 partial_W
        +2W partial_Z partial_W^2+6 partial_Z partial_W
        -tZ partial_Z^2-tW partial_Z partial_W-3t partial_Z
        +partial_t partial_Z+partial_Z^2.                   (11)

Now B=1 and C_rec remains one. The base images are 1,t,0 and again give
C_pure=1. The checker verifies both exact exponential-amplitude identities,
the corresponding recurrences, and a mutation where b(-d) gives a false
identity. It also checks that zero or non-lowering certificates are rejected.
No minimality of either b-function certificate is asserted.

## 6. Boundaries and consequences for the research programme

A finite sufficient certificate consists of (3)--(4), the finitely many
base polynomials with their exact origin-indicial zero certificates, and
the rational slope inequalities in (7)--(9). It proves the all-Q bound
by a written induction. Checking finitely many mixed moments alone does
not establish that induction or replace a base-series identity certificate.

The coefficient in (2) cannot be fixed across all degrees. The retained
Gaussian shears

    S_j=W+Z^j T-Z^(2j+1)/2

are safe and satisfy F_W(t)=-t^(2j+1)/2. Thus any valid C_pure for S_j
is at least 2j+1. This is the existing all-degree barrier, fully compatible
with the new existence of a linear cutoff for each fixed safe P.
No common coefficient for all safe degree-D polynomials, optimal C_pure,
or coefficient-locus constructibility is established here.

The [quartic pole-image obstruction](FINITE_MULTIPLIER_MIXED_DECISION.md)
still prevents safety: one base image in any valid finite b-function
test must be nonpolynomial. A rank-one rational image is not enough.
The new result does not classify full Gaussian images or their pole and
irregular defects. Nor does it transfer to an ordinary-Laplacian GVC
statement without preserving the operator's separate hypotheses.

The general qualitative-versus-linear distinction is now closed internally,
and the constants are theoretically computable. Remaining priorities are
structural coefficient classification, uniform bounds at fixed n,D,
optimal cutoff geometry, full-image defects, valid operator transfers and
independent review. Numerical cubic optimization stays secondary and no
new cubic number is computed.
[Source audit](linear_mixed_cutoff_source_audit.json).

