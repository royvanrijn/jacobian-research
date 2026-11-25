# Mixed cutoffs, finite moment jets and certified infinite tails

Let K be a characteristic-zero field, A a filtered polynomial algebra,
and L:A->K a linear functional. Fix a known pure phase P0 and a
nonnegative rational constant C such that, for every original marker Q,

    L(P0^n)=0 for n>=1,
    L(Q P0^n)=0 whenever n>C deg Q.                 (1)

The marker degree in (1) is the original spatial degree. A cutoff after
an uncontrolled change of marker filtration does not supply this hypothesis.
Purity and the unrestricted mixed cutoff are premises; neither follows
from the finite moment computations below.

## 1. The finite-jet theorem

Let R=K[u1,...,us], m=(u1,...,us), and let H in m R[X] have spatial
degree at most d. The coefficients of H may depend polynomially on u;
in particular H=sum ui Qi is allowed. Set

    P=P0+H, F_n(u)=L(P^n),
    I_N=(F_1,...,F_N), I_infinity=(F_n:n>=1),
    B_r=floor((Cd+1)r).                            (2)

The convention I_0=0 is used. For every integer r>=0,

    F_n in m^(r+1) for every n>B_r,
    I_infinity+m^(r+1)=I_(B_r)+m^(r+1).             (3)

The same statements hold after passing to a quotient local slice,
localizing at the chosen origin, completing, or pulling back along a
formal parameterization H in m K[[u]][X] of bounded spatial degree d.
No transversality or reducedness is required.

For a linear coefficient parameterization, the homogeneous perturbation
degree-j part of F_n is exactly

    binom(n,j)L(H^j P0^(n-j)).                     (4)

It vanishes if n-j>Cdj. Indeed deg_X H^j<=dj, and (1) applies to
every parameter coefficient. For a nonlinear parameterization, the
binomial term with H^j has parameter order at least j. Modulo m^(r+1)
only j<=r remain; if n>(Cd+1)r, every one of those terms vanishes by
the same inequality. The j=0 term vanishes by purity. This proves (3)
before any base change, and the same finite argument works coefficientwise
in the formal case. Extending the resulting identities proves the
quotient, local and completed versions.

In the ternary Gaussian and normalized two-sphere degree-at-most-three
settings, the retained [Gaussian Theorem `thm:cubic`, equation `eq:cutoff`](../../../results/gaussian-boundary/manuscript.tex)
and [sphere Theorem `thm:sphere`, equation `eq:cutoff`](../../../results/cubic-sphere-image/main.tex)
supply C=3. Degree-at-most-three
perturbations therefore give

    I_infinity+m^(r+1)=I_(10r)+m^(r+1).             (5)

These are local jets at a KNOWN pure phase, not a universal finite prefix
test for purity at arbitrary coefficient points. The Gaussian and sphere
functionals, and their spatial degree conventions, remain separate.

## 2. A finite certificate for an infinite tail

Work in a local polynomial slice or its completion at the chosen origin.
Suppose an exact certificate proves

    m^nu subset I_N, nu>=1.                        (6)

Put B=floor((Cd+1)(nu-1)). By (3), every F_n with n>B already belongs
to m^nu and hence to I_N. Consequently

    I_infinity=I_N
      iff F_n in I_N for N<n<=B.                   (7)

If N>=B the interval is empty. For C=d=3 and nu=7, B=60. This bound
is conditional on the actual nilpotence certificate (6); knowing the
radical, the reduced zero set, or tangent-space equations is insufficient.

There is a wholly finite way to certify (6). It is enough to verify

    m^nu subset I_N+m^(nu+1)                       (8)

in the nu-jet algebra. Let J=I_N, including any slice relations in the
ambient presentation, and let M=(m^nu+J)/J. This is finitely generated.
Condition (8) gives M=mM, so Nakayama's lemma yields M=0 in the local
ring. Thus (6) follows, and also holds in the completed ring.
Explicitly, for a finite generating vector w of M write w=Bw with
every entry of B in m. The adjugate identity gives det(I-B)w=0;
det(I-B) is one modulo m and hence a local unit. This proves M=0.

This jet-to-ideal step is exactly classical Nakayama, in the finite-module
form of [Stacks, Lemma 10.20.1(2), tag 00DV](https://stacks.math.columbia.edu/tag/00DV),
applied to the displayed M in the local ring. The adjugate proof above
is its determinant proof. Separately, if the full moment ideal in a
finite coefficient polynomial ring has radical m, then
[Hilbert Nullstellensatz, tag 00FS](https://stacks.math.columbia.edu/tag/00FS),
[Noetherianity, tag 00FM](https://stacks.math.columbia.edu/tag/00FM), and
[the radical-power lemma, tag 00IM](https://stacks.math.columbia.edu/tag/00IM)
give some finite prefix containing some power of m. They do not select
N or nu, certify a chosen jet identity, or provide the mixed-degree
bound (3). Those inputs and the intervening memberships remain the
specific content of the certificate here.

To verify (8), retain for EVERY parameter monomial u^alpha of total
degree nu a polynomial identity

    u^alpha=sum_i A_(alpha,i) F_i+sum_j B_(alpha,j) g_j
                                             modulo m^(nu+1),       (9)

where F_i are the first N actual moments and g_j the slice relations.
Multipliers may have parameter degree up to nu; using only constant
linear combinations would not compute the truncated ideal. For each
intervening F_n in (7), it then suffices to retain the analogous identity
modulo m^nu. Its omitted remainder lies in I_N by (6).

This certificate is LOCAL. For example J=(u^nu(1-u)) passes (8) at
u=0, but u^nu is not in J globally in K[u]. The factor 1-u is a unit
only after localization or completion. A global statement requires
actual global ideal-membership identities, not just these jets.

## 3. The exceptional shear: exact original images and sharp jet order

Use a Witt pair Z,W with Gaussian E[ZW]=1 and an independent ordinary
Gaussian Y with E[Y^2]=1. Let

    q=2ZW+Y^2, P*=W+ZY-Z^3/2.

For normalized sphere moments on q=1, a homogeneous polynomial of degree
2a has sphere mean equal to its Gaussian mean divided by (2a+1)!!;
odd means vanish. Both functionals are invariant under the polynomial
shear

    phi(Z)=Z, phi(Y)=Y-Z^2, phi(W)=W+ZY-Z^3/2,
    phi^-1(Z)=Z, phi^-1(Y)=Y+Z^2,
    phi^-1(W)=W-ZY-Z^3/2.                          (10)

For completeness, its generating derivation
delta=Z(-Z partial_Y+Y partial_W) kills q and has ordinary divergence
zero. It is locally nilpotent, and Gaussian integration by parts gives
E[delta A]=0. It raises homogeneous degree by one, so the displayed
Gaussian-to-sphere normalization gives L_sphere[delta A]=0 as well.
Exponentiating the finite derivation proves both invariances. Since
phi(W)=P*, both phases are pure and

    Phi_A(t)=L(A exp(tP*))=L(phi^-1(A) exp(tW)).      (11)

For a monomial Z^a W^b Y^(2c) on the right, its period is zero when
a<b. Otherwise it is the single monomial

    a!/(a-b)! (2c-1)!! t^(a-b)                     Gaussian;
    a!/(a-b)! (2c-1)!! t^(a-b)/(2a+2c+1)!!         sphere.            (12)

Odd Y powers give zero. Thus the original marked image of degree<=D
is contained in K[t]_(degree<=3D). For the reverse inclusion, the
inverse shear of Z^iY^jW^k has unique top spatial term

    (-1/2)^k Z^(i+2j+3k).

It gives a nonzero top period coefficient at degree n=i+2j+3k, equal
to (-1/2)^k in the Gaussian case and (-1/2)^k/(2n+1)!! on the sphere.
Every integer n from 0 to 3D is represented with i+j+k<=D. Selecting
one original marker for each n yields a triangular basis. Therefore,
for BOTH functionals independently,

    {Phi_A:deg A<=D}=K[t]_(degree<=3D).             (13)

This proves the sharp original mixed cutoff C=3 at P*. It also proves
that the moment Jacobian in all cubic coefficient directions has rank
10: dF_n(P*)[B]=n L(B P*^(n-1)), rows 1,...,10 are independent and
all subsequent rows vanish. Raw ternary polynomials of degree<=3 have
dimension 20, so their tangent kernel has dimension 10. For sphere
functions the quotient by (q-1)K[Z,W,Y]_(degree<=1) has dimension 16;
the Jacobian still has rank 10 and its tangent kernel has dimension 6.

The factor 10 in the MOMENT jet cutoff is sharp. For H=epsilon W^3,
the coefficient of epsilon^r in F_(10r) is

    (10r)!/r! (-1/2)^(3r)                          Gaussian;
    (10r)!/r! (-1/2)^(3r)/(18r+1)!!                sphere,           (14)

which is nonzero. This follows from the unique top term Z^(9r) in
phi^-1(W^(3r)). It does NOT prove that 10r is the least ideal-prefix
cutoff: lower moments can already imply that equation. Indeed in this
same one-parameter family F_2=-6 epsilon for Gaussian moments and
F_2=-2 epsilon/35 for sphere moments.

## 4. Exact controls and the verifier's scope

The genuine linear coefficient slice H=epsilon(Y-Z^2) supplies a small
nonreduced test. By (10), its moments equal those of W+epsilon Y:

    F_(2j)=(2j-1)!! epsilon^(2j)                    Gaussian;
    F_(2j)=epsilon^(2j)/(2j+1)                      sphere;
    F_(2j+1)=0.                                   (15)

Thus I_2=(epsilon^2), with exact nilpotence index two. Since d=2,
the bound in (7) is 7; all intervening checks succeed.

A negative control shows why those checks matter. Take
H=epsilon(Y^2-1) for Gaussian moments, or
H=epsilon(Y^2-1/3) for sphere moments. In either case I_2=(epsilon^2),
but F_5 has a nonzero linear term: 120 epsilon in the Gaussian case,
8 epsilon/63 in the sphere case. A verifier must reject I_2 stabilization.

The [companion implementation](moment_jet_tail_certificates.py) constructs exact shear moment jets using
(10)-(12), exports the finite identities (9), and verifies them from the
original phase input. Its executable coefficient domain is QQ; the
written theorems hold over every characteristic-zero field. The backend
is bound to the known pure P*, and is not a purity test for an arbitrary
base phase or a varying covariance. Its output proves local and completed ideal-tail
stabilization using (1)-(9). The unrestricted algebraic theorem and
filtration proof above are written arguments, not inferences from finite
controls. Independent replay, formal verification and external review
are not claimed.

## 5. A genuine index-seven sphere slice and its complete moment ideal

The [constant-profile calculation](../../../results/cubic-sphere-image/finite.tex),
Proposition `prop:finite-constant-profile` and equation
`eq:constant-profile-triangular`, gives the following polynomial graph
through P*. Fix its frame coefficients a=d=1,c=0 and put

    P(u,v,e)=W+(u+vZ)Y^2+(Z+eZ^2)Y
       -u/3-(2u^2+3v)Z/15
       -4u(31u^2+54v)Z^2/945
       -[1/2+2312u^4/14175+188u^2v/525+8v^2/75]Z^3.                 (16)

The first four normalized sphere moments vanish identically. In
QQ[u,v,e], the FULL moment ideal of this specified family is

    I_infinity=I_8
      =(u+e-1184e^5/155925,
        v+e^2/2-6022904e^6/638512875, e^7).                         (17)

In particular its quotient is exactly QQ[e]/(e^7), its length is seven,
and the maximal ideal has exact nilpotence index seven. The class e^6
is nonzero. These conclusions persist after localization, completion,
or characteristic-zero field extension.

At the parameter origin the least local prefix is even smaller:

    (I_infinity)_m=(I_7)_m,  F_7 not in (I_6)_m,                    (19)

and the same holds in the completed ring. Thus seven is the EXACT least
local prefix on this family. This does not assert that seven suffices
globally or that eight is the least global prefix.

This is a genuine closed immersion of the three-dimensional parameter
space into the cubic coefficient space: u,v,e are recovered as the
coefficients of Y^2, ZY^2 and Z^2Y. No parameter has been replaced by
its seventh power. It is not asserted to be transverse or versal to
the entire ambient cubic moment scheme. The length-seven conclusion
concerns this graph slice, rather than the full 16-dimensional sphere
coefficient space.

Here is the exact finite certificate for (17). Let N_5,...,N_8 be the
primitive integer numerators of the residual moments in the cited
constant-profile proof. Their nonzero ratios to the actual F_i are,
in order,

    9823275/8, -13408770375/8,
    273648375/32, -1465387048125/256.                               (18)

The [literal witness packet](sphere_pstar_constant_profile_tail_certificate.json)
retains seven intermediate generators, polynomial lifts from the N_i
to them, polynomial lifts from those generators to the three generators
in (17), and reverse lifts of every N_i from (17). Every identity is
checked over QQ against the original moment polynomials. Thus equality
of these finite ideals does not depend on a claimed Groebner flag,
dimension report, or set-theoretic radical calculation.

For a short structural explanation, the intermediate generators include

    e^3+2ve,  ve^2+2v^2,  v(u+e),
    4736v^2e-155925(u+e),
    96366464v^3+638512875e^2+1277025750v.

The first two give ve=-e^3/2 and v^2=e^4/4. The fourth gives the stated
expression for u; the fifth then gives the expression for v, and the
third forces e^7=0. Conversely substituting these u,v in all seven
intermediate generators gives zero modulo e^7. The retained lifts make
this an identification with I_8, and establish the asserted surviving
class e^6 without inferring nilpotence from the reduced locus.

Finally, the packet checks EVERY moment 1,...,60 in this quotient by
a separate inverse-shear/binomial contraction. All remainders vanish.
The verified global containment m^7 subset I_8 and the general tail
theorem in Section 2 then prove (17) at all orders. The local jet
certificate is also replayed separately by (9), using prefix seven:
its nilpotence minimality dual annihilates the entire six-jet ideal
and detects a degree-six monomial. It checks all 53 intervening moments
8,...,60 modulo m^7. This proves the equality in (19). A further retained
dual on the seven-jet algebra pairs actual F_7 to one and annihilates
every monomial multiple of each of F_1,...,F_6. It proves the nonmembership
in (19), including after localization or completion, and establishes
the least-prefix assertion without guessing from a failed computation.

There is a stronger independent bound specific to (16). Its perturbation
has no W; after the inverse shear it has Z-degree at most five. In the
length-seven quotient its seventh power is zero. A term Z^aY^b in its
j-th power contributes only to moment n=a+j<=6j, with j<=6. Hence
moments beyond 36 already vanish there. The general threshold 60 is
valid without this special structural optimization.

## 6. Replay and assurance

The [general local checker](verify_moment_jet_tail_certificates.py)
replays the [retained jet certificates](moment_jet_tail_certificate_examples.json)
and records its [receipt](moment_jet_tail_certificates_verification.json).
The [global slice checker](verify_sphere_pstar_constant_profile_tail.py)
checks the literal global lifts and original sphere contractions, with
its own [receipt](sphere_pstar_constant_profile_tail_verification.json).
From the repository root run

    .venv/bin/python research/moments/vanishing/verify_moment_jet_tail_certificates.py
    .venv/bin/python research/moments/vanishing/verify_sphere_pstar_constant_profile_tail.py

For another rational cubic perturbation at P*, the sparse instance format
in the retained examples can be supplied to the implementation with
`--produce --output certificate.json`; replay that output without
`--produce` to check every witness. A rejection at an intervening moment
is an actual finite obstruction to the requested prefix, not a failure
of the mixed-cutoff theorem.

The local checker includes deliberately false prefix-stabilization and
corrupted-certificate controls. Its local conclusion does not become
global merely because all finite jet identities pass. The global result
(17) has the additional exact polynomial ideal lifts that justify it.
These are internal proofs and exact computations, with no independent
replay, formal verification, external review, publication or priority
claim. General Gaussian safety implies purity and a least universal
pure-moment prefix for all ternary cubics remain open.
