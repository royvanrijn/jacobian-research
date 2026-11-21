# A mixed cubic primitive Image family in every balanced degree divisible by three

**Theorem OffcenterCubicEveryMultipleOfThree.** For every integer j>=1,
the literal integer original polynomial F_j below has balanced degree
D=3j+3, is radial-primitive, and is not a scalar proper power over any
characteristic-zero extension. For every positive integer m,
E2(F_j^m)=0 and E2(Z F_j^m)=(mD+2)!(-8)^m kappa_m is nonzero.
The SAME original multiplier Z works for every j and every m.

This note gives an explicit integer two-pair Image witness in every
balanced degree D=3j+3, j>=1, with one fixed original multiplier Z
surviving at EVERY positive order. It uses the conditional off-center
formal endpoint calculation, not a claim that all pure cubics satisfy
its extra invariant. Existing degree-six/even witnesses and older
odd-composite existence families are antecedents; no priority claim is made.

Let K have characteristic zero and put
A=K[x,y,t]/(t^2+2xy-1), with the rational interval/constant-term functional
Lambda as in the retained two-pair normalization. Set

    Q=1+x,  s=-1+x^j,  E=Q^4,  h=2/Q,  j>=1,
    P_j=Q^4(t-s)^3-4Q^2(t-s),
    f_j=P_j/x.                                         (N1)

## 1. Literal regularity and exact degree

At x=0, s=-1 and Q=1, so

    P_j(0,t)=(t+1)^3-4(t+1)
             =(t^2-1)(t+3).

Thus the following is an INTEGER polynomial representative in A:

    G_j=-2y(t+3)+[P_j-((t+1)^3-4(t+1))]/x.              (N2)

The bracket is coefficientwise divisible by x. If B_t=t+1, an
equivalent completely expanded regular formula is

    G_j=(x^3+4x^2+6x+4)B_t^3-4(x+2)B_t
        -3Q^4 x^(j-1)B_t^2+3Q^4 x^(2j-1)B_t
        -Q^4 x^(3j-1)+4Q^2 x^(j-1)-2y(t+3).

Thus no symbolic division remains in the original construction.
Its cubic numerator coefficients are

    C=-Q^4 s^3+4Q^2 s,   B=3Q^4 s^2-4Q^2,
    D=-3Q^4 s,          E=Q^4.                         (N3)

Hence C0=-3,B0=-1,D0=3,E0=1, and the two regularity conditions hold.
The exact normal-form degree formula gives

    deg_A f_j=max(3j+3,2j+4,j+5,6)=3j+3.               (N4)

For j=1 all four bounds equal6; for j>=2 the first is largest.
The x-only top coefficient is -1 in every case. No leading quadric
can remove the normal-form top, so these are least degrees, not
padded homogenization degrees. Also J^3 f_j=6x^2 Q^4 is nonzero
and J^4 f_j=0.

## 2. All pure moments and EVERY marked order

The numerator has the exact depressed form
E[(t-s)^3-h^2(t-s)]. It lies in the off-center I=0 branch because
C0=3B0 and E0 is nonzero. Use the formal square root h with h0=2.
Then

    A0=-E h^3=-8Q,    L=h/2=1/Q,
    L=1-x A0'/A0.                                    (N5)

For each m>=1, set

    kappa_m=integral[0,1] [v(1-v^2)]^m dv
            =2^m m! / product_(r=0)^m(m+2r+1).         (N6)

This is a nonzero rational number. The formal substitution v=(t-s)/h
has endpoints0 and1 modulo x. Since its integrand vanishes to order m
at both endpoints, the rational polynomial interval average obeys

    (1/2) integral[-1,1] P_j(x,t)^m dt
      =kappa_m (-8)^m Q^(m-1) modulo x^(m+1).           (N7)

This congruence includes both coefficient orders needed below:
it gives

    Lambda(f_j^m)
      =kappa_m (-8)^m [x^m]Q^(m-1)=0,

    Lambda(x f_j^m)
      =kappa_m (-8)^m [x^(m-1)]Q^(m-1)
      =(-8)^m kappa_m !=0.                            (N8)

Thus f_j is pure and unsafe, with the SAME sphere marker x for
every positive m and every parameter j. The all-order conclusion
comes from (N7), not from a vanishing finite prefix.

## 3. Literal primitive original lift and absolute nonpower

In the original polynomial ring K[xi1,xi2,z1,z2], define

    R=xi1*z1+xi2*z2,  Z=xi1*z2,
    W=2xi2*z1,       T=xi1*z1-xi2*z2.

Then T^2+2ZW=R^2. For D=3j+3, define the explicit homogeneous lift

    F_j=R^D G_j(Z/R,W/R,T/R).                          (N9)

Because the ordinary degree of G_j is D, this is an integer polynomial,
homogeneous of degreeD separately in xi and z. Its ordinary total
degree is2D.

At the ORIGINAL point (xi1,xi2,z1,z2)=(1,0,0,1), one has
R=T=W=0,Z=1. The top x-only coefficient in (N2) is-1, so

    F_j(1,0,0,1)=-1.                                  (N10)

Therefore R does not divide F_j, over K or any field extension.
This is direct radial primitivity; no padding or generic point is used.

For absolute scalar nonpower, take the ORIGINAL polynomial curve

    xi1=z2=1,   xi2=(1-tau)/2,  z1=(1+tau)/2.

It has R=1,Z=1,T=tau,W=(1-tau^2)/2. At x=1,
Q=2 and s=-1+1^j=0, for every j. Consequently

    F_j|curve=16tau^3-16tau.                           (N11)

The zero tau=0 has derivative-16. A nonzero scalar times a kth power,
k>=2, cannot specialize to a nonzero polynomial with a simple zero.
This rules out every scalar proper-power factorization in the FULL
original ring over EVERY coefficient-field extension.

## 4. The exact original Image contraction

The retained two-pair balanced contraction identity gives

    E2(F_j^m)=(mD+1)! Lambda(f_j^m)=0,

    E2(Z F_j^m)=(mD+2)! (-8)^m kappa_m !=0,
         for EVERY m>=1.                              (N12)

All factorials are invertible nonzero rational scalars in characteristic
zero. The original Image-kernel identity therefore puts every positive
power F_j^m in the two-pair Image sum, while the ONE fixed original
multiplier Z survives at every positive order. This is an actual
primitive absolutely nonproper-power Image witness in each
D=6,9,12,15,... .

The factorial in (N12) uses the actual new degreeD. One does not reuse
a previous seed's radial factorial. The distinction between angular
moments, literal original primitivity/nonpower, and the final original
Image identity is essential.

## 5. Exact scope

The new deduction strengthens the marker behavior in the displayed
multiple-of-three family. It does not settle balanced degree5 or7.
It does not assert that arbitrary mixed cubic phases satisfy I=0,
that arbitrary transported markers stay unchanged, or that every
balanced degree has an all-order fixed-Z witness. The degree9 and15
instances already had existence witnesses in older composite-degree
families; the explicit family and EVERY-order marker formula here are
the statements being proved.

The all-order endpoint argument from the conditional mixed-cubic branch
is reproduced directly in (N5)--(N8). The earlier mixed note is an
antecedent, not an indispensable theorem input for this specialization. The original balanced contraction
and Image-kernel identities are retained narrow inputs. Neither a
prime-density theorem nor an automorphism-group classification is used.
No external review, formal verification or historical priority is asserted.


## 6. Exact source and verification boundaries

The only retained theorem inputs are the arbitrary-balanced-degree
normalization (3)--(4) in TWO_PAIR_IMAGE_CUBIC_STRATIFICATION.md and
the direct characteristic-zero Image-kernel identity (2.1)--(2.2) in
TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md. The complete cubic safety
theorem and old existence families in those source files are NOT used.
The elementary beta/product formula, endpoint congruence, literal
normal form, radial point and nonpower curve are proved in this note.

The exact standard-library checker verify_offcenter_cubic_all_order.py
reconstructs the literal original four-variable lifts before contracting
the nine j=1,2,3 / m=1,2,3 controls. Its symbolic profile/divisibility,
coefficient and normalization controls authenticate the formulas; the
written endpoint argument proves the universal j,m quantifiers.
The checker also authenticates exactly three own roles and the two
raw source roles above. Neither424/all-Aut, FPYZ, the full mixed-cubic
theorem nor an old witness-existence principal is a premise.

```bash
python3 research/foundations/extended-geometry/verify_offcenter_cubic_all_order.py \
  --repo-root . --output /tmp/offcenter-cubic-all-order.json
```
