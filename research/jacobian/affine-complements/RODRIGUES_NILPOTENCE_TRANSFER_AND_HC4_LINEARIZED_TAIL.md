# Rodrigues nilpotence transfer and the HC4 linearized tail

PF2RULED152. 17 November 2026.

The completed Image local certificate was tested against the ACTUAL HC4
normal recurrence using its actual length17 multiplication algebra.
The result is a concrete obstruction to an automatic transfer:

    J^2=0 and m^7=0,
    but the HC4 trace deformation has infinitely many nonzero
    normal coefficients in the SAME fixed ideal J.

The coefficient sequence and its finite propagation operator are exact.
This is not an inference from checking a longer prefix.

This does NOT disprove a Rodrigues-style certificate on a more restricted
HC4 trace quotient. No such quotient ideal and increasing-power certificate
are supplied here. The test identifies a linearized condition which that
quotient must satisfy before a positive transfer can work.

## 1. What the Image certificate actually provides

Use the preserved
[Rodrigues local certificate](../../moments/vanishing/RODRIGUES_LENGTH17_AND_GLOBAL51.md),
not the weighted initial quotient. Its actual completed local algebra H
has dimension17, maximal ideal m=(a,b), and

    dim(m^j)=(17,16,14,11,7,4,1,0), m^7=0,
    eta=b^6!=0, m*eta=0, eta^2=0.                 (1)

The multiplication matrices M_a,M_b in the original certificate_data.json
are used literally. They commute, reproduce the stated17-element cyclic
unit basis and annihilate the actual prefix11--14 residuals. The original
finite-row/Nakayama argument identifies this algebra with the full completed
local quotient, rather than merely its order-seven truncation.

Set J=(eta). It is a NONZERO square-zero ideal of the actual algebra.

For the Image residuals, the additional input is the written weight bound
5i+4j>=n, together with F^13 H=0 and finite remaining memberships.
That makes residual indices n>24 vanish. The global length51/selected-locus
exhaustion problem is separate and remains open; it is not imported as a
completed full Image theorem.

The HC4 test below uses H only as a COEFFICIENT ALGEBRA. It does not assert
that H is an HC4 trace-parameter quotient or transplant the moment equations
into the Hessian recurrence.

## 2. The HC4 recurrence survives Artinian base change

The retained universal normal identity
[HC4MYGJ2](../../foundations/HC4_MENG_YANG_GRAPH_OBSTRUCTIONS.md)
is

    D_(R+x^k U)-D_R
      = -4*L*N^3*k*(k+1)*x^(k-1)*U modulo x^k.   (2)

It is a polynomial identity with rational coefficients. Consequently it
holds over EVERY commutative Q-algebra with L,N units, including H.
All normal multipliers are units there. The same induction gives existence
and uniqueness of the formal graph for every plane trace and target
constant. This is a valid transfer of formal solvability to the finite
algebra; it gives no increasing ideal power by itself.

The exact algebraic formal model from
[HC4MYGAP1](../../hessian/hc4-polynomial-graphs/PROOF.md) lets us test the missing step
without generating an unbounded prefix.

## 3. A fixed-determinant quadratic trace deformation

Put U=1+x*y and retain

    A0=U^3*p+3*x*U^2*q,
    B0=y^2*U*(4+3*x*y)*p
        +(y+3*x*y^2*(4+3*x*y))*q,
    Psi0=L*A0^2+M*A0+N*B0.

Choose the polynomial potential

    F_tau=L*p^2+M*p+4*N*y^2*p+N*y*q
              +kappa*x*p+tau*x*y^2.               (3)

DIRECT calculation gives

    det Hess(F_tau)=(N*kappa)^2,                 (4)

independent of tau. The added term is cubic, not an affine term being
discarded. Thus this is a trace deformation with the SAME target Hessian
constant and the SAME parent L,M,N; varying the determinant is not the cause.

Let

    C=N*(2-3*x*y)-M*x^2-2*L*x^2*A0,
    D_tau=(Psi0-F_tau)/x.

The equation Psi(graph R_tau)=F_tau is exactly

    L*x^5*R_tau^2+C*R_tau+D_tau=0.                 (5)

Its implicit derivative at x=0 is 2N, a unit. Let R_0 denote its solution
at tau=0 and let Delta_0=C^2-4*L*x^5*D_0.
Use the square root with leading value 2N. Since

    D_tau=D_0-tau*y^2,
    Delta_tau=Delta_0+4*L*tau*x^5*y^2,

implicit differentiation gives the EXACT linearized response

    dR_tau/dtau at0
      = y^2/(C+2*L*x^5*R_0)=y^2/sqrt(Delta_0).    (6)

The plane trace changes only by tau*y^2/(2N).
For v2 L=1,M=13,N=2,kappa=1 the body trace is

    T0= p/4
       -3*(2*p+13)*(y*p+q)/4
       -7*y^3*p/2-6*y^2*q.                       (7)

Thus tau varies an actual QUADRATIC trace coefficient, with determinant4 fixed.

## 4. Use the ACTUAL Rodrigues multiplication algebra

Substitute tau=eta=b^6 in H. Because eta^2=0, Taylor expansion truncates
exactly after its first parameter derivative:

    V=R_eta-R_0=eta*y^2/sqrt(Delta_0).             (8)

Every coefficient of V belongs to J. There are no ignored higher eta terms.

On y=p=0 the discriminant is

    Delta_0=(2N-M*x^2)^2-24*L*N*x^3*q.

Extracting y^2 and then p^0 q^0 from (8), at v2, gives

    [y^2 p^0 q^0]V = eta/(4-13*x^2)
       = sum_(j>=0) (13^j/4^(j+1))*eta*x^(2j).   (9)

EVERY displayed coefficient is nonzero. Multiplication by a nonzero
rational scalar is a unit operation on H, and eta!=0 by actual multiplication.

The exact propagation operator on this projected coefficient state is

    T=(13/4)*Id_H.                               (10)

It is a UNIT, not a nilpotent operator. The coefficient at x^(2j)
is T^j(eta/4). This finite-state calculation is an ALL-ORDER certificate
of failure of increasing ideal powers.

In particular there is no finite c>0 for which all normal coefficients
v_k=[x^k]V satisfy

    v_k in J^ceil(k/c).

For k>c the right side is zero; infinitely many v_(2j) have the nonzero
projection in (9). Nor can v_k be forced into m^ceil(k/c): for k>6c
the right side is m^7=0, with the same contradiction.

More sharply, if Tail_D(V) denotes the ideal in H generated by ALL plane
monomial coefficients of v_k for k>D, then

    Tail_D(V)=J for EVERY finite D>=0,
    Tail_D(V)^2=0, but V is NOT polynomial.       (11)

One late even coefficient already recovers eta by a rational unit;
all coefficients are in J by (8). Thus a stable nilpotent tail ideal
is not a termination certificate unless its powers actually increase.

## 5. The obstruction is not removed just by setting M=0

The same response has an all-order coefficient formula for every L,N unit:

    [x^(3j) y^2 p^0 q^j]V
      = eta/(2N) * binomial(-1/2,j)*(-6L/N)^j !=0
      for EVERY j>=0.                            (12)

This follows by expanding

    ((2N-M*x^2)^2-24LN*x^3*q)^(-1/2)

in q. Its first q^j term occurs at x order3j and has the displayed
coefficient. The rational binomial is nonzero in characteristic zero;
its product with the units L,N is a unit multiplying eta.

Hence M=0 removes the simple geometric sequence (9) at q=0, but does
not remove the infinite first-order normal response over the free
source variable q. There is no weight estimate on trace-ideal degree
alone analogous to the Rodrigues residual estimate on this family.

For v2 the obstruction is stronger still: the SAME plane monomial y^2
and the SAME trace-ideal layer occur at arbitrarily large normal order.

## 6. What this test does and does not rule out

The test disproves the inference

    actual finite local algebra + nilpotent maximal ideal
      + triangular HC4 formal recurrence
      => increasing ideal-power containment.

It also refutes any proposed universal HC4 certificate whose trace quotient
admits this parameter specialization and whose claimed power bound covers
this deformation. Even AFTER subtracting the body, the first-order tail
does not gain ideal order.

It does NOT refute the user's proposed strategy on a genuinely restricted
remaining trace quotient S/I. Such a certificate must name I and prove
that its actual relations either remove this base point or annihilate
the offending infinitesimal direction. The original Rodrigues ideal and
its nilpotence exponent cannot simply be declared to be HC4's I or J.

The applicability boundary is important. The body (7) is the retained
nonrational formal graph. At v2 its trace has degree4 and its first normal
coefficient has plane degree5, so the complete first jet has degree6.
It is ALREADY outside a degree-at-most-five polynomial-graph ansatz.
We therefore do not advertise (9) as a counterexample on that smaller
remaining quotient. It is an exact transfer guard for candidate certificates.

Also, the standard marked polynomial-graph collision route is already
closed by HC4MYGALL1. The unmarked polynomial-termination problem and
global minimal-counterexample coverage remain separate.

## 7. A concrete finite certificate for a positive transfer

Here is one SUFFICIENT way to prove the increasing-power statement.
It is not claimed necessary, and its hypotheses are not yet constructed
for the remaining HC4 quotient.

Let A be the actual candidate coefficient algebra, J^N=0.
First establish that the reduced graph is polynomial through normal
degree D. Represent the ENTIRE remaining recurrence, including its
actual multiplication, by a finite A-module E with A-linear transition T,
initial state v in J*E and A-linear output ell, so

    a_(D+k)=ell(T^(k-1)*v), k>=1.

A finite identity

    T^c(E) subset J*E                              (13)

then gives, by A-linearity,

    a_(D+k) in J^(1+floor((k-1)/c))
                 =J^ceil(k/c).                  (14)

Therefore a_(D+k)=0 for k>c*(N-1), with a uniform degree bound.
For a free module (13) is a finite matrix-membership test; equivalently,
the reduced transition has nilpotent c-th power. Every nonlinear product
and every coefficient-state closure used to construct E must be justified.
A finite truncation is not automatically an exact state module.

This separates the two missing certificates:
actual nilpotence of J, and nilpotence of propagation MODULO J.
Rodrigues has its independent weight argument for the latter role.
The HC4 response (10) instead has T mod J a nonzero scalar unit.

A positive control using the same ACTUAL H makes (13) concrete.
On E=H^2 take the block matrix

    T_plus = [0 M_b; Id 0].

It has T_plus^2=diag(M_b,M_b), so c=2 and N=7.
With initial state (b,0) and output the sum of the two components,
the sequence is

    b,b,b^2,b^2,...,b^6,b^6,0,...

and has precisely the bound (14), terminating after12 terms.
This is a CONTROL for the finite certificate, not the HC4 transition.

## 8. Next gate and preserved scope

The next useful calculation is a finite nilpotence test on the reduced
linearized propagation for a precisely identified HC4 trace quotient,
before generating more normal coefficients. If a unit mode such as (10)
survives with a nonzero observed tail state, the desired ideal-power bound
cannot hold there. If it vanishes,
one must still certify the full filtered recurrence or its weight bound,
not just its tangent operator.

No actual HC4 obstruction quotient I, uniform J^N certificate, or full
increasing-power containment on the remaining locus is claimed.
The positive Rodrigues-style HC4 application remains OPEN.

[PF151](NORMAL_FORM_COVERAGE_AND_COMMON_LEADING_DIVISORS.md)'s global
coverage obligations remain open. All full-normalization and ramified-prime
requirements, the dimension-three integrality controls and
[PF142](PARTIAL_KUMMER_PLANE_CHARTS_AND_CYCLIC_EQUIVARIANT_FACTORS.md)'s
mixed minimal retirement are preserved. No local coefficient algebra is
identified with the original finite normalization S or Lambda.

The written proof is (2)--(14), with actual multiplication in (1).
Finite symbolic controls are support, not an extrapolation of a prefix.
No new external theorem import or registered predecessor checker replay
is used. Independent whole-proof replay, formal verification, external
review and novelty priority are unclaimed. Full HC4/JC2 is OPEN.
