# Full complex quintic graph exclusion by a cubic-nilpotent residual algebra

PF2RULED156. 17 November 2026.

For the FIXED v2 Meng--Yang parent, no polynomial graph of total degree
at most five has nonzero constant pulled-back Hessian determinant over
any characteristic-zero field. This closes BOTH conjugate complex trace
classes retained by [PF155](QUINTIC_GRAPH_COVERAGE_AND_TWO_CONJUGATE_TRACE_CLASSES.md),
and therefore the entire degree-five graph branch for this parent.

The residual proof is a genuine application of actual coefficient algebra
and uniform nilpotence: a specified residual first/second-jet quotient is
finite free of rank three over its remaining coefficient base, its actual
obstruction ideal (c) satisfies (c)^3=0, and the next immutable coefficient
has an explicit finite inverse. An explicit three-generator polynomial
identity proves that the complete residual gate is the unit ideal.

This is not a normal-tail increasing-power theorem or a proof of HC4/JC2.
Degree five means degree of the GRAPH, not the degree of its potential.
Other graph degrees, other embeddings and global minimal-plane-counterexample
coverage retain their existing open scopes.

## 1. The complete residual trace, with all lower coefficients

PF155 proves that every v2 degree-at-most-five graph candidate has trace

    T=D y^3p+h y^3+i y^2p+e yp^2+f p^3+rho y^2q
       +a pq+b yq+c q^2+d p^2+g yp+j y^2+k p+l q+m y+n,             (1)
    D=(16rho+89)/2,
    Q(rho)=160rho^2+1968rho+6021=0.                               (2)

Both roots of Q are retained simultaneously. No coefficient among
h,i,e,f,a,b,c,d,g,j,k,l,m,n is discarded here. The nonzero desired
constant determinant is denoted C, distinct from c=[q^2]T.

The exact plane identity is

    D_R|x0=F_plane(T)-64S, S=R_x|x0.

For (1), its ENTIRE part above plane degree four is

    [degree>4] F_plane(T)=16Q(rho)y^4p.                           (3)

Thus the actual residual first-jet trace quotient has base

    K0=Q[rho]/(Q),
    B=K0[h,i,e,f,a,b,d,g,j,k,l,m,n,C,C^(-1)],

and is B[c]. In this algebra S=(F_plane(T)-C)/64 has degree at most four.
The polynomial Q is irreducible over Q, with discriminant 576*34;
K0 is the quadratic field Q(sqrt(34)). This construction is an explicitly
specified quotient of the actual first-jet coefficient ring. It is not
an identification of a field-valued radical with the unreduced full
candidate scheme.

## 2. Two complete second-jet coefficients

Set R=T+xS for this calculation and write D_R for its actual Hessian
determinant. Define

    Wc=[x y^4 q^2]D_(T+xS), Wp=[x y^5 p]D_(T+xS).

For the displayed polynomial representatives, take the degree-at-most-four
part of S. By (3) this represents the exact forced jet in B[c].

The universal k=2 normal multiplier is -192. An arbitrary x^2U allowed
by total degree five has deg U<=3, and changes [x]D by -192U.
Both displayed plane monomials have degree six, so they cannot be repaired.
Terms of normal order >=3 are invisible at this order. Therefore every
candidate degree-five graph must satisfy

    Wc=Wp=0.                                                     (4)

Exact extraction from the actual parent potential gives, with EVERY
coefficient in (1) still present,

    alpha=12582912=3*2^22,
    Wc=alpha*c^3,                                                (5)
    L(rho)=2898rho+16821,
    U(rho)=(96/25)L(rho),
    B1=196608(i-8b)(e-4a-rho-21/2),
    Tq=-64e+(704/5)rho+16704/25,
    Wp=U+c B1+Q*Tq.                                              (6)

In particular both witnesses are independent of f; no top-kernel
restriction on the residual classes was needed to remove it. The full
witness polynomials and their exact decomposition are in the
[residual witness certificate](sources/quintic-residual-closure/residual_witness_certificate.json).

For reproducibility the checker expands the parent only through x^3,
forms the two Hessian jets H0,H1, and extracts the two coefficients of
trace(adj(H0)H1). It retains full f and all lower parameters. The bounded
source-jet multiplication is an exact coefficient extraction, not an
extrapolation of normal recurrence coefficients.

## 3. Actual cubic nilpotence and a finite inverse

Before imposing Wp, the exact residual gate is

    A_pre=B[c]/(Wc)=B[c]/(c^3).                                  (7)

It is FREE of rank three over B, with basis 1,c,c^2. Actual multiplication
by c is

    M_c=[0 0 0; 1 0 0; 0 1 0],
    M_c^3=0, M_c^2!=0, J=(c), J^3=0.                            (8)

This uniform identity holds over the entire remaining parameter base.
There is no initial-algebra replacement or hidden passage to a reduced
quotient. The nilpotent direction and its square are genuinely nonzero
before the Wp equation is imposed.

The unit residue in (6) has an explicit inverse modulo Q. Put

    v(rho)=(25/96)((368/76545)rho+1328/42525),
    E=529/6075.

Direct polynomial multiplication gives

    vU=1+E Q.                                                    (9)

Consequently in A_pre, Wp=U+c B1 is a unit with inverse

    Wp^(-1)=v(1-t+t^2), t=v c B1.                               (10)

Indeed t^3=0 by ACTUAL multiplication (8). This inverse is finite and
uniform in every lower coefficient, including f and the determinant C.
It proves

    A_pre/(Wp)=0.                                                (11)

Here nilpotence certifies an inverse for a finite repair-invariant
obstruction; it does not establish increasing normal ideal powers.

## 4. Explicit polynomial unit-ideal certificate

The scheme-level contradiction is also a small exact identity BEFORE
quotienting. Equations (5)--(6) and (9) imply

    v(1-t+t^2) Wp
      -(v^3 B1^3/alpha) Wc
      -(E+v Tq)(1-t+t^2) Q=1.                                   (12)

All multipliers are polynomials with rational coefficients. The
[unit-ideal certificate](sources/quintic-residual-closure/residual_unit_ideal_certificate.json)
retains them literally, and the checker verifies their sum equals one.
For example (1-t+t^2)(1+t)=1+t^3, with

    t^3=(v^3 B1^3/alpha)Wc.

Thus (Q,Wc,Wp)=(1) over the full residual polynomial coefficient ring,
without any assumptions on lower coefficients or on which conjugate
root of Q is used. This excludes the two PF155 classes, including their
nilpotent thickenings in the specified residual quotient.

## 5. An actual finite algebra control with nonzero c^2

Fix all free coefficients except c to rational values with i=1,
a=b=e=f=h=d=g=j=k=l=m=n=0 and C=1. The quotient BEFORE Wp is exactly

    H=K0[c]/(c^3), dim_(K0) H=3, dim_Q H=6.                      (13)

This is a finite local algebra with residue field K0 and maximal ideal
(c), of nilpotence index three. In rational basis

    (1,rho,c,rho*c,c^2,rho*c^2),

actual multiplication by rho uses three copies of

    M_rho0=[0 -6021/160; 1 -123/10],

while multiplication by c shifts the three two-dimensional blocks.
The matrices commute, reproduce the six cyclic unit basis vectors,
annihilate Q(M_rho) and have M_c^2!=0, M_c^3=0. The checker evaluates
Wp and the finite inverse (10) using these ACTUAL matrices. Their product
is the identity; Wp is not set to its residue scalar during the check.
This retains the nonzero nilpotent and its square rather than checking
only the reduced field.

The algebra (13) is a quotient of the actual residual trace algebra,
not a transplanted Image algebra. The Image Rodrigues certificate
supplied the strategy. No Artin parameter algebra is identified with
the original finite normalization or its codifferent.

## 6. Full degree-five graph closure and its exact scope

Suppose a v2 degree-at-most-five graph with nonzero constant Hessian
determinant existed over a characteristic-zero field. PF155 puts its
trace in (1)--(2). Its first jet is the forced S by the universal normal
recurrence, and its possible x^2 repair cannot change Wc or Wp.
Equations (4) contradict (12). This proves nonexistence over every
characteristic-zero field, including C. It closes the previously retained
complex classes; it does not impose a dual-kernel hypothesis on them.

There is also a precise consequence for the FULL nonreduced polynomial
candidate scheme. Let G5 be the finitely presented Q-algebra generated
by the 126 coefficients of a four-variable polynomial R of degree<=5,
C and C^(-1), modulo ALL coefficients of D_R-C. Then

    G5=0.                                                        (14)

If it were nonzero, a maximal ideal would have a residue field finite
algebraic over Q by Zariski's lemma, hence admitting an embedding in C.
It would give a complex graph prohibited above. This is the standard
maximal-ideal argument applied to the written nonexistence theorem;
a 128-variable Groebner calculation is not claimed. In particular no
unaccounted nilpotent thickening rescues the full bounded-degree graph
scheme. The explicit small certificate is (12) on the residual quotient.

PF154's A_2 ALONE is not asserted zero: it does not contain the full
polynomial graph's necessary highest-Hessian kernel condition used in
PF155. Neither a power of every leading deviation in that larger
partial-jet ring nor an all-normal increasing-power containment is
claimed. The uniform J^3 certificate here belongs to (7), not to an
unspecified larger trace algebra.

The degree-five v2 graph branch should now be retired rather than
regenerated as an open search. Unmarked higher-degree graphs and other
normal forms retain their existing open scopes. The marked standard
graph collision route was already closed in all degrees separately.
Global [PF151 normal-form coverage](NORMAL_FORM_COVERAGE_AND_COMMON_LEADING_DIVISORS.md)
and generic first-derivative integrality on the FULL finite normalization
remain open and take priority over another coordinate-architecture
expansion. All ramified boundary primes, the 3D integrality controls and
PF142's mixed minimal retirement are preserved. Full HC4/JC2 and an
actual counterexample remain OPEN. The narrow checker passes 43 exact
controls, including actual six-dimensional multiplication and the explicit
unit identity. Independent whole-proof replay, formal verification and
external review remain unclaimed.
