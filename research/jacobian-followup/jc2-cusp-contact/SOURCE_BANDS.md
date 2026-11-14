# Selected actual source equations on the contact507 stratum

14 November 2026. Claim JC2CFSB1. The purpose is a small exact
coefficient obstruction on a surviving F2 stratum, rather than another
degree bound or a reconstruction of the historic nonlinear circuit.

## 1. Scope, coordinates, and the retained root-continuation input

Work on the ACTUAL lowest stratum L=39, cusp contact507 of JC2CFN1,
with q=y, v=xy^5, g=-(5/9)Q and J={P,g}!=0.
Normalize the Wronskian coefficients by J: dividing W and all its
coefficient functions by J reduces the identities below to J=1 without
changing any support, rank, or vanishing condition.

Write

    P=q^-15*C(q,v)^3,
    g=q^-25*C(q,v)^5 + terms beyond descent35,
    W/J=q^-39*E(q,v),
    C=c+q*k1+q^2*k2+q^3*k3+...,
    E=H+q*K1+q^2*K2+q^3*K3+...,
    H=c^7*d.

Here C is the formal cube root of q^15*P, with chosen leading coefficient
c=v(v-1)^2 R. It is initially only a series over C(v).
The absence of all seven integer centres on L39 means
g=P^(5/3) through relative q-order35. Thus the displayed substitutions
are valid through every coefficient used below, even after differentiating.

We use the retained PF2L1 primitive-residue theorem: the first seven
common-root coefficients k1,...,k7 are POLYNOMIAL in v.
This is a nonlinear conclusion, not the false claim that the first four
tangent rows force divisibility. In brief, at a simple prime of c the
second/third binomial coefficients force root continuation. At a double
prime, the potentially surviving valuation-three defect gives residues

    E4=a^4-9a^2*b+27b^2,
    E5=7a^4-60a^2*b+135b^2,
    Res_b(E4,E5)=1701a^8.

If the first defect is at spacing n<=7, both residues occur before
the relative target descent36 since 5n<36; they cannot both vanish for
a nonzero leading defect. Induction after subtracting the previous
polynomial root coefficients gives the stated continuation. On L39 there
are no intervening target centralizers through descent35. The new checker
reproduces the two primitive coefficients and the exact resultant.
The complete valuation/spacing argument remains the retained PF2L1 input.

## 2. A two-variable equation retaining the actual Keller forcing

Since {v,q}=q^5, the exact identity {g,W}=-5J*P^4 gives, to the
orders used here,

    (5C-q*C_q)*E_v - C_v*(39E-q*E_q) = -C^8.            (2.1)

The checker verifies the full rational-coordinate identity before
truncation. Its constant coefficient is

    5c*H'-39c'*H=-c^8,
    5c*d'-4c'*d=-c.                                    (2.2)

If earlier kn and Kn vanish, the relative order-n equation is the
homogeneous LINEAR equation

    5c*Kn'-(39-n)c'*Kn
       +(5-n)kn*H'-39kn'*H+8c^7*kn=0.                  (2.3)

Equivalently, putting Kn=c^6*e when such divisibility is justified,

    5c*e'+(n-9)c'*e
      =-(35+n)c*kn/5 +39(n-5)c'*kn*d/5+39c*kn'*d.

The checker uses the full polynomial Kn in (2.3); it does not assume
this optional stronger factorization.

For the squarefree case, vanishing at n=1 then n=2 makes (2.3)
the ACTUAL n=2 and n=3 equation, not just its tangent approximation.

## 3. Complete admissible coefficient spaces

The retained source envelope and translated terminal edge imply, for
A=P_(15-n),

    ord_0 A>=ceil((15-n)/5),
    deg_v A<=floor((90-n)/6),
    ord_1 A>=ceil((72-5n)/12).

When the earlier root coefficients vanish, A=3c^2*kn.
For n=1,2 this gives

    kn in v(v-1)^2 * C[v]_(<=1),

and at n=3 it gives

    k3 in v(v-1) * C[v]_(<=2).

The support of the ACTUAL W coefficient Kn at band m=39-n satisfies

    ord_0 Kn>=ceil(m/5),
    deg_v Kn<=floor((507+2m)/15),
    ord_1 Kn>=ceil((5m-15)/12).                         (3.1)

The first is source polynomiality, the second follows from the established
positive weight bound deg_(5,2) W=507, and the third is the retained
terminal height at most15. These are necessary support conditions,
sufficient to specify the full spaces used in this selected subsystem.
They are not asserted to be the complete remaining Keller equations.

For n=1,2 there are16 Kn coefficients and2 kn coefficients. At n=3
there are17 Kn coefficients and3 kn coefficients. Equating all
coefficients of the polynomial (2.3) gives small exact linear systems.

## 4. Squarefree carrier: two forced zero bands, then one parameter

Take

    R=(v^2-3v+3)/25,
    c=v(v-1)^2*R,
    d=v(v-1)*(v^2-3v+3)/3.

The exact systems have the following kernel dimensions:

| relative order n | unknowns | rank | kernel dimension |
| ---: | ---: | ---: | ---: |
| 1 | 18 | 18 | 0 |
| 2 | 18 | 18 | 0 |
| 3 | 20 | 19 | 1 |

Consequently every actual squarefree-carrier F2 pair on L39 has

    k1=k2=K1=K2=0,
    P_14=P_13=g_24=g_23=0.                             (4.1)

This excludes the entire sublocus on this stratum where either P_14
or P_13 is nonzero. It is a statement about actual source coefficients,
beyond the highest-face calculation.

The complete remaining third-order possibility in this subsystem is

    k3=tau*v(v-1)*(5v^2-15v+12)/5,
    K3=5tau*c^7*v*(13v^2-39v+31).                      (4.2)

Thus P_12=3c^2*k3 and g_22=5c^4*k3. For unnormalized W, its band36
coefficient is J*K3. The coefficient tau need not be nonzero and is
not proved to extend to a full Keller pair. The checker verifies both
the full ranks and that (4.2) spans the exact kernel.

This is the substantive actual-coefficient obstruction sought in the
scoped goal. It does not exclude the whole squarefree L39 stratum:
tau and deeper source coefficients still require the other equations.

## 4A. The squarefree third parameter also survives its quadratic return

The order6 equation is nonlinear even though orders1/2 vanish.
Its forcing from (4.2), removing tau from k3,K3, is

    tau^2 * (2*k3*K3' -36*k3'*K3 +28c^6*k3^2).

The complete coefficient matrix has43 rows and24 columns, rank22,
and augmented rank22. An especially small particular solution is

    k6_part=-3*(2v-3),
    K6_part=c^6*v(v-1)*(221v^5-1326v^4+3027v^3
                         -3090v^2+1137v+72).            (4A.1)

Multiplying both by tau^2 solves this forcing. Two homogeneous
coordinates remain. The known cubic-root contribution3c*k3^2 to P_9
and10c^3*k3^2 to g_19 already satisfies the respective source/terminal
bounds, so the full new-k6 space used here is legitimate.

In particular setting orders4,5,7 to zero and taking (4.2),(4A.1)
gives a nonzero-parameter solution of this SELECTED PDE subsystem through
relative order7. It is not a polynomial Keller solution or a point of
the complete historic circuit. It is an exact countercontrol: testing
the quadratic return alone cannot exclude tau.

A literal generic SymPy Matrix.rank call became too slow on this tiny
rational system and was terminated. Exact DomainMatrix arithmetic over QQ
computed ranks and the displayed witness promptly; this implementation
change is retained in the checker. No interrupted calculation is used
as evidence.

## 5. Double carrier: first parameter survives the nonlinear next row

Take rho^2-3rho+1=0, alpha=3(rho+1)/5,

    R=(v-rho)^2/[25(1-rho)^2],
    d=v(v-1)(v-rho)(v-alpha)/(rho*alpha).

The first-order system instead has nullity one. Normalize its parameter
as t and define

    beta=(54rho-261)/142,
    A=(-6438+2397rho)/142,
    B=(9265-3332rho)/142,
    D=(15rho-45)/2.

Then its complete kernel is

    k1=t*v(v-1)^2*(beta*v+1),
    K1=t*c^6*v^2*(v-1)^3*(v-rho)*(A*v^2+B*v+D).        (5.1)

The checker uses one embedding rho=(3+sqrt(5))/2 and exact arithmetic
over Q(sqrt(5)). Galois conjugation proves the identical result for the
other retained carrier. No floating-point tolerance is used.

At order2 the kn/Kn equations are no longer homogeneous. The nonlinear
forcing, with the t removed from the displayed first-order generators,
is

    t^2 * (4*k1*K1' -38*k1'*K1 +28c^6*k1^2).           (5.2)

The coefficient matrix has18 columns, rank17, and augmented rank17.
The checker obtains an exact particular solution by row reduction,
retains its coefficients in the receipt, and substitutes it back.
Thus this quadratic forcing lies in the image: it does NOT force
t^2=0 or t=0. One new homogeneous coordinate remains.

The use of the same k2 support as in Section3 is justified here despite
nonzero k1: the known contribution3c*k1^2 to P_13 already satisfies
every source-band degree and terminal-jet bound. Likewise the known
contribution10c^3*k1^2 to g_23 satisfies its corresponding bounds.
Subtracting them therefore leaves the same allowed k2 space.

This is an exact solvability countercontrol for the next nonlinear row.
It is neither a full-circuit point nor evidence that the double branch
exists globally. It records where this particular obstruction stops.

## 6. Remaining gate and reproducibility

The squarefree branch now starts with two vanished root coefficients
and one explicit next parameter whose quadratic return is solvable.
A next projection must retain the order4/5 freedoms and all couplings.
The cubic return at order9 meets the previously unresolved later
first-defect range: polynomial root continuation beyond order7 must
not be assumed without the required source divisibility equations.
The double branch retains its first parameter and a compatible
quadratic correction. A later cokernel projection, the missing global
attachment theorem, or a different obstruction is needed.

[verify_source_bands.py](verify_source_bands.py) checks the written
identities, exact source spaces, ranks, complete kernels and nonlinear
particular solution. Its receipt records the retained source hashes.
No historical large artifact is rebuilt. The universal deductions are
written proofs supported by exact finite linear algebra, without
formal proof-assistant verification, independent replay, external
review, or novelty-priority claims.

