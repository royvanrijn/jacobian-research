# Gao F7: all 118 normalized fivefold-root targets have seven-point fibres

14 November 2026. Registered claim **GAOF7FIVE1**.

Work over C with precisely Gao's F7 on the nonzero first-target-coordinate
chart c!=0, using normalized Y=(cA,c^2 B,c^3 D,c^4 Z).
Exactly118 distinct normalized targets Y have an eliminant root of
multiplicity at least five. Each has one root of multiplicity exactly
five and seven other distinct simple roots. Its F7 fibre has exactly
seven source points.

The [exact parameterization](parameterization.json) is the origin plus
117 values indexed by the distinct roots of an explicit degree 117
integer polynomial. The a,b coordinates and Y are rational functions of
the marked root. All denominator-zero cases are covered below. In the
original c!=0 target chart this is a disjoint union of 118 curves, each
isomorphic to C*, obtained by varying c in the displayed normalization.
This counts geometric points and curves; no assertion about a potentially
nonreduced scheme structure is made.

In particular, three-point nonsimple lists(5,4) and(5,2,2) are impossible.
Together with [GAOF7CUBE1](../gao-f7-three-point-certificate/PROOF.md),
[GAOF7MULT1](../gao-f7-linear-cube/HIGH_ROOT_PROOF.md) and the complete
c=0 boundary count, only(4,3,2) and(3,2,2,2) remain for a possible
three-point fibre. Their occurrence or exclusion, full F7 strata, F5
and JC2 remain **OPEN**.

The [separate checker](verify.py) passed [5490 checks](receipt.json) in
10.15 seconds. The [integer certificate](certificate.json) is 1345979 bytes;
the [finite witness](finite-witness.json) is a small additional file.
No Groebner basis, symbolic resultant algorithm, factorizer or rational
reconstruction is run by the checker. Source transport, all exceptional
identities, determinants and forward existence identities are replayed.
The written arguments have no independent end-to-end or formal review,
and novelty priority is unestablished.

## Full source coverage and the residual septic

Use [Gao 2608.00222v2](https://arxiv.org/html/2608.00222v2), with

    b(t)=Y4-3Y2*t^2+2Y1*t^3-t^4,
    a(t)=Y3-2Y2*t+Y1*t^2-t^3-t^2*b(t),
    R(t;Y)=a(t)^2+a(t)+t*b(t)+t^3+Y1*t-Y2.

Precisely the simple roots of R lift to source points on c!=0.
For a marked multiple root s put a=a(s), b=b(s), and

    G2=a^2+a+s*b+s^3,
    G3=2s*G2+a+s^3+s^2*b,
    G4=3s^2*G2+b+s^4,
    X1=-det(d(G2,G3,G4)/d(s,a,b)).

The definitions and root equation force Y=(xi,G2+s*xi,G3+s^2*xi,G4+s^3*xi),
with xi=Y1. Substitution into R'(s;Y), with Y fixed during differentiation,
gives xi-X1. Thus every multiple root has Y=X(s,a,b), where

    X=(X1,G2+s*X1,G3+s^2*X1,G4+s^3*X1).

The checker rederives the root a,b reconstruction, R(s;Y)=0 and this
first-derivative identity by exact integer polynomial arithmetic. This
is full coverage, with no variable divided out.

Let f2,f3,f4 be the three primitive derivative constraints stored in
[quintuple-input.json](quintuple-input.json). The checker reconstructs
R(s+x;X) and proves each coefficient of x^j is a nonzero rational
multiple of f_j for j=2,3,4; coefficients0 and1 vanish identically.
It also verifies the full monic quotient identity

    R(t;X)=(t-s)^5*Q7(t;s,a,b)
              +sum_(j=2..4) c_j*f_j*(t-s)^j,

with nonzero rational constants c_j. The degree-seven monic Q7 has 905
coefficient terms. Thus f2=f3=f4=0 is equivalent to a root of multiplicity
at least five, at any s,a,b. This equivalence supplies both directions
of the eventual classification.

## Reduce the generic incidence to a finite projection

Write f2=A*b^2+B*b+C, where A=4s^2*(3s^3-2)^2. On
D=s*(3s^3-2)!=0 define the integer linear polynomials L_i=P_i*b+Q_i by

    L1=(A*f3-[b^2]f3*f2)/(4s*(3s^3-2)),
    L2=(A*f4-[b^2]f4*f2)/4.

The checker replays these literal polynomial identities. The bivariate
necessary expressions are

    C12=P1*Q2-P2*Q1=P1*L2-P2*L1,
    E_i=A*Q_i^2-B*P_i*Q_i+C*P_i^2
       =P_i^2*f2+L_i*(A*Q_i-A*P_i*b-B*P_i), i=1,2.

After dividing only nonzero rational constants and powers of s and
3s^3-2, call the primitive expressions e0,e3,e4 respectively. They are
the retained [coefficient equations](../gao-f7-linear-cube/high-root-bivariate-input.json)
from the maximum-multiplicity proof. No P_i is divided at this stage.

The full integer Sylvester determinant

    U(s)=Res_a(e0,e3)
        =c*s^14*(3s^3-2)^9*p(s), c!=0 in Q,

has degree 341, with p of degree 300. The a-degrees are7,8 and both
s-degrees are27, so the determinant degree is at most405. At406 distinct
integer values the verifier computes the nominal15x15 determinant in
FLINT and independently by plain-Python integer Bareiss, matching the
stored integer polynomial every time. This proves the entire polynomial
identity in characteristic zero. Nominal a-degrees are retained when
leading coefficients vanish; no such point is discarded.

At a common finite a-root the Sylvester matrix M(s) annihilates the
nonzero vector(a^14,...,a,1). Let J_a,J_0 be the signed cofactors deleting
row 0 and columns13,14. At rank 14, every adjugate column is a scalar
multiple of that kernel vector; at rank at most13 all cofactors are zero.
In either case

    J_a(s)-a*J_0(s)=0.

The stored factorization is J_a=K*N and J_0=K*T, with integer polynomials
N,T of degree 243 and K of degree 64. Each original cofactor has degree at
most378. The checker verifies both signed14x14 determinants at 379
distinct integers, in both integer engines, proving both polynomial
identities. Dividing K is postponed until a unit certificate justifies it.

The certificate supplies exact integer identities

    H210^2=p*V120,
    H210=(3s^3-1)*(3s^3+1)*H204,
    H204=H117*H87.

Subscripts indicate degrees, not irreducibility. The checker multiplies
these integer polynomials exactly. The first identity ensures every root
of p is a root of H210, regardless of multiplicity. No squarefreeness
assumption is needed for this coverage. The three cubic exceptional loci
and s=0 are handled explicitly below.

## The finite-free good-reduction unit lemma

The following is the characteristic-zero justification of the finite
witness. It is essential; a general modular unit-ideal report would
not suffice.

Let ell=1000000007 and R0=Z_(ell), the rationals whose denominators are
not divisible by ell. For an integer polynomial H of degree d with
ell not dividing its leading coefficient L, the algebra

    A_R0=R0[s]/(H/L)

is FREE of rank d over R0, with basis1,s,...,s^(d-1). Its reduction is
exactly F_ell[s]/(H mod ell), with the same dimension. A polynomial H
need not be squarefree for this statement.

For z in A_R0, form the dxd multiplication matrix m_z in that basis.
If det(m_z) mod ell is nonzero, det(m_z) is a unit of R0. Its adjugate
therefore gives an inverse over R0, so z is a unit in A_R0. In particular
z cannot vanish under evaluation at ANY complex root of H. Previously
certified inverses and source polynomial evaluations remain in A_R0,
and reduction commutes with all these operations. Hence the same lemma
applies successively to rational recovery formulas.

Every H below is specified over Z before reduction. The checker proves
ell does not divide its actual leading coefficient, preserves degree,
and constructs the complete multiplication matrix. FLINT and a separate
plain-Python Gaussian elimination verify each determinant. Matrix entries
and source evaluations are also reconstructed using independent Python
coefficient-list polynomial arithmetic. The prime is proved by complete
trial division. Thus there is no unverified flatness, rank or leading-
coefficient assumption in passing from the finite witness to C.

## Recover a,b and prove residual squarefreeness

The exact nonzero determinant residues are:

| Algebra | Element whose multiplication is tested | Dimension | Residue |
|---|---|---:|---:|
| Q[s]/H204 | K |204|959718345|
| Q[s]/H204 | T |204|480911749|
| Q[s]/H87 | e4(s,N/T) |87|100160574|
| Q[s]/H117 | P1(s,N/T) |117|487226214|
| Q[s]/H117 | Res_t(Q7,Q7') |117|782152754|
| Q[s]/H117 | Q7(s) |117|76753123|
| Q[s]/H117 | s*(3s^3-2) |117|121992835|
| Q[s]/H117 | H117' |117|242907863|

In the last five rows substitute a=N/T and b=-Q1/P1. The finite-free
lemma makes all inverses legitimate over R0 and over C.

For an actual generic fivefold point outside3s^3=+/-1, projection gives
H204(s)=0. The first two rows make K,T nonzero there, so the adjugate
relation forces a=N/T. If H87(s)=0, the third row contradicts the
necessary e4=0. Thus H117(s)=0. The fourth row now makes P1 nonzero,
so L1=0 forces b=-Q1/P1. No zero of a b-linear coefficient is discarded;
its exclusion on this retained factor is a proved unit statement.

The coefficients of Q7 in this finite algebra are evaluated directly
from the characteristic-zero source formulas. The checker reproduces
them with FLINT polynomial arithmetic and plain Python coefficient
lists. Their canonical representatives have s-degree at most116.
The 13x13 Sylvester determinant with Q7' therefore has s-degree at most
13*116=1508. Its finite-field polynomial witness is checked at 1509
distinct field values, with both determinant engines; this proves the
entire identity over F_ell, not a sampled test. It correctly computes
the reduction of the abstract resultant in A_R0, since determinants
commute with specialization. The fifth norm then proves this resultant
a unit in characteristic zero. The sixth norm similarly proves
Q7(s)!=0. Consequently R has exactly a fivefold root at s and seven
other distinct simple roots on the generic piece.

The finite-field interpolation alone is not a characteristic-zero
coefficient reconstruction. The characteristic-zero conclusion comes
from the proved finite FREE presentation and its nonzero multiplication
determinants. This distinction prevents the component-loss problem of
unrestricted modular ideal calculations.

## Every retained root really occurs: exact count117+1

It remains to prove existence for every retained parameter, rather than
merely use it as an overestimate. Let nu=N, tau=T. For e0,e3 respectively,
the checker constructs the integer homogeneous substitutions

    tau^7*e0(s,nu/tau),
    tau^8*e3(s,nu/tau).

Their degrees are1728 and1968. They are exactly divisible by H117,
with integer quotient degrees1611 and1851. Both integer division and
full coefficient multiplication are checked, the latter also with
plain Python integer convolution. Thus e0=e3=0 in Q[s]/H117.

The certified units tau,P1,D justify the recovery formulas at every
root of H117. Set a=N/T and b=-Q1/P1. Then L1=0 by definition, and
E1=0 implies P1^2*f2=0, so f2=0. Also C12=0 gives P1*L2=0, so L2=0.
The two cancellation identities now give A*f3=A*f4=0. Since A=4D^2
is a unit, f3=f4=0. The full residual identity therefore constructs a
genuine fivefold target X(s,a,b) for EVERY root of H117.

The last norm proves H117' a unit modulo H117, so its117 complex roots
are distinct. Its constant coefficient is nonzero, so none is zero.
Each root has a unique recovered a,b. Two distinct roots cannot give
the same target Y: that eliminant would have two fivefold roots,
contradicting the already proved seven-simple-root residual property.
Thus these are117 distinct normalized targets.

## Every exceptional root value is accounted for

The four source specializations in the certificate are reconstructed
from f2,f3,f4. Only nonzero rational constants are cleared.

- At s=0, two exact rational identities express a and b in the ideal of
  the three specialized equations, with 30 multiplier terms each.
  Hence a=b=0 is the only possible point. It is an actual point:
  R=t^5*(t^7-2t^4+2t-1), and the septic has discriminant713 and value -1
  at zero, so its seven roots are simple and disjoint from the marked root.
- At3s^3=2, adjoin that equation to the three polynomial remainders.
  A399-term rational unit identity excludes all three complex s-values.
- At3s^3=1, a 301-term unit identity excludes all three values.
- At3s^3=-1, a 452-term unit identity excludes all three values.

All five identities are replayed after clearing a positive common
constant denominator, in FLINT and a separate sparse Python integer
convolution. The origin supplies one additional target, different from
the117 generic ones. The earlier projection identities and exceptional
split prove that there are no others. This proves the exact total118.
No numerical complex roots or finite-box extrapolation is used.

## Fibre counts and remaining open cases

Every target just classified has exactly seven simple eliminant roots,
so Gao's reconstruction gives exactly seven F7 source points. Varying
c!=0 with Y fixed gives original targets

    (c,Y1/c,Y2/c^2,Y3/c^3,Y4/c^4).

Distinct normalized Y-values give disjoint curves because their first
coordinate fixes c. There are exactly 118 such C* curves on the chart.
No conclusion about their closures at c=0 is needed.

The fivefold types(5,4) and(5,2,2) cannot occur in a three-point fibre.
The earlier cube obstruction and maximum-five bound removed the other
four types from the original eight. Only(4,3,2) and(3,2,2,2) remain.
The complete c=0 boundary result gives only one- and two-point fibres.
This does not settle global cardinality three, all multiplicity strata,
F5 or JC2.

## Reproduction, literature boundary and assurance

With [the pinned dependencies](requirements.txt), from the repository root:

```sh
timeout 180s .venv/bin/python research/jacobian/gao-f7-quintuple/verify.py
```

The receipt binds checker, integer certificate, finite witness and source
coefficient input by SHA256. Python3.12.3, SymPy1.14.0 and python-flint0.9.0
were used. Source derivation takes about 3 seconds; direct univariate
resultant and cofactor construction take fractions of a second. The
untrusted exceptional exporters each had20 second,1536 MiB address-space
and10 MiB output caps. Projection probes used120 second and2 GiB caps.
The final checker completed in 10.15 seconds within its180 second time cap.

An additional exact rational-univariate probe independently checks
f2=f3=f4=0 in Q[s]/H117, but recovering b by a dense rational inverse took
about 82 seconds. The homogeneous integer identities above avoid that
coefficient growth and are the canonical faster replay. Both attempts
are retained. No failed calculation was repeated unchanged or promoted
to a theorem from modular output. The earlier exclusion-only checker
and receipt are retained in[history](history/exclusion-receipt.json).

Gao v2 already proves the seven-point axis fibre and the existence of
an isolated algebraic quintuple-root point; those facts are prior work.
The complete locus and118 count here are supported by this certificate,
without a priority claim. arXiv version checks on 14 November still list
[Gao v2](https://arxiv.org/abs/2608.00222v2) and
[de Bondt v1](https://arxiv.org/abs/2610.24889v1). No new external theorem
is imported into the elimination or finite-free unit argument.

There is no independent end-to-end review, formal verification, external
review or author endorsement. MATH_STATUS.json is the status authority.
