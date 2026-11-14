# Universal exclusion of the gap-21 quartic conic-contact gate

14 November 2026. PF2R7QALL1 extends the two-class contact obstruction to
the **entire** birational polynomial degree-(6,10), gap-21 family.
PF2R7QEX1 applies it to the remaining quartic infinity escape of the
degree-six Keller, entire single ramified affine-prime (e,f)=(2,1)
programme. The immersed quartic family is excluded as well. The special
reciprocal cubic family and the entire-ramification hypothesis remain
open gates. No source finality or Theorem3.8 conclusion is used.

This is a written algebraic/geometric proof with literal multivariate
unit certificates. No formalization, external review or novelty priority
is claimed. The raw Gröbner-basis discoveries and capped unsuccessful
cofactor constructions are not used as proof.

## 1. Universal input and normalization

Use the retained normalization

    p=t^6+a*t^4+b*t^3+c*t^2+d*t,
    q=[p^(5/3)+B*p^(4/3)+D*p^(2/3)+E*p^(1/3)]_+,
    rho1=rho3=rho5=rho7=rho9=0.

The gap-21 family also has rho11!=0, but that inequality is not needed
below. The literal parity identity in PF2R7PAR2 shows b=0 forces d=0,
making both p and q even, contrary to birationality. Thus b!=0. Over C,
scale t and the two target coordinates to normalize b=1. These are
affine changes; they preserve conics, tangency and quartic polynomial
parameter maps. Free terms C*p+F in q are removed by an affine shear and
translation. Hence no coefficient stratum of the stated family is lost.

Suppose there is a quartic polynomial parameter

    Z(t)=t^4+u*t^3+v*t^2+w*t

whose four generic conjugate points on the branch curve lie on a common
conic tangent there. The [doubled-quartic contact matrix](../jc2-r7-quartic-contact/PROOF.md)
is necessary. Denote its four-row minors by Delta_I, with zero-based
row indices in the six-row matrix, and put

    F_j=[Z^j] Delta_(0,1,2,3),
    L=-21*B-20*a+30*v,
    K=-21*B+20*a-60*d.

Every selected minor coefficient below must vanish. All identities are
over Q in the displayed parameters before imposing the gap equations,
unless explicitly stated otherwise.

## 2. The cubic and quadratic quartic coefficients are forced

The universal coefficient identity is

    [Z^10] Delta_(0,2,3,4)=-80*u^5.

Therefore u=0. At u=0 and b=1,

    F10=-(4/9)*(3*w-2)^2*L.                          (A)

We first prove L=0. Otherwise (A) forces w=2/3. At w=2/3 the next
coefficient is

    F9=(4/81)*L*(-4*a+3*v+6*d)^2.

Thus v=4*a/3-2*d and L=K!=0. With this substitution,

    F8=-(16/81)*K*(2*a*d-c-3*d^2)^2.

So c=2*a*d-3*d^2. At these same substitutions the literal gap/contact
identity is

    F7-(8/27)*rho1 = 4*K/6561.

Both F7 and rho1 vanish, a contradiction. Consequently

    L=0,  v=2*a/3+7*B/10.                            (B)

## 3. An unused minor forces the linear quartic coefficient

At (B), define

    R=-3321*B^3-6300*B*a^2+18900*B*c+32400*D
      +4000*a^3-18000*a*c-6000-22500*w+27000*w^2.

Exact coefficient extraction gives

    F9=(3*w-2)^2*R/6075.                             (C)

If w!=2/3, then R=0 and D is eliminated by its nonzero constant
coefficient32400:

    D=(3321*B^3+6300*B*a^2-18900*B*c-4000*a^3
       +18000*a*c+6000+22500*w-27000*w^2)/32400.       (D)

In the b=1 chart rho1 is affine-linear in E with coefficient1/3, so
rho1=0 uniquely eliminates E as a polynomial in the other parameters.
After these substitutions, a different necessary minor has the exact
coefficient

    [Z^8] Delta_(0,1,3,5)=-(40/81)*(3*w-2)^4.         (E)

No further gap equation or ideal reduction is required for (E). It
contradicts w!=2/3. Therefore every possible quartic contact must satisfy

    u=0, v=2*a/3+7*B/10, w=2/3.                      (F)

The additional minor is essential to the short proof. Attempting to
expand a single global cofactor from only the principal minor led to
large expressions and capped computations; none of those outputs is
promoted to a certificate.

## 4. Two four-variable unit certificates finish the exclusion

Put w=2/3 and v as in (F), and let R0 be R at w=2/3. A further exact
identity is

    F8+(4*K/15)*rho1=-K^2*R0/5467500.                (G)

Thus either K=0 or R0=0. In each branch eliminate E using rho1, whose
coefficient1/3 is a unit; use only the equations listed below.

### 4.1 K=0

Set B=(20*a-60*d)/21 and v=4*a/3-2*d. The six equations

    rho3,rho5,rho7,rho9,F7,F6

are polynomials in D,c,d,a after eliminating E. Clear rational
denominators and divide each by its rational coefficient content.
The first certificate in [certificates.json](certificates.json) gives

    sum_(i=1)^6 U_i(D,c,d,a)*g_i(D,c,d,a)=1.          (H)

The checker reconstructs every generator from the original gap recurrence
and universal minors, then replays (H) by separate sparse Fraction
polynomial multiplication. Hence this branch has no complex solution.

### 4.2 R0=0

Set D by (D) with w=2/3, and eliminate E by rho1. The five equations

    rho3,rho5,rho7,rho9,F7

are polynomials in B,c,d,a. The second certificate gives

    sum_(i=1)^5 V_i(B,c,d,a)*h_i(B,c,d,a)=1.          (I)

It is reconstructed and replayed in the same way. This branch too has
no complex solution. Equations (G)--(I) include the overlap K=R0=0;
no division by either K or R0 has silently discarded it.

All alternatives from (A) through (I) are exhausted. There is therefore
no conic with four generic contacts over any quartic polynomial parameter
map for any pair in the stated birational gap-21 family. This proves
PF2R7QALL1, including every immersed specialization and multiple-root
derivative stratum. No assertion that p' is squarefree is used.

## 5. Exact coefficient extraction and its boundary

[weighted_contact.py](weighted_contact.py) extracts the necessary universal
minor coefficients without expanding entire large determinants. This
uses the exact grading

    wt(t)=1, wt(Z)=4,
    wt(a,b,c,d,B,D,E,u,v,w)=(2,3,4,5,2,6,8,1,2,3).

The column weights for q,p^2,p*q,q^2 are10,12,16,20. In row corresponding
to t^r, an entry's coefficient of Z^k has weight column_weight-r-4*k.
For a selected four-row determinant coefficient of Z^j, the total
parameter weight is therefore the sum of column weights minus the sum
of row powers minus4*j. All parameter weights are nonnegative. An entry
coefficient with larger weight cannot contribute and is omitted exactly.

This filtering is performed before substitutions such as b=1 or (D),
which need not preserve the grading. Coefficient extraction commutes
with those parameter substitutions because none involves Z. The result
is an exact polynomial identity, not an asymptotic truncation, modular
calculation, or numerical sample. The monic doubled-quartic remainders
and pivot formula are the same as in the earlier directly checked
contact matrix. Additional exact specializations compare the weighted
implementation to that full remainder-and-determinant construction.

The two multivariate certificates are frozen polynomial arrays. Their
units are checked by a separate Fraction implementation; it neither
computes a Gröbner basis nor trusts Singular's printed basis. Mutating
one multiplier makes each unit check fail. The discovery logs retain
the unsuccessful computational routes separately.

## 6. Geometric consequence and the remaining frontier

PF2R7QBASE1 proves that every remaining quartic infinity escape under
PF2R7INF1 hypotheses must have exactly the polynomial base and four
contacts excluded above. Therefore that entire quartic escape is
impossible. This is PF2R7QEX1. Its inherited geometric inputs are written
proofs; the present coefficient certificates do not formalize them.

Combined with the earlier balanced-infinity and varying-direction cubic
exclusions, only the reciprocal case can remain:

    infinity indices(1,5), fibre degrees(1,0),
    deg(E -> pencil base)=3,
    a finite u with ord_u(p')-ord_u(q')=2.

For an immersed image that finite point is an exact double root of p'.
PF2R7G1's finite/infinite coefficient profiles still apply. Both exact
cusp classes have squarefree p', so they remain excluded independently
of finality; now every squarefree-derivative immersed pair is excluded
as well. The multiple-root immersed family is not thereby closed.

The entire unique ramified-prime hypothesis is still explicit and must
be established in intended terminal-passport applications. Other
ramification patterns, normal r5, nonnormal slices, (75,125), and JC2
remain unresolved. No whole normal-r7 or JC2 closure is claimed.
