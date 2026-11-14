# Generic single-branch polar obstruction and the two-root frontier

14 November 2026. PF2RULED21. This closes the two higher-polar profiles
left as necessary, unresolved data in [PF2RULED20](QUARTIC_GENERIC_POLAR_GAP.md).
The full finite-map and actual compatible-normalization goal stays active.

## Statement and hypotheses

Work in characteristic zero, with fixed local boundary coordinates z,w,
B=(z=0), and a pencil F_c=A(z,w)-c*z^m, where A is independent of c.
Generic means over the algebraic closure of C(c). The whole reduced local
fibre germ, rather than just one branch in a multibranch germ, is assumed
unibranch when specified below. Put l=I(F_c,z) and I=I(F_c,A_w).

1. If the generic whole germ is unibranch and ord_t(w)=2, then, for m>=3,
   I <= (m/2+1)*l-2. In the finite-map setting with volume divisor E-3B,
   the local companion pole is consequently at most
   (3-m/2)*l-2-I(F_c,E). Thus this type is impossible for every m>=6.
2. A generic unibranch pencil A-c*z^6 with A(0,w)=w^8*unit cannot have
   polar intersection33. Its conjugate-difference data would force
   ord_t(w)=2, whereas the preceding bound gives I<=30.
3. In an actual finite quartic pole(6,9) pair, every away common-nonic
   root of multiplicity4 costs at least two units of Q pole degree.
   This is only a lower bound; sharpness or two local branches is not proved.
4. Noncontact432 and contact3/away42 each cost at least five Q poles,
   exceeding degree4. The only remaining noncontact types are72 and54,
   with the sixteen marked leading candidates of PF2RULED19. Every
   currently retained common nonic has support at most two points. Its
   auxiliary cubic-cover normalization has rational components, and its
   leading-period test already passes. These are not finite-map existence
   or target-intermediate-field assertions.

## The critical-value proof in the original coordinates

Suppose ord_t(w)=2. The intersection formula gives ord_z F_c(z,0)=2.
Since m>=3, A_zz(0,0) is nonzero and A_z(0,0)=0. The formal implicit
function theorem gives a unique critical curve z=phi_c(w), phi_c(0)=0,
solving (F_c)_z=0. Define the critical value K_c(w)=F_c(phi_c(w),w).
All these coefficients are algebraic or rational functions of c; indeed
recursive division is only by the fixed nonzero A_zz(0,0).

Taylor expansion and a formal square root give

    F_c(z,w)=K_c(w)+(z-phi_c(w))^2 U_c(z,w)=K_c(w)+Z^2,
    U_c(0,0)!=0, Z=(z-phi_c(w))*sqrt(U_c(z,w)).

This coordinate is used to compute valuations only. Parameter dependence
is differentiated in the original coordinates, before using this square
completion. Write N=ord_w K_c. K_c is nonzero: otherwise F_c is not reduced.
The reduced unibranch assumption makes N odd, and on the normalization
ord_t Z=N. Therefore ord_t (F_c)_z=N, since 2Z*Z_z differs from Z by a unit.
Differentiating along the branch gives

    (F_c)_w w'(t)+(F_c)_z z'(t)=0,
    ord_t w'=1, ord_t z'=l-1,
    I=N+l-2.

If N<=l then N<=m*l/2 immediately. If N>l, the equality
ord_t(z-phi_c(w))=N forces the leading z term to come from phi_c(w).
Thus l is even and ord_w phi_c=l/2. At fixed w the exact chain rule is

    partial_c K_c = (F_c)_z(phi_c(w),w)*partial_c phi_c - phi_c(w)^m
                  = -phi_c(w)^m.

Differentiating coefficients cannot create a lower w power. Hence
N<=ord_w(partial_c K_c)=m*l/2. In either case,

    I=N+l-2 <= (m/2+1)*l-2.

The [generic-polar identity](POLAR_CUSP_LEADING_SECTORS.md) now gives

    nu = I-(m-2)*l-I(F_c,E)
       <= (3-m/2)*l-2-I(F_c,E).

For m>=6 the right side is negative, contradicting the strictly positive
pole required by finiteness. This argument is all-series and all-height;
it is not a bounded jet computation. It does not exclude multibranch germs
whose individual branches have w order2. In particular it does not remove
the two-cusp quartic model in PF2RULED20.

## Exhaustive degree-eight conjugate differences

For a reduced unibranch degree-eight germ normalize the parameter so
z=t^8. Its other coordinate w(t) generates the degree-eight extension.
The roots of its monic Weierstrass polynomial are w(epsilon*t), epsilon^8=1;
a Weierstrass unit does not affect the derivative valuation on F_c=0.
Therefore

    I=sum_{epsilon!=1} ord_t(w(t)-w(epsilon*t)).

Let j_1<...<j_s be precisely the exponents that lower the running gcd,
e_0=8, e_i=gcd(e_(i-1),j_i), e_s=1. Counting the roots of unity first
separated at j_i gives the exact formula

    I=sum_i (e_(i-1)-e_i)*j_i.

For I=33 each j_i<=33, so an exhaustive finite search is justified by
this formula itself. The only chain is (j_1,e_1)=(2,2), (j_2,e_2)=(21,1).
For a direct check: e_1=1 would give 7*j_1=33; e_1=2 forces j_1=2 and
j_2=21; e_1=4 forces j_1=4, and its remaining contribution17 cannot be
3*j_2, or 2*j_2+j_3 with j_2>=6,j_3>=7. Larger j_1 already exceed33.
Thus six conjugates separate at order2 and one at order21, giving12+21.
No exponent divisible by8 can precede order2, so the original w itself
has order2. No parameter-dependent ambient coordinate change is needed.

For m=6,l=8,I=33 this would require N=27 and ord_w phi_c=4. But the
critical-value identity gives ord_w partial_c K_c=24, contradicting N=27.

## A selected fibre really can have polar33

The family

    A=(z-w^4)^2-w^27, F_c=A-c*z^6

has at c=0 the reduced single branch w=t^2,z=t^8+t^27. Direct substitution
gives F_w=-8*t^33-27*t^52 and F_z=2*t^27, so its polar intersection is33.
For generic c substitute w=t,z=t^4+t^12*X:

    F_c/t^24=X^2-t^3-c*(1+t^8*X)^6.

The two simple initial roots of X^2-c lift to two smooth branches. Each
has polar order15, giving generic total30 and Q pole15-4*4=-1. This
pencil is not a finite-map example. Its critical curve begins
phi_c=w^4+3*c*w^20+O(c^2), and K_c=-w^27-c*w^24+O(c^2).
This both replays the obstruction and rules out an all-selected-fibres
interpretation. No generic multibranch polar33 exclusion is asserted.

## Charging poles and the current root frontier

At an away common-nonic multiplicity4 root, P has pole6 and leading
multiplicity8. If its total local Q pole degree were1, finiteness would
force exactly one normalized branch and the polar identity would force
I=4*8+1=33. The preceding theorem excludes this. The local cost is
therefore at least2. From PF2RULED20 an away common multiplicity2 root
has exactly two cusp branches and cost2. Other distinct roots cost at
least1, including q when present.

| Preceding root type | New minimum Q pole degree | Result |
| --- | ---: | --- |
| noncontact432 | 2+1+2=5 | excluded |
| contact3, away42 | 1+2+2=5 | excluded |
| noncontact72 | 1+2=3 | remains necessary |
| noncontact54 | 1+2=3 | remains necessary |

The sixteen exact marked points in the last two sectors retain their
fixed-coordinate meaning; no new quotient by source automorphisms is taken.
The three numerical groups remain96/144 at contact0;80/128 at contacts
3,4,5,6,7,9;48/128 at contacts7,9. For contact, every remaining common
nonic has divisor (9-k)*p+k*q, k in {3,4,5,6,7}, or9*q. The contact3
away6 case is a cube. These are necessary forms, with no claim that
all positions or full sections produce a finite pair.

## What the auxiliary period gate now does and does not do

For noncontact72 or54, PF2RULED18 already proves exactness for every two
distinct roots. For a contact noncube there is one finite root when q
is infinity, so its one-root criterion applies. The cube cases6p+3q,
3p+6q,9q pass the retained rational criterion. Thus every presently
retained pole69 root form passes the leading-period test. This conclusion
concerns this filtered frontier, not every degree-nine period locus.

For clarity one can explicitly parametrize the auxiliary cover. For a
two-support noncube, an abstract Mobius coordinate x sends its two support
points to0,infinity. Since the common section has degree9, its defining
rational function differs by a cube from x^k, k in {2,4,5,7}. Write

    theta^3=x^k, 3u+kv=1, T=x^u*theta^v.

Then T^3=x and theta=T^k, so the function field is C(T). If the section
is a cube, its cover splits into three rational components. This coordinate
is used only on the auxiliary curve. It is neither a source automorphism
preserving the marked model nor an intermediate field of an actual quartic
target extension. The global target cusp equations, source Jacobian,
finiteness and actual compatible completion remain unproved.

The next useful gate is full cusp/Jacobian/finiteness on the sixteen marked
two-root leading candidates and the remaining contact cases. Three/four-root
period-position campaigns cannot advance the current pole69 frontier.

## Reproduction and limits

The [checker](verify_unibranch_polar33.py), [certificates](unibranch_polar33_certificates.json)
and [receipt](unibranch_polar33_receipt.json) replay the exhaustive conjugate
gcd chains, critical-flow identities, exact special/generic countercontrol,
pole costs, predecessor seals and rational auxiliary parametrizations.
The all-series coverage is supplied by the written proof above, not by
SymPy or a finite expansion. No formal verification, independent full
reproof, external review, novelty claim, global P,Q construction, full
pole69 exclusion, height-nine exclusion or completion of the original
goal is asserted.
