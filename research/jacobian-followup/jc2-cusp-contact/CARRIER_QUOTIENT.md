# Affine carrier gauge and the finite degree-seven quintic parameter map

14 November 2026. Claims JC2CFGA1 and JC2CF7F1. This is a scoped continuation
of [the leading contact spectrum](PROOF.md) and
[the contact507 branch](BRANCH507.md). It does not solve the full F2
coefficient equations or the global carrier-to-nonproper-curve attachment.

## 1. Exact target action on seven carrier centres

Use g=-(5/9)Q, pi=P^3/g^2 and h=P^5/g^3. At the generic retained F2
carrier, pi has source order5. There are seven integer-power centres

    h=1+h1*pi+h2*pi^2+...+h7*pi^7+remainder,

where the remainder has order strictly greater than7 in pi. It need not
have order8: the first nonintegral correction can have order36/5.
The retained carrier descent supplies these seven centres. All calculations
below are in the quotient by terms of order greater than7; they therefore
apply with either remainder convention.

The determinant-one affine target transformations

    T(mu,eta,nu): (P,g) -> (P+mu,g+eta*P+nu)

act by the exact rational identities

    A=1+mu*pi^3/h^2,
    B=1+eta*pi^2/h+nu*pi^5/h^3,
    pi'=pi*A^3/B^2,
    h'=h*A^5/B^3.                                      (1.1)

This action is tangent to the identity in pi. Substitution cannot bring a
remainder of order greater than7 into one of the first seven centres.
The checker recomputes the full truncated series in (1.1), independently
of the retained literal forward-jet formulas, and verifies each coefficient.

The first triangular pivot equations are

    h1'=h1,
    h2'=h2-3eta,
    h3'=h3+2h1*eta+5mu,

and the coefficient of nu in h5' is -3. Thus there is a unique
transformation taking ANY seven-jet to h2'=h3'=h5'=0:

    eta0=h2/3,
    mu0=-2h1*h2/15-h3/5,
    nu0=-4h1^3*h2/15-11h1^2*h3/15-h1*h2^2/9
          +2h2*h3/3+h5/3.                              (1.2)

There are no divisions by parameters: these are global polynomial formulas
over characteristic zero, including all vanishing and nonreduced branches.

For clarity the group law, with the second transformation applied last, is

    (mu,eta,nu) then (mu2,eta2,nu2)
      = (mu+mu2,eta+eta2,nu+nu2+eta2*mu).

The action is not assumed commutative. Each orbit meets the polynomial
section h2=h3=h5=0 exactly once. Uniqueness follows already from the
three triangular pivots. The normalization and its polynomial inverse
(the inverse target transformation applied to a point of the section)
give a global product parametrization by the section and the group.
Consequently the invariant ring of this seven-jet action is a polynomial
ring in the remaining normalized coefficients I1,I4,I6,I7.

Explicitly, with the numerator signs included,

    I1=h1,
    I4=(6h1^2*h2+24h1*h3+10h2^2+15h4)/15,

    I6=-(126h1^4*h2+414h1^3*h3+96h1^2*h2^2
          -612h1*h2*h3-360h1*h5-100h2^3-360h2*h4
          -189h3^2-135h6)/135,

    I7=(540h1^5*h2+1710h1^4*h3+342h1^3*h2^2
         -2214h1^2*h2*h3-1125h1^2*h5-160h1*h2^3
         -240h1*h2*h4-153h1*h3^2+885h2^2*h3
         +750h2*h5+765h3*h4+225h7)/225.                 (1.3)

Their weights, assigning hj weight j, are 1,4,6,7.

## 2. Five normal forms for the actual F2 contact strata — JC2CFGA1

These transformations preserve the retained coordinate degrees75/125,
both support envelopes, the carrier leading faces, terminal heights3/5,
and the nonzero constant Jacobian. In particular, adding eta*P to g
has smaller ordinary degree, lower eta degree, and lower terminal height
than g. Constants have still lower values. Hence they act within the
same F2 problem, including its actual nonlinear Keller equations.

In the representative (1.2), the first nonzero integer centre is one of
1,4,6,7. If all seven vanish, JC2CFN1 forces the residual band39.
Since the highest eta band is L=75-5j when centre j first survives,
the complete necessary list UP TO this target equivalence is

| first normalized centre | highest band L | cusp contact13L | geometric-degree cap |
| --- | ---: | ---: | ---: |
| I1 | 70 | 910 | 2019 |
| I4 | 55 | 715 | 1899 |
| I6 | 45 | 585 | 1554 |
| I7 | 40 | 520 | 1384 |
| none through7 | 39 | 507 | 1347 |

The degree caps are those already proved in JC2CFMV1. No transformation
of source coordinates and no reconstruction of the missing historical
954-variable circuit is involved.

The original contact values650,780,845 are NOT proved impossible for an
unnormalized Keller map. Their representatives are covered by one of
these five rows. The contact polynomial P^5-g^3 changes under the target
action, so its degree is not an invariant of the target orbit.
The smallest contact507 row remains open and has the forced Wronskian
endpoint established in JC2CFN1.

## 3. Quotient of the retained rational quintic family — JC2CF7F1

For a normalized coefficient presentation

    p=t^3+a*t,
    q=t^5+b*t^4+c*t^2+d*t,

compute the seven carrier centres using pi=p^3/q^2,h=p^5/q^3, and then
apply (1.3). The resulting polynomial map N:A^4 -> A^4 is

    N1=-3b,
    N4=-(50a^2+48a*b^2-27b^4-36b*c+45d)/15,

    N6=(125a^3+9144a^2*b^2+8640a*b^4-216a*b*c+135a*d
         -5049b^6-6480b^3*c+8100b^2*d+81c^2)/135,

    N7=(650a^3*b+39118a^2*b^3+25a^2*c+36816a*b^5
         -1152a*b^2*c+690a*b*d-21627b^7-27612b^4*c
         +34515b^3*d+432b*c^2+45c*d)/75.               (3.1)

The checker independently verifies the input carrier graph through order7
by substitution in z=1/t, and then verifies all four quotient formulas.
Domain weights(a,b,c,d)=(2,1,3,4) make N weighted homogeneous of
degrees(1,4,6,7).

The exact zero-fibre Groebner basis, lex order b>d>c>a, is

    b,
    d+10a^2/9,
    c^2-25a^3/81,
    a^2*c,
    a^5.                                               (3.2)

Thus the complete standard monomials are

    1,a,a^2,a^3,a^4,c,a*c,                               (3.3)

and the zero fibre is a length7 scheme supported only at the origin.
This is the same eligibility scheme found by the full transported
zero-seven-jet calculation JC2CFJ1, now as the special fibre of a global
map with four invariant coordinates.

Here is the universal finiteness and flatness argument, beyond the finite
Groebner calculation. Let A=C[a,b,c,d] and B=C[N1,N4,N6,N7].
Because A/(N1,N4,N6,N7) is finite dimensional and every weight is
positive, induction on weight expresses A as a finite B-module generated
by homogeneous lifts of (3.3). Equivalently this is the usual graded
Nakayama argument, proved by subtracting the residue and reducing weight
in each coefficient of an Ni. The finite extension has dimension4,
so the four Ni are algebraically independent and B is a polynomial ring.

The surjection from the free graded B-module on (3.3) to A is an
isomorphism. Indeed its Hilbert series equals that of A, since

    H_A/H_B
      =(1-t^6)(1-t^7)/((1-t^2)(1-t^3))
      =1+t^2+t^3+t^4+t^5+t^6+t^8.                       (3.4)

The basis weights in (3.3) are 0,2,4,6,8,3,5, exactly that numerator.
A nonzero graded kernel would have a nonzero homogeneous component,
contradicting equality of Hilbert series. This proves that N is finite
FREE of rank7. In particular it is surjective and flat, every
scheme-theoretic fibre has length7, and a generic fibre has seven
distinct geometric points in characteristic zero.

These are seven normalized coefficient presentations over invariant
carrier data. No assertion that they are seven inequivalent polynomial
curves, seven nonproper components, or seven actual Keller solutions
is made.

## 4. Exact ramification locus and an explicit seven-row fibre

For p,q as above the determinant with columns(a,b,c,d) is

    det(dN)=-Res_t(p',q')/15
      =-(25a^4+48a^3*b^2-144a^2*b*c+90a^2*d
           +108a*c^2+81d^2)/15.                        (4.1)

Thus N is etale exactly when the parametrized curve has no finite
nonimmersed parameter. This is a condition on the parametrization:
self-intersection of distinct immersed parameters does not itself force
this determinant to vanish. The same resultant controls the retained
raw-jet transport determinant, with a different nonzero scalar.

For the invariant fibre (N1,N4,N6,N7)=(0,0,0,w), w!=0, (3.1) gives

    b=0, d=-10a^2/9, c^2=25a^3/81, -a^2*c/3=w,
    a^7=729w^2/25, c=-3w/a^2.                          (4.2)

There are exactly seven distinct rows. At each row

    Res_t(p',q')=175a^4/3 !=0.

In contrast the zero fibre w=0 is the nonreduced length7 monomial
endpoint. This explicit contrast prevents replacing the whole finite
map by its special zero-fibre support. In particular N is surjective:
carrier invariant values alone do not exclude the k1 coefficient family.

## 5. What this does and does not close

The exact affine gauge removes redundant cases before deeper Keller
elimination. The quotient map supplies a bounded algebraic interface:
after assigning the four carrier invariants, all k1 coefficient
presentations form a finite flat length7 scheme. This can replace an
unrestricted four-parameter search in a conditional attachment analysis.

It does not prove that any nonproper component follows the generic
carrier profile. PF2PPA1 keeps that global attachment gate open. The
nonzero seventh-centre fibre demonstrates surviving genuine coefficient
freedom rather than an exclusion. A global target-attachment theorem
or a contradiction using actual deeper source Keller coefficients is
still required, especially on the lowest contact507 row.

The exact checker imports only a pinned arithmetic literal list from the
retained source, validates it by independent truncated-series identities,
and certifies the small polynomial ideals and Hilbert identity. The
universal group, graded-module, and geometric implications are proved
in this note. There is no formal proof-assistant verification,
independent replay, external review, or novelty-priority claim.

