# Incidence, a finite cusp gate, and the nonimmersed singleton obstruction

14 November 2026. Three refinements of PF2R7B5REC1 are established here.

* PF2R7B5GEN1: on a nonempty open in each first-image family, all seven
  high-contact quintics have four ordinary nodes, avoid singularities
  of C, and meet C in twelve transverse affine points. C also avoids
  their singularities. These local tests permit the immersed rows C
  and A with a shared unramified image, rather than excluding them.
* PF2R7B5CUSP1: a nonimmersive high-contact quintic occurs only over
  exactly160 normalized first-image coefficient points, at most54
  cube-root orbits. This is a necessary finite gate for excess-one
  row B, not a realization of any of those points.
* PF2R7B5N1: the nonimmersed oriented row N has no singleton fibre if
  its additional image has same-coordinate polynomial normalization
  degrees(3,5) and infinity contact at least38.

All statements retain those explicit degree and contact conditions on B.
They are authored written geometry with exact finite-algebra checks,
without formalization or independent external review. No global Keller
cover is constructed, and JC2 remains open.

## 1. The seven-choice algebra and the nonimmersion function

Use the explicit rank-seven coefficient algebra and unique normalized
parametrizations of PF2R7B5REC1:

    P(s)=s^3+A*s+M,
    Q(s)=s^5+b*s^4+c*s^3+e*s^2+f*s+n,
    r9=r8=0.

The five lower coefficients are fixed polynomials in A,M and the first
parameter d, or z=x^3. Its monic weighted Groebner basis is valid over
the entire parameter polynomial ring; all parameter fibres have length7.
This counts normalization-coefficient multiplicity, not distinct curves
at every parameter.

The exact critical resultant is

    Delta=Res_s(P',Q')
      =25*A^4+48*A^3*b^2-144*A^2*b*e
       +90*A^2*(f-c*A)+108*A*e^2+81*(f-c*A)^2.       (1)

The checker recomputes the resultant. Nonimmersion is exactly Delta0.
Multiplication by Delta in the rank-seven algebra has an explicit
seven-by-seven matrix over QQ[d] or QQ[z]. Its determinant is zero
precisely when some closed normalization-coefficient point has Delta0:
an element of a finite algebra is a unit if and only if it is nonzero
at every closed point. Nilpotent fibres cause no exception to this test.

## 2. The exact finite excess-one gate

Exact matrix determinants give

    Norm_L(Delta)=k_L*d^2*F12(d^3),
    Norm_K(Delta)=k_K*z^15*G41(z),                   (2)

where k_L,k_K are nonzero rationals. The primitive integer polynomials
F12 and G41 are supplied in the two nonimmersion-parameters JSON files;
the matrices, scalars and literal determinant polynomials are in
[incidence.json](incidence.json). The checker reconstructs the matrices,
recomputes their determinants and verifies both factorizations exactly.

F12 and G41 are squarefree, have degrees12 and41, and are nonzero at0.
Thus L has exactly37 coefficient parameters with a nonimmersive choice:
d0 and the36 roots of F12(d^3). K has exactly41 nonzero z values, or123
x parameters. The three shared L/K coefficient points have d^3=x^3=-2/3;
the factors in(2) are nonzero there. The sets therefore do not overlap:
there are160 coefficient points in their normalized union.

The cube-root action identifies the36 nonzero L points in12 orbits,
with d0 one further orbit. The123 K points give41 orbits. Hence at
most54 corresponding first-image types are needed here. No assertion
of54 distinct affine moduli or complement types is made.

Exact gcds show F12(d^3) is coprime to every exceptional affine factor
in PF2R7COLL1, and G41 to the cusp, triple and tacnode factors in
PF2R7KCOLL1. Therefore159 of these160 first images have twelve nodes.
The remaining one is d0, with three triple points and three nodes.
All160 first images are immersive.

In residual row B, PF2OMIT6/PF2OMIT7 require the additional image to
have one ordinary cusp. Thus any row-B image in this degree/contact
range must lie over this finite set. Conversely norm zero provides a
nonimmersive polynomial image, but does not prove it has exactly one
ordinary cusp, the required incidences, or a Keller realization. The
finite gate is necessary, not sufficient. This proves PF2R7B5CUSP1.

## 3. Exact simultaneous witnesses for the generic incidence statement

Use d1 and z1. At each, PF2R7B5REC1 gives a squarefree degree-seven
eliminant h(A) and a polynomial formula M=m(A). Work in E=QQ[A]/h;
E includes all seven candidates without choosing numerical roots.

For each candidate the unordered self-collision quartic is

    R(v)=v^4+b*v^3+A*v^2+(2*A*b-e)*v-(A^2+f-c*A).

Its pair parameters satisfy s1+s2=v and s1*s2=v^2+A. A diagonal pair
is detected by 3*v^2+4*A. Its common target coordinates are

    U_B=M-v*(v^2+A),
    V_B=(v^2+A)*(v^3+2*A*v+A*b-e)+c*(U_B-M)+n.       (3)

The exact matrices below are multiplication matrices in explicitly
monic finite-free quotients. Their ranks are recomputed modulo101.
Every listed matrix is square and has full rank, for each of d1,z1.

| Element tested for invertibility | Algebra dimension | Consequence |
|---|---:|---|
| h'(A) | 7 | Seven distinct normalization choices |
| R'(v), separately 3v^2+4A | 28 each | Four simple off-diagonal collision pairs |
| A^2+A*b^2-b*e-A*c+f | 7 | No triple value of B |
| F_B(p_C(t),q_C(t)) modulo the first-image collision polynomial D_C(t) | 168 | B avoids Sing(C) |
| q_C(t)-V_B modulo h,R(v),p_C(t)-U_B | 168 | C avoids Sing(B) |
| Leading t^12 coefficient of F_B(p_C,q_C) | 7 | All seven finite intersection polynomials have degree12 |
| t derivative of the monic finite intersection polynomial | 84 | All twelve intersections are simple |

Here D_C is recomputed as the monic resultant in the other variable of
the two divided differences. It has degree24 and is squarefree at both
witnesses; the normalization is immersive. The already proved first-image
strata identify these members as twelve-node curves.

The finite-free hypothesis is essential to the modular inference.
All input rational denominators are nonzero modulo101. The displayed
relations are monic in successive variables: h has degree7 in A, R
degree4 in v, and p_C-U_B degree6 in t; similarly D_C is monic degree24.
Thus these quotients have the stated rational vector-space bases, and
the reduced matrices really are reductions of their characteristic-zero
multiplication matrices. A nonzero determinant modulo101 is a nonzero
rational determinant. It excludes every common zero over C.

For the degree12 intersection algebra, its leading coefficient is first
proved a unit in E with determinant nonzero modulo101. Its rational
inverse exists and has reduction equal to the modular inverse used by
the checker. Dividing by that unit gives a monic degree12 relation and
a rank84 finite-free algebra. This accounts for the only normalization
step involving a nonconstant unit. It is not an inference from an
arbitrary modular unit ideal, where solutions could escape to infinity.

Simple roots of R give transverse normalization pairs, since the
unordered-pair coordinate change is etale off the diagonal. Two different
pairs with the same target value would give three roots of the cubic
P(s)-U_B. Reducing Q modulo s^3+A*s-k shows this happens only if
A^2+A*b^2-b*e-A*c+f=0. Therefore each B has four ordinary nodes.
Mutual singular-locus avoidance and simple intersection roots give
twelve distinct transverse affine intersection values.

The conditions are algebraically open. The coefficient correspondence
is finite over the first-image base, so an open property holding for
all seven points of a fibre holds over a base neighbourhood after
removing the finite projected bad locus. The collision and intersection
algebras above remain finite after restricting to the open where their
leading units are nonzero. Hence both exact witnesses establish a
nonempty Zariski open on their entire respective first-image families,
on which all seven candidates have the stated properties.

## 4. Interpret the generic witnesses in the omitted-prime ledger

On these opens, C has twelve nodes, B four nodes, and they meet at
twelve smooth transverse points, with infinity contact38. There are
no other affine branch strata. The complete retained counts are:

| Target stratum | Row C | Row A with H=C | Row A with H=B |
|---|---:|---:|---:|
| Off both curves | 6 | 6 | 6 |
| Smooth C only | 4 | 3 | 4 |
| C node | 2 | 0 | 2 |
| Smooth B only | 3 | 4 | 3 |
| B node | 0 | 2 | 0 |
| Smooth C/B intersection | 1 | 1 | 1 |

These follow directly from PF2OMIT6/PF2OMIT11, counting omitted
normalization branches per prime. Empty source fibres are allowed.
In each of these three proposed rows there would be exactly twelve
singleton target values, all at the displayed intersections. Every
count is nonnegative, so these local branch tests do not exclude the
generic candidates. Row B, however, requires a cusp and is absent on
this open. H distinct still requires its actual equation and branches.
This proves PF2R7B5GEN1; it is a local compatibility statement, not a
finite covering or source-A2 construction.

## 5. The complete nonimmersed first-image classes

Every nonimmersed C coefficient point is represented, after the known
cube-root scaling, in the scaled K chart by one of

    z=-2/3, t*=1,  C cusp=(-1/9,-98/59049),
    z= 4/15,t*=-1, C cusp=(1/675,10991408/184528125).

These are the class-I and class-II points of PF2R7FAM2/PF2OMIT13.
The checker recomputes the coordinates, the unique common derivative
root, the nonzero second/third jet determinant and the clean double
normalization fibre. No cusp corank is counted as two fold branches.

At each class form the element F_B(U*,V*) in the same rank-seven
coefficient algebra, without assuming that its fibre is reduced.
Its exact seven-by-seven multiplication determinant is nonzero in
both cases; matrices and rational determinants are retained in
incidence.json. Hence every one of the seven candidate images avoids
that clean cusp target. The nonimmersion norms in Section2 are also
nonzero here, so these B choices are immersive as row N requires.

PF2OMIT13's oriented row N has two omitted primes and no H. Its retained
source cardinality is

    6-2*c(y)-2*b(y)-epsilon(y),

with epsilon1 only at the known clean cusp of C. At that value c1,b0,
so the fibre has three source points. Everywhere else the count is
even. Consequently no singleton fibre exists anywhere in row N for
any of these high-contact degree(3,5) B choices. Scaling back to the
six normalized coefficient points and any original target coordinates
preserves incidence and fibre cardinality. This proves PF2R7B5N1.

It excludes these configurations as explanations of a singleton fibre;
it does not exclude a Keller map merely because it has no singleton.
The actual extra image is not forced to have these degrees or contact.

## 6. Remaining global work and failed discovery route

The generic local incidence gate passes; a global S6 complement
representation, all omitted image equations including any H, boundary
attachment and identification of the source A2 remain necessary.
The160-point cusp gate still needs its actual ordinary-cusp and
incidence classification. Other contact ranges and target charts remain.
These results do not close unrestricted singleton-zero properness,
the whole r7 family, other terminal slices, or JC2.

Two direct characteristic-zero reverse-incidence Groebner invocations
had45-second bounds and returned terminal halt output without a
completed basis. They are retained as failed discovery. The canonical
proof uses the explicit finite-free algebras and matrix ranks instead;
neither capped run nor printed solver unit is a premise.
