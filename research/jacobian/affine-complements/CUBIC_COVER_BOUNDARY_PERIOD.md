# Cubic-cover boundary periods and the pole(6,9) frontier

14 November 2026. PF2RULED18 extends the leading-period method to every
pole pair of gcd3. On the auxiliary normalized cover zeta^3=h of B,
the intrinsic one-form e(w)*dw/zeta^2 must be exact. This is an exactness
condition on that cover, not a rational-primitive condition on B.
For noncube h it is equivalent to a small polynomial linear system.

The subsequent [PF2RULED19 generic-polar and marked cusp constraints](POLAR_CUSP_LEADING_SECTORS.md)
exclude simple common roots, contact1/contact8 and the pole(6,9)72,108
case. The five numerical rows below are the earlier PF2RULED18 handoff;
its cover-period classification and positive meromorphic controls remain
valid. Current source groups are96/144 at contact0,80/128 at contacts
3,4,5,6,7,9, and48/128 at contacts7,9.

At pole(6,9) the result also gives the complete necessary ordinary-degree
cases below and excludes contact multiplicity2 of the common nonic.
It does not close the whole pair or produce a finite map.

| Contact multiplicity k of common nonic h | Ordinary source degrees(P,Q) | Cusp multiplier slope |
| --- | --- | --- |
| 0 | (96,144) | 3/4 |
| 1 | (88,132) | 8/11 |
| 3,4,5,6,7,8,9 | (80,128) | 3/4 |
| 7,8,9 | (48,128) | 3/4 |
| 3 | (72,108) | 2/3 |

Each row must additionally satisfy the period restrictions below. No
row is asserted realizable, even as a full nonlinear common-symbol pair.
The first open reduced pairs remain(5,9),(6,9),(7,9),(8,9), with the
other three narrowed by [PF2RULED17](HEIGHT_NINE_RESONANT_FRONTIER.md).

## The correct cover and the formal critical coefficient

Use the compatible F1 geometry, L=O_B(B)=O_P1(3), and volume divisor
E-3B, with E|B=4q. A gcd3 pair of poles m,n has leading symbols
h^(m/3),h^(n/3), up to nonzero scalars, for h in H0(B,L^3). Write
locally B=(z=0), Omega=(e(w)+O(z))*z^-3 dz wedge dw.

Assume first that h is not a cube in C(B). Then the normalization of
zeta^3=h is connected. A rational homogeneous term of pole j commuting
with P has coefficient proportional to h^(j/3); it is rational only
when3 divides j. The formal root W=P^(3/m) has pole3 and rational
leading coefficient h. Subtract the finitely many permitted terms
W^(j/3), including negative powers, down through all j>2-m. This is
the same formal descent as PF2RULED11, with the divisibility changed
to3; the auxiliary root need not be a global regular function.

At the first unforced order one obtains a rational rho in C(B) with

 Q-F=z^(m-2)*rho+O(z^(m-1)),
 e=-m*h^(m/3-1)*(h*rho'+((m-2)/3)*h'*rho).

On the cover, set G=zeta^(m-2)*rho. Then

 dG=-(1/m)*e(w)*dw/zeta^2.                              (1)

Thus eta=e*dw/zeta^2 is an exact meromorphic differential. Under a
normal-frame change h_new=u^3*h, zeta_new=u*zeta and e_new=u^2*e,
the differential is unchanged. The construction is intrinsic to the
marked boundary data, up to the harmless constant cube-root choice.

Conversely, if eta has a meromorphic primitive on the cover, project
that primitive to the deck character zeta->epsilon*zeta. The projected
primitive still differentiates to eta and has character epsilon,
because eta has character epsilon^-2=epsilon. Dividing it by
zeta^(m-2) gives a rational function on B, since m-2=1 mod3. Therefore
cover exactness is equivalent to solvability of this first critical
coefficient equation. It is necessary for P,Q, not sufficient for the
remaining coefficients or any global finite map.

This auxiliary degree-three cover is not an intermediate field of the
actual quartic map. No root is put in the target algebra, and no field
divisibility or S4 contradiction is inferred.

## A complete polynomial test for a proposed noncube nonic

Choose a coordinate on B with q at infinity. The leading volume becomes
a nonzero constant and h becomes a polynomial of degree D<=9. Its
contact multiplicity at q is k=9-D. This is an intrinsic coordinate
choice for the period test, not an asserted automorphism of the marked
surface or source. Scaling the constant differential does not affect
exactness. Let s=rad(h) be monic and r=deg(s), the number of distinct
finite roots, with multiplicities k_i.

Normalize an eigencharacter primitive as3*zeta^4*rho_6, with
rho_6 rational. At a finite root of multiplicity k_i, the rational
critical equation forces rho_6 to have pole order exactly2k_i-1;
there can be no other finite pole. The derivative term at an additional
pole would dominate and could not equal the nonzero constant. At a
root, the forced pole has coefficient proportional to3-2k_i, which
never vanishes for integer k_i. A possible homogeneous order4k_i/3
is either nonintegral or smaller than2k_i-1 and cannot supply the
constant term. Therefore write

 rho_6=A*s/h^2, A a polynomial.

Direct differentiation proves that exactness is equivalent to

 3*s*A'+B*A=1,
 B=3*s'-2*s*h'/h.                                      (2)

B is polynomial. If(2) holds, an explicit primitive of dw/zeta^2 is

 3*zeta^4*A*s/h^2.

Conversely, every primitive of the indicated character yields such A.
Thus(2) is a necessary and sufficient exactness test for the cover
period, including positive-genus covers. It is not a map-realization
test. For noncube h the solution, if it exists, is unique: the difference
would be a rational multiple of h^-4/3, impossible unless h is a cube.

For r=1, equation(2) has A=1/(3-2D), after h=(w-a)^D is normalized.
For r>=2, let M=deg(A). The highest coefficient of the left side has
factor3M+3r-2D. Since M+r-1>0, this must vanish. Hence

 D in{3,6,9}, M=2D/3-r>=0.                             (3)

After this top cancellation, the left side has degree at most
2D/3-2. For D=9, the whole exactness test is thus at most five
coefficient equations, linear in at most five coefficients of A.
For D=6 it is at most three equations; for D=3 it is a constant
identity. For a fixed h this is a finite exact linear decision problem,
not a numerical period integration.

## A sharper root budget on the cover

For r>=2 and D divisible by3, the three points over infinity are
unramified. Equation(3) makes the primitive finite and nonzero there;
its three values differ by cube roots of unity. At each, the primitive
minus its value vanishes to order

 ell=2D/3-1.

At a root of h of multiplicity k_i, put g_i=gcd(3,k_i). There are g_i
points over it, with base ramification3/g_i. The primitive has order
(3-2k_i)/g_i: it vanishes at a simple root, and for k_i>=2 contributes
total pole degree2k_i-3 across these points. Thus its total pole degree is

 d_G=2D-3r+s1,

where s1 is the number of simple roots. A zero of G-G(infinity) cannot
exceed the entire pole degree. Consequently

 2D-3r+s1 >=2D/3-1.                                    (4)

This excludes more partitions than(3) alone. At D=3,r2 it leaves only
(2,1); at D=6,r4 only(3,1,1,1); at D=9,r6 only(4,1,1,1,1,1).
For D9,r5 it excludes(2,2,2,2,1). These are partition obstructions,
not an assertion that every remaining root configuration passes(2).

## Complete low-degree root-position classifications

The following statements concern noncube h; cube sections are treated
separately below. Positions may be normalized by affine changes on the
abstract B with q=infinity. All constraints are intrinsic to the root
configuration; the actual marked point of the ruled surface is not
moved or forgotten in a source comparison.

For D3, the only multi-root case is(2,1), and every two distinct positions
pass(2). For D6, the two-root partitions(5,1),(4,2) pass at every pair
of distinct positions. The full three-root cases, with roots0,1,t and
multiplicities in that displayed order, are exactly:

| Multiplicities | Necessary and sufficient condition, t!=0,1 |
| --- | --- |
| (4,1,1) | 2t^2+t+2=0 |
| (3,2,1) | 2t+1=0 |
| (2,2,2) | t^2-t+1=0 |

To prove completeness, put A=a1*w+a0 in(2). The constant equation fixes
a0, the quadratic coefficient fixes a1, and the remaining linear
coefficient gives precisely the displayed polynomial. The possible
denominator exception t=-1 in the first and third row contradicts
the quadratic coefficient, so it loses no solutions. All actual roots
of the displayed criteria satisfy the distinctness and denominator
guards. The checker retains every solved coefficient and residual.

For D6,r4, the root budget leaves(3,1,1,1). With the triple root at a,
write x=w-a and h=x^3*S(x), S monic cubic with S(0)!=0. Here A is
constant and B=x*S'-3*S. Thus B can be constant only when

 h=x^3*(x^3+c), c!=0.                                  (5)

This family passes with A=-1/(3c). Together with the single-root
cases and the cube criterion below, this classifies every degree-six
period configuration.

For D9,r6, the same argument gives the complete extremal family

 h=x^4*(x^5+c), c!=0, A=-1/(5c).                        (6)

For D9,r2, every noncube two-root partition(8,1),(7,2),(5,4) passes.
For the remaining degree-nine root-position loci, equations(2)-(4)
are retained as the complete decision test; their full parameter
classification is still open here.

Positive controls are essential. For example h=(w^3-1)^2 passes with
A=-w/3. On t^3=w^3-1, zeta=t^2, its primitive is -w/t, an exact
differential on the nonsingular cubic cover. Family(6) is a genuinely
noncontact nonic with six finite roots and an exact cover period.
Thus neither positive genus nor multiple noncontact roots alone exclude
the gcd3 layer. These controls construct no P,Q or finite quartic map.

## Cube sections and the actual finite-quartic input

If h=g^3, its cubic g is governed by [PF2RULED11](BOUNDARY_LEADING_PERIOD.md)'s rational period
criterion, applied to P leading g^m. There is at most one root of g
away from q. For an actual finite quartic target, [PF2RULED14](GLOBAL_CUSP_LEADING_OBSTRUCTION.md) excludes
the one-root noncontact case. The cube branch can therefore only have

 div(h)=6p+3q, 3p+6q, or9q, p!=q.                      (7)

This step uses the actual target-algebra cusp theorem. It must not be
extended to noncube h: the positive controls above show why it would fail
as a leading-period argument.

## The pole(6,9) ordinary-degree cases

For k=0, the source-degree theorem is exact at contact order zero and
gives(96,144). Cube h is already excluded there by(7); a noncube h
must pass the D9 test and root budgets. For k=1, the two component
symbol orders are2 and3. The exact contact-sensitive transfer gives
(88,132). The period condition has D8, so h has a single finite root:
div(h)=8p+q. Their chart-leading forms are v^6*(u+v)^16 and
v^9*(u+v)^24. The cusp-multiplier slope is8/11.

For k>=2, the source bounds are deg(P)<=80, deg(Q)<=128. The pinned
2022 alternative gives either exact Q degree128 or the normal72,108
exception; the reverse is impossible because P has degree at most80.
The complete Q source basis contains no degree125,126,127 row.

For Q128, the full contact-order>=6 subspace of W_(9,128), dimension139,
has highest chart form c*T^8, T=v*(u+v)^3. The complete W_(6,80) has
chart degree<=20. The leading Jacobian therefore forces P to lead as
T^j,1<=j<=5, with source degree16j. The source-gap and target-shear
argument of PF2RULED17 excludes j1,2,4, leaving P48 or80. Complete
constrained kernels retaining the actual pole-six symbol give:

| Source degree P | Chart leading form | Dimension | Minimum symbol order |
| --- | --- | --- | --- |
| 48 | T^3 | 43 | 13 |
| 80 | T^5 | 65 | 5 |
| 72, exceptional branch | R^6, R=v*(u+v)^2 | 61 | 6 |

An actual P symbol is h^2, with order2k. Therefore the128 branch
permits P48 only for k>=7, and P80 only for k>=3. In particular k2
is excluded in this branch.

For the72,108 exception, every W_(9,108) symbol has order>=8. Since
Q has order3k, this already gives k>=3. Its full contact-order>=9
high projection is the quadratic one of PF2RULED17. The same factor
argument, now with partner chart degree18 or19, retains only R^9,R^6.
In the entire partner kernel one has

 [r^6]symbol(P)=[v^18]P0!=0.

Thus2k=6, so k=3 exactly. This excludes k2 in every numeric alternative.
The surviving exception has slope2/3; the128 branch has slope3/4.
For k0 the slope is3/4. These slopes use the [actual PF2OMIT9 cusp
identity](../jc2-omitted-divisors/WRIGHT_IMAGE_RESTRICTION.md) and retain its nonzero value at the cusp; they are not sufficient
Jacobian or finiteness conditions.

The [checker](verify_cubic_cover_period.py),
[certificates](cubic_cover_period_certificates.json) and
[receipt](cubic_cover_period_receipt.json) verify the exact polynomial
identities, full low-degree cross-ratio conditions, every partition
in degrees1..9, positive/negative period controls and complete section
kernels for the numerical table. Formal descent, normalized-cover
interpretation and the divisor degree budget are written proofs;
the source-degree frontier also uses the already retained2022 theorem.
No formal verification, independent full reproof, conditional2026 upgrade,
external review or novelty priority is claimed.

The remaining degree-nine period-position loci, the full nonlinear
Jacobian/cusp-unit/finiteness problem on every retained row, higher heights,
all-height pole-four coordinates and actual compatible embedding of an
arbitrary quartic normalization remain open. The original full goal
stays active; neither full quartic packet is excluded.

The exact noncontact six-root controls above cannot be finite quartic
realizations under PF2RULED19's generic-polar criterion: their simple
common roots violate2*k>d. They remain positive period controls and
make clear why a period solution is not a finite map. The current
finite-map handoff and sixteen marked two-root leading candidates are
in [PF2RULED19](POLAR_CUSP_LEADING_SECTORS.md).
