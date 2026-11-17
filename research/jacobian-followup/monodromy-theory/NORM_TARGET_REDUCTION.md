# A norm double cover in every degree-eleven marked descent

This investigates the rational-source relative construction of KMFI80,
rather than repeating its actual-component unit obstruction. For the
Klein degree-eleven input, every dominant base change selecting an
eleven-sheet marked component must contain one explicit double-cover
field. A candidate top form on that field is regular away from the
two finite branch fibers, including at infinity. The finite-fiber gate is now completed by
[KMFI92--93](NORM_TARGET_FINITE_FIBERS.md), which prove geometric genus3
and exclude the specified dominant rational-target descent.
The narrower local proof and its checkpoint are retained here.
The exceptional action itself remains open.

Work over algebraically closed characteristic-zero k. Fix a monic
degree11 polynomial f with f(0)=0, f'(0)=1, geometric group
G=PSL2(11) in either exceptional action, and passport
\(2^4 1^3|3^3 1^2|11\), with finite branch values0,beta, with beta nonzero.
Choose an unramified point of the first fiber as0; two scalings can
make f monic and f'(0)=1 over k, while the second branch value beta
is retained. These three normalizations do not generally also allow
beta=1. Let c denote the T^10 coefficient.
The arguments below apply to any f satisfying these hypotheses.

## 1. The rational relative model and its full colored group

Put K=k(b,q,r), let t satisfy f(t)=b, and mark it by

\[
 a^{12}=f'(t)Q(t),\qquad
 Q(T)=1+rT^2+q(T^3-cT^2).
 \tag{1}
\]

Solving q gives the rational source field k(t,a,r). As shown in KMFI80,
this is an actual component of a relative marked construction with
n=14, m=12. Its original affine component has a nonconstant unit;
rationality of the field does not make it an affine-space source.

Let M/K split f(T)-b, label its roots t_i, and set

\[
 N=M(a_1,\ldots,a_{11}),\qquad a_i^{12}=f'(t_i)Q(t_i).
\]

Over M, the eleven cofactor linear polynomials Q(t_i) in q,r are
pairwise distinct, with constant term1. Each has an individual
value-one Kummer divisor not shared by another radicand.
Consequently

\[
       \operatorname{Gal}(N/K)=A\rtimes G,\qquad
                 A=(\mathbf Z/12)^{11}.
 \tag{2}
\]

G permutes the roots and the a_i; independent Kummer monomials show
this is a genuine section. This reproduces the relative group
calculation, not a new ordinary-affine occurrence theorem.

## 2. A compulsory norm subfield in every eleven-sheet descent

**KMFI90.** Let K'/K be any characteristic-zero field extension.
Suppose a connected component of the base-changed marked extension
(1) has degree11 and its underlying root action is still the full
exceptional G action. Then K' contains the nontrivial quadratic field

\[
 E=K(Y),\qquad
 \boxed{Y^2=R(b,q,r)=\operatorname{Res}_T(f(T)-b,Q(T)).}
 \tag{3}
\]

Here T is the polynomial variable. All complements
of A in (2) are A-conjugate; there are exactly12^10 of them.
The same embedded field E is fixed by every complement.

To see the subgroup necessity, embed N into NK' and write
J=Gal(NK'/K'), viewed as a subgroup of A semidirect G.
Its projection to G is surjective by the root-action hypothesis.
An eleven-point marked orbit contains exactly one color over every
root label. An element of J intersect A must fix that chosen color
at every label, hence all its coordinates are zero.
Thus J intersect A=1 and J is a complement isomorphic to G.
The fixed field N^J=N intersect K' is contained in K'.

Let U=product a_i. Since f has the stated passport,

\[
 \prod_i f'(t_i)=\delta(b)^2,\qquad
                 \delta(b)=\lambda b^2(b-\beta)^3,\quad\lambda\in k^*.
\]

Equivalently, the discriminant is square since G is contained in A11.
Then \(Y=U^6/\delta(b)\) satisfies (3).
For any complement, its color-vector sum is a homomorphism
G -> Z/12. The group G is perfect, so that sum is zero.
It fixes U and Y. Therefore E is contained in N^J for **every**
complement, not only the section with zero colors.
R is not square: over M it is a product of distinct cofactor
linear forms, each with valuation1.

For completeness, complement conjugacy also has an elementary proof.
Use colors as functions on the eleven labels; a complement has
a cocycle z(g) in A, with z(gh)=z(g)+g*z(h).
Fix a label0 and its stabilizer H=A5. The coordinate z(h)_0 is a
homomorphism H -> Z/12, zero since H is perfect.
Choose g_i with g_i(0)=i and put v_i=z(g_i)_i.
Writing \(g g_i=g_{g(i)}h\) shows

\[
                      z(g)_{g(i)}=v_{g(i)}-v_i.
\]

Thus z(g)=v-g*v. Translation by v conjugates the complement to the
zero-color section. Its stabilizer in A consists of the12 constant
vectors, giving12^11/12=12^10 complements.
Both exceptional A5 actions are included.

This theorem concerns dominant field base changes of (1), which
contain K. It does not cover a non-dominant specialization replacing
K by a quotient. Forgetting the mark and keeping the unmarked
eleven-sheet root field is also a different operation.

If E has a nonzero regular three-form on a smooth projective model,
then E is not unirational: a dominant rational map from projective
space would pull that form to a nonzero regular three-form on a
resolution of projective space. No such form exists.
A rational extension field K' containing E is therefore impossible,
even with additional independent variables. Establishing this
premise is supplied by KMFI92 in the follow-up note; (3) alone does
not prove it.

## 3. A top form with only two unresolved fibers

Homogenize Q as

\[
 Q_h(T)=h+rT^2+q(T^3-cT^2),\qquad
 R_h=\operatorname{Res}_T(f(T)-b,Q_h(T)).
\]

The resultant has b-degree3 and (h,q,r)-degree11.
Compactify E as the double cover of
P1_b times P2_[h:q:r] with fiber equation

\[
                       Y_h^2=h R_h .
 \tag{4}
\]

More precisely, use [b0:b1] on P1 and
Rbar=b0^3 R_h(b1/b0). The normal double cover has equation
Z^2=b0*h*Rbar, of bidegree(4,12), with Z in O(2,6).
Its dualizing sheaf is the pullback of O(0,3).
A dualizing section on this singular model need not extend to a
resolution at noncanonical singularities.

In the affine chart h=b0=1 consider

\[
                 \boxed{\eta=q^3\,db\wedge dq\wedge dr/Y.}
 \tag{5}
\]

**KMFI91.** This nonzero rational three-form extends over a resolution
of (4) above P1 minus {0,\beta}. In particular **infinity is resolved**.
The assertion does not include the fibers over0,beta; it is not a
global differential-form or non-unirationality theorem.

First consider finite b away from0,beta. Over an etale splitting
neighborhood, the branch divisor in P2 consists of h=0 and the
eleven distinct lines Q_h(t_i)=0. A generic fiber has no triple
intersections, as follows.
For three finite roots the determinant, apart from nonzero pair
differences, is
\(t_i t_j+t_i t_k+t_j t_k\). If it vanished, the three reciprocals
would sum to zero in the splitting field. Averaging its G-conjugates
gives (3/11) times the sum of all reciprocals equal to zero.
But logarithmic differentiation of f(T)-b at0 gives
\(\sum_i 1/t_i=1/b\), a contradiction.
No two root lines meet on h=0, because the corresponding determinant
is \(t_i^2t_j^2(t_j-t_i)\ne0\).

This generic argument must not be used at every specialized b:
ordinary triple intersections can occur at additional finite values.
There can never be four root lines through one point, because a
nonzero polynomial Q_h of degree at most3 cannot vanish at four
distinct roots. Away from b=0 no root is zero, and h=0 meets the root
lines at distinct points by the determinant just computed.
Each fiber away from0,beta therefore has only ordinary double or
triple branch-line intersections.

At a double intersection the transverse surface section is A1;
at a triple intersection it is D4, the double cover of a plane
branched along three distinct concurrent lines. The total threefold
is a compound Du Val deformation of that section and is canonical.
We use [Corti--Haskins--Nordstrom--Pacini, Definition3.25 and
Proposition3.26](https://msp.org/gt/2013/17-4/gt-v17-n4-p02-p.pdf),
which state precisely this cDV-to-canonical implication.
Consequently a regular dualizing top form on these charts extends
to a resolution. This attributed singularity input is necessary
for the possible ordinary triple fibers; their existence is not
silently excluded.

For double intersections one can also check extension directly.
At a line
intersection the double cover has local equation z^2=xy with one
extra smooth parameter. The form dx wedge dy/z extends to its
resolution: in a chart x=u, y=u v^2, z=u v it is2 du wedge dv.
The other chart is identical. Across the projective fiber charts,
a degree3 numerator is exactly regular: the branch degree is12,
the covering coordinate has weight6, and the base plane top form
has weight-3. This verifies (5) on the whole proper fiber, not only
on the h=1 affine patch.

At infinity put b=epsilon^(-11), q=epsilon^3 Q, r=epsilon^2 R.
The roots are t_i=epsilon^(-1) z_i(epsilon), with
z_i(0) the eleven distinct eleventh roots of unity; this follows
formally from monicity and Hensel lifting. Then

\[
 R(b,q,r)=\prod_i
       [1+R z_i^2+Q(z_i^3-c\,\epsilon z_i^2)].
\]

Use the corresponding projective P2 family in [h:Q:R], and keep
the affine covering coordinate Y. Direct differentiation gives

\[
 \eta=-11\,\epsilon^2 Q^3\,
           d\epsilon\wedge dQ\wedge dR/Y .
 \tag{6}
\]

At epsilon=0 its twelve branch lines are h=0 and
\(h+R\zeta_i^2+Q\zeta_i^3=0\).
A triple of finite lines could be concurrent only if three distinct
eleventh roots of unity summed to zero after taking reciprocals.
That is impossible: a rational polynomial with three unit
coefficients and degree at most10 cannot be divisible by
Phi_11=1+X+...+X^10. Intersections with h=0 are again distinct.
Thus the new projective family has only the resolved z^2=xy charts,
and (6) is regular everywhere above infinity, with a factor epsilon^2.

Finally regularity descends under the finite base extension and the
birational change of proper model. A pole of integer order a>=1
at a divisor pulls back with order -e*a+(e-1)<0 in characteristic
zero. It cannot become regular after a generically finite cover.
A common resolution of the proper models supplies the needed
divisors. This proves the infinity assertion.

### Retained finite-fiber gate, now completed

At b=0 four pairs of root lines coalesce. The simple tail t=0
also meets the extra branch line h=0 with higher contact.
At b=beta three triples coalesce. Their intersections with each
other and with the remaining tails were the final resolution gate.
KMFI92 now resolves them using a finite ordered-root cover and
crepant complete-intersection blowups. Fiberwise adjunction, the epsilon^2 calculation,
or regularity at a generic point of a finite fiber is insufficient.

A useful choice of the tail normalized to0 can ensure that the
cubic critical factor of the order3 fiber has nonzero derivative
there: it has only two derivative roots, whereas the first fiber
has three simple tails. This prevents Q_h from vanishing at all
three triple-critical points at once. It does not by itself
complete the resolution by itself. The new proof does not need
this restriction: it includes all three triple-critical points.

The local checkpoint model records only the KMFI90--91 proof scope.
Current global extension, geometric genus3 and dominant-descent
exclusion are in KMFI92--93. Exceptional affine occurrence and the
full classification remain **open**.

## Reproduction and assurance

Run

    .venv/bin/python research/jacobian-followup/monodromy-theory/verify_norm_target.py

The checker controls both exceptional actions and perfect stabilizers,
colored complement orbits, norm/resultant identities, projective
form transitions, the full infinity scaling and all cyclotomic line
determinants. An algebraic degree11 fixture is used only to test the
resultant identity; it is not called a Klein polynomial or a PSL
realization. The all-field arguments are written above.
The earlier direct m9 derivative-only exploration is retained
separately and is not confused with this rational-source m12 model.
No formal verification or independent whole-proof review is claimed.
