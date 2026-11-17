# No missing branches over a tautological sweep's ramification surface

The earlier finite-sweep obstruction left a precise escape: a nonfinite
polynomial sweep might have additional ramified normalization sheets
over the same surface as its visible ramification. This note closes
that escape without a degree bound on any coefficient, primitive pair
or scalar, and without a multiplier-chart assumption.

Work over an algebraically closed characteristic-zero field. Let
\[
 S(\gamma,w,z)=X(w,z)+\gamma(1,w,z),\qquad
 X\in k[w,z]^3,\quad \det DS=\gamma^2,\quad \operatorname{rank}DX=2 .
                                                               \tag{1}
\]
Put D=closure(X(A2)), and delta=[k(w,z):k(X)].
KMFI32 proves delta<=2 and shows that Delta=(1,w,z) is an
asymptotic tangent direction to D at X(w,z).

## 1. Horizontal line rulings force the finite polynomial chart

Call a line in the affine target horizontal if its Y1 coordinate
is constant. Saying D is covered by such lines concerns actual
straight lines, not birational ruledness.

**KMFI100.** If D is covered by horizontal lines, then, after an
affine source change in w,z and the matching affine target change,
the nilpotent normalization multiplier is a nonconstant polynomial
in one source coordinate. In particular S is finite by KMFI44.

Here is a proof covering reducible generic plane sections and
nonlinear initial parametrizations.

D is not a plane: its constant normal would have to annihilate
Delta=(1,w,z) for all independent w,z. Its second fundamental
quadratic is therefore nonzero at the generic point, as in KMFI32.
Every horizontal ruling direction is an asymptotic direction with
first coordinate0. Each actual Delta has first coordinate1, so it
is the other root of that quadratic. Over a general point of D
there is at most one such root. Distinct preimages under X give
distinct Delta directions. Thus **delta=1**.

Y1 is nonconstant on D, since Delta is tangent and its first
coordinate is1. Let E be the relative algebraic closure of k(X1)
in k(D)=k(w,z). It is a finite extension of k(X1).
Normalize A1_(X1) in E. Its coordinate ring embeds in k[w,z]:
an integral element is integral over that normal polynomial ring.
A general source affine line maps nonconstantly to this normalized
curve. Its projective completion has genus0, and every omitted
point must pull back to the one point at infinity of that line.
The normalized affine curve is therefore A1. Consequently
\[
                     E=k(t),\qquad X1=f(t),\quad
                  t\in k[w,z],\ f\in k[t]\setminus k .       \tag{2}
\]
This is the curve-normalization argument also used in KMFI38.

The geometric generic fibers of D->A1_(Y1) are lines.
After a constant linear change of Y2,Y3, write their equations as
\[
             X2=v,\qquad X3=r(t)v+s(t),\qquad r,s\in k(t).
                                                               \tag{3}
\]
The field of a line's coefficients is contained in E; taking the
relative algebraic closure accounts for all geometric components.
Thus k(D)=k(t,v). The change in Y2,Y3 corresponds to an invertible
affine change of the two source directions, so preserves (1).

Use primes for t derivatives and set
\[
 A=f'r''-f''r',\qquad B=f's''-f''s'.
\]
A normal to (f(t),v,r(t)v+s(t)) is
\[
                n=(-(r'v+s')/f',-r,1).
\]
Tangency and the asymptotic condition for Delta give, respectively,
\[
 z=rw+(r'v+s')/f',\qquad
                 Av+B+2r'(f')^2w=0.                       \tag{4}
\]
The identities follow by differentiating the surface: in coordinates
(t,v), the Delta direction is (1/f',w), the vv second derivative
vanishes, and the mixed second derivative has normal component r'.

If r'=0, (4) makes B=0, hence s'/f' is constant and D is a plane,
a contradiction. Thus r'!=0.

If A=0, the second equation of (4) puts w in k(t).
Then t is algebraic over k(w) and lies in k(w,z), whose relative
algebraic closure of k(w) is k(w) itself. Hence t in k[w].
Since w in k(t), k(t)=k(w), and t is affine linear in w.
The normal slope r(t) depends only on w.

If A!=0, eliminate v in (4):
\[
 z=p(t)w+q(t),\quad
 p=r-2(r')^2 f'/A,\quad
 q=s'/f'-r'B/(f'A).                                        \tag{5}
\]
Since t is polynomial in w,z, (5) gives
k(w,z)=k(w,t). Over k(w), the rational function wp(t)+q(t) must
therefore have degree1 in t.
Write p=P/R,q=Q/R with coprime P,Q,R in k[t].
There is no cancellation between wP+Q and R over k(w), so all
three polynomials have degree at most1. Write their coefficients
with subscripts1,0. Solving (5) yields
\[
 t=\frac{R0\,z-P0\,w-Q0}{P1\,w+Q1-R1\,z}.                  \tag{6}
\]
Polynomiality of t forces the denominator to be a nonzero constant.
Otherwise a nonconstant linear polynomial would divide another
linear polynomial, making t constant. Thus t is affine linear
in w,z, again with r depending on one affine source coordinate.

In either case the normal slope rho=r(t) becomes, after the
permitted affine direction change, a rational function of w alone.
The slope transforms by a constant Mobius transformation; it does
not acquire dependence on z. The new normalization matrix remains
polynomial. KMFI33 then makes rho a polynomial h(w): the top
z coefficient of lambda has second derivative0, whereas the
square of any denominator of h must divide it. A constant h would
force rank DX<=1 by KMFI38, so h is nonconstant.
KMFI44 now proves finiteness for every polynomial lambda satisfying
the normalization wave equation. This proves KMFI100.

## 2. Compactify the entire family of lines

**KMFI101.** Every sweep (1), including nonfinite sweeps, is finite
over a target neighborhood of the generic point of D. In particular
every prime of its finite normalization above D is present in the
polynomial source, and its **full** inertia is
\[
                         3^\delta 1^{N-3\delta},
                    \qquad \delta=1\ \text{or}\ 2,          \tag{7}
\]
where N=[k(gamma,w,z):k(S)].
No extra ramified or unramified normalization branches are omitted
over the generic point of this particular surface.

To see all possible missing branches, associate to (w,z) the
projective line through
\[
                   (1,X1,X2,X3),\qquad(0,1,w,z)
\]
in P3_[Y0:Y1:Y2:Y3]. The map A2->Gr(1,3) is regular.
Take the graph closure in P2_[a:w:z] times Gr(1,3), where a=1 is
the original affine direction chart, and normalize it to a
projective surface C. It contains A2 as an open subset; the
boundary is a finite union of curves and points over a=0.

Pull back the universal projective line to C, obtaining a normal
projective threefold T with a proper evaluation map e:T->P3.
Over the original A2, its points over Y0!=0 are exactly (1):
the Y1 coordinate along each line is X1+gamma, so gamma is the
unique affine line parameter. Thus A3_(gamma,w,z) is an open of T.

For any boundary point, its line contains the limit direction
(0,0,w,z). This follows from the closed incidence relation between
the line and its marked point (0,a,w,z).
If that line has an affine part, its unique point at infinity has
Y1=0. The affine line is horizontal.
It follows that every surface component of e(T-A3) meeting the
affine target is covered by horizontal straight lines.
This conclusion retains all boundary curves of the graph
normalization, including curves contracted by its projection to P2.
It does not assume that X itself is finite or that the graph
boundary is a single divisor.

If D were contained in this boundary image, it would be one of
these horizontal-line surfaces. KMFI100 would make S finite.
But then the inclusion A3_X->e^-1(A3_Y) is proper, because A3_X
is proper over A3_Y and the latter target is separated over A3_Y.
It is also an open immersion. It is therefore closed, and, being
dense in the irreducible e^-1(A3_Y), it is the whole space.
There would be no affine boundary image at all, a contradiction.

Hence D is not contained in the closed image of the projective
boundary. There is a target open V containing its generic point
over which T_V=A3_X x_Y V and S is proper.
The generic D fiber is finite: gamma0 maps generically finitely
to D with degree delta, while every gamma!=0 point is etale.
A positive-dimensional fiber component cannot meet the etale
locus and would have to lie in a fiber of X on gamma0.
Shrinking V again makes S quasi-finite, hence finite there.

Equivalently, factor e by its finite normalization
([Stacks, Tag03H0](https://stacks.math.columbia.edu/tag/03H0)).
The proper birational map to that normalization is an isomorphism
at its codimension-one generic points
([Stacks, Tag0BFP](https://stacks.math.columbia.edu/tag/0BFP)).
Thus the compactification argument really controls every
normalization prime over D, not just a selected proper model.

On this finite neighborhood the only ramification divisor in the
source is gamma0. Its different exponent2 gives index3 in
characteristic0. Its geometric residue degree is delta, so it
contributes exactly delta three-cycles. Every other point over D
is unramified. This proves (7). KMFI32 supplies delta<=2.

For a concrete nonfinite control put kappa^2+kappa+1=0,
a=-kappa-2 and
\[
 X=(awz,(a+1)w^2z,-wz^2).
\]
It collapses w0 on gamma0. Its branch surface is
d(Y)=Y2Y3-Y1^3/(3a), with d(S)=-gamma^3/(3a):
the whole three-sheet ramification over D is present.
For e->0 and v0!=0, the escaping arc
\[
 w=e,\quad z=v0/e,\quad \gamma=v0+eB/v0
\]
has target limit ((a+1)v0,0,B). These horizontal lines sweep the
separate nonproperness plane Y2=0. Thus completeness over D
genuinely holds even when there are missing branches elsewhere.

The result does not assert that S is globally finite. Its
nonproperness locus can have other surface components, and its
finite normalization can have omitted primes over them.

## 3. Consequences for the group, without an affine-source shortcut

**KMFI102.** Neither exceptional degree11 PSL2(11) action is the
generic field action of any normalized polynomial tautological
sweep S as in (1), regardless of finiteness, multiplier chart,
primitive-pair degree or scalar degree. The statement survives
common independent stabilization and arbitrary birational
recharts of the same field pair.

Indeed every nonidentity order3 element in either exceptional
action has type 3^3*1^2, contradicting (7).
For rank DX<=1, KMFI38 instead gives a natural symmetric group.
Thus both exceptional actions are excluded throughout the whole
three-variable normalized polynomial tautological sweep class,
and consequently throughout its KMFI59 branch-centered affine4
completions. Higher primitive pairs and nonconstant cubic scalars
cannot reopen this route to exceptional11.

More generally every rank-two sweep group contains a permutation
of type 3*1^(N-3) or 3^2*1^(N-6).
A primitive action with neither type is excluded as a sweep field.

The natural projective actions PSL2(q) on q+1 points have a complete
answer within this class, also allowing rank<=1:
\[
                 \text{they occur exactly for }q=2,4.      \tag{8}
\]
For every prime power q>=8, an order3 element fixes at most two
projective points. It moves at least q-1>=7 points, hence at
least nine because all nontrivial cycles have length3.
This contradicts the support bound6 in (7).
The q5 action on six points has order3 type3^2 with no fixed sheet.
Its point stabilizer normally generates the simple group, so
KMFI59 would complete such a rank-two field to a Keller map.
The divisor in (7) persists in that completion, contradicting
KMFI1. The q3 and q7 actions are excluded by the same completion
criterion and the absolute central-cover obstructions KMFI66--67.
Rank<=1 is symmetric and cannot give those actions.

For q2, rank-one X=(z^2,0,2z^3/3) gives the generic cubic pencil,
with field group S3=PSL2(2) on3. For q4, h=w^2,lambda1 in KMFI44
gives natural A5=PSL2(4) on5. These are actual sweep-field witnesses;
the corresponding group actions already have ordinary-affine Keller
witnesses in KMFI11 and12. Equation (8) does not exclude these
abstract groups in different actions, or natural projective actions
under unrelated Keller constructions.

The same cycle test excludes every natural odd dihedral action of
degree>=5: an order3 rotation, when it exists, has no fixed points
and moves at least nine letters; otherwise the group has no
order3 element. F20 has no order3 element at all. Rank<=1 gives
symmetric actions instead. Thus these still-open Keller actions
also require a different architecture.

The original objective remains unchanged. This closes the entire
normalized rank-two tautological route to exceptional11; it does
not prove that every Keller map has such a presentation.
An unrelated polynomial affine-space realization or an absolute
action-level obstruction is still needed. The full faithful
transitive action spectrum is open.

## Exact controls and assurance

    .venv/bin/python research/jacobian-followup/monodromy-theory/verify_sweep_boundary.py

[Control helper](sweep_boundary.py) · [model](sweep_boundary_model.json) ·
[receipt](sweep_boundary_receipt.json) ·
[primary pins](sources/sweep-boundary-manifest.json).

The compactification and horizontal-ruling arguments are written
universal proofs, not conclusions from finite tests.
The checker verifies Plucker/incidence identities, the complete
horizontal-ruling differential/elimination identities, polynomial
inverse constraints, affine covariance, exact finite and nonfinite
controls, and both exceptional action cycle dictionaries.
Formal verification, independent whole-proof review and external
review are not claimed.
