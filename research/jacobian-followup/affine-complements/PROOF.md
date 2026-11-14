# Affine complement recognition and the quartic ruled surface gate

14 November 2026. Work over C. PF2RULED1 is a conditional application to the
quartic finite normalization, with the recognition theorem imported from the
corrected paper. The application narrows compatible ruled completions to one
surface and two boundary classes. Existence of such a completion for a
hypothetical quartic Keller map remains unknown.

The subsequent [pole-reduction theorem](POLE_REDUCTION.md), PF2RULED2,
excludes every compatible quartic map of target pole height at most three.
The [complete pole-four atlas](POLE_FOUR_EXCLUSION.md), PF2RULED4, raises
that lower bound to five. The numerical m=2 cluster below is retained as
a necessary-condition calculation now excluded by the stronger argument.

## The corrected recognition theorem

Use Remy van Dobben de Bruyn,
[Divisors in projective bundles over the projective line whose complement is affine space, v2](https://arxiv.org/html/2608.27341v2),
14 September 2026, Theorem 3.1 and the discussion preceding it.

For a rank n bundle V over P1, use the quotient convention
X=P(V), f=pi*O(1), h=O_X(1). For distinct irreducible boundary divisors
D_i~h+d_i f, the complement is A^n precisely when there are two divisors,
their d_i differ by one, and their fibre hyperplanes coincide over exactly
one geometric point. Version 1 required torsion length one; v2 requires only
one point in its support. Higher contact is allowed.

The proof reduces by linear projection to a ruled surface, then changes
its ruling by boundary blowups and contractions. The initial projection
is an affine-bundle torsor; it becomes trivial over the recognized affine
plane. Lemmas 2.2, 2.4 and Corollary 2.7 supply the bundle descriptions.
The appendix applies the criterion to the cubic marked-root construction;
it supplies no plane quartic exclusion.

The next arguments are deductions for our retained quartic problem.

## What a compatible completion must preserve

Write pi:N -> A2 for the canonical affine finite normalization and
U=A2_source inside N. QPK4 leaves either one omitted prime R of transverse
index two, or R and a second unramified omitted prime. QCL4 gives the
Weil canonical/different class [K_N]=[R].

A **compatible ruled completion** here means an open embedding N -> S=P(V)
with V of rank two over P1, such that every component of S\U is an
irreducible divisor of relative degree one. In particular the whole affine
normalization N, its finite map, and its omitted primes are preserved.
Since S is smooth, this hypothesis also requires N smooth. Merely finding
a ruled surface birational to N does not establish it.

The complement of an affine open in a smooth surface has pure codimension
one. One way to see this is to localize near a hypothetical isolated
boundary point: normal Hartogs extension leaves the ring of functions
unchanged, whereas an affine open immersion cannot remove only that
codimension-two point. Thus S\N is a nonempty union of curves. It is nonempty
also because a positive-dimensional projective variety is not affine.

The localization sequence, using U^*/C^*=0 and Pic(U)=0, identifies the free
group on the components of S\U with Pic(S)=Z^2. There are exactly two
components. The closures of different omitted primes are different curves,
and at least one further curve lies at infinity.

Consequently the two-omitted-prime packet has no compatible ruled
completion. This excludes this form of completion, without excluding that
quartic packet. In the one-prime packet write

\[
 S\setminus N=B,\qquad \overline R=E,\qquad
 S\setminus U=B\cup E.
\]

Both B and E are sections. The recognition theorem then gives adjacent
classes and one intersection point, with arbitrary positive intersection
multiplicity. Because the completion preserves N, the notation R denotes
E intersected with N, not the target image of R.

## The different class forces a unique ruled surface

Normalize S as the Hirzebruch surface

\[
 \mathbb F_a=P(\mathcal O\oplus\mathcal O(a)),\quad a\geq0,\qquad
 f^2=0,\quad f h=1,\quad h^2=a.
\]

Set B=h+b f and E=h+e f. The canonical formula is

\[
 K_S=-2h+(a-2)f.
\]

Since N=S\B is smooth, Pic(N)=Pic(S)/Z[B]=Z[f]. Put

\[
 k=B^2=a+2b,\qquad \delta=e-b\in\{-1,1\}.
\]

Restriction to N gives

\[
 [K_N]=(k-2)[f],\qquad [R]=\delta[f].
\]

Simple ramification along the sole omitted prime means K_N=R as a divisor
class: the rational two-form pi*(du wedge dv) has exactly R as its zero
divisor on N. Therefore k-2=delta. If delta=-1, then k=1 and

\[
 B E=k+\delta=0,
\]

contradicting the required nonempty intersection. The only possibility is

\[
 \boxed{\delta=1,\quad B^2=3,\quad E^2=5,\quad B E=4.}
\]

All four intersections are supported at the one coincidence point. They
need not be transverse; imposing reduced intersection would discard this
entire possible application.

For a>0 let C0=h-a f be the negative section. An irreducible section
different from C0 has b>=0, since B.C0=b. The case B=C0 has B^2<=0 and is
already excluded. The equation a+2b=3 thus gives only (a,b)=(1,1),(3,0);
a=0 would require nonintegral b. On F3 the choice b=0 leaves C0 disjoint
from B, hence leaves a complete P1 inside the affine N. This is impossible.
We obtain the conditional reduction

\[
 \boxed{S=\mathbb F_1,\quad B=h+f=C_0+2f,\quad
 E=h+2f=C_0+3f,\quad B\cap E\text{ has length }4
 \text{ at one point}.}
\]

Here B is ample: B.f=1 and B.C0=1, so its complement is affine. This proves
PF2RULED1's surface statement for every compatible completion, without
claiming every quartic packet admits one.

## A positive control with fourth order contact

The surviving boundary geometry exists. With base coordinates [s:t] on P1,
take the line subbundle injections into O direct-sum O(1)

\[
 \mathcal O(-1)\longrightarrow\mathcal O\oplus\mathcal O(1):
       (s,t^2),\qquad
 \mathcal O(-2)\longrightarrow\mathcal O\oplus\mathcal O(1):
       (st,s^3+t^3).
\]

The two entries in either column have no common projective zero. Their
quotient bundles define smooth irreducible sections B and E in the
displayed classes. The determinant is

\[
 s(s^3+t^3)-t^2(st)=s^4.
\]

Thus their coincidence scheme is the length-four point s=0. The corrected
criterion gives S\(B union E)=A2. Also N=S\B is affine and
[K_N]=[E intersect N]. This control realizes the entire ruled boundary and
canonical-class package.

The subsequent [period obstruction](VOLUME_PERIOD_OBSTRUCTION.md) rules
out the required exact volume form on this entire marked-contact family,
including this positive control, at every pole height. Thus this control
has no finite quartic realization with the prescribed ramification. Its
abstract A2 complement and canonical-class package remain valid. The
[zero-period pair in the other marked position](WRIGHT_MARKED_MODEL.md)
remains open.

The blowdown F1 -> P2 gives an alternative small description of the gate.
If L is the pullback of a line and C0 the exceptional curve, h=L, f=L-C0.
Then B=2L-C0 is the strict transform of a smooth conic through the blowdown
point, and E=3L-2C0 is the smooth strict transform of an irreducible cubic
with a double point there. The conic and cubic have their remaining four
intersections concentrated at one point on F1.

## A finite quartic map must have base points at infinity

Assume a finite degree-four pi:N -> A2 exists on this compatible model,
with ramification divisor R. Regard its two coordinates P,Q as rational
functions on S, regular on N. Their nonzero pole orders m_P,m_Q along B
are positive integers: absence of a pole would make that function regular
on the projective S, hence constant.

The pulled-back two-form has divisor

\[
 \operatorname{div}_S(dP\wedge dQ)=E-3B.
\]

Indeed its only affine divisor is E and
K_S=E-3B=-2h-f. In local coordinates z=0 for B and w along B, write
P=z^{-m_P}(p(w)+O(z)), Q=z^{-m_Q}(q(w)+O(z)). The leading coefficient of
the wedge at order -(m_P+m_Q+1) is

\[
 m_Q p'(w)q(w)-m_P p(w)q'(w).
\]

When this coefficient is nonzero, the pole order three forces
m_P+m_Q=2, so both orders are one. But two members of |B| intersect with
total multiplicity B^2=3. They cannot contain the four affine intersection
points of a generic finite quartic fibre. Thus any such quartic map
requires both a pole order at least two and cancellation of its initial
Jacobian term.

Moreover the rational map [P:Q:1]:S -> P2 contracts B generically to a
point. If B instead dominated the target line at infinity, the two pole
orders would agree, say m, and the leading ratio p/q would be nonconstant.
The displayed coefficient would then be nonzero in characteristic zero,
forcing m=1 and the same degree-at-most-three contradiction.

In particular the rational map cannot extend to a morphism S -> P2.
If it did, the pullback of the target line at infinity would be mB and
the projection formula would give

\[
 4=(mB)^2=3m^2,
\]

which has no integer solution. This conclusion does not assume that the
extension would be finite; a generically finite proper extension suffices.

For an actual resolution of the rational map, write r_i for all proper
and infinitely near base-point multiplicities of the three-dimensional
linear system generated by P,Q,1 in |mB|, m=max(m_P,m_Q). Every center
lies over B. The resolved pullback of a line has square four, so

\[
 \boxed{m\geq2,\qquad \sum_i r_i^2=3m^2-4.}
\]

Since B is contracted, its intersection with that resolved linear system
is zero. B remains smooth under successive blowups, giving the additional
identity

\[
 \boxed{\sum_{i:\,\text{center on the successive strict transform of }B}
 r_i=3m.}
\]

Both sums include infinitely near centers; counting only proper points
would be invalid. At the first possible pole bound m=2, the arithmetic has
only two patterns: eight simple centers, six on the successive strict
transform of B; or one double center and four simple centers, all on that
strict transform.

The second pattern is actually impossible. At any proper base point x,
a multiplicity r=m would make the new exceptional coefficient in the
pullback of target infinity equal to m-r=0. Its generic point would map
into A2. However the pulled-back form with divisor E-3B has order
epsilon-2 there, where epsilon is 1 if x lies on E and 0 otherwise.
Both orders are negative, contradicting regularity of the pullback of
du wedge dv over A2. Thus every proper base point has r<m.

At m=2 the lone double center in the second pattern must be proper:
proximity requires a double predecessor for every infinitely near double
center, and no second double is available. Hence only eight simple
centers, six on successive strict B, remain at this pole bound. This
pattern still requires proximity, actual sections, the ramification ledger,
degree and finiteness on N.

## A computable coincidence test

For a specified splitting V=sum O(a_j), a divisor h+d f is represented by
a column v whose j-th entry is a homogeneous form of degree a_j+d.
First check each column is nowhere zero on P1; otherwise it is not an
irreducible relative hyperplane divisor. For two such columns v,w, form
all minors v_i*w_j-v_j*w_i. Their common homogeneous divisor g cuts out
the coincidence scheme W. All-zero minors mean the divisors are identical
and must be rejected.

Over C the number of support points of W is

\[
 r=\deg\left(g/\gcd(g,\partial_sg,\partial_tg)\right).
\]

Thus r=1 allows g to be any positive power of one homogeneous linear form.
Do this homogeneously: dehomogenizing can lose the point at infinity.
For g constant, r=0. This test also works with the differing entry degrees
of a nontrivial split bundle.

The checker includes the fourth-order F1 control, a length-three
coincidence on P1 times P1, a rank-three length-two control, two distinct
coincidence points, no coincidences, coincident columns, a nonsaturated
column, and nonadjacent classes with one coincidence point. It checks the
symbolic intersection/canonical identities and the pole/base-point
arithmetic. It does not certify a recognition theorem or resolve a
hypothetical finite map.

## Application to the retained programme

| Retained object | Consequence of this audit |
| --- | --- |
| QPK4 one omitted prime | A compatible ruled completion must be the F1 pair above. A singular N cannot embed in a smooth ruled completion. The finite map must contract B generically and have infinity base points satisfying both displayed budgets. |
| QPK4 two omitted primes, either target grouping | At least three total boundary curves would be needed. No compatible completion in this theorem's class exists. |
| QLD4 completed cusp and 2+2 charts | Local conductor compatibility remains intact. Local contact of ramification and companion branches is different from the global B,E coincidence; the number four above does not bound those local contacts. |
| QBX4 threefold W with W_s=A2 times Gm | This open has a nonconstant unit s and is not A3. Neither W nor its normalization comparison has been identified with a rank-three projective bundle over P1 with the specified boundary. |
| QGM4 and DZE4 global section and endpoint programme | A boundary-only birational model must prove that all exceptional curves and all omitted valuations are tracked. Recognition of an abstract A2 complement does not provide its finite normalization map, companion pairing or degree bound. |
| Pointed Keller pencil follow-up | Its generic curve can have positive genus. A section and a smooth pencil do not establish a P1-bundle; the recognition theorem supplies no missing rationality theorem. |

The immediate next gate is precise: either construct a compatible embedding
of the actual N into this F1 pair and then control the infinity base-point
cluster, or prove that the quartic finite-normalization valuations cannot
admit such an embedding. If N needs boundary contractions or a resolution
that changes N, record that as a different model and prove what it
preserves before applying recognition.

This settles the application audit, not the geometric-degree-four plane
Jacobian problem. The recognition theorem is an external input; the
surface, divisor, differential and base-point arguments here are written
deductions supported by exact symbolic controls, without formal
verification or independent external review.

