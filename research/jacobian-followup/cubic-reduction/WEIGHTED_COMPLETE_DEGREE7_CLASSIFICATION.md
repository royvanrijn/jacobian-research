# Complete degree-seven Keller class of weights (1,-1,-2)

**Theorem (KCUB15).** Over a characteristic-zero field, a polynomial Keller
map in three variables, of ordinary degree at most seven, with source and
output weights (1,-1,-2), is either a polynomial automorphism with linear
first component, or a diagonal input/output transform of the established
13-term representative. Thus every counterexample in this entire weighted
degree range has **exactly 13 nonlinear monomial occurrences**. The minimum
ordinary degree of a counterexample in these weights is **seven**, attained
by the known map. No collision or coefficient-nonvanishing assumption is made.

The 13-term sparsity minimum has the explicit degree-at-most-seven cutoff.
The degree minimum does not need an upper cutoff, since all smaller degrees
are covered. Neither assertion concerns other weights or arbitrary polynomial
coordinates, and neither closes the global simplicity goal.

## Complete coefficient cases

A weighted map has no constant terms and has diagonal linear part. A nonzero
constant determinant makes each diagonal coefficient nonzero. Dividing each
output by its linear coefficient preserves degree, nonlinear support and
injectivity, and gives identity linear part and determinant one.

Put u=ab,w=a^2c and F=(af,a^{-1}g,a^{-2}h). All compatible source atoms
through degree seven give

\[
f=1+\mu u+\nu w+R u^3+P u^2+Q u w+T w^2. \tag{1}
\]

There are six possible nonlinear first-output atoms, eight in the second
output and nine in the third. The determinant equation is

\[
J=g(f_u h_w-f_w h_u)+2h(g_u f_w-g_w f_u)
       +f(g_u h_w-g_w h_u)=1. \tag{2}
\]

[KCUB14](WEIGHTED_UNMARKED_DEGREE7_CLASSIFICATION.md) classifies the case
mu*nu!=0: diagonal conjugacy gives the unique normalized 13-term map.
Here we close all coefficient-zero boundaries. Diagonal scalings normalize
the three remaining cases to

\[
(\mu,\nu)=(0,-1),\quad(-3,0),\quad(0,0).
\]

Each model retains all 21 other coefficient parameters, with no collision or
sparsity equations. Zero designated atoms are recorded rather than silently
removed from the common 23-slot source enumeration.

The smooth-cubic boundary lemma from
[KCUB13](WEIGHTED_FULL_DEGREE7_CLASSIFICATION.md) gives RT=0 in all three cases.
Indeed, if RT!=0, completing the square with y=2Tw+Qu+nu gives

\[
y^2=(Qu+\nu)^2-4T(Ru^3+Pu^2+\mu u+1).
\]

The right side has degree three. A singular boundary contradicts (2);
a smooth one has the same single-infinity-point differential and pole-order
contradiction as before. Neither argument depends on mu,nu or on a collision.
The general completed-square identity is independently replayed.

## A univariate boundary obstruction

Two rotations of the source use the following elementary argument. Suppose

\[
g=A(L)t+B(L),\qquad h=C(L)t^2+D(L)t+E(L),
\]

with deg A<=2,deg C<=1 and B(0)=E(0)=0, and f=phi(L) is nonconstant with
phi(0)=1. Let epsilon be either sign, according to the coordinate rotation.
Equation (2) becomes

\[
J=\epsilon\{\phi'(g h_t-2h g_t)+\phi(g_Lh_t-g_th_L)\}=1. \tag{3}
\]

If A=0 identically, then J=epsilon*(phi B)'*h_t. The product could be one
only if (phi B)' were a nonzero constant. Since B(0)=0, nonzero B and
nonconstant phi make deg(phi B)>=2, impossible.

Otherwise the t^2 coefficient gives AC'-2A'C=0, so C=cA^2 for a constant c.
The degree bounds allow either C=0 or constant A,C. Replace h by h-cg^2;
this preserves J, B(0)=E(0)=0 and D(0), and removes C. No later degree bound
on the modified h is needed. The t coefficient now gives

\[
\phi D=kA
\]

for a constant k. If k!=0, the constant coefficient is

\[
J=\epsilon D\,\bigl(\phi B-\phi^2E/k\bigr)'=1.
\]

Hence D is a nonzero constant. Integrating and evaluating at L=0 makes a
nonconstant phi divide L, impossible since phi(0)=1.

If k=0, then D=0 and (3) reduces to
J=-epsilon*A*(2phi'E+phi E')=1. This forces A to be constant. The remaining
polynomial equation is impossible when deg phi>=2 by its leading degree.
When deg phi=1 it forces E constant, contrary to E(0)=0. Thus (3) admits
no such completion, provided either D(0)!=0 (which already rules out k=0)
or D(0)=0 and the source normalization makes A(0) nonzero.

For a squarefree boundary with at least two roots, the degree restrictions
needed above follow directly from source supports:

- For f=phi(u), take L=u,t=w. The coefficients of above-threshold t powers
  in g,h have degree less than two in L. Their vanishing at every root of phi
  makes them zero identically. Three complete linear equations give the
  displayed forms and D(0)=1. Consequently every nonconstant phi(u) of
  degree at least two is excluded.
- For f=phi(sigma*u+w) of degree two in its argument, take t=u,L=sigma*u+w.
  The same coefficient-degree observation gives seven complete linear
  equations. Here D(0)=-sigma and A(0)=1-sigma*g8. If sigma!=0 then
  D(0)!=0. If sigma=0 then A(0)=1. Thus this entire parallel-line branch
  is excluded as well.

A multiple root would make f=f_u=f_w=0 and already contradict (2).
The independent audit reconstructs both complete linear systems and their
source bounds; the degree/divisibility argument is written mathematics,
not proof-assistant formalization.

## The mu=0,nu=-1 branch has no Keller completion

Nine of the new identities below are cleared generic-field identities;
every target has explicit rational polynomial multipliers. A nonzero scalar
factor in a displayed target does not affect its vanishing consequence.

| Branch | Necessary target or treatment |
|---|---|
| T=0 | Q^16 forces Q=0 |
| R=0 | T^2(PT-Q^2/4)^7 forces T=0 or 4PT=Q^2 |
| T=Q=0, polynomial-boundary constraints | P R^3 forces P=0 or R=0 |
| T=Q=P=0,R!=0 | Reduced full-ring unit with 1-inv*R |
| T=Q=R=0 | P^4 forces P=0 |
| R=P=Q=T=0 | Exact unit for f=1-w, all 17 tail atoms allowed |
| R=0,T!=0,4PT=Q^2 | Write P=sigma^2*T,Q=2sigma*T |
| Previous row,sigma!=0 | sigma^5*T^9 forces a contradiction |
| Previous row,sigma=0 | f=1-w+T*w^2: the parallel-line obstruction |

The polynomial charts used for the cubic/quadratic rows are u=t,
w=1+Pt^2+Rt^3. Their boundary Wronskian has nonzero constant right side,
forcing deg G<=1,deg H<=2. For the conic chart with sigma!=0 use
u=(-1+t-Tt^2)/sigma,w=t-sigma*u. All denominator branches are displayed.

## The nu=0,mu=-3 branch has no Keller completion

The T=0 graph identity is Q^5. For R=0, the conic identity is
(PT-Q^2/4)^7. Thus the remaining cases are:

1. T=Q=0: f=1-3u+Pu^2+Ru^3. Any nonzero P or R gives the univariate
   obstruction above. The only remaining quotient is f=1-3u.
2. R=0,T!=0,4PT=Q^2: set P=sigma^2*T,Q=2sigma*T. The coefficient of u
   after using L=sigma*u+w is always -3, so there is no parallel exception.
   The polynomial chart u=(1+Tt^2)/3,w=t-sigma*u supplies 12 boundary
   linear equations. A reduced full-ring unit, with saturation 1-inv*T^2,
   excludes this entire branch, including roots of the extra factor in
   the unused generic certificate. Its only linear-solution denominator
   is the already nonzero T.

For f=1-3u, use v=1-3u and x=w. Write G(v,x),H(v,x) for the tails. The
boundary Wronskian forces G2=K*v and H3=0, where subscripts denote powers
of x. The source supports give deg G1<=2,deg H1<=3,deg H2<=1 in v, and
H1(1)=1. Equation (2) is -3 times

\[
\mathcal E=(G+vG_v)H_x-G_x(2H+vH_v).
\]

Its two highest coefficients are

\[
[x^3]\mathcal E=-2Kv^2 H_2',\qquad
[x^2]\mathcal E=v(2G_1'H_2-G_1H_2'-2KH_1-2KvH_1').
\]

If K=0, g is globally linear in w and h globally quadratic in w; the
univariate obstruction applies even though phi has degree one.
If K!=0, the first identity makes H2 a constant C. If C=0, the second
identity would make (vH1)'=0, contradicting H1(1)=1. Thus C!=0. It then
forces deg H1<=1. These are complete necessary linear constraints on
this branch. An exact unit with saturation 1-inv*K*C, containing only
13 multiplier monomials, excludes all its coefficients. No numerical K
or C value is fixed in this proof. A normalized auxiliary calculation
is retained but unnecessary.

The original full 17-tail solve timed out at 20 seconds. This structural
split and exact unit, rather than that timeout, settle its mathematical class.

## The mu=nu=0 branch is an automorphism sector

The graph and conic cleared targets are Q^5 and (PT-Q^2/4)^3 respectively.
After RT=0, the graph remainder is f=1+Pu^2+Ru^3, excluded by the univariate
obstruction unless P=R=0. The conic remainder with T!=0 is a quadratic
phi(sigma*u+w), excluded by the parallel-line obstruction. Thus f=1.

Now F1=a and (2) is simply det D(g,h)=1. The quotient polynomials g,h
have ordinary degree at most four. We explicitly import the plane result
from [Nakai-Baba, Theorem3(1)-(2), printed p408](https://ir.library.osaka-u.ac.jp/repo/ouka/all/8846/ojm14_02_08.pdf):
the prime-degree and degree-four cases, together with the elementary
linear case, cover every possible component degree here. Hence (g,h)
has a polynomial inverse (U,V). This import is not re-proved by our audit.
The source PDF, theorem/page identification, checksum and visual-read
receipt are pinned in [the import record](weighted-zero7/sources/plane-quartic-import.json).

There is also a direct polynomial inverse for the full source map:

\[
a=A,\quad b=A^{-1}U(AB,A^2C),\quad c=A^{-2}V(AB,A^2C).
\]

The expressions are polynomial. Both inverse polynomials have zero constant
term. Since g=u+lambda*w+higher terms and h=w+higher terms, the linear part
of V is the second target variable, with no first-variable term. Every U
monomial therefore has weight at least one under weights(1,2), and every V
monomial has weight at least two. The displayed powers of A never become
negative. The inverse identities hold for A!=0 and extend as polynomial
identities to A=0. Characteristic-zero fields follow by base extension and
uniqueness of the formal inverse, which descends its coefficients.

## Conclusion and exact audit

The zero-coefficient branches supply only this automorphism sector.
Combining them with KCUB14 proves the theorem and the matching weighted
bounds. The known representative has degree seven, 13 nonlinear terms
and a verified collision. No lower-degree weighted counterexample exists,
and no other weighted degree-at-most-seven counterexample has a different
normalized coefficient vector.

Run

```sh
.venv/bin/python research/jacobian-followup/cubic-reduction/verify_weighted_zero7.py
```

[Validation](weighted-zero7/validation.json) reconstructs the three full
universal source-variable determinants, every atom and boundary linear
system, all allowed-monomial row reductions, the high-w identities and
source degree restrictions, the general completed square, and all thirteen
new certificates. The latter contain 12,342 multiplier monomials, with
maximum multiplier degree 15: nine cleared parameter identities and four
units. It also checks the degree-four plane reduction and seals the
attributed primary theorem import. The written curve and univariate
arguments and the imported plane theorem are not formal machine proofs.

The generic nu-zero parabola certificate has an extraneous denominator
factor and is not used; its reduced full-ring unit covers every T!=0.
The direct nu-zero linear-tail timeout is retained. Initial replay logs,
including a rejected sign in an abstract coefficient identity, remain
stored; the corrected independent replay passes. No timeout is relabeled
UNSAT, and no generic result is used without its exceptional branches.

The next raw-sparsity route must leave these weights, exceed degree seven,
or change polynomial coordinates into a different class. For example,
weights(1,-1,-3) change the determinant reduction and exclude the first
w^2 atom by the degree cutoff; prior work must be checked before a new
bounded census. Ordinary cubic12/36, cubic-linear89 and normalized3D[10,48]
remain separate, unchanged construction bounds. Their global minima and
the full simplicity goal remain open.

The subsequent [KCUB16 result](WEIGHTED_123_DEGREE7_EXCLUSION.md) excludes
the entire degree-at-most-seven weight(1,-2,-3) class. The first bounded
(1,-1,-3) and(2,-1,-3) searches remain UNKNOWN, with source models and
all timeout checkpoints retained; unchanged caps are not the next step.
