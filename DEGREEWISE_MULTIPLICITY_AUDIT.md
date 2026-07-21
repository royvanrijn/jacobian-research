# Degreewise stable multiplicity: five-lemma audit

This is the verification companion to the canonical standalone paper
[`Marked-Root Keller Maps and Degreewise Stable Multiplicity`](papers/marked-root-multiplicity/main.tex).
It splits that proof into independently checkable obligations; it is not a
second canonical theorem source.  Fix an integer `N>=4`, put `n=N-1`, and
work over `C`.  The audited conclusion is

\[
 \#\{\text{stable left--right classes of degree-}N\text{ Keller maps}
       \}\ge \tau(N-1).
\]

Here a *Keller map* is a polynomial map `A^3 -> A^3` with nonzero constant
Jacobian, and stable left--right equivalence allows polynomial automorphisms
of source and target after adjoining identity variables.  Each lemma below
has its own input list; in particular, the construction lemmas do not use the
boundary calculation, and the functoriality lemma does not use either family.

## Lemma 1 — the displayed weighted seed is admissible

For every `N>=4`, set

\[
 H_N(W)=W^2(1-W)(1+W^{N-3}).                              \tag{1}
\]

Then `H_N` has degree `N`, has a double zero at `0`, a simple zero at `1`,
and `N-3` further distinct simple zeros, none equal to `0` or `1`.  If
`c=-H_N'(1)`, then

\[
 c=2,\qquad \kappa=H_N''(1)/c=-(N+1)\ne-2.               \tag{2}
\]

Thus (1) satisfies the weighted-seed hypotheses.  Taking `b_0=1` and

\[
 a_0=-\frac{1+\kappa}{2+\kappa}
\]

in the weighted construction gives a polynomial map `G_N:A^3->A^3` with
constant nonzero Jacobian and generic fiber degree `N`.

### Check

The roots of `1+W^{N-3}` are simple in characteristic zero, are nonzero, and
cannot equal `1`.  Differentiating (1) at `1` gives (2).  For completeness,
the weighted construction is

\[
 v=xy,\quad S=x^2z,\quad u=1+v,\quad
 \gamma=1+a_0v+b_0S,\quad W=u\gamma,
\]

Put `p=H'`, `q=(Wp-H)/c` and use

\[
 C=x\gamma,\qquad B=\frac{c+p(W)/\gamma}{x},\qquad
 A=\frac{u+q(W)/\gamma^2}{x^2}.                            \tag{3}
\]

The conditions at `0` and `1`, together with the displayed value of `a_0`,
make both quotients in (3) polynomial; direct triangular differentiation
gives `det DG_N=b_0c=2`.  The inverse equation

\[
 H_N(W)-BCW+cAC^2=0                                      \tag{4}
\]

has degree `N` and reconstructs the source at every generic simple root, so
the generic degree is exactly `N`.

**Uses only:** elementary differentiation and the general weighted
divisibility/Jacobian identity.  It does not use Lemmas 2--5.

## Lemma 2 — every proper divisor gives a cancellation map

For every proper positive divisor `r|n`, put

\[
 m=n/r-1.                                                  \tag{5}
\]

Then `m,r>=1`.  There is a polynomial cancellation map
`F_{m,r}:A^3->A^3` with nonzero constant Jacobian, generic fiber degree

\[
 r(m+1)+1=N,                                               \tag{6}
\]

and an explicit fiber of `N` distinct points.

### Check

The cancellation map starts from

\[
 A=1+xy^m,\quad B=A^{r+1}z+y^{m+1}h(A),\quad
 P=AB,\quad Q=y+xB,
\]

\[
 R=C\int_0^{x/A}\{1-t(Q-Pt)^m\}^r\,dt.                  \tag{7}
\]

Its localized Jacobian is identically `-C`.  Polynomiality is equivalent to
a finite order-`r` cancellation condition.  Writing `q=h(0)`, its constant
condition is the monic degree-`mr` polynomial

\[
 M_{m,r}(q)=\sum_{j=0}^{mr}(-1)^j
 \binom{mr+r+1}{j}q^{mr-j}.                               \tag{8}
\]

It is useful to audit squarefreeness before invoking algebraic closure.  Up
to a nonzero scalar, (8) is

\[
 y(q)=\int_0^1u^r(1-q+qu)^{mr}\,du
     =\frac1{r+1}\,{}_2F_1(-mr,1;r+2;q).
\]

The hypergeometric differential equation and uniqueness show that a common
zero of `y` and `y'` away from `0,1` would force `y=0` identically.  At
`q=0,1` the integral is nonzero.  Hence (8) is squarefree.  It has a complex
root because `mr>=1`; at each root the nonzero derivative recursively and
uniquely determines `h_1,...,h_r` in
`h(A)=q+h_1A+...+h_rA^r`.  Formula (7) is then polynomial.

The inverse polynomial

\[
 \Psi(T)=C\int_0^T\{1-t(Q-Pt)^m\}^r\,dt-R               \tag{9}
\]

has degree (6), and its regular reconstruction proves exact generic degree.
At `(P,Q,R)=(1,0,0)` the standard beta-integral argument shows that (9) has
`N` simple roots, all reconstructing to distinct source points.

**Uses only:** the cancellation operator, the elementary ODE squarefreeness
test, and reconstruction.  It does not use Lemmas 1 or 3--5.

## Lemma 3 — the boundary vertices and labels are exhaustive

Let `U=A^3`, let `Y=A^3`, and for either map let

\[
 \bar X_F=\operatorname{Norm}_Y k(U),\qquad
 \partial_F=(\bar X_F\setminus U)_{\mathrm{red}}.         \tag{10}
\]

For the split weighted map, the complete target boundary list is the
discriminant `Delta_H` and `C=0`.  There is one boundary prime over
`Delta_H` with `(e,f)=(2,1)` and `N-3` geometric boundary primes over `C=0`,
each with `(e,f)=(1,1)`.

For the cancellation map, the complete target boundary list is the
discriminant `Delta_{m,r}` and `P=0`.  There is one boundary prime over the
discriminant with `(e,f)=(r+1,1)` and `mr-1` geometric boundary primes over
`P=0`, each with `(e,f)=(1,1)`.  Since `N>=4`, `mr-1>0`, so both listed target
vertices actually receive boundary.

### Check

Regular reconstruction off the two displayed divisors confines every
boundary image to them.  Completeness then follows at the generic points
from the DVR identity

\[
 [k(U):k(Y)]=\sum_{E\mid Z}e(E/Z)f(E/Z).                  \tag{11}
\]

For the weighted map the two sums, including affine primes, are

\[
 2+(N-2)=N,\qquad 2+1+(N-3)=N.                            \tag{12}
\]

For the cancellation map they are

\[
 (r+1)+N-(r+1)=N,\qquad (r+1)+1+(mr-1)=N.                \tag{13}
\]

Every omitted term in (11) would contribute positively, so (12)--(13) rule
out an unseen normalization prime.  The ramified prime cannot lie in the
affine source because the Keller map is etale there.  The canonical finite
normalization is affine and normal; its affine open `U` has pure divisorial
complement, and a height-one prime in a finite integral cover lies over a
height-one target prime.  Thus no higher-codimension boundary image is being
silently discarded.

The discriminant vertex is intrinsically marked as the unique vertex
receiving ramification.  The second vertex receives only unramified boundary
primes, so the two vertices cannot be swapped.

**Uses only:** the two reconstruction charts, normality, and (11).  It does
not use the intersection calculation or stable invariance.

## Lemma 4 — `(e_Delta,mu)` is a stable left--right invariant

For a map satisfying Lemma 3, let `Z_Delta` be the unique target boundary
vertex receiving a ramified boundary prime, and `Z_0` the other vertex.  Put

\[
 e_\Delta=e(E_\Delta/Z_\Delta),\qquad
 \mu=\min\{a\ge1:\operatorname{Nil}
 \mathcal O_{Z_\Delta\cap Z_0}^{,a}=0\}.                 \tag{14}
\]

For a reduced intersection the convention in (14) gives `mu=1`.  The pair
`(e_Delta,mu)` is preserved by polynomial left--right equivalence and by
adjoining identity variables.

### Check

A left--right equivalence identifies the two function-field extensions over
the corresponding targets.  Functoriality and uniqueness of normalization
therefore identify the finite covers, their distinguished affine opens,
their boundary primes, and the target scheme intersections.  Ramification
indices and nilradicals are intrinsic under this isomorphism.

After stabilization by `A^s`, the normalization is
`bar X_F x A^s`, and every intersection ring `R` becomes
`R[t_1,...,t_s]`.  Polynomial extension is faithfully flat,
`Nil(R[t])=Nil(R)[t]`, and the exact nilpotency index is unchanged.
Ramification indices are unchanged by the same regular base change.  Hence
(14) is stable.

**Uses only:** the canonical boundary object furnished by Lemma 3 and
standard functoriality of finite normalization.  It does not use the family
formulas of Lemmas 1, 2, or 5.

## Lemma 5 — weighted contact is reduced; cancellation contact is thick

For the split weighted map of Lemma 1, the two canonical target boundary
divisors meet in

\[
 C=0,\qquad B^2-4h_2cA=0,                                 \tag{15}
\]

where `H(W)=h_2W^2+O(W^3)` and `h_2c!=0`.  Its coordinate ring is
isomorphic to `C[B]`, so it is reduced and `mu=1`.

For every cancellation pair supplied by Lemma 2, the two canonical target
boundary divisors meet in a scheme with ring

\[
 \frac{\mathbb C[Q,R]}
 {\bigl(Q^{mr(m+1)}((r+1)RQ^m-C)\bigr)}
 \simeq
 \frac{\mathbb C[Q,R]}{(Q^{mr(m+1)})}
 \times\mathbb C[Q,Q^{-1}].                              \tag{16}
\]

Its nilradical has exact index

\[
 \mu=mr(m+1)=m(N-1)>1.                                   \tag{17}
\]

### Check

For (15), write `H(W)=W^2J(W)`, scale `W=Cu` near the double-zero chart,
and divide out the forced power of `C`.  The limiting quadratic is
`h_2u^2-Bu+cA`; its discriminant is (15), which is linear in `A`.

For (16), normalize the critical divisor of (9) with

\[
 Y=Q-PT,\qquad T=Y^{-m},\qquad P=(Q-Y)Y^m.
\]

Eliminating `Y` and setting `P=0` gives, up to a unit,
`Q^{mr(m+1)}((r+1)RQ^m-C)`.  The two factors are comaximal because `C!=0`,
giving the product decomposition.  Its first factor has nilradical generated
by `Q` with exact index (17).  The restriction `N>=4` excludes `(m,r)=(1,1)`,
the sole case in which `P=0` is only a coordinate trace rather than a
boundary vertex.

**Uses only:** the local equations of the two normalized critical divisors.
It does not use the stable-invariance argument.

## Assembly check

There are `tau(n)-1` proper positive divisors of `n=N-1`.  Lemma 2 supplies
one cancellation map for each.  Distinct divisors `r` have distinct
`e_Delta=r+1` (and also distinct values
`mu=n(n/r-1)`), so Lemmas 3--5 make these maps pairwise stably inequivalent.
Lemma 1 supplies one weighted map; Lemma 5 separates it from every
cancellation map by reducedness.  Therefore there are at least

\[
 (\tau(N-1)-1)+1=\boxed{\tau(N-1)}
\]

stable left--right classes in generic degree `N`.  Cite the canonical paper,
not this audit, for the theorem statement and proof.

## Reproduction hooks

The abstract degree sums and family charts are independently checked by

```bash
python3 scripts/audit_boundary_exhaustion_independent.py
```

The all-parameter thick-intersection formula has a bounded symbolic
regression in

```bash
.venv/bin/python scripts/verify_scheme_boundary_all_parameters.py
```

These scripts are regressions, not substitutes for Lemmas 3 and 5.  The
standalone paper is in
[`papers/marked-root-multiplicity/`](papers/marked-root-multiplicity/).
