# Conductor does not eliminate rerooting ambiguity

## Result

The answer is **no**.  Node pairing, the conductor morphism, and all
Fitting data of the finite discriminant normalization do not distinguish
rerootings.  This is not merely a collision in a coarse numerical invariant:
every rerooting gives an explicit isomorphism of the entire decorated finite
cover.

Let \(H\) be a normalized seed and let \(a\ne0\) be a simple root.  Put

\[
 \kappa_a=-\frac1{aH'(a)},\qquad G_a(w)=\kappa_aH(aw).
\]

Then \(G_a\) is normalized and

\[
 E_{G_a}(w;s,t)
 =\kappa_a E_H\left(aw;\frac{s}{\kappa_a a},
                         \frac{t}{\kappa_a}\right).       \tag{1}
\]

Consequently the finite incidence covers are isomorphic over the affine
target change

\[
 (s,t)\longmapsto
 \left(\frac{s}{\kappa_a a},\frac{t}{\kappa_a}\right),
 \qquad w\longmapsto aw.                                \tag{2}
\]

The discriminant-normalization square is the restriction of (2):

\[
\begin{array}{ccc}
\mathbb A^1_w&\xrightarrow{\ w\mapsto aw\ }&\mathbb A^1_r\\
\nu_{G_a}\downarrow&&\downarrow\nu_H\\
\mathbb A^2_{s,t}&\xrightarrow{\ (s,t)\mapsto
(s/(\kappa_a a),t/\kappa_a)\ }&\mathbb A^2_{S,T}.
\end{array}                                               \tag{3}
\]

Here

\[
 \nu_H(r)=(H'(r),rH'(r)-H(r)).
\]

Equation (3) identifies the two finite extensions

\[
 k[H',rH'-H]\subset k[r],
 \qquad
 k[G_a',wG_a'-G_a]\subset k[w].                         \tag{4}
\]

Everything in the proposed richer package is functorially constructed from
(4), so it is transported rather than separated:

* \(G_a''(w)=\kappa_a a^2H''(aw)\), with every Fitting multiplicity retained;
* the ordered equal-image equations transform by
  \[
   D_s^{G_a}(w,u)=\kappa_a a^2D_s^H(aw,au),\qquad
   D_t^{G_a}(w,u)=\kappa_a aD_t^H(aw,au);
  \]
  diagonal saturation and the unordered quotient therefore agree;
* under the ring isomorphism (4),
  \[
   \operatorname{Ann}_{A_{G_a}}(B_{G_a}/A_{G_a})
   \longleftrightarrow
   \operatorname{Ann}_{A_H}(B_H/A_H),
  \]
  so both the upstairs conductor ideal and the finite downstairs conductor
  quotient agree.

The same argument transports zero, infinity, and the second-boundary center.
An ordering of Fitting points also does not help if it is intrinsic: the
isomorphism transports that ordering.  A labeling anchored to untransported
coordinate values is an additional mark, not an invariant of the finite
cover.

## Smallest tests

### Generic quartic: two indistinguishable rerootings

Write the normalized quartic line as

\[
 H_q(w)=-\frac12w^2(w-1)\bigl(qw-q+w+1\bigr).
\]

Its second nonzero root is

\[
 a=\frac{q-1}{q+1},
\]

and rerooting there gives the particularly simple involution

\[
 G_a=H_{-q}.
\]

For the rational clean point \(q=2\),

\[
 H_2=-\frac12w^2(w-1)(3w-1),\qquad
 H_{-2}=\frac12w^2(w-1)(w-3).
\]

They are distinct normalized seeds.  The rerooting has
\(a=1/3\), \(\kappa=-27\), hence

\[
 w_H=\frac{w_G}{3},\qquad
 s_H=-\frac{s_G}{9},\qquad
 t_H=-\frac{t_G}{27}.                                  \tag{5}
\]

Both discriminant curves are ordinary.  Each has two simple Hessian points,
one ordinary node with its two paired normalization branches, and conductor
degree six.  Substitution (5) identifies the saturated ordered node ideal,
its unordered quotient, the full Fitting polynomial, the upstairs conductor
polynomial, the downstairs conductor ideal, and the finite conductor map.

### Generic quintic: three indistinguishable rerootings

Take the fully split rational seed

\[
 H(w)=-\frac12w^2(w-1)(w-2)(w-3).
\]

The three rerootings, at \(a=1,2,3\), are

\[
\begin{aligned}
 G_1(w)&=-\frac12w^2(w-1)(w-2)(w-3),\\
 G_2(w)&=w^2(w-1)(2w-1)(2w-3),\\
 G_3(w)&=-\frac12w^2(w-1)(3w-1)(3w-2).
\end{aligned}
\]

All three are ordinary and pairwise distinct.  Their cubic Hessians are
squarefree; the original Hessian discriminant is \(39960\).  Each decoration
has three Fitting points, three nodes, six ordered node branches, and
conductor degree twelve.  For \(a=2\) and \(a=3\), substitution
\[
 w_H=aw_G,\quad
 s_H=s_G/(\kappa_a a),\quad
 t_H=t_G/\kappa_a
\]
identifies the complete decorated conductor squares.

There is also a stronger automorphism example already recorded in
`INTRINSIC_SELECTOR_ATTACK.md`: for

\[
 H=(w^2-w^5)/3
\]

the group \(\mu_3\) acts on one normalized incidence cover, preserves the
Fitting, node, and conductor constructions, and cyclically permutes all
three primitive-root candidates.  Thus the failure is visible both between
distinct generic rerooted seeds and as a symmetry of a single special
decorated cover.

## First Hessian-collision divisor

The quartic Hessian discriminant is

\[
 \operatorname{disc}(H_q'')=12(q^2+2).                 \tag{6}
\]

Work over \(K=\mathbb Q(q)/(q^2+2)\).  Put

\[
 m=\frac{q}{2(q+1)}.
\]

On (6),

\[
 H_q''(w)=-6(q+1)(w-m)^2.                              \tag{7}
\]

An exact implicit-equation/adjunction computation gives:

* the saturated off-diagonal equal-image ideal is the unit ideal, so the
  ordinary node has merged into the higher cusp and no separate node pair
  remains;
* the full Fitting divisor is \(2[m]\);
* the conductor upstairs is
  \[
   \mathfrak c_HK[w]=(w-m)^6.                          \tag{8}
  \]

The two primitive roots remain simple on this divisor.  Rerooting at
\(a=(q-1)/(q+1)\) again sends \(H_q\) to \(H_{-q}\).  It maps the doubled
Fitting point and the length-six conductor scheme by \(w_H=aw_G\), and the
two conductor squares are isomorphic.  Thus keeping nonreduced Fitting and
conductor structure does not rescue the selector at the first collision.

## Consequence for F2

F2 cannot be strengthened by replacing the affine root sheet with node
pairing, conductor, and full Fitting data.  Those decorations are invariants
of the unmarked finite normalization and therefore remain constant on the
rerooting groupoid.

What does eliminate the ambiguity is the extra Zariski--Main incidence datum:
the distinguished affine open, equivalently the unique unramified root sheet
whose divisorial center meets the regular reconstruction locus.  Among the
rerooted candidates, (3) sends the affine sheet of \(G_a\) to the root-\(a\)
sheet of \(H\); only \(a=1\) is affine.  The affine mark is therefore not
merely convenient.  It records information absent from the decorated finite
cover proposed in the question.

## Reproduction

Run

```bash
.venv/bin/python scripts/verify_conductor_rerooting_ambiguity.py
```

The checker verifies the universal transport identities, both generic
low-degree examples, the upstairs and downstairs conductor maps, and the
nonordinary quartic collision calculation over
\(\mathbb Q(\sqrt{-2})\).
