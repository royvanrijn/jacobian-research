# Plane boundary exclusion: residue immersion and the conductor gap

## Result

This note isolates a plane theorem suggested by the repository's
three-dimensional boundary-cancellation constructions.  It gives a complete
exclusion for every puncture number once "saturated rank-one link with no
additional divisors" is expressed as an exact generic sheet budget.

Work over an algebraically closed field \(k\) of characteristic zero.  Let

\[
 F:\mathbb A^2\longrightarrow\mathbb A^2
\]

be a Keller map, let

\[
 \pi:\overline X=\operatorname{Norm}_{\mathbb A^2} k(x,y)
 \longrightarrow\mathbb A^2
\]

be its finite Zariski--Main normalization, and identify \(\mathbb A^2\) with
its distinguished open subset \(U\subset\overline X\).

Assume:

1. the reduced boundary \(\overline X\setminus U\) has one prime divisor
   \(D\);
2. \(\overline X\) is regular along \(D\), and \(D\) is smooth;
3. the smooth completion of \(D\) is \(\mathbb P^1\);
4. no residual ramification divisor occurs in addition to the tame
   transverse ramification of \(D\);
5. over the generic point of \(B=\pi(D)\), the exhaustive divisor
   decomposition consists of the boundary prime \(D\), with ramification and
   residue degrees \((e,1)\), and one primitive affine linked prime, with
   degrees \((1,1)\).

Condition 5 is the precise strengthened meaning of "saturated rank-one
valuation link, with no additional divisorial components" used in the
theorem.  It includes residue primitivity, not only primitivity of the normal
valuation or of the unit lattice.

Condition 4 is automatic from 1--2 and the Keller condition if the divisor
audit is exhaustive: any additional zero of the residual Jacobian would
give a curve meeting \(U\), where \(F\) is etale.

Put \(B=\pi(D)\).  Then:

> **Residue-immersion theorem.**
> The finite residue map \(D\to B\) is immersive everywhere.  If
> \(\widetilde B\) is the normalization of \(B\), the induced finite map
> \[
> g:D\longrightarrow\widetilde B
> \]
> is unramified everywhere.

The nonproper-value theorem for plane nonsingular polynomial maps says that
\(B\) has one place at infinity.  Thus
\(\widetilde B\simeq\mathbb A^1\).  Riemann--Hurwitz now gives:

> **Puncture exclusion.**
>
> If \(D\) has two or more punctures, it is impossible.  If
> \(D\simeq\mathbb A^1\), then \(g\) has degree one.  Hence
> \(D\to B\) is the normalization map and is immersive.

Consequently the only possible one-puncture survivor identifies two or more
points of \(D\simeq\mathbb A^1\) in the target.  Equivalently, \(B\) must
have nontrivial conductor gluing.  Let \(a\) be the total generic degree
contributed by affine primes over \(B\).  The generic DVR degree formula is

\[
 d=e+a.
\]

A conductor identification forces \(d\ge2e\), hence \(a\ge e\).  The case
\(e=1\) is impossible independently: conditions 1--4 then make \(\pi\) a
connected finite etale cover of \(\mathbb A^2\), so it has degree one and
cannot have a boundary.  Thus \(e\ge2\).

Condition 5 has \(a=1\).  It therefore contradicts \(a\ge e\).
Equivalently, its specialized degree formula is

\[
 d=e+1.
\]

Every conductor identification puts at least two distinct points
\(p,q\in D\) in one finite fiber.  Flatness makes that fiber have length
\(d\), while the local intersection multiplicity at each of \(p,q\) is at
least \(e\).  Hence

\[
 2e\le d=e+1.
\]

Since \(e\ge2\), this is impossible.

If no conductor identification occurs, \(B\simeq\mathbb A^1\), which is
impossible for a nonsingular plane polynomial map by Nguyen Van Chau's
line-component theorem.

We therefore obtain the following precise restricted theorem.

### Plane boundary-exclusion theorem

There is no nonproper plane Keller map satisfying conditions 1--5.

More sharply:

- conditions 1--4 already exclude every case with at least two punctures;
- in the one-puncture case they force a degree-one immersive normalization;
- the primitive minimal-sheet condition 5 excludes its only possible
  conductor gluing.

## 1. Why generic tameness is not enough

In characteristic zero, the extension of the generic DVR at \(D\) is
automatically tame.  Thus the phrase "tame generic ramification" alone adds
no restriction.  The useful hypothesis is **logarithmic purity along the
whole boundary curve**:

- \(\overline X\) is regular along \(D\);
- the transverse different has order \(e-1\) at every point of \(D\); and
- after removing that factor, the Jacobian is a unit.

An exhaustive "no additional divisors" condition on a regular graph model
implies this.  A condition checked only at the generic valuation does not.
Exceptional divisors over special boundary points are precisely where the
Riemann--Hurwitz deficit could otherwise hide.

## 2. Proof of residue immersion

Fix \(p\in D\).  Since \(\overline X\) is regular at \(p\), let \(z=0\) be a
local equation of \(D\), and choose a parameter \(t\) along \(D\).  Let \(e\)
be the transverse ramification index.  Tameness gives

\[
 \operatorname{ord}_D(\operatorname{Jac}\pi)=e-1.
\]

Write

\[
 \operatorname{Jac}\pi=z^{e-1}h.
\]

The function \(h\) is a unit at every point of \(D\).  Indeed, if \(h(p)=0\),
the zero set of \(h\) contains a curve through \(p\).  If that curve is
\(D\), the coefficient of \(D\) in the ramification divisor is larger than
\(e-1\).  Otherwise its generic point lies in \(U\), contradicting the
Keller condition.

The logarithmic differential at \(D\) separates into a nonzero normal map
and the tangential differential of

\[
 \varphi=\pi|_D:D\longrightarrow B.
\]

Its determinant is \(h|_D\).  Since \(h|_D\) is a unit and the normal factor
is nonzero, \(d\varphi\) is nonzero at \(p\).  Hence \(\varphi\) is immersive
everywhere.

Let \(\nu:\widetilde B\to B\) be the normalization.  Finiteness of
\(\varphi\) gives a unique factorization

\[
 D\xrightarrow{g}\widetilde B\xrightarrow{\nu}B.
\]

If \(dg\) vanished, then \(d\varphi=d\nu\circ dg\) would vanish.  Therefore
\(g\) is unramified.

## 3. Riemann--Hurwitz leaves no two-puncture cover

The curve \(B\) is an irreducible component of the nonproperness set.  Since
there is only one boundary prime, it is the whole nonproperness set.
The plane nonproper-value theorem gives one point at infinity, so
\(\widetilde B\simeq\mathbb A^1\).

Extend \(g\) to

\[
 \overline g:\mathbb P^1\longrightarrow\mathbb P^1
\]

of degree \(n\).

### One puncture

If \(D=\mathbb A^1\), the unique point at infinity is the unique pole of
\(\overline g\), of order \(n\).  It contributes \(n-1\) to the ramification
divisor.  There is no ramification in \(D\).  Riemann--Hurwitz requires total
ramification \(2n-2\), hence

\[
 n-1=2n-2,\qquad n=1.
\]

### Two punctures

If \(D=\mathbb G_m\), finiteness over \(\mathbb A^1\) forces both punctures
to lie over infinity.  Let their pole orders be \(a,b>0\), so \(a+b=n\).
Their total ramification contribution is

\[
 (a-1)+(b-1)=n-2.
\]

There is no ramification in \(D\), whereas Riemann--Hurwitz requires
\(2n-2\).  This would give \(n-2=2n-2\), hence \(n=0\), a contradiction.

This is the divisor-freedom obstruction sought in the plane: a
two-puncture residue cover necessarily spends \(n\) units of ramification
inside \(\mathbb G_m\), but the Keller boundary has nowhere to put them.

## 4. Closing the one-puncture conductor case

Degree one does not by itself make \(B\) normal.  The normalization can
identify distinct points while remaining immersive at each branch.  For
example,

\[
 t\longmapsto (t^2-1,\;t(t^2-1))
\]

has image

\[
 y^2=x^2(x+1).
\]

It has normalization \(\mathbb A^1\), one place at infinity, and nonzero
derivative everywhere, but it identifies \(t=1\) and \(t=-1\) at the node.
Its puncture and rank-one valuation lattices are the same as those of a
normal affine line.  Thus divisorial saturation alone does not detect this
escape.

This example is not a Keller map.  It is a countermodel to an attempted
proof from puncture and valuation data alone.  It shows that the remaining
hypothesis must control either the conductor or the complete sheet budget,
not merely the boundary's generic DVR.

If conductor gluing is absent, \(B\simeq\mathbb A^1\).  The line-component
theorem says that a polynomial map of \(\mathbb C^2\) whose nonproperness set
has such a component must have a critical point.  That contradicts the
Keller condition.

It remains to exclude conductor gluing without assuming the target curve is
normal.

### 4.1 Flatness and the generic sheet budget

The surface \(\overline X\) is smooth: \(U=\mathbb A^2\) is smooth and
condition 2 covers its complement \(D\).  A finite morphism from the
Cohen--Macaulay surface \(\overline X\) to the regular surface
\(\mathbb A^2\) is flat by miracle flatness.  It is therefore finite locally
free of constant degree \(d\), and every closed fiber has scheme length
\(d\).

At the generic point of \(B\), the degree formula for the finite extension
is

\[
 d=\sum_i e_i f_i.
\]

There must be an affine term.  If \(h=0\) defines \(B\), then
\(h\circ F\) is a nonconstant polynomial and therefore has a divisor in
\(U\).  Since \(F\) is etale, every such affine prime has ramification index
one.  Condition 5 says that the saturated link selects exactly one of them,
with primitive residue degree one, and that the divisor list is exhaustive.

The boundary term is \(e\cdot1\).  Every affine prime is unramified because
it lies in the Keller open.  Write

\[
 a=\sum_{\text{affine }E_i} f_i.
\]

Then \(d=e+a\).  Condition 5 gives \(a=1\), so

\[
 \boxed{d=e+1.}
\]

This is the operational content of a **primitive saturated rank-one minimal
link**.  A rank-one unit lattice alone does not imply this formula.

### 4.2 A conductor collision costs two boundary packets

Suppose the normalization identifies distinct points \(p,q\in D\) over
\(b\in B\).  Work in completed local rings.  The image of the germ
\((D,p)\) is a smooth analytic branch because the residue map is immersive.
Choose regular target parameters \((u,v)\) such that this branch is \(v=0\)
and \(u\) restricts to a parameter on it.  At \(p\), parameters \((t,z)\)
with \(D=(z=0)\) satisfy

\[
 \pi^*u=t+\text{terms divisible by }z
\]

after rescaling \(t\), while

\[
 \pi^*v=z^e A(t,z)
\]

with \(A\) not divisible by \(z\).  Eliminating \(t\) shows that the
local fiber algebra

\[
 \operatorname{length}
 \mathcal O_{\overline X,p}/(\pi^*u,\pi^*v)
\]

has length at least \(e\).  The same holds at \(q\), using coordinates
adapted to its branch.  Since both local algebras are
summands of the finite flat fiber over \(b\),

\[
 d=\operatorname{length}\pi^{-1}(b)\ge 2e.
\]

In general, \(d=e+a\) turns this into the necessary inequality
\[
 \boxed{a\ge e.}
\]

### 4.3 Why total conductor reduces to two points under immersion

For a rational plane curve of degree \(c\), the normalization conductor has
total degree

\[
 (c-1)(c-2).
\]

This number alone does not force many distinct points: it can all be
concentrated at one unibranch cusp.  Such a cusp is unavailable here,
however, because the normalization map of a singular plane branch has zero
differential at its unique preimage, whereas Section 2 proves that the
boundary residue map is immersive everywhere.

The least complicated remaining conductor profile is therefore one
singular target point with two smooth normalization branches.  Arbitrarily
large conductor can still be stored in their contact order, so no
three-point bound follows.  But the two distinct normalization points are
exactly the packet used in Section 4.2.  Hence the inequality \(a\ge e\)
applies to the entire immersion-compatible Pareto-minimal conductor face,
not merely to an ordinary node.

The executable enumeration and the direct degree-independent Pareto proof
are in
[`cas/target_conductor_atlas.py`](cas/target_conductor_atlas.py).

If \(e=1\), the finite map is etale on \(U\), and the residue-Jacobian
argument proves it etale along \(D\) as well.  Hence \(\pi\) is a connected
finite etale cover of \(\mathbb A^2\), so it has degree one; this is
incompatible with a nonempty boundary.  Therefore \(e\ge2\).

The primitive minimal link has \(a=1<e\), a contradiction.  Equivalently,
using \(d=e+1\), one gets \(2e\le e+1\), impossible for \(e\ge2\).

This closes the conductor case and proves the theorem.

### 4.3 Conductor packets

The two-point estimate is the first case of a local finite-flat statement.
Let \(b\in B\), and let

\[
 p_1,\ldots,p_r\in \overline X\setminus U
\]

be distinct boundary points over \(b\).  At \(p_i\), assume:

1. the source and boundary germ are regular;
2. the corresponding analytic branch of \(B\) is smooth;
3. the residue map is immersive; and
4. the pullback of a transverse equation of that branch has order \(e_i\).

Then the local fiber algebra at \(p_i\) has length at least \(e_i\).
The local algebras are distinct summands of the finite flat fiber over
\(b\), and therefore

\[
\boxed{d\ge \sum_{i=1}^r e_i.}
\tag{9}
\]

This is the **conductor-packet inequality**.  It is more accurately a
finite-flat packet inequality: the conductor is what supplies several
normalization points in one target fiber, while flatness and transverse
multiplicity supply the estimate.

The generic and special-fiber data must be kept separate.  If the complete
generic boundary contribution over \(B\) is

\[
 E_B=\sum_D e_Df_D
\]

and the total affine contribution is \(a_B\), then

\[
 d=E_B+a_B,\qquad
 a_B\ge \sum_{p_i\mid b}e_i-E_B.
\tag{10}
\]

For one boundary divisor of generic index \(e\) and residue degree \(f\),
with a constant-index \(r\)-point packet, this becomes

\[
 d=ef+a_B,\qquad a_B\ge (r-f)e.
\tag{11}
\]

After residue immersion and the one-puncture Riemann--Hurwitz argument force
\(f=1\), the useful form is

\[
\boxed{a_B\ge(r-1)e.}
\tag{12}
\]

Here \(r\) is the number of normalization points over one *closed* singular
value; it is not the generic residue degree \(f\).  For \(r=2\), (12)
recovers the current estimate \(a_B\ge e\).

No logarithmic-purity hypothesis is logically needed once assumptions
1--4 above are supplied point by point.  In the boundary theorem,
logarithmic purity is the mechanism that proves those local assumptions
uniformly.  Without residue immersion the estimate can fail in the proposed
form: a singular unibranch image may have only one normalization point, and
the tangential differential can spend ramification at that point.

### 4.4 Arbitrary-puncture residue budget

The canonical finite-normalization programme supplies an additional input:
every normalization-boundary curve has rational projective normalization.
Suppose one such curve has residue degree \(f\) over a parametric target
nonproperness curve and has \(s\) punctures over target infinity.  If the
positive pole orders are \(m_1,\ldots,m_s\), then

\[
\sum_jm_j=f.
\]

Riemann--Hurwitz forces ramification on the affine residue curve of total
degree

\[
(2f-2)-\sum_j(m_j-1)=f+s-2.
\tag{13}
\]

Thus affine residue immersion is compatible only with \(f=s=1\).
The former one- and two-puncture calculations are the first two cases of
(13); no bounded puncture-rank assumption is needed once rationality and
residue immersion have been transferred to the canonical normalization.
The exact checker now enumerates every ordered positive pole profile through
`puncture_profile_budgets`.

## 5. Smallest numerical test

The intrinsic affine-plane boundary gate already contains the first
numerically admissible one-dicritical package.  Starting with
\((\mathbb P^2,L)\), make three successive free boundary blowups.  In
boundary order

\[
 (L,E_1,E_2,E_3)
\]

the canonical and target-pole vectors are

\[
 k=(-3,-2,-1,0),\qquad p=(3,2,1,0).
\]

The exact gate gives

\[
 Qp=(2,0,0,1),\qquad p^tQp=6.
\]

Thus:

- degree \(6\) is the first numerically possible free-depth package;
- \(E_3\) is its unique dicritical;
- \(E_3\) has one puncture;
- its target-line degree is \(1\).

The last equality forces both the residue degree and the projective image
degree to be one.  Hence the image is a line, so the line-component theorem
excludes this package.  This is a genuine strengthening of the existing
intersection-matrix gate: that gate accepts the degree-six divisor data,
while the residue/conductor gate rejects it.

The executable arithmetic audit is
[`cas/plane_boundary_exclusion.py`](cas/plane_boundary_exclusion.py), with
regression
[`cas/test_plane_boundary_exclusion.py`](cas/test_plane_boundary_exclusion.py).

Actual coefficient searches in generic degree \(3\)--\(8\) are unnecessary:
the published Newton reductions already place any plane counterexample far
above this range.  Low degree remains useful here as a structural regression
for the new boundary obstruction.  Under condition 5, a conductor collision
in generic degree \(d\) would require

\[
 2(d-1)\le d.
\]

Thus every degree \(3\)--\(8\) row is rejected before coefficient solving.

## 6. Relation to the proposed statement

The proposed statement becomes a theorem if its terms are read as follows:

| proposed phrase | operational requirement |
| --- | --- |
| one irreducible rational dicritical component | the Zariski--Main boundary has one prime \(D\), and \(B=\pi(D)\) is the unique nonproperness component |
| rational normalization with finitely many punctures | \(D\) is smooth with completion \(\mathbb P^1\); the canonical \(\mathbb A^2\) boundary theorem supplies this rationality |
| tame generic ramification | logarithmic purity along every point of \(D\), not only generic tameness |
| saturated rank-one valuation link | the linked affine prime is primitive in both normal and residue degree: \((e,f)_D=(e,1)\), \((e,f)_E=(1,1)\) |
| no additional divisorial components | the generic pullback of \(B\) has exactly the linked pair \(D,E\), and the regular graph audit has no exceptional residual-ramification divisor |

If "saturated" means only the rank-one unit or valuation lattice, the
one-puncture statement is still one condition short: the nodal conductor
model above is invisible to that lattice, and the degree formula \(d=e+1\)
does not follow.  The full operational dichotomy is

\[
\boxed{
\begin{array}{c}
\text{two-or-more-puncture boundary}\\
\text{excluded by residue immersion + Riemann--Hurwitz}
\end{array}
}
\qquad
\boxed{
\begin{array}{c}
\mathbb A^1\text{ boundary}\\
\text{degree one; conductor gluing violates }d=e+1
\end{array}
}.
\]

The theorem therefore proves the proposed statement under a checkable
meaning of "saturated" and "no additional divisors."  If the intended
meaning omits primitive residue degree or permits further affine primes over
\(B\), the remaining escape is equally precise: a nodal or multi-branch
one-place curve together with enough extra generic sheets to satisfy
\(d\ge2e\).

### 6.1 Typed finite-normalization gate

The executable audit now separates source-resolution data from the
finite-normalization hypotheses.  A
`OneDicriticalNormalizationCertificate` records

\[
(d,e,f,a;\ \text{punctures},\text{one normalization boundary},
  \text{log purity},\text{exhaustiveness},\text{target transfer}),
\]

and first checks the generic degree identity

\[
d=ef+a.
\]

It refuses to apply the theorem when a source dicritical has not been
transported through the resolved target graph to the finite Zariski--Main
normalization over the original affine target.  This is intentional: a
primitive source pole vector certifies neither that the normalization
boundary has one prime nor the affine contribution \(a\).

Once the target transfer, exhaustive pullback, one-place theorem, and log
purity are certified, the gate applies Riemann--Hurwitz.  Every profile with
two or more punctures is excluded; one puncture forces \(f=1\).  The image is then either normal, in
which case the affine-line component theorem excludes it, or nonnormal, in
which case a conductor collision requires

\[
2e\le d=e+a.
\]

Thus \(a<e\) is the exact sheet-deficient kill condition.

The source-selected `(72,108)` packages are included only as typed previews.
Their ledgers are

\[
29=3+26\quad\text{and}\quad29=5+24.
\]

The audit reports the source records as incomplete because the target
transfer is not yet certified.  Even if that transfer is supplied for the
locally log-pure Case 1 row, its affine contribution satisfies \(a\ge e\),
so it survives the conductor budget.  Case 2 additionally has a residual
special-point ramification factor and fails log purity.

Run the combined regression with

```bash
.venv/bin/python plane-jc/cas/test_plane_boundary_exclusion.py
```

## 7. Sources and status

- Nguyen Van Chau,
  [*Non-proper value set and the Jacobian condition*](https://arxiv.org/abs/math/0305088):
  a nonempty nonproperness set of a nonsingular plane polynomial map has one
  point at infinity.
- Nguyen Van Chau,
  [*Plane Jacobian conjecture for simple polynomials*](https://arxiv.org/abs/0711.3894):
  in particular, a nonsingular plane polynomial map cannot have a
  nonproperness component isomorphic to the affine line.
- Nguyen Van Chau,
  [*A note on singularity and non-proper value set of polynomial maps of
  \(\mathbb C^2\)*](https://arxiv.org/abs/0710.5212):
  the direct Newton--Puiseux/Abhyankar--Moh form of the affine-line
  component obstruction.
- Jorge A. Guccione, Juan J. Guccione, Rodrigo Horruitiner, and Christian
  Valqui,
  [*Increasing the degree of a possible counterexample to the Jacobian
  Conjecture from 100 to 108*](https://arxiv.org/abs/2204.14178):
  the published Newton reduction leaves only the \((72,108)\) row below
  larger degree \(125\); the repository separately audits the subsequent
  exact certificate for that row.
- Alexander Borisov,
  [*Frameworks for two-dimensional Keller maps*](https://arxiv.org/abs/1901.04073):
  compactification dual graphs and Picard-group constraints for hypothetical
  plane counterexamples.
- The
  [Stacks Project miracle-flatness lemma](https://stacks.math.columbia.edu/tag/00R4)
  and
  [finite locally free degree formalism](https://stacks.math.columbia.edu/tag/02K9)
  supply the constant fiber-length step.

As of July 23, 2026, the unrestricted complex plane Jacobian problem remains
open.  The recently announced three-dimensional counterexample and the
classification of graded plane Keller maps do not settle this ungraded
two-dimensional conductor problem.

## 8. What a boundary-leaf reduction would have to prove

The conductor packet gives a precise target for a global reduction, but the
weighted source tree and pole vector do not by themselves contain all of its
inputs.  For a target nonproperness component \(B\), define the packet
deficiency at \(b\in B\) by

\[
 \Delta(B,b)=
 \sum_{p\in(\overline X\setminus U)\cap\pi^{-1}(b)}e_p-d.
\tag{13}
\]

A positive value is an immediate contradiction.  Equivalently, using the
generic decomposition, the packet demands more affine degree than the
available budget:

\[
 \sum_{p\mid b}e_p-E_B>a_B.
\tag{14}
\]

Thus a useful boundary-leaf theorem for a minimal counterexample would need
to force one component \(B\) for which:

1. the relevant source dicriticals are transported to the finite
   Zariski--Main normalization and grouped over \(B\);
2. the residue curve has one puncture and the pointwise logarithmic
   differential is pure, so its normalization map is immersive;
3. the generic factors \(e_D,f_D\) and the affine remainder
   \(a_B=d-E_B\) are exhaustive; and
4. some singular value has packet weight greater than \(d\).

The complete source intersection matrix \(Q\) and pole vector \(p\) do
provide important parts of this ledger:

\[
 d=p^tQp,\qquad
 R_f=K_X+3f^*L,
\]

and, on a genuine resolved morphism, the coefficient of a dicritical in
\(R_f\) recovers its generic transverse index minus one.  They also identify
dicritical candidates by \(p_i=0\), \((Qp)_i>0\).

Three required pieces are not numerical invariants of \((Q,p)\):

- the factorization of the restriction degree into residue degree and
  projective image degree;
- the grouping of dicriticals over the same target component and of their
  special points over the same \(b\); and
- pointwise residual ramification, which is a section-level condition rather
  than only a divisor coefficient at the generic point.

Consequently the next compiler layer should be a **target-transfer ledger**,
not another source-lattice inequality.  It should record for each
dicritical:

\[
(B;\ e_D,f_D;\ \text{punctures};\
 \text{residual-Jacobian divisor};\
 \{(b,p,e_p)\}),
\]

together with the affine primes above \(B\).  The packet checker can then
evaluate (13) exactly.

Minimality alone presently gives no known reason for \(a_B=1\), nor for a
packet large enough to violate (14).  In particular, rationality and the
one-place-at-infinity theorem force a singularity when the image is not an
affine line, and immersion forces that singularity to have at least two
normalization branches; they do not force more than the already proved
two-packet bound.  The genuinely new global lemma must therefore constrain
the affine remainder \(a_B\), force a larger coincident packet, or relate the
two through monodromy or the boundary intersection pairing.

### 8.1 Conductor length does not enlarge the packet

It is tempting to replace the number of normalization points by the
conductor length or the \(\delta\)-invariant.  That does not strengthen the
finite-fiber estimate without an additional surface-level argument.

For every \(m\ge2\) and \(c\in k\), consider

\[
 t\longmapsto
 \left(
 x=t^2-1,\quad
 y=c(t^2-1)+t(t^2-1)^m
 \right).
\tag{15}
\]

Its image has equation

\[
 (y-cx)^2=x^{2m}(x+1).
\tag{16}
\]

The parametrization is birational: away from \(x=0\),

\[
 t=\frac{y-cx}{x^m}.
\]

It is also immersive everywhere.  Since \(x'=2t\), only \(t=0\) needs
checking, and there \(y'=(-1)^m\ne0\).  The two points \(t=1,-1\) map to
the origin.  Their smooth image branches have common tangent \(y=cx\) and
intersection multiplicity \(m\).  Thus
\(\delta=\operatorname{length}(\widetilde{\mathcal O}/\mathcal O)=m\);
because this is a plane-curve singularity, the conductor has colength
\(2m\) in the normalization.  Both grow arbitrarily while the conductor
packet still has only two points.

The coordinate degrees \(2\) and \(2m+1\) are coprime, so the projective
closure has one place at infinity.  This family therefore has exactly the
qualitative curve features available in the plane-Jacobian reduction:
rational normalization, one point at infinity, residue immersion, and a
large conductor.  Nevertheless its packet bound remains only

\[
 d\ge e_{+}+e_{-},
\]

or \(d\ge2e\) at constant transverse index.

Consequently neither the arithmetic genus of the image curve nor total
conductor length can by itself produce the desired sheet deficit.  Any
stronger inequality would have to show that tangency of target branches
forces extra length in the *surface-map fiber*.  The local normal form used
in Section 4.3 shows that this is not a formal consequence of residue
immersion: after choosing a parameter along either branch, its local fiber
contribution can still be exactly \(e_i\).

### 8.2 Executable target-transfer ledger

The checker now implements the proposed target-side object as a
`TargetComponentLedger`.  It keeps four kinds of data separate:

1. generic boundary primes \((D,e_D,f_D)\);
2. generic affine primes and their residue degrees;
3. named, distinct boundary points grouped into certified closed target
   fibers; and
4. the finite-flatness, target-transfer, exhaustiveness, and pointwise
   residue-immersion certificates required to apply the inequality.

For each component it checks

\[
d=\sum_D e_Df_D+\sum_A f_A
\]

and, independently for every certified packet over \(b\),

\[
\sum_{p\mid b}e_p\le d.
\]

It deliberately reports `incomplete` when a source-only packet has not been
transported, when the points are not known to share one target fiber, or
when residue immersion is missing at even one listed point.  The regression
contains an excluded primitive \(d=3,e=2,a=1\) ledger, a paid multi-boundary
packet, and an uncertified source-only preview.
