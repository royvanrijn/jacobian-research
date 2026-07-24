# JC(2) quartic finite-normalization packets

## Status

Orevkov's topology at infinity reduces every hypothetical geometric-degree
four plane Keller map to two global boundary packets.  It also excludes the
smallest clean one-boundary subcase:

> A one-boundary quartic packet consisting of one ordinary \(3+1\) cusp
> fiber and no \(2+2\) boundary self-collision is impossible.

Thus a surviving clean one-boundary packet must contain at least one target
value with two distinct simple ramification points and full fiber partition
\(2+2\).  The other global possibility has one ramified and one unramified
dicritical component, with no local-multiplicity jumps.

The exact arithmetic and finite permutation audit are implemented in
[`cas/plane_boundary_exclusion.py`](cas/plane_boundary_exclusion.py) and
regressed by
[`cas/test_plane_boundary_exclusion.py`](cas/test_plane_boundary_exclusion.py).

## 1. The global Orevkov budget

Let

\[
F:\mathbb C^2\longrightarrow\mathbb C^2
\]

be a polynomial local biholomorphism of geometric degree \(N\).  Use the
resolved and contracted topology-at-infinity model from Section 12 of
[`JC2_FINITE_NORMALIZATION_FRONTIER.md`](JC2_FINITE_NORMALIZATION_FRONTIER.md).
For every dicritical component \(\ell\), let \(\mu_\ell\) be its generic
local multiplicity and let \(\mu_x\) be the multiplicity at a finite special
point.  Orevkov's identity is

\[
\sum_{\ell}
\left(
 \mu_\ell+\sum_{x\in\ell^\circ}(\mu_x-\mu_\ell)
\right)=N-1,
\qquad
\mu_x-\mu_\ell\ge0.
\tag{1.1}
\]

Orevkov also remarks that a polynomial local biholomorphism cannot have a
dicritical component with

\[
\mu_\ell=N-1.
\tag{1.2}
\]

Now put \(N=4\).  Purity forces a ramified dicritical.  The positive affine
sheet budget rules out transverse index four, while (1.2) rules out index
three.  Hence there is exactly one ramified dicritical \(E\), with

\[
\mu_E=e_E=2.
\tag{1.3}
\]

Its residue degree is one: \(2f_E\le N-1=3\).  The remaining unit on the
right of (1.1) has exactly two possible uses.

| packet | dicritical multiplicities | exceptional jumps | Orevkov cost |
| --- | --- | --- | ---: |
| \(Q_{\mathrm{jump}}\) | \(2\) | one \(2\to3\) jump | \(2+1=3\) |
| \(Q_{\mathrm{two}}\) | \(2,1\) | none | \(2+1=3\) |

This proves:

### Proposition 1.1 -- quartic Orevkov dichotomy

Every hypothetical geometric-degree-four plane Keller map has precisely one
ramified dicritical \(E\) of transverse index two and belongs to exactly one
of the following rows.

1. \(E\) is the only dicritical and the total exceptional multiplicity jump
   on \(E\) is one.
2. There is exactly one further unramified dicritical \(D\), and every
   finite point of \(E\cup D\) has its generic local multiplicity.

In particular there cannot be two ramified dicriticals, a transverse
index-three dicritical, a jump of size at least two, or both a jump and an
extra dicritical.

This is a global statement.  It is stronger than enumerating the four
ramified degree-four signatures over one chosen target component, because it
simultaneously accounts for dicriticals mapping to different target curves.

## 2. The clean one-boundary row

Assume \(Q_{\mathrm{jump}}\) and suppose the unique exceptional point
\(p\in E\) is clean.  Proposition 9.1 of the cubic audit applies without
change.  If the first target order along \(E\) is \(m\), then the local fiber
length is \(m+1\).  The only available jump is

\[
\mu_p-\mu_E=1,
\]

so

\[
m=2,\qquad \mu_p=3.
\tag{2.1}
\]

The image branch is an ordinary cusp and finite flatness forces the complete
fiber partition

\[
\boxed{3+1.}
\tag{2.2}
\]

At the boundary point the length-three summand is curvilinear,
\(k[t]/(t^3)\); the other point is reduced and affine.  Thus the packet is
the cubic cusp collision plus one spectator sheet.

Over the generic point of the target curve \(C=(g=0)\), the affine
contribution has total residue degree two.  It can be one degree-two prime
or two degree-one primes.  Accordingly the global divisor has one of the
forms

\[
\operatorname{div}_B(g)=2E+A
\quad(f(A/C)=2),
\tag{2.3}
\]

or

\[
\operatorname{div}_B(g)=2E+A_1+A_2
\quad(f(A_i/C)=1).
\tag{2.4}
\]

Locally at the cusp, one affine sheet enters the curvilinear triple point
and the other is the reduced spectator.  If (2.3) holds, these two local
sheets belong to the same global affine prime.

There is one further special-fiber mechanism which costs no Orevkov jump.
Two distinct points \(x,y\in E\), each of local multiplicity two, can map to
the same target value.  Degree four then forces

\[
\boxed{\pi^{-1}(\pi(x))=2x+2y}
\tag{2.5}
\]

scheme-theoretically at the level of local lengths.  No affine point remains
over that value.  Each image branch is smooth by Orevkov's local link lemma,
although the two branches may have nontransverse contact.  We call (2.5) a
\(2+2\) boundary self-collision.

## 3. A lone cusp cannot connect the spectator

### Proposition 3.1 -- no-self-collision exclusion

In the clean one-boundary row, at least one \(2+2\) boundary
self-collision occurs.

#### Proof

Suppose none occurs.  The affine curve \(E^\circ\) is isomorphic to
\(\mathbb A^1\), the residue degree is one, and the normalization map

\[
E^\circ\longrightarrow C
\]

is injective away from the unique cusp.  An ordinary cusp is unibranch and
does not identify two normalization points.  Hence \(C\) is homeomorphic to
\(\mathbb A^1\), so it is a topologically contractible irreducible plane
curve.  The Lin--Zaidenberg theorem carries it by a polynomial
automorphism to

\[
C_0=(y^2-x^3=0).
\tag{3.1}
\]

The complement group has the braid presentation

\[
\pi_1(\mathbb C^2\setminus C_0)
=
\langle a,b\mid aba=bab\rangle,
\tag{3.2}
\]

where \(a,b\) are meridians.  Over the complement of \(C\), the finite
normalization is a connected degree-four étale cover: its source is a
nonempty open subset of the irreducible normalization surface.  Its
monodromy is therefore a transitive homomorphism

\[
\rho:\pi_1(\mathbb C^2\setminus C)\longrightarrow S_4.
\tag{3.3}
\]

Simple transverse ramification forces both \(\rho(a)\) and \(\rho(b)\) to
be transpositions.  Two transpositions satisfying the braid relation are
either equal or share exactly one letter.  In the first case they generate
\(S_2\); in the second they generate \(S_3\).  In either case they have an
orbit of size at most three on four letters.  This contradicts the
transitivity required in (3.3).

The checker enumerates all \(6^2\) ordered pairs of transpositions in
\(S_4\).  Exactly 30 satisfy the braid relation, their maximum orbit size is
three, and none is transitive.

### Proposition 3.2 -- one \(2+2\) collision saturates transitivity

The sheet-level monodromy obstruction goes no further.  A nondegenerate cusp
pair consists of two distinct transpositions satisfying the braid relation;
it generates \(S_3\) on three sheets and fixes the spectator.  A \(2+2\)
collision supplies two disjoint transpositions, hence one of the three
perfect matchings of four sheets.  Every perfect matching contains one edge
from the spectator to the cusp orbit.  Therefore the four transpositions
generate \(S_4\).

The exact enumeration has:

\[
24\ \text{ordered nondegenerate cusp pairs},\qquad
3\ \text{perfect matchings},\qquad
72\ \text{combined packets}.
\]

All 72 packets are transitive and all generate a group of order 24.  Thus
one \(2+2\) collision is both necessary and sufficient for transitivity at
the level of the local permutation generators.  Any further exclusion must
use the global relations among those meridians, or conductor/class-group
geometry; the sheet graph alone is exhausted.

## 4. Why Euler characteristic stops here

Let the clean one-boundary target curve have one cusp and \(r\)
double-branch self-collisions.  Its normalization is \(\mathbb A^1\).
The cusp does not change the topological Euler characteristic, while each
identification of two normalization points lowers it by one:

\[
\chi(C)=1-r.
\tag{4.1}
\]

Put

\[
\delta(y)=4-\#F^{-1}(y).
\tag{4.2}
\]

On the smooth part of \(C\), \(\delta=2\); at the cusp, \(\delta=3\); and at
a \(2+2\) collision, \(\delta=4\).  Since

\[
\chi\bigl(C\setminus\{\text{cusp and }r\text{ collisions}\}\bigr)
=-2r,
\]

Euler integration gives

\[
\int_C\delta\,d\chi
=2(-2r)+3+4r
=3.
\tag{4.3}
\]

But globally

\[
\int_{\mathbb A^2}\delta\,d\chi
=4\chi(\mathbb A^2)-\chi(\mathbb A^2)
=3.
\tag{4.4}
\]

Thus every number \(r\) is Euler-neutral.  A direct Suzuki comparison for
\(H=g(P,Q)\) repackages the same equality and cannot by itself exclude the
remaining \(2+2\) packets.  The next obstruction must retain monodromy,
conductor pairing, or the boundary class lattice.

## 5. The two-boundary row

In \(Q_{\mathrm{two}}\), the ramified component \(E\) has multiplicity two,
the unramified component \(D\) has multiplicity one, and neither has a local
multiplicity jump.  If they map to the same target component, the generic
ledger is forced to be

\[
(2,1)_E+(1,1)_D+(1,1)_{\mathrm{affine}}.
\tag{5.1}
\]

If their target images are different, the target component under \(E\) has
boundary contribution two and affine contribution two, while the residue
degree of \(D\) over its own target component can be one, two, or three,
subject to the positive affine remainder.

Orevkov's local link lemma makes every individual image branch smooth.
It does not prevent distinct boundary points from sharing a target value.
The target-transfer ledger must therefore retain at least:

- \(E\)-\(E\) coincidences, which have full packet \(2+2\);
- \(E\)-\(D\) coincidences, which consume length three and leave one affine
  point;
- whether \(E\) and \(D\) have the same target image; and
- self-identifications of the unramified residue curve.

This row belongs naturally to the two-generator free boundary class group
\(\mathbb Z[E]\oplus\mathbb Z[D]\).

## 6. Next closure target

The highest-value remaining statement is now:

\[
\boxed{
\substack{\text{No one-boundary quartic Keller normalization admits one
clean \(3+1\) cusp packet}\\
\text{together with any finite collection of \(2+2\) boundary
self-collisions.}}}
\tag{6.1}
\]

The first attack should encode the \(2+2\) pairings as transpositions on four
sheets.  Connectedness requires their edge graph to connect the cusp's
three-sheet \(S_3\)-orbit to the spectator.  At the same time every
self-collision consumes the full fiber and contributes conductor degree two
without spending an Orevkov jump.  The desired contradiction must compare
that conductor-pairing graph with the rank-one class relation (2.3) or
(2.4), rather than with Euler characteristic alone.

If (6.1) is closed, the only quartic survivor is \(Q_{\mathrm{two}}\), which
should then be compiled into the two-generator unimodular boundary lattice.

References:

- S. Yu. Orevkov,
  [*On three-sheeted polynomial mappings of
  \(\mathbb C^2\)*](https://doi.org/10.1070/IM1987v029n03ABEH000984),
  especially Lemma 4.2, Lemma 5.2, the proof of Theorem 1.1, and the closing
  remark.
- M. G. Zaidenberg and V. Ya. Lin,
  *An irreducible simply connected algebraic curve in
  \(\mathbb C^2\) is equivalent to a quasihomogeneous curve*,
  Soviet Math. Dokl. **28** (1983), 200--204.
