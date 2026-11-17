# Atomic spectrum and non-generation of the quartic Keller map

Let \(k\) be a characteristic-zero field and write

\[
 \mathcal K_d(k)=
 \{F:\mathbb A_k^d\longrightarrow\mathbb A_k^d:
   \det DF\in k^\times\}
\]

for the Keller, equivalently étale polynomial, endomorphism monoid under
composition.  Its units are the polynomial automorphisms.  A nonunit is
**atomic** if every factorization in \(\mathcal K_d(k)\) has a unit factor.
It is **stably atomic** if every identity stabilization is atomic.

This note answers the generation question raised in the post-counterexample
[survey](https://jacobianconjectures.com/jacobian/note/) for the meaning of
*quasi-torus lift* used there: a Keller map admitting a nontrivial
positive-dimensional algebraic-torus source--target equivariance after
polynomial left--right changes.  It does not use *quasi-torus* to include a
purely finite diagonalizable symmetry.  That distinction is necessary
because the quartic below has a residual \(\mu _5\)-symmetry.

## 1. The atomic spectrum

Define the degreewise atomic spectrum

\[
 \operatorname{ASpec}_d(k)=
 \{N>1:\text{some atomic }F\in\mathcal K_d(k)
                 \text{ has geometric degree }N\}.
 \tag{1}
\]

For stable left--right classes, let

\[
 \mathfrak A_{d,N}^{\mathrm{st}}(k)
 =
 \{\text{stably atomic degree-\(N\) Keller maps}\}/
       \text{stable polynomial left--right equivalence}.
 \tag{2}
\]

> **Atomic-spectrum theorem in dimension three.**  Over every
> characteristic-zero field,
> \[
> \boxed{\operatorname{ASpec}_3(k)=\{3,4,5,\ldots\}.}
> \tag{3}
> \]
> More precisely:
>
> 1. every admissible weighted or root-engineered quadratic-gauge map of
>    geometric degree \(N\ge3\) has geometric monodromy \(S_N\), and is
>    absolutely and stably atomic;
> 2. for every \(N\ge4\), the boundary-clean weighted family supplies an
>    explicit \((N-3)\)-dimensional family in
>    \(\mathfrak A_{3,N}^{\mathrm{st}}(k)\) after geometric base change;
> 3. in degree three, the weighted, cancellation, and diagonal
>    quadratic-gauge mechanisms give the same foundational stable class,
>    while fiber-invisible cubic gauge lifts give infinitely many stable
>    atomic classes outside those three minimal mechanisms.

The exclusion of degree two is the
[geometric-degree spectrum theorem](GEOMETRIC_DEGREE_SPECTRUM.md).
The all-degree existence and stable atomicity follow from symmetric
monodromy and
[primitive-monodromy atomicity](PRIMITIVE_MONODROMY_ATOMICITY.md).
The \((N-3)\)-dimensional lower bound is the
[degreewise stable-moduli theorem](SAME_DEGREE_STABLE_INEQUIVALENCE.md).
All members of that weighted family have \(S_N\) monodromy, so the separated
stable classes are stably atomic classes, not merely noninvertible classes.

Every degree-three nonunit is atomic for the simpler numerical reason that
geometric degree is multiplicative and \(3\) is prime.  This proves
atomicity, not uniqueness of its stable class.  The
[low-rank boundary theorem](LOW_RANK_MULTIPLICITY_BOUNDARIES.md) proves only
that the three currently established construction mechanisms collapse to
the foundational class.  The later
[fiber-invisible cubic gauge theorem](UNIVERSAL_CUBIC_GAUGE_MULTIPLICITY.md)
proves that unrestricted cubic uniqueness is false.

Likewise, an admissible quadratic-gauge map does not factor into smaller
quadratic-gauge atoms: its \(S_N\) monodromy makes the map itself atomic.
Compositions of quadratic-gauge atoms give explicit decomposable maps of
product degree.  Recognizing every arbitrary map that admits such a word
remains open; the exact current algebraization criterion is the polynomial
sandwich criterion in
[imprimitive Keller factorization](IMPRIMITIVE_KELLER_FACTORIZATION.md).

### 1.1 Degrees that force every map to be atomic

The absence of degree two gives a stronger degree-only statement.  Put

\[
 \mathcal F=
 \{p:p\ge3\text{ is prime}\}
 \cup\{2p:p\text{ is prime}\}
 \cup\{8\}.
 \tag{3a}
\]

Thus

\[
 \mathcal F=\{3,4,5,6,7,8,10,11,13,14,17,19,22,\ldots\}.
\]

> **Forced-atomic/decomposable degree dichotomy.**
>
> 1. If \(N\in\mathcal F\), every characteristic-zero Keller map of
>    geometric degree \(N\), in every dimension, is absolutely and stably
>    atomic.
> 2. If \(N\ge3\) and \(N\notin\mathcal F\), there are integers
>    \(a,b\ge3\) with \(N=ab\).  Composing the explicit degree-\(a\) and
>    degree-\(b\) quadratic-gauge atoms gives a decomposable Keller map of
>    degree \(N\) on \(\mathbb A^3\), and hence in every dimension at least
>    three by identity stabilization.
>
> Consequently, in dimension three, the geometric degrees admitting a
> decomposable characteristic-zero Keller nonunit are exactly
> \[
> \boxed{\{ab:a,b\ge3\}.}
> \tag{3b}
> \]
> The first is \(9\).  Thus every degree-three through degree-eight Keller
> nonunit is atomic, while degree nine is the first degree in which atomic
> and decomposable maps coexist.

**Proof.**  In a nontrivial factorization \(F=G\circ H\), both factors are
Keller and

\[
 \operatorname{gdeg}F
 =\operatorname{gdeg}G\operatorname{gdeg}H.
 \tag{3c}
\]

Degree-one Keller maps are automorphisms, and the
Campbell--Razar--Wright Galois case excludes degree two.  Thus a
factorization into two nonunits forces both degrees to be at least three.

It remains to classify the integers with no factorization \(N=ab\),
\(a,b\ge3\).  An odd composite has such a factorization.  Write an even
integer as \(N=2m\).  If \(m\) is prime, this is the listed family \(2p\).
If \(m=4\), then \(N=8\).  Every other composite \(m\) can be written
\(m=rs\) with \(r\ge2\) and \(s\ge3\), and then
\(N=(2r)s\) is an allowed factorization.  This proves that the exceptional
integers are exactly \(\mathcal F\).

The argument survives every characteristic-zero field extension and every
identity stabilization because geometric degree is unchanged and the
degree-two Galois exclusion holds in every dimension.  This proves (1).
For (2), choose the explicit degree-\(a\) and degree-\(b\)
root-engineered quadratic-gauge maps.  Each has \(S_a\) or \(S_b\)
monodromy and is an atom, while their composition is visibly decomposable
and has degree \(ab=N\). \(\square\)

More generally, every ordered multiplicative partition

\[
 N=n_r\cdots n_1,\qquad n_i\ge3,
 \tag{3d}
\]

is realized by an explicit word of \(r\) quadratic-gauge atoms.  This is an
existence classification of quadratic-gauge factor words by degree.  It is
not yet an intrinsic recognition theorem for deciding whether a supplied
Keller map belongs to that generated submonoid.

## 2. A general atomic generation obstruction

The following elementary monoid lemma is the bridge from atomicity to
generation.

> **Lemma 2.1 (one-nonunit word lemma).**  Let \(M\) be a monoid equipped
> with a multiplicative degree
> \(\delta:M\to\mathbb Z_{\ge1}\) such that \(\delta(g)=1\) exactly for the
> units \(U\), and let \(a\in M\) be an atom.  If
> \[
> a=g_m\cdots g_1
> \tag{4}
> \]
> is a finite word, then exactly one \(g_i\) is a nonunit.  Consequently
> \[
> a=u_Lg_i u_R
> \tag{5}
> \]
> for units \(u_L,u_R\).

**Proof.**  Split (4) before its last factor.  Atomicity says that the last
factor or the preceding subword is a unit.  If the subword is a unit,
multiplicativity gives
\(\prod_{i<m}\delta(g_i)=1\), so every factor in it is a unit.  If the last
factor is a unit, remove it and repeat.  Since \(a\) is a nonunit, exactly
one nonunit remains. \(\square\)

Let \(\mathcal Q_d(k)\) be any left--right saturated collection of Keller
maps whose nonunits admit nontrivial positive-dimensional algebraic-torus
equivariance.  Lemma 2.1 gives

\[
 \boxed{
 a\in\langle\operatorname{Aut}(\mathbb A^d),\mathcal Q_d(k)\rangle
 \Longrightarrow
 a\text{ is polynomially left--right equivalent to one member of }
 \mathcal Q_d(k)
 }
 \tag{6}
\]

for every atomic \(a\).  Thus no putative composition-stable numerical
torus invariant is needed.  Primitive monodromy first collapses the whole
generator word to one nonunit left--right class; an intrinsic invariant can
then separate that single class.

## 3. The explicit non-quasi-torus atom

Over \(k\), put

\[
 t=1+xy,\qquad q=t^2z-y^2(1+3t)
\]

and define \(F=(F_1,F_2,F_3):\mathbb A^3\to\mathbb A^3\) by

\[
\boxed{
\begin{aligned}
F_1&=-\frac12tq,\\
F_2&=y-3xq-tq+2t^2x^2q^4,\\
F_3&=x(5-3t)+x^3z-(xq)^4.
\end{aligned}}
\tag{7}
\]

The root-engineered quadratic-gauge theorem and its exact certificates give

\[
 \det DF=1,\qquad
 \operatorname{gdeg}(F)=4,\qquad
 \operatorname{Mon}_{\mathrm{geom}}(F)=S_4.
 \tag{8}
\]

The map also sends the four distinct rational points

\[
 (0,1,5),\quad(-1,2,-9),\quad
 (1/3,-4,-27),\quad(2/3,-1,45)
 \tag{9}
\]

to \((-1/2,0,0)\).  Thus it is a noninvertible Keller map.  Since the natural
\(S_4\)-action is primitive, the primitive-monodromy theorem makes \(F\)
absolutely and stably atomic.

The
[intrinsic algebraic-torus exclusion](../cancellation/NO_ALGEBRAIC_TORUS_EQUIVARIANCE.md)
proves that no polynomial left--right representative of \(F\) admits a
nontrivial positive-dimensional algebraic-torus equivariance.  Apply (6):

> **Non-generation theorem.**
> \[
> \boxed{
> F\notin
> \left\langle
> \operatorname{Aut}(\mathbb A_k^3),\,
> \{\text{positive-dimensional quasi-torus Keller maps}\}
> \right\rangle .
> }
> \tag{10}
> \]

In particular \(F\) is not generated by automorphisms and any collection of
the currently known quasi-torus families.  The statement is stronger than
family-by-family separation: it excludes every noninvertible Keller
generator having such a torus-equivariant left--right representative.

## 4. What is and is not composition-preserved

Algebraic-torus equivariance itself is not safely compositional: two factors
may carry incompatible source, intermediate, and target actions.  The
composition obstruction used here is instead the intermediate-field
structure:

\[
 F=G\circ H,\qquad \operatorname{gdeg}G,\operatorname{gdeg}H>1
 \quad\Longrightarrow\quad
 k(F)\subsetneq k(H)\subsetneq k(x).
 \tag{11}
\]

This gives a nontrivial block system in the generic monodromy action.
Therefore a primitive-monodromy map cannot be a composition of two
nonunits.  Once the word has been reduced to one nonunit, the canonical
decorated-normalization torus obstruction is a genuine polynomial
left--right invariant and excludes that final generator.

The resulting atomic-spectrum programme has the following present boundary:

- the degree set is exactly known;
- the degrees forcing every map to be atomic are exactly
  \(\mathcal F\), and the first decomposable degree is nine;
- positive-dimensional families of stable atomic classes are known in every
  degree at least four;
- infinitely many stable atomic classes are known in degree three: the
  [fiber-invisible cubic lifts](UNIVERSAL_CUBIC_GAUGE_MULTIPLICITY.md) have
  prime geometric degree and pairwise distinct canonical boundary counts;
- all admissible quadratic-gauge maps are atoms;
- every ordered factorization of an integer into factors at least three is
  realized by a word of quadratic-gauge atoms;
- the foundational cubic is the unique class only inside the three
  established marked-root mechanisms;
- (7) is an explicit non-quasi-torus atom and proves non-generation by all
  positive-dimensional quasi-torus generators;
- the remaining degreewise classification problem is the geometry of the
  [candidate stable Keller-moduli object](../extended-geometry/STABLE_KELLER_MODULI_PROBLEM.md):
  its dimension, irreducible components, and componentwise generic
  monodromy.  The positive-dimensional weighted loci and the infinite cubic
  set do not determine those invariants;
- recognition of the full quadratic-gauge-generated submonoid remains open.

## 5. Reproduction

The explicit polynomial, collision, and torus-decoration calculations are
checked by

```bash
.venv/bin/python scripts/verify_root_engineered_quadratic_gauge.py
.venv/bin/python scripts/verify_universal_cubic_gauge_multiplicity.py
.venv/bin/python scripts/verify_quartic_monodromy.py
.venv/bin/python scripts/verify_quartic_algebraic_torus_obstruction.py
```

The monoid reduction in Lemma 2.1 and the passage from primitive monodromy to
atomicity are exact arguments, not bounded searches.
