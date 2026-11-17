# Homological filling obstruction for the `A4+A2+A1` cubic row

> **Status.** Exact computational-topological exclusion, conditional only on
> the displayed Zariski--van Kampen presentation and the one-component
> hypothesis.  On the remaining severe `k=1` stratum with affine packet
> `A4+A2+A1`, there is exactly one transitive degree-six action in which every
> geometric meridian has cycle type `3+1+1+1`.  Its image is the natural
> `A6`.  The associated six-sheet cover has first homology `Z^2`.  Filling
> the three fixed-sheet peripheral lifts of any geometric meridian leaves
> first homology `Z`; adding all nine conjugate representatives leaves the
> same class.  If this curve were the complete nonproperness set of a Keller
> map, the filling would recover `A2`, whose first homology vanishes.
> Therefore the row is impossible.  Together with the cyclic-complement,
> Chau-tangency, and permutation exclusions, this closes the one-component
> degree-six `k=1` target atlas.  It does not exclude additional
> nonproperness components, larger geometric degree, or `k=2,...,24`.

The integral cellular calculation is replayed by
[`verify_f2_affine_k1_a4_a2_a1_filling_obstruction.py`](../scripts/verify_f2_affine_k1_a4_a2_a1_filling_obstruction.py).
The source presentation and the complete action enumeration can be
recomputed with SageMath and SIROCCO using
[`research_f2_affine_k1_severe_complements.py`](../scripts/research_f2_affine_k1_severe_complements.py).

## 1. The last degree-six severe stratum

Use the normalized `k=1` parametrization

\[
 p=t^3+at,
 \qquad q=t^5+bt^4+ct^2+dt.                    \tag{1.1}
\]

The exact witness

\[
 (a,b,c,d)=\left(-3,\frac52,-5,-5\right)       \tag{1.2}
\]

has affine singularity packet `A4+A2+A1`.  Elimination gives

\[
\begin{aligned}
F={}&P^5-\frac{55}{8}P^4+\frac{15}{2}P^3Q
 +\frac{15}{4}P^3-\frac{45}{2}P^2Q+15PQ^2-Q^3\\
 &+\frac{115}{4}P^2-45PQ+15Q^2
 +\frac{55}{4}P-\frac{33}{4}Q.                \tag{1.3}
\end{aligned}
\]

An exact Zariski--van Kampen computation, followed by Tietze reduction
without eliminating the three geometric meridian generators, gives

\[
 \pi_1(\mathbb A^2-C)
 =\langle x_1,x_2,x_3\mid r_1,r_2,r_3\rangle,  \tag{1.4}
\]

where, writing a signed integer `i` for `x_i` and `-i` for its inverse,

\[
\begin{aligned}
r_1={}&(3,-1,-3,-1,3,1),\\
r_2={}&(3,1,2,3,-1,-3,-2,-1),\\
r_3={}&(1,2,1,2,1,-2,-1,-2,-1,-2).            \tag{1.5}
\end{aligned}
\]

Its abelianized relation matrix has Smith diagonal `(1,1,0)`, as required
for an irreducible affine plane curve.

Enumerating triples of three-cycles in `S6`, imposing (1.5), and requiring
transitivity gives `18` labeled solutions after fixing the first cycle.
They form one simultaneous-conjugacy class.  A representative is

\[
\begin{aligned}
X_1&=(1\ 2\ 3),\\
X_2&=(1\ 4\ 5),\\
X_3&=(1\ 3\ 6).                                \tag{1.6}
\end{aligned}
\]

The generated group has order `360`; it is the natural `A6`.  Each
meridian has one ramified orbit of length three and three fixed sheets.
The other three noncyclic severe packets (`A4+A3`, `A6+A1`, and `D5+A2`)
have no such transitive degree-six action, so (1.6) is the sole degree-six
cubic-inertia survivor of the severe atlas.

## 2. Lifted cellular complex

Use the presentation two-complex of (1.4).  Its connected six-sheet cover
has six vertices and `3*6=18` oriented one-cells.  A spanning tree leaves

\[
 18-6+1=13                                      \tag{2.1}
\]

fundamental one-cycles.  Lift each of the three relators from every sheet.
The resulting `13 x 18` cellular relation matrix has Smith diagonal

\[
 \boxed{(1,1,1,1,1,1,1,1,1,1,1,0,0).}         \tag{2.2}
\]

There is no torsion, and hence the unfilled cover `V` satisfies

\[
 H_1(V;\mathbb Z)\simeq\mathbb Z^2.             \tag{2.3}
\]

Choose any one of the three geometric meridians.  Append its three closed
lifts based at the fixed sheets.  The `13 x 21` matrix has diagonal

\[
 \boxed{(1,1,1,1,1,1,1,1,1,1,1,1,0),}         \tag{2.4}
\]

so

\[
 H_1(V;\mathbb Z)/\langle\text{fixed-sheet meridians}\rangle
 \simeq\mathbb Z.                              \tag{2.5}
\]

The checker repeats (2.4) for each meridian generator.  It also appends all
nine fixed lifts simultaneously; the diagonal remains (2.4).  Thus the
surviving class is not an artifact of which conjugate meridian
representative is used.

## 3. Filling contradiction

Suppose that `C` is the complete nonproperness set of a degree-six Keller
map `f:A2->A2`.  Then

\[
 V=f^{-1}(\mathbb A^2-C)\longrightarrow\mathbb A^2-C             \tag{3.1}
\]

is the connected unramified cover defined by (1.6).  The length-three
peripheral orbit is the missing ramified boundary component.  The three
fixed generic sheets are precisely the affine pullback of `C`.  Adding
those affine curve components to `V` normally kills their meridians; it
recovers the affine source.  Consequently

\[
 H_1(\mathbb A^2;\mathbb Z)
 \simeq H_1(V;\mathbb Z)/
 \langle\text{fixed-sheet meridians}\rangle.    \tag{3.2}
\]

The left side is zero, whereas (2.5) is `Z`.  This contradiction excludes
the unique row (1.6).

This argument is insensitive to how the three fixed sheets are permuted by
a peripheral longitude.  Such a permutation only identifies conjugate
component meridians; the calculation with all nine representatives already
imposes the larger set of relations and still leaves `Z`.

## 4. Exact scope

The result completes the degree-six, one-component `k=1` monodromy audit:

- cyclic complement strata cannot support a connected cover with fixed
  affine sheets;
- all `a=0` and pure-quintic strata are forbidden as complete exceptional
  sets by Chau's tangency theorem;
- three noncyclic severe packets have no cubic degree-six action; and
- the only remaining natural-`A6` packet fails the filling test above.

The one-component and degree-six qualifications are essential.  Extra
nonproperness components add base relators and filling relations, while
larger degrees admit other permutation actions.  Neither is classified by
this certificate.

## Reproduction

```bash
.venv/bin/python scripts/verify_f2_affine_k1_a4_a2_a1_filling_obstruction.py
```
