# Intrinsic numerical boundary gates for `A^2`

## Result

Let `X` be a smooth projective surface over an algebraically closed field of
characteristic zero, let

\[
 D=D_1+\cdots+D_r
\]

be a reduced SNC divisor, and suppose

\[
 X\setminus D\simeq\mathbb A^2.
\]

Write

\[
 Q=(D_i\cdot D_j)_{i,j}.
\]

Then the intrinsic numerical boundary package satisfies all of the following.

1. The `D_i` are smooth rational curves and their dual graph is a tree.
2. Their classes form an integral basis of `Pic(X)`.
3. The intersection matrix is unimodular and has Hodge inertia
   \[
   \det Q=(-1)^{r-1},\qquad
   \operatorname{inertia}(Q)=(1,r-1,0).
   \]
4. Adjunction reconstructs the canonical class uniquely.  If
   \[
   K_X=\sum_i k_iD_i,
   \]
   then
   \[
   \boxed{
   Qk=(-2-D_1^2,\ldots,-2-D_r^2)^t.}
   \tag{1}
   \]
   In particular `k` is integral.
5. The reconstructed class obeys the rational-surface Noether identity
   \[
   \boxed{k^tQk=K_X^2=10-r.}
   \tag{2}
   \]

Conditions 1--5 are necessary, not sufficient.  They are nevertheless
strictly stronger than the existing Smith-normal-form gate.  For example,
the rational star with central weight `-1` and leaf weights `-5,-3,-2` has
determinant `-1` and inertia `(1,3,0)`, but (1) gives

\[
 K_X^2=-2
\]

instead of `10-4=6`.  It cannot be the boundary of `A^2`, despite passing the
unimodularity and Hodge-signature tests.

The exact implementation is
[`cas/intrinsic_a2_boundary.py`](cas/intrinsic_a2_boundary.py).

Here “intrinsic” means that, once the complete resolved source--target graph
is fixed (preferably its minimal smooth SNC model), `K_X`, ramification, and
the dicritical labels are reconstructed from its divisor and pullback data.
The weighted tree by itself is not invariant under optional extra boundary
blowups, and an unresolved list of Newton corners is not yet such a package.

## 1. Proof of the affine-plane gate

The standard compactification theorem transforms a smooth normal-crossing
completion of the complex affine plane to `(P^2,L)` by blowups and blowdowns
on the boundary.  In particular `X` is rational and `D` is a rational tree.
This boundary transformation statement is recalled in Section 1.4 of
M. Suzuki, *Affine plane curves with one place at infinity*,
[doi:10.5802/aif.1678](https://doi.org/10.5802/aif.1678), citing
Ramanujam and Morrow.

Divisor localization gives

\[
 0\longrightarrow
 \mathcal O(\mathbb A^2)^*/k^*
 \longrightarrow \mathbb Z^r
 \longrightarrow\operatorname{Pic}(X)
 \longrightarrow\operatorname{Pic}(\mathbb A^2)
 \longrightarrow0.
\]

The end groups vanish.  Hence the boundary classes are an integral basis of
`Pic(X)`.  The intersection form of a smooth rational projective surface is
unimodular with inertia `(1,r-1,0)`, proving items 2--3.

Since every `D_i` is rational, adjunction gives

\[
 K_X\cdot D_i=-2-D_i^2.
\]

Writing `K_X` in the boundary basis gives (1).  Finally a smooth rational
projective surface satisfies

\[
 K_X^2+\rho(X)=10.
\]

Here `rho(X)=r`, which proves (2).

There is also a useful log form.  If `v_i` is the valency of `D_i`, then

\[
 (K_X+D)\cdot D_i=v_i-2.
\]

Thus the entire numerical log-canonical class is recoverable from the
weighted dual tree; it is not an additional marking.

## 2. The Keller pole-vector gate

Now let

\[
 F:\mathbb A^2\longrightarrow\mathbb A^2
\]

be a polynomial Keller map.  Resolve the rational extension to obtain a
morphism

\[
 f:X\longrightarrow\mathbb P^2
\]

without changing the affine source, and let `L` be the target line at
infinity.  Because all base points lie on the source boundary,

\[
 f^*L=\sum_i p_iD_i,\qquad p_i\in\mathbb Z_{\ge0}.
\tag{3}
\]

Put `p=(p_i)` and retain the canonical vector `k` from (1).  Four further
intrinsic tests follow.

First, `f^*L` is nef and its square is the geometric degree:

\[
 \boxed{Qp\ge0,\qquad p^tQp=\deg(f)>0.}
\tag{4}
\]

Second, the Keller condition says that the ramification divisor has no
component on the affine source.  The canonical formula

\[
 K_X=f^*K_{\mathbb P^2}+R_f
\]

therefore gives the boundary coefficient vector

\[
 \boxed{R_f=k+3p\ge0.}
\tag{5}
\]

Third, the determinant of the logarithmic differential map

\[
 f^*\Omega_{\mathbb P^2}^2(\log L)
 \longrightarrow\Omega_X^2(\log D)
\]

has effective divisor, so

\[
 \boxed{R_{\log}=k+\mathbf1+2p\ge0.}
\tag{6}
\]

Finally, `D_i` maps to a curve meeting the affine target exactly when

\[
 \boxed{p_i=0,\qquad (Qp)_i>0.}
\tag{7}
\]

The first equality says that `D_i` is not contained in the inverse image of
`L`; the second says it is not contracted.  Such an index is the intrinsic
dicritical label.  A nonproper generically finite polynomial map must have at
least one.

For a Keller map this last assertion can be seen without adding a separate
purity hypothesis.  Zariski Main factors the quasi-finite map through an open
immersion

\[
 \mathbb A^2\hookrightarrow\overline X
\]

with `\overline X` normal affine and finite over the target.  If the
complement were nonempty of codimension at least two, normal Hartogs
extension would identify

\[
 \Gamma(\overline X,\mathcal O_{\overline X})
 =\Gamma(\mathbb A^2,\mathcal O_{\mathbb A^2}),
\]

so the open immersion would be an equality.  Hence nonproperness produces a
boundary divisor.  A finite morphism preserves its dimension, so its target
image is a curve; on the resolved projective graph it gives an index
satisfying (7).

Combining (5) and (7) gives the simple coordinate-free obstruction

\[
 \boxed{\text{every dicritical boundary prime has }k_i\ge0.}
\tag{8}
\]

Consequently:

> If the canonical vector reconstructed from the complete boundary
> intersection package is strictly negative in every coordinate, that
> package cannot support a nonproper plane Keller map.

This conclusion uses neither a preferred Newton coordinate nor a guessed
dicritical marking.

## 3. Free-point depth

The canonical labels also have a transparent blowup recursion.  Begin with
`(P^2,L)`, where `k_L=-3`.

- Blowing up a smooth point of one boundary component `D_i` creates an
  exceptional component with
  \[
  k_E=k_i+1.
  \]
- Blowing up an SNC crossing `D_i\cap D_j` creates one with
  \[
  k_E=k_i+k_j+1.
  \]

Before the first nonnegative label appears, a crossing blowup cannot create
one: if both parent labels are at most `-1`, their new label is again at most
`-1`.  Starting from `-3`, at least three one-parent boundary blowups must
occur along an ancestry path before a dicritical prime can occur; intervening
crossing blowups do not increase this free depth.

This is a small but unconditional classification result for the one-pair
frontier:

> Every complete Newton boundary package built on the standard
> `(P^2,L)` source resolution of a hypothetical plane Keller counterexample
> must contain a boundary prime of canonical free depth at least three.
> Packages with only shallower or crossing-created primes are impossible
> before any coefficient elimination.

After three free blowups the labels are `(-3,-2,-1,0)`, so the bound is
sharp at the numerical level.  The pole vector

\[
 p=(3,2,1,0)
\]

passes (4)--(7), has `Qp=(2,0,0,1)`, and gives numerical geometric degree
six.  This is only a consistency witness; it does not construct a Keller
map.

## 4. Executable use

Run

```bash
.venv/bin/python plane-jc/cas/test_intrinsic_a2_boundary.py
```

The regression:

1. checks exact inertia even for isotropic matrices such as the standard
   `P^1 x P^1` form;
2. verifies (1)--(2) after several smooth and crossing boundary blowups;
3. rejects the unimodular fake star above;
4. recovers zero ramification for the identity `P^2 -> P^2`;
5. verifies the canonical labels `-3,-2,-1,0` and the sharp numerical pole
   witness.

For a compiled Newton candidate the intended workflow is now

```text
complete source boundary and intersection matrix
    -> A2 tree/unimodularity/Hodge gate
    -> adjunction reconstruction of K_X
    -> K_X^2 + rho = 10
    -> target pole vector p
    -> nefness, degree, ramification, log ramification
    -> intrinsic dicritical test
    -> only then coefficient elimination
```

The remaining frontier gap is still the global one: the published Newton
case tree must be glued on compatible source and target completions before
its `Q` and `p` are complete inputs.  The new gate makes precise what should
be computed once that gluing is available.
