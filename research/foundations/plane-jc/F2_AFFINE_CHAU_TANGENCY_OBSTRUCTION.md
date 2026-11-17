# Chau tangency obstruction on the F2 affine target atlas

> **Status.**  Exact literature reduction plus exact symbolic
> classification.  Nguyen Van Chau proved that the exceptional-value set of
> a polynomial map of `C^2` with nonzero constant Jacobian cannot, in one
> affine target coordinate system, be a union of polynomial curves
> parametrized as `t -> (t^n,q(t))`.  On the normalized F2 `k=1` chart
> `p=t^3+a*t`, `q=t^5+b*t^4+c*t^2+d*t`, a linear target direction restricts
> to a pure power after affine reparametrization exactly on one of two loci:
>
> 1. `a=0`, where `p=t^3`; or
> 2. `25*c=2*b^3` and `125*d+50*a*b^2=b^4`, where
>    `q+(2*b^2/5)*p` is a translated fifth power.
>
> Consequently either locus is impossible when this curve is the complete
> exceptional-value set.  In particular the one-component `E_6+A_1` face,
> the monomial `E_8` endpoint, and every other `a=0` degeneration are
> excluded in every geometric degree.  This supersedes the apparent
> one-component `E_6` monodromy escape and makes the degree-six E8 cellular
> filling calculation an independent corroboration rather than the needed
> exclusion.  The theorem does not exclude an `a!=0` component off the pure
> quintic locus, a union containing another component not monomial in the
> same target coordinate, or any `k>1` chart.

The coefficient comparisons and first uncovered image-merger witness are
replayed by
[`verify_f2_affine_chau_tangency_obstruction.py`](../scripts/verify_f2_affine_chau_tangency_obstruction.py).

## 1. The external theorem and its exact scope

For a Keller map

\[
 F=(P,Q):\mathbb A^2_{\mathbb C}\longrightarrow\mathbb A^2_{\mathbb C},
\]

the critical-value set is empty, so its exceptional-value set `E_F` equals
its nonproperness set `A_F`.  Chau's Theorem 2 states:

\[
 \boxed{
 E_F\text{ cannot be a curve composed of images }
 t\longmapsto(t^n,q(t))
 \text{ in one target coordinate system}.}       \tag{1.1}
\]

His Corollary 3 separately says that `E_F` cannot be simply connected.  The
proof of (1.1) is stronger than that corollary.  It uses the preceding
tangency theorem to put all exceptional values of the first coordinate
polynomial in one value, Suzuki's equality for its fibres, and then the
Abhyankar--Moh and Gwozdziewicz injectivity theorems.

There are two scope points.

- Equation (1.1) is a statement about the **complete** exceptional set, not
  about an isolated component in an otherwise unknown union.
- For a reducible set, all components must have the displayed form in the
  same target coordinate.  Choosing an unrelated good coordinate on each
  component is insufficient.

These are exactly the qualifications retained below.

## 2. Complete linear-projection classification for `k=1`

Use the normalized chart

\[
 p=t^3+at,qquad q=t^5+bt^4+ct^2+dt.             \tag{2.1}
\]

A nonconstant linear target coordinate has degree three or five on the
normalization.

### 2.1 Cubic direction

The degree-three directions are the nonzero multiples of `p`.  Since

\[
 p'(t)=3t^2+a,                                  \tag{2.2}
\]

an affine reparametrization turns `p` into a pure cube exactly when

\[
 \boxed{a=0.}                                   \tag{2.3}
\]

On this whole hypersurface, (2.1) already has Chau form

\[
 t\longmapsto(t^3,t^5+bt^4+ct^2+dt).            \tag{2.4}
\]

### 2.2 Quintic direction

After scaling, every degree-five direction is `q+lambda*p`.  Compare it
with `(t-h)^5+h^5`, whose constant term is zero.  Coefficients of
`t^4,t^3,t^2,t` successively give

\[
 h=-\frac b5,qquad
 \lambda=\frac{2b^2}{5},qquad
 c=\frac{2b^3}{25},qquad
 d=\frac{b^4}{125}-\frac{2ab^2}{5}.             \tag{2.5}
\]

Thus the second and only other locus is

\[
 \boxed{
 25c=2b^3,qquad
 125d+50ab^2=b^4.}                              \tag{2.6}
\]

On (2.6), the invertible affine target change

\[
 (P,Q)\longmapsto
 \left(Q+\frac{2b^2}{5}P,,P\right)             \tag{2.7}
\]

and the parameter translation `s=t+b/5` put the normalization in the form
`(s^5, q_2(s))`, up to a harmless target translation.  Therefore Chau's
theorem applies.

Combining (2.3) and (2.6) gives a complete classification of the `k=1`
curves reached by a **linear** monomial projection.  Nonlinear target
automorphisms may enlarge the locus, but are not assumed here.

## 3. Named F2 faces closed at once

The former first noncyclic escape is

\[
 (p,q)=(t^3,t^5+t^4),                            \tag{3.1}
\]

with affine packet `E_6+A_1`.  It lies on `a=0`, so it cannot be the whole
exceptional set.  Its transitive degree-six `S_4` permutation action remains
a correct complement-group calculation; it merely does not evade (1.1).

The concentrated endpoint

\[
 (p,q)=(t^3,t^5),qquad C:P^5-Q^3=0,             \tag{3.2}
\]

is excluded twice: by Theorem 2 in the form (1.1), and by Corollary 3 since
the cusp is homeomorphic to `C`.  Hence

\[
 \boxed{
 \text{no one-component E8 exceptional set occurs at any degree}.} \tag{3.3}
\]

This includes the simple-inertia orbifold atlas and the cubic-inertia
natural `A_6` equality row.  Their local logarithmic modules and monodromy
atlases remain valid local classifications, but none can be a global
one-component Keller completion.

The conductor-conservation witness

\[
 (p,q)=(t^3,t^5+t^2),                            \tag{3.4}
\]

whose affine singularities are an ordinary cusp and a three-branch point,
is likewise on `a=0` and is excluded as a complete exceptional set.

## 4. A merger beyond the theorem

The collision formulas also isolate what Chau's obstruction does **not**
settle.  Reduce `q` modulo `p-v`, where `v` is a target value:

\[
 q\equiv
 (v-ab+c)t^2+(a^2+bv+d)t-av
 \pmod{p-v}.                                    \tag{4.1}
\]

All three points of a `p`-fibre have one `q`-value precisely when

\[
 v=ab-c,qquad
 \boxed{a^2+ab^2-bc+d=0.}                       \tag{4.2}
\]

The exact witness

\[
 a=1,quad b=c=0,quad d=-1,qquad
 (p,q)=(t^3+t,t^5-t)                            \tag{4.3}
\]

maps `t=0,i,-i` to `(0,0)`.  Its implicit equation factors suggestively as

\[
 \boxed{
 P^5=(Q+P)(Q+2P)^2,}                            \tag{4.4}
\]

and its `Q`-discriminant is

\[
 -P^8(27P^2+4).                                 \tag{4.5}
\]

It has `a!=0` and fails the second equation of (2.6), so it lies outside
both certified Chau loci.  It is an honest `k=1` image-merger target beyond
the theorem.  The subsequent
[`complete singularity atlas`](F2_AFFINE_K1_COMPLETE_SINGULARITY_ATLAS.md)
identifies it as a `D6` packet.  Its affine complement is cyclic, so the
separate complement-monodromy obstruction excludes it as the only ramified
component.  Chau's theorem itself does not supply that later exclusion.

## 5. Consequence for the `(75,125)` classification

The one-component `k=1` branch now splits as follows.

\[
\begin{array}{c|c}
\text{target stratum}&\text{current disposition}\\ \hline
\text{seven immersed/generic cusp strata with }\pi_1=\mathbb Z
  &\text{excluded as the only ramified component}\\
a=0\text{, including }E_6+A_1\text{ and }E_8
  &\text{excluded as the complete exceptional set by Chau}\\
\text{pure-quintic locus }(2.6)
  &\text{excluded as the complete exceptional set by Chau}\\
a\ne0\text{ severe nonimmersion/image mergers off }(2.6)
  &\text{classified by the complete singularity atlas}.
\end{array}                                      \tag{5.1}
\]

Multiple exceptional components remain open unless they are all monomial
in one common target coordinate.  The other 23 normalization charts also
remain open.  Thus this is a strict closure of the E6/E8 branch and a finite
reclassification of `k=1`, not an exclusion of `(75,125)` or a proof of
`JC(2)`.

## Literature

- Nguyen Van Chau,
  [*Two remarks on non-zero constant Jacobian polynomial maps of C2*](https://arxiv.org/abs/math/0408048),
  Theorem 1, Theorem 2, and Corollary 3.
- Nguyen Van Chau,
  [*Non-proper value set and the Jacobian condition*](https://arxiv.org/abs/math/0305088),
  for the one-point-at-infinity and `(3k,5k)` normalization input.
- S. Yu. Orevkov,
  [*On three-sheeted polynomial mappings of C2*](https://www.math.univ-toulouse.fr/~orevkov/jc86.pdf),
  for the finite-cover completion and Euler-multiplicity identity used by
  the surrounding F2 program.

## Reproduction

```bash
.venv/bin/python scripts/verify_f2_affine_chau_tangency_obstruction.py
```
