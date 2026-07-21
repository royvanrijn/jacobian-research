# Exact image, fibers, and nonproperness set

Write a target as `(a,b,c)` and put

\[
P(T)=cT^3-2T^2+bT-2a,\qquad
Q=27a^2c^2-18abc+16a+b^3c-b^2.
\]

Then the exact statements are

\[
F(\mathbb C^3)=\mathbb C^3\setminus\Gamma,
\qquad S_F=V(Q),
\qquad
\Gamma=V(3bc-4,12a-b^2).
\]

Here `S_F` denotes the set of target limits of sequences escaping every compact
subset of the source.  Since `3bc=4` on `Gamma`, it can also be parametrized as

\[
\Gamma=\{(4/(27c^2),4/(3c),c):c\ne0\}.
\]

The fiber cardinalities (without multiplicity) are

\[
\#F^{-1}(a,b,c)=
\begin{cases}
3,&Q\ne0,\\
1,&Q=0\text{ and }(a,b,c)\notin\Gamma,\\
0,&(a,b,c)\in\Gamma.
\end{cases}
\]

## Proof of the fiber statement

On `x != 0`, put `t=y+1/x`.  The identities in `scripts/cubic_model.py`
show that `P(t)=0`, `P'(t)=2/x`, and that every simple root reconstructs one
and only one source point:

\[
x=2/P'(t),\quad y=t-P'(t)/2,\quad
z=5P'(t)^2/4-3tP'(t)/2-cP'(t)^3/8.
\]

Thus the `x != 0` points are in bijection with the simple roots of `P`.  On
`x=0` the map is

\[
F(0,y,z)=(z+4y^2,y,0),
\]

so a target has exactly one additional point `(0,b,a-4b^2)` precisely when
`c=0`.

Because `Disc(P)=-4Q`, off `V(Q)` the cubic has three simple roots when
`c != 0`.  When `c=0`, it has two simple finite roots and the additional
`x=0` point.  On `V(Q)`, away from `Gamma`, it has exactly one simple root
(with the `c=0` case counted using the extra point).  It has no simple root
exactly when it is a scalar multiple of a cube.  Comparing its coefficients
gives `3bc=4` and `12a=b^2`, hence exactly `Gamma`; here `c` is automatically
nonzero, so there is no `x=0` point.

## Proof of the nonproperness statement

Every point of `V(Q)` has a finite repeated root `t_0` of `P`.  Keep `b,c`
fixed, let `t` tend to `t_0`, and set

\[
a(t)=(ct^3-2t^2+bt)/2.
\]

Then `P(t)=0`, the targets tend to the given point, and `P'(t)` tends to zero.
Choosing `t` away from the finitely many zeros of `P'` gives reconstructed
points with `x=2/P'(t)` tending to infinity.  Hence `V(Q) subset S_F`.

Conversely, near a target outside `V(Q)`, every finite root stays simple, so
the reconstruction formulas keep its source branch bounded.  The only issue
is a root tending to infinity when `c` tends to zero.  Set `s=1/t`.  Its root
equation and reconstruction simplify to

\[
c=2s-bs^2+2as^3,
\quad x={s\over1-bs+3as^2},
\quad y=b-3as.
\]

The first equation has derivative `2` in `s` at `s=c=0`, so this branch is
regular there.  Substitution in the third coordinate gives
`z -> a-4b^2`; therefore it extends to the finite point
`(0,b,a-4b^2)`.  All local branches over a target outside `V(Q)` are bounded,
which proves `S_F subset V(Q)`.

Finally, direct Groebner reduction gives

\[
(Q_a,Q_b,Q_c)=(16a-b^3c,(3bc-4)^2)
\]

after replacing the generators by a lexicographic Groebner basis.  Its radical
is `(3bc-4,12a-b^2)`, so `Gamma` is exactly the singular locus of `V(Q)`.

Run `.venv/bin/python scripts/image_nonproperness.py` for the exact symbolic
certificates used above. The complementary
`.venv/bin/python scripts/verify_exceptional_fibers.py` makes the complete
degeneration case split executable: `c!=0` double/triple roots, `c=0` finite
quadratic roots, the simple projective root at infinity, and every denominator
in both reconstruction charts.

## Denominator-safe inclusion audit

The complementary script
`.venv/bin/python scripts/verify_image_nonproperness_inclusions.py` verifies
both image inclusions and both nonproperness inclusions using the following
exhaustive chart split:

| Source/root chart | Condition | Outcome |
|---|---|---|
| affine `x=0` | exactly over `c=0` | the finite point `(0,b,a-4b^2)` |
| finite root `t` | `P'(t)!=0` | one finite `x!=0` source point |
| finite root `t` | `P'(t)=0` | the escaping boundary divisor over `V(Q)` |
| projective root `t=infinity` | necessarily `c=0` | regular `s=1/t` branch returning to `x=0` |

This split checks the strata that raw denominator clearing could mishandle.
The elimination of the finite boundary is `(Q)`, but that fact is not used
alone: off `Gamma` its rational inverse is checked modulo `Q`; the failed
inverse denominator is handled explicitly by the triple-root parameter on
`Gamma`; and the leading-coefficient slice `c=0`, `a=b^2/16` is reached by
`t=b/4`. Conversely, every boundary parameter satisfies `Q=0`, so the
elimination introduces no extraneous target stratum.

For nonproperness, every repeated-root parameter supplies an exact punctured
escape path with `x=2/P'(t)`. Outside `V(Q)`, finite roots have `P'` invertible,
while the only possible root at infinity lies in the regular `s` chart above.
Thus no unexamined projective-root direction remains.

## Boundary normalization and lost sheets

Put `r=P'(t)=2/x`.  In `(t,r,c)` coordinates the graph becomes polynomial:

\[
a=t^2-ct^3+rt/2,\qquad b=4t+r-3ct^2,\qquad c=c.
\]

Consequently `r=0` is the dicritical boundary divisor.  Eliminating `t` from
its graph gives the principal ideal `(Q)`.  Its map to the discriminant is

\[
\nu:\mathbb A^2_{t,c}\longrightarrow V(Q),\qquad
(t,c)\longmapsto(t^2-ct^3,4t-3ct^2,c).
\]

This is the normalization map.  Indeed `t` is integral over the discriminant
ring because

\[
t^2-bt+3a=0,
\]

and away from `Gamma` it is recovered rationally as

\[
t={b-9ac\over4-3bc}.
\]

The exceptional locus of this rational inverse is `ct=2/3`; its image is
exactly `Gamma`.  The normalization is nevertheless one-to-one there: the
singularity is unibranch, not a crossing of boundary components.

Along a general transverse parameter `r -> 0`, the source valuations are

\[
\operatorname{ord}_r(x)=-1,\qquad
\operatorname{ord}_r(y)=0,\qquad
\operatorname{ord}_r(z)\ge1.
\]

Thus this single dicritical divisor accounts for the escaping sheets.  In the
projective cubic fiber, a smooth point of `V(Q)` has one simple root and one
double root: the double root is the length-two part lost at this boundary.
On `Gamma` the cubic has a triple root, so the entire length-three fiber lies
at the boundary.  For `c=0`, the remaining projective root is the regular
`s=1/t` branch and becomes the finite `x=0` source point described above.

## Monodromy

Local monodromy around a smooth point of `V(Q)` exchanges the two roots that
form its double root, hence is a transposition.  These transpositions generate
the full `S_3`: for example, in the slice `c=1`, take roots

\[
m+\sqrt\epsilon,\quad m-\sqrt\epsilon,\quad2-2m.
\]

Their sum is fixed at `2`, so they define targets in this family.  A loop of
`epsilon` around zero exchanges the first pair.  Choosing different collision
centers supplies adjacent transpositions, and therefore the monodromy group is
`S_3`.  At the triple-root curve, the transverse discriminant has cusp type;
the corresponding braid products give three-cycles.
