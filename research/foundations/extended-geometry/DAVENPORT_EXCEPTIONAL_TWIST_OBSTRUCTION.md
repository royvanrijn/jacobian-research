# No exceptional degree-seven twist of the Davenport action

The degree-seven Davenport point/line monodromy cannot be converted into an
exceptional permutation cover merely by an arithmetic twist.  The reason is
finite-group rigidity: the natural degree-seven copy of
`GL_3(F_2)` is self-normalizing in `S_7`.

Thus the permutation-reduction part of the Cox programme must change the
geometric permutation action, not just the field of definition or the
orientation character.

## 1. Point and line actions

Let

\[
 G=GL_3(\mathbb F_2),\qquad |G|=168.
\]

It acts transitively on the seven nonzero vectors of `F_2^3`, equivalently
on the seven points of the Fano plane.  It also acts on the seven nonzero
covectors, equivalently on the seven Fano lines, through

\[
 M\longmapsto(M^{-1})^{\mathsf T}.                   \tag{1}
\]

The point and line actions have the same image as an abstract subset of
`S_7`, because (1) is an automorphism of `G`.  Their difference is the
identification of an abstract element with that image: (1) is the outer
duality automorphism and is not implemented by a permutation conjugation
of the seven sheets.

Exact enumeration gives

\[
 \boxed{N_{S_7}(G)=G.}                               \tag{2}
\]

Equivalently, the homomorphism

\[
 N_{S_7}(G)/G\longrightarrow\operatorname{Out}(G)
\]

has trivial image in this permutation representation.  Fano duality does
not extend the degree-seven point action inside `S_7`.

## 2. Arithmetic-monodromy consequence

Let

\[
 \phi:X\longrightarrow Y
\]

be a separable degree-seven cover whose geometric monodromy is this point
action of `G`.  Write `A<=S_7` for its arithmetic monodromy.  Geometric
monodromy is normal in arithmetic monodromy, so

\[
 G\mathrel{\triangleleft}A\leq S_7.
\]

Every element of `A` therefore normalizes `G`, and (2) forces

\[
 \boxed{A=G.}                                        \tag{3}
\]

There is no nontrivial constant-field quotient `A/G` and hence no
Frobenius coset which could turn the point action into an exceptional one.
In the equivalent fixed-point test, elements of `G` have

\[
 0,\ 1,\ 3,\ \text{or }7
\]

fixed points, rather than exactly one throughout an arithmetic coset.

The same conclusion holds for the line action.

### Theorem 2.1

No arithmetic twist which preserves the degree-seven geometric
`GL_3(F_2)` point or line action can be an exceptional cover.  In
particular, the oriented Davenport Cox maps cannot acquire permutation
reductions at infinitely many good primes by changing only constants,
orientation, or the point/line descent datum.

This is stronger than checking the particular Davenport equations: it
rules out every degree-seven cover with the same embedded geometric action.

## 3. Relation to the toric specialization

At `T=0`, the Davenport polynomial degenerates to a translated and scaled
seventh-power map.  Such a map can permute `F_q` when

\[
 \gcd(7,q-1)=1.
\]

There is no contradiction with Theorem 2.1, because the specialization
changes the geometric monodromy from `GL_3(F_2)` to a cyclic toric action.
It does not retain the generic Davenport--Sunada cover.

Consequently the remaining arithmetic construction must alter at least one
of:

1. the degree;
2. the geometric monodromy group;
3. the permutation representation; or
4. the cover through a degeneration which does not preserve the generic
   Davenport geometry.

## 4. Reproduction

Run

```bash
.venv/bin/python scripts/verify_davenport_exceptional_twist_obstruction.py
```

The checker enumerates `GL_3(F_2)`, its point and line images, all `7!`
sheet permutations, the full normalizer in `S_7`, and all point-action
fixed-point counts.
