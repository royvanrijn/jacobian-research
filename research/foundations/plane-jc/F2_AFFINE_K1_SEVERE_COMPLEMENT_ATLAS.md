# Severe complement atlas for the F2 `k=1` target

> **Status.** Exact computational-topological classification.  On the
> twelve severe rows left after the generic complement theorem, eight
> connected equisingular strata have affine complement group `Z`:
> `A2+A3+A1`, `A4+2A1`, `D4+A1`, `D4+A2`, `D5+A1`, `D6`, and the two `E7`
> realizations.  The four noncyclic rows are `A4+A3`, `A6+A1`,
> `A4+A2+A1`, and `D5+A2`.  At geometric degree six, with cubic local
> inertia, the first, second, and fourth admit no transitive action.  The
> third admits one simultaneous-conjugacy class, the natural degree-six
> action of `A6`.  Its fixed-sheet filling is excluded by the separate
> homology certificate.  Combined with the generic complement theorem and
> Chau's theorem on `E6/E8`, this closes every one-component degree-six
> `k=1` target.  It does not close larger geometric degree, extra target
> components, or `k>1`.

The implicit equations, braid continuation, Zariski--van Kampen groups, and
complete three-cycle enumeration are recomputed by
[`verify_f2_affine_k1_severe_complement_atlas.py`](../scripts/verify_f2_affine_k1_severe_complement_atlas.py).
The integral filling of the unique survivor is checked by
[`verify_f2_affine_k1_a4_a2_a1_filling_obstruction.py`](../scripts/verify_f2_affine_k1_a4_a2_a1_filling_obstruction.py).

## 1. Why one witness per row suffices

The
[`complete singularity atlas`](F2_AFFINE_K1_COMPLETE_SINGULARITY_ATLAS.md)
describes each row by factors of the collision quartic.  Each severe row is
connected over `C`, after identifying the two choices of a marked critical
point by `t -> -t`:

- `A2+A3+A1` is the open part of one linear-in-`c` discriminant factor;
- `A4+2A1` is an open subset of the line `c=6b+20`;
- `A4+A3`, `A6+A1`, and `A4+A2+A1` are single scaling orbits;
- `D4+A1`, `D4+A2`, and `D6` are nonempty opens in irreducible loci of the
  line--cubic merger factorization;
- the double-cubic locus is parametrized irreducibly by `(rho,b)`, giving
  `D5+A1`, with `D5+A2` and the two `E7` realizations as its three connected
  special sections.

On each row the projective quintic, its line at infinity, and its fixed
`(2,5)` infinity branch have constant embedded topological type.  Isolate
the affine singular points in disjoint Milnor balls.  Local equisingular
triviality and Thom's first isotopy lemma on the proper complement glue to
an ambient topological trivialization.  The affine complement group is
therefore constant along the connected row, so an exact computation at one
rational witness determines the row.

## 2. Exact complement groups

Sage/SIROCCO elimination and braid continuation give the following table.
The parameters refer to

\[
 p=t^3+at,\qquad q=t^5+bt^4+ct^2+dt.            \tag{2.1}
\]

\[
\begin{array}{c|c|c}
\text{packet}&(a,b,c,d)&\pi_1(\mathbb A^2-C)\\ \hline
A_2+A_3+A_1&(-3,-1/2,3,-9)&\mathbb Z\\
A_4+2A_1&(-3,-9/4,13/2,-9)&\mathbb Z\\
D_4+A_1&(1,1,0,-2)&\mathbb Z\\
D_4+A_2&(-3,2,0,3)&\mathbb Z\\
D_5+A_1&(-3,0,-2,-9)&\mathbb Z\\
D_6&(1,1,-1,-3)&\mathbb Z\\
E_7\text{ (intersection-three)}&(-3,-1,1,-7)&\mathbb Z\\
E_7\text{ (A4-cusp)}&(-3,2,-8,-13)&\mathbb Z\\ \hline
A_4+A_3&(-3,10,80,-205)&\text{noncyclic}\\
A_6+A_1&(-3,-7/2,-1,11)&\text{noncyclic}\\
A_4+A_2+A_1&(-3,5/2,-5,-5)&\text{noncyclic}\\
D_5+A_2&(-3,-2,4,-5)&\text{noncyclic}.
\end{array}                                                     \tag{2.2}
\]

In every cyclic row, the simplified presentation has one generator and no
relation.  This is exactly `Z`, not merely a cyclic quotient: the total
linking-number map of an irreducible affine plane curve surjects onto `Z`.

If such a curve were the complete ramified nonproperness set, its connected
finite cover over the complement would be a transitive action of `Z`.
The local meridian nevertheless fixes every affine sheet lying generically
over the curve.  A transitive cyclic permutation has no fixed point, so all
eight rows are impossible as the only ramified component, in every degree.

## 3. Complete cubic degree-six enumeration

For the four noncyclic groups, retain the three geometric meridians in the
raw Zariski--van Kampen presentation.  Cubic boundary inertia requires each
one to act in `S6` as a three-cycle, with three fixed affine sheets.
Fixing the first of the forty three-cycles and enumerating the other two
gives

\[
\begin{array}{c|c|c}
\text{packet}&\text{labeled solutions}&\text{conjugacy classes}\\ \hline
A_4+A_3&0&0\\
A_6+A_1&0&0\\
A_4+A_2+A_1&18&1\\
D_5+A_2&0&0.
\end{array}                                                     \tag{3.1}
\]

Every one of the eighteen labeled solutions in the survivor row generates
a group of order `360`.  They are the single natural `A6` class.  The
homological filling theorem computes its lifted presentation complex:

\[
 H_1(V;\mathbb Z)=\mathbb Z^2,
 \qquad
 H_1(V;\mathbb Z)/\langle\text{all fixed meridian lifts}\rangle
 =\mathbb Z.                                    \tag{3.2}
\]

Filling the affine sheets of a one-component Keller cover would recover
`A2`, whose first homology is zero.  Equation (3.2) excludes the sole
survivor.

## 4. Exact closure obtained

The complete `k=1` singularity table now has no uncomputed severe row.
Combining the independent inputs gives:

\[
\boxed{
 \text{No degree-six Keller cover has a one-component `k=1`
 nonproperness set.}}                           \tag{4.1}
\]

The qualifications in (4.1) are the remaining global frontier.  The F2
boundary calculation has not yet proved geometric degree six, has not
proved that the complete exceptional set is irreducible, and has not
excluded normalization contact `k=2,...,24`.

## Reproduction

With SageMath and SIROCCO installed:

```bash
sage scripts/verify_f2_affine_k1_severe_complement_atlas.py
.venv/bin/python scripts/verify_f2_affine_k1_a4_a2_a1_filling_obstruction.py
```
