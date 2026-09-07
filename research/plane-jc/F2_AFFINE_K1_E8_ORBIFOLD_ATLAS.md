# F2 `k=1` `E_8` simple-inertia orbifold atlas

> **Status.**  Exact finite computation.  Imposing simple meridian inertia
> on the `(3,5)` cusp group produces an orbifold quotient of order `240`.
> All of its `30` subgroup-conjugacy classes are enumerated.  Exactly `13`
> transitive coset actions in the F2 range have a fixed sheet; their degrees
> are `6,10,12,15,20,24,30,40,60,120`.  The central order-four action and
> preferred longitude split their ramified peripheral orbits into rows
> `(e,f)=(2,f)` with `f` in `{1,2,4}`.  This exhausts every E8 action whose
> meridian is an involution.  It does not cover inertia greater than two and
> does not construct a Keller cover or exclude `(75,125)`.

The exact Todd--Coxeter and subgroup replay is
[`verify_f2_affine_k1_e8_orbifold_atlas.py`](../scripts/verify_f2_affine_k1_e8_orbifold_atlas.py).

## 1. The finite universal simple-inertia quotient

For the monomial cusp `P^5-Q^3=0`, put

\[
 G=\langle a,b\mid a^3=b^5\rangle,
 \qquad m=a^{-1}b^2,
 \qquad z=a^3=b^5.                               \tag{1.1}
\]

If every geometric-meridian cycle has length one or two, its permutation
image satisfies `M^2=1`.  Hence every such action factors through

\[
 \boxed{
 \Gamma=\langle a,b\mid a^3=b^5,\ (a^{-1}b^2)^2=1\rangle.}     \tag{1.2}
\]

Exact coset enumeration at the trivial subgroup gives

\[
 |\Gamma|=240.                                    \tag{1.3}
\]

In the regular action, `z` is central of exact order four.  The quotient by
`<z>` is the order-sixty icosahedral quotient from the degree-six theorem.
The earlier `A_5` atlas is therefore precisely the part of the present atlas
on which `z` acts trivially.

Every transitive permutation action of a finite group is its action on the
cosets of a subgroup.  The checker starts from the regular action, constructs
the exact multiplication table, enumerates subgroups by adjoining generators,
and canonicalizes each subgroup under all `240` conjugations.  It obtains
`30` subgroup-conjugacy classes, with subgroup-order census

\[
\begin{array}{c|rrrrrrrrrrrrrrrr}
|H|&1&2&3&4&5&6&8&10&12&16&20&24&40&48&120&240\\
\#&1&2&1&3&1&3&3&3&3&1&3&2&1&1&1&1.
\end{array}                                                    \tag{1.4}
\]

This makes the classification exhaustive without a degree cutoff search in
symmetric groups.

## 2. Peripheral normalization data

The preferred longitude is

\[
 \ell=zm^{-15}.
\]

Since `M^2=1` and `z` commutes with `M`, its image is

\[
 L=zM,
 \qquad \langle M,L\rangle=\langle M,z\rangle.                 \tag{2.1}
\]

Thus the central action is exactly the missing normalization datum.  A
ramified peripheral orbit has size `2f`: one meridian transposition and its
orbit under `z` give residue degree `f`.  Since `z^4=1`, only

\[
 \boxed{(e,f)=(2,1),(2,2),(2,4)}                \tag{2.2}
\]

occur.  Write `q_f` for the number of rows of residue degree `f`, put

\[
 q=q_1+q_2+q_4,
 \qquad R=q_1+2q_2+4q_4.                        \tag{2.3}
\]

Here `q` is the number of distinct source-boundary divisors and `R` is the
number of meridian transpositions, equivalently the total ramified residue
degree.  If `u` is the fixed-sheet count, then

\[
 d=u+2R.                                         \tag{2.4}
\]

## 3. Complete fixed-sheet spectrum

The positive affine-sheet remainder requires `u>0`.  After removing degrees
one and five below the F2 floor, the complete simple-inertia spectrum is

\[
 \boxed{d\in\{6,10,12,15,20,24,30,40,60,120\}.} \tag{3.1}
\]

There are `13`, rather than ten, action classes: degrees `24` and `40` each
have two subgroup classes with the same peripheral signature, while degree
`30` has two different signatures.

\[
\begin{array}{c|c|c|c|c|c|c|c}
d&\#&u&q_1&q_2&q_4&q&R\\ \hline
6&1&2&2&0&0&2&2\\
10&1&2&4&0&0&4&4\\
12&1&4&0&2&0&2&4\\
15&1&3&6&0&0&6&6\\
20&1&4&0&4&0&4&8\\
24&2&4&0&1&2&3&10\\
30&1&4&1&6&0&7&13\\
30&1&2&14&0&0&14&14\\
40&2&4&0&1&4&5&18\\
60&1&4&0&14&0&14&28\\
120&1&4&0&1&14&15&58
\end{array}                                                    \tag{3.2}
\]

The rows `d=6,10,15` and the second `d=30` row have trivial central action
and recover the fixed-sheet part of the `A_5` atlas.  All other rows are
genuine normalization-gluing packets invisible to determinant and generic
branch Smith data.

## 4. Uniform logarithmic-Chern sieve

Let the `q` source rows have squares `E_i^2=-n_i` and put

\[
 N=\sum_{i=1}^q n_i.
\]

Suppose they traverse the same `b` smooth carrier centers, and let `s_X`
count the further smooth-boundary blowups in the completed source.  The
divisorial term is weighted by the total residue degree `R`, while the
self-intersection term is summed over the actual `q` divisors.  After
subtracting the multiplicity-three cusp lower charge `2R`, the doubled
residuals are

\[
 \boxed{
 2\ell_{\rm rest}^{\rm sf}
 =7d-62+4N-4R(b-6)-s_X\ge0,}                    \tag{4.1}
\]

and

\[
 \boxed{
 2\ell_{\rm rest}^{\rm dbl}
 =7d-67+4N-4R(b-6)-s_X\ge0.}                   \tag{4.2}
\]

Their parity gates remain
`s_X \equiv d \pmod 2` and `s_X \equiv d+1 \pmod 2`, respectively.  The
conditional source-component
floors are now `27+q` and `48+q`, not `27+R` and `48+R`: the latter would
count normalization embeddings rather than divisors.

At minimal negativity `N=q`, maximal contact `b=8`, and `s_X=0`, (4.1)
has the following values:

\[
\begin{array}{c|rrrrrrrrrrr}
d&6&10&12&15&20&24&30_a&30_b&40&60&120\\ \hline
2\ell_{\rm rest}^{\rm sf}
&-28&-8&-2&19&30&38&72&92&94&190&374.
\end{array}                                                    \tag{4.3}
\]

The double-row values are lower by five.  Therefore minimal maximal-contact
packets are excluded in squarefree degrees `6,10,12` and in double-row degree
`12`; every higher row survives even this strongest minimal specialization.
Extra negativity can also repair the three negative squarefree values.

This closes the full simple-inertia group theory but sharply identifies the
remaining geometric gap.  To exclude the E8 endpoint one must now do at
least one of the following:

1. prove that its meridian inertia is simple and rule out every surviving
   row of (3.2) by the global boundary filtration;
2. bound the stable divisorial charge replacing the model-dependent raw
   negativity `N`;
3. constrain the actually attained contact `b`; or
4. treat meridian cycles of length greater than two.

More Laurent coefficients confined to the terminal/carrier packets do none
of these.

The subsequent
[`complete-chain budget`](F2_AFFINE_K1_COMPLETE_CHAIN_BUDGET.md) removes the
raw-negativity gap.  Once the full Cartier determinant cycle and all its
node matching are retained, the point budget is exactly `u-1`; every row of
(3.2) has E8 cusp lower `2R>u-1`.  Thus all simple-inertia rows are excluded
for an effective cyclic completion.  A survivor must instead supply the
explicitly tabulated negative normalization/`Fitt_1` class at the unresolved
global attachment.

<!-- status-consumer: PF2K1CB1 5cc386dba344a867 -->

Finally, isolated `Fitt_1` corrections cannot supply that negative class:
the cyclic-submodule positivity theorem proves that they are effective
quotients of the complete Cartier packet.  This excludes all `13`
one-component simple-inertia actions in the atlas.

<!-- status-consumer: LCSP1 8658eebeb1d65671 -->

## Reproduction

```bash
.venv/bin/python scripts/verify_f2_affine_k1_e8_orbifold_atlas.py
```

The command requires SymPy for exact Todd--Coxeter enumeration.  It takes
about forty seconds on the reference workstation.  It checks the group
order and central order, all `30` subgroup classes, every peripheral orbit
in (3.2), and every value in (4.3).
