# Why the components use only twos and threes

## Contact-order semigroup

The reconstruction theorem fails exactly at repeated inverse roots.  Thus an
omitted value has local root multiplicities in

\[
S_2=\{2,3,4,\ldots\}.
\]

Adjoin zero and regard `S_2` as an additive numerical semigroup.  Collision
adds contact orders.  Its indecomposable nonzero elements are exactly

\[
\operatorname{Atoms}(S_2)=\{2,3\}.
\]

Indeed, neither two nor three is a sum of two allowed contact orders, whereas

\[
m=2+(m-2)\qquad(m\ge4).
\]

Consequently every full-contact partition refines to a partition using only
twos and threes.  The collision-closure theorem turns this elementary
semigroup statement into the component classification.

## Contact-atom theorem

For the normalized degree-`n` exceptional seed locus in characteristic zero:

1. primitive contact types are the factorizations

   \[
   n=2a+3b;
   \]

2. their strata have dimension

   \[
   a+b-1;
   \]

3. every contact order at least four is a collision boundary of a refinement
   by these primitive types.

Thus twos and threes appear before Mason--Stothers is used.  They are selected
by the local reconstruction threshold and the additive collision law.

Geometrically, multiplicity two is ordinary tangency or simple ramification,
while multiplicity three is the first unsplittable higher contact, represented
on the dual discriminant by a cusp.  Higher contacts are degenerations inside
the no-simple-root locus:

\[
4\rightsquigarrow2+2,\qquad
5\rightsquigarrow3+2,\qquad
6\rightsquigarrow2+2+2\ \text{or}\ 3+3.
\]

## Dimension optimization

For a full-contact partition `lambda`,

\[
\dim E_\lambda=\ell(\lambda)-1.
\]

Splitting a part increases the support length and hence the stratum dimension.
The unsplittable partitions therefore index maximal closures even though the
resulting components need not be equidimensional.

For `lambda=2^a3^b`,

\[
\dim C_{a,b}=a+b-1={n-b\over2}-1.
\]

The unique top-dimensional component minimizes `b`: it is all-double for even
`n`, and has one triple part with all remaining parts double for odd `n`.
Here the triple contact is the smallest permitted parity correction.

## Omitted-value phase diagram

The unique omitted-value theorem identifies the nonsurjective seed locus
`N_n` with the disjoint union of the full-contact strata.  Maximizing
`a+b-1` subject to `2a+3b=n` gives

\[
 \boxed{\dim N_n=\left\lfloor\frac n2\right\rfloor-1}.          \tag{1}
\]

Since the normalized admissible degree-`n` seed space `A_n` has dimension
`n-3`, subtraction yields

\[
 \boxed{\operatorname{codim}_{A_n}N_n
 =n-3-\left(\left\lfloor\frac n2\right\rfloor-1\right)
 =\left\lfloor\frac{n-3}{2}\right\rfloor}.                     \tag{2}
\]

Thus degree four is the transition case: `dim A_4=dim N_4=1`, so the
nonsurjective locus is dense and a generic normalized quartic seed is
exceptional.  For every `n>=5`, (2) is positive; the nonsurjective locus lies
in a proper closed subset, and a generic normalized degree-`n` seed is
surjective.  In short,

\[
 \boxed{n=4:\ \text{generic nonsurjectivity};\qquad
        n\ge5:\ \text{generic surjectivity}.}                  \tag{3}
\]

## What Mason--Stothers does

Mason does not select the atoms.  It prevents different contact types from
coexisting away from collision.

Suppose two exact degree-`n` full-contact polynomials satisfy

\[
P-Q=L,\qquad \deg L\le1.
\]

Away from collision, `P,Q,L` are pairwise coprime.  Polynomial
Mason--Stothers gives

\[
n\le\ell(\lambda)+\ell(\mu).
\]

Define the excess above double contact by

\[
\epsilon(\lambda)
=\sum_i(\lambda_i-2)
=n-2\ell(\lambda).
\]

Then

\[
\ell(\lambda)+\ell(\mu)
=n-{\epsilon(\lambda)+\epsilon(\mu)\over2}.
\]

The excess is nonnegative.  The unique zero-excess partition is the all-double
partition, so two distinct partitions have positive total excess.  Hence

\[
\ell(\lambda)+\ell(\mu)<n,
\]

contradicting Mason.  Therefore:

> Two distinct exact full-contact multiplicity types never support distinct
> omitted values away from collision.

This strengthens the maximal-component separation theorem: it applies to
every pair of distinct full-contact partitions.

## Threshold-`r` generalization

The atom mechanism is universal, but the particular atoms two and three are
not.  Suppose a different reconstruction problem allowed only multiplicities
at least `r`.  Its additive contact semigroup would be

\[
S_r=\{r,r+1,r+2,\ldots\},
\]

with

\[
\operatorname{Atoms}(S_r)=\{r,r+1,\ldots,2r-1\}.
\]

Every `m>=2r` splits as `r+(m-r)`, while no smaller allowed multiplicity can
split into two parts at least `r`.

ABC separation depends only on radical support.  Two threshold-`r` contact
polynomials satisfy

\[
\ell(\lambda)+\ell(\mu)\le {2n\over r}.
\]

For `r>=3` this is automatically less than `n`.  The present `r=2` theory is
the critical borderline case: equality in the support bound is possible only
for the unique all-double type, and the excess calculation separates every
distinct pair.

The general structural principle is therefore:

> The atoms of the allowed contact semigroup label the primitive collision
> types; the dimension function selects the dominant ones; and abc rigidity
> controls whether different types can coexist away from collision.

Turning this principle into a component theorem for another reconstruction
problem would still require the geometric inputs proved here: nonempty
incidences, collision deformation, irreducibility, and finite recovery.

## Executable certificate

Run:

```bash
python scripts/verify_contact_atom_principle.py
```

The script verifies the atoms of `S_r` for thresholds two through eight,
constructs atomic refinements, checks the excess/Mason identity for every
distinct full-contact pair through degree twenty-four, and audits the strict
threshold-`r` support bound.
