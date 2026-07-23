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

The [unique omitted-value lemma](UNIQUE_OMITTED_VALUE.md) identifies the
nonsurjective seed locus `N_n` with the disjoint union of the full-contact
strata.  This is the separation input used here: it rules out two distinct
omitted values uniformly, including the all-double borderline case and pairs
having the same multiplicity type.  Maximizing
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

## Optional comparison: the old Mason separation

The following argument is not an additional input to the contact-atom theorem
or the phase diagram.  It is retained as a shorter alternative proof for the
special case of two distinct multiplicity types, and as motivation for the
bounded-difference estimate in the abstract threshold theorem.

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

This is strictly weaker than the unique omitted-value lemma used above.  As
recorded, it is only a distinct-type statement; the same excess estimate also
handles equal non-all-double types, but remains borderline for a pair of
all-double polynomials.  The square-factorization step in the unique
omitted-value proof closes that case.  Thus the Mason calculation here is an
optional type-separation proof, not a second dependency.

## Threshold-`r` generalization

The direct `r=2` argument above uses the contact semigroup, the collision
theorem, the stratum dimension formula, and the unique omitted-value lemma; it
does not invoke a threshold theorem.  The standalone
[threshold contact-atom theorem schema](THRESHOLD_CONTACT_ATOM_SCHEMA.md)
abstracts the same mechanism.  Under explicit reconstruction and collision
hypotheses it proves that the
atoms are `r,...,2r-1`, primitive components are their restricted partitions,
dominant dimensions are a one-row integer program, component counts and
codimensions are quasipolynomial, and Mason separation for bounded-degree
differences is automatic in the precise range

\[
n-2\left\lfloor n/r\right\rfloor-d+1>0.
\]

The old type-separation calculation explains why that abstract Mason bound is
strict for `r>=3` but borderline at `r=2`.  In the actual inverse-pencil
problem, the unique omitted-value lemma supplies the stronger `r=2`
separation independently of the schema.

## Executable certificate

Run:

```bash
python scripts/verify_contact_atom_principle.py
python scripts/verify_unique_omitted_value.py
```

The first script verifies the atoms of `S_r` for thresholds two through eight,
constructs atomic refinements, and audits the dimension optimization.  It
also retains the excess/Mason identity for distinct `r=2` types as an optional
regression and checks the strict threshold-`r` support bound.  The second
script verifies the load-bearing unique omitted-value lemma, including the
all-double equality case.
