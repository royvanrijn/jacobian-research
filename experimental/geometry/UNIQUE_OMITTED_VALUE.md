# Unique omitted-value theorem

Fix a normalized admissible seed of inverse degree `n>=3` over an
algebraically closed field of characteristic zero.  Omitted value below means
an omitted value `(s,t)` of its normalized inverse pencil

\[
H(W)-sW+t.
\]

## Theorem

Every normalized admissible seed has at most one omitted inverse-pencil value.

Consequently, if `E_lambda` denotes the exact full-contact stratum of
multiplicity type `lambda`, then

\[
\boxed{
\mathcal N_n=
\bigsqcup_{\substack{\lambda\vdash n\\\lambda_i\ge2}}
\mathcal E_\lambda.}
\]

Thus the contact strata form a genuine disjoint stratification of the
nonsurjective locus, not merely a covering by incidence images.

## Proof

Suppose one seed omitted two distinct pencil values.  After division by their
common leading coefficient, the two inverse polynomials become monic
full-contact polynomials `P,Q` satisfying

\[
P-Q=L,
\qquad 0\ne L,
\qquad \deg L\le1.
\]

Let their exact multiplicity partitions be `lambda,mu`.  The three
polynomials `P,Q,L` are pairwise coprime.  A common root of `P` and `Q` would
have order at least two in their difference, which is impossible for a
nonzero affine polynomial.  The remaining pairwise gcd statements follow
immediately.

If `L` is linear, polynomial Mason--Stothers gives

\[
n\le\ell(\lambda)+\ell(\mu).
\]

Every partition part is at least two, so both support lengths are at most
`n/2`.  Equality is possible only when `n` is even and

\[
\lambda=\mu=2^{n/2}.
\]

If `L` is constant, its radical contributes degree zero rather than one, and
Mason gives the strictly stronger inequality

\[
n\le\ell(\lambda)+\ell(\mu)-1,
\]

which is impossible.  It remains only to eliminate the linear all-double
case.

Write

\[
P=A^2,
\qquad Q=B^2,
\]

where `A,B` are monic of degree `m=n/2>=2`.  Since the omitted values are
distinct, `A!=B`, and

\[
L=P-Q=(A-B)(A+B).
\]

In characteristic zero, `A+B` has degree `m`, while the nonzero polynomial
`A-B` has degree at least zero.  Therefore

\[
\deg L=\deg(A-B)+m\ge m\ge2,
\]

contradicting `deg L<=1`.  Two distinct omitted values cannot exist.

Finally, one omitted polynomial has one exact multiplicity partition by
unique factorization.  Hence a seed belongs to exactly one `E_lambda`, proving
the disjoint-union formula.

## Relation to the component theorem

The distinction between strata and components is now clean:

- `N_n` is the disjoint union of the exact strata `E_lambda`;
- closures introduce the collision order;
- maximal 2/3 partitions index the irreducible components of `closure(N_n)`;
- component intersections consist of the same unique omitted polynomial
  acquiring a common coarser multiplicity partition.

Thus intersections of component closures never represent two different
omitted pencil values.  They represent different atomic descriptions
colliding to one exact polynomial.

For the three-dimensional weighted map, a single omitted pencil value may
lift to a positive-dimensional locus of target coordinates because of the
weighted target parameterization.  The theorem concerns the normalized
inverse-pencil value `(s,t)`.

## Executable certificate

Run:

```bash
python scripts/verify_unique_omitted_value.py
```

The script checks all Mason support alternatives through degree twenty-four
and realizes every possible degree of `A-B` in the all-double square case.
