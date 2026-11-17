# Unique omitted-value lemma

Fix a normalized admissible seed of inverse degree `n>=3` over an
algebraically closed field of characteristic zero.  By the
[omitted-value classifier](OMITTED_VALUE_CLASSIFICATION.md), an inverse-pencil
value `(s,t)` is omitted exactly when every root of

\[
H(W)-sW+t
\]

has multiplicity at least two.

## Foundational corollary

Every normalized admissible seed has at most one omitted inverse-pencil
value.

Equivalently, two distinct members of one affine inverse pencil cannot both
have no simple root.

## Proof

Suppose two distinct omitted values give

\[
P(W)=H(W)-sW+t,
\qquad
Q(W)=H(W)-s'W+t'.
\]

Then every root of `P` and `Q` has multiplicity at least two, while

\[
L(W)=P(W)-Q(W)=(s'-s)W+(t-t')
\]

is a nonzero polynomial of degree at most one.

The polynomials `P,Q,L` are pairwise coprime.  Indeed, if `r` were a common
root of `P` and `Q`, then

\[
P(r)=Q(r)=P'(r)=Q'(r)=0.
\]

Thus `L(r)=L'(r)=0`, which forces the affine polynomial `L` to vanish
identically, contrary to `P!=Q`.  Coprimality with `L=P-Q` follows as well.

Let `ell(P)` and `ell(Q)` denote the numbers of distinct roots.  If `L` is
linear, polynomial Mason--Stothers applied to `Q+L=P` gives

\[
n\leq \ell(P)+\ell(Q).
\]

If `L` is constant, its radical contributes degree zero and the stronger
bound is

\[
n\leq \ell(P)+\ell(Q)-1,
\]

which is impossible.  Since every root multiplicity of `P` and `Q` is at
least two,

\[
\ell(P)\leq {n\over2},
\qquad
\ell(Q)\leq {n\over2}.
\]

The linear case therefore forces equality throughout.  In particular, `n`
is even and both polynomials are all-double.  Since they have the same
leading coefficient `c`, write

\[
P=cA^2,
\qquad
Q=cB^2,
\]

with `A,B` monic of degree `m=n/2>=2`.  Then

\[
L=c(A-B)(A+B).
\]

If `A!=B`, characteristic zero gives `deg(A+B)=m`, so

\[
\deg L=\deg(A-B)+m\geq m\geq2,
\]

contrary to `deg L<=1`.  Hence `A=B`, so `P=Q`, again a contradiction.
Therefore at most one omitted pencil value exists.

## Structural consequences

If `E_lambda` denotes the exact full-contact stratum of multiplicity type
`lambda`, unique factorization assigns the unique omitted polynomial one
exact multiplicity partition.  Hence the nonsurjective locus is the genuine
disjoint stratification

\[
\boxed{
\mathcal N_n=
\bigsqcup_{\substack{\lambda\vdash n\\\lambda_i\geq2}}
\mathcal E_\lambda.}
\]

In particular:

- closure-order limits use the same omitted polynomial with roots collided;
- component intersections cannot encode two different omitted values;
- coincident-root factorizations reconstruct one canonical omitted
  polynomial; and
- the omitted-value phase diagram may optimize over disjoint exact contact
  strata.

For the three-dimensional weighted map, this statement concerns the
normalized inverse-pencil value `(s,t)`.  Its lift to weighted target
coordinates can still be positive-dimensional.

## Executable certificate

Run:

```bash
.venv/bin/python scripts/verify_unique_omitted_value.py
```

The checker audits every Mason support alternative through degree twenty-four
and realizes every possible degree of `A-B` in the all-double equality case.
