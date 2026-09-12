# Four coherent harmonic channels in degree eight

## Status

This note records one exact boundary theorem and one unpromoted modular
calculation for the balanced degree-eight family

\[
 P=\sum_{i=0}^3 a_i\rho^{(8-\ell_i)/2}\langle v_i,z\rangle^{\ell_i},
 \qquad (\ell_0,\ell_1,\ell_2,\ell_3)=(2,4,6,8),
\]

where every \(v_i\) is isotropic.  Coefficient-zero faces reduce to the exact
two/three-channel theorem in
[GVC3_ISOTROPIC_HARMONIC_CHANNEL_OBSTRUCTION.md](GVC3_ISOTROPIC_HARMONIC_CHANNEL_OBSTRUCTION.md).

The exact result proved here is:

> If two direction points coincide, or if the coefficient \(B\) below
> vanishes on the four-distinct-direction chart, the first six pure Reynolds
> moments exclude the coefficient torus.  Hence these strata contain no
> \(\Delta^4\) GVC witness.

There are thirteen nonterminal direction-collision partitions.  Each has a
literal characteristic-zero msolve basis \([1]\), at a cutoff no larger than
six.  The four-distinct \(B=0\) boundary also has a literal basis \([1]\) over
\(\mathbb Q\) through moment six.

On the complementary four-distinct \(B\ne0\) chart, the saturated moment ideal
through moment seven is \([1]\) modulo each of \(101,103,107\).  Several
characteristic-zero attempts reached their declared time bounds without a
certificate.  That chart is therefore **strong modular evidence, not a
characteristic-zero theorem**.

## Invariant moment compiler

Normalize the distinct direction points in \(\mathbb P^1\) to

\[
 \infty,\quad 0,\quad 1,\quad \lambda,
 \qquad \lambda(\lambda-1)\ne0,
\]

and put \(a_0=1\).  If \(g_{ij}\) is the squared bracket of the corresponding
points, define

\[
 W(s)=[u^s]\exp\!\left(\sum_{i<j}g_{ij}u_i u_j\right)
 =\sum_{\deg(e)=s}\prod_{i<j}\frac{g_{ij}^{e_{ij}}}{e_{ij}!}.
\]

For \(|n|=m\), put \(s_i=\ell_i n_i\) and \(S=\sum_i s_i\).  Up to one fixed
nonzero normalization, the pure Reynolds moment is

\[
 \mu_m=
 \sum_{|n|=m}
 \binom{m}{n}\prod_i a_i^{n_i}
 \frac{\prod_i s_i!}{(S+1)!!}\,W(s).
\]

This is a finite proper-hypergeometric sum: shifting any occupation or edge
index changes a summand by a rational function of the summation indices.
Thus the sequence satisfies the requested promotion gate before any
Groebner candidate is considered.  The equivalent recurrence used by the
checker is

\[
 s_iW(s)=\sum_{j\ne i}g_{ij}W(s-e_i-e_j),\qquad W(0)=1.
\]

## The four-distinct chart

The nonzero compiled moments have term counts

\[
 \#\mu_3=27,\quad \#\mu_4=53,\quad \#\mu_5=112,
 \quad \#\mu_6=260,\quad \#\mu_7=417.
\]

The third moment is linear in \(a_2\):

\[
 \mu_3=A+B a_2,
\]

with

\[
 A=19a_1(112a_1a_3\lambda^{16}+2431)
\]

and

\[
 B=4\left(
 1344a_1a_3\lambda^6(\lambda-1)^{10}
 +4845a_1+1064a_3(\lambda-1)^{12}
 \right).
\]

### The exact \(B=0\) boundary

On the coefficient torus, \(B=0\) and \(\mu_3=0\) imply \(A=0\), hence

\[
 a_3=-\frac{2431}{112a_1\lambda^{16}}.
\]

After this substitution, the equations coming from \(B,\mu_4,\mu_5,\mu_6\)
have respectively \(25,51,93,199\) terms in \(a_1,a_2,\lambda\).  Saturating
by \(a_1a_2\lambda(\lambda-1)\) gives \([1]\) modulo all three discovery
primes.  The exact Rabinowitsch ideal over \(\mathbb Q\) has the literal
reduced msolve basis \([1]\).

### The open \(B\ne0\) frontier

Keep the compact equations in \(a_1,a_2,a_3,\lambda\), rather than expanding
the substitution \(a_2=-A/B\), and saturate by

\[
 a_1a_2a_3\lambda(\lambda-1)B.
\]

Through \(\mu_7\), Singular returns the unit basis at all three primes:

| characteristic | basis | saturation exponent | measured checker time |
|---:|:---:|---:|---:|
| 101 | \([1]\) | 1 | 17.4 s |
| 103 | \([1]\) | 1 | 17.9 s |
| 107 | \([1]\) | 1 | 18.4 s |

These are full polynomial-ring saturations, so they exclude algebraic, not
merely prime-field-rational, cross-ratios in those characteristics.  They do
not by themselves prove the characteristic-zero ideal is a unit.  The
exponent-one result says that the open-set product itself belongs to each
finite-field moment ideal; this is the smallest possible target for a future
rational Nullstellensatz lift.

## Direction collisions

The coefficient torus has thirteen nonterminal set partitions of the four
channels: six one-pair/three-direction partitions, three two-pair partitions,
and four triple--single partitions.  Modular discovery at \(101,103,107\) and
exact characteristic-zero replay give:

- every three-direction collision chart is a unit by moment six;
- every two-pair or triple--single chart is a unit by moment four, with two
  triple--single charts already units at moment three;
- the all-coincident partition is the one-direction terminal stratum, killed
  eventually against every fixed multiplier by phase weight.

The exact output and the partition-by-partition cutoffs are stored in
`artifacts/generated-results/gvc3_four_coherent_channels.json`.

## What remains

The only unclosed stratum of this degree-eight four-channel family is the
four-distinct \(B\ne0\) chart over characteristic zero.  The most useful next
calculation is not another finite-field scan; it is one of:

1. reconstruct and exactly verify a rational Nullstellensatz identity from
   modular lifts;
2. find a smaller affine cover on which exact msolve terminates; or
3. derive the parameter eliminant in \(\lambda\) and certify its finitely many
   fibers.

Only after this chart is exact should the result be promoted to a complete
four-channel obstruction.  The next coherent families are then the five
four-of-five degree-ten profiles and the single five-channel degree-ten
profile.

## Reproduction

```bash
.venv/bin/python scripts/verify_gvc3_four_coherent_channels.py
```

The optional command below attempts the unresolved characteristic-zero chart
and deliberately fails rather than silently promoting a timeout:

```bash
.venv/bin/python scripts/verify_gvc3_four_coherent_channels.py --exact-open
```
