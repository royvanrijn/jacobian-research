# Degree-ten four-channel coherent frontier

## Status

This note records two results of different strength.  The direction-collision
boundary is excluded exactly in characteristic zero.  The pairwise-distinct
chart has only an exact finite-field computation and is **not** promoted to a
characteristic-zero obstruction.  In balanced degree ten, choose four of the
five positive even harmonic summands and put

\[
 P=\sum_{i=0}^3 a_i\rho^{(10-\ell_i)/2}
       \langle v_i,z\rangle^{\ell_i},
 \qquad
 \{\ell_i\}\subset\{2,4,6,8,10\}.
\]

For the four genuinely new profiles containing \(\mathcal H_{10}\), every
nonterminal direction-collision coefficient-torus chart has an exact
characteristic-zero unit certificate.  There are thirteen set partitions per
profile and hence 52 new charts.  Together with the radially equivalent
\((2,4,6,8)\) collision theorem in
[GVC3_FOUR_COHERENT_CHANNEL_FRONTIER.md](GVC3_FOUR_COHERENT_CHANNEL_FRONTIER.md),
this closes the direction-collision boundary for all five four-of-five
profiles.

On the chart where the four isotropic direction points are pairwise
distinct and every coefficient is nonzero, the invariant moment ideals
through \(\mu_7\) are the unit ideal over each of
\(\mathbb F_{101},\mathbb F_{103},\mathbb F_{107}\), for all five
four-of-five profiles.  Thus none of these charts supplies a finite-field
pure-moment survivor at the declared primes.  This does **not** exclude a
characteristic-zero point, a direction-collision stratum, or the
five-channel profile.  The exact five-channel rational-slice result recorded
below is a separate boundary theorem and does not promote this chart.

## Exact direction-collision certificate

For each of

\[
 (2,4,6,10),\quad(2,4,8,10),\quad(2,6,8,10),\quad(4,6,8,10),
\]

enumerate all set partitions of the four channel labels into two or three
direction blocks.  The all-singleton partition is the pairwise-distinct chart,
and the one-block partition is terminal by one-sided phase weight.  After
normalizing \(a_0=1\), exact elimination over \(\mathbb Q\) gives

\[
 \bigl(\mu_2,\ldots,\mu_c,
        z a_1a_2a_3-1\bigr)=(1)
\]

for every one of the remaining 52 charts.  Thus the localized moment ideal is
the unit ideal, not merely radical or zero-dimensional.  Modular quotient
saturation at \(101,103,107\) discovers and replays the cutoffs; msolve then
returns the literal reduced basis \([1]\) in characteristic zero.

| last required moment \(c\) | charts |
|---:|---:|
| 3 | 21 |
| 4 | 3 |
| 5 | 5 |
| 6 | 22 |
| 9 | 1 |

The unique depth-nine chart is profile \((2,4,8,10)\) with collision blocks
\(\{0,2\}\mid\{1,3\}\).  Only moments divisible by three are nonzero on
this phase configuration, so the decisive equation is \(\mu_9\).  This is why
compilation through moment eight would have produced a false survivor.  A
second arithmetic subtlety occurs for \((2,6,8,10)\) with blocks
\(\{0,3\}\mid\{1,2\}\): the exact rational cutoff is five and primes 101 and
103 agree, whereas characteristic 107 needs moment six.  The rational
certificate, rather than unanimity of modular cutoffs, determines the theorem.

A zero coefficient reduces to the exact two/three-channel theorem GVC3IHC.
When every direction coincides, the form has one-sided phase support, so every
fixed multiplier channel eventually vanishes.  These facts complete the
coefficient and collision boundary; they do not address the pairwise-distinct
coefficient torus.

## Moment system

Normalize the direction points to

\[
 \infty,\quad 0,\quad 1,\quad \lambda,
 \qquad \lambda(\lambda-1)\ne0,
\]

and set \(a_0=1\).  The Wick/occupation compiler is the same as in
[GVC3_FOUR_COHERENT_CHANNEL_FRONTIER.md](GVC3_FOUR_COHERENT_CHANNEL_FRONTIER.md).
In particular, its finite proper-hypergeometric formula supplies the
declared moment-sequence promotion gate independently of the modular
elimination.

Radial powers do not change the pure spherical moment equations: on the
unit quadric, only the selected harmonic degrees remain.  Consequently the
\((2,4,6,8)\) profile has the same pure moment system in balanced degrees
eight and ten.

The compiled term counts are:

| harmonic degrees | \(\mu_3\) | \(\mu_4\) | \(\mu_5\) | \(\mu_6\) | \(\mu_7\) |
|---|---:|---:|---:|---:|---:|
| \(2,4,6,8\) | 27 | 53 | 112 | 260 | 417 |
| \(2,4,6,10\) | 15 | 79 | 101 | 236 | 428 |
| \(2,4,8,10\) | 34 | 46 | 129 | 295 | 441 |
| \(2,6,8,10\) | 31 | 62 | 129 | 276 | 489 |
| \(4,6,8,10\) | 31 | 60 | 178 | 339 | 509 |

In every row \(\mu_2=0\).  Saturating the compact ideal

\[
 (\mu_3,\mu_4,\mu_5,\mu_6,\mu_7)
 :\bigl(a_1a_2a_3\lambda(\lambda-1)\bigr)^\infty
\]

by successive ideal quotients gives the literal basis \([1]\) at all
fifteen profile/prime pairs.  The saturation exponent is one in every
case.  This quotient representation is materially smaller than adjoining
a Rabinowitsch variable: the corresponding raw and linear-pivot
Rabinowitsch runs reached their declared time bounds without results.

## What is and is not learned

The calculation gives a coherent modular pattern across every
pairwise-distinct four-channel profile for \(\Delta^5\).  No multiplier
channel is tested, because no pure-moment survivor reaches the promotion
stage.  The finite proper-hypergeometric formula is recorded so that a
future survivor would already satisfy the sequence gate.

The remaining exact tasks are:

1. promote the degree-eight \((2,4,6,8)\) system by a homogeneous rational
   certificate;
2. promote the pairwise-distinct degree-ten units over \(\mathbb Q\); and
3. attack the degree-ten five-channel chart in its two cross-ratio
   parameters.

Coefficient-zero faces and direction collisions are now exact.  Only the
pairwise-distinct coefficient torus remains within each four-channel profile.

### Five-channel pivot and rational-slice frontier

For the full \((2,4,6,8,10)\) profile, normalize five directions to
\(\infty,0,1,\lambda,\mu\).  A compile-only replay gives

\[
 \#\mu_3=271,\qquad \#\mu_4=1142,\qquad \#\mu_5=3686.
\]

There are six chart variables after normalizing \(a_0=1\), so dimension
counting puts the first plausible zero-dimensional cutoff at \(\mu_8\), not
at \(\mu_5\).  The rapid term growth makes a raw Gröbner run the wrong next
step.

Write the first nonzero moment as a polynomial in one chart variable.  Only
\(a_2,a_3,a_4\) are linear pivots for \(\mu_3\); the degrees in
\(a_1,\lambda,\mu\) are respectively \(2,16,20\).  Solving
\(\mu_3=A+B a_i=0\) and clearing the minimal power of \(B\) gives the
following exact term counts for the transformed fourth moment:

| pivot | terms in transformed \(\mu_4\) |
|---|---:|
| \(a_2\) | 29,002 |
| \(a_3\) | 32,148 |
| \(a_4\) | 7,333 |

Thus \(a_4\) is the unique practical pivot on this normalization.  Its
transformed fifth moment already has 58,971 terms, so the main
\(A B\ne0\) localization has not been sent into a raw higher-moment
elimination.  On the exceptional boundary \(B=0\), nonzero \(a_4\) and
\(\mu_3=0\) force \(A=B=0\).  Here \(A\) and \(B\) have only 27 and 244
terms, but generic two-cross-ratio quotient saturation of \((A,B)\) timed
out at 60 seconds at all three discovery primes; the ideal through \(\mu_5\)
also timed out at 101 after 300 seconds.  These are recorded time bounds, not
survivors or characteristic-zero evidence.  They show that the next generic
calculation must stratify the \((A,B)\) incidence before saturation.

There is a canonical nested stratum.  The polynomial \(A\) is the irreducible
four-channel \((2,4,6,8)\) third moment: it is independent of \(\mu\), has
gcd one with \(B\), and is linear in both \(a_2\) and \(a_3\).  Solving
\(A=0\) on the two possible main charts gives:

| nested pivot | transformed \(B\) | transformed \(\mu_4\) | transformed \(\mu_5\) |
|---|---:|---:|---:|
| \(a_2\) | 958 | 6,173 | 34,782 |
| \(a_3\) | 648 | 7,340 | 15,300 |

The \(a_3\) chart is therefore the declared continuation.  Its first
transformed equation is much smaller, but quotient saturation by the complete
localization product still timed out after 60 seconds at each of
101, 103, and 107.  The next implementation should saturate successively by
the coefficient, configuration-discriminant, and nested-pivot factors, and
record the first factor that causes growth.  The successive calculation below
does this.  A longer run with the same single product would add no structural
information.  The nested denominator-zero boundary must be treated as a
separate projective stratum.

#### Successive saturation and the cutoff-four component cycle

Rechart by \(a_1=1\), retain the \(a_4\)-boundary and nested \(a_3\)-pivot,
and factor the complete localization denominator over \(\mathbb Q\).  At
cutoff three, factor-by-factor saturation completes at 101, 103, and 107.  In
the original \(a_0=1\) chart only the factor \(a_1\) is removed, with
valuation one; after recharting by \(a_1=1\), no localization factor divides
the sole cutoff-three generator.  The remaining nine factors are certified
coprime by exact polynomial division or a principal gcd.  Thus the earlier
cutoff-three timeout was algorithmic, not a survivor.

Cutoff four is different.  The two transformed moment equations form a
codimension-two complete intersection in
\(\mathbb F_{101}[a_0,a_2,a_4,\lambda,\mu]\).  On \(a_0=0\), their common
boundary divisor has degree 28.  The first exact colon generator has 10,132
terms and total degree 67.  Its residual boundary gcd has degree 17 and is
exactly the restriction of the final irreducible nested-denominator factor
\(D\).  Writing the three current generators as

\[
 P_i=a_0 A_i+D_0B_i,\qquad D_0=D|_{a_0=0},
\]

gives certified residual-intersection minors
\(A_iB_j-A_jB_i\): multiplication by both \(a_0\) and \(D_0\) returns an
explicit combination of the \(P_i\).  Two sparse minors, with 16,661 and
19,039 terms modulo 101, reduce the \(a_0\)-boundary gcd to one.  Hence the
\(a_0\)-factor is fully saturated at exponent two in this modular chart.

The same unmixed-boundary test then certifies unit boundary gcds for

\[
 a_2,\quad a_4,\quad \lambda,\quad \mu,\quad \lambda-1,\quad \mu-1.
\]

Together with \(a_0\), seven of the ten irreducible localization factors are
therefore scheme-theoretically complete at all three declared primes at
cutoff four.  The two remaining linear
factors

\[
 143a_0+60a_2,\qquad \lambda-\mu
\]

both have degree-17 boundary gcd equal to the corresponding restriction of
\(D\).  The three minors from the original three-generator submatrix were
adjoined exactly on each boundary (term counts 22,635, 25,260, 47,943 and
26,192, 31,190, 129,969), but the boundary gcd remains degree 17.  The final
\(D\)-hypersurface dimension test and direct saturation each exceed the
declared five-second diagnostic bound.  The resulting eleven-generator
record is a bounded modular component analysis, not a unit certificate and
not characteristic-zero evidence.

The reciprocal calculation now closes the *support* of this cycle at all
three declared primes.  The factor \(D\) is linear in \(a_0\).  Write
\(D=C a_0+R\), solve \(a_0=-R/C\), and for a polynomial of \(a_0\)-degree
\(n\) use the cleared restriction \(C^nP(-R/C)\).  The coefficient
\(C=38a_2(\lambda-1)^{12}\) is already invertible on this chart.  Singular
verifies, for every one of the five current generators, the exact identity

\[
 C^nP-C^nP(-R/C)\in(D)
\]

by a zero remainder.  The gcd of the five cleared restrictions has total
degree 29.  After stripping all seven completed localization factors it has
precisely the valuations

\[
 v_{143a_0+60a_2}=1,
 \qquad
 v_{\lambda-\mu}=12,
\]

and unit residual gcd.  The individual valuations are respectively

\[
 (1,2,1,1,1),
 \qquad
 (12,20,12,12,12).
\]

The degrees, valuations, and five zero-remainder identities agree over
\(\mathbb F_{101},\mathbb F_{103},\mathbb F_{107}\).  Thus every divisorial
piece of the cutoff-four ideal lying on \(D=0\) is supported on the two known
linear boundaries, and the localization dependency graph is set-theoretically
closed for all ten factors.  This does **not** construct the saturated ideal:
the \((D,\lambda-\mu)\) piece is visibly thick of order 12.  A preliminary
restricted-minor shortcut conflicted with these verified valuations and was
therefore retired rather than promoted.

The same normal calculation now promotes generically over characteristic
zero, without a Gröbner basis over a rational-function field.  For the two
original transformed equations, containing 648 and 7,340 terms, the exact
cleared restrictions to \(D=0\) contain 3,284 and 7,581 terms and have
\((\lambda-\mu)\)-valuations 12 and 20.  More explicitly,

\[
 C=38a_2(\lambda-1)^{12},\qquad
 R=\lambda^6\bigl(48a_2(\lambda-1)^{10}+19\lambda^{10}\bigr),
\]

and, up to a nonzero rational scalar, their characteristic-zero gcd is

\[
 R\,\bigl(60a_2C-143R\bigr)(\lambda-\mu)^{12}.
\]

The first two factors are precisely the cleared restrictions of \(a_0\) and
\(143a_0+60a_2\).  The gcd has total degree 46 and 728 terms; division by
these two factors and \((\lambda-\mu)^{12}\) leaves a constant.  Thus the
support calculation itself is exact over \(\mathbb Q\), while the earlier
degree-29 five-generator calculation remains an independent three-prime
replay after the seven completed saturation steps.

There is also an exact generic multiplicity certificate.  Work over

\[
 K=\mathbb Q(a_2,a_4,\lambda),\qquad d=D,\qquad e=\lambda-\mu.
\]

The two normal derivatives have 25 and 394 terms.  At order \(e^{12}\), the
restriction coefficients have 152 and zero terms respectively, because the
second restriction starts at order 20.  The resulting two-by-two normal-jet
determinant is a nonzero polynomial with 959 terms and total degree 154.
Exact division removes \(C^3\), the square of the cleared \(a_0\)
restriction, and the square of the cleared \(143a_0+60a_2\) restriction.
The only remaining exceptional factor is a 115-term polynomial \(H\) of
total degree 47.  Thus, on the declared localization chart, the determinant
can vanish only on \(H=0\).  Consequently one equation is a regular
\(d\)-parameter and the other has first nonzero image of order \(e^{12}\);
in the completed local ring at the generic \((D,\lambda-\mu)\) point the
ideal is equivalent to

\[
 (d,e^{12}),
\]

and has length exactly 12.  This closes the generic normal profile in
characteristic zero.

The exceptional branch is exact as well.  The 115-term factor \(H\) is
irreducible and squarefree over \(\mathbb Q\).  Modulo \(H\), normal jets one
through five vanish and jet six is nonzero, while the lower collision initial
and lower normal derivative remain nonzero.  Hence the generic completed
ideal on \(H=0\) has length \(12+6=18\).  The next possible exceptional
factor reduces, after removing chart units, to an irreducible 154-term
polynomial \(K\) of total degree 57.  Its resultant with \(H\) in \(a_2\)
is

\[
 \lambda^{120}(\lambda-1)^{120}W(\lambda),
\]

where \(W\) is irreducible of degree 32.  This same \(W\) is the projection
factor obtained from the already inverted boundary
\(143a_0+60a_2=0\).  Over \(\mathbb Q[\lambda]/(W)\), the gcd of \(H\) and
\(K\) is linear in \(a_2\), and the companion boundary reduces to zero
modulo that gcd.  Thus every putative deeper \(H=K=0\) point lies outside the
declared chart.  The complete \(D\)-supported normal stratification at cutoff
four therefore has generic lengths 12 and 18, with no deeper chart stratum.
This still does not construct the ambient cutoff-four colon.

One rational parameter fiber gives a literal standard-basis replay.  Put

\[
 (a_2,a_4,\lambda)=(2,3,5),\qquad
 d=D,\qquad e=\lambda-\mu,
\]

and regard \((d,e)\) as local coordinates.  The two exact transformed
equations, originally containing 648 and 7,340 terms, specialize and combine
to bivariate polynomials with 40 and 82 terms.  A characteristic-zero local
standard basis is literally

\[
 (d,e^{12}),
\]

so the quotient length is exactly 12.  The modular local bases on the same
fiber have length 12 at 101, 103, and 107 as well.

The fifth moment removes the need for an ambient colon on the \(D=0\) fibre.
After the nested pivot it has 15,300 terms, and its exact cleared restriction
to \(D=0\) has 38,015 terms and total degree 152.  The common gcd of all
three restrictions now has only 56 terms and total degree 34; up to a
nonzero rational scalar it is

\[
 R\bigl(60a_2C-143R\bigr),
\]

with one copy of each factor and no residual divisor.  In particular the
factor \((\lambda-\mu)^{12}\) has disappeared.  After all already inverted
chart factors are removed from the fifth restriction, its specialization to
\(\mu=\lambda\) is, up to scalar,

\[
 a_2^3\lambda^8(\lambda-1)^{36},
\]

which is a unit on the chart.  The second restriction is a chart unit times
\((\lambda-\mu)^{20}\).  The finite geometric-series identity for
\(u+(\lambda-\mu)g\) modulo \((\lambda-\mu)^{20}\) therefore gives the unit
ideal on the localized \(D=0\) fibre.  Thus moment five closes the entire
\(D\)-supported branch, including \(H=0\), in characteristic zero.

For discovery, changing coordinates from \(a_0\) to \(d=D\) makes the
normal saturation linear.  The corresponding residual-intersection minor
has 64,467 terms and total degree 132 modulo 101; its boundary gcd is again
exactly the two transformed chart factors above.  Full coordinate-colon and
component-quotient runs still time out at 600 seconds, so they are retained
as algorithmic diagnostics rather than certificates.  The remaining target
is now the residual \(D\ne0\) chart, not an exceptional normal branch.

One rational fiber can already be closed exactly.  Set
\((\lambda,\mu)=(2,3)\), which retains five distinct directions, and compile
with the numeric Gram matrix before expanding.  The specialized moment term
counts through order eight are

\[
 0,9,22,42,79,127,202.
\]

For

\[
 I_c=(A,B,\mu_4,\ldots,\mu_c):(a_1a_2a_3a_4)^\infty,
\]

the dimensions modulo 101 at \(c=3,4,5\) are \(2,1,0\); at \(c=5\) the
quotient length is 36.  The next moment makes \(I_6=(1)\) modulo
101, 103, and 107, always with saturation exponent one.  Exact Rabinowitsch
elimination over \(\mathbb Q\) then returns the literal reduced basis
\([1]\).  Hence the \(a_4\)-pivot boundary contains no pure-moment point on
this rational cross-ratio fiber.  This is an exact slice theorem, but it does
not address the main \(A B\ne0\) localization, the generic \((\lambda,\mu)\)
boundary, or any other cross-ratio fiber.

## Reproduction

```bash
.venv/bin/python \
  scripts/research_gvc3_degree10_distinct_four_channels.py \
  --workers 11 --timeout 900
```

The generated record is
`artifacts/generated-results/gvc3_degree10_distinct_four_channels_modular.json`.

The exact direction-collision boundary is reproduced by

```bash
.venv/bin/python \
  scripts/verify_gvc3_degree10_four_channel_collisions.py \
  --workers 8 --modular-timeout 180 --exact-timeout 300
```

It writes
`artifacts/generated-results/gvc3_degree10_four_channel_collisions.json`.
This replay requires Singular and msolve.  Its three-prime computations are
discovery and consistency checks; the promoted result is the 52 exact
characteristic-zero unit bases.

The preliminary five-channel term counts are reproduced by

```bash
.venv/bin/python \
  scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 5 --compile-only \
  --output \
  artifacts/generated-results/gvc3_degree10_five_channel_compile5.json
```

The linear-pivot scan and the selected main-chart compilation are

```bash
.venv/bin/python scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 4 --compile-only \
  --scan-linear-pivots --scan-pivot-max-order 4 \
  --output \
  artifacts/generated-results/gvc3_degree10_five_channel_pivot_scan4.json
.venv/bin/python scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 5 --compile-only \
  --linear-pivot a4 --output \
  artifacts/generated-results/gvc3_degree10_five_channel_a4_pivot_compile5.json
```

The exact \((\lambda,\mu)=(2,3)\) boundary slice is reproduced by

```bash
.venv/bin/python \
  scripts/research_gvc3_degree10_five_channel_slice.py \
  --lam 2 --mu 3 --max-order 8 \
  --modular-timeout 180 --exact-timeout 300
```

It writes
`artifacts/generated-results/gvc3_degree10_five_channel_a4_boundary_slice.json`.
The characteristic-zero \([1]\) basis at moment six is the certificate; the
three finite fields are discovery and replay only.

The selected nested generic-boundary compilation and its bounded first
saturation attempt are

```bash
.venv/bin/python scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 5 --compile-only \
  --linear-pivot-boundary a4 --boundary-linear-pivot a2 \
  --output \
  artifacts/generated-results/gvc3_degree10_five_channel_a4_boundary_a2_pivot_compile5.json
.venv/bin/python scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 5 --compile-only \
  --linear-pivot-boundary a4 --boundary-linear-pivot a3 \
  --output \
  artifacts/generated-results/gvc3_degree10_five_channel_a4_boundary_a3_pivot_compile5.json
.venv/bin/python scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 3 \
  --linear-pivot-boundary a4 --boundary-linear-pivot a3 \
  --saturation-method quotient --report-quotient-dimension \
  --primes 101 103 107 --timeout 60 \
  --output \
  artifacts/generated-results/gvc3_degree10_five_channel_a4_boundary_a3_pivot_modular3.json
.venv/bin/python scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 3 \
  --linear-pivot-boundary a4 --boundary-linear-pivot a3 \
  --saturation-method successive --report-quotient-dimension \
  --primes 101 103 107 --timeout 60 \
  --output \
  artifacts/generated-results/gvc3_degree10_five_channel_a4_boundary_a3_pivot_successive3.json
.venv/bin/python scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 4 \
  --normalize-coefficient a1 \
  --linear-pivot-boundary a4 --boundary-linear-pivot a3 \
  --saturation-method successive --report-quotient-dimension \
  --primes 101 --timeout 5 \
  --output \
  artifacts/generated-results/gvc3_degree10_five_channel_a1norm_a4_boundary_a3_pivot_successive4_p101.json
.venv/bin/python scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 4 \
  --normalize-coefficient a1 \
  --linear-pivot-boundary a4 --boundary-linear-pivot a3 \
  --saturation-method successive --report-quotient-dimension \
  --primes 101 103 107 --timeout 90 \
  --output \
  artifacts/generated-results/gvc3_degree10_five_channel_a1norm_a4_boundary_a3_pivot_component_cycle_modular3.json
```

The quotient record contains three timeouts and makes no emptiness or survivor
claim.  The first successive record completes every factor at cutoff three.
The first cutoff-four record is the historical seven-of-ten component
certificate; its final status is deliberately `timeout`, not `completed`.
The second replays the reciprocal solved-\(D\) support classification.  Its
three modular results have status `component_support_classified`, with
`factors_component_classified=10`, `support_closed=1`, and
`scheme_closed=0`.  These fields distinguish the closed support graph from
the still-open ambient saturated ideal.  The same command records
`exact_generic_primary_result`: over `Q(a2,a4,lam)` the two exact
restrictions have valuations 12 and 20, their support gcd has valuations one,
one, and twelve on the cleared `a0`, `143*a0+60*a2`, and `lam-mu` factors,
and the nonzero normal-jet determinant proves generic normal length 12.  It
also records `exact_local_primary_result`: on the rational normal fiber
\((a_2,a_4,\lambda)=(2,3,5)\), its characteristic-zero basis is exactly
`dd,ee^12` and its local length is 12.

The pre-nested generic-boundary time bounds are reproduced by

```bash
.venv/bin/python scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 3 \
  --linear-pivot-boundary a4 --saturation-method quotient \
  --primes 101 103 107 --timeout 60 --output \
  artifacts/generated-results/gvc3_degree10_five_channel_a4_boundary_modular3.json
.venv/bin/python scripts/research_gvc3_many_coherent_channels.py \
  --degrees 2,4,6,8,10 --max-order 5 \
  --linear-pivot-boundary a4 --saturation-method quotient \
  --primes 101 --timeout 300 --output \
  artifacts/generated-results/gvc3_degree10_five_channel_a4_boundary_modular5.json
```

Both records contain only timeouts.
