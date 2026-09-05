# Small bad primes do not expose the exceptional quotient

The three retained high/low comparisons now have exact local Kummer
calculations at \(S=\{3,5,7,11,13\}\). **The marked generic subgroup fills
the complete local point image at all thirty curve-prime cases, and fills
the product of the five images on every curve.** Twenty-one cases have bad
reduction. The known \(+9,+9,+7\) directions open no additional support
in this dictionary, either individually by prime or jointly.

This is an **incidence diagnostic**, specifically an exclusion of a proposed
local-support explanation. It is not a point-visibility measurement or a
positive rank predictor. The five-prime result is not a full Selmer
computation and does not include 2, the real place, or other bad primes.

The [protocol](BAD_PRIME_PROTOCOL.json), portable
[local arithmetic inputs](../../artifacts/generated-results/elliptic-curves/rank_jump_bad_prime_inputs_v1.json)
and [report](../../artifacts/generated-results/elliptic-curves/rank_jump_bad_prime_support_v1.json)
are independent of the active search. Replay:

    python3 elliptic-curves/rank-jump/bad_prime_support.py check
    sage -python elliptic-curves/rank-jump/bad_prime_support.py verify

The first command checks the exact binary accounting. The second recomputes
all local arithmetic and checks the 25 explicit generic corrections using
a separate local-power interface.

## Exact local computation, with bounded setup

For each curve, only the previously certified independent point subgroup is
used. The first \(m\) points are its marked generic basis. Clearing
denominators by \(x'=d^2x,\ y'=d^3y\) gives an integral cubic and multiplies
each Kummer representative by the rational square \(d^2\).

The computation imports the existing
[LocalSquareclasses implementation](../cas/research_runtime/local_kummer.py)
without a cache adapter. PARI constructs an order maximal at the five
specified primes; global maximality and global discriminant factorization
are not requested. The guarantee for local prime decomposition in this
setup is documented in [PARI's number-field manual](https://pari.math.u-bordeaux.fr/dochtml/html-stable/General_number_fields.html#nfinit).
At each prime ideal, valuation parity and the square character of the
remaining unit give the full odd-adic squareclass.

Local minimal reduction is computed separately with
[elllocalred](https://pari.math.u-bordeaux.fr/dochtml/ref-stable/Elliptic_curves.html#elllocalred).
This avoids treating denominators or a nonminimal discriminant as evidence
of bad reduction. Each sequential curve worker had a 30-second bound and
its own checkpoint. All six completed; no class-group or point search was
launched.

For odd \(p\),
\[
\dim_{\mathbf F_2} E(\mathbf Q_p)/2E(\mathbf Q_p)
=\dim_{\mathbf F_2}E(\mathbf Q_p)[2]=k_p-1,
\]
where \(k_p\) is the number of cubic local factors. Take a sufficiently
small formal subgroup on which multiplication by 2 is invertible, and apply
the finite-group kernel/cokernel count. The three cubic factorization types
then give 0, 1 or 2 rational two-torsion dimensions.

All generic image ranks attain these dimensions. Completeness is proved by
a known finite dimension, not by failing to find another local point.

## Three paired results

The fourth column lists the complete local image dimensions at
\(3,5,7,11,13\). The generic image in the **product** has dimension \(s\),
equal to the sum of those dimensions.

| Family and parameter | Generic \(m\) | Observed quotient \(q\) | Local dimensions | Joint \(s\) | Added quotient support |
|---|---:|---:|---|---:|---:|
| MW16-05, \(307/206\) | 16 | 9 | \(1,1,2,1,1\) | 6 | 0 |
| MW16-05, \(-3158/1291\) | 16 | 0 | \(1,1,2,1,1\) | 6 | 0 |
| MW16-04, \(-1647/91\) | 16 | 9 | \(1,2,1,1,0\) | 5 | 0 |
| MW16-04, \(-2177/2397\) | 16 | 0 | \(1,1,1,1,0\) | 4 | 0 |
| published-R17, \(-2300/843\) | 17 | 7 | \(1,2,2,1,1\) | 7 | 0 |
| published-R17, \(-1561/3133\) | 17 | 0 | \(1,2,1,1,1\) | 6 | 0 |

The original matching qualifications remain: the first pair has a search
box/parameter-scale mismatch, the second shares its H4096 cohort but not
coefficient size, and the third matches coefficient size closely within
its cohort. Zero gain is observed and censored, not a rank upper bound.

**MW16-05.** Nine independent quotient directions accompany exactly the same
five-prime dimension profile as the control. No extra local component is
available beyond the generic product image.

**MW16-04.** The high fibre has an extra local dimension at 5: it has
\(I_0^*\) reduction with Tamagawa number 4, while the control has good
reduction. Its generic subgroup already spans both dimensions. This local
change supplies no new quotient support.

**Published R17.** At 7 the high fibre has \(I_0^*\) reduction and a
two-dimensional local image; the control has \(I_4\) reduction and a
one-dimensional image. Both have Tamagawa number 4 there. Generic
subgroups fill both images; Tamagawa number alone does not distinguish
even their local Kummer dimensions.

The curves have different cubic fields. These comparisons do not subtract
their full Selmer dimensions or identify their global classes.

## Explicit simultaneous generic corrections

Surjectivity separately at each place would not suffice: different generic
combinations might be needed at different primes. The joint rank proves
that **one** combination can match all five local images at once.

The report gives a correction for each of the 25 selected exceptional
directions. Let \(G_i\) be the generic points in their retained order and
\(Q_1\) the first selected quotient point. With one-based indices, examples
of points locally divisible by 2 at every \(p\in S\) are
\[
\begin{array}{c|l}
\text{MW16-05 high}&Q_1+G_1+G_8\\
\text{MW16-04 high}&Q_1+G_3+G_6\\
\text{published-R17 high}&Q_1+G_3+G_7+G_8.
\end{array}
\]

The cross-check multiplies the corresponding cubic Kummer representatives
and verifies a square in **every** completion above each selected prime,
using [nfislocalpower](https://pari.math.u-bordeaux.fr/dochtml/html-stable/General_number_fields.html#nfislocalpower).
Integral-basis denominators are cleared by a rational square. This avoids
a nonintegral-input issue in PARI 2.17.3 without changing the squareclass.

These corrections preserve quotient directions. They do not construct new
independent points or prove global divisibility by two. Rational solubility
of a genus-one two-cover and a rational half of a specified point remain
different questions.

## Where the quotient classes must reside

Let \(V=\operatorname{Sel}_2(E)\), \(G=\delta(M)\), and
\[
T_S=\bigoplus_{p\in S}\delta_p(E(\mathbf Q_p)),\qquad
\ell_S:V\longrightarrow T_S.
\]
All six curves satisfy \(\ell_S(G)=T_S\). For
\(K_S=\ker\ell_S\), elementary linear algebra gives
\[
\boxed{V/G\ \cong\ K_S/(G\cap K_S).}
\]
Subtract a generic class with the same local image from any class of \(V\).
The resulting surjection \(K_S\to V/G\) has kernel \(G\cap K_S\).

Every residual **Selmer** class, including an unobserved Sha class, can
therefore be represented with trivial localization throughout this
dictionary. This uses only generic surjectivity and does not require
exceptional points.

For the specified rational witness subgroup \(A\), the dimensions are:

| Fibre | \(\dim(G\cap K_S)=m-s\) | \(\dim(\delta(A)\cap K_S)=m+q-s\) |
|---|---:|---:|
| MW16-05 high | 10 | 19 |
| MW16-05 control | 10 | 10 |
| MW16-04 high | 11 | 20 |
| MW16-04 control | 12 | 12 |
| published-R17 high | 10 | 17 |
| published-R17 control | 11 | 11 |

The full dimension of \(K_S\) remains UNKNOWN. The observed gains are
extra independent rationally soluble classes in this global localization
kernel modulo its generic part. Constructing that kernel before supplying
points requires further arithmetic; the identity does not compute it.

## Mechanisms and the next bounded step

This rules out new quotient support at these five small odd primes as the
event explaining the three high gains. It does not rule out global effects
of changing local conditions, other primes, ideal-class constraints or
simultaneous higher-descent solubility.

The strongest remaining mechanism is a global construction or secondary
descent condition forcing several independent classes in
\(K_S/(G\cap K_S)\) to represent rational points. The fixed-field
experiments show why large admissible spaces and first norm equations alone
do not suffice.

Before expanding a bad-prime list or attempting a full class group, the next
small comparison should handle **2 and the real place** on these six
fibres. They are intrinsic to two-descent and deliberately absent from
this protocol. Test their generic images and then the joint image with
the current five primes. Two-adic units and real components may behave
differently; the present result makes no claim about them.

Agent 1 could use generic local surjectivity as an **incidence diagnostic
eliminating a proposed support feature**. It should not score low rank or
discard these high-gain fibres. No search policy changes follow.
