# All in-field radical operations on the retained norm pool are exhausted

The retained 4134-element norm dictionary cannot produce an additional
Selmer class by **any sequence of multiplication, inversion, rational
rescaling, norm projection or integer-root extraction inside its cubic
field**. This excludes arbitrarily many rounds, not just the single root
previously tested.

The exact calculation is:

| Quantity | Dimension or index |
|---|---:|
| Original norm-one multiplicative group, free rank | 4150 |
| Its original squareclass image | 4149 |
| Index in its complete 2-saturation | 2 |
| Squareclass image after the one required root | 4150 |
| Independent good-prime ramification directions in that image | 4134 |
| Intersection with the elliptic 2-Selmer group | exactly the generic 16-dimensional subgroup |
| Additional Selmer classes from the full in-field radical closure | **0** |

This strengthens
[the single-root affine obstruction](RELATION_ROOTS_CREATE_RAMIFICATION_NOT_THE_MISSING_BLOCK.md)
to a complete saturation theorem for its frozen input group. It gives a
rigorous stopping rule for that constructor. The extra classes required
by the large-jump problem must use a genuinely new nonrational generator
outside this radical closure, or a construction outside the allowed
operations.

The experiment remains on the MW16-05, t=3/17 generic control. It supplies
no new class basis for the fresh high/low pairs and does not bound the
full Selmer or cubic class group. No exceptional points or classes enter.

## Normalize before taking roots

Write K for the fixed cubic field and N for its rational norm. Define

\[
 D(a)=a^3/N(a),\qquad \pi(a)=N(a)a.
\]

Then D is a multiplicative map into the norm-one group K^1, with

\[
 N(D(a))=1,\qquad
 D(a)/\pi(a)=(a/N(a))^2,\qquad D(ca)=D(a)\ (c\in\mathbb Q^*).
\tag{1}
\]

Thus D and pi give the same cubic squareclass. Normalization removes
rational scalars exactly, rather than treating them as an untracked
source of square roots. For any norm-square element gamma, [D(gamma)]
also equals [gamma].

The norm-one group is torsion-free. Any root of unity in a cubic field
has degree dividing three; the only possibilities are 1 and -1, and
N(-1)=-1. In particular two norm-one elements with the same square are
equal.

Let a_j be the 4134 raw dictionary elements and gamma_i the sixteen
generic Kummer representatives. Put

    u_j=D(a_j),      v_i=D(gamma_i),
    Γ0=<u_j, v_i> inside K^1.

The independent calculation uses the original input polynomials and
complete retained ideal-parity certificates. It does not factor a new
integer or initialize a number field.

## The sole root gives an exact index-two extension

The unique dictionary element with zero retained outside-S parity is
at primitive address (20941,464), zero-based dictionary index **1769**.
Call its index j0. The previous exact relation is

\[
 \pi(a_{j0})\prod_{i\in I}\gamma_i=w^2,
 \qquad I=\{2,5,8,12,13,15\},
\]

with one-based generic indices and mask 22674. Applying D to this exact
equality, and using D(pi(a))=D(a), gives

\[
 u_{j0}\prod_{i\in I}v_i=h^2,
 \qquad h=D(w).
\tag{2}
\]

The certificate supplies the exact norm-one coefficients of h, u_j0
and all sixteen v_i, and verifies their norms and equation (2).

Define the replacement group

    Γ1=<h, u_j (j!=j0), v_i>.

It has 1+4133+16=4150 displayed generators. Their squareclasses are
independent, for two separately certified reasons:

1. On the retained good-prime ideal columns, the 4133 remaining u_j
   have rank 4133. The h row raises that rank to **4134**. These are the
   same parity rows as pi(a_j) and pi(w), by (1).
2. Every generic v_i has even valuation outside S: its squareclass is a
   rational elliptic Kummer class and hence unramified at good odd
   places. Consequently a square relation among the replacement
   generators first forces all h and u_j coefficients to vanish modulo
   two. The sixteen remaining generic classes are independently
   certified by finite-field square characters of rank **16**.

The new character replay uses 38 split-prime blocks with p<=2003.
These generic classes are the only points used; no exceptional points
or outcome labels select a block. The good-prime parity ranks are
recomputed from the fixed 10606 prime-ideal columns. A restricted set of
columns suffices to prove independence; no assertion about untested
prime columns is needed.

Squareclass independence also proves that the 4150 replacement
generators are multiplicatively independent. Any nonzero integer
relation can be divided by its largest common power of two in the
torsion-free group, leaving a relation with an odd exponent, which would
contradict squareclass independence.

Thus Γ1 is free on this displayed basis. Equation (2) expresses the old
generator as

    u_j0=h^2*product(v_i : i in I)^(-1).

The inclusion Γ0⊂Γ1 has determinant two in these bases. The verifier
checks the 17-by-17 affected block; the other 4133 coordinates are
unchanged. Hence

\[
 [\Gamma_1:\Gamma_0]=2.
\tag{3}
\]

## Why no second root can add a direction

Suppose z∈K^1 and z^2∈Γ1. Write z^2 as a product of the 4150 basis
generators with integer exponents. Its class is zero in K*/K*2.
Independence forces every exponent to be even, so z^2=g^2 for some
g∈Γ1. Torsion-freeness gives z=g. Therefore

\[
 \boxed{\Gamma_1\text{ is 2-saturated in }K^1,
 \qquad\operatorname{Sat}_2(\Gamma_0)=\Gamma_1.}
\tag{4}
\]

This is an infinite statement proved from the finite independence
certificate. It does not rely on failing to find a second root.

Allow roots of arbitrary positive integer order, and define

    Rad_K(Γ0)={z in K^1 : z^n belongs to Γ0 for some n>=1}.

If n=2^e*m with m odd, (4) implies z^m∈Γ1. Since m is odd, z and z^m
have the same squareclass. Conversely Γ1⊂Rad_K(Γ0) by (2). Therefore

\[
 \boxed{\operatorname{image}(\operatorname{Rad}_K(\Gamma_0)
               \to K^*/K^{*2})
       =\operatorname{image}(\Gamma_1\to K^*/K^{*2}).}
\tag{5}
\]

Odd-order roots could enlarge the abstract group, but cannot enlarge
this squareclass image. The statement also covers products and roots
interleaved in any finite sequence, since the radical closure is a
subgroup.

There is no rational-rescaling loophole. If a raw element z satisfies
z^n=c*product(a_j^e_j)*product(gamma_i^f_i), with arbitrary c∈Q*, then
D(z)^n lies in Γ0 by (1). The resulting norm projection pi(z) has the
same squareclass as D(z), so is already covered by (5). Norm projection
itself does not escape the normalized group, since D(pi(z))=D(z).
If the final raw element already has square norm, its class likewise
equals its D-normalization.

## The entire closure still has no extra Selmer class

Let V1 denote the squareclass image in (5), and G the generic subgroup.
The outside-S parity argument above shows exactly

    kernel(V1 -> retained good-prime valuation parities) = G.

Every elliptic Selmer class is unramified at these places, while all of
G consists of rational Kummer classes and is in Selmer. Hence

\[
 \boxed{V_1\cap\operatorname{Sel}_2(E/\mathbb Q)=G.}
\tag{6}
\]

Intersecting with the strict group U also gives V1∩U=G∩U. All possible
new directions in this radical closure fail incidence already at good
primes; rational solubility and point-search visibility never enter.

This does not exclude multiplying a representative by an arbitrary
square z^2 with z imported from outside the raw radical closure and then
extracting its root. Such a step would import the new class pi(z)
itself. It is a new generator, not an operation on the frozen dictionary.
Nor does the theorem cover adjoining fields and subsequently applying
norms or other constructions back to K.

## What this changes in the research plan

The previous bounded relation-root experiment left open repeated root
extraction. That route is now closed completely on the retained pool,
including arbitrary rational rescalings and odd roots. Enlarging the
number of root rounds has no mathematical prospect of supplying the
missing strict block.

The next constructor must change the nonrational generator source:
for example, an independently obtained ideal class and its square
principalization, with the full unramified/local-square certificates.
The preceding
[fixed-resolvent obstruction](FIXING_THE_RESOLVENT_POLYNOMIAL_TESTS_SOLUBILITY.md)
means that such a construction must allow arbitrary cubic generators.
Neither result supplies those new classes or a condition on the original
family parameter t; the class creation of the largest jumps remains open.

The current repository audit still finds incomplete class-group stages
in the retained fresh/historic descent records. Their missing bases were
not treated as zero-dimensional and no previously timed-out class-group
campaign was restarted. The separate public upgrade of inventory188
remains an outcome/visibility control, outside this arithmetic input.

## Reproduction

The [protocol](RETAINED_RADICAL_CLOSURE_PROTOCOL.json) fixes the existing
pool and the finite generic-character check. The
[main certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_retained_radical_closure_v1.json)
contains the norm-one root relation and the replacement ranks. The
[independent verifier](../../artifacts/generated-results/elliptic-curves/rank_jump_retained_radical_closure_verification_v1.json)
rebuilds valuation parity using 3v_P(a)-e_P*v_p(Na), eliminates with the
opposite pivot order, checks norm-one identities with rational companion
matrices, and recomputes finite characters using Sage finite fields.
The retained lattice-based ideal valuation verifications are hash-bound
as upstream evidence; they are not rerun.

```sh
timeout 30 python3 elliptic-curves/rank-jump/retained_radical_closure.py check
timeout 30 sage -python elliptic-curves/rank-jump/verify_retained_radical_closure.py check
```

Only new analysis files and immutable rank-jump artifacts are added.
Search protocols, candidate populations, scoring, workers and
mathematical-status entries remain untouched.
