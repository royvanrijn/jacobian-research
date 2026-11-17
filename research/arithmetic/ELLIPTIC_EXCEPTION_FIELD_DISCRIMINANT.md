# Maximal-order discriminant of the certified elliptic Hasse example

This note certifies the cubic **field** discriminant in the
[elliptic Hasse example](ELLIPTIC_EXCEPTION_HASSE_EXAMPLE.md). Its original
coordinate and local-point certificate is preserved. A small integral basis,
four Eisenstein congruences, and the discriminant congruence modulo four
establish maximality without using PARI/GP in verification.

## Theorem EFD1: exact cubic field discriminant

Let `L=Q[W]/h(W)` be the cubic component of Theorem EEH1. Then

\[
 \boxed{D_L=-470381786514785447886543801463724621078040}
\]

and its exact factorization is

\[
 D_L=-2^3\cdot3\cdot5\cdot13^2\cdot1321^2\cdot9601^2
                \cdot2709173\cdot230703199^2.                 \tag{1}
\]

The quadratic resolvent has fundamental discriminant

\[
 d_K=-325100760=4(-81275190),
\]

and

\[
 D_L=d_K f^2,\qquad
 f=38037875815735627=13\cdot1321\cdot9601\cdot230703199.          \tag{2}
\]

**Proof: a reduced defining polynomial and an order.** Put

\[
 R(Z)=Z^3-Z^2-3321974487907578091Z
                     +228828238540688308908196032871.
\]

It is irreducible modulo seventeen. If `z` is its root, the integral-model
root `T` used in EEH1 is exactly

\[
 \alpha=\frac{2359787545600z^2-3111490635502052736000z
                   +59131265742732523839186831769600}{288929893}. \tag{3}
\]

Reduction modulo `R` gives `H(alpha)=0`, for the unchanged monic integral
cubic `H` of EEH1. Both polynomials have degree three and are irreducible,
so this gives an isomorphism of their fields. The relation between `H`
and `h` is the supplied rational scaling in the existing certificate.

Set

\[
 I=1733579358,\quad a=451757256,\quad B=2214649658454466309,
 \qquad w=\frac{z^2-az-B}{I}.
\]

The lattice `O=Z+Zz+Zw` is closed under multiplication, with table

\[
 \begin{aligned}
 z^2={}&2214649658454466309+451757256z+1733579358w,\\
 zw={}&-132574676509950914027+521026169z-451757255w,\\
 w^2={}&69761464271595910271-76610311059z-521026169w.
 \end{aligned}                                                 \tag{4}
\]

A finite free ring over `Z` is integral, so `O` is an order in `L`.
The determinant of the trace pairing on `(1,z,w)` is (1). Independently,
it equals `Disc(R)/I^2`; both exact computations are in the verifier.

**Proof: maximality.** At an odd prime with valuation one in (1), the
index of `O` in the maximal order cannot be divisible by that prime,
because the index squared divides the order discriminant.

At the four remaining odd primes, the following exact congruences show
that `R(Z+r)` is Eisenstein. Its quadratic and linear coefficients vanish
modulo `p`, and the table gives its constant coefficient divided by `p`,
which is nonzero modulo `p`.

| `p` | `r` | `R(r)/p mod p` |
|---|---:|---:|
| 13 | 9 | 8 |
| 1321 | 881 | 1315 |
| 9601 | 6401 | 4547 |
| 230703199 | 153802133 | 105958963 |

Consequently `Z_p[z]` is the full local ring of integers: an Eisenstein
root is a uniformizer in a totally ramified extension and the uniformizer
generates that ring. Since the integral order `O` contains `Z[z]`, it is
also maximal at these primes.

At two, the order discriminant has valuation three. An even index would
force the field discriminant to have valuation one, since its valuation
is nonnegative and differs by an even integer. This is impossible: an
integer of valuation one is two modulo four, whereas every number-field
discriminant is zero or one modulo four. Thus `O` is maximal at two as
well. No other prime divides its discriminant, so `O` is the maximal order.

The exact imported facts are Milne, *Algebraic Number Theory*, v3.00,
[Proposition 7.55 and Remark 7.56, printed pp. 122–123, PDF pp. 130–131](https://www.jmilne.org/math/CourseNotes/ANT300.pdf),
for the Eisenstein generator and its local ring, and Proposition 2.40(b),
printed p. 32, PDF p. 40, for the discriminant congruence. Their statements
and proofs were inspected in the [retained text](sources/milne-ant300.pdf).

Finally the squarefree part of (1) is `-81275190`, which is two modulo
four. Its quadratic fundamental discriminant is therefore four times that
integer. Division gives (2). All primes and complete factorizations in this
argument are checked by exact arithmetic and trial division. \(\square\)

The index of the original integral model in the maximal order is

\[
 [\mathcal O_L:\mathbb Z[T]/H]
       =2^{31}3^3 5^7\cdot7\cdot353\cdot2591\cdot9601^7.        \tag{5}
\]

Here `Z[T]/H` denotes its image in the identified field. Squaring (5)
and multiplying by (1) recovers precisely the original polynomial
discriminant. In particular the many extra primes and valuations in that
polynomial discriminant are not additional ramification conditions on `L`.

## Corollary EFD2: the field-discriminant Hasse criterion agrees

For this field, every prime dividing `f` in (2) is one modulo three, and
three does not divide `f`. The
[rank-five field-discriminant criterion](RANK_FIVE_DISCRIMINANT_CLASSIFICATION.md)
therefore proves that `K x L` is Hasse-failing. The quadratic fundamental
discriminant in (2) is the same field already identified in EEH1.

This is a second arithmetic verification route using the maximal order.
The existing whole-fibre identification and target certificate still
supply the realization by the fixed map. No statement about infinitely
many Hasse targets on the elliptic slice follows from this single example.

## Exact verification and scope

Run from the repository root:

```sh
python3 research/arithmetic/verify_elliptic_exception_field_discriminant.py --write-receipt
```

The data are
[elliptic_exception_field_discriminant_certificate.json](elliptic_exception_field_discriminant_certificate.json).
The [verifier](verify_elliptic_exception_field_discriminant.py) uses only
Python's standard library. It verifies (3), the multiplication table,
trace discriminant and index, the complete prime factorizations, and all
four Eisenstein rows. It reads the original certificate without editing it.
The [receipt](elliptic_exception_field_discriminant_checks.json) records
exact source identities and scope. No PARI/GP result or numerical embedding
is an assumption of the verifier. General local-ring, discriminant, and
Hasse-classification theorems remain explicitly imported written inputs;
this does not establish formal verification or external review.
