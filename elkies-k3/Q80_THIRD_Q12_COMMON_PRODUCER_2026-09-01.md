# Q80 third-q12 common producer checkpoint

## Status

The fixed `u=-2` construction now has two distinct layers.

First, `p=19` is complete as an exact finite-field child: the resolved
genus-one pencil, generic Jacobian, birational maps in both directions,
minimal `I6+I4+3I2+8I1` model, transported `A5+A3+3A1` component marking,
and Frobenius-invariant coefficient encoding all have independent literal
replay certificates.

Second, one common polynomial-closure producer selects the same horizontal
orbit and connected genus-one divisor at `p=19` and `p=61`.  The p=61
realization is now complete as well: 72 exact mapped fibres, the generic long
Jacobian, birational maps in both directions, the minimal short child, and the
transported component/zero marking all pass literal replay.

## Prime-independent producer

The producer starts from one audited good-reduction surface at rational `u`.
It exports the six-variable polynomial-section closure, decodes every
degree-at-most-two RUR support point in one quadratic finite field, tests all
pairwise differences, and quotients accepted results by section sign and
Frobenius.  Its acceptance rule is exactly:

- `P.O=2`;
- canonical height `8`;
- identity component at the finite `I1*` (`D5`) fibre;
- identity component at the infinite `I3*` (`D7`) fibre.

At each of 19 and 61 the closure has squarefree geometric degree 12.  The
producer tests 156 relative-sign pairs, retains two unsigned candidates, and
finds one Frobenius/sign orbit.  The RUR happens to split differently at the
two primes; this is local factorization data and is not part of the alignment.

The generic connected compiler then gives, at both primes:

| invariant | value |
| --- | --- |
| Smith degrees | `(0,0,6)` |
| saturated ambient dimension | `7` |
| complete D7 gate rank | `4` |
| combined D7+D5 gate rank | `5` |
| pencil dimension | `2` |
| moving degrees `(new base,W,x)` | `(2,9,3)` |
| divisor square / old-fibre degree | `0 / 3` |
| generic member | primitive, separable, irreducible, genus one |

The resolved-pencil alignment certificate is
`artifacts/generated-results/q80-third-q12-um2-p19-p61-common-producer-alignment.json`.
The stronger completed-child certificate is
`artifacts/generated-results/q80-third-q12-um2-p19-p61-full-child-alignment.json`.
It aligns 1,947 ordered Frobenius-invariant coefficient slots and deliberately
does not identify the two local quadratic generators.

## First-marking field test

The stretch hypothesis has a decisive local answer.  Write the exact
first-marking genus-two cover as `w^2=f(t)`, with
`t=u-u_CM24`, and specialize at `u=-2`.  The value `f(-2-u_CM24)` is a
nonzero nonsquare rational number, so it defines a quadratic number field.
Its reduction is:

- `16` at 19, a square: the number field splits at 19;
- `2` at 61, a nonsquare: the number field is inert at 61.

Therefore the quadratic p=19 horizontal field cannot be the reduction of
this first-marking number field.  At p=61 the residue fields are isomorphic
(the radicand-to-target-discriminant ratio is the square `12`), but that
abstract local agreement does not identify the horizontal as the reduction
of a characteristic-zero section.  Descending this particular genus-two
cover cannot explain the mandatory p=19 control.

The replay certificate is
`artifacts/generated-results/q80-first-marking-field-um2-local-behavior.json`.

## Normalization breakthrough and next reconstruction path

The original p=61 mapped-fibre worker spent more than eight minutes and over
1.1 GB in Sage's non-prime-constant-field maximal-order routine.  A second
fibre showed the same behavior.  The cause is algorithmic: that routine uses
a characteristic-`p` power map, so its intermediate degrees grow with `p`.

`sample_q80_third_q12_weierstrass_modp2.py` now uses Singular's
Grauert--Remmert module normalization and the same reversed-Hermite module
reduction as Sage.  The formerly blocked p=61 fibre completes in about one
second.  A deterministic batch then retained all 72 attempted fibres (64
training and eight held out), after which generic interpolation and literal
map replay recovered:

- long coefficient degrees `(2/2,4/4,4/4,6/6,8/8)`;
- forward-map block degrees identical to p=19;
- inverse weighted bounds `4` and `10`;
- minimal degrees `(8,12)` and fibres `I6+I4+3I2+8I1`;
- the `R5` zero and old-zero alignment with the unique `I6` factor.

The coefficient interface uses trace, norm, and coefficient discriminant as
generator-free CRT values.  The anti-invariant coefficient remains local
until a characteristic-zero quadratic field is reconstructed.  Section sign,
base `PGL2`, and Weierstrass scaling have separate pinned ledgers.  The next
step is therefore additional-prime collection with held-out primes, not
another horizontal search.

That collection has begun.  The complete common pipeline also passes at
`p=67`: one horizontal orbit, the same connected pencil, 72/72 mapped fibres,
the same generic map degree profile, minimal `I6+I4+3I2+8I1`, the transported
`R5` marking, and 1,947 invariant slots.  The CRT interface now accumulates
the generator-free values at `19`, `61`, and `67`, with combined modulus
`77653`; `p=71` is certified good and reserved for held-out replay.  The
interface artifact is
`artifacts/generated-results/q80-third-q12-um2-frobenius-crt-interface.json`.
This modulus is intentionally reported as too small for rational
reconstruction.  No equation over `QQ` is asserted at this checkpoint.
