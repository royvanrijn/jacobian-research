# An explicit third carrier direction and its parameter images

The auxiliary Jacobian was proved to have rank three, but its recorded
descent supplied only two explicit points. The pointed quartic construction
already gives the x-coordinate of its opposite marked rational point.
Its two y-values follow from a quadratic that factors as y(y+a1*x+a3).
Transporting either sign to the recorded minimal Jacobian supplies a third
independent point, with one choice

```
x = 49005891476829078312687044788/314316088321
y = -4649890735583771052832412990871048905192266/176217857440197119
```

Both modulo3 and modulo5 exact finite-quotient calculations prove rank three
for the two old points and this point. Separate good-reduction witnesses
exclude rational3- and5-torsion; rational2-torsion causes no problem.
An independent ordinary-Python replay reconstructs the opposite points,
model transports and finite certificates. This gives an explicit full-rank
auxiliary subspace, without claiming saturation or an integral generator
basis. It does not increase the rank of the original known fibre.

The [auxiliary proof](../../artifacts/generated-results/elliptic-curves/native_rank3_carrier_marked_point_v1.json)
and [independent replay](../../artifacts/generated-results/elliptic-curves/native_rank3_carrier_marked_point_replay_v1.json)
precede a fixed125-word cube[-2,2]^3 in these three points. All words were
frozen before parameter inversion. The Sage build and replay pass, giving
124 finite images and one recorded birational exception. Initially79 images
were specialized while45 exceeded the declared512-bit parameter cap.

The [independent full height audit](../../artifacts/generated-results/elliptic-curves/native_rank3_fullspan_heights_v1.json)
now checks all124 finite images, including those45, by homogeneous integer
arithmetic. It verifies both native point identities and exact model
transports without an original-curve point search. There are75 distinct
parameters and75 j-invariants: the known origin and74 other parameters.
Sixty-two parameters are beyond the original twelve-image sample and the
origin. Every nonanchor image exceeds the400-bit normalized integral model
gate; the smallest proved coefficient lower bound remains519 bits. Large
invariants are retained exactly in hexadecimal to avoid decimal conversion
limits. This closes only the frozen cube's height question, not the full
rank-three carrier or larger point-exposure budgets.

The original twelve-image subset has separate rank-at-least19 proofs.
Specialized independence and catalogue comparisons for the additional62
parameters have not been claimed. No new near-record result follows from
this auxiliary construction alone.

Sources: `../cas/audit_native_rank3_carrier_marked_point.sage`,
`../cas/verify_native_rank3_carrier_marked_point.py`,
`../cas/audit_native_rank3_fullspan_images.sage` and
`../cas/verify_native_rank3_fullspan_heights.py`. Frozen ledgers are under
`native-rank3-carrier-marked-point-v1` and `native-rank3-fullspan-images-v1`
in `artifacts/local/elliptic-curves`.


## Exact halving gap closed

The three explicit independent auxiliary points together with the rational
point of order two generate a **2-saturated subgroup**. The fixed one-layer
audit checks all 15 nonzero parity combinations. None has a rational half.

For a short model y^2=x^3+Ax+B and S=(s_x,s_y), a half must have x-coordinate
among the roots of

\[
x^4-4s_xx^3-2Ax^2+(-8B-4As_x)x+A^2-4Bs_x.
\]

The certificate gives exact factorizations for all fifteen polynomials and
for the 2-torsion cubic. Every nonlinear factor has a good finite-prime
reduction with no root. The independent Python verifier multiplies all
factors back, checks those finite obstructions, classifies rational square
ordinates and verifies doubling on every rational lift. The torsion cubic
has exactly one rational root. This proves completeness of this halving
layer without a point search or a numerical height test.

The result rules out missing halves as the explanation for the earlier
coefficient cube's large images. It does not prove a full integral
Mordell-Weil basis, bound the whole carrier, or change any specialized rank.
The sources are `../cas/audit_native_carrier_halving.sage` and
`../cas/verify_native_carrier_halving.py`; compact certificates are
`native_carrier_halving_v1.json` and `native_carrier_halving_replay_v1.json`.
Both bounded build and independent replay pass under the fixed 120-second
stage limits, with logs in `native-carrier-halving-v1`.
