# A coefficient-only obstruction excludes the horizontal mechanism on A1/MW16-05

The exact equal-ordinate mechanism from the
[two-direction control](ONE_SQUARE_CONDITION_TWO_RATIONAL_DIRECTIONS.md)
cannot occur on the A1/MW16-05 fibre at 307/206, whose retained quotient
gain is +9. This is an equation-only exclusion, valid for **every** rational
point pair on the curve, not merely the supplied point list.

In a completed short model y²=x³+Ax+B, this fibre has

\[
 \boxed{v_{23}(-A)=1,\qquad23\equiv2\pmod3.}
\]

Two distinct rational abscissae with the same cubic value would require
−A to be a norm from Q(√−3), which this odd valuation forbids.
The obstruction is invariant under rational Weierstrass changes.

It also excludes equal-ordinate rational multisections over Q(t) in the
A1/MW16-05 family, of any rational-function degree: a solution would
specialize to a forbidden norm representation at this smooth rational
fibre. The precise specialization argument is below.

This rules out one specific simultaneous-solubility mechanism. It does
not rule out high rank, a shared quadratic cover with unequal ordinate
multipliers, or a different low-degree construction.

## Completed ordinates are invariant; raw ordinates are not

For a general Weierstrass equation, use

\[
 \eta(P)=2y(P)+a_1x(P)+a_3.
\]

Under x=u²x′+r and y=u³y′+u²sx′+t, the transformed coefficients give
η=u³η′. Hence η(P)=±η(Q) is invariant. This is the usual Weierstrass
coordinate convention documented by
[Sage's Weierstrass morphism reference](https://doc.sagemath.org/html/en/reference/arithmetic_curves/sage/schemes/elliptic_curves/weierstrass_morphism.html).
Completing the square replaces y by η/2, and removing the quadratic
x term leaves that ordinate unchanged.

In contrast, any two known points with distinct x can be made to have
equal **raw** y′ by choosing s to be their secant slope. The equation
then acquires a y′x′ term. This operation does not make their completed
ordinates equal. A purported common-value signal based on raw ordinates
would therefore be a coordinate artefact.

The exact experiment checks both statements on all six retained paired
production controls: a fixed nontrivial Weierstrass change preserves the
completed-ordinate groups, while a deliberately chosen shear manufactures
raw-y equality and the invariant detector correctly rejects it.

## One horizontal triple has at most two directions

For a short cubic F(x)=x³+Ax+B, suppose F(a)=F(b)=d with a≠b. Then

\[
 F(x)-d=(x-a)(x-b)(x+a+b).
\]

Over the shared cover w²=d, the three points with ordinate w satisfy

\[
 (a,w)+(b,w)+(-a-b,w)=O.
\]

Thus one horizontal triple supplies at most two independent directions,
and at most two modulo an existing generic subgroup. A repeated third
intersection or a further relation can lower this bound.

The small mechanism achieves this bound generically, but its n=c=1
specialization has three displayed points in a single horizontal group
and only one retained direction, already generic. Merely counting the
three points would give the wrong answer.

If a gain q were entirely explained by such triples, at least ⌈q/2⌉
triples would be needed: 4,5,7 for gains +8,+10,+14. This is **not** a
bound on the number of auxiliary covers. One quadratic cover can carry
several distinct triples or other multisections.

## The equation-only norm obstruction

Put z=(a+b)/2 and ℓ=(b−a)/2. Subtracting F(a)=F(b) gives

\[
 a^2+ab+b^2+A=0,
 \qquad\boxed{\ell^2+3z^2=-A.}
\]

The common value, if the conic is solved, is

\[
 d=B-2Az-8z^3.
\]

This separates two gates even within the horizontal mechanism:

1. **Incidence of rational abscissae:** the norm conic must have a
   rational solution with ℓ≠0.
2. **Rational solubility at a fibre:** its common value d must be square.

Solving the conic alone does not prove rational points. Failing its
necessary local norm condition rules out the mechanism altogether.

For an odd prime p≡2 mod 3, −3 is a nonsquare modulo p. If ℓ,z are
rational, factor out the smaller of their p-adic valuations. The remaining
two integral coordinates are not both divisible by p, so their norm is
nonzero modulo p. Therefore

\[
 v_p(\ell^2+3z^2)=2\min(v_p(\ell),v_p(z))
\]

is even. The conclusion also holds if one coordinate is zero. This proves
the obstruction at p=23 without a global factorization or a norm solver.
The independent check enumerates all 529 residue pairs modulo 23 and
finds that the only zero of ℓ²+3z² is (0,0).

Since A scales by u⁴ under a short Weierstrass change, its valuation
parity is invariant. Equivalently c₄=−48A; at p=23 the same obstruction
is visible directly in c₄. The retained model has

```
c4 = 8572015517453275854364236700326495593653786816
v23(c4) = 1
```

The [independent certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_horizontal_norm_gate_verification_v1.json)
recomputes c₄ directly from the original five Weierstrass coefficients
with Sage and confirms the valuation, separately from the producer's
short-model rational arithmetic.

## Why this excludes the exact mechanism over the original parameter line

Take a short Weierstrass model regular near a smooth rational fibre t=t₀.
Suppose ℓ,z∈Q(t) solved ℓ²+3z²=−A(t). If either had a pole at t₀,
let N>0 be the largest pole order. After multiplication by (t−t₀)²ᴺ,
specialization would give

\[
 \ell_0(t_0)^2+3z_0(t_0)^2=0
\]

with at least one nonzero rational coordinate. That is impossible over Q.
So both functions are regular and specialize to a rational norm solution.

At the A1/MW16-05 fibre 307/206 no such solution exists. Consequently
there is no rational-function pair of distinct abscissae giving identical
completed cubic values over Q(t) for this family. A rational change to a
regular local Weierstrass model preserves the completed-ordinate property,
so a denominator in some other global presentation cannot evade the
argument. No degree bound on ℓ or z is needed.

The small family E_(m,c) has an explicit norm representation for every
rational m,c:

\[
 \ell=3/2,\qquad z=1/2+m/3,\qquad
 -A=\ell^2+3z^2.
\]

Thus the obstructed A1 fibre cannot be a rational Weierstrass presentation
of a member of that exact small family either. The constructive control
remains valid; its direct transfer to this production family is excluded.

## Three paired retrospective comparisons

The [point-list audit](../../artifacts/generated-results/elliptic-curves/rank_jump_completed_ordinate_blocks_v1.json)
uses only the already retained generic and witness lists, deduplicated up
to sign and repeated x. The separate norm detector reads a masked file
containing **only equations and case numbers**.

| Retained family and parameter | Generic / independent witness / observed quotient | Horizontal groups in supplied lists | Coefficient-only gate |
|---|---:|---:|---|
| A1/MW16-05, 307/206 | 16 / 25 / +9 | 0 | EXCLUDED at 23 |
| A1/MW16-05, −3158/1291 | 16 / 16 / observed 0 | 0 | UNKNOWN |
| A1/MW16-04, −1647/91 | 16 / 25 / +9 | 0 | UNKNOWN |
| A1/MW16-04, −2177/2397 | 16 / 16 / observed 0 | 0 | UNKNOWN |
| Published R17, −2300/843 | 17 / 24 / +7 | 0 | UNKNOWN |
| Published R17, −1561/3133 | 17 / 17 / observed 0 | 0 | UNKNOWN |

For UNKNOWN rows, absence from the supplied list does not exclude hidden
combinations or other multisections. The fixed dictionary consists only
of primes 5≤p≤1999 with p≡2 mod 3. It gives no complete norm decision
when it finds no obstruction. No larger factorization or prime search
was launched. Observed zero gains remain censored.

The small controls (c,n)=(1,5),(3,7),(1,1) each give one detected
horizontal group. Their previously certified quotient ranks are 2,2,0,
respectively. The two positive controls also provide exact norm witnesses,
and neither is falsely excluded by the coefficient gate.

## Reproducibility and consequence for mechanism selection

The [point-list protocol](COMPLETED_ORDINATE_PROTOCOL.json) caps each input
at 64 points and forbids new subgroup enumeration. The
[norm protocol](HORIZONTAL_NORM_GATE_PROTOCOL.json) freezes the six
equations and small-prime dictionary. The
[masked input](../../artifacts/generated-results/elliptic-curves/rank_jump_horizontal_norm_gate_inputs_v1.json)
contains no points or ranks; the
[result](../../artifacts/generated-results/elliptic-curves/rank_jump_horizontal_norm_gate_v1.json)
retains every nonzero tested valuation. Five narrow tests cover negative
coordinate scaling, misleading shears, dependent triples and rational
valuation edge cases.

```sh
python3 elliptic-curves/rank-jump/completed_ordinate_blocks.py check
python3 elliptic-curves/rank-jump/horizontal_norm_gate.py check
sage -python elliptic-curves/rank-jump/verify_horizontal_norm_gate.py check
python3 -m unittest discover -s elliptic-curves/rank-jump -p test_completed_ordinate_blocks.py
```

The ranked interpretation is:

1. **Still viable:** several multisections sharing a nontrivial quadratic
   cover, with a proven independent rational image after specialization.
2. **Excluded here:** the exact horizontal-triple version on A1/MW16-05.
   The obstruction applies before exceptional points are supplied and
   prevents this route from explaining that family's +9 fibre.
3. **Weak or misleading:** raw-y equality and counting points in a triple.
   Moreover, after rational points are supplied, their cubic values are
   already squares; a shared squareclass observed only at that fibre is
   tautological. The meaningful shared class must be established over
   the original parameter field before specialization.
4. **Next missing computation:** test a broader identity
   F_t(x_i(t))=d(t)h_i(t)², with unequal h_i, rather than identical
   values. The norm obstruction above does not address that identity.
   Its generic rank and exceptional quotient accounting must still be
   proved. No production instance is currently supplied by this audit.
5. **Potential use for Agent 1:** this coefficient-only gate can exclude
   a proposed horizontal construction, but cannot exclude a high-rank
   fibre or rank candidates. The active search remains unchanged.
