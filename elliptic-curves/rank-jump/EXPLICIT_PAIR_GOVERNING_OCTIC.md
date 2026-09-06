# An explicit octic for a pair's governing obstruction

Follow-up: [a cubic norm witness constructs the cochain without elliptic
points](UNPOINTED_NORM_COCHAINS_AND_THE_DYADIC_BLOCK_SWITCH.md). The retained
class-selected control gives an order-24 governing field and an exact
dyadic explanation of its known rational/Sha block switch. The rational-pair
construction and production replays below remain unchanged.

For two independent rational Kummer classes on an S3 elliptic curve, an
explicit even octic encodes the extra central bit needed to evaluate their
governing cochain. Its splitting field has degree 192 on all three frozen
production controls. At an unramified inert cubic prime, factor degrees
`1,1,3,3` give governing bit zero; `2,6` give bit one.

This closes the explicit-cochain endpoint proposed in
[the governing theorem application](GOVERNING_COCHAINS_SEPARATE_INCIDENCE_FROM_OBSTRUCTION.md).
It does **not** construct a carrier whose rational points produce a new jump.
The two supplied points are retrospective inputs. No original parameter,
twist, point search, class group, or full Selmer computation is performed.

## Formula and conventions

Fix the rational short model

\[
E:y^2=x^3+Ax+B,
\qquad P=(x_P,y_P),\quad Q=(x_Q,y_Q),\quad c=x_Q-x_P.
\]

Assume the cubic has Galois group S3 and the Kummer classes of P and Q are
independent. In particular, c and both y coordinates are nonzero. Choose roots
theta_i and square roots

\[
 a_i^2=x_P-\theta_i,\quad b_i^2=x_Q-\theta_i,
 \qquad\prod_i a_i=y_P,\quad\prod_i b_i=y_Q.
\]

Put delta_i=a_i+b_i and F=product(delta_i). All delta_i are nonzero since
their two squares differ by c. Then T=sqrt(F) satisfies

\[
\boxed{
h_{P,Q}(T)=T^8-4(y_P+y_Q)T^6
 +6c^2(x_P+x_Q)T^4-4c^3(y_Q-y_P)T^2+c^6=0.
}
\tag{1}
\]

The model, points and signs are part of this cochain choice. This is not an
invariant polynomial depending only on E or on the unrepresented classes.
For example, scaling x by d^2 and y by d^3 multiplies F by d^3, which can
change the extra quadratic character. The experiment uses precisely the
integral models already frozen in the bad-prime panel.

## Why this is the required cochain

Let V be the even-weight masks in F2^3. Its dot product is the Weil pairing
in additive notation, and S3 permutes its coordinates. Let L contain the
cubic roots and all a_i,b_i. Independence, the simple S3 module V, and
H1(S3,V)=0 give

\[
\operatorname{Gal}(L/\mathbf Q)=(V\times V)\rtimes S_3,
\qquad [L:\mathbf Q]=96.
\]

This uses the joint-class-field assertion in Morgan, Proposition 4.3;
the S3 hypotheses were checked in the preceding note. Write sigma=(e,f,g)
with target-indexed sign masks: sigma(a_i)=(-1)^{e_{g(i)}}a_{g(i)}, and
similarly for b and f. Let D=e+f and, for any v in V, put

\[
q_v=\frac{c^{|v|/2}}{\prod_{i:v_i=1}\delta_i}.
\]

Flipping only a_i replaces delta_i by c/delta_i; flipping only b_i replaces
it by -c/delta_i; flipping both replaces it by -delta_i. Because f has even
weight,

\[
\sigma(F)/F=q_D^2,
\qquad
\sigma(q_v)q_D=(-1)^{f\cdot gv}q_{gv+D}.
\tag{2}
\]

Consequently N=L(sqrt(F)) is normal. Define kappa by
sigma(T)=(-1)^kappa(sigma) q_D T. Equation (2) gives

\[
d\kappa(\sigma,\tau)=f_\sigma\cdot g_\sigma(e_\tau+f_\tau).
\]

The quadratic function q(v)=|v|/2 mod 2 on V is zero at zero and one
elsewhere; its polarization is the dot product. Therefore the correction

\[
\gamma(\sigma)=\kappa(\sigma)+e_\sigma\cdot f_\sigma+q(f_\sigma)
\]

satisfies exactly

\[
d\gamma(\sigma,\tau)=e_\sigma\cdot g_\sigma f_\tau=e\cup f.
\tag{3}
\]

There really is an extra central extension. If F were a square in L,
kappa would be a cochain on the order-96 group. Its coboundary would be
symmetric on any commuting pair. But on the commuting translations
(0,f,1) and (e,0,1), with f dot e=1, the displayed coboundary has values
one and zero in the two orders. This is impossible. Thus [N:L]=2 and
[N:Q]=192; this argument does not require a production class-group computation.

The eight displayed conjugates have labels (v,s) in V x F2 and values
(-1)^s q_v T. Their action is

\[
(v,s)\longmapsto(gv+e+f,\ s+\kappa+f\cdot gv).
\tag{4}
\]

The finite certificate checks (3) for all 96^2 pairs and (4) for all
4*96^2 pairs of lifted elements. There are 192 distinct permutations,
acting transitively on eight labels. On each production case, the exact
nonzero discriminant of (1) verifies that these labels give distinct
algebraic roots. Hence (1) is irreducible and its splitting field is N.
These conclusions use the written group argument and discriminant checks,
not a black-box production Galois-group command.

Degree eight is minimal for a **single polynomial encoding this full
governing extension**: a faithful permutation representation of a group of
order 192 cannot have degree at most seven, since 192 does not divide 7!.
This is not a minimality theorem about a simultaneous rational-point carrier.

## Exact octic identity and Frobenius evaluation

The four conjugates of F under a-sign changes are F and
c^2 delta_k/(delta_i delta_j). Their elementary symmetric functions are

\[
4(y_P+y_Q),\quad 6c^2(x_P+x_Q),\quad
4c^3(y_Q-y_P),\quad c^6.
\]

To check this without numerical roots, the script works in the free
rank-eight algebra with b_i^2=a_i^2+c. It verifies that all eight
coefficients of

\[
F^4-4(p+q)F^3+c^2(4\sum a_i^2+6c)F^2-4c^3(q-p)F+c^6
\]

vanish, where p=product(a_i), q=product(b_i). For the short cubic,
sum a_i^2=3x_P, yielding (1).

At a fixed-point-free Frobenius g, let u=(g-1)^(-1)e. The governing value is

\[
\psi=e_2(u,f)+\gamma.
\tag{5}
\]

Conjugate away the two translation parts. With e=f=0, gamma=kappa, and g
fixes zero while cycling the three nonzero elements of V. Formula (4)
then has cycle type (1,1,3,3) for psi=0 and (2,6) for psi=1. The finite
certificate checks all 64 lifts over the two three-cycles: 32 of each type.

An independent arithmetic evaluation avoids the octic. Over F_(p^3), let
theta be a cubic root and choose square roots of x_P-theta and x_Q-theta
with norms y_P and y_Q. Such roots exist because their norm squareclasses
are trivial; the prescribed norm selects the sign in this odd-degree
extension. Then

\[
\psi=0\quad\Longleftrightarrow\quad
\operatorname{Norm}_{\mathbf F_{p^3}/\mathbf F_p}(a+b)
\text{ is a square in }\mathbf F_p.
\tag{6}
\]

This follows by choosing Frobenius-compatible triples of roots, so e=f=0
and (5) becomes the quadratic Frobenius action on sqrt(F).

## Frozen production replay

The [protocol](EXPLICIT_GOVERNING_OCTIC_PROTOCOL.json) fixes the first two
exceptional directions in two old high-gain controls and the first two
generic directions in an observed-zero-gain control. Gain labels and point
orders are those of the immutable original panel, not a current rank census.
Each high pair is independently verified modulo its full marked generic
subgroup. The low pair has dimension two in E(Q)/2 and lies in its generic
subgroup. Its zero observed gain remains censored; no exact rank is asserted.

| Frozen fibre | Supplied pair | Pair dimension modulo generic | Eligible inert primes <=199 | psi=0 | psi=1 |
|---|---|---:|---:|---:|---:|
| A1/MW16-05 307/206, observed +9 | exceptional basis indices 16,17 | 2 | 15 | 5 | 10 |
| R17 -2300/843, observed +7 | exceptional basis indices 17,18 | 2 | 13 | 7 | 6 |
| R17 -1561/3133, observed +0 | generic basis indices 0,1 | 0 | 18 | 9 | 9 |

All 46 octic factorizations agree with (6). The
[arithmetic certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_explicit_governing_octic_v1.json)
retains the exact points, models, signatures, octics, discriminants, excluded
primes and values. The
[independent replay](../../artifacts/generated-results/elliptic-curves/rank_jump_explicit_governing_octic_verification_v1.json)
uses pure Python: Bareiss determinants of the three Sylvester matrices,
polynomial gcds with Frobenius powers for factor degrees, and Tonelli--Shanks
in the cubic finite fields. It reproduces every discriminant and all 46
values without importing Sage or the generating script.

This is a construction check, not a fitted discrimination experiment.
The full pair group and its two inert cycle types are universal under the
stated independence hypotheses. Different cochain gauges also preclude
treating the raw bit counts as comparable rank features.

## A concrete condition on a twist parameter

Clear denominators by choosing the recorded positive integer d_h and put
H(U)=d_h^8 h(U/d_h), a monic integral octic. Its nonzero discriminant D_H
is an explicit integer, independently verified. Every finite ramified prime
of its splitting field divides D_H. This is a certified support superset;
we have not factored it or claimed an exact field discriminant.

On the frozen integral model, let D_f be the cubic discriminant and put
M=8|D_f D_H|. This deliberately large integer requires no factorization.
For a positive prime twist parameter ell satisfying

\[
\ell\equiv1\pmod M,\qquad
f\bmod\ell\text{ irreducible},\qquad
H\bmod\ell\text{ has factor degrees }(2,6),
\tag{7}
\]

the quadratic character of Q(sqrt(ell)) is locally trivial at 2, infinity,
all bad places and governing ramification. It ramifies only at a new inert
prime. Morgan's
[Proposition 3.3](https://arxiv.org/pdf/2309.02374v2) then gives
CT_(E^(ell))(P,Q)=1, since the original rational pair has CT=0. Lemma 2.8
preserves every local Kummer image and the full Sel2 group. Consequently
every nonzero class in this two-dimensional block is nonrational on the
twist: the nonsingular restricted pairing excludes it from the full radical.

Replacing (2,6) by (1,1,3,3) in (7) forces this one CT entry to remain zero.
It does not force either class to have a rational point, or the block to
annihilate the rest of Sel2. The preceding governing theorem and its
abelianization argument give existence of primes of both types compatible
with the congruence; none is computed here. The small-prime replay table
was not required to meet (7) and is not a twist experiment.

Thus we now have an explicit **rational-to-Sha pair switch at fixed Selmer
incidence**. The parameter in (7) is a twist parameter. It is not t in the
original R17/A1 families and twisting need not preserve their rational
generic subgroups. This distinction is essential.

## Lessons and next goals

1. **Solubility obstruction:** governing cochains are now computable for a
   retained rational pair, including their extra central bit and a certified
   ramification bound. This is the strongest new mechanism-level result.
   A family-wide degeneration of the relevant obstruction matrix remains a
   candidate explanation, not an observed or proved cause of the jumps.
2. **Weak discrimination:** merely having a degree-192 pair field, an octic,
   or both inert Frobenius types cannot distinguish large gains. These
   properties also hold for independent generic classes in the low control.
   Further enumeration of such fields from exceptional point pairs is a
   lower-priority goal.
3. **Missing computation:** construct governing data for a class specified
   before exceptional points are supplied. The formula above exploits the
   special representatives x_P-theta and x_Q-theta, whose difference is a
   rational scalar. Arbitrary Selmer representatives need additional descent
   data; this experiment has not removed that requirement. A generic-point
   pair gives an available calibration, but no new-rank signal by itself.
4. **Missing implications:** a zero restricted CT block must first be tested
   against the entire Selmer space. Even its full-radical membership permits
   higher-divisible Sha. A rational construction or higher descent is still
   needed to reach simultaneous rational solubility, followed by exact
   independence modulo the specialized generic subgroup to reach a jump.
5. **Information Agent 1 could eventually use:** an inexpensive exact nonzero
   obstruction for independently specified prospective classes could reject
   their proposed rational block. A zero bit is not a rank score. No current
   search policy change follows from this retrospective octic construction.

The next bounded endpoint should be a **non-point-supplied class with a
certified obstruction calculation**, rather than another known-point pair
or a larger prime table. This separates constructing the arithmetic
obstruction from observing that a supplied rational point makes it vanish.

## Replay

```bash
/home/royvanrijn/.local/bin/sage -python elliptic-curves/rank-jump/explicit_governing_octic.py check
python3 elliptic-curves/rank-jump/verify_explicit_governing_octic.py check
```

Workers are bounded at 60 seconds per frozen pair and checkpointed in the
ignored local directory. Both replays complete in seconds. Source hashes
bind the immutable inputs and scripts. Concurrent search, theorem navigation,
MATH_STATUS and STATUS are untouched by this change.
