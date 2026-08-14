# Prime-power discriminant engineering by CRT and rational lattices

## Exact local-to-global core

Let \(F\in\mathbb Z[t]\) be a squarefree-base factor of a family
discriminant.  Roots modulo \(p\) are lifted one digit at a time: from a root
\(r\bmod p^e\), test

\[
r+c p^e,\qquad 0\le c<p,
\]

modulo \(p^{e+1}\).  Exhausting every digit handles simple and singular roots
uniformly.  Several compatible congruences combine to

\[
t\equiv r\pmod M.
\]

For a rational parameter \(t=a/b\), with \(p\nmid b\) for every \(p\mid M\),
this becomes

\[
a-rb\equiv0\pmod M.
\]

The solutions form the determinant-\(M\) lattice

\[
L_{r,M}=\langle(M,0),(r,1)\rangle.
\]

This is a homogeneous two-dimensional shortest-vector problem.  Exact
Lagrange--Gauss reduction is enough; LLL or BKZ becomes relevant only after
there are several parameters or genuinely higher-dimensional coupled
constraints.  The normal size scale remains \(\sqrt M\): reduction can expose
exceptional residues, but it does not make height independent of \(M\).

## Primitive-pair gate

The shortest lattice vector need not define a valid rational residue.  If
\(d=\gcd(a,b)\) shares a factor with \(M\), dividing by \(d\) can destroy the
congruence.  Every proposed pair is therefore normalized and then rechecked:

1. divide by \(\gcd(|a|,|b|)\) and choose \(b>0\);
2. require \(\gcd(b,M)=1\);
3. recheck \(a-rb\equiv0\pmod M\).

For example, the constraints

\[
t\equiv1\pmod{5^3},\qquad t\equiv-1\pmod{7^2}
\]

give \(r=1126\bmod6125\).  A Gauss basis is
\((-49,-49),(38,-87)\).  Its shortest vector reduces to \(t=1\) and is
invalid, while \((-38,87)\) is valid and satisfies

\[
a^2-b^2=-6125=-5^3 7^2.
\]

The tests exhaust the entire height box rather than trusting a fixed window in
reduced-basis coordinates.

## Homogenization and valuations

If \(F\) has degree \(d\), use

\[
F^h(a,b)=b^dF(a/b).
\]

The affine root condition plus \(p\nmid b\) implies
\(p^k\mid F^h(a,b)\).  It promises only \(v_p(F^h)\ge k\); exact valuations
are computed on the global pair.  A zero value is rejected because it gives a
singular specialization.

## From raw discriminant to conductor

Large powers in a raw discriminant are useful only if they survive
minimalization.  For clean multiplicative reduction at \(p\ge5\), the minimal
discriminant may contain \(p^k\) while the conductor pays only \(p\).  The
following cases must be audited separately:

- primes 2 and 3 and their wild terms;
- primes dividing contents, denominators, factor discriminants, or resultants;
- collisions between discriminant factors;
- additive reduction and nonminimal specializations;
- uncontrolled cofactors, which may dominate the conductor.

Accordingly, every candidate promoted beyond the local-seed stage must pass an
exact global minimal-model and local-reduction computation.  The calibration
manifest does this with PARI/GP's `ellglobalred`; the Fermigier seed has not yet
crossed that promotion gate.

## Pinned calibration result

For the rank-at-least-two calibration family

\[
E_t:y^2=x^3-t^2x+t^2,
\]

the independent sections are \(Q=(t,t)\) and \(R=(-t,t)\).  PARI's exact
rank bounds at \(t=5\) are `[2,2]`, with specialized points `(5,5)` and
`(-5,5)`, so any generic relation would contradict their independence.  The
eight-class search at \(23^3,47^2,73^2\) finds

\[
t=-110627/84367,\qquad
27b^2-4a^2=23^3 47^2 73^2.
\]

The integral specialization is global minimal.  Its shaped minimal-
discriminant valuations are `(3,2,2)`, while all three conductor exponents are
one; PARI reports trivial torsion and rank bounds `[3,3]`.  This is EC-CRT1,
not a record candidate.

## Local constraint optimization

The constraint layer should keep distinct feature tables:

```text
shaping prime p: root/lift, v_p(Delta), reduction type, local Euler term,
                 Tamagawa and root-number features
good prime q:    t mod q (and denominator square class if needed), a_q,
                 point count, information/cost score
```

A beam, meet-in-the-middle search, or rare-event splitting scheme can choose
local symbols.  CRT assembles them, and exact rational reconstruction supplies
the height penalty.  Scoring must be trained or calibrated on held-out
specializations; an additive Nagao score is a heuristic ranking, not a rank
certificate.

The residue-only local table resolves good and clean multiplicative fibers.
When both `c4` and the discriminant vanish modulo a prime, it records
`unresolved_bad`; additive versus nonminimal or globally singular behavior is
decided only after an actual rational specialization is minimized.

## Candidate evidence ladder

The search stages are deliberately asymmetric:

1. cheap exact congruence and height checks;
2. exact minimal conductor and root number;
3. local-score and point-search heuristics;
4. verified rational points and a height-pairing independence certificate;
5. descent/Selmer upper bounds if exact rank is claimed;
6. independent replay before a record claim.

The low-rank artifact reaches the full conductor-calibration level.  The first
Fermigier-family seed reaches exact high-family local shaping but stops before
factoring its uncontrolled 469-bit cofactor, computing a global conductor, or
certifying specialized points.  Neither result establishes that prime-power
shaping raises the probability of exceptional rank.  The next statistical
experiment is a height-matched comparison against unshaped specializations,
with good-prime scores and cofactor smoothness recorded separately.

Concretely, EC-FSEED1 uses the canonical adapter coordinate \(u=s/2\), lifts
the two simple split roots of its primitive degree-20 discriminant factor at
each of \(89^2,131^2,137^2\), and exhausts all eight CRT classes.  The first
height is attained by

\[
u=673709/29965,\qquad M=2551312982089.
\]

The three exact factor valuations are two.  Residue classification and PARI's
local Euler coefficients give split \(I_2\), Tamagawa number two, and conductor
exponent one at each shaped prime.  Its stored uncontrolled cofactor has 469
bits and is not factored, so the global conductor and rank remain open.
