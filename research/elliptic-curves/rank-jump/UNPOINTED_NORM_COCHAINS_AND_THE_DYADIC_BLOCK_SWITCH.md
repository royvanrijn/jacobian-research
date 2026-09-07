# A cochain without elliptic points, and its exact dyadic block switch

Follow-up: [the two cross-cochains now complete the full three-dimensional
Selmer pairing](FULL_GOVERNING_BLOCK_SEPARATES_SELMER_INCIDENCE_AND_CT.md).
Twists 41 and 113 have zero full CT at unchanged Selmer incidence; twist 97
has CT rank two and rank at most one. Rationality of the zero-pairing
classes remains UNKNOWN.

A cubic norm equation constructs the governing cochain for two strict
Selmer classes without requiring rational points on their elliptic covers.
The frozen small class-group control gives an explicit example:

\[
g(T)=T^8-18T^6+92T^4-112T^2+16,
\quad\operatorname{Gal}(g)=\mathrm{SL}_2(\mathbf F_3),
\quad\operatorname{disc}(g)=2^{36}163^4.
\]

The two input classes are known Sha classes on one elliptic curve and
known rational classes on its minus twist. The **same cochain and octic**
serve both. Their CT difference is a single dyadic Hilbert symbol, computed
below from the norm witness without elliptic points. This is a new direct
cochain derivation of an already certified switch, not a new rank result.

It advances the endpoint in
[the rational-pair octic note](EXPLICIT_PAIR_GOVERNING_OCTIC.md): general
strict classes no longer require representatives of the form x_P-theta.
Solving a norm equation remains a computational cost, and the construction
does not by itself detect whether a class is rational.

## General construction

Let K be a cubic field with labelled embeddings i=1,2,3. Let alpha,beta
be nonzero elements of square rational norm. Suppose a witness in K satisfies

\[
X^2-\alpha Y^2=\beta,\qquad Y\ne0.
\tag{1}
\]

Choose A_i^2=alpha_i and B_i^2=beta_i, with prescribed rational products
a_0 and b_0. Let

\[
C_i=X_i+B_i,\quad F=\prod_i C_i,\quad T=\sqrt F,
\quad q_v=\prod_{i:v_i=1}\frac{A_iY_i}{C_i}
\]

for even masks v in V subset F2^3. The C_i are nonzero since
C_i(X_i-B_i)=alpha_i Y_i^2 is nonzero. As usual, V is the two-torsion
module, with dot product equal to the additive Weil pairing.

Write sigma=(e,f,g) for its two sign masks and permutation of the embeddings.
Flipping B_i replaces C_i by alpha_i Y_i^2/C_i. It follows that

\[
\sigma(F)/F=q_f^2,\qquad
\sigma(q_v)q_f=(-1)^{e\cdot gv}q_{gv+f}.
\]

Thus, defining gamma by sigma(T)=(-1)^gamma(sigma) q_f T gives directly

\[
\boxed{d\gamma(\sigma,\tau)=e_\sigma\cdot g_\sigma f_\tau=e\cup f.}
\tag{2}
\]

There is no correction term in this version. The labelled eight roots
(-1)^s q_v T transform by

\[
(v,s)\longmapsto(gv+f,\ s+\gamma+e\cdot gv).
\tag{3}
\]

The independent verifier checks this cup-extension action on all 192^2
pairs in the universal S3 group. The formula remains valid when the actual
Galois image is smaller; the group must then be determined separately.

Put

\[
P=N_{K/\mathbf Q}(X),\quad
D=N_{K/\mathbf Q}(X^2-\beta),\quad
C=2D\operatorname{Tr}_{K/\mathbf Q}
\left(\frac{X^2+\beta}{X^2-\beta}\right).
\]

The resulting octic is

\[
\boxed{h(T)=T^8-4(P+b_0)T^6+CT^4-4D(P-b_0)T^2+D^2.}
\tag{4}
\]

Indeed, the four even-sign conjugates of F have elementary symmetric
functions 4(P+b_0), C, 4D(P-b_0), D^2. The verifier expands all four
identities as sparse integer polynomials in the six independent variables
X_i,B_i. No numerical roots or elliptic points enter this check.

This gives a sufficient construction whenever (1) is soluble. It is not
claimed that every arbitrary pair of non-strict Selmer representatives
passes this particular norm equation.

## Why the norm gate always opens for strict classes

Let S contain the places above 2 and all archimedean places. Suppose both
classes are square at every place in S and have even valuations elsewhere.
At a place outside S, divide off squares to make both classes units. Over
a local field of odd residue characteristic the Hilbert symbol of two
units is trivial. At S one argument is already square. Hence
(alpha,beta)_v=1 everywhere over K.

The quadratic extension K(sqrt(alpha))/K is cyclic, so the
[Hasse norm theorem, Milne VIII.3.1](https://www.jmilne.org/math/CourseNotes/CFT.pdf)
gives (1). This is an existence statement; a bounded norm solver can still
return UNKNOWN. The equation is therefore a universal construction gate
for these strict classes, **not a rational-solubility discriminator**.

For a strict class paired with an S-unramified norm-square class, the same
argument works: the strict argument is square at S and both valuations
are even outside S. This supplies a route to cross-pairings against a
larger Selmer subspace, not only the strict block's internal matrix.

## The fixed class-selected control

The [protocol](UNPOINTED_GOVERNING_NORM_PROTOCOL.json) uses

\[
\begin{aligned}
K&=\mathbf Q(\theta),&\theta^3-11\theta^2-14\theta-1&=0,\\
\alpha&=\theta^2-10\theta+1,&
\beta&=\theta^2-13\theta+12.
\end{aligned}
\]

These are the existing two class-group generators of the strict space for
S={2,163,infinity}, each of norm 625. Their selection preceded the earlier
rational-point test. Their outcomes were already known when this protocol
was frozen, so this is a confirmatory construction experiment.

One norm computation, with a 60-second worker cap and fixed seed, gives

\[
X=\frac{37}{2}+21\theta-2\theta^2,
\qquad
Y=-\frac{66}{5}-\frac{51}{5}\theta+\frac9{10}\theta^2.
\tag{5}
\]

The exact identity (1) is independently replayed by rational arithmetic
modulo the cubic. The generator reads no elliptic points. Here

\[
P=19375/8,\quad b_0=25,\quad D=131675625/64,
\]

and (4) becomes

\[
h(T)=T^8-\frac{19575}{2}T^6+\frac{859156875}{32}T^4
-\frac{2524880109375}{128}T^2+\frac{17338470219140625}{4096}.
\]

If z satisfies g from the beginning of this note, the exact substitution

\[
T=-\frac{45}{32}z^7+\frac{195}{8}z^5
  -\frac{435}{4}z^3+\frac{135}{2}z
\]

satisfies h(T)=0. Its first eight powers span the entire quotient algebra,
so this is an isomorphism of the two degree-eight algebras. The independent
verifier checks the substitution, spanning determinant, and discriminants.
Consequently all finite ramification of the splitting field is contained
in {2,163}. No maximal-order or exact field-discriminant claim is needed.

## The smaller Galois block

The automorphism

\[
\tau(\theta)=\theta^2-12\theta-2
\]

has order three, sends alpha to beta, and satisfies

\[
\tau(\beta)=\alpha\beta r^2,\qquad
r=(-29-23\theta+2\theta^2)/25,\qquad \alpha\beta r=-25.
\]

The three conjugates of alpha therefore span just the two classes
alpha,beta. Their independence over K is certified by two independent
quadratic characters at the three roots modulo 37. The joint class field
has group V semidirect C3=A4 and degree 12, not the S3 control's degree 96.

Choose the B roots to be the cyclic shift of the A roots. On the kernel
over K, the sign masks obey f_i=e_(i+1). Restricting (2) to this subgroup
gives a nontrivial central extension: every nonzero translation has lifts
of order four, and two suitable translations have nontrivial commutator.
The kernel over K is Q8. The order-three action cycles its three pairs of
noncentral elements, giving Q8 semidirect C3=SL2(F3), of order 24.

The verifier constructs this subgroup inside (3), proves its faithful
transitive action on eight roots, and counts orders as
1,1,8,6,8 elements of orders 1,2,3,4,6, respectively. The nonzero octic
discriminant makes the labelled roots distinct, proving irreducibility and
the asserted splitting group. This independently confirms PARI's group
identification without applying the S3 independence theorem to a cyclic cubic.

At all 29 eligible inert primes at most 199, the octic's factor degrees
agree with the Legendre symbol of N(X+sqrt(beta)), where sqrt(beta) in
F_(p^3) has norm 25. There are 14 zero bits and 15 one bits. These are
replay counts, not rank frequencies.

## A local formula when the classes are strict

Let E and E^(d) share the labelled two-torsion module and let the two
classes belong to both Selmer groups. At every exceptional place where
their local cohomology classes are zero, choose local primitives for their
cocycles. After choosing compatible local square roots, both sign cocycles
are zero and the local quadratic character of the cochain is represented by

\[
F_v=N_{K\otimes\mathbf Q_v/\mathbf Q_v}(X+\sqrt\beta),
\quad N(\sqrt\beta)=b_0.
\]

Suppose outside a finite set S the curves have good reduction, the twist
and governing cochain are unramified, and the two classes are unramified.
Then

\[
\boxed{\operatorname{CT}_{E}(\alpha,\beta)+
\operatorname{CT}_{E^{(d)}}(\alpha,\beta)
=\sum_{v\in S}\operatorname{inv}_v(d,F_v).}
\tag{6}
\]

This is a deduction from the cochain calculation in the
[proof of Morgan Proposition 3.3, equation (3.6)](https://arxiv.org/pdf/2309.02374v2),
not an application of that proposition's stated locally-trivial-twist
hypothesis. If a=dP_v locally, use the local lift d(1 tensor P_v). It is a
coboundary and hence satisfies the required local Kummer condition. The
summand becomes inv_v(chi cup (P_v cup b+gamma)); in local zero-cocycle
coordinates this is the Hilbert symbol above. Outside S the cocycles and
local lift are unramified, so their degree-two unramified cohomology class
vanishes. This argument allows the twist to be nontrivial at S precisely
because the selected classes are locally zero there.

## The entire difference is dyadic in this control

Take

\[
E_0:y^2=x^3-11x^2-14x-1,
\qquad E_+=E_0^{(-1)}:y^2=x^3+11x^2-14x+1.
\]

The common strict block U=<alpha,beta> has already been proved entirely
Sha on E0 and entirely rational on E+ in
[the exact quotient comparison](NORM_LIFTS_CAN_BE_ENTIRELY_SHA.md).
Their full Selmer spaces have dimension three but different local boundary
lines; this experiment preserves U, not the entire Selmer group.

Evaluate (6) with d=-1 and S={2,163,infinity}:

* **At 2:** the cubic is irreducible and unramified. The unique residue
  square root of beta in F8 is 1+theta^2. For any actual local square root,
  2X+2sqrt(beta) is therefore congruent modulo 4 to
  39+42theta-2theta^2. Its norm is -79985, congruent to 3 modulo 4.
  Thus F_2 has valuation -3 and odd unit 3 mod 4, giving (-1,F_2)_2=-1.
  The higher-precision choice of root does not affect this residue.
* **At 163:** X and sqrt(beta) are integral, while X^2-beta has unit norm.
  Hence X+sqrt(beta) is a unit at every place above 163. Its norm is a
  unit, so (-1,F_163)_163=1.
* **At infinity:** alpha is positive at all three roots. Equation (1)
  implies |X_i|>|sqrt(beta_i)|. Consequently sign(F_infinity)=sign(NX)=+1,
  giving Hilbert symbol +1.
* **Elsewhere:** the good/unramified argument in (6) gives zero contribution.

The [dyadic certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_strict_cochain_dyadic_switch_v1.json)
checks these exact residues, norm valuations and real signs. The result is

\[
\boxed{\operatorname{CT}_{E_0}|_U+
\operatorname{CT}_{E_+}|_U=
\begin{pmatrix}0&1\\1&0\end{pmatrix}.}
\]

It agrees with the previously retained independent ideal-cup calculation.
The new contribution is an explicit cochain and a concrete local explanation
of the simultaneous two-class obstruction change. The previously known
rational witnesses on E+ still supply the rationality endpoint; a zero
pairing alone would not do that.

## Consequences for the rank-jump programme

1. **Solubility obstruction, strongest new result:** strict classes admit a
   point-free cochain construction by a cubic norm equation. Its local
   quadratic characters can explain an entire alternating obstruction block.
   The dyadic unit here supplies an explicit specialization-dependent change
   for the fixed minus twist, with no new point search.
2. **Weak explanation:** the same classes, norm witness, cochain field,
   ramification support and good-prime Frobenius table work on both the Sha
   and rational sides. Those objects alone are insufficient predictors.
   Local geometry and the pairing's baseline must also enter.
3. **Incidence versus obstruction:** norm solvability for strict classes is
   automatic; the new feature is the local evaluation of its cochain, not
   whether the norm solver succeeds. Solver runtime is not a rank feature.
4. **Missing production computation:** apply this construction to a small
   fixed set of class-selected production representatives, without reading
   exceptional points. Large cubic arithmetic may still be costly. A
   complete class group is not part of this pair-level endpoint once valid
   representatives and their strictness are supplied.
5. **Next bounded test:** compute the two cross-cochains from this strict
   basis to the retained non-strict Selmer generator. This tests the full
   matrix rather than only its restriction. Subsequently compare the local
   matrices under a predeclared deformation. The original R17/A1 parameter
   families are not scalar twists, so formula (6) must not be applied to
   their t variation without the appropriate comparison theorem.
6. **Remaining solubility theorem:** full CT-radical membership still permits
   higher-divisible Sha. An independent rational construction or higher
   descent is necessary before claiming a soluble block, and independence
   modulo the specialized generic subgroup is necessary before claiming a
   jump. No current Agent 1 selection-policy change follows.

## Reproducibility

The [norm artifact](../../artifacts/generated-results/elliptic-curves/rank_jump_unpointed_governing_norm_v1.json)
retains the deterministic PARI witness, exact octic and 29 prime values.
The [independent verifier artifact](../../artifacts/generated-results/elliptic-curves/rank_jump_unpointed_governing_norm_verification_v1.json)
binds the old class-selected source and checks rational identities, Galois
relations, polynomial reduction, ramification and all modular values using
only Python's standard library. No elliptic points are read by these scripts.

```bash
python3 elliptic-curves/rank-jump/verify_unpointed_governing_norm.py check
python3 elliptic-curves/rank-jump/strict_cochain_dyadic_switch.py check
```

The one norm worker completed in under a second. Both independent replays
complete in seconds. All outputs are immutable, with local checkpoints in
the ignored rank-jump directory. The active search, its outputs, shared
navigation, MATH_STATUS and STATUS are untouched.
