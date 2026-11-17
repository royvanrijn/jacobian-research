# Uniform candidate primes for finite obstructions of a given mass

The nonrational numerator branch contributes a third independent
point over a fixed quadratic field. A rank-three collision argument
therefore gives a stronger estimate at primes where a numerator
branch is actually supported. It also yields a finite, computable
set containing every new prime that can kill a smaller-prime prefix
with any of the four fixed core tests and a prescribed positive mass,
uniformly over these prefixes and all four classes. Common numerator
returns for two independent Haar vectors also have a quantitative
all-primes tail. This does not control the Hasse tail for one vector
or rule out finite obstructions with arbitrarily small prefix mass.

## 1. A rational height bound for three independent directions

Use the short model of the existing elliptic curve

\[
 E:y^2=X^3-10758382X+6600759744,\qquad X=x-3552,
\]

and the independent rational points `P,Q` of ETC3a, translated from
the model in
[the local classification](ELLIPTIC_EXCEPTION_COMPLETE_THREE_ADIC_CLASSIFICATION.md).
Write `Gamma=E(Q)` and
`R_p=im(Gamma -> E(F_p))`, `r_p=|R_p|`. Put

\[
 d=51730,\qquad L=\mathbb Q(\sqrt d),\qquad
 R=\left(\frac{25865}{128}-3552,
                    \frac{639735\sqrt d}{2048}\right).       \tag{1.1}
\]

Here `R` is the branch point called `R_+` in EKX, not the rational
point `P+T` used in other notes. EBT1 in
[the branch-twist note](ELLIPTIC_EXCEPTION_BRANCH_TWIST.md)
proves that `R` is nontorsion and conjugates to `-R` over `L/Q`.
Thus `P,Q,R` are independent modulo torsion: conjugating any
relation first makes its `R` coefficient zero, and independence
of `P,Q` handles the remaining two coefficients.

For a rational number `t`, let `h(t)` be its logarithmic projective
height, with `h(0)=0`. The rational canonical-height comparison
retained in
[EEC1](ELLIPTIC_EXCEPTION_COUNT.md) gives

\[
 h(X(uP+vQ))+h(y(uP+vQ))
                     \ll_E 1+u^2+v^2                         \tag{1.2}
\]

whenever the point is not the identity. The bound on `y` follows
also from its monic cubic equation, since
`2h(y)<=3h(X)+O_E(1)`.

The existing rational twist point is

\[
 Q_d=\left(d\left(\frac{25865}{128}-3552\right),
                      \frac{639735d^2}{2048}\right)
 \in E^{(d)}(\mathbb Q),
\]
\[
 E^{(d)}:Y^2=Z^3-10758382d^2Z+6600759744d^3.
                                                               \tag{1.3}
\]

Its twist isomorphism to `E/L` is
`(Z,Y) -> (Z/d,Y/(d sqrt(d)))`, sending `Q_d` to `R`.
For `w!=0`, write

\[
 [w]R=(B,b\sqrt d),\qquad
 B=Z([w]Q_d)/d,\qquad b=Y([w]Q_d)/d^2.
\]

The same rational height theorem, now on this fixed rational twist,
gives

\[
 h(B)+h(b)\ll_E 1+w^2.                            \tag{1.4}
\]

Constants include the fixed curve, points and twist; no height
theorem over a number field is imported.

For each nonzero integer triple `(u,v,w)`, set

\[
 D_{u,v,w}=uP+vQ+wR,\qquad
 X(D_{u,v,w})=r_{u,v,w}+s_{u,v,w}\sqrt d,
 \quad r_{u,v,w},s_{u,v,w}\in\mathbb Q.
\]

If `w=0`, or if `u=v=0`, the coordinate is rational, and
(1.2) or (1.4) bounds its height. Otherwise write
`U=uP+vQ=(A,c)` and `[w]R=(B,b sqrt(d))`. Both points are
nonidentity, and `A!=B`: equal coordinates would make
`[w]R=U` or `-U` rational, contradicting anti-invariance and
nontorsion. The short-model addition formula gives exactly

\[
 r=\frac{db^2+c^2}{(B-A)^2}-A-B,\qquad
 s=-\frac{2bc}{(B-A)^2}.                          \tag{1.5}
\]

Using `h(ab)<=h(a)+h(b)`, `h(a+b)<=h(a)+h(b)+log2`
and `h(1/a)=h(a)` in this fixed formula proves

\[
 h(r_{u,v,w})+h(s_{u,v,w})
                       \ll_E 1+u^2+v^2+w^2.                  \tag{1.6}
\]

Let the positive integer

\[
 j_{u,v,w}=\operatorname{lcm}
       (\operatorname{den}r_{u,v,w},
        \operatorname{den}s_{u,v,w})
\]

use reduced rational denominators, with `den(0)=1`. It satisfies

\[
 j_{u,v,w}X(D_{u,v,w})\in\mathbb Z[\sqrt d],\qquad
 \log j_{u,v,w}\ll_E 1+u^2+v^2+w^2.              \tag{1.7}
\]

In particular, if `D_(u,v,w)` reduces to the identity at a prime
of `L` above a good split rational prime `p`, then `p` divides
`j_(u,v,w)`. Otherwise both rational coefficients are integral
at `p`, and `sqrt(d)` is integral, contradicting the pole of `X`
at that reduction. This only uses the integral good-reduction model,
as in ESR1; it is unaffected by affine denominators of the
individual generating points.

## 2. The rank-three small-image product

For a good rational prime split in `L`, choose either square root
of `d` in `F_p` and put

\[
 H_p=\langle\bar P,\bar Q,\bar R\rangle
                   \subseteq E(\mathbb F_p),\qquad t_p=|H_p|.
\]

The other choice changes `bar R` to its negative and hence gives
the same subgroup. For `z>=2` define

\[
 M_z=\lfloor z^{1/3}\rfloor,\qquad
 J_z=\prod_{\substack{|u|,|v|,|w|\le M_z\\
                       (u,v,w)\ne(0,0,0)}}j_{u,v,w}.
                                                               \tag{2.1}
\]

This is an explicit positive integer; all its factors are determined
by rational group arithmetic on `E` and its fixed twist and by
(1.5).

**Proposition EUP1 (rank-three bound for supported branch images).**

\[
 \log J_z\ll_E z^{5/3},\qquad
 \{p\text{ good, split in }L:t_p<z\}
                         \subseteq\{p:p\mid J_z\}.            \tag{2.2}
\]

Consequently

\[
 \sum_{\substack{p\text{ good, split in }L\\t_p<z}}\log p
                 \ll_E z^{5/3},\qquad
 \#\{p\text{ good, split in }L:t_p<z\}
                 \ll_E\frac{z^{5/3}}{\log z}.                 \tag{2.3}
\]

The same weighted and cardinality bounds hold for the entire,
possibly larger set of prime divisors of `J_z`.

**Proof.** The `(M_z+1)^3>z` points
`uP+vQ+wR` with coefficients between zero and `M_z` are distinct,
by independence. If `t_p<z`, two have the same reduction.
Their difference is a nonzero triple in (2.1), whose denominator
clearing integer is divisible by `p`, by (1.7).
There are `O(z)` such triples and each has
`log j_(u,v,w)=O_E(z^(2/3))`. This proves (2.2) and bounds
`log rad(J_z)` by `O_E(z^(5/3))`. For the cardinality, separate
prime divisors at most `z`, of which there are at most `z`, from
those greater than `z`, each contributing at least `log z`.
The result is
`z+O_E(z^(5/3)/log z)=O_E(z^(5/3)/log z)`.
No full basis of `E(L)` is needed. \(\square\)

## 3. A stronger density-one bound for numerator returns

Let `V_p` be the intrinsic Haar event of positive valuation
`v_p(F)>0`, with `F=128x-25865`, on the rational closure used in
EFA2. Finitely many excluded base points have measure zero.
For a good active prime `p=2 mod3` outside the fixed constants,
`V_p` can occur only in the two finite zero residues. If either
one is supported by the rational reduction image, then both are,
and

\[
 H_p\subseteq R_p,\qquad t_p\le r_p.              \tag{3.1}
\]

This follows because the zero residues are `bar R` and `-bar R`,
and `R_p` is a subgroup containing `bar P,bar Q`. If `p` is
nonsplit in `L`, neither zero is rational over `F_p`, so `V_p`
is empty. Thus in every supported case

\[
 \nu_\infty(V_p)\le\frac2{r_p}\le\frac2{t_p}.     \tag{3.2}
\]

**Corollary EUP2 (numerator losses with exponent three fifths).**
For every fixed `A>0`, outside a set of rational primes satisfying

\[
 \#\{p\le X\text{ exceptional}\}
                    \ll_{E,A}\frac{X}{(\log X)^{1+5A/3}},      \tag{3.3}
\]

all sufficiently large good active primes satisfy

\[
 \boxed{\nu_\infty(V_p)
                     \le\frac{2(\log p)^A}{p^{3/5}}.}         \tag{3.4}
\]

The same bound applies to the mass removed from any finite-place
class event by forbidding positive valuations, uniformly over that
event. In particular it bounds the numerator contribution to local
failure. The total local failure can also have a pole contribution;
no exponent-three-fifths bound for that total is asserted.

**Proof.** Apply (2.3) with
`z=X^(3/5)/(log X)^A`. The threshold function is eventually
increasing, and
`z^(5/3)=X/(log X)^(5A/3)`, `log z~(3/5)log X`.
Thus (3.3) counts all split primes below their corresponding
threshold. Away from this exceptional set, either no finite zero
is supported, or (3.1)--(3.2) give (3.4). Intersecting `V_p` with
any event cannot increase its absolute mass. As in ESR2, the
exceptional set has density zero among rational primes, using the
elementary prime-count lower bound proved there. \(\square\)

## 4. Effective uniform candidates with the explicit four-prime core

Use the fixed test core

\[
 C_0=\{2,3,5,1889\}
\]

of EDC1--EDC2 in
[the dependency-core note](ELLIPTIC_EXCEPTION_DEPENDENCY_CLOSED_CORE.md).
For `0<eta<=1/4` put

\[
 z_\eta=1+2/\eta,\qquad
 \mathcal C_\eta=C_0\cup\{p:p\mid J_{z_\eta}\}.                \tag{4.1}
\]

For supplied rational `eta`, this is an explicit finite computation:
(2.1) uses rational group arithmetic and denominators, followed by
integer factorization. One may use any supplied rational lower bound
on a positive prefix coefficient. No unspecified numerical cutoff
or full Mordell--Weil basis is an input.

**Theorem EUP3 (effective candidates uniform in the prefix).** For
every class `a`, prime `p`, and finite prime set satisfying

\[
 S\subseteq C_0\cup\{q\text{ prime}:q<p\},
\]

one has

\[
 \boxed{\mu_{S,a}\ge\eta,\qquad
        \mu_{S\cup\{p\},a}=0
          \quad\Longrightarrow\quad p\in\mathcal C_\eta.}    \tag{4.2}
\]

The candidate set is independent of `S,a` and satisfies

\[
 \sum_{p\in\mathcal C_\eta}\log p\ll_E\eta^{-5/3},\qquad
 |\mathcal C_\eta|\ll_E
                     \frac{\eta^{-5/3}}{\log(2/\eta)}.        \tag{4.3}
\]

**Proof.** Every finite-place set is already feasible in classes
`1,10q_0`, by ECO1 and the certified-class theorem. Thus neither
class can satisfy the killing hypothesis. For either remaining
class `2,5q_0`, primes one modulo three and the bad prime `q_0`
are automatic, by EDC1. The four core primes are included in
`C_eta`. It remains to consider a good active prime `p>=11`,
`p!=1889,q_0`.

The actual prefix is measurable in the non-`p` primary factors.
Indeed the core dependencies are only `2,3,5,7,1889`, and none
of them equals this `p`. At a remaining smaller good prime `q`,
the reduction order has all prime divisors smaller than `q`, and
its formal coordinate is `q`-primary. EDC2 therefore supplies the
same fresh formal `p`-coordinate and exact residue-hazard kernel.
This applies to any subset of the core tests: we do not impose
missing core conditions or reduce the given prefix mass.

In particular, a larger core place such as 1889 may already be in
`S` when `p<1889`. It causes no problem after conditioning on all
non-`p` factors: its finite reduction quotient has order
`1890=2*3^3*5*7`, and its formal factor is pro-1889. This is the
core-adjusted conditioning of EDC2, rather than the original
numerically lower-prime filtration with an unspecified cutoff.

If the new test kills this positive prefix, all its reduction
support must lie on the finite zero residues. Every other supported
residue has positive conditional local survival. Uniform reduction
Haar measure on `R_p` then gives

\[
 \eta\le\mu_{S,a}\le\frac2{r_p}.                 \tag{4.4}
\]

At least one zero is supported, so `p` splits in `L` and (3.1)
holds. Consequently

\[
 t_p\le r_p\le2/\eta<z_\eta.
\]

EUP1 makes `p` divide `J_(z_eta)`, proving (4.2). The entire
divisor set, including any irrelevant extra primes, satisfies
(2.2)--(2.3). Adding the fixed four primes changes only the
curve-dependent constant. Since `z_eta` is comparable to
`1/eta` on the indicated range, this proves (4.3). \(\square\)

The permitted prefix may include a numerically larger core prime.
Every other nonautomatic prime must be smaller than the new one;
arbitrary later tests can constrain its primary factor, so this
hypothesis is essential. The theorem replaces a coset-dependent
exceptional set by one effective set uniform over all permitted
prefixes with a common mass lower bound. Membership in `C_eta` is
only necessary: no candidate is asserted to kill a prefix or form
a minimal obstruction. As `eta` decreases the set can grow, so
there is no mass-independent finite candidate set or proof of
all-finite-set positivity in either unresolved class.

The proof uses rational heights on the two fixed rational curves,
the exact addition identity (1.5), an explicit product of integer
denominators, and the dependency-core terminal-support argument.
It needs no finite-prime extrapolation, new moving family, full
Mordell--Weil basis or artificial checker.
## 5. Square-summable numerator marginals over all primes

Away from the fixed finite set of bad-reduction and branch constants,
let `A_p` be the supported finite-zero reduction event in the intrinsic
rational closure, and put `n_p=nu_infinity(A_p)`. Set `n_p=0`
when the finite zeros are not rational or are unsupported. Reduction
is uniform on `R_p`, and support of one zero implies support of its
negative, so

\[
 n_p=\begin{cases}2/r_p,&\bar R\in R_p,\\
                  0,&\text{otherwise}.
       \end{cases}                                             \tag{5.1}
\]

This definition may include all good primes, not just primes two
modulo three. The finitely many omitted primes do not affect the
following conclusions. At the other primes, `A_p` and the event
of positive valuation `F` differ at most on the measure-zero base
points.

**Corollary EUP4 (all-primes numerator summability).** For `0<t<=1`,

\[
 \#\{p:n_p>t\}\ll_E
                       \frac{t^{-5/3}}{\log(2/t)}.             \tag{5.2}
\]

Consequently

\[
 \sum_p n_p^s<\infty\quad\text{for every }s>5/3,
 \qquad\boxed{\sum_p n_p^2<\infty.}             \tag{5.3}
\]

More quantitatively, for sufficiently large `Y`,

\[
 \boxed{\sum_{p>Y}n_p^2\ll_E
                 \frac{1}{(\log Y)^{1/5}\log\log Y}.}         \tag{5.4}
\]

Two independent intrinsic-Haar vectors therefore have simultaneous
numerator returns at only finitely many primes, almost surely. The
probability of any such common return above `Y` is at most the
right-hand side of (5.4), with its fixed implicit constant.

**Proof.** If `n_p>t`, the branch is supported and
`t_p<=r_p<2/t`. EUP1 with `z=2/t` proves (5.2); its logarithm
is comparable to `log(2/t)` throughout the indicated range.
For `s>5/3`, nonnegative integration, first for finite prime sets
and then by monotone limits, gives

\[
 \sum_p n_p^s
   =s\int_0^1 t^{s-1}\#\{p:n_p>t\}\,dt<\infty.
\]

For the quantitative tail, if `n_p>0`, use EUP1 with `z=t_p+1`
to obtain `log p<=log J_z<<_E(t_p+1)^(5/3)`. Since `t_p>=1`
and `n_p<=2/t_p`, this implies
`n_p<<_E(log p)^(-3/5)`. Also, for small positive `t`,

\[
 \sum_{p:n_p\le t}n_p^2
  \le 2\int_0^t u\,\#\{p:n_p>u\}\,du
  \ll_E\frac{t^{1/3}}{\log(2/t)},
\]

by (5.2). All primes exceeding `Y` have
`n_p<=C_E(log Y)^(-3/5)`. Substituting this value of `t` in
the last estimate gives (5.4).

For independent vectors `omega_1,omega_2`, the probability of a
common return at `p` is `n_p^2`. The probability of a common
return at any prime exceeding `Y` is at most
`sum_(p>Y)n_p^2`, which tends to zero. Thus infinitely many
common returns have probability zero. No independence between
different primes is used. \(\square\)

This controls overlap between two independent Haar vectors. It
neither gives a bound on same-vector correlations across distinct
primes nor asserts convergence of `sum_p n_p`. Its uniform all-prime
pair tail is not a single-vector Hasse tail. In particular it
supplies no all-primes numerator-survival measure or asymptotic
for rational-target pairs. The critical exponent `s=5/3` is not
covered by this proof.
