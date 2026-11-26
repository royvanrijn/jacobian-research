# Random interpolation for joint shallow residue compression

A near-full family of genuine abc triples cannot remain compressed
in two coordinates across a large logarithmic weight of primes.
Each prime can compress a different subset of the rows. The
conclusion also controls an averaged collision statistic over prime
powers. The proof constructs a low-degree polynomial through random
anchors and uses the complete-letter curve bound with all three
original radical constraints.

[Shallow residue registers on actual charts (SRL1–SRL5)](ABC_SHALLOW_RESIDUE_REGISTER_SATURATION.md)
· [Complete-letter curve count (ACG)](ABC_WIDE_COFACTOR_CONGRUENCES.md#an-absolute-curve-saving-with-arbitrary-fixed-budget-degrees)

## Interpolation against proper congruence hits

Use the polylogarithmic-degree version of (ACG). It provides
absolute constants \(c_*,\sigma_*,L_*>0\) such that for every
fixed \(0<L<L_*\), every sufficiently large \(X\), and every
nonzero rational polynomial in two complete coordinates of
degree at most \((\log X)^{c_*}\), the number of positive
primitive triples on its zero set with \(A+B=C\le X\) and
\(\operatorname{rad}(ABC)<C^L\) is at most
\[
         O_{L,\varepsilon}\left(X^{(1/3-\sigma_*)L+\varepsilon}\right).
\]
The constant is uniform in the polynomial, including
polynomials moving with the height. Substituting \(C=A+B\)
and permuting complete roles allows any two different roles.
We may decrease \(\sigma_*\) to have \(0<\sigma_*<1/3\).

Put \(n=X^L\) and let \(S\subseteq[1,X]^2\cap\mathbb Z^2\)
be the pairs \((A,B)\) of any actual family of these triples,
with
\[
                     N=|S|\ge n^{1/3-\sigma_*/4}.
\]
All original eligibility restrictions can be retained in
choosing \(S\). In particular, one can start on a fixed
shallow-register chart, with the original radical cap
and primitive additive law intact. No bound on multiplicity
in an individual complete coordinate is needed.

For a modulus \(q\), write \(\mu_q(z)\) for the proportion of
\(S\) in the pair residue class \(z\in(\mathbb Z/q\mathbb Z)^2\).
Choose an integer
\[
        4\le r\le (\log X)^{c_*},\qquad
        V_r=\binom{r+2}{2},\qquad m=\lfloor V_r/4\rfloor .
\]
An anchor sample consists of \(m\) independent uniform draws
from \(S\). For a query \(x\in S\), a proper hit at \(q\)
means that some sampled anchor \(a\) satisfies \(a\ne x\)
and \(a\equiv x\bmod q\) in both coordinates. Its probability
is exactly
\[
     P_{q,m}(x)=1-\left(1-\bigl(\mu_q(x\bmod q)-1/N\bigr)\right)^m .
\]
This counts only distinct pairs and remains valid when
the query itself occurs among the sampled anchors.

**JSI1.** Let \(\mathcal Q\) be any finite set of distinct
prime powers. Then
\[
 \sum_{p^j\in\mathcal Q}\log p\,
           \frac1N\sum_{x\in S}P_{p^j,m}(x)
       \le \left(3r+O_L(m n^{-\sigma_*/2})\right)\log X .
 \tag{JSI1}
\]
The same bound holds for a finite set of pairwise coprime
integer moduli \(q\ge2\) using weight \(\log q\) at \(q\).
Prime powers can be taken at arbitrary depths.
Every proper hit on a modulus larger than \(X\) has
zero probability.

To prove this, fix all \(m\) anchors. The space of polynomials
of total degree at most \(r\) has \(V_r\) monomials. Take
coefficient vectors from \(\{0,\ldots,\lfloor X^r\rfloor\}^{V_r}\)
and map them to their integer evaluations at the anchors.
There are at most \((1+V_r X^{2r})^m\) output vectors. At
sufficiently large heights this is less than \(X^{r V_r}\),
since \(4m\le V_r\) and \(1+V_r X^{2r}\le X^{3r}\).
The domain has more than \(X^{r V_r}\) elements.
Subtract two coefficient vectors with the same evaluations.
This gives a nonzero integer polynomial \(F\) vanishing on
the anchors, with \(|F(x)|\le V_r X^{2r}\le X^{3r}\)
throughout the height box. These bounds are uniform over
the allowed \(r\) and include samples with repeated anchors.

For a fixed query \(x\), every hit at \(p^j\) makes
\(p^j\mid F(x)\). The sum of the weights \(\log p\) on all
hit depths of one prime is at most the deepest hit depth
times \(\log p\). Therefore the exponential of the total hit
weight divides \(F(x)\), even when different anchors supply
the hits at different powers. For the pairwise coprime
version use the product of the hit moduli instead.
In either setting the total weight for a fixed query
is at most \(m\log X\). Indeed each distinct anchor–query
pair contributes weight at most the log of the positive gcd
of its two coordinate differences, which is at most \(\log X\).
If one difference is zero, use the absolute value of the
other in this gcd. Summing over anchors can overcount hits
and hence is an upper bound.

If the total hit weight exceeds \(3r\log X\), the polynomial
vanishes at \(x\). Taking the \(\varepsilon\) of (ACG) to be
\(L\sigma_*/4\), the fraction of such queries is
\(O_L(n^{-\sigma_*/2})\). The mean hit weight at the fixed
anchor sample is at most \(3r\log X\) plus this fraction
times \(m\log X\). Average over the anchor samples to obtain
(JSI1). The bound is a statement at each fixed \(L\) before
any further small-budget limit.
[Polylogarithmic-degree form of (ACG)](ABC_WIDE_COFACTOR_CONGRUENCES.md#degrees-bounded-by-a-power-of-the-logarithm-of-height)

## Large subsets may depend on the modulus

**JSI2.** For every \(p^j\in\mathcal Q\), suppose \(q=p^j\ge m\)
and suppose a subset \(S_q\subseteq S\) with \(|S_q|\ge\eta_q N\)
occupies at most \(K_q q\) pair classes modulo \(q\), with
\[
         K_q\ge1,\qquad 0<\eta_q\le1,\qquad K_q q\le \eta_q N/2.
\]
Then at sufficiently large heights at each fixed budget,
\[
       \sum_{p^j\in\mathcal Q}
          \frac{\eta_{p^j}^{\,2}}{K_{p^j}}\frac{\log p}{p^j}
                       \le \frac{4096\log X}{r}.
 \tag{JSI2}
\]
For pairwise coprime integer moduli the same assertion uses
\(\log q/q\) at \(q\). The subsets and the classes can vary
arbitrarily with \(q\).

For a proof, omit singleton classes of the full set \(S\)
from the at most \(K_q q\) classes occupied by \(S_q\).
They account for at most \(K_q q\le\eta_q N/2\) rows.
The remaining classes have total \(\mu_q\)-mass at least
\(\eta_q/2\) and off-diagonal mass
\(\mu_q(z)-1/N\ge\mu_q(z)/2\) in each class. The probability
of a proper hit at \(q\) is consequently at least
\(\frac12\min\{1,m\mu_q(z)/2\}\) on these classes, using
\(1-(1-t)^m\ge1-e^{-mt}\ge\frac12\min\{1,mt\}\).
If classes of mass at least \(2/m\) account for at least
\(\eta_q/4\), their contribution to the overall query
probability is at least \(\eta_q/8\). Otherwise the lighter
classes carry at least \(\eta_q/4\), and their contribution
is at least \(m/4\) times the sum of their squared masses.
Cauchy–Schwarz bounds this from below by
\(m\eta_q^2/(64K_q q)\). The heavier case gives at least
the same bound since \(q\ge m\).

Insert these lower bounds into (JSI1).
For \(r\ge4\) we have \(m\ge r^2/16\), and the error
\(O_L(m n^{-\sigma_*/2})\) is at most \(r\) at sufficiently
large heights. Thus \(64(4r)/m\le4096/r\), as claimed.

In particular fix \(0<\beta<1/3-\sigma_*/4\) and
\(Q=n^\beta\). For the primes \(p\le Q\) where an
\(\eta N\)-row subset occupies at most \(Kp\) classes,
assume \(KQ\le\eta N/2\). Their total weight is at most
\[
          \frac{4096 K\log X}{\eta^2 r}
                            +O((\log(2m))^2).
\]
The last term bounds the primes \(p<m\) by summing over all
integers up to \(m\). Thus a positive fraction
\(\kappa\log Q\) of weighted primes cannot have this property
when \(r\ge 8192K/(\eta^2\kappa\beta L)\) and the small-prime
term is smaller than the remaining margin. Equivalently
take \(r=\lfloor(\log X)^{c_*}\rfloor\). If
\(K/\eta^2\le(\log X)^{c_*-\xi}\) for fixed \(\xi>0\),
the proportion of direct prime weight goes to zero at a
power of \(\log X\) for fixed \(\beta,L\).
This bound allows slowly growing \(K\) and decreasing \(\eta\).

## A distinct-pair collision consequence

**JSI3.** For any selected family of prime powers at which
\(\max_z(\mu_q(z)-1/N)\le1/m\), one has
\[
    \sum_{p^j\in\mathcal Q}\log p\,
       \left(\sum_z\mu_{p^j}(z)^2-\frac1N\right)
                   \ \ll\ \frac{\log X}{r}.
 \tag{JSI3}
\]
At each such modulus the proper-hit probability is at least
\(m/2\) times the off-diagonal residue mass seen by the query.
Its mean is at least \(m/2\) times the parenthetical
quantity. Apply (JSI1) and \(m\ge r^2/16\).
The term in parentheses subtracts the diagonal exactly.
For example, when every prime on a specified interval meets
this occupancy hypothesis, the mean logarithmic product of
the distinct interval primes dividing both coordinate
differences of an ordered distinct pair is
\(O(\log X/r)\). Conditioning on distinct pairs only
rescales the average by \(N/(N-1)\). The same conclusion
over nested powers counts their successive depths in the gcd.

### Deep secant gcds near zero complete multiplicity

**JSI5.** On a regular near-balanced chart of exponent
\(d=1/3\), take any fixed \(0<\zeta<d/2\) and a common
complete-degree exponent \(u<\zeta/4\). Choose the chart,
raw-cap and counting allowances smaller than the
strict margin. For two distinct uniformly sampled
rows let
\[
              g=\gcd(A-A',B-B')
\]
be the positive gcd. Then for every allowed interpolation
degree \(r\) in JSI1,
\[
   \mathbb E_{\rm distinct}
     \sum_{\substack{p^j\mid g\\p^j\ge n^\zeta}}\log p
                          \ll_L \frac{\log X}{r}.
\]
In particular the large-prime part supported on
\(p\ge n^\zeta\) has logarithmic size \(o(\log n)\)
on average when \(r\) grows as a small fixed power of
\(\log X\). The statement also limits deep powers of
smaller primes, with only their layers above \(n^\zeta\)
counted.

To verify the occupancy hypothesis in JSI3, fix a pair
residue class at \(p^j\ge n^\zeta\). On every nonempty
class at least one of the two coordinates is a unit at
\(p\), by original primitivity. Its complete residue is
specified. In its exact-layer box the SMI1 list count
with \(K=1\) bounds its dictionary by
\(n^{d-\zeta/2+o(1)}\), since the moment gain has already
reached \(\zeta\) at that half depth and is nondecreasing.
This estimate holds uniformly even for macroscopic
layers and when the other coordinate is a nonunit.
The complete-degree upper bound multiplies it by at
most \(n^{u+o(1)}\). Thus every pair residue class has
mass at most \(n^{u-\zeta/2+o(1)}\), power-below \(1/m\)
for the polylogarithmic anchor size. Sum JSI3 over
all the prime powers in the displayed range, and
condition on distinct rows. The modulus upper bound
\(q\le X\) is sufficient since larger depths have
zero proper-hit probability. The chart is taken
within the raw budget of JSI1, or with its fixed
arbitrarily small budget enlargement and the
corresponding allowances. This bound controls
unconditioned secants of the regular chart. Short
cells selected for particularly high gcds would
need a new count.
[Half-depth moment saturation for arbitrary profiles (SMI2)](ABC_SHALLOW_MOMENT_INTERPOLATION.md#full-half-depth-occupancy-and-a-deeper-bound)

## Comparison and the remaining arithmetic hypothesis

Walsh's Theorem 2.4 starts from a single complete point set
whose image modulo each prime in a set of positive
logarithmic weight has at most a fixed multiple of \(p^k\)
classes. It obtains a bounded-degree polynomial on almost
all the points. The argument above needs only an averaged
proper-hit statistic in the complete-letter plane and
produces enough polynomial zeros to use (ACG). This allows
a different dense concentrated subset at every modulus and
makes the degree budget explicit when the class multiple
or the density varies with height. The method is a
random-anchor interpolation argument in the setting of
Walsh's characteristic sets. Its arithmetic consequence
uses the uniform curve exclusion for genuine primitive
abc triples. Related fixed-degree work in this repository
uses large least common multiples of anchored secant gcds
to force exact vanishing. Here a random sample of
polylogarithmically many anchors controls the average
proper-hit weight across a prime family.
[Walsh, Proposition 2.2, Theorem 2.4 and §4](https://arxiv.org/pdf/1105.1551)
· [Bounded-degree curves and local secant pools (BC2–BC3)](ABC_BOUNDED_DEGREE_CURVES_AND_SECANT_GCD.md#local-prime-pools-for-common-secant-gcds)

The direct marked-prime events alone have too little weight.
If \(\eta_p=\Pr_{(A,B)\in S}(p\mid ABC)\), then
\(\sum_p\eta_p\log p\le L\log X\) by the original radical cap.
At each \(p\) these rows lie on the union of
\(A=0,B=0,A+B=0\bmod p\), so their joint image has at most
\(3p\) classes. Yet splitting at \(p=\log X\) gives
\[
      \sum_p \eta_p^2\frac{\log p}{p}\ll\log\log X+L .
\]
For the initial segment use the standard prime weight
\(\sum_{p\le y}(\log p)/p=\log y+O(1)\). Above the split
use \(\eta_p^2\le\eta_p\) and the raw radical sum.
A single set of fixed locks has the same limitation.
Indeed \(\sum_{p\mid D}(\log p)/p\ll_B\log\log X\) when
\(D\le X^B\), by the same split. Shallow truncations at
finitely many actual locks therefore do not directly
supply enough weighted joint compression. A route to a
strict unrestricted bound would still have to transfer
their actual joint lifting constraints to a much larger
family of prime tests, or force a different consequence of
the joint congruences.
[Actual shallow Latin charts (SRL5)](ABC_SHALLOW_RESIDUE_REGISTER_SATURATION.md#a-three-register-partial-latin-square)

## Limits of the direct p-adic projection transfer

Shen's restricted-graph theorem asks for \(|\mathcal A|\le
\delta^{-\alpha}\) with \(\alpha<1\) in
\(\mathbb Z/p^N\mathbb Z\), for a
\((\delta,\beta,\delta^{-\chi})\)-set \(\mathcal B\), and for
a sufficiently nonconcentrated unit coefficient set. For a
selected coefficient \(c\), it guarantees expansion for
every graph of density at least \(\delta^\chi\) in
\(\mathcal A\times\mathcal B\). Here \(\delta=p^{-N}\)
and \(0<\beta\le\alpha<1\). The following elementary obstruction
bounds the scope of the density threshold in that statement.

**JSI4.** The exponent in Shen's Theorem 1.8 necessarily
satisfies
\[
                          \chi\le \frac{\beta(1-\alpha)}2 .
\]
To see this, fix any \(\chi>\beta(1-\alpha)/2\) and choose
\[
   \max(0,\beta(1-\alpha)-\chi)<\chi'
                             <\min(\chi,\beta(1-\alpha)).
\]
Set \(k=\lceil(1-\alpha)N\rceil\).
Take \(\mathcal A=1+p^k(\mathbb Z/p^N\mathbb Z)\).
This has at most \(\delta^{-\alpha}\) entries.
Build \(\mathcal B\) by freeing \(f(j)\) of the first \(j\)
base-\(p\) digits, starting with the lowest-order digit.
For \(j\le k\) let \(f(j)=\max(0,\lfloor\beta j-\chi'N\rfloor)\).
From \(k\) to \(N\), take the integer part of the linear
interpolation from \(f(k)\) to \(f(N)=\lfloor\beta N\rfloor\).
Successive increments are zero or one for large \(N\),
since \(\chi'<\beta(1-\alpha)\le\alpha(1-\beta)\).
Give every free digit all \(p\) options and prescribe each
remaining digit. Prescribe the first digit to be one,
so \(\mathcal B\) consists of units. All occupied
\(p^j\) classes have relative mass \(p^{-f(j)}\), and
\(f(j)\ge\beta j-\chi'N-O(1)\). Thus \(\mathcal B\)
satisfies the required upper ball bound with factor
\(\delta^{-\chi}\) at all sufficiently large \(N\).
All units form one eligible coefficient set, since their
relative mass in a \(p^j\) class is \(O_p(p^{-j})\).

Choose one occupied \(p^k\) class \(\mathcal B_0\).
Its relative mass is \(p^{-f(k)}\ge\delta^\chi\).
For every unit \(c\), the graph
\(\mathcal A\times\mathcal B_0\) passes the density
test, and its projected values \(a+cb\) form a single
coset of \(p^k(\mathbb Z/p^N\mathbb Z)\) of size
\(|\mathcal A|\). This contradicts the proposed expansion
for that \(\chi\).
[Shen, Definition 1.5 and Theorem 1.8](https://arxiv.org/pdf/2607.12251v1#page=5)

For the direct complete-letter residue graph at zero
multiplicity, put \(d=1/3\), and give an external unit probe
\(p^N\) depth \(\sigma=\log_n p^N\).
In a fine-layer chart (SRL4) gives both original free
marginal size exponents \(\min(d,\sigma)\).
At \(\sigma\le d\) this rules out a fixed \(\alpha<1\)
for the full first marginal. At \(\sigma>d\), put
\(t=d/\sigma>0\). The two marginal sizes have exponents
\(t+o(1)\) in base \(\delta^{-1}\), and the graph of
projected original pairs has at most \(\delta^{-t-o(1)}\)
entries. Its density is at most \(\delta^{t-o(1)}\),
so it would need \(\chi\ge t-o(1)\).
The terminal scale of the ball condition implies
\(\beta-\chi\le t+o(1)\), even if one chooses \(\beta\)
larger than the natural marginal exponent. Hence (JSI4)
gives
\[
        \chi\le \frac{1-\alpha}{1+\alpha}\,t+o(1),
\]
which is incompatible with that density when
\(\alpha\ge t-o(1)\). This comparison concerns the
direct graph on the two full marginal alphabets.
A transfer based on a new choice of variables or a
dense arithmetic subgraph would require its own density
and image bounds.

There is a separate regularity issue at the sparse
terminal depth \(\sigma>d\). After the rank pruning
in (SRL4), a \(p^j\) ball contains at most
\(n^{d-\min(d,\log_n p^j)+o(1)}\) terminal
representatives. A subfamily of \(n^{d-\eta}\)
terminal representatives therefore occupies at least
\(n^{\min(d,\log_n p^j)-\eta-o(1)}\) of the \(p^j\)
balls. At depth \(\log_n p^j=d+o(1)\), this is
\(n^{d-\eta-o(1)}\). Upper covering regularity of
dimension \(t=d/\sigma\), with factor \(n^\varepsilon\)
in Shen's Definition 1.6(2), would bound this count
by \(n^{td+\varepsilon+o(1)}\), taking the parent ball
of radius one. Thus removing a subpower fraction in
the exponents cannot make these representatives upper
regular at their natural terminal dimension. The
upper pointwise ball estimate of (SRL4) itself
continues to hold.
[Shen, Definitions 1.5–1.6](https://arxiv.org/pdf/2607.12251v1#page=4)
· [Terminal shallow residues with rank pruning (SRL4)](ABC_SHALLOW_RESIDUE_REGISTER_SATURATION.md#several-actual-probes-at-zero-multiplicity)
