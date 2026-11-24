# Modulus size for the Eisenstein divisor mark

At \(N^{1/5}\) roughness, the positive
[value-normalized proper-divisor mark](EISENSTEIN_ROUGH_FACTOR_DISTRIBUTION.md#what-the-normalization-changes)
can require testing a prime as large as \(N^{4/5}\).
A replacement based on one- and two-prime divisibility tests would
be easier to distribute if it only used moduli within the existing
square-root level. The best square-root mark in the class below
has a negative margin on almost all targets. A correction
using pair products up to \(a^{0.620}\) gives a small positive
margin. Its extra credit is available when the unmarked
cofactor has to be prime by roughness. The signed refinement
below shortens the modulus range further and locates the
roughness where square-root marks become viable.

## Divisibility marks on one summand

Use ordered pairs \(a+b=N\) with \(N\equiv2\bmod6\), \(a,b>1\)
and \(ab\) squarefree, with every prime factor \(1\bmod6\) and
strictly greater than \(N^{1/5}\). Count each norm-value pair once.
Write \(R_N\) for the total and \(G_N\) for the prime pairs.
Each summand has at most four prime factors.
For any function \(B(a)\) that vanishes on primes and is at
least one on composites in this family,
\[
 G_N\ \ge\ R_N-2\sum_{a+b=N}B(a).                             \tag{1}
\]
Each pair with at least one composite is charged at least one
by \(B(a)+B(b)\), and symmetry gives (1).

Use the exponents \(x_p=\log p/\log a\) for \(p\mid a\).
They sum to one and are strictly greater than \(1/5\).
Let \(K(a)=\#\{p\mid a:x_p<1/2\}\), and let \(V_\theta(a)\)
sum the pair weights
\[
 h_\theta(s)=
 \begin{cases}
 1,&s<1-\theta\ \text{or}\ 1/2<s\le\theta,\\
 0,&\text{otherwise}
 \end{cases}
 \quad\text{at }s=x_p+x_q
 \quad(p<q,\ x_p,x_q<1/2),
 \qquad \tfrac12\le\theta<\tfrac23.
\]
The equality \(x_p+x_q=1/2\) is impossible for a squarefree
value, and otherwise \(h_\theta(s)+h_\theta(1-s)=1\).
Define \(B_\theta=K-V_\theta\).
The single marks have \(p<\sqrt a\) and pair marks have
\(pq\le a^\theta\le N^\theta\).

This majorizes the composite indicator. A prime has no marks.
A semiprime has one small prime and no pair mark.
A triple with a prime larger than \(\sqrt a\) has two small
primes and at most one pair mark, so \(B_\theta\ge1\).
For a triple with all three primes below \(\sqrt a\), each
pair has \(s>1/2\), and a mark is allowed exactly when the
omitted exponent is \(\ge1-\theta>1/3\). At most two marks
are possible. For a four-factor value, every prime lies
below \(\sqrt a\) because the other three exponents sum to
more than \(3/5\). The six pairs form three complementary
pairs of pairs, yielding \(V_\theta=3\) and \(B_\theta=1\).

There is a second mark \(\widetilde B_\theta\) for
\(3/5\le\theta<2/3\). Give one extra credit to a marked
pair with
\[
 3/5 < x_p+x_q\le\theta,\qquad x_p,x_q<1-\theta.              \tag{2}
\]
Namely subtract that credit in the definition of \(B_\theta\).
The cofactor \(a/(pq)\) is \(<a^{2/5}\), so it cannot have
two prime factors \(>N^{1/5}\). It is greater than one, as
\(p,q<\sqrt a\), and therefore prime. There are three
prime factors, all below \(\sqrt a\). The cofactor is the
unique factor whose exponent is at least \(1-\theta\).
The pair in (2) is consequently the only mark of \(V_\theta\)
on this triple. The extra credit lowers its mark from two to
one. Thus \(\widetilde B_\theta\) is still a majorant using
only the same two-prime moduli. For instance on \(a=7\cdot13
\cdot19=1729\), which can occur at \(N=1760\), \(B_{1/2}=3\).
At \(\theta=31/50\), the pair \(7,13\) meets (2), giving
\(\widetilde B_{31/50}=1\).

## The square-root obstruction and the longer correction

**Proposition (optimal square-root mark).** Consider a mark
formed by summing a function of one marked exponent \(x_p<1/2\)
and a function of two marked exponents with \(x_p+x_q<1/2\).
The mark must vanish on prime profiles and dominate one on
every composite exponent profile with entries \(>1/5\), total
one and no proper subsum equal to \(1/2\). Squarefreeness
forces this last restriction. The optimum in the class is
\(B_{1/2}\), pointwise on all such profiles.

Indeed, any \(x\in(1/5,1/2)\) can be the small exponent on
the semiprime profile \((x,1-x)\). Thus the corresponding
single coefficient is at least one. On a balanced triple,
meaning all three exponents below \(1/2\), no pair lies in
the modulus range. Its charge is at least three. On all
other composite profiles a permitted majorant must give
at least one. The explicit \(B_{1/2}\) attains all these
lower bounds simultaneously. The same argument at any
fixed outer scale holds with \(1/5\) replaced by the
actual lower exponent \((\log N)/(5\log a)\) when testing
profiles for a first summand near that scale. The
semiprime complement \(1-x\) still meets this cutoff.
Values \(a<X/(\log X)^4\) give only \(O(X/(\log X)^4)\)
pairs at a prescribed target \(N\in(X,2X]\), so omitting
them cannot repair the limiting margin. In particular,
moduli \(d\le N^{1/2-\varepsilon}\) with a fixed
\(\varepsilon>0\) fall in the square-root range on the
macroscopic summands.

**Theorem (margins for the pair marks).** For almost all
\(N\equiv2\bmod6\), the ratio of the right side of (1) for
the marks above to \(R_N\) tends to the respective constants
in the following table. Each row is a rational enclosure
of the corresponding constant. For any fixed tolerance
and any \(A>0\), the ratios are within that tolerance
outside \(O_A(X/\log^A X)\) targets in \((X,2X]\), with the
constant also depending on the tolerance and the mark.

| Mark | Largest pair exponent | Limiting margin divided by \(R_N\) |
| --- | ---: | ---: |
| \(B_{1/2}\) | \(0.500\) | \([-0.064637,-0.063883]\) |
| \(B_{629/1000}\) | \(0.629\) | \([-0.001198,-0.000427]\) |
| \(B_{631/1000}\) | \(0.631\) | \([0.000628,0.001403]\) |
| \(\widetilde B_{617/1000}\) | \(0.617\) | \([-0.001316,-0.000563]\) |
| \(\widetilde B_{31/50}\) | \(0.620\) | \([0.002523,0.003283]\) |
| \(\widetilde B_{5/8}\) | \(0.625\) | \([0.008593,0.009364]\) |
| \(\widetilde B_{16/25}\) | \(0.640\) | \([0.024435,0.025246]\) |

In particular, exact counts at the square-root level already give
the wrong sign for every mark covered by the proposition. The
obstruction is specific to one-sided divisor majorants on a
single summand with coefficients based on the marked sizes.
The longer correction uses information about the two marked
prime sizes jointly and gives a target main term of the right
sign.

### Factor shapes and the constants

For an exponent vector of degree \(j\) use the finite measure
\[
 d\mu_j(\mathbf x)=\frac{dx_1\cdots dx_{j-1}}
                       {j!\,x_1\cdots x_j},
 \quad x_j=1-\sum_{i<j}x_i,\quad x_i>1/5 .
\]
Its mass is \(c_j=c_j(1/5)\) from
[the factor-count distribution](EISENSTEIN_ROUGH_FACTOR_DISTRIBUTION.md#the-factor-count-distribution),
with the prime profile having unit mass. For a symmetric
bounded shape function with finitely many hyperplane cuts, its
first-summand limiting mass relative to total value mass on the
good targets is
\[
 \frac1T\sum_{j=1}^4 2^{1-j}\int g(\mathbf x)\,d\mu_j,
 \qquad T=1+\frac{c_2}{2}+\frac{c_3}{4}+\frac{c_4}{8}.         \tag{3}
\]
Here the full sum over the second summand cancels on division
by the total. The earlier fixed-factor proof gives (3)
by refinement to finitely many fixed boxes in \(\log p/\log X\).
To make that refinement, specify one box for the \(j\)
labelled prime factors. Split it into at most \(j!\) orderings
for the minor arcs so that the least factor is the free
prime and the product of the other factors is the grouped
variable. The boxes restrict the grouped variable and
intersect the interval tests for the free prime. The
same Cauchy–Schwarz estimate therefore applies. On major
arcs use labelled prime tuples with the symmetry factor
\(1/j!\). Siegel–Walfisz applies in the fixed prime intervals
with the product cutoff, giving the same shared local
multiplier. After the change to the simplex, both real
convolution arguments can be restricted to
\([X/(\log X)^4,2X]\) at cost \(O(X/(\log X)^5)\) in the
continuous model. On this range \(\log p/\log X\) and
\(\log p/\log a\) differ by \(O(\log\log X/\log X)\).
Inner and outer unions of fixed boxes therefore bracket
each of the shape sets. The integral of boxes near a cut
tends to zero with the mesh, since the density is bounded
on the rough simplex. For any fixed tolerance choose the
mesh first, then take the union of the finitely many
exceptional sets. Finally the two frozen roughness cutoffs
and deletion of common primes work as in the factor-count
proof. This supplies the claimed exceptional-set form
for (3) without a new uniform estimate as the mesh shrinks.

For triples, put
\[
 g(s)=\frac{\log((1-s-1/5)/(1/5))}{s(1-s)},\qquad
 L=\int_{1/2}^{3/5}g(s)\,ds .
\]
The marginal density of one particular exponent in \(\mu_3\)
is \(g(s)/3\). This follows by integrating
\((6s x(1-s-x))^{-1}\) for \(1/5<x<1-s-1/5\).
Thus \(L\) is the triple mass with one large prime. For
\(b=1-\theta>1/3\), write
\[
 L_\theta=\int_{1/2}^{\min(\theta,3/5)}g(s)\,ds,\qquad
 K_\theta=\int_b^{1/2}g(s)\,ds.
\]
For \(b<s<1/2\), the other two exponents automatically lie
below \(1/2\), so \(K_\theta\) is the mass of eligible
omitted exponents on balanced triples, counted with
multiplicity. The triple contribution to \(B_\theta\) is
\(3c_3-2L+L_\theta-K_\theta\). The \(L_\theta\) term
charges a triple with large exponent \(1/2<s\le\theta\)
two units instead of one. The \(K_\theta\) term removes
one credit per qualifying pair on a balanced triple.
For \(\theta\ge3/5\) the extra credit (2) has mass
\[
 J_\theta=\int_b^{2/5}
               \frac{\log(b/(1-s-b))}{s(1-s)}\,ds .
\]
Indeed the cofactor exponent is \(s\in[b,2/5)\).
Both marked exponents must be below \(b\), so one lies
in \((1-s-b,b)\), with the other fixed by the sum.
The lower endpoint is at least \(1/5\) in this range.
Integrating the other two denominators as above gives
the displayed mass, as at most one cofactor can qualify.

Combining these triple masses with the charges on primes,
semiprimes and four-factor values gives the margin
\[
 \frac{1-c_2/2-5c_3/4-c_4/8
         +L-L_\theta/2+K_\theta/2}{T}                         \tag{4}
\]
for \(B_\theta\), with \(J_\theta/(2T)\) added for
\(\widetilde B_\theta\). For \(B_\theta\) the numerator is
strictly increasing in \(\theta>1/2\) up to \(2/3\).
Below \(3/5\), \(K_\theta-L_\theta\) has derivative
\(g(1-\theta)-g(\theta)>0\): the denominators agree and the
first logarithm is larger. Above \(3/5\), \(L_\theta=L\)
while \(K_\theta\) increases. The simple family's unique
zero lies between \(0.629\) and \(0.631\). For the improved
family the displayed rows certify a sign change between
\(0.617\) and \(0.620\).

## Implication for the pointwise estimates

The [Eisenstein joint coset estimate](EISENSTEIN_BOUNDED_PRIME_FACTORS.md#all-target-valuations-and-joint-cofactor-levels)
has a summed level approaching \(N^{1/2}\) at macroscopic
width for its stated tests. The square-root obstruction says
that reaching that level with these single-summand marks
would still leave a negative main term at \(N^{1/5}\)
roughness.

For the improved \(\theta=31/50\) choice, write \(L_N\) for
the single-mark total \(\sum K(a)\) and \(V_N\) for the
two-prime total including the extra credit in (2).
The sufficient pointwise bound is \(R_N-2L_N+2V_N>0\).
This needs a lower count for pair products up to \(N^{31/50}\)
on the same value-normalized rough shell. All error terms
must fit a combined main-term margin of about \(0.3\%\)
of the total. At \(\theta=16/25\) the error budget rises
to about \(2.5\%\), with a longer modulus range.
This is a quantitative trade-off between modulus and
allowable error. Distribution of primes in arithmetic
progressions to large moduli with factorable weights
concerns different coefficients and does not supply
these marked rough-shell counts.
[Ford–Maynard on prime-producing sieves](https://arxiv.org/abs/2407.14368) ·
[Lichtman on distribution beyond the square-root barrier](https://arxiv.org/abs/2309.08522)

## Reproduction

Run `python3 research/arithmetic/verify_eisenstein_divisor_modulus.py`
from the repository root with assertions enabled. It reuses
the exact simplex enclosures for \(c_2,c_3,c_4\).
For each one-dimensional integral it bounds the log argument
by its endpoints in 512 intervals, bounds \(s(1-s)\) by
its extrema there, and uses the positive atanh series with
a rational tail bound. All rounding for the integrals and
the printed table is outwards. At \(\theta=31/50\) the
verified margin is between \(0.002523\) and \(0.003283\).
The computation certifies the signs of the limiting
constants, given (3).

## Signed pair weights and the roughness trade-off

The preceding pair weights were nonnegative. Allowing a negative
weight on a short pair lets its complementary longer pair receive
more credit. This yields a lower modulus range at \(N^{1/5}\).
A separate comparison locates the roughness threshold at which
the optimal square-root mark itself has positive margin.

### Signed complementary weights at one fifth

For \(1/2\le\theta\le3/5\), retain \(K(a)\) and the form \(B=K-V\)
from (1), with \(V\) summed over pairs of small primes. Replace
the weight \(h_\theta\) by
\[
 h^\pm_\theta(s)=
 \begin{cases}
 1,&s<1-\theta,\\
 -1,&1-\theta\le s<1/2,\\
 2,&1/2<s\le\theta,\\
 0,&s>\theta.
 \end{cases}                                                \tag{5}
\]
The half-exponent equality is excluded by squarefreeness.
This still gives a composite majorant on the \(N^{1/5}\)-rough
family. The semiprime mark remains one. A triple with a large
factor has two small prime marks and either gets one unit of
credit or pays one unit for its short pair. A balanced triple
has three small primes and at most one eligible pair, because
two omitted exponents at least \(1-\theta\ge2/5\) would leave
the third exponent at most \(1/5\). Its mark is one or three.
For four prime factors, \(h^\pm_\theta(s)+h^\pm_\theta(1-s)=1\).
The three complementary pairings still leave a mark of one.

Relative to the square-root rule, (5) gains two units on a
balanced eligible triple and loses two units on each unbalanced
triple with large exponent in \((1/2,\theta]\).
Using \(L_\theta,K_\theta,T\) from (4), its limiting margin
as a fraction of \(R_N\) is
\[
 \frac{H_{1/2}+K_\theta-L_\theta}{T},\qquad
 H_{1/2}=1-c_2/2-5c_3/4-c_4/8+L .                           \tag{6}
\]
The argument using (3) gives the same fixed-tolerance
\(O_A(X/\log^A X)\) exceptional set. Since
\(g(1-\theta)>g(\theta)\) for \(1/2<\theta\le3/5\),
the numerator in (6) strictly increases.

**Proposition (optimal signed complementary weight).** Among
bounded piecewise-continuous pair weights depending only on \(s=\log(pq)/\log a\), supported
on \(s\le\theta\), and satisfying the complementary identity
\(h(s)+h(1-s)=1\) on the four-factor range \(2/5<s<3/5\),
rule (5) maximizes the limiting margin with the fixed single
mark \(K\).

For \(\tfrac12<s\le\theta\) the selected pair can occur on a
balanced triple with exponents close to
\((s/2,s/2,1-s)\). The other two pairs are outside the selected
range, so the composite constraint gives \(h(s)\le2\).
Each increase in this weight requires the opposite decrease
at \(1-s\). The resulting increase in the mean credit is
proportional to \(g(1-s)-g(s)>0\). The weights at
\(s>\theta\) and their complements are already forced. These
facts give the optimizer (5), with the endpoint cases
understood by continuity of the limiting integrals.
The scope is the stated complementary identity and the
fixed single mark.

The limiting sign turns positive between \(0.591\) and \(0.592\),
below the \(0.620\) unsigned correction. At \(0.592\) the
margin is small. At \(3/5\), the short positive band of
(5) is empty because two rough factors have product
\(>a^{2/5}\). The mark then uses a lower count for longer
pairs and an upper count for the shorter pairs that carry a
negative coefficient.

### When square-root marks suffice in the limiting model

Now use any fixed \(1/5\le\eta<1/2\) and retain the optimal
square-root mark \(B_{1/2}\). The shape proof (3) applies with
\(\eta\) in place of \(1/5\): the same box refinement, frozen
cutoffs, and deletion of shared factors work at every fixed
\(\eta>0\). Let \(c_{3,\mathrm{bal}}(\eta)\) and
\(c_{3,\mathrm{unbal}}(\eta)\) be the two triple masses under
the corresponding simplex measure. Each composite has
mark one except a balanced triple, which costs three.
The limiting margin numerator is therefore
\[
 1-\frac{c_2(\eta)}2
   -\frac{5c_{3,\mathrm{bal}}(\eta)
                 +c_{3,\mathrm{unbal}}(\eta)}4
   -\frac{c_4(\eta)}8,                                      \tag{7}
\]
and the positive normalization is
\(1+c_2(\eta)/2+c_3(\eta)/4+c_4(\eta)/8\).
For the numerical calculation,
\[
 c_{3,\mathrm{unbal}}(\eta)=
 \int_{1/2}^{1-2\eta}
     \frac{\log((1-s-\eta)/\eta)}{s(1-s)}\,ds
\]
when \(1-2\eta>1/2\), and is zero otherwise.
The coefficient masses of both triple shapes and of the
other composite types are continuous and nonincreasing as \(\eta\) rises,
by nesting of the cut simplices. The \(c_2\) term is strictly
decreasing. Thus (7) has at most one zero. In particular
the square-root obstruction from one-fifth roughness changes
sign between \(\eta=0.214\) and \(\eta=0.215\).

**Theorem (certified trade-off).** On the almost-all sets of
the factor-shape comparison, the limiting margins in
(6) and (7) have the following enclosures.

| Roughness exponent | Pair-product bound | Mark | Limiting margin divided by \(R_N\) |
| ---: | ---: | --- | ---: |
| \(0.200\) | \(a^{0.591}\) | Signed complementary | \([-0.001262,-0.000255]\) |
| \(0.200\) | \(a^{0.592}\) | Signed complementary | \([0.000186,0.001200]\) |
| \(0.200\) | \(a^{0.595}\) | Signed complementary | \([0.004639,0.005672]\) |
| \(0.200\) | \(a^{0.600}\) | Signed complementary | \([0.012428,0.013495]\) |
| \(0.214\) | \(a^{1/2}\) | Square-root | \([-0.002465,-0.001944]\) |
| \(0.215\) | \(a^{1/2}\) | Square-root | \([0.002057,0.002563]\) |
| \(0.220\) | \(a^{1/2}\) | Square-root | \([0.024834,0.025276]\) |

Run `python3 research/arithmetic/verify_eisenstein_modulus_tradeoff.py`
to verify the rational enclosures. It reuses the exact simplex
bounds and the interval integration from the previous checker.
Assertions also check two sample triple profiles against their
integer-power thresholds. The strict sign comparison between
the two roughness exponents concerns the limiting criterion.
It entails no estimate for an individual exceptional target.

At \(\eta=0.220\) the mark itself needs only square-root
moduli. The positive limiting sign moves the problem to
a pointwise lower estimate for the normalized rough shell
and upper and lower estimates for its marks. At this roughness
the unrestricted affine-pair lower sieve is within the
range \(\eta<1/4.267\) from
[the random-semigroup proof](RANDOM_SPLIT_PRIME_SUMSETS.md#rough-decompositions-with-controlled-supports).
It counts arbitrary rough pairs. Applying that dimension-two
lower-bound route directly to the current norm shell would
require \(D>N^{4.267(0.220)}\) for the stated sufficient
sifting ratio, well above its known summed Hecke level.
At \(\eta=1/5\), the signed rule keeps the original rough
support but asks for pair-modulus information out to
\(N^{0.592}\) or further for a larger margin. Thus changing
the mark controls which distribution input is needed, and
positive shape margins coexist with a pointwise analytic gap.

