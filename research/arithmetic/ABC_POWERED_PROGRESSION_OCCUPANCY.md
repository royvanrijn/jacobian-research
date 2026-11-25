# Small-radical integers in short progression blocks

Write \(\operatorname{rad}(v)\) for the product of the distinct prime
divisors of the positive integer \(v\). We count such values in a block
of a reduced arithmetic progression. The difference of two values
produces a primitive \(abc\) triple. In the derived radical budget the
cost of the original progression step is its radical, even when the
step is large. This gives a one-letter estimate independent of how the
values arose as sums. A second argument uses an adaptively chosen full
divisor and gives subpower occupancies for targets with a positive
diffuse valuation component.

## A close pair gives a primitive triple

We use two available small-budget counts. For each fixed
\(a\in(1/3,1/2)\), there is \(\lambda_0(a)>0\) such that
\(N_\lambda(Z)\ll_{a,\lambda,\epsilon} Z^{a\lambda+\epsilon}\)
for \(0<\lambda<\lambda_0(a)\). There are also absolute
\(h_*,\lambda_*,\sigma>0\), with \(0<\sigma<1/3\), for which the
same count holds with \(a=1/3-\sigma\) in the sector with a short
summand at most \(c^{h_*}\), for \(0<\lambda<\lambda_*\).
We decrease \(\sigma\) and \(h_*\) if needed so that \(a>1/4\) and \(h_*\le1\).
[Unrestricted small-budget bound](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)
· [Short-addend bound](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-strict-saving-when-the-short-term-is-small)

Fix nonnegative \(\rho,\mu,\tau\), put \(Y=X^\tau\), and require
\(B=2\rho+\mu+\tau\) to lie strictly below the available budget
threshold and below \(1/2\). Let \(M\ge1\) be an integer with
\(\operatorname{rad}(M)\le X^\mu\). Its size may vary with \(X\).
If \(a<1/3\), require also \(MY\le X^{h_*}\).
Then a block of parameter width at most \(Y\) in any unit class
modulo \(M\) contains at most
\[
                      O_{\rho,\mu,\tau,a,\epsilon}
                                (X^{a B/2+\epsilon})                   \tag{PI1}
\]
integers \(X<u\le2X\) with \(\operatorname{rad}(u)\le X^\rho\).
The bound is uniform in \(M\), the class and the block.

Indeed, for any pair \(u<v\) of eligible elements in the block put
\(g=\gcd(u,v)\) and \(j=(v-u)/M\). Since the class is a unit,
\((g,M)=1\), hence \(g\mid j\) and \(g\le Y\). The normalized
positive primitive triple is
\[
                 \frac{u}{g}+\frac{M j}{g}=\frac{v}{g}.
\]
Its radical is at most \(X^{2\rho}\operatorname{rad}(M)Y/g
\le X^B/g\le (X/g)^B\). If \(MY\le X^{h_*}\), its second
summand is at most \(X^{h_*}/g\le(v/g)^{h_*}\).
For \(M=1\) these facts just use the ordinary size of the gap. For a
general progression the gap may be tall and only its radical is
charged to the unrestricted estimate.

The derived triple determines \(u/g\) and \(v/g\). All eligible
scalings \(g\) lie in one class modulo \(M\), since
\(g(u/g)\) is in the fixed unit class. If two scalings occur, their
first values differ by at least \(M(u/g)>MX/Y\), which exceeds
\(MY\) because \(\tau<1/2\). Thus the map is injective on pairs
inside one block. The absolute height of each derived triple is
at most \(2X\). Use \(N_{B+\delta}(2X)\), with a sufficiently
small fixed \(\delta>0\) absorbed in \(\epsilon\), to bound the
number of pairs. This proves (PI1). It also includes a possible
singleton. Bounds for a finite collection of fixed \(B\)'s, with
each positive margin below its threshold, hold simultaneously.

For instance at \(M=1,\tau=\rho>0\), sufficiently small and fixed,
the short-addend bound gives \(O(X^{(1/2-3\sigma/2)\rho+\epsilon})\)
in every interval of length \(X^\rho\). This is a power saving over
the interval length and the ordinary total radical-cap count,
without a hypothesis on the size of its prime factors. The same
calculation in a reduced progression of parameter width \(X^\rho\)
and step with \(\mu\le\rho\) gives \(O(X^{(2/3+2\eta)\rho+\epsilon})\)
using \(a=1/3+\eta\), even when its physical span exceeds the range
of the short-addend estimate. The gain in the latter example
requires \(\eta<1/6\).

## Counting nearby pairs across all blocks

Let \(\mathcal A\) be any set of integers in \((X,2X]\) coprime to
\(M\), each of radical at most \(X^\rho\). Count all unordered pairs
in \(\mathcal A\) in the same residue class modulo \(M\) with
\(0<(v-u)/M\le Y\), over the whole height interval. Under the same
parameter assumptions the count \(E_M\) satisfies
\[
       E_M\ll X^{aB+(1-3a)_+\min(\rho,\tau)+\epsilon},
       \qquad (t)_+=\max(t,0).                                      \tag{PI2}
\]
Here distinct scalings of one derived triple can occur. If, in
addition, each \(u\in\mathcal A\) has
\(\prod_{p:\,0<v_p(u)<k}p\le X^\xi\) for a fixed integer \(k\)
with \(1/k\le a\), then
\[
                           E_M\ll X^{a B+2\xi+\epsilon}.              \tag{PI3}
\]
In particular there is no extra principal exponent for \(k\)-full
inputs under this condition. The constants allow the fixed
parameters and \(k\), and are uniform in the set and the step.

For the proof group \(s=\log_X g\) and the exponents
\(r_1=\log_X\operatorname{rad}(u/g)\) and
\(r_2=\log_X\operatorname{rad}(v/g)\) in bins of sufficiently small
fixed width. We write the cost at an exact scale, with arbitrarily
small errors from the endpoints of the bins. The reduced height
is of order \(X^{1-s}\), and its radical cap in base \(X\) is
\(r_1+r_2+\mu+\tau-s\). This numerator is nonnegative and is at
most \(B-s\), so the relative budget is at most \(B\).
Use the original count at each bin's height, with a fixed budget
rounded up slightly. Thus the exponent for primitive triples in
that bin is at most \(a(r_1+r_2+\mu+\tau-s)\). If the short-addend
condition is used, the comparison \(Mj/g\le(v/g)^{h_*}\) holds
before any binning. All available budgets and the number of bins
are fixed in the height limit.

For a fixed derived triple there are at most \(O(X^s)\) possible
scalings in its height bin, up to the same exponent errors.
Also the part of \(g\) whose support lies outside
\(\operatorname{rad}((u/g)(v/g))\) has radical at most
\(X^{\rho-\max(r_1,r_2)}\). This follows from both original radical
caps. On the complementary prime support there are subpower
choices for the full part of \(g\). The one-variable radical count
bounds the number of choices for the outside part by its radical cap
times \(X^\epsilon\), uniformly. Hence the scaling exponent is at
most \(\min(s,\rho-r)\), with \(r=\max(r_1,r_2)\).
[Uniform fixed-support and radical counts](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

Since replacing \(r_1+r_2\) by \(2r\) only enlarges the exponent,
we maximize
\[
                     a(2r+\mu+\tau-s)+\min(s,\rho-r)
\]
over \(0\le r\le\rho\) and \(0\le s\le\tau\).
For \(0\le s\le\rho\) the expression is increasing in \(r\) up
to \(r=\rho-s\) and decreasing thereafter, because \(2a<1\).
Its maximum at that \(s\) is \(aB+(1-3a)s\). For \(s\ge\rho\)
the maximum is \(a(\mu+\tau-s)+\rho\), decreasing in \(s\).
This proves (PI2). For (PI3), any prime with a positive valuation
less than \(k\) in \(g\) has that valuation in at least one
endpoint. Hence \(\operatorname{rad}(g)\le X^{s/k+2\xi}\).
Use this scaling cap in place of \(\min(s,\rho-r)\) and use
\(s/k-as\le0\). This gives (PI3) after resizing the bin errors
and the subpower exponents.

For example, at \(M=1\), \(\tau=\rho\), and with a \(k\)-full
dictionary for \(k\ge 4\), the global number of nearby pairs is
\(O(X^{(1-3\sigma)\rho+\epsilon})\). A dictionary of nearly full
radical-cap size \(X^{\rho-o(1)}\) can therefore be made
\(X^\rho\)-separated by removing a power-smaller set. This
conclusion remains true with a small positive low-valuation
radical cap \(\xi\) when \(2\xi<3\sigma\rho\), with strict room
for the count allowances. More generally partition \((X,2X]\)
into intervals of width \(MY\), and also into unit classes modulo
\(M\). For a set of size \(S\), the number \(J\) of occupied cells
obeys \(J\ge S^2/(S+2E_M)\) by Cauchy–Schwarz. This is a spreading
bound on any such set.

## Stopping a full divisor above the block width

The following bound does not use the \(abc\) counts. Fix an integer
\(k\ge2\) and \(P\ge2\), and in a reduced progression block of
parameter width \(Y\ge1\) count integers \(v\) for which
\[
                  \prod_{\substack{p\le P\\v_p(v)\ge k}}
                                 p^{v_p(v)}>Y .                      \tag{PI4}
\]
Then the count is
\[
                                O_k(Y^{1/k}P) .                     \tag{PI5}
\]
The physical step can be arbitrarily large. The power threshold
only applies to the selected part of \(v\).

To see this, build a divisor of the product in (PI4) starting at
one. At each eligible prime insert \(p^k\) as a whole block,
then insert the remaining factors \(p\) one by one. Stop at the
first product \(D>Y\). It is \(k\)-full throughout, and
\(Y<D\le YP^k\). The reduced progression class makes \((D,M)=1\),
so each possible \(D\) divides at most one value in the block:
the parameters of its multiples are in one class modulo \(D>Y\).
The number of \(k\)-full integers up to \(Z\) is \(O_k(Z^{1/k})\).
Indeed write each as \(w^k\prod_{j=k+1}^{2k-1}b_j^j\): each
prime with a nonzero residue exponent modulo \(k\) contributes
one copy of \(p^{k+r}\) for that residue \(r\), and the remaining
exponent goes in \(w^k\). The sum over the \(b_j\) after counting
\(w\) converges, bounded by the product of \(\zeta(j/k)\).
Using \(Z=YP^k\) proves (PI5). Stopping below \(Y\) and bounding
multiples by \(Y/D\) instead would give the greater overshoot
factor \(P^{k-1}\). The upper stopping rule keeps the loss linear
in \(P\).

For \(v\ge X^\nu\) with radical at most \(X^\rho\) and all repeated
primes at most \(P\), condition (PI4) follows from
\(\nu>(k-1)\rho+\log_X Y\). More generally use just a positive
numerical weight of small primes in a fixed positive scaled
valuation band. For example, along the iterated small-budget
limit \(L\to0\) let \(n=X^L\), fix \(h>0\), and suppose one band
of orders \(e\ge s/L\), for fixed \(s>0\), contributes numerical
weight at least \(w>0\) on each retained value, using primes
\(p\le n^{b_L}\) with \(b_L\to0\). For \(Y\le n^h\), take
\(k=\max(2,\lfloor s/(2L)\rfloor)\). Then (PI4) holds for all
small enough \(L\), and (PI5) gives at most
\(O_L(n^{h/k+b_L})\) values per block. A positive nonatomic mass
in a compact positive scaled-degree band supplies such a block
of primes after fixing anatomic boxes and removing large root
coordinates. Each prime in a retained coordinate is bounded by
its root. The nonatomic mass allows retaining a fixed positive
weight among coordinates whose own normalized caps tend to zero.
This is sufficient even if other primes of \(v\) are macroscopic.
All divisor counts take \(k\) fixed at their own budget before
the inner height limit.

## Use on a distinct-target record and comparison with short-interval work

Suppose an affordable pair record of normalized cost \(t\)
has remaining allowance \(r=1/3-t\) in base \(n=X^L\).
Assume its actual target dictionary has radical cap \(n^z\).
If it is covered by at most \(n^{s+\epsilon}\) unit progression
blocks with parameter width \(n^h\) and step radical at most
\(n^\mu\), (PI1) bounds its exponent by
\[
                         s+\frac{a}{2}(2z+\mu+h)
\]
up to the budget and exponent allowances. A strict margin
below \(r\) suffices for a saving after restoring the pair
records. If all targets instead have the positive small-prime
weight used in (PI4) and (PI5), the bound is \(s+o_L(1)\)
for fixed \(h\) and variable reduced steps. None of these
results bounds the number of progression blocks required by
arbitrary fibres. Even zero nearby pairs permits one target in
each of many different blocks.
[Affordable pair records](ABC_WIDE_COFACTOR_CONGRUENCES.md#affordable-whole-coordinate-pair-locks-at-every-balanced-endpoint)
· [Short affine cells with actual target divisors](ABC_SPARSE_AFFINE_CONTENT.md)

Chan's short-interval problem asks for \(O_{\ell,\epsilon}
(Y^{1/\ell+\epsilon})\) values with \(\operatorname{rad}(v)\le
v^{1/\ell}\) uniformly over \(1\le Y\le X\). His conditional
upper bound for all such intervals uses the \(abc\) conjecture
and excludes certain three-term progressions. We use established
small-budget counts and two-term gaps. Inequality (PI1) gives
an unconditional power saving in a restricted range of widths
and does not give Chan's conjectured exponent there. Chan's
unconditional estimate in the sequel uses a smoothness
condition on repeated primes but no lower height constraint of
the form (PI4). In its January 2026 revision, Theorem 1 gives
\(O_k(Y^{(3k+1)/(4k)}+Y^{1-\delta}/\log(Y+1))\) for
\(k\)-powered values whose repeated primes are at most
\(Y^{1-\delta}\). Proposition 1 of that paper instead selects a
smooth divisor without maintaining its fullness. The upper stopping argument gives a different
range when a large smooth full divisor is available.
[Chan, *Powered Numbers in Short Intervals*, Theorem 1 and Conjecture 1](https://math.colgate.edu/~integers/a6Proc23/a6Proc23.pdf)
· [Chan, *Powered numbers in short intervals II*, Theorem 1](https://arxiv.org/html/2406.10014v2#S1)

Zeng's count of distinct exceptional targets in Theorem 11.7(ii)
uses only \(\operatorname{rad}(c)\le\operatorname{rad}(abc)\)
and therefore has exponent \(\lambda\) when the full triple
budget is \(\lambda\). It supplies no replacement for the
balanced \(1/3\) coefficient. The density deficit for sums of
three cube-full integers in Beyer de Ryke's Theorem 1.1 is a
relative density result in a progression chosen after its
error allowance. It does not bound a power-size population
of occupied progression blocks with an additional radical cap.
[Zeng, *A lattice-theoretic approach to the abc conjecture*, Theorem 11.7(ii)](https://hexagonmath.org/pdf/2610.00630v1)
· [Beyer de Ryke, Theorem 1.1](https://arxiv.org/html/2609.35772v1)
