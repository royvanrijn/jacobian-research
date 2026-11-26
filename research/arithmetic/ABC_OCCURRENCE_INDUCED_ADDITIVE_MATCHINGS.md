# Induced additive matchings at balanced abc occurrence

The occurrence endpoint can be reduced to an induced
matching of exact primitive sums with all three original
radical caps. A short real filter reduces the size of a
modulus needed to preserve the matching, with the real
filter still imposed on the model. Two radical caps alone
permit induced matchings of the full normalized exponent,
even in the presence of fixed full divisors.

## Extracting an exact induced matching under all three caps

Take \(n=X^L\). The family consists of positive primitive
rows \(A+B=C\), \(X/2<C\le X\), with
\(\operatorname{rad}(W)\le n^{1/3+\tau}\) for all three
complete roles. Keep each stipulated height restriction
and exact least-repair cutoff. At fixed budget use the
usual allowances on the exponent. Assume for the first
extraction that the maximum degree in each complete role
is at most \(n^{u+\tau}\) in the family. Take all
asymptotics in \(X\) with \(L\) fixed.

**IM1.** If the original row exponent is \(d\), there
is a complete-label matching of exponent
\(t\ge d-u-O(\tau)-o(1)\). Here the \(o(1)\) is an
arbitrarily small counting allowance at sufficiently
small budgets. It has a subset that is *induced* in
its three complete dictionaries. In other words
\[
                A_i+B_j=C_k\quad\Longleftrightarrow\quad i=j=k
\]
on that subset. Its exponent can be taken at least
\[
                 t-\frac{(h-t)_+}{2}-\eta
                                                               \tag{IM1}
\]
for every fixed positive allowance \(\eta\), where
\(n^{h+o(1)}\) is any upper bound on the number of
true cross equalities from the three matching
dictionaries. In particular, \(h=1/3+O(\tau)+\eta\)
is available at sufficiently small fixed budgets.
The complete occurrence endpoint \(d=1/3,u=0\)
therefore admits an induced matching still of exponent
\(1/3\) in the outer limit.

First colour the primes with three colours and retain
each row whose prime support on every letter belongs
to that letter's colour. The probability of retention
is \(3^{-\omega(ABC)}=X^{-o_X(1)}\) uniformly, so some
colouring retains the exponent. All cross triples from
the surviving complete dictionaries are now pairwise
coprime. A greedy selection gives a matching of size
at least the retained row count divided by three times
the maximum degree. Every mixed true equality still
has all three radical caps and a target in the original
height shell. Its individual height restrictions and
any target-only repair cutoff are unchanged. When
bounding cross equalities by the unrestricted abc
bound, further eligibility restrictions may be
dropped. Use \(c\asymp X\) and a slightly enlarged
fixed budget \((1+O(\tau)+\eta)L\) to obtain the
stated choice of \(h\). An ambient bound on each
target's representation degree within these same
caps gives an alternative value of \(h\).
[Unrestricted small-budget bound](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)

In a complete matching, if two indices in a true
cross equality coincide, all three must coincide.
A true off-diagonal equality therefore requires
three distinct indices. Keep each index independently
with probability \(p\). The expected size of the
retained family minus the number of off-diagonal
equalities is at least \(pM-p^3 n^{h+o(1)}\), where
\(M\) is the matching size. Choose
\(p\) of order \(n^{-(h-t)_+/2-\eta}\) and delete
one index per remaining off-diagonal equality.
This proves (IM1). The conclusion is about these
actual complete dictionaries, with no modular or
additive inverse theorem assumed.
One can also make the three sets of exact radical
labels injective within their own roles. The fixed-support
count bounds the number of complete values with a
prescribed radical by \(X^{o_X(1)}\), so greedy thinning
for the three sets loses only an inner subpower.
Inducedness is preserved on the remaining complete
representatives and every target keeps its own
least-repaired height.
[Support occurrence under the exact repair cutoff](ABC_MINIMAL_REPAIR_SUPPORT_OCCURRENCE.md#exact-reduction-to-occurring-support-unions)

## Preserving an induced matching with a random modulus

Let an exact induced matching from (IM1) have size
\(M=n^{t+o(1)}\). Choose a prime in \([n^\beta,2n^\beta]\)
with fixed \(\beta>0\). One may instead choose a fixed
number of primes in fixed power-size intervals,
condition the primes to be distinct, and use their
product \(Q\), with \(\beta\) the sum of the exponents.
Every nonzero cross error \(A_i+B_j-C_k\) is of
polynomial height and is divisible by at most a
bounded number of primes in any one such interval.
The prime number theorem gives
\[
 \Pr_Q\{A_i+B_j\equiv C_k\bmod Q\}\le n^{-\beta+o_X(1)}
                           \qquad (i,j,k)\text{ off-diagonal}.
\]
The same bound holds for any nonzero label difference.

**IM2.** Some \(Q\) from the indicated scales has a
subset of indices of size
\[
              \min(M,n^{\beta/2})\,n^{-o(1)}
\]
which is an induced matching in the additive group
\(\mathbb Z/Q\mathbb Z\). For a fixed budget the small
loss can be replaced by any positive normalized
allowance. This size matches the full exponent of
the input when \(\beta\ge2t\).

For the proof again sample indices independently,
then discard one index for each label collision or
off-diagonal modular equality. Cross equalities whose
indices have exactly two distinct values account for
at most \(O(M^2)\) possible conflicts, and all other
conflicts for at most \(M^3\). The expected loss is
\[
            \ll n^{-\beta+o_X(1)}
                                 (p^2 M^2+p^3 M^3).
\]
Take \(pM\) a sufficiently small fixed-power factor
below \(\min(M,n^{\beta/2})\), then send the factor's
exponent to zero. Diagonal equalities persist.
This argument controls expectation over the prime
choice and the row choice. It is not a claim that every
prime in the interval has this property.
[Average prime tests for true support occurrence](ABC_MINIMAL_REPAIR_SUPPORT_OCCURRENCE.md#an-average-finite-field-test-equivalent-to-true-occurrence)

This scale does not yield a generic induced-matching
exclusion. At the critical relation
\(M\asymp Q^{1/2}\), even a group \(H\times H\) of order
\(Q\) has the induced matching
\(((x,0),(0,x),(x,x))_{x\in H}\). Moreover the additive
image of the ordinary integer labels modulo a list of
primes is cyclic, by the Chinese remainder theorem.
The power bounds for fixed-prime vector groups do not
apply to that image. A cyclic group of order \(Q\)
admits induced matchings of size \(Q^{1-o(1)}\):
choose a progression-free subset \(S\) of
\([1,Q/3]\) by Behrend's construction and use
\((s,s,2s)\) for \(s\in S\). This is an arithmetic group
comparison and imposes no radical restrictions.
[Kleinberg–Sawin–Speyer and the bounded-torsion setting](https://discreteanalysisjournal.com/article/3734-the-growth-rate-of-tri-colored-sum-free-sets)
· [Behrend's construction](https://www.pnas.org/doi/10.1073/pnas.32.12.331)

## A short real filter cuts the required modulus exponent in half

Assume all three roles are at least \(X^{h_0}\) for
some fixed \(h_0>0\). For an arbitrarily small fixed
count allowance \(\eta>0\), choose \(\theta>0\) so
that every interval of length \(O(X^\theta)\) contains
at most \(n^\eta X^\varepsilon\) eligible values from
any of the three role dictionaries. This uses (SI2)
for the original one-letter radical caps and the
lower bound on the numerical height. Choose the
inner \(\varepsilon/L\) smaller than the desired
normalized allowance. At a putative one-third
saturation point the fixed positive lower height
bound can also be obtained by the strict short-term
sector theorem.
[Low radical counts in short intervals (SI2)](ABC_SHORT_INTERVAL_RADICAL_FILTER.md#arbitrarily-small-normalized-counts-in-shorter-intervals)
· [Strict short-height sector](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-strict-saving-when-the-short-term-is-small)

Put \(H=\lfloor X^\theta\rfloor\) and
\(b(x)=\lfloor x/H\rfloor\). Choose one of the two
carry values \(\kappa=0,1\) to keep at least half
the diagonal rows with
\(b(A_i)+b(B_i)+\kappa=b(C_i)\).
A mixed triple satisfying this coarse identity has
\(\lvert A_i+B_j-C_k\rvert<2H\).

**IM3.** Some modulus \(Q\) of exponent \(\beta\) and
some choice of the carry value preserve an induced
matching on
\[
             \min\{M,\sqrt{M n^\beta}\}\,n^{-O(\eta)}
\]
indices for the joint test
\[
 b(A_i)+b(B_j)+\kappa=b(C_k),\qquad
                      A_i+B_j\equiv C_k\bmod Q .
                                                               \tag{IM3}
\]
The first relation is an equality of ordinary
integers. In particular \(\beta>t\) preserves the full
exponent with sufficiently small fixed allowance.
For an off-diagonal candidate with all indices
distinct, fixing \(i,j\) leaves at most
\(n^{\eta+o(1)}\) possible \(C_k\) by the short
interval count. Thus there are
\(O(M^2 n^{\eta+o(1)})\) possible false cross equalities
in the first test. If two indices coincide, one
complete value differs by \(O(H)\) from another
on its own role. This gives
\(O(M n^{\eta+o(1)})\) possible conflicts by the same
count, now using all three height and radical caps.
Apply the prime-test probability and the deletion
argument with \(p^3M^2\) and \(p^2M\) in place of
\(p^3M^3\) and \(p^2M^2\), paying the \(n^\eta\)
allowance. This proves (IM3).

Equivalently the matching is induced in the
additive group \(\mathbb Z\times\mathbb Z/Q\mathbb Z\)
on labels \((b(A_i),A_i\bmod Q)\),
\((b(B_i),B_i\bmod Q)\), and
\((b(C_i)-\kappa,C_i\bmod Q)\). The numerical
coordinate can be large. Removing it would reintroduce
the discarded cross incidences. This leaves a
recognizable quantitative obligation for the
occurrence endpoint. A small-model or incidence
bound for this mixed relation must use more than the
existence of an induced matching, and must retain the
target's actual radical condition.

## A primitive two-cap induced matching with diffuse degrees

**IM4.** There are actual primitive two-cap families
at arbitrarily small fixed \(L\) with both prescribed
radicals at most \(n^{1/3}\), all three numerical
heights one, and an induced complete-label matching
of exponent \(1/3\). Each affordable joint slot
record has the same residual entropy exponent as in
the simultaneous affine-record counterexample.
The target radical is deliberately unrestricted.
The two free limiting profiles can both be wholly
nonatomic with no escaped numerical weight.
Alternatively both free letters may carry prescribed
positive-weight common full factors with vanishing
normalized radical cost.

Let \(m\) tend to infinity through odd integers, put
\(L=m^{-2}\), \(n=X^L\), \(J=2m-1\), and
\(\xi=1/(3J)\). Set \(e_j=m(2m+j)\) for \(1\le j\le J\).
Choose \(2J\) pairwise disjoint equal-size pools of
primes in \([n^\xi/4,n^\xi/2]\) away from 2 and 3,
each of size \(P\asymp_m n^\xi/\log X\).
Let \(\kappa_A\in[2/m,2/3)\) and
\(\kappa_B\in[0,2/3)\) be obtained by dividing
nonnegative integers \(s_A,s_B\) by \(3m^2\).
Require \(s_B\) divisible by \(m\).
Let \(F\) be a power of 2 in \([X^{\kappa_A},2 X^{\kappa_A})\).
Let \(G\) be an \(m\)-th power of 3 in
\([X^{\kappa_B},3^m X^{\kappa_B})\), taking \(G=1\)
when \(\kappa_B=0\). Form
\[
 A=F\prod_j q_j^{\,e_j-s_A},\qquad
 B=G\prod_j p_j^{\,e_j-s_B},\qquad C=A+B ,
\]
with \(q_j,p_j\) from their respective pools.
Take margins below \(2/3\) while fixing the two
\(\kappa\)'s. The exponents are positive for large
\(m\), and their means are respectively
\(3m^2-s_A,3m^2-s_B\). The free numerical height
exponents are one by construction. The prime
interval bounds show \(A+B\le X\) for sufficiently
large \(X\) at each fixed \(m\) and the three letters
are comparable to \(X\) with constants depending
on \(m\). The radical caps hold since each pool
contributes at most \(n^\xi/2\), and the fixed
prime contributes just 2 or 3. Primitivity follows
from the disjoint prime supports.

Write \(B=z^m\). The positive integer \(z\) is
smaller than \(F\), because \(z\le X^{1/m}\) and
\(\kappa_A\ge2/m\). Since \(z\) is odd and \(m\)
is odd, raising to the \(m\)-th power is a bijection
of the unit group modulo \(F\), a finite 2-group.
Hence the complete values of \(B\) have distinct
residues modulo \(F\). Since \(F\mid A\), the entire
Cartesian sum map is injective, before sampling.
Consequently any complete matching between the two
dictionaries is induced on all three complete roles.

For simultaneous record bounds take the random
sample of Cartesian edges with probability
\(1/(P^J\log X)\), as in the existing record
construction. Its mean size is \(P^J/\log X\).
Expected collisions in either free dictionary are
\(o(P^J/\log X)\); targets never collide. By
injectivity and the fixed-support count, there are
at most \(n^{1/3+o_X(1)}\) candidate edges whose
target radical is at most \(n^{1/3}\). Their
expected sampled count is \(n^{o_X(1)}\).
Delete collisions and these targets. For every
one of the finitely many slot subsets \(I\) the
pre-sampling fibre size is \(P^{2J-|I|}\);
a union bound on binomial upper tails leaves all
sampled fibres of size at most
\(n^{(1/3-\xi|I|)_++o_X(1)}\) simultaneously.
This property survives deletion. Complementing
the fibre bound by the size bound for the slot
dictionary gives the claimed entropy for every
overlapping union. If both \(\kappa\)'s tend to zero,
the two limiting diffuse measures are
\(\tfrac13\operatorname{Unif}[2,4]\). For
fixed limiting positive \(\kappa_A,\kappa_B<2/3\),
the respective uniform bands are shifted by
\(3\kappa_A,3\kappa_B\) toward zero, and the
fixed numerical weights appear as escaped moments.
[Original simultaneous slot-record sampling](ABC_OVERLAPPING_FULL_RECORD_COUNTEREXAMPLES.md#simultaneous-entropy-saturation)

This example gives no three-cap exceptional set.
It shows why inducedness, numerical height in
the original equation, and even a pair of common
full divisors do not by themselves supply the
additional sparse-incidence bound in (IM3).
