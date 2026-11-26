# Near-maximal signed sumsets in a rational set or its reciprocal

## Standalone height theorem

**OS1.** Fix \(M\ge1\), a finite list of integer weight vectors
of lengths at most \(M\) all of whose weights are nonzero,
and \(\varepsilon>0\). Let
\(\mathcal R\) be a nonempty finite set of distinct positive
rationals, written in lowest terms with numerator and denominator
at most \(H\ge 2\), and put \(N=|\mathcal R|\). There is one choice
\(\mathcal T=\mathcal R\) or \(\mathcal R^{-1}\) for which, for every
weight vector \((e_1,\ldots,e_j)\) on the list,
\[
             |e_1\mathcal T+\cdots+e_j\mathcal T|
                 \gg_{M,\mathbf e,\varepsilon} N^j H^{-\varepsilon}.
\]
In fact, for each \(j\le M\) there is a selection of
\(\gg_M N^j\) ordered \(j\)-tuples in that same \(\mathcal T^j\)
on which each corresponding signed-sum map has fibres
\(O_{M,\mathbf e,\varepsilon}(H^\varepsilon)\). The selections
can be made simultaneously for the finite list. Rational weights
are allowed by clearing their fixed denominators.

### Prescribed support and private factors

Write \(P(x)\) for the prime support of the reduced numerator
times denominator of \(x\in\mathcal R\). Let
\(\pi(p)=\Pr_{x\in\mathcal R}(p\in P(x))\). For
\(k=2M-1\), choose \(\eta>0\) sufficiently small in terms of
\(k,\varepsilon\). Put in a prescribed prime pool \(S_0\) all
primes with \(\pi(p)\ge\eta\) and all primes of the specified
nonzero coefficients. Its product is at most a fixed power
of \(H\), since
\(\sum_p\pi(p)\log p\le2\log H\).
For \(k\) independent uniform samples, let \(S_1\) be the set
of primes outside \(S_0\) occurring in at least two positions.
With repetitions allowed, the union bound gives
\[
 {\mathbb E}\sum_{p\in S_1}\log p
   \le {k\choose2}\sum_{\pi(p)<\eta}\pi(p)^2\log p
   \le 2{k\choose2}\eta\log H .
\]
Thus at least half of the tuples have
\(\prod_{p\in S_1}p\le H^{\varepsilon/4}\), for example by
choosing \(\eta\le\varepsilon/[16(1+{k\choose2})]\).
Write \(S=S_0\cup S_1\) for their associated support.

In each position independently orient the fraction as \(U_i/W_i\)
so its private numerator \(u_i=U_i^{S^c}\) is no larger than its
private denominator \(d_i=W_i^{S^c}\). Here \(J^{S^c}\) retains the
exact powers of all primes outside \(S\), and \(J^S\) retains the
complementary exact powers. Choose a fixed convention in a tie.
At least \(M\) positions use the same orientation. Pigeonhole
their position set and common orientation, then project to those
\(M\) positions. Each projected tuple has at most \(N^{k-M}\)
preimages. This leaves \(\gg_M N^M\) distinct tuples in
\(\mathcal T^M\) with a witness \(S\) on which all those private
inequalities hold. The private prime supports in the retained
positions are disjoint because they already were disjoint before
projection. Keeping a prescribed \(j\)-position sublist gives
\(\gg_M N^j\) tuples with the same property. Use any one witness
per surviving tuple.

### Recovery from one exact weighted sum

Fix a reduced output \(A/B=\sum_{i=1}^j e_i U_i/W_i\) with
\(B>0\). Enumerate the auxiliary shared support \(S_1\), at
most \(H^{\varepsilon/4}\) possible squarefree products.
On the fixed support \(S\), enumerate all \(U_i^S,W_i^S\)
at \(H^{o_H(1)}\) cost. This standard bound is uniform even
for varying polynomial-height support pools. For any supported
integer at most \(H\), its support is chosen from a pool
of \(O(\log H/\log\log H)\) primes since the pool product is
polynomial in \(H\). Each exponent sum is \(O(\log H)\), so
the number of exponent vectors is \(H^{o_H(1)}\) by
stars and bars.

Every prime of \(d_i\) is a unit on all other denominators
and on \(e_i U_i\). Thus its exact power in \(d_i\) survives
in \(B\). In particular \(d_i\mid B\), and \(B\le H^j\).
Enumerating the \(d_i\) at divisor cost fixes all \(W_i\).
For \(D=\operatorname{lcm}(W_1,\ldots,W_j)\), multiply the
sum identity by \(D\). Modulo \(d_i\) all terms except the
\(i\)-th vanish and \(e_i(D/W_i)U_i^S\) is a unit.
This gives a single class for \(u_i\) modulo \(d_i\).
Since \(1\le u_i\le d_i\), it gives \(u_i\) uniquely, including
\(d_i=1\). Thus there are at most \(H^{\varepsilon/4+o_H(1)}\)
tuples over this exact output. This proves the stated bounds,
with constants for the finitely many small heights. It applies
whether the sum is zero or nonzero.

## Two rational maps with different poles

The same conclusion on signed-sum cardinalities holds for
one of \(f({\cal R})\) or \(g({\cal R})\), where \(f,g\)
are fixed rational Möbius maps with different poles.
Omit their poles and at most one additional exceptional
input. After applying \(f\), the second map has the form
\(g\circ f^{-1}(x)=u+v/(x-w)\) with \(v\ne0\).
Split \(f({\cal R})-w\) by sign and apply the theorem
to the larger part, making it positive by a fixed sign
if necessary. The height of each image is \(O(H)\), and
an affine change preserves every fixed weighted-sum
cardinality. Bounded-size sets satisfy the claimed
cardinality bound directly whenever the images under
consideration are nonempty.

For original primitive odds the three choices of
denominator role have different poles, represented by
\(\infty,0,1\) in the parameter \(A/C\). The two possible
numerator roles over a given denominator differ by an
affine change with fixed coefficients. Therefore on any
retained row subfamily, at least two of the three
denominator roles have near-maximal *whole-dictionary*
weighted-sum growth simultaneously for the stated finite
list. In particular two compressed whole-dictionary odds
projections on the same large row subfamily require
additional restrictions on the tuples.

The reciprocal choice matters even for unrestricted tuple
sets. For prime \(H\), the reduced fractions \(a/H\) with
\(H/3\le a\le H/2\) form a set of size \(\asymp H\) whose
signed sums in that fixed orientation have only \(O(H)\)
possible outputs for each fixed weight list and length.

## Genuine three-cap abc application with a power-small exceptional set

**OS2.** Put \(n=X^L\) and consider \(N\ge n^{1/3-\tau}\) distinct positive
primitive shell rows \(A+B=C\), \(X/2<C\le X\), with genuine
total radical cap \(n\) and individual caps
\(n^{1/3+\tau}\). Fix two distinct roles to make their
primitive odds, for example \(A/C\) and \(C/A\). Either odds
determines the other letter and the complete original row.
There is no numerical lower bound on the smaller of these roles.

For any fixed tuple length \(k\), a fixed prescribed pool of
polynomial-height product, and any \(\sigma>0\), discard tuples
where the product of shared original primes outside the pool
exceeds \(n^\sigma\). At sufficiently small fixed \(\tau,L\)
these form only \(O(N^k n^{-c})\) tuples for some \(c>0\).
Indeed, use the union over pairs of positions, with threshold
\(\sigma/{k\choose2}\) for a pair if \(k\ge2\). Fix the first
row of a pair. Its prime support is a prescribed zero-price
pool in EB5, and a second row with outside overlap of that
mass has exponent at most \(1/3-\sigma/(3{k\choose2})+\eta\).
Choose \(\eta,\tau\) and the fixed-height-count \(X^\epsilon\)
allowance smaller than the gap. The first row and other
positions are unrestricted. Any concern about distinctness can
also be addressed by deleting \(O_k(N^{k-1})\) repeated tuples.

On the surviving tuples let \(S\) consist of the prescribed
pool, all coefficient primes, and the actual shared original
support. The same private-size orientation and exact weighted
recovery just proved apply, but the auxiliary support costs
\(n^\sigma\), with no Markov loss and with only height-subpower
additional cost. Therefore in any selected good sector with
at least \(n^{k/3-\rho}\) tuples, its naturally oriented
exact signed-sum image has size at least
\(n^{k/3-\rho-\sigma-o_X(1)}\). The same statement uses a
fixed pattern whenever it agrees with the chosen orientations,
or one can pigeonhole the \(2^k\) patterns before taking the
image. The \(2M-1\) projection gives a single pure orientation
for all specified fixed lengths at most \(M\), retaining a
constant-proportion tuple sector for each length.

This argument retains all three original caps through EB5.
It applies to every balanced dictionary profile, including
occurrence. A compressed exact odds-shift record on a
power-large tuple sector must either fail the stated
orientation or pay the full tuple count up to the explicit
retention and overlap losses. An extraction argument
would have to supply a different arithmetic record
or extra structure on the actual tuple selection.

### Prescribed tuple sectors and pure orientations

For a prescribed sector of pairs or triples with small
shared support, one can use a pure orientation without
forming longer tuples. In each tuple take a majority
of the private-size orientation choices. At most one slot
then has private numerator larger than private denominator.
From one exact weighted sum, CRT recovers all favorable
slots as before. The sum recovers the last fraction
directly. Pigeonhole the majority orientation over the
sector. The total fibre cost is still only the actual
shared support price and a height-subpower factor.

More generally for \(j\) slots take
\(t=\max(1,\lfloor j/2\rfloor)\) exact weighted sums
whose weight matrix has full column rank on every
sublist of at most \(t\) columns, with all entries of
its first row nonzero. A fixed Vandermonde matrix with
distinct integer column parameters works. The first sum
gives the private denominators and reconstructs the
favorable slots. After subtracting their contributions,
the other fractions are uniquely recovered by the
linear system. This gives the same subpower completion
after the shared support price.

One pure orientation and one sum do not have this property
on arbitrary prescribed sectors of length four. For a
power of two \(H\), use primes \(p\) between \(H/3\) and
\(H/2\) and the set containing both \(p/H\) and \(H/p\).
Select tuples with two instances of the first form in
the first two slots and two of the second form in the
last two, with all four primes distinct. The sector has
\(H^{4-o(1)}\) tuples, with shared support just \(\{2\}\)
or the empty set. In either pure orientation a pair of
slots with common denominator \(H\) gives only \(O(H)\)
weighted partial sums, and the other two slots give at
most \(O(H^2)\). Either pure signed-sum image therefore
has size \(O(H^3)\). This separates arbitrary retained
tuple sectors from the whole-dictionary conclusion.

## Exact shifts across actual divisor fibres

**OS4.** Put \(n=X^L\) with fixed positive \(L\).
Suppose \(N\ge n^d\) genuine primitive rows of height at
most \(X\) in one signed orientation are assigned actual
divisors \(D\) of the third role, with total exact divisor
dictionary at most \(n^c\) and \(d>c\) by a fixed margin. Write \(N_D\)
for the partition sizes and use the odds of the other
two original roles. For any fixed weight vector of
length \(j\ge2\) with nonzero coefficients of sum zero,
the global exact weighted-sum dictionary, allowing either
odds orientation for each whole fibre, has exponent at
least
\[
                          jd-(j-1)c
\]
up to arbitrarily small fixed allowances.

Indeed the odds and its inverse are units congruent to
the same fixed sign modulo the actual \(D\). For a
nonzero exact weighted sum \(A/B\) in lowest terms,
\(B\) is a unit at \(D\) and \(D\mid A\). The same output
therefore belongs to at most divisor-many actual \(D\)
fibres, at polynomial height. Apply OS1 to the rows
of each \(D\) fibre, with \(H=X\) and a height allowance
arbitrarily small compared with the fixed \(L\).
Discard fibres of size less than \(n^{(d-c)/2}\);
they contain at most \(n^{(d+c)/2}\) rows altogether.
The output zero accounts for one dictionary entry and
does not affect the lower bound in each remaining large
fibre. Sum the surviving image sizes and use
\(\sum_D N_D^j\ge(\sum_D N_D)^j/(\#D)^{j-1}\).
The proof also allows an orientation tag at a bounded
factor cost. It keeps all original caps and the actual
divisors throughout. In APC11 one may take
\(c=b_0+C(r)\) up to allowances. Thus an improving
APC7–APC8 certificate would force a strict increase
\(d+(j-1)(d-c)\) in these exact shift dictionary exponents.
[Actual divisor catalogue and fibre bounds (APC11)](ABC_ADAPTIVE_POWER_COEFFICIENT_PRICES.md#actual-divisor-multiplicities-in-a-joint-pool-extraction)

## Fresh ordinary radicals and power defects of summed odds

**OS3.** Keep the genuine three-cap hypotheses and fixed nonzero weights
\(e_1,\ldots,e_j\), with \(j\ge2\). Let \(Z=A_*/B_*\) be
their naturally oriented reduced exact signed sum as above,
where \(B_*>0\). Fix a finite list of squarefree integral
homogeneous binary forms \(F\), each of degree at least two
and satisfying \(F(1,0)\ne0\). Fix any prescribed
polynomial-height prime pool. Define fresh support by omitting
this pool and all original primes of the \(j\) rows. For
every \(\delta>0\), outside \(O(N^j n^{-c})\) tuples at
sufficiently small fixed \(\tau,L\), every \(F(A_*,B_*)\)
is nonzero and has fresh radical greater than
\(n^{j/3-\delta}\). For forms of degree at least five the same
bound holds for their fresh \(m\)-defect radicals, for each
fixed \(m\ge2\). All choices of weights, forms, orders and
positive margins precede the small-budget cutoff and the
height limit. As in GP2, finite order lists at each fixed
budget are permitted with constants depending on the list,
if only the inner height limit is taken afterward.

To remove old support, first remove the power-small shared
support exception with allowance \(\sigma\). Enlarge the
prescribed pool by primes dividing coefficients and endpoint
values. Use raw numerator \(T=\sum_i e_iU_i\prod_{\ell\ne i}W_\ell\)
and denominator \(D=\prod_i W_i\). For a private prime on
\(W_i\), \(T\) is a unit and \(D=0\bmod p\), so \(F(T,D)\)
is a unit. For a private prime on \(U_i\), \(U_i/W_i=0\bmod p\).
For one on the third role, this odds is \(1\) or \(-1\bmod p\),
according to the additive orientation. In the latter two cases
\(F(T,D)\) at \(p\) is a unit multiple of the same form
evaluated at the raw anchored sum from the \(j-1\) other
rows, shifted by \(0\) or \(\pm e_i\). With the others fixed
this is a value of polynomial height. If it is nonzero, its
support is a prescribed pool in EB5 for the varying \(i\)-th
row. That estimate discards a power-small fraction on which
this overlap exceeds an arbitrarily small fixed positive
normalized threshold. If the anchor value is zero, fix all
but one of its rows. A nonzero univariate polynomial of
degree at most \(\deg F\) evaluated at that remaining
original odds vanishes. Its possible exact odds determine
its primitive row, so all such zero anchors account for
\(O(N^{j-1})\) full tuples. These arguments are applied over
the finitely many fixed orientation patterns and then
restricted to the naturally chosen pattern, so there is
no conditioning on how orientations were chosen.
Reduction from \((T,D)\) to \((A_*,B_*)\) cannot create a
new prime dividing \(F\). Hence its old support outside the
prescribed pool has arbitrarily small normalized mass off a
power-small exception. The same one-variable argument
discards the tuples with \(F(A_*,B_*)=0\).

Now suppose the fresh radical of one form costs at most
\(n^{j/3-\delta}\). Combine it with the small old support
and fix exact powers on the resulting support together with
the prescribed pool at height-subpower cost. For a
nondegenerate quadratic, completing the square and the
uniform quadratic norm estimate give height-subpower
primitive completions over the exact nonzero value.
For degree at least three, use Evertse's Thue–Mahler bound
as in GP1. For a fresh \(m\)-defect upper cap use a prime
divisor \(p\) of \(m\), enumerate the signed \(p\)-power-free
coefficient of the form, and use BT1 on each such
coefficient as in GP2. This has only height-subpower
primitive completions at fixed \(p,F\). Each completion
gives an exact rational output and so reconstructs the
original tuple at \(n^{\sigma+o_X(1)}\) cost. Choose
\(\sigma\), the overlap and cap tolerances and \(\tau\)
sufficiently small against \(\delta\). Comparing with
\(N^j\ge n^{j/3-j\tau}\) gives the claimed power-small
exception. Pairwise geometrically coprime forms have
disjoint fresh supports after deleting their resultant
primes and enjoy the corresponding independent defect
classes, exactly as in GP2.

## Comparison with rational energy bounds

Bloom and Maynard's Theorem 2 bounds the full even additive
energy using a denominator bound and the maximum count
at a single denominator. The height theorem above has
a tuple-selection conclusion after an optional reciprocal
change, with no per-denominator cap or assumption that the
available denominators are nearly filled. Bloom and
Kuperberg's Theorems 4–5 count odd rational equations with
interval or numerator restrictions. Their relative-gcd
decomposition includes the private denominator isolation
used here. The additional step in OS1 is to choose the
orientation using the *exact private factors* and then
complete the favorable numerators by CRT.
[Bloom–Maynard, Theorem 2](https://arxiv.org/pdf/2011.13266)
· [Bloom–Kuperberg, Theorems 4–5 and Appendix A](https://arxiv.org/html/2312.09021v2)
· [Earlier primitive comparison bounds (GP1–GP2)](ABC_GENERIC_PAIR_FRESH_BINARY_FORMS.md)
