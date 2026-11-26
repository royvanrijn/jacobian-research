# Three-role divisor records up to the projective height threshold

Primitive integer triples can have saturated dictionaries for every union of selected
exact valuation slots on all three roles. The prescribed prime powers have joint
numerical weight arbitrarily close to two in the projective height normalization.
The construction retains the full small-radical cap on one letter and gives explicit
lower bounds for the unrecorded radical on both other letters. This isolates the
scope of an argument that charges only its selected divisor records.

## Statement

Fix \(0<\theta<1/4\). For each sufficiently large integer \(m\), fix
\(L=1/m\), \(d=1/3\), \(J=2m-1\), and \(\beta=d/J\). Put \(n=X^L\).
For all sufficiently large \(X\) at this fixed budget, there is a family \(\mathcal H\)
of positive primitive integer triples \(a+b=c\) with \(X/2<c\le X\).
Each complete letter lies in \([X/8,X]\), and each of the three complete-value
projections is injective on \(\mathcal H\). Each row has specified divisors
\(F\mid a,\ G\mid b,\ M\mid c\), made of \(J\) exact valuation layers per
letter, with the following properties. Each selected layer consists of one
prime, and no unmarked prime has a valuation in that layer on its letter.

- The marked prime in each slot belongs to its own pool, disjoint from every
  other slot pool, including the pools on the other letters. Every pool has
  size \(P=n^{\beta+o_X(1)}\). Each marked prime has normalized radical mass
  \(\beta+o_X(1)\). Thus every marked divisor has radical at most \(n^d\)
  and normalized radical mass tending to \(d\).
- With \(w_m=\lfloor 3m(1-\theta)/2\rfloor/(3m)\), the numerical weights of
  \(F,G,M\), in base \(X\), are respectively \(1,w_m,w_m\) up to \(o_X(1)\).
  In particular their joint weight approaches \(2-\theta\) as \(m\to\infty\)
  after the height limit.
- In fact \(\operatorname{rad}(a)\le n^d\). On both other letters the
  unmarked cofactors have radical greater than \(X^{\theta/2}\), supported
  on primes outside their own marked part and with all valuations bounded
  in terms of \(\theta\) alone.
- The size of \(\mathcal H\) is \(n^{d+o_X(1)}\). For each subset \(I\) of the
  \(3J\) marked slots, put \(t_I=\beta |I|\). The number of rows with any
  fixed assignment on \(I\) is at most \(n^{(d-t_I)_++o_X(1)}\), uniformly
  for all \(I\) and assignments. When \(t_I<d\), at least one assignment
  has \(n^{d-t_I+o_X(1)}\) rows.

The marked valuation measures, using scaled orders \(Lv_p\), converge first
in height and then in \(m\) to
\[
 \frac13\operatorname{Unif}[2,4],
 \qquad
 \frac13\operatorname{Unif}\!
          \left[\frac32(1-\theta)-1,\frac32(1-\theta)+1\right],
 \qquad
 \frac13\operatorname{Unif}\!
          \left[\frac32(1-\theta)-1,\frac32(1-\theta)+1\right].
\]
Here \(\operatorname{Unif}\) denotes the probability measure on the indicated
interval. All three marked measures are diffuse in the limit. These also
give the complete-letter limits tested against any continuous function
with compact support in positive scaled orders, since the unmarked orders
on \(b,c\) are \(O_\theta(L)\) after scaling and the unmarked prime two
on \(a\) has vanishing normalized mass at each fixed budget.

## Primitive realization of every assignment

For \(1\le j\le J\) the orders of the marked powers on the three letters are
\[
 e_{A,j}=2m+j,\qquad
 e_{B,j}=e_{C,j}=\left\lfloor\frac{3m(1-\theta)}2\right\rfloor-m+j.
\]
They are positive for sufficiently large \(m\). The prime number theorem
provides \(3J\) pairwise disjoint pools of the same size
\(P\asymp_m n^\beta/\log X\) in \([n^\beta/4,n^\beta/2]\).
Take \(X\) sufficiently large that all pool primes exceed two.
Choose one prime from each pool and form the three products \(F,G,M\)
of the powers in their designated slots. The orders and the prime intervals
give
\[
              F\asymp_m X,\qquad G,M\asymp_{m,\theta} X^{w_m}.
\]
For example \(\beta L\sum_j e_{A,j}=1\), and the same sum for either
of the other two letters is \(w_m\). The logarithmic errors are
\(O_{m,\theta}(1/\log X)\). Their prime supports are disjoint.

The upper endpoint of the prime interval gives \(F<X/8\) for
large \(m\). Choose the least \(k\ge0\) with \(a=2^kF\ge X/8\),
so \(a<X/4\) and \(k=O_m(1)\). The same bound
\(F\le 2^{-\sum_j e_{A,j}}X\) gives \(k\ge \sum_j e_{A,j}-3>\max_j
e_{A,j}\). Thus none of the chosen \(a\)-orders contains the unmarked
prime two. Also \(\operatorname{rad}(a)\le 2(n^\beta/2)^J\le n^d\).

Seek \(b=Gy\) in \([X/2,5X/8]\). The condition \(M\mid a+Gy\)
prescribes one residue class \(y=y_0+M v\), since \(G\) is a unit
modulo \(M\). There are
\[
       K\asymp_{m,\theta} X/(GM)\asymp_{m,\theta}X^{1-2w_m}
\]
consecutive indices \(v\) in that interval, where \(1-2w_m\ge\theta\).
Take \(y\) odd and coprime to \(FG\). At each odd prime of \(FG\)
this rules out one class of \(v\), as \(M\) is a unit there. For
each prime of \(M\) also require \(p\nmid(a+Gy)/M\).
On the progression this last quotient is
\((a+Gy_0)/M+Gv\), with slope a unit at those primes.
Thus each of the \(3J\) marked primes forbids at most one class
of \(v\), and parity forbids at most \(\lceil K/2\rceil\) indices.
The number of excluded indices for the marked primes is
\(O_m(K/n^\beta+1)=o_{m,\theta}(K)\) as \(X\to\infty\).
The remaining cases are primitive because \(a\) is supported on
\(2F\), \(G\) is coprime to \(2F\), and \(y\) avoids \(2F\).
They have exact marked orders on \(b\) and \(c=a+b\) by the two
other groups of exclusions.

Fix an integer \(E\ge3\) with \(1/E<\theta\). Take \(m\) large enough
that \(E\le\min_j e_{B,j}=\min_j e_{C,j}\). On the retained progression,
both cofactors are units at the primes of \(2FGM\). In addition to
the preceding exclusions, the condition \(M\mid a+Gy\) makes \(y\)
a unit at each prime of \(M\), and \(\gcd(a,G)=1\) makes
\((a+Gy)/M\) a unit at each prime of \(G\). Primitivity gives the
unit condition on the quotient at the primes of \(F\), and the
chosen parity makes the quotient odd.

For \(p\nmid 2FGM\), rule out \(p^E\) dividing either cofactor.
Each cofactor is affine in \(v\) with slope a unit at \(p\).
The two exclusions cost at most \(2(K/p^E+1)\) indices per prime
\(p\le X^{1/E}\), since both cofactors are at most \(X\). In total
this is at most
\[
        2K\sum_{j\ge3}j^{-E}+O(X^{1/E})<K/4+o(K).
\]
The strict inequality in the main term follows from \(E\ge3\).
This leaves a positive proportion of the original progression.
All unmarked valuations on \(b,c\) are now less than \(E\), so
the selected full valuation layers consist of the claimed primes.

Both \(y\) and \((a+Gy)/M\) lie in \([1,X]\) and are injective
as functions of \(v\). The number of integers in \([1,X]\) of
radical at most \(X^{\theta/2}\) is \(O_{\theta,\epsilon}
(X^{\theta/2+\epsilon})\) for each \(\epsilon>0\), by enumerating
the radical followed by its supported powers. This is \(o(K)\)
if \(\epsilon<\theta/2\). Choose one of the remaining indices
for each assignment. The cofactors' radicals have the required
lower bounds and are disjoint from their own marked supports.
The target has \(X/2<c<7X/8\).
[Uniform supported-power enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

## Thinning all record fibres at once

In the assignment space of size \(P^{3J}\), a fixed complete letter
has at most \(O_{m,\theta}(P^{2J})\) assigned rows. Indeed each eligible
prime factor of that letter is at least \(n^\beta/4\), so an integer
at most \(X\) has only \(O_m(1)\) such prime factors. The marked
prime tuple for that role has \(O_m(1)\) possibilities. The two
other roles give at most \(P^{2J}\) possibilities. This bound
also handles different assignments that yield an identical complete
triple.

Select each assignment independently with probability
\(\pi=P^{-2J}/\log X\). The expected sample size is
\(S=P^J/\log X\). There are \(O_{m,\theta}(P^{5J})\) unordered
pairs of assignments sharing any complete letter, so the expected
number of such selected pairs is \(O_{m,\theta}(S/\log X)\).
The binomial lower-tail bound and Markov's inequality give, with
probability tending to one, a sample of size at least \(S/2\) with
\(o(S)\) colliding pairs.

For a fixed assignment on a subset \(I\), the fibre sample is binomial
with mean \(\mu_I=P^{J-|I|}/\log X\). Its exponential moment at one
is at most \(\exp((\mathrm e-1)\mu_I)\). By Markov's inequality the
probability of exceeding \(4\mu_I+C_m\log X\) is at most \(X^{-C_m}\).
There are at most \((1+P)^{3J}\) records over all subsets \(I\).
Fixing \(C_m\) large enough makes the simultaneous bound hold
with probability tending to one. This is exactly
\(n^{(d-t_I)_++o_X(1)}\) as an upper bound on the fibre sizes.
All of these upper bounds persist after deleting one endpoint
per collision. The remaining size is \(n^{d+o_X(1)}\), since
the deletion costs \(o(S)\) and the injective projection on
\(a\) has at most \(P^J\) values.

The pools are disjoint, so the product of any selected marked prime
powers and their tuple of roots determine one another. The same
fibre bounds therefore apply to unions of overlapping slot records,
using the union's cost and counting each occupied slot only once.
If \(t_I<d\), pigeonhole among the \(P^{|I|}\) root assignments
to obtain the stated lower bound on one fibre. Equivalently, for
a uniform row of \(\mathcal H\) the entropy of the record on \(I\)
is \(\min(d,t_I)\log n+o_X(\log n)\), and its conditional row entropy
is \((d-t_I)_+\log n+o_X(\log n)\). This follows from the upper
bounds on the number of records and the largest fibre, whose
entropy bounds add to the entropy of a uniform row.
[Two-letter simultaneous record construction](ABC_OVERLAPPING_FULL_RECORD_COUNTEREXAMPLES.md#simultaneous-entropy-saturation)

## The role of the unrecorded radicals

For two distinct primitive rows with a common record \(F,G,M\),
the three moduli divide the nonzero integer \(a b'-a'b\), by
using \(a+b=c\) for the modulus on the third letter. This explains
the numerical threshold two for direct projective uniqueness.
The construction attains any fixed joint numerical weight
\(2-\theta\) in the height-first limit, but its three full packet
records together cost one in base \(n\), above the affordable
exponent. Every affordable subset still exhibits the full residual
entropy \(d-t_I\). In particular this tests a simultaneous-lock
inference based on the assigned exact valuations, the complete
heights, primitivity, and the costs of all unions of the records.

There are affordable records involving all three letters whose numerical
weight exceeds one by a fixed margin. Start with all marked slots on \(a\).
In the limiting measures remove mass \(2q+\epsilon\) from its lowest orders
and add mass \(q\) from the highest orders of each of the other two letters.
The combined cost is \(d-\epsilon\). The removed weight is
\(2(2q+\epsilon)+O((2q+\epsilon)^2)\), and the added weight is
\(2q(5/2-3\theta/2)+O(q^2)\). Since \(5/2-3\theta/2>2\), their
difference is positive for a sufficiently small fixed \(q>0\) and then
sufficiently small fixed \(\epsilon>0\). Choose slot counts nearest
these masses at sufficiently large fixed \(m\). Both the affordability
and the weight margin persist, by the scaled-order interval calculation.
The selected three-role record still has a fibre with exponent equal to
the residual budget \(d-t_I\).

The two unmarked radicals are essential to the scope. On each of
\(b,c\) they add at least \(\theta/(2L)\) in base \(n\) to the
marked normalized mass. An argument for the exceptional-set count
must keep a constraint on these remainders, or derive one from
some other property of the actual family. The construction therefore carries three diffuse marked profiles and
one complete radical cap. The complete radicals of \(b,c\) exceed
the balanced caps.
[Three-factor projective locks](ABC_FIXED_FIBRE_VALUATION_WIDTH.md#discrete-rounding-on-the-full-three-letter-family)
