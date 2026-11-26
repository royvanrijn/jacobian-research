# Cyclic projective transfer for pairs of primitive triples

The three cyclic determinants of a pair of primitive triples give
another additive triple. At a prime that cancels between the pair
in a relation of the projective fractions, the full prime power
appears squared in a specified determinant. The resulting primitive
triple retains those squared factors. Eisenstein factorization
recovers the input pair from the exact determinants with only
subpower multiplicity.

## Cyclic determinants and squared factors

Take two positive ordered primitive triples \(a+b=c\) and
\(u+v=w\) with \(X/2<c,w\le X\). Put
\[
 x=(a,b,-c),\qquad y=(u,v,-w),\qquad
 D_x=abc,\qquad D_y=uvw .
\]
Indices of signed letters are read modulo three. Define
\[
 \Delta_r=a y_{r+1}-b y_r\qquad (r=0,1,2).
\]
Then \(\Delta_0+\Delta_1+\Delta_2=0\) because the signed
letters of \(y\) add to zero. Explicitly,
\[
 \Delta_0=av-bu,\qquad
 \Delta_1=-aw-bv,\qquad
 \Delta_2=au+bw .
\]
In particular \(\Delta_1<0<\Delta_2\). The largest absolute
determinant is at most \(2X^2\) and at least \(cw/2>X^2/8\):
if \(a\ge c/2\), use \(\lvert\Delta_1\rvert\ge aw\), and
otherwise use \(\Delta_2\ge bw\).

Use \(F(z)=z-1/z+1/(1-z)\), and let \(z_x=a/c,\ z_y=u/w\).
The two nonidentity cyclic permutations of the signed letters
send \(z\) to \(1-1/z\) and \(1/(1-z)\). Direct substitution
shows that \(F\) is unchanged by both permutations. Exchanging
the two positive summands sends \(F(z)\) to \(1-F(z)\).

Let \(\mathcal P\) be any set of primes for which
\[
       v_p(D_x)=v_p(D_y)=e_p>0,\qquad
       F(z_x)-F(z_y)\in\mathbb Z_p .
\]
Since the triples are primitive, each such prime divides a
unique letter of each triple. Let \(r(p)\) be the unique cyclic
shift putting its pole in \(y\) in the same role as its pole
in \(x\). Put
\[
             Q_r=\prod_{\substack{p\in\mathcal P\\ r(p)=r}}p^{e_p},
             \qquad Q=Q_0Q_1Q_2.
\]
Then
\[
 p^{2e_p}\mid\Delta_{r(p)},\qquad
 p\nmid\Delta_{r'}\quad(r'\ne r(p)).
 \tag{CP1}
\]
To prove the first assertion, work with the shifted signed
triple \(y\), using its own value of \(c\) as minus its third
signed letter. The expressions for the principal parts of
\(F\) at the roles \(a,b,c\) are respectively \(-c/a\),
\(c/b\) and \(a/c\). This also holds for signed nonzero
triples, since the other two terms in \(F\) are integral at
the prime. In the chosen shift, both pole terms have the
same role \(W\) and order \(e_p\), so integrality of their
difference implies
\[
 p^{2e_p}\ \bigm|\ \ell_W(x)W(y)-\ell_W(y)W(x),
\]
where \((\ell_a,\ell_b,\ell_c)=(c,c,a)\). The linear map
\((a,b)\mapsto(W,\ell_W)\) is unimodular. Its determinant
on the two rows is therefore \(\pm\Delta_{r(p)}\).
For another shift, the two rows modulo \(p\) lie on
different projective pole lines. These are distinct over
every prime field, and each row has a unit coordinate.
Their determinant is nonzero modulo \(p\).
[Projective principal parts](ABC_THREE_CAP_RATIO_SUMSETS.md#full-congruences-from-the-last-occurrence-of-a-denominator-prime)

This applies to a pair of active rows in a fixed linear
relation in \(F\) when the two weights are equal up to sign.
Take the primes outside the coefficient and target pool
that occur on exactly these two rows. Their pole orders
must agree, since unequal orders leave a unique term of
least valuation. If the weights are opposite, their
contribution is an integral difference at each such prime.
If the weights have the same sign, first exchange the
positive summands of the second row. The local condition
then has the same difference form, since the exchange
subtracts the original fraction from one.
[Equal orders at exclusive pair primes](ABC_THREE_CAP_RATIO_SUMSETS.md#radical-surplus-and-separation-by-letter)

If the oriented rows are distinct, \(\Delta_0\) is
nonzero by primitivity, and all three determinants are
nonzero. Write \(g=\gcd(\Delta_0,\Delta_1,\Delta_2)>0\).
Rearrange the absolute values of \(\Delta_r/g\) so that
the largest is the sum of the other two. This is a
positive primitive triple, with each \(Q_r^2\) dividing its
corresponding term. Indeed (CP1) gives \(\gcd(g,Q)=1\).
An identical oriented pair has \(\Delta_0=0\). Before the
optional exchange this is the diagonal or the reflection
pair in the corresponding sign convention.

## Cofactors and recovery

The exact product identity is
\[
   \Delta_0\Delta_1\Delta_2
                  =-D_xD_y\bigl(F(z_x)-F(z_y)\bigr).
   \tag{CP2}
\]
For a derivation, write \(F(a/c)=N_x/D_x\), where
\(N_x=a^2b+(a-b)c^2\), and use the analogous expression
for \(y\). The form \(N_xD_y-D_xN_y\) vanishes when
\(x\) is proportional to a cyclic permutation of \(y\).
For each \(r\) the condition is the irreducible bilinear
determinant \(\Delta_r=0\). Indeed, changing \(y\) by the
cyclic permutation reduces to the determinant of two
independent pairs of variables. Hence \(\Delta_r\)
divides the form. The three determinants are pairwise
nonassociate and their product has the same bidegree
\((3,3)\). Evaluate at \((a,b,u,v)=(1,2,1,3)\), where
\((\Delta_0,\Delta_1,\Delta_2)=(1,-10,9)\) and the form
equals \(90\), to obtain (CP2).

If \(M_r=\Delta_r/Q_r^2\) for distinct oriented rows, then
\[
       g\mid M_r,\qquad
       \prod_r |M_r|
          =\frac{D_xD_y}{Q^2}\bigl|F(z_x)-F(z_y)\bigr|.
       \tag{CP3}
\]
For example, if \(Q=D_x=D_y\) and the two oriented
fractions lie in one interval of width less than one,
then the right side is less than one. The nonzero integral
cofactors exclude a distinct pair. In particular a bounded
interval for \(F\) can be split into finitely many smaller
intervals for this test. In the same-sign case the oriented
second fraction is \(1-F(z_y)\), so the interval condition
must hold after that exchange.

For recovery let \(\omega^2+\omega+1=0\), with
\(\omega\ne1\). The identity
\[
   (a-b\omega)(u-v\omega^2)=-\Delta_1+\Delta_0\omega
   \tag{CP4}
\]
follows by multiplication. Its norm is
\[
 \Delta_0^2+\Delta_0\Delta_1+\Delta_1^2
                  =(a^2+ab+b^2)(u^2+uv+v^2)\ll X^4.
\]
The ring \(\mathbb Z[\omega]\) is norm-Euclidean. Rounding
the two coefficients in the basis \(1,\omega\) leaves
error of norm \(s^2-st+t^2\le3/4\) for
\(\lvert s\rvert,\lvert t\rvert\le1/2\). Thus it is a UFD.
A nonzero element of norm \(N\) has at most
\(6\tau(N)^2\) divisors, since there are at most two prime
factors in the ring over each rational prime and each
exponent is bounded by the corresponding exponent in \(N\).
The exact \(\Delta_r\) therefore give \(O_\varepsilon(X^\varepsilon)\)
ordered input pairs. Given the primitive output and \(g\),
recover the exact \(\Delta_r\) up to the bounded number of
sign and role choices. The same multiplicity bound applies
to those data.

## A two-full-letter saving for pairs with a small cofactor product

Fix \(B>0\). There exist \(\theta,\tau,\sigma,L_0>0\),
depending on \(B\), with the following property. For every
fixed \(0<L<L_0\), count pairs in the shell above with
\(\operatorname{rad}(D_x),\operatorname{rad}(D_y)\le X^{BL}\).
Require distinct oriented rows and a selection of primes
\(\mathcal P\) as in (CP1) such that
\[
 \frac{D_xD_y}{Q^2}\max\{1,|F(z_x)-F(z_y)|\}\le X^{\theta L}.
 \tag{CP5}
\]
Suppose there is an integer \(p_*\) with \(p_*BL\ge3-\tau\)
and two of the three \(Q_r\) are \(p_*\)-full. Here an
integer is \(p_*\)-full if every nonzero prime valuation is
at least \(p_*\), with \(1\) allowed. Then the number of
these ordered pairs is
\[
                   O_{B,L,\varepsilon}
                          (X^{(B/3-\sigma)L+\varepsilon}).
 \tag{CP6}
\]
The choices of prime packets and their two full roles may
vary with the pair. The assertion counts pairs admitting
such a choice, with no multiplier for the number of choices.

Indeed (CP3) gives \(\prod_r|M_r|\le X^{\theta L}\) and
\(g\le X^{\theta L/3}\). The output height \(Y\) lies
between \(X^{2-\theta L/3}/8\) and \(2X^2\), and its
radical is at most \(X^{(B+\theta)L}\), since it divides
the radical of \(Q\prod_r|M_r|\). Set
\(\lambda'=(B+2\theta)L/2\). For sufficiently small
fixed \(L\) and large \(X\), this radical is below
\(Y^{\lambda'}\) by a fixed budget margin. On the two
specified roles, primes coming from \(Q_r\) have
valuation at least \(2p_*\), so their scaled degree
at budget \(\lambda'\) is at least \(p_*(B+2\theta)L
\ge 3-\tau\). The remaining primes on these roles
divide \(\prod_r|M_r|\) and have joint radical at most
\(X^{\theta L}\).

Use the two-full-letter saving with the small low-degree
mass allowance. In this application the normalized
radical mass of the low-degree primes is bounded by
\(\theta/(B+\theta)\) using the lower bound on \(Y\).
For clarity, the robust form follows from the same
compactness proof as the two-full-letter corollary.
If no positive allowance and exponent gap worked, take a
sequence with the allowance and the threshold deficit
tending to zero, and a count exponent tending to one third
at budgets tending to zero. The total anatomic root
budget tends to \(\lambda'\) in relative units. The
two prescribed limiting measures are then supported
on \([3,\infty)\), because the low-degree mass tends
to zero and the fixed-budget truncation can have
vanishing normalized cost. This contradicts the necessary
positive mass below three on at least two letters.
Thus choose \(\theta,\tau\) below the resulting fixed
allowances. Dyadic height bins absorb the use of \(Y\)
in the radical cap and in the normalized low mass. The
output count is \(O_{B,L,\varepsilon}
(X^{(1/3-\delta)(B+2\theta)L+\varepsilon})\) for some
absolute \(\delta>0\) from that corollary. Recovery by
(CP4) costs at most \(X^{\theta L/3+\varepsilon}\) over
all choices of \(g\). Choosing additionally
\(\theta<\delta B/2\) gives (CP6), taking
\(\sigma=\delta B/2\).
[Two-full-letter bound and its low-mass allowance](POWERFUL_TRIPLES_ABC_ENDPOINT_TRANSFER.md#corollary-a-small-budget-saving-with-two-full-letters)

Without (CP5), a pair can have a large \(Q\) in normalized
radical mass while the omitted prime powers still have
large numerical height. A compactly supported valuation
profile with unpaired raw mass \(o(1)\) only gives numerical
defect \(o(1)\) in base \(X\), whereas (CP5) needs a defect
of order \(L\). Likewise, fullness on two cyclic packets
is a restriction on where the high orders occur, distinct
from high fullness of two original letters. These remain
separate inputs for a counting application to the even
relations.
[Matching alternative for near-maximal relations](ABC_THREE_CAP_RATIO_SUMSETS.md#near-maximal-relations-pair-the-repeated-prime-supports)

## Raw imbalance and concentrated role transitions

There are two additional counting consequences of (CP5) that
require no fullness. Fix \(B>0\) and allow arbitrary small
fixed normalized exponent losses for all sufficiently small
fixed \(L\).

First, if \(\operatorname{rad}(Q_r)\ge X^{(B/3+\beta)L}\)
for some \(r\) with fixed \(\beta>0\), the pair count has
normalized exponent at most
\[
                        B/3-\beta/2+5\theta/6 .
 \tag{CP7}
\]
As above, the output is primitive and has normalized total
radical at most \(B+\theta\) in base \(X^L\). Its largest
radical cap is at least \(B/3+\beta\), because the output
retains \(Q_r^2\). In a fixed dyadic height shell for
the output, bin the three actual radical masses on a finite
mesh. Apply the two cheapest letter bound (GD3) with prime
prices equal to \(\log_X p\). The sum of those two masses is
at most \(2B/3+\theta-\beta\) up to the rounding loss.
The output shell heights have \(\log_X Y=2+O(L)+o_X(1)\),
so their relative budgets tend to zero with \(L\), uniformly.
The output exponent is at most half that sum. Recover the
pair by (CP4), paying at most \(\theta/3\) for \(g\).
[Two cheapest letter bound](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#group-dictionary-count)

Second, let \(\pi\) be a permutation of the three signed
roles for which the residues \(\pi(i)-i\) are all distinct
modulo three. Suppose the primes of \(D_x/Q\), together with
the marked primes whose role transitions \(i\mapsto j\)
have \(j\ne\pi(i)\), have product at most \(X^{\xi L}\).
Then the pair count is
\[
                            O_{B,L,\varepsilon}
                                (X^{(\xi+\theta)L+\varepsilon}).
 \tag{CP8}
\]
Indeed, fix this exceptional radical and the full powers
of its primes on \(x\), recording the transition of each
marked prime and indicating the unmarked ones. The support
has at most \(X^{\xi L}\) choices. For any one support,
allocations to the bounded set of roles and transitions,
and powers on it bounded by \(X\), have subpower count.
Let \(U_i\) be the complete supported part of \(x_i\) on
this fixed support. For \(r_i=\pi(i)-i\), write
\[
                     Q_{r_i}= E_{r_i}\,\frac{|x_i|}{U_i},
\]
where \(E_{r_i}\) is the fixed product of exceptional marked
powers assigned to that cyclic index. Every unrecorded
prime appears on exactly the indicated cyclic index, by
the hypothesis on \(\pi\). Fix the three signed nonzero
\(M_r\) at cost \(X^{\theta L+o_X(1)}\), using
\(\prod_r|M_r|\le X^{\theta L}\) and the divisor bound.
The equation for the cyclic determinants now becomes
\[
       \sum_{i=0}^2
            M_{r_i} E_{r_i}^{\,2} U_i^{-2} x_i^2=0 .
\]
In \(x_0=a,\ x_1=b,\ x_2=-(a+b)\) this is a nonzero
homogeneous binary quadratic. The coefficient of \(ab\)
is twice the nonzero coefficient of \(x_2^2\).
There are at most two rational projective solutions
and hence at most two positive primitive first rows.
The \(Q_r\) and \(M_r\) then determine the exact determinants,
so (CP4) gives only subpower further cost. There are three
choices of \(\pi\), since they are precisely the permutations
\(i\mapsto j_0-i\) modulo three, so one can also allow
\(\pi\) to vary across pairs.

## Replacing numerical cofactor size by cofactor radical

The two-full-letter saving (CP6) remains valid if (CP5) is
replaced by the weaker condition
\[
                    \operatorname{rad}\!\left(\prod_r |M_r|\right)
                                     \le X^{\theta L}.          \tag{CP9}
\]
The three cofactors may be numerically large. The constants
\(\theta,\tau,\sigma,L_0>0\) may be decreased in this form.

Write \(Y=\max_r|\Delta_r|/g\) for the positive primitive
output height. Since \(g\mid M_r\) and \(g\le 2X^2\), the
possible scales \(g\) number \(O_{\varepsilon}
(X^{\theta L+\varepsilon})\) by radical enumeration followed
by fixed-support enumeration. This pays for recovery in
(CP4) even if the output height is much shorter than \(X^2\).

For outputs with \(Y\ge X^{1/2}\), use a dyadic output shell
and choose a relative budget slightly above
\((B+2\theta)L\log X/\log Y\). It is at most
\(2(B+3\theta)L\), up to harmless dyadic losses. At each fixed
\(L\) a finite upward mesh can preserve both the radical
margin and the global bound \((B+3\theta)L\) on the
corresponding budget exponent in base \(X\). Since
\(Y\le2X^2\), the threshold \(2p_*\lambda'\ge3-2\tau\)
holds for large height. The output radical is at most
\(X^{(B+\theta)L}\). The two specified roles have scaled
orders at least \(2p_*\lambda'\) away from primes supported
on \(\prod_r |M_r|\). By (CP9) the low-order mass is at most
\(\theta/(B+\theta)\) in the output budget units.
Apply the same robust two-full-letter count as in (CP6).
The exponent for outputs is at most
\((1/3-\delta)(B+3\theta)L+\varepsilon\), for some fixed
\(\delta>0\). Adding the scale cost \(\theta L\) preserves
a strict gap when \(\theta\) is small enough in terms of
\(B,\delta\).

For outputs with \(Y<X^{1/2}\), write \(Q_{r_1},Q_{r_2}\)
for the two full packets. Each satisfies \(Q_{r_i}^2\le Y\).
Their combined radical therefore has normalized mass at most
\[
               \frac{\log_X Y}{p_*L}\le \frac{B}{2(3-\tau)} .
\]
The two complete output values on these roles have combined
radical at most \(X^{\{B/[2(3-\tau)]+\theta\}L}\).
Enumerating their supports and powers fixes both values,
hence the entire additive output up to a bounded choice of
role and signs. After scale recovery this costs at most
\(B/[2(3-\tau)]+2\theta\) in base \(n\), which stays below
\(B/3\). This proves the claimed extension.

An input formulation follows from (CP3). If
\[
  \operatorname{rad}(D_xD_y/Q^2)\,
  \operatorname{rad}\!\left(
    \left|\operatorname{num}\bigl(F(z_x)-F(z_y)\bigr)\right|
                         \right)
                 \le X^{\theta L},
\]
with numerator taken in reduced form on the nonzero oriented
difference, then (CP9) holds. Division by the denominator
can only delete prime factors of the integral cofactor
product. The first radical is small for a pair with almost
all its raw mass exclusively repeated at matching orders.
The additional condition is on the numerator support of
the oriented difference, together with the two required
full cyclic packets. Smallness of the raw defect alone
does not control that numerator.
[Uniform supported-part enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

