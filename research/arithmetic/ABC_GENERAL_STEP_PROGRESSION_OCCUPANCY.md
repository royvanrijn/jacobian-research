# Progression occupancy for arbitrary steps and low-degree roots

Fix \(\eta,\varepsilon>0\) and fixed nonnegative \(\rho,\tau\)
with \(0<2\rho+\tau=\lambda\) sufficiently small. Let
\(M\) be any positive integer and take a block of parameter
width at most \(X^\tau\) in a unit residue class modulo
\(M\). Its number of values \(X<u\le2X\) with
\(\operatorname{rad}(u)\le X^\rho\) is
\[
 O_{\eta,\rho,\tau,\varepsilon}\left(
    X^{\,\min\{\tau,\,3\rho/5,\,(2\rho+\tau)/5\}
                         +\eta\lambda+\varepsilon}\right).
 \tag{GO1}
\]
The constant is uniform in the step and the block.
For example the exponent is \(3\rho/5\) before allowances
at \(\tau=\rho\). The estimate permits arbitrary prime
supports and valuations on the values. When the step
has small radical or the block has a fixed physical
height gap, the stronger two-letter progression count
applies directly.
[Counting with a small step radical or a numerical gap](ABC_POWERED_PROGRESSION_OCCUPANCY.md#discounting-the-step-in-a-two-letter-bound)

The two ingredients in (GO1) are a conditional pair count
that charges only half the low-position cost on one side
and a direct count when the low part has large radical.
The pair count treats the progression step as a fixed
coefficient and the variable gap as one small integer.

## A grouped count with low positions on the two variable terms

Use the three-letter fixed-coefficient setting of (GD3),
with singleton root dictionaries and raw cap \(O(\lambda)\).
Designate one complete letter \(L=C R\) with arbitrary
fixed coefficient \(C\), and \(R\le X^{C_0\lambda}\).
Treat \(R\) as one degree-one variable, of charged cost
\(W_L\). The other two complete letters \(P,Q\) have
positive limiting variable numerical height, with
arbitrary fixed coefficients. Their original root
positions can have any degrees. Fix a scaled cutoff
\(a>0\), let \(B_P,B_Q\) be the charged costs of the free
positions of degrees below \(a/\lambda\), and put
\(B=B_P+B_Q\), \(B_L=0\). Write \(W=W_P+W_Q+W_L\)
for the summed charges and choose a most expensive
letter \(S\), so \(W_S=\max_T W_T\). Then
\[
          D\le\max\left\{B,\ \tfrac12(W-W_S+B_S)\right\}
                        +o(\lambda).
 \tag{GO2}
\]
The proof uses the same height-first convention, with
fixed polynomial height bounds and fixed \(a,C_0\).
If one can choose among tied maximizers, any one of
them gives the stated bound. As usual one can use
arbitrarily small fixed exponent margins in taking
limiting costs. If some positions are fixed by another
record first, their fixing cost is paid before applying
this conditional estimate.
[All-degree group count (GD3)](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#the-two-letter-effective-cost)
· [High roots with a determined low root (SHF5)](ABC_FIXED_SHIFT_STRICT_HALF_SAVING.md#releasing-the-high-roots-with-a-determined-low-root)

We prove (GO2) along a hypothetical violating sequence
of boxes and use units of \(\lambda\) for limiting
dictionary costs \(w_T,b_T\), total costs \(w,b\), and
count exponent \(d\). Thus \(d>b\). A letter is called
supported if \(w_T>w-2d\); this forces its variable
numerical weight into boundedly many original positions,
as in (GD3). This use of concentration does not bound
any fixed coefficient: it only uses large full factors
from the variable roots and complete unit residues,
absorbing the fixed coefficients into the target.
All variable weights are measured using the actual
dyadic raw root caps. The weight of \(L\)'s variable
part vanishes by its hypothesis.

Suppose \(w_C<d\) for some letter \(C\), and another
letter \(A\) is supported; let \(B'\) be the remaining
letter. Fix \(C\) completely. Under every positive-weight
high lock on \(B'\), the release of \(A\) by (GW) has
vanishing normalized cost. Thus \(B'\) concentrates as
well, since these locks must cost at least \(d-w_C\).
If \(C\) is \(P\) or \(Q\), it has positive variable
numerical height and supplies a full lock at cost
\(w_C\). Release both \(A,B'\) by (GW), a contradiction.
If \(C=L\), both tall variable lists are concentrated.
Freeze all their low positions at cost at most \(b\),
and all high positions outside their finite lists at
vanishing cost because the high degree lower bound
converts vanishing weight to vanishing raw mass.
Apply (SHF5) with \(R\) as the one determined low
variable. The resulting \(d\le b\) is again impossible.
We have proved that \(w_C<d\) precludes any other
supported letter whenever \(d>b\). Taking a cheapest
\(C\) gives the quarter inequality and hence \(d\le w/3\)
in this case, by the same algebra as in (GD3).

The violation of (GO2) makes \(S\) supported. There cannot
be three supported letters since \(d\le w/3\), and with
two supported letters the third has cost below
\(4d-w\le d\), which was just ruled out. Thus \(S\) is the
unique supported letter. If neither of the unsupported
letters is \(L\), they have positive variable numerical
weight, and the two (GW) switches in the proof of (GD3)
give a contradiction using only \(d>(w-w_S)/2\).

In the remaining case \(S\) is a tall letter and \(L\)
is one of the unsupported letters; name the other
tall letter \(Z\). If \(w_Z<d-b_S\), fix \(Z\) and all
low positions on \(S\), at combined cost below \(d\).
The other high roots on \(S\) outside its bounded list
have vanishing cost, and (SHF5) releases the list and
\(R\). Otherwise \(w_Z\ge d-b_S\). The strict inequality
\(2d>w_Z+w_L+b_S\) forces \(w_L<d\). Fix \(L\)
completely for the preliminary lock argument. A
positive-weight high lock on \(Z\) costs at least
\(d-w_L\), by releasing \(S\). As in (GD3), greedy
packing shows there is an original position on \(Z\)
of positive numerical weight and cost
\(\beta_Z\ge d-w_L-\xi\), for arbitrarily small fixed
\(\xi>0\). Its degree is \(\gg1/\lambda\) regardless
of whether it is above the specified cutoff \(a/\lambda\).
Freeze its complement and the low positions of \(S\).
With \(\xi\) sufficiently small this costs strictly less
than \(d\), because
\[
       w_Z-\beta_Z+b_S\le w_Z+w_L+b_S-d+\xi<d .
\]
The retained high roots on \(S,Z\) and the low root \(R\)
are again covered by (SHF5) at vanishing normalized cost.
This proves (GO2).

## Pair counts at a fixed low radical scale

Set \(k=\lceil a/\lambda\rceil\) with a fixed sufficiently
small \(a>0\), say \(a=1/4\). For an original endpoint
\(u\), call the radical of its part with positive
valuations below \(k\) its low radical. Restrict to endpoints
whose low radical lies in a fixed narrow exponent bin
\([b,b+\Delta]\) in base \(X\), where \(0\le b\le\rho\).
It suffices for occupancy to count one bin at a time.
We prove an upper bound on the exponent of the number
of pairs in this bin. All errors from the bin endpoints,
truncation and small-budget estimates can be absorbed
in arbitrarily small fixed multiples of \(\lambda\)
and arbitrary inner height-exponent allowances.

For two endpoints \(u<v\) in the block define
\(j=(v-u)/M\), \(g=\gcd(u,v)\), \(s=\log_X g\), and
\[
 P=u/g,\qquad Q=v/g,\qquad L=M(j/g).
\]
Here \(g\mid j\) since the block is in a unit class.
The triple is primitive. The first two terms are of size
\(X^{1-s+o_X(1)}\), with \(s\le\tau\). Apply the anatomic
factorization to \(P,Q\) with exceptional coefficient
loss arbitrarily small relative to \(\lambda\). Its
packing index can be chosen at least \(k\), so the raw
caps on the positions of degree below \(k\) are exactly
the low radical scales up to the dyadic errors. The
total raw cap is at most \(2\rho+\tau+o(\lambda)\)
after fixing the coefficient \(M\) on \(L\) and taking
\(R=j/g\) as a degree-one singleton. No radical of \(M\)
enters. The other fixed coefficients are small, so
the positive limiting numerical height requirements of
(GO2) hold. The low classification agrees with the
cutoff \(a/\lambda\), since \(e<k=\lceil a/\lambda\rceil\)
is equivalent to \(e<a/\lambda\) for integer degrees.
[BBLT, Lemma 2.2](https://arxiv.org/html/2410.12234v2#S2)

Write \(r_P,r_Q\) for the radicals of the primitive
endpoints in base \(X\), and \(B\) for their combined
low radical exponent. For a prescribed primitive triple,
let \(v_0\) denote the smallest of the following scaling
caps on the \(X\)-exponent, up to the bin allowances.
The part of \(g\) supported on the two primitive endpoints
has only subpower-many possibilities. The off-support
part is at most \(X^s\) in height and, upon attaching it
to each primitive endpoint, has radical cap
\(X^{\rho-r_P}\) and \(X^{\rho-r_Q}\), respectively.
Thus \(v_0\le s\) and \(v_0+r_T\le\rho\) for \(T=P,Q\).
These bounds follow by counting the off-support integer
at its smallest size or radical cap. A further cap is
\[
                  B+v_0\le 2b+O(\Delta+\tau/k).
 \tag{GO3}
\]
Indeed the low primes of \(P,Q\) and the primes of \(g\)
outside both primitive supports are disjoint sets, each
contained in the union of the original two low prime
supports and the primes of the high part of \(g\).
For a low prime of a quotient which was originally
high and disappears from the low support after multiplying
back by \(g\), membership in the gcd is necessary. If
its gcd valuation is below \(k\), then the minimum of
the two original valuations is below \(k\) and it belongs
to an original low support. Otherwise it belongs to
the high part of \(g\). The same dichotomy applies
directly to primes outside both primitive supports.
The high part of \(g\) is \(k\)-full and has radical at
most \(X^{\tau/k}\). The combined radical cap, and
the unrestricted one-variable radical count on the
off-support part of \(g\), give (GO3).

For a tall position \(S=P\) or \(Q\), let \(b_S\) be
its primitive low radical exponent, and let \(Z\) be
the other tall term. At least
\[
                       \ell=\max\{0,b_S-b-\Delta\}
\]
of \(b_S\) consists of primes that were originally high
on that side and hence occur in \(g\). They divide
exactly one primitive endpoint and therefore lie outside
the support of \(Z\). Together with the disjoint
off-support part of \(g\), they consume a part of the
original \(Z\)-side radical. Thus the cap in (GO3)
may additionally be restricted by
\[
                          v_0+r_Z+\ell\le\rho .
 \tag{GO4}
\]
No enumeration of the high-to-low primes is needed:
their mass lower bound holds for every admissible
scaling of the fixed primitive triple. All fixed-support
choices for the complementary part of \(g\) are uniform.

Write \(m_*=\min\{\rho,(\rho+\tau)/2\}\).
If \(L\) is a most expensive charged letter, use
\(B_L=0\) and \(r_P+v_0,r_Q+v_0\le\rho\)
and \(v_0\le s\) to bound the two-letter term of (GO2)
after scaling by \(m_*\). If \(S\) is tall, the two-letter
term in (GO2) involves \(Z,L\), whose combined cap is
at most both \(r_Z+\tau-s\) and \(r_Z+r_S\), up to the
bin allowances. The latter holds because \(W_L\le W_S\).
Using (GO4), \(b_S-\ell\le b+\Delta\), and
\(v_0+r_S\le\rho\), the tall-pair cap gives at most
\(\rho+b/2+O(\Delta)\) after scaling. The \(Z,L\) cap
gives at most \((\rho+\tau+b)/2+O(\Delta)\), since
\(v_0\le s\). The first term of (GO2) plus scaling
is at most \(2b+O(\Delta+\tau/k)\) by (GO3).
Consequently the pair exponent is bounded by
\[
                 \max\{2b,\ m_*+b/2\}
                                   +o(\lambda)+\varepsilon.
 \tag{GO5}
\]
This argument also permits counting across all residue
classes for the fixed step. Dividing the pair exponent
by two bounds the occupancy of a single block, apart
from the possible singleton. The term \(\tau/k\) is
\(O(\lambda^2)\).

## Counting high factors and finishing the occupancy bound

If \(M<X^{1/2}\), then \(\tau<1/4\) at sufficiently
small budget gives the fixed physical gap needed in
(SHF4). The two-letter progression bound then gives
occupancy at most \(X^{m_*/2+o(\lambda)+\varepsilon}\),
already strong enough for (GO1). Suppose instead
\(M\ge X^{1/2}\). For each allowed \(u\) write
\(u=H L_0\), where \(L_0\) is its complete part supported
on primes of valuations less than \(k\), and \(H\) the
remaining high part. In an original low radical bin
of lower exponent \(b\), the high part has radical
at most \(X^{\rho-b}\). There are at most
\(X^{\rho-b+\varepsilon}\) choices for \(H\) by the
uniform radical count. For each \(H\), invert it in
the block's unit class. Its low part \(L_0\) then
lies in a prescribed residue class modulo \(M\).
The upper bound \(L_0\le X^{k\rho}\) and the choice
\(a=1/4\), \(\rho\le\lambda/2\), give \(L_0<M\)
at sufficiently small budget. Thus \(H\) determines
the original value in the block, giving occupancy
exponent at most \(\rho-b\) in this bin.
[One-variable radical count](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

Combine this with (GO5) and the trivial parameter
bound. If \(\tau<\rho/2\), the trivial bound already
gives (GO1). Otherwise \(m_*\ge3\rho/4\). Put
\(b_*=(4\rho-2m_*)/5\), so
\[
       \rho-b_*=m_*/2+b_*/4,\qquad
       b_*\le\rho-b_* .
\]
For \(b\le b_*\), the square-root upper exponent
\(\max(b,m_*/2+b/4)\) is at most \(\rho-b_*\).
For \(b\ge b_*\), the high-part upper exponent
\(\rho-b\) gives the same conclusion. Therefore
the occupancy exponent is at most
\(\min\{\tau,(\rho+2m_*)/5\}\), which is exactly (GO1).
If \(\rho=0\) there are no positive integers in this
height range with radical one for large \(X\), so the
argument can assume \(\rho>0\). Narrow fixed bins
can be summed before taking \(X\) to infinity.
