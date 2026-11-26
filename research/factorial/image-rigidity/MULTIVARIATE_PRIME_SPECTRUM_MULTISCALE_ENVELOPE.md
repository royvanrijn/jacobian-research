# Prime-spectrum separation for multiscale factorial envelopes

Fix an ordinary two-monomial core and finitely many appended row types
with smaller exponent \(c_iT^{e_i}+a_i\) for positive integers \(c_i,e_i\)
and \(a_i\ge0\), positive gap \(d_i\), and specified orientation. Allow
independent nonnegative multiplicities with at least one row present.
The following test permits equal \(2\)-adic valuations of the scale
exponents and permits common polynomial factors among ramps belonging
to different scales.

Use the absolute cutoff \(J\) and the residual adjacent orders
\(n,N=n+1\) from the
[unpaired-envelope theorem](MULTIVARIATE_UNPAIRED_ENVELOPE_MULTISCALE_MARKS.md).
At one order \(j\) of such a pair, choose a prime \(\pi\mid j\).
Choose one value of \(v_2(e_i)\) occurring among present types.
Among the types in this class, maximize
\[
                    \tau=\frac{v_\pi(j)+v_\pi(c_i)}{e_i}.
\]
Choose the largest scale exponent \(e\) of a maximizing type.
In its degree block let \(\mathcal T\) be the maximizing types,
and put \(c_*=\gcd_{i\in\mathcal T}c_i\).
Let \(\mathcal L\) be the set of odd primes dividing one of
\(e/\gcd(e,f)\) for the other scale exponents \(f\) of maximizing
types in the chosen class. Set
\[
         L=2 \quad (\mathcal L=\varnothing),\qquad
         L=4(|\mathcal L|+2)\quad(\mathcal L\ne\varnothing).
\]
Say this choice passes the envelope test if every interior grid index
has an attaining row of \(d_i/c_i\ge L/c_*\) in
\[
                    f(k)=\max_{i\in\mathcal T}
                     \frac{j a_i+d_i k_i}{c_i},\qquad
 k_i=k\text{ for increasing rows},\quad k_i=j-k\text{ otherwise}.
\]

**Theorem.** At each residual pair, one passing choice with \(e>1\)
suffices for coprimality above a computable height. If both orders
have passing choices with \(e=1\), they also suffice. The thresholds
are uniform in the multiplicities. The same conclusions hold when
the threshold \(L/c_*\) is replaced by the exact forward-distance
test at the selected mark below.

In particular, write \(\omega_{\rm odd}(e)\) for the number of odd
prime divisors of \(e\), and set \(L_e=2\) when this number is zero,
\(L_e=4(\omega_{\rm odd}(e)+2)\) otherwise. If every type has
\(d_i\ge c_i L_{e_i}\), the whole template is regular in every
degree at all sufficiently large heights, without any restriction
on its slopes or scale exponents. This stronger numerical assumption
ensures the grid test whichever type is present. When all present
slopes are \(\pi\)-units, the smallest exponent within the chosen
parity class is a unique maximizing block; then \(L=2\), regardless
of how many different exponents have the same \(2\)-adic
valuation.

## An admissible mark in a common mesh

Let \(F=\min_k f(k)\) on \(0\le k\le j\), and let \(Q=j c_*\).
Choose an integer \(u\) among the first \(L\) integers strictly
above \(c_*F\), with \(\pi\nmid u\) and such that \(u/Q\) is
not a rational \(\ell\)-th power for any \(\ell\in\mathcal L\).
Such an integer exists. If \(\mathcal L\) is empty, choose a
\(\pi\)-unit from the two successive integers. Otherwise, for an
odd prime \(\ell\), the integers with \(u/Q\) an \(\ell\)-th power
are exactly \(u=Q_\ell b^\ell\), \(b\ge1\), where \(Q_\ell\) is
the part of \(Q\) with prime exponents reduced modulo \(\ell\).
Their spacing is at least \(2^\ell-1\). Thus the fraction of an
interval of \(L\) integers lost to all these families is bounded
by
\[
             \sum_{\ell\in\mathcal L}\frac1{2^\ell-1}
                         +\frac{|\mathcal L|}{L}
               \le \frac4{21}+\frac{|\mathcal L|}{L}.
\]
The bound on the sum follows from
\((2^\ell-1)^{-1}\le(8/7)2^{-\ell}\) and the odd geometric series
starting at three. At most \(L/2+1\) additional integers are
multiples of \(\pi\). Their union has size less than \(L\) for
the stated value of \(L\).

Set \(x=u/c_*\), so \(0<x-F\le L/c_*\). Every \(c_i x\) is
an integer for \(i\in\mathcal T\). Its candidate offset
\(h_i=c_i x-j a_i\) is positive, since \(j a_i/c_i\le F\) for
each tier row by evaluating at a minimizing grid index. It is
a \(\pi\)-unit, as is \(c_i/c_*\), since \(\pi\mid j\) and \(\pi\nmid u\).
Thus a grid index samples the mark in a tier row precisely when
its factorial argument normalized by \(c_i\) is at least \(x\).
At an interior point attaining \(f(k)\), advancing one index
in the increasing direction of the attaining argument reaches
\(x\) by the stated sufficient test. One may instead verify this
reach directly at every interior index, which is the exact
forward-distance test. In either case at least one tier row
samples the mark somewhere on the grid. The mark-free indices
\(\{k:f(k)<x\}\) form a nonempty interval by convexity.
Three consecutive such indices would give a middle attaining
row reaching \(x\) at a neighboring index, so its width is at
most one.

## Isolation across scales

A literal realizing this mark is proportional to
\[
                       T^e+\frac{u}{j c_*}.
\]
Every one of its roots has \(\pi\)-valuation \(-\tau\), since
\(v_\pi(u)=0\) and \(c_*\) has the common valuation of the
maximizing slopes in block \(e\).
For a possible companion-order \(q\) literal in a scale \(f\)
with slope \(c_l\), the valuation of its root is at least
\(-v_\pi(c_l)/f\), strictly larger than \(-\tau\) if \(f\)
lies in the chosen parity class. This follows from \(\pi\nmid q\)
and from maximality of \(\tau\). For a literal of the marked
order in the same parity class, its root has valuation at least
\(-(v_\pi(j)+v_\pi(c_l))/f\), so only maximizing types can share
a root.

A literal from a different parity class cannot share any root:
raising the two negative real defining powers to the exponents
obtained after dividing by their gcd gives opposite signs. For
a maximizing type of a different degree \(f\) in the same
class, \(f<e\) by the choice of \(e\). If two monic literals
\(T^e+\alpha\) and \(T^f+\beta\) with positive rational constants
share a root and \(g=\gcd(e,f)\), then
\(\alpha^{f/g}=\beta^{e/g}\). The two powers \(f/g,e/g\)
are coprime and odd. It follows prime by prime that
\(\alpha\) is a rational \(\ell\)-th power for any prime
\(\ell\mid e/g\), and at least one such prime exists.
Here \(\alpha=u/(j c_*)\); this contradicts the selection of
\(u\). Only same-degree literals in \(\mathcal T\) with exact
normalized offset \(x\) can share the mark.

Include the primes of all nonzero literal resultants, the fixed
core ramps, order factorials and leading constants in a finite
set \(S\). Any outside prime of a selected mark now isolates
exactly the allowed marked-order factors, so the companion Newton
polygon is flat and the marked polygon has horizontal width at
most one. As in the
[unpaired test](MULTIVARIATE_UNPAIRED_ENVELOPE_MULTISCALE_MARKS.md#good-primes-in-a-fixed-template),
the effective polynomial \(S\)-part bound supplies an outside
prime for all large \(T\) if \(e>1\). For marks \(e=1\) in both
orders, companion separation makes the marks nonproportional.
They cannot both be \(S\)-smooth beyond the effective \(S\)-unit
height bound. These inputs are respectively
[Bugeaud, Evertse and Győry, Theorem 2.2](https://arxiv.org/pdf/1708.08290)
and
[Bérczes, Evertse and Győry, Theorem 2.1](https://math.leidenuniv.nl/reports/files/2008-16.pdf).
There are finitely many literal choices over the fixed template
and residual pairs, so the bounds are uniform. Even real
positivity and the Newton polygon product rule then preclude
any rational common factor; rational Bézout identities transfer
coprimality to every characteristic-zero field.
