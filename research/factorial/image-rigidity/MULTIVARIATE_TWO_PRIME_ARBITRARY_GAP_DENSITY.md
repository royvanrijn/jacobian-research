# Two-prime coprimality with one arbitrary increasing gap

Fix integers \(d\ge2\), \(B\ge0\) with \(B\equiv-2\bmod d\).
Consider one ordinary row \((dB,dB+d)\) and \(R\ge1\)
unit-gap rows \((a_i,a_i+1)\) with \(0\le a_i\le B\).
Equal rows are immaterial. The factorial moment polynomial
is irreducible at a set of orders of relative lower density
at least \(\log 2\) within
\[
                        n\equiv0\bmod \operatorname{lcm}(2,d).
\]
At each such order, both consecutive pairs containing it
are coprime. This is an arithmetic criterion for the
unaugmented monomials, for every fixed gap \(d\).

The proof uses the same
[two-prime mechanism as for unit gaps](MULTIVARIATE_TWO_PRIME_UNIT_GAP_DENSITY.md).
The congruence condition makes the nearby-prime edge's
denominator equal its full width even when its height is \(d\).

## The prime dividing the order

Write \(n=Kq\) with \(q\) prime and choose it so that
\(q>C\sqrt n\), where \(C\) is fixed sufficiently large
in terms of \(d,B,R\). For all large orders this ensures
\[
         q>d(B+1)K,\qquad q>d+R,\qquad
         q> d(d+R)/R,\qquad q\ge4d .
\]
All factorial arguments lie below \(q^2\). At index
\(k=hq+s\), \(0\le s<q\), the valuation relative to the
constant coefficient is
\[
        (d+R)h+\lfloor ds/q\rfloor+\mathbf1_{s>0}.
\]
The binomial term follows by counting the carry at the
\(q\) level, and the factorial terms follow since all starting
arguments are multiples of \(q\).

Subtract the supporting slope \((d+R)/q\) from the ordinates.
Within a block the unique minimum is at \(s=q-1\).
Indeed for \(s=q-e>0\) with \(e<q/d\) the shifted ordinate
is \(-R+(d+R)e/q\). For \(e\ge q/d\) it is greater than or
equal to \(-R+R/d\), which is strictly larger than the
value at \(e=1\) by the prime-size bound. The point
\(s=0\) lies higher as well. Thus the middle edge has slope
\((d+R)/q\) and width \((K-1)q\), while the final edge
has width one and slope \(R\).

Each earlier edge, in the first block, joins the right ends
of two different constant plateaus, or starts at index zero.
Every earlier slope is positive, so the polygon cannot touch
an earlier point of a constant plateau and then run to
its right end. The plateau endpoints have consecutive
distances at least \(q/d-2\), and each edge has integer
rise at most \(d\). Therefore its reduced denominator is
at least \(q/(2d^2)\). The middle slope has denominator \(q\).
The Newton polygon product rule consequently allows only
a linear factor among factors of positive degree less than
\(q/(2d^2)\).

## The nearby prime

Choose a prime \(p\le n\) with \(d<t=n-p\le\sqrt n\).
For all sufficiently large \(n\), there are no factorial
valuation contributions above the \(p\) level.
The binomial valuation is one precisely on \(t<k<p\).
A unit row makes its jump at \(k=p-a_i t\ge p-Bt\).
In the gap-\(d\) row the factorial floor increases from its
value at zero by
\[
                   \left\lfloor\frac{dBt+dk}{p}\right\rfloor.
\]
Each ordinal \(h\) before the final value \(d\) has the last
coefficient of total relative valuation \(h\) at
\[
                    k_h=\left\lfloor
                         \frac{hp-dBt-1}{d}\right\rfloor,
                       \qquad 1\le h\le d .
\]
At large \(n\) these points are strictly between \(t\) and \(p\).
There are no unit-row jumps before \(k_d=p-Bt-1\).
After this index all ordinates are at least \(d+1\), since
either the binomial term remains present or all the unit
rows have made their jumps at \(k\ge p\).

The plateau at zero ends at \(t\). Put
\[
                       w=k_d-t=p-(B+1)t-1 .
\]
For each \(h<d\),
\[
 h w-d(k_h-t)\ \ge\ (d-h)(B+1)t-h+1>0.
\]
Thus the earlier positive plateaus lie strictly above the
line from \((t,0)\) to \((k_d,d)\).
The later points also lie above it, because their slopes
from \((t,0)\) are at least \((d+1)/p>d/w\) for large \(n\).
This gives an edge of height \(d\) and width \(w\).
Since \(n\equiv0\bmod d\) and \(B+2\equiv0\bmod d\),
\[
                       w=n-(B+2)t-1\equiv-1\bmod d.
\]
The slope \(d/w\) is reduced. Some rational irreducible
factor therefore has degree at least \(w\).

Choose \(C>2d^2(B+3)\) in the first prime condition.
The complement of this large factor has degree
at most \(n-w=(B+2)t+1<q/(2d^2)\). The polygon at \(q\)
leaves only a possible linear complement. At even \(n\)
the factorial moment polynomial has no real zero by the
positive Gamma integral, so the complement is constant.
The nondivisibility argument for two monomials then gives
coprimality with the following order; the preceding order
has insufficient degree to contain the irreducible factor.

## Density

For fixed modulus \(D=\operatorname{lcm}(2,d)\), write
\(n=Dm\) in this progression. A prime \(q>C\sqrt{Dm}\)
does not divide \(D\) for large \(m\), and hence divides
\(n\) precisely when it divides \(m\).
The prime-divisor density is therefore \(\log 2\) on this
progression, by
[Filaseta and Rocks, equation (3.2)](https://people.math.sc.edu/filaseta/papers/ColloquiumMathPaper2016.pdf).
The near-prime condition has relative density one by their
equation (3.3). Excluding \(0\le n-p\le d\) removes a subset
of finitely many translates of the primes, of density zero.
The intersection of these conditions proves the claimed
density.

