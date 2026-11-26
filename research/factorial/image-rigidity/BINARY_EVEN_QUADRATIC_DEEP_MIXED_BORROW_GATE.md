# A binary borrow constraint in the deep mixed quadratic sector

Normalize the constant coefficient to one and consider
\(f=ax^2+bxy+cy^2+1\), with \(ab c\ne0\) and
the usual factorial functional \(\Lambda_2\).
Suppose the moments at \(N-3,N-2,N-1,N\) all vanish.
Write \(v=v_2\) with \(v(2)=1\).

**Theorem.** In the strict mixed sector
\(v(4a)\le-3,\ v(c)>-1\), put \(\tau=v_2(N)\).
Then \(\tau\ge4\) and \(v(4a)=1-\tau\).
Set \(L=2^\tau\), \(W=(\tau-3)L\) and \(Y=N/L-1\).
The two-state binary minimum
\[
       \min_{0\le q\le Y}
            \bigl[Wq+2s_2(q)-\operatorname{carr}_2(q,Y-q)\bigr]
                                                               \tag{1}
\]
must be attained at least twice. Here \(\operatorname{carr}_2\)
counts carries in binary addition. The minimum and whether it
ties are decidable in time linear in the binary length of \(Y\),
with integer operations on that length.

In particular, for some \(q\ge1\), starting at bit
\(h=\lfloor\log_2 q\rfloor+1\), the binary expansion of \(Y\)
must contain a block of at least
\(Wq+2s_2(q)-h\) zeros before a set bit.
For fixed \(\tau\), among the possible even \(Y\) with
this form the upper proportion of starting indices that
can still allow such a four-zero block is at most
\[
                  2\sum_{q\ge1}2^{-(Wq+2s_2(q)-h)}
                   \ \le\ \frac{2^{-W}}{1-2^{-(W-1)}} .
\]
The assertions also hold with \(a,c\) exchanged. The
argument does not exclude genuine ties in (1).

## The plane moments and their exact reverse ordinates

Write \(r=v(4a)\le-3\) and let \(M_j^0\) be the moment
on the comparison plane \(b=2a\) with the same \(a,c\).
By the [coupled sector calculation, (7)–(10)](BINARY_EVEN_QUADRATIC_COUPLED_TWO_ADIC_SECTORS.md#a-dominant-axis-and-approximation-by-its-single-face-plane),
the two boundary zeros at \(N-1,N\) force \(N\) even,
\(r=1-\tau\), and
\[
        v(M_j-M_j^0)\ \ge\ (j-N)r+s_2(N)-s_2(j).          \tag{3}
\]
Here \(s_2\) denotes binary digit sum. The equality
\(r=1-\tau\) uses \(v(c)>-1\) to make both relevant
\(B_j(c)\) axis values units. In particular
\(\tau\ge4\), and (3) is nonnegative at \(j=N-1,N-2\).
For \(j=N-1\) its right side is zero; for \(j=N-2\)
it is \(\tau\).

The exact degree-\(q\) contribution of the quadratic
part on the plane \(b=2a\) is
\[
             (q!)^2\sum_{s=0}^q
                       \binom{2s}{s}(4a)^{q-s}c^s .
\]
Since \(v(c)>r\), the term with \(s=0\) is uniquely
lowest in this sum. Consequently, dividing the
moment \(M_j^0\) by \((j!)^2(4a)^j\), its coefficient
group of reverse degree \(d=j-q\) has exact valuation
\[
 V_j(d)= w d+s_2(d)+s_2(j)-s_2(j-d),\qquad 0\le d\le j.
                                                               \tag{4}
\]
Indeed the ratio of the dominant degree-\(j-d\) term to
the dominant degree-\(j\) term is
\((j-d)!/[j!\,d!\,(4a)^d]\), whose valuation is (4).
The leading normalizer has valuation
\(-wj-2s_2(j)<0\), whereas \(V_j(0)=0\).
If the minimum of (4) were unique for either of
\(j=N-1,N-2\), the plane moment would therefore
have strictly negative valuation. The error (3)
could not cancel it. Thus both minima have ties.

## The two low blocks have the same minimum

Set \(w=\tau-3\), so \(W=wL\).
For \(j=N-1=LY+L-1\), write \(d=Lq+b\), \(0\le b<L\).
There is no borrow within the top part when subtracting \(b\).
Thus \(V_j(Lq+b)-V_j(Lq)=wb+2s_2(b)>0\) for \(b>0\).
The minimum is at an exact multiple of \(L\).

For \(j=N-2=LY+L-2\), the same conclusion holds.
If \(0<b\le L-2\), split \(b=2z\) or \(2z+1\).
The increments of (4) from \(Lq\) are respectively
\(wb+2s_2(z)\) and \(wb+s_2(z)+s_2(z+1)\), both positive.
For the borrowing case \(b=L-1\), write \(d=Lq-1\)
with \(1\le q\le Y\). Comparing to \(N-1\) at \(Lq\)
gives
\[
 V_{N-2}(Lq-1)=V_{N-1}(Lq)+1+v_2(q).
\]
This follows from \(s_2(Lq)-s_2(Lq-1)=1-\tau-v_2(q)\)
and \((N-2)-(Lq-1)=(N-1)-Lq\). Also
\(V_{N-1}(Lq)=V_{N-2}(Lq)\), since the low digit sums of
\(j\) and \(j-Lq\) cancel in either case. Thus the borrowing
case is strictly above the term at \(Lq\).

At exact multiples of \(L\) the ordinates for both orders
are \(Wq+s_2(q)+s_2(Y)-s_2(Y-q)\).
Each binary carry decreases the digit sum in addition by
one, giving exactly (1). Its minimum must tie by (3)–(4).
A positive minimizer exists since \(q=0\) has value zero.

## A binary dynamic program and the sparse-index bound

Process the digits of \(Y\) from low to high. If its bit
at \(i\) is \(Y_i\), the two addend digits \(x,y\) and
the incoming and outgoing carries \(c_i,c_{i+1}\)
must satisfy \(x+y+c_i=Y_i+2c_{i+1}\).
The additive cost at that bit is
\[
                       W\,2^i x +2x-c_{i+1}.
\]
Start and end with carry zero. Retain the minimum cost
and the number of tied paths, truncated at two, at
each carry state. A pair of addends is determined by
the first addend and \(Y\), so this detects exactly the
stated ties. The number of states and transitions per
bit is bounded independently of the moment order.

For a positive minimizer \(q\), the value in (1) is
nonpositive. There must be at least \(Wq+2s_2(q)\)
carries. At most \(h\) of these occur below position \(h\).
At and above \(h\) all bits of \(q\) are zero. There an
incoming carry propagates as long as \(Y\) has bit zero,
and it stops at the first set bit. Once stopped it cannot
restart. This proves the claimed consecutive-zero condition.

To count possible indices, among \(Y\le X\) there are
at most \(1+\log_2 X\) carries, so only
\(q\le(1+\log_2 X)/W\) need be considered.
For each \(q\), the required zero block is a congruence
condition of density \(2^{-(Wq+2s_2(q)-h)}\). Since
\(h\ge1\), imposing even \(Y\) does not change this
density within the even integers. The counting error
for \(Y\le X\) is at most \(2^h\), so the total error
over eligible \(q\) is \(O((1+\log X)^2)\).
Taking the union bound gives the displayed upper
proportion. For \(q=1\) the exponent is \(W+1\),
and for \(q\ge2\) it is at least \((W-1)q+2\);
their geometric sum is below the stated bound.

There are exact ties in this optimization. For
example \(Y=2^{W+2}\) gives value zero at both
\(q=0\) and \(q=1\). For \(q>1\) the cost is
positive: a carry count is at most \(W+2\),
whereas \(Wq+2s_2(q)>W+2\).
Thus the valuation-only criterion still leaves a
substantive cancellation problem.
