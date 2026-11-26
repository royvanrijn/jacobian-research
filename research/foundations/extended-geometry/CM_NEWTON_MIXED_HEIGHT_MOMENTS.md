# A mixed-height certificate for moments divisible by the residue prime

At every fixed odd prime \(p\), every positive moment order \(m\)
divisible by \(p\) has supporting-face certificates with unbounded
connected cyclic-cover conductor. The top angular face supplies
a directed coherent Frobenius coefficient in order \(m-1\).
A single lower radial term moves that coefficient onto a large
interval denominator in order \(m\). Its multiplicity is \(m\),
and the certificate divides out its exact \(p\)-valuation.
The primitive factor is nonordinary at the certifying prime, and
the character-one return map is zero.

The Frobenius input is the ordinary Kummer matrix construction in
[Bouw, §§4–5](https://doi.org/10.1023/A:1017513122376) and
[Lin–Mantovan–Singhal, §2.4](https://arxiv.org/html/2406.09632v1#S2.SS4).
The [prescribed prime-to-\(p\) construction](CM_NEWTON_PRESCRIBED_MOMENT_EIGENLINES.md)
supplies the directed coefficient and the zero map in its return
orbit. The added input here is a higher interval-residue identity
and a choice of subleading radial height.

## Division by the \(p\)-part of the power

Work in the literal Laurent chart of \(t^2+2xy=1\), with coefficients
integral at a prime above an odd \(p\). For \(m\ge1\) write
\(d=v_p(m)\), and write \(\partial_E=t\,\partial/\partial t\) for
the derivative in that chart with \(x\) fixed. If \(p\nmid j\),
the coefficient identity
\[
 j[x^0t^j]f^m
  =m[x^0t^j](\partial_E f) f^{m-1}
\]
shows that the coefficient on the left is divisible by \(p^d\).
Its quotient by \(p^d\) is governed by the residue of a derivative
times the power \(m-1\), even when every first-residue coefficient
at a pole denominator vanishes by Frobenius.

More generally, suppose that all angular-constant heights
satisfy \(j+1<p^{R+1}\) with \(R>d\). The normalized moment is
integral at \(p\), and
\[
 p^{R-d}\mathcal A(f^m)\equiv
 -\frac{m}{p^d}
   \sum_{\substack{1\le h<p\\ h\ \mathrm{odd}}}
   \frac{[x^0t^{h p^R-1}](\partial_E f)f^{m-1}}{h}
 \pmod p .
\]
Indeed integration divides by \(j+1\) at even height and kills
odd height. If a denominator is divisible by \(p\), then \(j\)
is a unit and the displayed coefficient identity applies.
After normalization, all terms with denominator valuation below
\(R\) vanish modulo \(p\); unit denominators vanish since \(R>d\).
At \(j=h p^R-1\), the factor \(j^{-1}\) reduces to \(-1\).

## A subleading radial height and the highest crossing face

Let \(p\mid m\), \(k=m-1\), and \(d=v_p(m)\). Fix a polynomial
\(P\) with \(p\)-integral coefficients and the root pattern from the
[prescribed directed coefficient theorem](CM_NEWTON_PRESCRIBED_MOMENT_EIGENLINES.md#every-prime-to-p-moment-at-a-fixed-odd-prime)
at order \(k\). Put \(\deg P=b+1\). Its coefficient
\[
                         C=[z^{bk}]P(z)^k
\]
is a \(p\)-unit. The same construction provides a nonempty open
interval \(I\subset(0,1)\) of character ratios where the reduced
one-step Frobenius map between the coherent lines is zero.

Choose positive \(\lambda\) and an odd \(s\equiv1\bmod 2kp\) with
\[
                  k p^\lambda/s-k\in I .
\]
Existence follows from the increasing length of the corresponding
interval for \(s\) as \(\lambda\) increases. Take large \(r\)
with \(p^r\equiv1\bmod k\), and put
\[
 n=\frac{s p^r-1}{k},\qquad R=r+\lambda,\qquad
 J=p^R-s p^r,\qquad \delta=n-J=mn+1-p^R .
\]
The ratio \(J/n\) tends to \(k p^\lambda/s-k\in I\), hence
\(0<J<n\) for large \(r\). Also \(J\) is even and
\(\delta\equiv2\bmod p\), because \(kn\equiv-1\bmod p\) and
\(k\equiv-1\bmod p\). In particular \(\delta\) is a unit.

Give \(x^it^j\) weight \(j+2pi\). Let \(g\) have highest weight
\(n\) with angular polynomial \(z^{-b}P(z)\), and arbitrary
strictly lower terms integral at \(p\). Leave out the distinguished
free coefficient of \(t^J\), and set \(f=g+c t^J\). Its highest
crossing face is unchanged. For large \(n\) all equality
monomials are regular liftable. At a negative angle \(-v\),
\[
 (-1)^v x^{-v}t^n(1-t^{2p})^v
 =
 (-2)^v y^v t^n(1+t^2+\cdots+t^{2p-2})^v
\]
has highest weight \(n\) with the prescribed leading coefficient.
All its literal weights are congruent to \(n\bmod p\).
Positive angles have nonnegative heights because \(n\ge2p\).
The equality lift has least degree \(n+(2p-1)b\), which is also
the phase degree if the lower normal terms have degree at most this.

The largest angular-constant height in \(f^m\) is at most \(mn\).
Since \(mn+1=p^R+\delta<2p^R\), the higher interval-residue
identity leaves only height \(p^R-1\). If at least two factors of
a term come from \(c t^J\), the possible angular-constant heights
are at most \(mn-2\delta<p^R-1\). With exactly one radial
factor the target height is achieved by the top angular constant
coefficient of \(g^k\). This is \(C\), since
\(J+kn=p^R-1\). Thus the exact coefficient at that height is
\[
 [x^0t^{p^R-1}]f^m
   = [x^0t^{p^R-1}]g^m+ m c C .
\]
Using the derivative identity on the first term proves
\[
 p^{R-d}\mathcal A(f^m)
 \equiv \frac{m}{p^d}
 \left(c C-
 [x^0t^{p^R-1}](\partial_E g)g^{m-1}\right)\pmod p .
\]
Consequently **every lower-weight choice in \(g\)** admits
all but one residue class of the distinguished radial coefficient
\(c\). Since \(p\) is odd, one of the two units \(1,2\) always
certifies the \(m\)-th moment. Only products whose total
deficiency from the top weight is \(\delta\) can affect the
exceptional class. In particular terms of weight strictly less
than \(J\) cannot affect it.

A simple condition makes every unit \(c\) work. If the reduction
of all literal support of \(g\) has weight congruent to \(n\bmod p\),
the derivative term in the last congruence is zero: the derivative
preserves weight, so the angular-constant height in that product
would be congruent to \(mn\equiv0\bmod p\), rather than to
\(p^R-1\equiv-1\bmod p\). One may still add arbitrary strict
lower terms with \(p\)-divisible coefficients, since they disappear
from this derivative product on reduction. The precision-one
allowance here follows from dividing out \(p^d\) by the derivative
identity; it does not require the lower terms to vanish to
precision \(d+1\).

The connected Kummer cover \(w^n=z^{n-b}P(z)\) belongs to the
actual highest crossing face. At the certifying prime, the
\(r\)-fold coherent Frobenius sends the overlap class of
\(w^s/z^s\) to \(C[w/z]\). The exponent relation is
\(s p^r=kn+1\), as in the
[transition formula](CM_NEWTON_SELECTED_EIGENLINES.md#a-moment-detecting-a-transition-rather-than-a-return).
On the other hand the primitive character represented by
\(p^R\bmod n\) has ratio \(p^R/n-k=(J+1)/n\) in \(I\) for
sufficiently large \(r\). It lies
in the character-one Frobenius orbit and has a zero outgoing
map between coherent lines. The return on character one is
therefore zero, and the primitive factor is nonordinary.
The selected incoming map still certifies the moment through
the mixed-height coefficient.

## A dominant radial variant

There is also an isolation with a freely prescribed unit radial
coefficient and no weight congruence or residue exception on \(g\).
Start with any of the prime-to-\(p\) crossing faces just used,
with highest weight \(n\) and unit coefficient \(C\) at order \(k\).
Put \(S=nk+1\), choose \(R\) with \(p^{R-d}>mS\), and set
\(J=p^R-S\). The radial term in \(f=g+c t^J\) is now dominant.

In a summand with \(l\) radial factors and angular-constant
height \(h\) from the other factors, the interval denominator is
\[
 j+1=l p^R+(h+1-l S).
\]
We have \(0\le h\le(m-l)n\). The offset in parentheses vanishes
only for \(l=1,h=nk\), and its absolute value is less than
\(p^{R-d}\) otherwise. Thus \(p^{R-d}\) annihilates all
other terms modulo \(p\), regardless of their coefficients.
The surviving term is again \((m/p^d)c C\).
This variant preserves the unrestricted lower coefficients and
the connected auxiliary-cover certificate, but its crossing face
is subleading rather than topmost.

At \(p=2\) all the denominators from the even interval heights
are odd. The pole-isolation arguments here concern odd primes
and give no analogous two-adic certificate.

## Exact replay

The [checker](verify_cm_newton_mixed_height_moments.py) evaluates
the selected Frobenius coefficients, confirms a zero map in each
return orbit, and expands the mixed moments over \(\mathbb Q\).
It passes the leading-face examples \((p,m)=(3,3),(5,5),(3,9)\)
at connected orders \(289,18401,21115\), all with normalized
residue 1. The third tests \(v_p(m)=2\).

For an unrestricted lower choice at order 289 it computes the
derivative residue and finds that radial coefficient 1 gives
normalized residue zero while radial coefficient 2 gives one.
The dominant-radial variant of the same unrestricted lower choice
gives residue one with radial coefficient 1. The replay uses the
auxiliary polynomial and character routines in
[the prescribed-moment checker](verify_cm_newton_prescribed_moment_eigenlines.py).

    python3 -I verify_cm_newton_mixed_height_moments.py
