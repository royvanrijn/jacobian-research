# Pure winding-one cubic phases with zero leading center coefficient

Let \(K\) be a field of characteristic zero. In the sphere algebra
\(A=K[x,y,t]/(t^2+2xy-1)\), consider a regular winding-one cubic
\[
 f=\frac{C(x)+B(x)t+D(x)t^2+E(x)t^3}{x},\qquad E\ne0,
\]
with \(C,B,D,E\in K[x]\). Regularity requires \(C(0)+D(0)=B(0)+E(0)=0\).
Define \(\Lambda(f)\) by substituting \(y=(1-t^2)/(2x)\), taking the
constant term in \(x\) and applying the interval functional
\(t^j\mapsto 1/(j+1)\) for even \(j\), and \(0\) for odd \(j\).

**Theorem.** Suppose \(E(0)=0\). The positive pure moments
\(\Lambda(f^m)\) vanish eventually if and only if all four coefficient
polynomials are divisible by \(x^2\). In that case every positive pure
moment vanishes, and \(\Lambda(g f^m)=0\) for every fixed \(g\in A\)
and all sufficiently large \(m\). There is no bound on the sphere degree.

## The upper face with negative angular exponent

Write \(S\) for the support of the literal Laurent expansion
\(\sum_{i,j} f_{ij} x^it^j\). The exponent \(i\) is the angular exponent
and \(j\) is the height. Suppose \(D(0)\ne0\). Regularity makes the
negative angular part exactly \(x^{-1}D(0)(t^2-1)\). At all nonnegative
angular exponents the height is at most three, and height three occurs
because \(E\ne0\) and \(x\mid E\).

If \(S\) has a strictly positive angular exponent, let \((0,H)\) be the
uppermost point of \(\operatorname{conv}(S)\) at angular exponent zero.
The hull has a supporting line \(j+s i=H\) there, since it has points on
both sides of that axis. Also \(H>0\) by the segment joining \((-1,2)\)
to any point of positive angular exponent. Only the occupied nodes on
this line can contribute to \([t^{mH}]\operatorname{CT}_x f^m\) when
\(mH\) is integral, and the angular constant term has height at most
\(mH\).

If the occupied nodes on the line have both angular signs, the negative
node is \((-1,2)\), giving \(H=2-s\). Since an occupied point of
height three has nonnegative angle, \(s\ge0\) would give a weight of
at least three at that point, above the supporting line. Thus \(s<0\).
The equation \(j=2-s(1+i)>2\) forces height three at every other
occupied node of the line, because \(i\ge0\) there. A nonhorizontal
line contains at most one such node. The entire occupied face therefore
consists of \((-1,2)\) and \((a,3)\) for an integer \(a\ge1\), with
\(H=(2a+3)/(a+1)\).

If the occupied face has only one strict angular sign, or has none, the
node \((0,H)\) is occupied and only its coefficient contributes to
the top angular-constant coefficient. The same conclusion holds when
there are no positive angular exponents in \(S\). In that case the
height-three node is at angular exponent zero, and the constant term of
each power is the power of this angular core.

## The prime test

First work over a number field. Exclude the finitely many rational primes
at which a Laurent coefficient is nonintegral or a specified leading
coefficient fails to be a unit. In the two-node case put \(N=2a+3\)
and take good primes \(p\equiv1\pmod{2N}\). For \(k=(p-1)/N\) and
\(m=(a+1)k\), the top possible angular-constant height is \(p-1\).
Its coefficient is
\[
 \binom{(a+1)k}{k}\, f_{-1,2}^{ak} f_{a,3}^{k}.
\]
The factorial indices are smaller than \(p\), because
\((a+1)k\le Nk=p-1\). This coefficient is a unit at every prime ideal
above each such good \(p\).

In the core case \(H\) is a positive integer and
\(m=(p-1)/H\) for good primes \(p\equiv1\pmod{2H}\) gives the same
top possible height \(p-1\), now with coefficient \(f_{0,H}^m\).
In both cases \(p\) is odd, so the term of height \(p-1\) survives
interval integration. For all smaller heights the interval denominators
are prime to \(p\). Consequently \(p\Lambda(f^m)\) reduces to the
nonzero top coefficient at any prime ideal above \(p\). The needed primes
exist in arbitrarily large size by the
[prime progression argument](SPHERICAL_LEADING_BINARY_FORM_MOMENT_OBSTRUCTIONS.md#3-a-leading-coefficient-detected-at-good-primes).

Over arbitrary \(K\), generate a finitely generated \(\mathbb Q\)-subalgebra
from the Laurent coefficients and invert every occupied coefficient. One
maximal-ideal quotient is a number field by Zariski's lemma. This single
specialization preserves the support and all power moments, since each
moment is a rational polynomial in the coefficients. The specialized
nonvanishing at infinitely many orders implies nonvanishing at the same
orders over \(K\). This excludes every case with \(D(0)\ne0\).

## The zero center and safety

When \(D(0)=E(0)=0\), regularity gives \(C(0)=B(0)=0\), so
\(f\in K[x,t]\). The existing
[zero-center criterion](MIXED_CUBIC_ZERO_CENTER_AND_CONDITIONAL_OFFCENTER.md#1-classical-interval-rigidity-and-the-complete-zero-center-branch)
identifies purity with divisibility by \(x\) in that ring. The same
criterion applies to eventual purity. Indeed all pure moments are the
interval moments of \(h(t)=f(0,t)\). A nonzero constant \(h\) gives
nonzero moments at every order. For nonconstant \(h\), the core prime
test applied to its highest coefficient gives infinitely many nonzero
moments.

Thus eventual purity forces \(x\mid f\) in \(K[x,t]\), or equivalently
\(x^2\mid C,B,D,E\). Conversely every such phase has strictly positive
angular support. Every fixed marker has a lower bound on its Laurent
angular exponents, so its product with \(f^m\) has no angular constant
term for sufficiently large \(m\). The proof of the nonzero
\(D(0)\) exclusion imposes no assumptions on the cubic invariant.
