# An exact lower pole for supercritical sphere moments

There are fixed phases with a supercritical weighted angular root
for which the top interval-pole residue vanishes and the bounded-phase
volume diverges. In the family below the unmarked moment at the same
order still has a pole of depth \(r/2\) at primes from a set of
density \(1/2\). A polynomial change of sphere coordinates reduces
the entire moment to a beta integral, so its exact depth can be
calculated from base-\(p\) carries.

## The phase and its exact moments

Let \(d\ge2\), \(B\ge1\), \(b=dB\) and \(n=(2d-1)B=2b-B\).
Take an integer \(1\le a<b\), and take \(c,\alpha,\beta\ne0\) with
\(\alpha\ne\beta\) in a characteristic-zero field. Put \(w=x^d-ct\)
and define the phase on \(2xy+t^2=1\) by its literal Laurent form
\[
  f=x^{-b}(w-\alpha x-c)^a(w-\alpha x+c)^{b-a}
            (w-\beta x-c)^{b-a}(w-\beta x+c)^a.
\]
For the numerator, the constant term in \(x\) of each linear factor
is a nonzero scalar times \(t+1\) or \(t-1\). There are \(b\) of
each sign. A contribution to the coefficient of \(x^j\) uses at
most \(j\) nonconstant factors, so for \(j<b\) this coefficient is
divisible by \((1-t^2)^{b-j}\). This is the factor needed to replace
the negative angular power \(x^{j-b}\) by a polynomial in \(y,t\).
Thus \(f\) is regular on the sphere.

For the weight \(j+i/d\), the support lies in
\([n-2b,n]\). The top part comes from choosing \(w\) in all
\(2b\) factors. Its angular polynomial is
\[
                       z^{-b}(z^d-c)^{2b},
                     \qquad z=x/t^{1/d}.
\]
Each of the \(d\) angular roots has multiplicity \(2b=n+B>n\).
The coefficient of \(x^{(2d-1)b}\) is one. The degree bound
\((2d-1)b\) also holds for the polynomial in the normal sphere
basis, since each numerator factor has ordinary degree at most
\(d\) and replacing a regular negative angle by powers of \(y\)
preserves its degree bound. Hence the least sphere degree is
\((2d-1)b=dn\). For \(a=1\) the phase is an absolute nonpower, as
it has a simple zero on the prime divisor \(w-\alpha x=c\)
at its generic point in the chart \(x\ne0\).

Let \(\Lambda\) be angular constant term followed by the normalized
interval integral.

**Proposition (Exact moments).** For every \(M\ge1\) one has
\[
 \boxed{\displaystyle
 \Lambda(f^M)=
 \frac{(\alpha-\beta)^{bM}(-1)^{aM}(2c)^{bM}}
 {(bM+1)\binom{bM}{aM}} .}
\]
To prove this, use the polynomial change of sphere coordinates
\[
 T=t+\phi(x),\quad Y=y-\frac{\phi(x)}x t-\frac{\phi(x)^2}{2x},
 \qquad \phi(x)=\frac{\alpha x-x^d}{c}.
\]
It fixes \(x\), preserves the sphere equation and preserves its
volume form \(dx/x\wedge dt\) on the Laurent chart. Scaling
\(\phi\) from zero to its stated value gives an isotopy of
automorphisms. After a coefficient-field embedding into
\(\mathbb C\), the real sphere cycle and its transforms are
homologous, and the holomorphic volume forms multiplied by powers
of a regular phase are closed. Hence the change preserves all
sphere moments. The identity then holds over the original
characteristic-zero field. [Sphere residue formula](GLOBAL_SPHERE_PURITY_RECOGNITION.md#1-a-closed-rational-period)

In the new coordinates, \(W=w-\alpha x=-cT\). The last two
factors in \(f\) become \((W+(\alpha-\beta)x-c)^{b-a}\) and
\((W+(\alpha-\beta)x+c)^a\). Their combined degree in \(x\) is
exactly \(b\), so after the original \(x^{-b}\) all angular
exponents are nonpositive. Only the constant angular coefficient
\((\alpha-\beta)^b(W-c)^a(W+c)^{b-a}\) contributes to the
constant term in any power. The normalized integral over
\(W\) from \(-c\) to \(c\) gives the stated beta formula.

## The vanished top residue and the exact lower depth

Work at an odd rational prime \(p\) with \(p\nmid dn(d-1)\)
where \(c,\alpha,\beta\) are integral and \(c(\alpha-\beta)\) is a
unit. Fix \(r\ge1\) such that \(q=p^r\equiv1\pmod{dn}\) and put
\(M=(q-1)/n\). In the angular constant term of \(f^M\), the
height is at most \(q-1\) and lies within \(2bM<2q\) of this
upper height. As in the [critical-face calculation](SPHERICAL_CRITICAL_FACE_P_POWER_MOMENTS.md),
the only denominator divisible by \(q\) is the top denominator
\(q\). Its numerator reduces to
\[
             [Z^{BM}](Z-c)^{2bM}
                 =[Z^{BM}](Z-c)^{q+BM-1}\equiv0\pmod p.
\]
Indeed \(BM<q\), and \((Z-c)^q\) contributes only its constant
term while the other factor has degree \(BM-1\). Thus
\(q\Lambda(f^M)\equiv0\pmod p\).

Write \([u]_n\) for the least nonnegative residue of an integer
modulo \(n\).

**Theorem (Secondary pole).** At the same moment order the exact
valuation is
\[
 \boxed{\displaystyle
 v_p\bigl(\Lambda(f^M)\bigr)=-\#\{0\le j<r:
              [a p^j]_n+[(b-a)p^j]_n\ge n\}.}
\]
The factor \(bM+1=(d q+d-1)/(2d-1)\) is a \(p\)-unit. For
\(0<A<n\) and \(1\le i\le r\), one has
\[
                \left\lfloor\frac{AM}{p^i}\right\rfloor
                     =\left\lfloor\frac{A p^{r-i}}{n}\right\rfloor .
\]
In fact \(AM/p^i=A p^{r-i}/n-A/(np^i)\). The residue
\(R=[A p^{r-i}]_n\) satisfies \(R p^i\equiv A\pmod n\) and
\(R p^i\ge A\), since \(0<A<n\) and \(R>0\). The subtraction
therefore leaves the floor unchanged. Counting multiples of each
\(p^i\) in the three factorials of \(\binom{bM}{aM}\) gives the
formula. Its summand is one exactly when the two residues cross
\(n\), the carry behind [Kummer's binomial valuation](https://doi.org/10.1515/crll.1852.44.93).

When \(-1\in\langle p\rangle\subset(\mathbb Z/n)^\times\), exactly
half of the terms are carries. Pair an orbit element with its
negative. Each residue is nonzero and the two residues cannot
sum to \(n\), as \(b<n\) and the orbit elements are units. Thus
exactly one of each pair crosses \(n\). For every admissible \(r\)
in this prime class, the same original moment order has
\[
               v_p\bigl(\Lambda(f^{(p^r-1)/n})\bigr)=-r/2.
\]
**Corollary (Prime density).** Taking \(d=2\), \(B=3^{h-1}\)
and \(a=1\) gives an unbounded family of absolute nonpowers with
\(n=3^h\) having exact pole depth \(r/2\) on a density-one-half
set of primes. The group of units
modulo \(3^h\) is cyclic of order \(2\cdot3^{h-1}\), so
\(-1\in\langle p\rangle\) exactly when \(p\equiv2\pmod3\).
These primes have natural density \(1/2\), and each good prime
has infinitely many admissible positive exponents \(r\). Over a
number field the valuations are normalized by \(v_p(p)=1\) at
the good coefficient places; only finitely many rational primes
are excluded by the fixed coefficients.

## Two divergent clusters

The preceding moment conclusion applies even though the bounded
phase volume diverges. In the original shear chart \((x,w)\),
the support of \(f\) has both extreme angles \(-b,b\). On the
region where \(x\) tends to infinity and \(W=w-\alpha x\) stays
in a fixed disk, the first group of factors is bounded and the
second has growth at most \(|x|^b\), canceled by the denominator.
In these coordinates the volume form is a nonzero constant times
\(dx/x\wedge dW\), whose squared integral diverges there. The
same reasoning applies separately to \(w-\beta x\) bounded.
Both regions approach the nonzero supercritical angular roots
\(z^d=c\). The exact lower pole therefore applies at a face
where the [Puiseux volume test](SPHERICAL_WEIGHTED_PUISEUX_VOLUME.md)
detects divergence.

For an exact replay, take \(c=1,\alpha=1,\beta=2,a=1,d=2\).
At \(B=1,p=5,r=2\) the order is \(8\), and literal Laurent
multiplication gives \(5\Lambda(f^8)\equiv2\pmod5\) while
\(25\Lambda(f^8)\equiv0\pmod5\). The replay also checks the
valuation formula and the density class using finite residue orbits.
[Exact replay](verify_spherical_supercritical_second_pole.py)
