# Local tests for the odd-start quadratic Gaussian window

For a positive order \(q\), write
\[
 G_q(A,B,C)=\mathbb E(A+B X_q+C X_q^2)^q,\qquad
                        X_q\sim\mathcal N(0,1/q).
\]
These are rational homogeneous polynomials of degree \(q\).
Nonvanishing of the three consecutive forms on projective
coefficient space gives the
[large-height exclusion and fixed-support decision](SFC_INDEPENDENT_COMMON_HEIGHT_DECISION.md#3-three-bivariate-columns-of-equal-degree)
for homogeneous bivariate trinomials at the same start. The
[even-start theorem](SFC_COLLINEAR_TERNARY_EVEN_START_AND_DYADIC_REDUCTION.md)
handles all even starts. Here \(n\) is even and the window
consists of \(n-1,n,n+1\).

## 1. Degree and the odd-pair resultant

**Theorem.** For every positive odd \(j\), the projective curves
\(G_j=0\) and \(G_{j+2}=0\) have no common component. In the
chart \(B=1\), the resultant in \(x=A/B\) of the two forms,
as a polynomial in \(h=C/B\), has a zero of exact order
three at \(h=0\). Every closed point of the odd-pair
intersection has degree at most \(j(j+2)\) over \(\mathbb Q\).

Write \(F_j(x,h)=\mathbb E(X_j+x+hX_j^2)^j\). It is monic in
\(x\). At \(h=0\) its roots are \(\mathrm i\) times the scaled
Hermite nodes. Each is simple, and the two odd orders share
only the zero node by
[scaled Hermite isolation](SFC_SCALED_HERMITE_ROOT_ISOLATION.md).
At that node the derivative in \(x\) is
\(j\mathbb E X_j^{j-1}\ne0\).
The ratio of the derivatives in \(h\) and \(x\) is one,
so the simple root has an expansion
\[
                x_j(h)=-h+\frac{2(j-1)}{3j}h^3+O(h^5).
\]
Indeed the parity identity \(F_j(-x,-h)=-F_j(x,h)\) rules
out even powers. For the cubic coefficient substitute
\(x=-h+\gamma h^3\) and compare the terms of order three
after dividing by \(j\mathbb E X_j^{j-1}\). The cubic term
from \(h(X_j^2-1)\) is
\[
 \frac{\binom j3}{j}
 \frac{\mathbb E X_j^{j-3}(X_j^2-1)^3}{\mathbb E X_j^{j-1}}
                  =-\frac{2(j-1)}{3j}.
\]
For \(j\ge3\), under the tilt by \(X_j^{j-3}\) the consecutive
first three moments of \(X_j^2\) are
\((j-2)/j,(j-2)/j,(j^2-4)/j^2\), giving the displayed value.
For \(j=1\) the cubic term is absent and the root is \(-h\).

The cubic coefficients at \(j\) and \(j+2\) differ by
\(4/(3j(j+2))\). The monic resultant is the product of the
cross-order root differences. All pairs other than the two
zero branches have nonzero differences at \(h=0\).
This proves exact order three and excludes a common component
meeting \(B\ne0\). The line \(B=0\) is not a component because
the first form is nonzero at \([A:B:C]=[1:0:0]\).
The degree bound is plane Bézout, applied to the proper
intersection and then to the full Galois orbit of any point.

## 2. The pure quadratic direction at intermediate primes

**Theorem.** At every prime \(n+1<p<2n-1\), a common zero
of \(G_{n-1},G_n,G_{n+1}\) over \(\overline{\mathbb Q}_p\)
cannot reduce to \([0:0:1]\).

More generally put \(P=(p-1)/2\), and consider at least three
distinct orders in \(P<q<p\). If a common zero reduces to
\([0:0:1]\), all the orders other than the least one must
take the same value under \(q\mapsto 4q^2+q\bmod p\).
Consequently four orders in this interval already rule out
that reduction, and for three orders the two upper orders
\(q_1,q_2\) must satisfy \(4(q_1+q_2)+1=0\bmod p\).

Normalize \(C=1\), put \(T=B/2\), \(R=A-B^2/4\), and use
\(v(p)=1\). Then \(\alpha=v(R)>0,\ \beta=v(T)>0\), with
infinity allowed. Completing the square gives the coefficient
of \(R^aT^{2b}\) at order \(q\) as
\[
 c_{ab,q}=\binom qa\binom{2(q-a)}{2b}
              \frac{(2(q-a-b)-1)!!}{q^{q-a-b}} .
\]
For \(a+b\le q\), its valuation is
\(\mathbf 1_{\{q-a>P\}}-\mathbf 1_{\{b>P\}}\).
Indeed the binomial times the double factorial equals
\((2(q-a))!/((2b)!\,2^{q-a-b}(q-a-b)!)\), and all other
factors are units. The two indicators cannot give a
negative valuation.

Set \(D_q=q-P\). The only candidates for minimum valuation
are \(1,R^{D_q},T^{p+1}\) with valuations
\(1,D_q\alpha,(p+1)\beta\), ignoring zero terms.
The constant term is uniquely smallest among the
coefficients of valuation one. In the valuation zero region
either \(a\ge D_q\) or \(b\ge P+1\), giving the two remaining
candidates. The two constant bounds \(1,(p+1)\beta\) must
agree. Otherwise, a tie with \(D_q\alpha\) is possible
for at most one \(q\), and no tie is possible if \(R=0\).
Let \(q_0\) be the least order. We obtain
\[
 \beta=\frac{1}{p+1},\qquad D_{q_0}\alpha\ge1 .
 \tag{1}
\]
In particular \(\alpha>2\beta\), since \(D_{q_0}<P+1\).
At an order \(q>q_0\), the terms with \(a\ge D_q\)
have valuation at least \(D_q\alpha\ge1+\alpha\)
when \(R\ne0\), so they lie strictly above \(1+2\beta\).
All other terms involving \(R\) also lie above
\(1+2\beta\): their coefficient has valuation one or
they have \(b\ge P+1\).

Define \(V=T^{p+1}/p\), a unit, and put
\(k_q=p c_{0,P+1,q}/c_{00,q}\).
At all \(q>q_0\), the terms of valuation at most \(1+2\beta\)
in \(G_q\), after dividing by \(c_{00,q}\), give
\[
 0=1+k_q V+T^2(a_q+k_q V b_q)+E_q,\qquad v(E_q)>2\beta ,
 \tag{2}
\]
where \(a_q=c_{01,q}/c_{00,q}\) and
\(b_q=c_{0,P+2,q}/c_{0,P+1,q}\). These are successive
coefficient ratios, since
\[
 \frac{c_{0,b+1,q}}{c_{0,b,q}}
         =\frac{2q(q-b)}{(2b+2)(2b+1)}.
\]
Hence \(a_q=q^2\) and
\(b_q\equiv q(2q-1)/6\bmod p\). There are no terms
omitted at the valuation of the correction in (2):
for \(a=0\) the next candidates have \(b\ge2\) in
the valuation one region or \(b\ge P+3\) in the
valuation zero region.

The constant residue in (2) gives \(\overline{k_q V}=-1\)
for all upper orders. Thus any difference of \(k_q\)'s
has valuation at least one, since these are rational
\(p\)-adic units. Subtract (2) at two upper orders, divide
by \(T^2\), and reduce. The coefficient
\(\overline{a_q-b_q}=(4q^2+q)/6\) must be the same for all
upper orders. There are at most two distinct upper orders
with the same quadratic value. For two such orders this is
equivalent to \(4(q_1+q_2)+1=0\bmod p\).

For the consecutive window we can also use the constant
residue in (2) in its explicit form. Wilson's congruence
gives
\[
 \overline{k_q}=-(2q)^{P+1}\,q(q-1)\cdots(q-P).
\]
Equality at \(q=n,n+1\) requires
\[
 1=\left(\frac{n+1}{n}\right)^{P+1}\frac{n+1}{n-P}
   =\chi\left(\frac{n+1}{n}\right)
                  \frac{2(n+1)^2}{n(2n+1)} ,
\]
with \(\chi\) the quadratic character. The two signs imply
\(p\mid (3n+2)(4n^2+5n+2)\). The second residue condition
requires \(8n+5=0\bmod p\). At \(n=-5/8\) the two factors
reduce to \(1/8\) and \(7/16\). The only possible prime is
therefore \(7\), which would require \(n=5\) from the interval
and fails \(8n+5=0\bmod7\). This proves the theorem.

## 3. Exact tests over finite fields

Multiplying \(G_q\) by \(q^q\) gives an integer polynomial,
with the coefficient of \(A^a B^{2b} C^{q-a-2b}\) equal to
\[
 \binom{q}{a,\,2b,\,q-a-2b}\,
      (2(q-a-b)-1)!!\;q^{a+b}.
\]
A common projective zero over \(\overline{\mathbb Q}\) would
reduce to a common projective zero over the algebraic closure
of every residue field. A modular test must therefore keep
the homogeneous directions at infinity.

There are two tests. At a prime \(p>2n+2\), check the chart
\(C=1\) by a Gröbner basis in \(A,B^2\); over an algebraic
closure replacing \(B\) by \(B^2\) loses no coefficient
directions. On \(C=0\), the point \(B=0\) is excluded by the
pure \(A\) term. Test \(C=0,\ B=1\) by a univariate gcd.
Unit ideals in both charts give an exact exclusion.

At an intermediate prime of Section 2, it suffices that the
reductions in \(A=1\) and in \(A=0,\ B=1\) have unit ideals.
The only direction outside those charts is \([0:0:1]\),
which cannot lift by Section 2. This combines exact residue
calculations with a common-lift obstruction. The pure
quadratic direction always survives in the residue system
since \(2(n-1)-1\ge p\).

### Lifting test at a non-pure residue point

In \(A=1\), put \(U=B^2\), \(V=C\), and
\(F_q(U,V)=q^q G_q(1,B,C)\). For a rational residue point,
choose integer representatives \(u,v\). The three rows
\[
 \left(\partial_U F_q(u,v),\ \partial_V F_q(u,v),\
                         F_q(u,v)/p\right)\bmod p
 \tag{3}
\]
give an exclusion if their determinant is nonzero.
The first two columns then contain an invertible two-row minor.
Any lift to the algebraic closure has both coordinate corrections
of valuation at least one. Otherwise their minimum valuation
would lie in \((0,1)\). For the two rows of the invertible minor,
some linear term would have that valuation, while the constant
term has valuation at least one and all nonlinear terms have
strictly higher valuation. Now divide the coordinate corrections
by \(p\), apply the three first-order equations and reduce. The
result is a linear system whose augmented matrix is (3), so it
has no solution. The same argument applies at a point over a
finite extension of the residue field using representatives
in its unramified extension of \(\mathbb Q_p\).

For example, at start \(71\) and prime \(109\), the
\(A=1\) residue ideal has reduced basis \(V-43,U+17\).
Thus the non-pure point \((U,V)=(92,43)\) survives the
modular test, and the \(A=0,B=1\) chart has unit ideal.
The rows in (3) are \((67,68,105),(34,32,90),(98,39,56)\),
with determinant \(106\bmod109\). This supplies a
certificate at \(109\) by combining Section 2 with the
first-order lift test at the other point.

### The constant-zero chart at fixed distance from the prime

For \(q=p-d\) in the intermediate interval calculation, put
\(D=q-(p-1)/2\). The restriction of \(q^qG_q\) modulo \(p\)
to \(A=0,\ B=1\) is a nonzero scalar times
\[
 J_d(C)=\sum_{t=0}^{\lfloor(d-1)/2\rfloor}
  \frac{(d-1)_{\underline{2t}}}{(2d-2)_{\underline{2t}}}
         \frac{(d/2)^t}{t!}\,C^{d-1-2t},
 \tag{2}
\]
where the underlined exponents are falling factorials and
the displayed formula is valid when \(p>2d\). To see this,
the surviving indices \(b\) of the \(B^{2b}\) terms run from
\(D\) to \(\lfloor q/2\rfloor\), because the double factorial
index is \(q-b\). Set \(b=D+t\) and divide by the coefficient
at \(t=0\). The multinomial factor gives
\((d-1)_{\underline{2t}}/(2d-2)_{\underline{2t}}\) after
reducing \(2D\equiv 1-2d\bmod p\). The other factors give
\[
 q^t\,\frac{(2(P-t)-1)!!}{(2P-1)!!}
             \equiv \frac{(d/2)^t}{t!}\pmod p.
\]
This proves (2).

In particular, take \(n=p-3\) for a prime \(p\ge11\).
The two upper orders have \(d=3,2\), so the polynomials are
\(C^2+1/4\) and \(C\), with no common root.
The residue direction \([0:0:1]\) cannot lift by
Section 2, including at \(p=23\). Thus for these starts any
algebraic common zero must have \(A\) a unit at every place
over \(p\) in projective integral coordinates. To prove
exclusion at such a start by this route, a unit-ideal test
in \(A=1\) is enough.

## 4. A verified prefix of odd windows

**Theorem.** The quadratic Gaussian window has no common
nonzero projective coefficient direction for every odd
start from 1 through 155 over any characteristic-zero field.

The [exact checker](verify_sfc_odd_gaussian_local_windows.py)
uses the two large-prime charts of Section 3 at \(65537\)
for all 64 odd starts through 127. The remaining 14 starts
are checked at the intermediate prime \(163\), where the
new local comparison applies even at start \(155\), when
the constant residues alone agree. At orders 156 and 157
both residues of \(k_q\) are 92, and the two residues of
\(a_q-b_q\) are 113 and 23. The commands, run from
the repository root with Python 3 and Singular, are

```sh
python3 research/factorial/image-rigidity/verify_sfc_odd_gaussian_local_windows.py
python3 research/factorial/image-rigidity/verify_sfc_odd_gaussian_local_windows.py \
    --first 129 --last 155 --prime 163 --intermediate
python3 research/factorial/image-rigidity/verify_sfc_odd_gaussian_local_windows.py \
    --first 71 --last 71 --prime 109 --intermediate --residue-point 92 43
```

Combining the odd prefix with the even-start gate proves
the large-height exclusion and the finite independent-height
decision for every fixed homogeneous bivariate trinomial
support and every prescribed start through 156. The cutoff
can depend on the prescribed start. The unresolved modular
condition is to find, for every remaining odd start, an
intermediate prime where all non-pure residues are absent
or have an obstruction to lifting.

An independent sufficient test uses the odd-pair resultant
of Section 1. The
[dyadic reduction](SFC_COLLINEAR_TERNARY_EVEN_START_AND_DYADIC_REDUCTION.md#odd-starts-and-the-remaining-valuation-pair)
forces \(v_2(C/B)=v_2(n)/2\) at every genuine common zero
of the triple. Hence it suffices to exclude a Newton polygon
slope \(-v_2(n)/2\) in the resultant in \(h=C/B\).
If that slope occurs, the resultant supplies an odd-pair zero
at the specified valuation. The even-order equation must
still vanish at the same coefficient direction.
