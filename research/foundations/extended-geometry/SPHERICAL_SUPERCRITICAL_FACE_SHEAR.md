# Sphere moments above the critical weighted root multiplicity

Write \(A=k[x,y,t]/(2xy+t^2-1)\) over a characteristic-zero field and use
the literal Laurent chart \(y=(1-t^2)/(2x)\). The normalized functional
is \(\Lambda(f)=\frac12\int_{-1}^1\operatorname{CT}_x f(x,t)\,dt\).
The [critical multiplicity theorem](SPHERICAL_CRITICAL_FACE_P_POWER_MOMENTS.md)
can certify nonvanishing from a face whose bounded-phase volume is
infinite. Above that multiplicity the leading face alone cannot decide
purity. There are regular absolute nonpowers with a common leading
angular polynomial, common strip and divergent volume, with opposite
answers to the same moment question.

## A purity criterion after an angular shear

Let \(P(x)\in xk[x]\). Introduce \(v=t-P(x)\) and
\(Y=y+(P(x)/x)t-P(x)^2/(2x)\), so \(2xY+v^2=1\).
The substitution preserves \(\Lambda\). Indeed the derivation
\(J=x\partial_t-t\partial_y\) acts as \(x\partial_t\) in the
Laurent chart, and this shear is its exponential with multiplier
\(P(x)/x\), up to the direction of the substitution. For any regular
\(g\), the angular constant term of \(P(x)g\) vanishes at both
interval endpoints: every negative angular term \(x^{-j}\) of a regular
phase is divisible by \((1-t^2)^j\). Thus
\(\Lambda((P(x)/x)Jg)=0\) by integration of a derivative. The same
argument applied to the higher derivatives in the finite exponential
proves invariance. This is the
[moment-preserving sphere shear](HOPF_LND_ORBIT_DEGREE_PROPAGATION_AND_POWER_DESCENT.md#3-a-completely-explicit-eventual-degree-bound).

**Theorem.** Suppose a fixed regular phase has Laurent expansion in
\((x,v)\) supported on nonpositive angles. Put
\(R(v)=\operatorname{CT}_x f(x,v+P(x))\). Then all its unmarked moments
are the normalized interval moments of \(R\). The following conditions
are equivalent: every positive moment vanishes, all sufficiently high
moments vanish, and \(R=0\). In this case every fixed regular marker
has zero marked moments at sufficiently high orders. More explicitly,
if the largest angle of a nonzero pure phase is \(-a<0\) in this chart,
and the largest angle of the transformed marker is \(L\), all orders
\(M>L/a\) vanish.

The nonpositive support makes the angular constant term of \(f^M\)
equal to \(R^M\). To test a nonzero \(R\) of degree \(k>0\), choose an
odd prime \(p\nmid k\) at which its coefficients are integral and its
leading coefficient is a unit. This is possible after one number-field
specialization preserving the leading coefficient if necessary. Take
any \(r>0\) with \(q=p^r\equiv1\pmod{k}\), and put \(M=(q-1)/k\).
Then the leading height is the even integer \(q-1\), and all lower
interval denominators have \(p\)-valuation below \(r\). Consequently
\[
 q\Lambda(f^M)\equiv\operatorname{lc}(R)^M\pmod p.
\]
This proves nonvanishing at infinitely many orders. A nonzero constant
\(R\) has every moment nonzero. When \(R=0\), the Laurent support of
the marked power lies at angles at most \(L-Ma\), proving the bound.

## The same supercritical face with two moment behaviours

Take positive integers \(d,B,n\), write \(b=dB\), and assume
\(n\ge2b\) and \(1\le\ell\le B\). Take \(c\ne0\) and
\(U\in k[Z]\) with \(\deg U\le B-\ell\) and \(U(0)U(c)\ne0\).
Use \(w=x^d-ct\) and let
\[
 f_0=x^{-b} t^{B-\ell}w^{n+\ell-2b}(w^2-c^2)^b U(x^d/t),
 \qquad v=t-x^d/c.
\]
The factor \(t^{B-\ell} U(x^d/t)\) is interpreted as a polynomial.

**Corollary.** The phase \(f_0\) is regular and has top angular
polynomial
\[
   Q(z)=z^{-b}(z^d-c)^{n+\ell}U(z^d)
       \quad\text{at weight } j+i/d=n.
\]
The full literal weights lie between \(n-2b\) and \(n\). For any
polynomial \(R(v)\) with support in degrees \([n-2b,n-1]\), the
phase \(f_0+R(v)\) has the same top face and strip. The squared
volume of a bounded phase region diverges at each root \(z^d=c\),
including for \(R=0\). Its moments are exactly the interval moments
of \(R\). The phases are absolute nonpowers if
\(\gcd(b,n+\ell)=1\) or \(U\) has a simple nonzero root other than
\(c\). More generally, within regular phases in the given face and
strip whose Laurent support is nonpositive after this shear, purity
is exactly the linear condition that the angle-zero polynomial vanish.

For regularity, put \(u=x^d\). At \(u^j\), the coefficient of
\((w^2-c^2)^b=(u^2-2ctu+c^2(t^2-1))^b\) is divisible by
\((t^2-1)^{b-j}\) when \(j\le b\). In the full numerator the
coefficient of \(u^i\) for \(i<B\) is therefore divisible by
\((t^2-1)^{b-i}\), enough for the negative angle \(di-b\) since
\(b-i\ge b-di\). The leading part uses \(w^2\) in all factors of
\((w^2-c^2)^b\), and the deficits are the even integers at most
\(2b\). Each monomial of \(v^k\) has weight \(k\).

After the shear, \(w=-cv\), and the factor
\((v+x^d/c)^{B-\ell}U(x^d/(v+x^d/c))\) is a polynomial in \(x^d,v\)
of degree at most \(B-\ell\) in \(x^d\). After division by \(x^b\),
all angles are at most \(-d\ell\). The theorem gives the moment
identity and purity criterion. On a disk in \(v\), with \(|x|\)
large, \(f_0+R(v)\) is uniformly bounded and the volume density is
a nonzero constant times \(|x|^{-2}d^2x\,d^2v\), which diverges.
For a specified root \(z^d=c\) one may restrict \(x\) to an angular
sector in this region and use the corresponding branch of \(t^{1/d}\).
The leading angular function of a proper power would itself be a
proper power over \(\overline k(z)\). Its orders at \(z=0\) and at
the critical roots are \(-b\) and \(n+\ell\), and any simple
nonzero root of \(U(z^d)\) other than these gives order one.

These two types have the same isolated top-pole residue, namely
zero at the orders in the critical theorem. In fact with \(q=p^r\)
and \(nM=sq-1\), one reduces modulo \(Z^q\) in the coefficient
\([Z^{BM}](Z-c)^{(n+\ell)M}U(Z)^M\). The factor
\((Z-c)^{sq}\) contributes only its constant term because \(BM<q\).
The remaining polynomial \((Z-c)^{\ell M-1}U(Z)^M\) has degree
at most \(BM-1\). The zero residue is thus compatible with a
nonzero moment from a lower denominator. No uniform nonvanishing
theorem on arbitrary thickening terms can include these faces.

For example take \(d=2,B=2,n=37,\ell=1,c=1,U(Z)=Z-2\).
The order-two moments of \(f_0\) and \(f_0+v^{33}\) are \(0\)
and \(1/67\), respectively. The top denominator at \(p=5,r=2,s=3\)
is \(75\) and its normalized residues vanish for both phases.
[Exact rational replay](verify_spherical_supercritical_face_shear.py)

## A complete test for a family with separated root clusters

Fix \(d>K\ge1\), \(b\ge1\), \(c\ne0\), distinct \(\alpha_i\), and
nonnegative integers \(a_i,b_i-a_i\) with
\(\sum_i b_i=2b\), \(\sum_i a_i=b\), and \(b_i>0\).
Set \(w_i=x^d-ct-\alpha_i x^K\) and consider
\[
 f=x^{-b}\prod_i (w_i-c)^{a_i}(w_i+c)^{b_i-a_i}.
 \tag{1}
\]
The original weight \(1/d\) sees the angular face
\(z^{-b}(z^d-c)^{2b}\), of height \(2b-b/d\). It therefore
sees a supercritical multiplicity at every nonzero angular root.

**Proposition.** Put \(T=2b-b/K\). If \(\max_i b_i>T\),
the phase is pure and every fixed marker vanishes eventually.
If \(\max_i b_i<T\), the bounded-phase squared volume is finite
and infinitely many moments are nonzero. At equality the volume
diverges and the exact moment sequence, for any index \(i\) with
\(b_i=T\), is
\[
 \Lambda(f^M)=C_i^M(-1)^{a_iM}(2c)^{b_iM}
     \frac{(a_iM)!((b_i-a_i)M)!}{(b_iM+1)!},
 \qquad C_i=\prod_{j\ne i}(\alpha_i-\alpha_j)^{b_j}.
 \tag{2}
\]
In particular every moment at equality is nonzero. If one
nonzero exponent in (1) is one, the phase is an absolute nonpower.

To see regularity, at \(x=0\) exactly \(b\) factors vanish at
each of \(t=1\) and \(t=-1\), counting multiplicities.
The coefficient of \(x^s\) in the numerator, for \(s<b\),
can use positive \(x\)-degrees in at most \(s\) factors.
It is thus divisible by \((1-t^2)^{b-s}\). The last assertion
about powers follows from the multiplicity-one factor in the
localization \(k[x,x^{-1},t]\), since the displayed degree-one
polynomials in \(t\) are pairwise nonassociate there.

For a specified cluster \(i\) shear to
\(v=t-(x^d-\alpha_i x^K)/c\).
The cluster contributes \((-cv-c)^{a_i}(-cv+c)^{b_i-a_i}\).
The angular degree of the numerator from the other factors is
\(K(2b-b_i)\), with leading coefficient \(C_i\). If \(b_i>T\)
all angles are strictly negative. At equality their maximum
is zero, giving (2) by the interval beta integral. In either
case, \(f\) stays bounded on a disk in \(v\) as \(|x|\to\infty\)
and the same logarithmic volume density proves divergence.
This gives a divergent end for each cluster satisfying \(b_i\ge T\).

If \(b_i<T\) for every \(i\), shear instead to a centre
\((x^d-\alpha x^K)/c\) with \(\alpha\ne\alpha_i\) for all
\(i\). At weight \(1/K\) in this chart the height is \(T\)
and the leading angular polynomial is
\[
          z^{-b}\prod_i\bigl((\alpha-\alpha_i)z^K-c\bigr)^{b_i}.
\]
Its extremes are exactly those of the full support, and its
nonzero roots have multiplicities \(b_i<T\). Apply the
[bilateral weighted criterion](SPHERICAL_WEIGHTED_FACE_PROPER_VOLUME_OBSTRUCTION.md#the-bilateral-weighted-criterion)
in this chart. The shear preserves both the moments and the
volume form, since \(dP(x)\wedge dx=0\).

When \(d\mid b\), write \(B=b/d\) and \(n=2b-B\).
For \(q=p^r\equiv1\bmod n\) at an odd good prime and
\(M=(q-1)/n\), the unique highest interval pole at \(q\)
has zero normalized residue for every phase (1). Indeed its
coefficient is \([Z^{BM}](Z-c)^{2bM}\), and
\(2bM=q+BM-1\). Below degree \(q\) the factor
\((Z-c)^q\) contributes only its constant coefficient.
This pole gives the same zero for all three cases in the
proposition.

**Corollary (Secondary depths for higher separation).** At equality
let \(b_i=T\), \(0<a_i<T\), and \(d\mid b\). At an odd good prime
\(p\nmid n(n-T)\), with \(c\) and the differences in \(C_i\) units,
the valuation at the same order \(M=(p^r-1)/n\) is
\[
 v_p(\Lambda(f^M))=
 -\#\{0\le j<r:[a_i p^j]_n+[(T-a_i)p^j]_n\ge n\}.
\]
It equals \(-r/2\) whenever \(-1\in\langle p\bmod n\rangle\).

Since \(K<d\), one has \(T<n\). The factor \(TM+1\) in the beta
denominator is a unit because it reduces to \((n-T)/n\).
The [factorial floor argument for a secondary pole](SPHERICAL_SUPERCRITICAL_SECOND_POLE.md#the-vanished-top-residue-and-the-exact-lower-depth)
applies with \(T,a_i,T-a_i\) in place of \(b,a,b-a\).
For clarity, for \(0<A<n\) and \(1\le j\le r\), subtracting
\(A/(np^j)\) from \(A p^{r-j}/n\) leaves the floor unchanged:
the nonzero residue \([A p^{r-j}]_n\), after multiplication by
\(p^j\), is a positive integer congruent to \(A\) modulo \(n\).
Legendre's factorial formula then yields the displayed count.
Pairing an orbit element with its negative gives exactly one
carry, since the residue sum cannot equal \(n\) when \(T<n\).

For an unbounded family with nonlinear separation, take
\(d=14\), \(K=2\), \(b=14\cdot3^{h-3}\), \(h\ge3\), and two
clusters of sizes \(3b/2\) and \(b/2\). Assign \(a_1=b-1\)
and \(a_2=1\) and take \(c=1,\alpha_1=1,\alpha_2=2\).
The small cluster supplies a simple factor, so the phase is an
absolute nonpower. The original face has \(n=3^h\) and the large
cluster gives divergent volume. The primes \(p\equiv2\bmod 3\)
have density one half and, after excluding \(2\), give exact
secondary pole depth \(r/2\) at infinitely many original orders.

## A second pole after a lower change of centre

**Example.** Let \(n>15\), \(n\equiv1\pmod6\), and \(k=n-2\).
For nonzero \(a,c\) set \(w=x^2-ct-cax\) and consider
\[
            f=x^{-2}w^{n-3}(w^2-c^2)^2+w^k.
\]
The top face at weight \(1/2\) is \(z^{-2}(z^2-c)^{n+1}\).
In the chart \(t-x^2/c\), which subtracts only the centre
prescribed by the face, this phase still has both positive and
negative angles. The largest angle is \(n-1\), from the term
with every factor \(w\) replaced by \(-cax\) in \(x^{-2}w^{n+1}\).
The angle \(-2\) has a nonzero polynomial coefficient. Regularity
follows from the first example by the additional shear
\(t\mapsto t+ax\).

The full literal weights lie between \((n-5)/2\) and \(n\).
For the first summand the minimum arises by taking the
constant from \((w^2-c^2)^2\) and the term in \(x\) from every
remaining \(w\), and the second summand has weights at least
\(k/2\). Along \(w\) bounded with \(|x|\) large, the phase is
bounded and the squared volume diverges. With the shifted
centre \(P(x)=x^2/c-ax\), the first summand has angle \(-2\)
and the angle-zero polynomial is \((-cv)^k\). Hence
\[
 \Lambda(f^M)=
   \begin{cases}
      (-c)^{kM}/(kM+1)&\text{if \(kM\) is even},\\
      0&\text{if \(kM\) is odd}.
   \end{cases}
\]
In particular put \(q=p^r\) at odd good primes, with
\(3q\equiv1\pmod{\operatorname{lcm}(2n,k)}\).
At \(M=(3q-1)/n\) the strip has width less than \(2q\).
The isolated top residue at \(3q\) vanishes by the same degree
test as above: it is a coefficient at degree \(M<q\) of
\((Z-c)^{3q}(Z-c)^{M-1}\). At the separate order
\((3q-1)/k\), the lower-weight polynomial supplies exactly
the denominator \(3q\) with unit numerator, giving a nonzero
normalized residue. The chosen prime powers exist at every good
prime \(p\equiv3^{-1}\pmod{\operatorname{lcm}(2n,k)}\) since
the modulus is prime to 3.
