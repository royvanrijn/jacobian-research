# Auxiliary unit tests for repaired target occurrence

Two support-prime layers, all three balanced radical caps, and the
original and repaired height conditions admit false candidates at the
full repaired-target exponent. We prove that these candidates can also
pass prescribed disjoint auxiliary prime-power tests. The count is
uniform in the auxiliary modulus, including moduli larger than the
height of the two integer parameters.

The underlying positive triples and the independent random-prime
formulation (SO1) are in
[True support occurrence for minimal repair](ABC_MINIMAL_REPAIR_SUPPORT_OCCURRENCE.md).
For fixed \(M\), write
\[
    H_M(n)=\prod_{p\mid n}p^{\max(M,v_p(n))}.
\]

## AU1. Uniform supply of false candidates

Fix an integer \(e\ge5\) with \(\gcd(e,6)=1\), an integer \(M\ge2e\),
and \(0<\sigma<4/3\). Let \(H\) tend to infinity. For every integer
\(1\le Q\le H^{4/3-\sigma}\) whose prime factors are \(1\bmod12\),
there are \(\gg_\sigma H^2\varphi(Q)/Q^2\) distinct pairs
\[
    3H\le u\le4H,\qquad H\le v\le2H
\]
such that the following conclusions hold. Put
\[
    \ell=u^2-v^2,\qquad m=2uv,\qquad n=u^2+v^2,\qquad
    (A,B,C)=(\ell^{2e},m^{2e},n^{2e}),\qquad
    \lambda=\frac{3}{2e}.
\]

- The three coordinates are positive and pairwise coprime.
  Each is comparable to \(H^{4e}\), and \(A+B<C\).
- Their radicals are \(\ell,m/2,n\), respectively. Each is comparable
  to \(C^{\lambda/3}\) and
  \(\operatorname{rad}(ABC)<C^\lambda\).
- The actual minimal repair of the target is \(H_M(C)=n^M\).
- \(ABC\) is coprime to \(Q\), and
  \[
       Q\operatorname{rad}(ABC)^2\mid A+B-C .
  \]

Uniformly in this range of \(Q\), the construction gives
\(\gg_{\sigma,\varepsilon} H^{2-\varepsilon}/Q\) distinct
target radicals for every \(\varepsilon>0\). At a repaired cutoff
\(Y=20^M H^{2M}\), the lower bound is \(Y^{1/M-o(1)}/Q\),
with implied constants at fixed \(e,M,\sigma\). In particular
subpolynomial \(Q\) preserves the full \(1/M\) exponent. The prime
powers in \(Q\) are arbitrary within the stated modulus range.

### Prescribing the auxiliary classes

In each \(\mathbb Z/q^j\mathbb Z\) with \(q^j\parallel Q\), choose
a primitive twelfth root \(\eta\) by lifting a root of order twelve
modulo \(q\). Such a root exists because \(q\equiv1\bmod12\), and
the derivative of \(T^{12}-1\) is a unit at the root. Put
\[
 I=\eta^3,\qquad s=\eta+\eta^{-1},\qquad
 t=\frac{(1+I)(1+s)}2 .
\]
The identities \(I^2=-1\), \(s^2=3\) and \(2\eta=s+I\) follow
by evaluating the twelfth cyclotomic polynomial. Direct substitution
gives
\[
 t^2=I(2+s)=\frac{1+\eta}{1-\eta},\qquad
 \frac{t^2-1}{t^2+1}=\eta,\qquad
 D(t)=0,\qquad D(T)=T^4-2T^3+2T^2+2T+1.
\]
Here \(t,t^2-1,t^2+1\) are units. For the quartic relation one may
multiply the two quadratics
\(T^2-(1+s)T+(2+s)\) and \(T^2-(1-s)T+(2-s)\).
All these identities hold at each prime power. Form \(t_Q\) and
\(I_Q\) by CRT. When \(Q=1\) the congruence and unit conditions are
vacuous.

Impose
\[
             u\equiv t_Q v\bmod Q,\qquad \gcd(v,Q)=1.       \tag{AU2}
\]
This makes \(u,v,u-v,u+v,n\) units modulo \(Q\). With
\(X=(\ell/n)^2\) in the residue ring, the identity
\(\ell^2+m^2=n^2\) gives \((m/n)^2=1-X\). On a component modulo
\(q^j\), \(X=\eta^2\) and \(1-X=X^{-1}\) because
\(X^2-X+1=0\). Thus
\(X^e+(1-X)^e=1\) since \(e\equiv\pm1\bmod6\).
This gives \(A+B-C\equiv0\bmod Q\).

### The parameter lattice and the small-prime sieve

The congruence lattice
\[
          L_Q=\{(u,v)\in\mathbb Z^2:u\equiv t_Qv\bmod Q\}
\]
has determinant \(Q\). On any nonzero vector the homogenization
\(v^4D(u/v)\) is a nonzero integer divisible by \(Q\). For
\(v\ne0\) it is nonzero by the rational root test for \(D\)
(the only possibilities are \(1,-1\), neither a root). For
\(v=0\) it is \(u^4\). Its absolute value is bounded by a constant
times the fourth power of the vector length, so a shortest
nonzero lattice vector \(b_1\) has \(|b_1|\gg Q^{1/4}\).

This shortest vector is primitive in \(L_Q\) and extends to a
lattice basis. Subtract an integer multiple of \(b_1\) from its
companion \(b_2\) so that the projection of \(b_2\) along \(b_1\)
has absolute length at most \(|b_1|/2\). The perpendicular length
is \(Q/|b_1|\), and \(|b_2|\ge|b_1|\) by shortestness. This gives
\[
           |b_1|\,|b_2|\asymp Q,\qquad |b_2|\ll Q^{3/4},
\]
with the sine of their angle bounded below. In the coordinates
\((x,y)\) of this basis, the image of the prescribed parameter
box lies in \(|x|\ll N_1,|y|\ll N_2\) for \(N_i=H/|b_i|\).
It contains a translated axis-parallel rectangle of side lengths
comparable to \(N_1,N_2\), by taking a small rectangle around the
inverse image of the center. Both \(N_i\gg H^{3\sigma/4}\), and
\(N_1N_2\asymp H^2/Q\).

Take \(P_0=\lfloor\sqrt{\log H}\rfloor\). Sieve \(u\) odd and
\(v\equiv2\bmod4\), and at every odd \(p\le P_0\) with \(p\nmid Q\)
exclude the conditions
\[
 p\mid\gcd(u,v),\qquad
 p^2\mid u,\ v,\ u-v,\ u+v,\ \text{or }u^2+v^2 .
\]
For each odd \(p\) the forbidden pairs modulo \(p^2\) have
proportion \(O(p^{-2})\). For the norm, outside common
divisibility, \(p^2\mid u^2+v^2\) forces \(v\) to be a unit and
\(u/v\) to be one of at most two roots of \(T^2+1\bmod p^2\).
Each prime admits, for instance, \(u\equiv1,v\equiv p\bmod p^2\).
Hence the product of allowed proportions, including the parity
condition, is bounded below absolutely. These residue conditions
have collective modulus
\[
          B=4\mathop{\prod_{2<p\le P_0}}_{p\nmid Q}p^2
           =H^{o(1)}
\]
and transfer to basis coordinates with unchanged proportions
since \(B\) is coprime to \(Q\).

One can also impose the unit condition of (AU2) with growing
\(\omega(Q)\). In basis coordinates \(v=\alpha x+\beta y\) with
\(\gcd(\alpha,\beta)=1\), since \(L_Q\) contains \((t_Q,1)\).
For a squarefree \(d\mid Q\), let \(d_1=\gcd(\alpha,d)\) and
\(d_2=d/d_1\). The restriction \(d\mid v\) first puts \(y\) in
the zero class modulo \(d_1\) and then \(x\) in one class modulo
\(d_2\) for each \(y\). Within fixed \(x,y\) classes modulo \(B\),
both progressions combine by CRT. On the inner rectangle, the
count for this fixed class and \(d\) is its area divided by
\(B^2d\), with error \(O(N_1+N_2+1)\). Inclusion-exclusion for
\(d\mid Q\), summed over the small-sieve classes, gives at least
\[
      c N_1N_2\frac{\varphi(Q)}Q
           -O\!\left(B^2 2^{\omega(Q)}(N_1+N_2+1)\right)
\]
pairs. All logarithms in this argument are natural. The bound
\((\omega(Q)+1)!\le Q\) shows \(2^{\omega(Q)}=H^{o(1)}\).
Ordering the auxiliary primes as \(q_1<q_2<\cdots\), one has
\(q_j\ge12j+1\). Thus
\[
 \frac{\varphi(Q)}Q
    \ge\prod_{j=1}^{\omega(Q)}\frac{12j}{12j+1}
    \gg (1+\omega(Q))^{-1/12}
    \gg(\log H)^{-1/12}.
\]
Consequently the error is lower order uniformly. Here and below,
lower order is relative to \(H^2\varphi(Q)/Q^2\).

### The large-prime deletions

First discard \(y=0\) at a cost \(O(N_1+1)\). For \(p>P_0\)
dividing both \(u\) and \(v\), one has \(p\nmid Q\) by (AU2), so
both basis coefficients are divisible by \(p\). As \(y\ne0\),
such primes are at most \(O(N_2)\). Using the containing rectangle
and summing \((N_1/p+1)(N_2/p+1)\) over these primes bounds this
loss by
\[
 O\!\left(\frac{H^2}{QP_0}+(N_1+N_2)\log(2H)+N_2\right).
\]
This is lower order. We may now use primitivity of each remaining
parameter vector.

We use the following direct observation. A collection of primitive
integer vectors of length at most \(R\), all in one projective
class modulo a modulus \(D_0\), has \(O(1+R^2/D_0)\) elements.
Here a projective class means a line generated by a vector
primitive modulo \(D_0\). Two vectors in one class have
determinant divisible by \(D_0\). Two distinct noncollinear
vectors therefore have angular distance at least \(D_0/R^2\)
within any sector of angular width at most \(\pi/2\). Collinear
primitive vectors differ by a sign, which can be separated into
different sectors.

If \(p>P_0\) and \(p^2\) divides one of \(u,v,u-v,u+v\),
the two parameters are in one projective class modulo \(p^2\),
and \(p\nmid Q\) by (AU2). Together with the class modulo \(Q\)
this is one class modulo \(Qp^2\). Since \(p\ll H^{1/2}\),
the observation bounds the loss by
\[
                      O\!\left(H^{1/2}+\frac{H^2}{QP_0}\right),
\]
again lower order.

It remains to exclude \(p^2\mid n\) for \(p>P_0\). For
\(P\le p<2P\), fixing \(p\) first gives at most two classes modulo
\(p^2\) from \(u^2+v^2=0\bmod p^2\) and primitivity. The same
observation and CRT give the dyadic upper bound
\[
                      O\!\left(P+\frac{H^2}{QP}\right).     \tag{AU3}
\]
A complementary count fixes the remaining Gaussian factor first.
Indeed \(\mathbb Z[i]\) is Euclidean for the norm by rounding each
coordinate of a complex quotient. If \(p\) were prime in this
ring, \(p\mid (u+iv)(u-iv)\) would give \(p\mid u,v\). Thus \(p\)
splits as two conjugate nonassociate primes of norm \(p\): a
proper factorization of \(p\) has factors of norm \(p\), and
an odd prime norm rules out a Gaussian prime being associated
to its conjugate. Only one of them can divide \(u+iv\) by
primitivity. Hence
\[
          w=u+iv=\pi^2 a,\qquad
          N(\pi)=p,\qquad 0<|a|\ll H/P .
\]
Here \(N\) denotes the Gaussian norm. There are \(O(H^2/P^2)\)
choices for \(a\) since \(P\ll H\).

Evaluate Gaussian integers modulo \(Q\) at \(i=I_Q\), writing
\(\bar w,\bar\pi,\bar a\) for the respective conjugates
before evaluation. All six evaluations are units because
\(N(w)\) is coprime to \(Q\). For a fixed \(a\), (AU2) specifies
\[
      \left(\frac{\pi}{\bar\pi}\right)^2
           =\frac{t_Q+I_Q}{t_Q-I_Q}\frac{\bar a}{a}
              \pmod Q .
\]
There are at most \(k=2^{\omega(Q)}\) unit square roots of the
right side. Each possible ratio \(r=\pi/\bar\pi\bmod Q\)
specifies one projective class for the integer coordinates
\((z_1,z_2)\) of \(\pi\). In fact it imposes
\((1-r)z_1+I_Q(1+r)z_2=0\bmod Q\), with at least one coefficient
a unit at each prime of \(Q\). These classes combine by CRT.
The coordinates are primitive since \(N(\pi)\) is prime, and
their length is \(O(P^{1/2})\). The angular bound now gives the
complementary estimate
\[
                     O\!\left(k\frac{H^2}{P^2}
                                 +k\frac{H^2}{QP}\right).  \tag{AU4}
\]
For \(Q=1\) the same bound follows with a single unconstrained
class. An arbitrary choice of the factorization for each \(w\)
suffices for the upper bound.

Use (AU3) for dyadic \(P\le H^{2/3}k^{1/3}\) and (AU4) above
that threshold. Summing the geometric series gives at most
\[
                     O\!\left(\frac{H^2}{QP_0}
                                 +H^{2/3+o(1)}\right)
\]
deleted parameters. In the high range, the second term of (AU4)
sums to \(O((H^2/Q)H^{-2/3}k^{2/3})\). These errors are lower order
at \(Q\le H^{4/3-\sigma}\). This completes the squarefree-value
count in the auxiliary classes. The unrestricted binary-form
sieve for a fixed form is classical. In particular Helfgott's
Theorem 5.2 applies to the homogeneous degree-six product
\(uv(u-v)(u+v)(u^2+v^2)\) in a full box. Its stated error for a
fixed form alone does not yield the uniform lower bound in
these growing-modulus classes. The two Gaussian counts above
supply that uniformity directly.
[Helfgott, *On the square-free sieve*, Theorem 5.2](https://www.impan.pl/shop/publication/transaction/download/product/83119)

### Radicals, repair and target multiplicity

The retained pairs have \(u\) odd, \(v\equiv2\bmod4\), and
\(\gcd(u,v)=1\), with \(u,v,u-v,u+v,n\) free of odd square
divisors. This gives pairwise coprime \(\ell,m,n\), with
\(\ell,n\) squarefree and \(\operatorname{rad}(m)=m/2\).
Since \(\ell,m<n\), it follows that
\(\operatorname{rad}(ABC)=\ell mn/2<n^3=C^\lambda\).
The squarefreeness of \(n\) gives \(H_M(C)=n^M\). With
\(U=\ell^2,\ V=m^2,\ W=n^2=U+V\), strict convexity gives
\(U^e+V^e<W^e\). The difference is divisible by \(UVW\) by
reducing separately modulo \(U,V,W\), using oddness of \(e\)
in the last reduction. Consequently every support prime
divides the error to depth at least two. Combine this with
the independent congruence modulo \(Q\) proved above.

Each \(n\) has at most \(4\tau(n)=H^{o(1)}\) representations
as a sum of two squares. This follows directly from Gaussian
unique factorization by allocating the prime factors above
each rational prime between \(u+iv\) and \(u-iv\).
Together with the parameter lower bound and \(n\le20H^2\)
this proves the claimed lower bound on distinct targets.

## Scope

The lower bound concerns prescribed candidate representatives.
It demonstrates the exponent retained on false representatives
if a relaxation uses the two support-prime layers and auxiliary
unit tests with total modulus \(H^{o(1)}\), even with all three
radical caps and actual repair heights. For polynomial \(Q\),
the displayed lower bound has the explicit \(1/Q\) cost,
up to subpowers. The original valuations are all equal to the
common exponent 2e except at 2 on the middle coordinate. Thus
this family has a single valuation degree carrying asymptotically
all of each radical at fixed budget. The construction leaves the
genuine occurrence count and the average test (SO1) unresolved. In
(SO1), the tested primes are sampled across an interval and
their product can be much larger than the range above.
