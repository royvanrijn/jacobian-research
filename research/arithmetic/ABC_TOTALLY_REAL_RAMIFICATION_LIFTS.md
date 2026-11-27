# Totally real fields for the abc radical

The abc radical can be transferred to the root discriminant of an
explicit totally real cyclic field, including when the normalized
radical on the lift is required to be bounded. An exact degree
formula controls this construction. On infinitely many inputs, any
abelian field of degree polynomial in the logarithmic height leaves
a height-sized radical on the lift. Totally real covers of the
three-marked line give a different constraint at the two real
endpoint marks, with Chebyshev covers attaining large ramification
at the third mark.

## Exact degree for prescribed odd-prime inertia

Take \(t=a/c\in(0,1)\) with \(a+b=c\) primitive, and put
\(h=\log c\), \(R=\operatorname{rad}(abc)\) and \(A=\log R\).
Write \(r_L(t)\) for the normalized logarithmic cusp radical over
a number field \(L\), as in the Legendre comparison. For abelian
\(L/\mathbb Q\), its odd-prime part is
\[
  r_L^{\rm odd}(t)=\sum_{\substack{p\mid R\\p>2}}
                         \frac{\log p}{e_p(L/\mathbb Q)}.
\]
Indeed every prime over \(p\mid R\) has positive cusp depth and
the sum of their residue degrees is \([L:\mathbb Q]/e_p\).
The contribution above \(2\) is at most \(\log 2\).

For each odd prime \(p\) in a finite set \(S\), prescribe a positive
inertia order \(e_p=p^{k_p}d_p\) with \(k_p\ge0\) and \(d_p\mid p-1\).
There is a totally real cyclic field of exact degree
\[
                   D=\operatorname{lcm}_{p\in S}e_p           \tag{TR1}
\]
with these inertia orders, unramified outside \(S\cup\{2\}\).
For empty \(S\) set \(D=1\). This degree is minimal even among
Galois fields realizing the orders, since every inertia order
divides the field degree.

Use the cyclotomic field of conductor
\(4\prod_{p\in S}p^{k_p+1}\). Map each odd-prime cyclotomic
factor onto the subgroup of order \(e_p\) of \(C_D\).
The images generate \(C_D\) by the definition of \(D\).
The image of the odd-prime components of \(-1\) has order
dividing two. Send the generator of the conductor-four factor
to that image's inverse. Then \(-1\) lies in the kernel, so the
fixed field is totally real. The inertia at each odd prime is
the image of its own factor. At \(2\), only the conductor-four
factor contributes to inertia, and the normalized logarithmic
discriminant contribution is at most \(\log 2\).

Conversely, by Kronecker–Weber the inertia at odd \(p\) in every
abelian extension of \(\mathbb Q\) is a quotient of some
\((\mathbb Z/p^n)^\times\), which is cyclic of order
\(p^{n-1}(p-1)\).
[Milne, Chapter I, Theorem 4.16](sources/milne-cft403.pdf).
Thus, under the divisibility budget \([L:\mathbb Q]\mid B\), its
order divides \(p^{v_p(B)}\gcd(p-1,B)\). Taking all these maximal
orders in (TR1) simultaneously shows the exact optimum
\[
 \min_{\substack{L/\mathbb Q\ {\rm abelian}\\[L:\mathbb Q]\mid B}}
      r_L^{\rm odd}(t)
   =\sum_{\substack{p\mid R\\p>2}}
       \frac{\log p}{p^{v_p(B)}\gcd(p-1,B)}.                  \tag{TR2}
\]
A totally real cyclic field attains the minimum. For a size
budget \([L:\mathbb Q]\le B\), take the minimum of (TR2) over
positive integers at most \(B\) in place of the divisibility
budget.

## Discriminant cost and bounded radical

For an abelian field with odd-prime inertia order \(e=p^kd\)
as above, its contribution to \(\log\operatorname{rd}(L)\) is
\[
  \left(1-\frac1e+k-\sum_{j=1}^{k}p^{-j}\right)\log p.       \tag{TR3}
\]
For the proof, pass locally to a finite unramified extension
containing the residue-field contribution from the other
cyclotomic factors. The normalized different is unchanged.
The totally ramified part is then the unique degree-\(e\)
subfield of the local \(p^{k+1}\)-cyclotomic extension. In the
full cyclotomic extension the uniformizer \(\zeta_{p^{k+1}}-1\)
gives the different via \(\Phi'_{p^{k+1}}(\zeta_{p^{k+1}})\).
Using \(\Phi_{p^{k+1}}(x)=(x^{p^{k+1}}-1)/(x^{p^k}-1)\) gives
the normalized different \(k+1-1/(p-1)\), in units of
\(\log p\). The extension over the degree-\(e\) subfield has
tame degree \((p-1)/d\). Subtracting its contribution by the
different tower formula gives (TR3).

First take \(S=\{p>2:p\mid R\}\) and \(e_p=p-1\). This gives a
tame odd-prime lift with
\[
     r_L(t)=O(\log(A+2)),\qquad
     \log\operatorname{rd}(L)=A+O(\log(A+2)).                 \tag{TR4}
\]
Split the sum of \(\log p/(p-1)\) at \(A\). Below \(A\) the
Chebyshev bound for \(\sum_{p\le x}\log p\) and partial summation
give \(O(\log(A+2))\). Above \(A\) use
\(\sum_{p\mid R}\log p=A\). The degree in (TR1) is at most \(R\).

For any fixed \(0<\eta<1\), there is also a totally real cyclic
lift satisfying
\[
 \begin{aligned}
 r_L(t)&=O_\eta(1),\\
 A&\le \log\operatorname{rd}(L)+r_L(t)
                     \le (1+\eta)A+O(1),\\
 \log[L:\mathbb Q]&\le(1+\eta)A .
 \end{aligned}                                               \tag{TR5}
\]
Take \(e_p=p(p-1)\) for \(p\le x\) in \(S\) and \(e_p=p-1\)
for other \(p\in S\), with \(x=\eta A/C\) where \(C\) is an
absolute Chebyshev constant such that
\(\sum_{p\le x}\log p\le Cx\). The extra different in (TR3)
sums to at most \(\sum_{p\le x}\log p\le\eta A\). Formula (TR1)
and \(\log(p-1)\le\log p\) give the degree bound. The radical
below \(x\) is bounded by the full sum of \(\log p/[p(p-1)]\).
Between \(x\) and \(A\), partial summation gives
\(O(1+\log(C/\eta))\), and above \(A\) the contribution is
bounded. Bounded \(A\) is absorbed in the constant. The lower
bound in (TR5) also follows directly at every \(p\mid R\) from
\(d_p^{\mathrm{loc}}\ge e_p-1\) for the local different.

In particular abc is equivalent to its pointwise height bound
with \(A\) replaced by \(\log\operatorname{rd}(L)\) for the
explicit tame fields in (TR4). The error in (TR4) is
\(O(\log(h+2))\) since \(A\le 3h\). One can use (TR5) in the
same equivalence with \(\eta\) chosen small in terms of the
target height tolerance. The base change keeps \(h(t)\)
unchanged. At odd primes the primitive Frey curve
\(y^2=x(x-a)(x+b)\) has multiplicative reduction precisely on
the cusp support, and base change preserves this assertion.
Its normalized odd-prime log conductor over \(L\) therefore
equals \(r_L^{\rm odd}(t)\).
[Primitive Frey reduction](ABC_FREY_REPAIR_PROJECTIVE_INERTIA.md#projective-representations-and-inertia-orders).

## An abelian degree obstruction on unbounded inputs

There are infinitely many primes \(p\equiv3\bmod4\) such that
every odd prime factor of \(p-1\) exceeds \(p^{1/10}\).
Apply the lower linear sieve to \(p-1\) for primes
\(X<p\le2X\) with \(p\equiv3\bmod4\), sifting only odd primes.
For odd squarefree \(d\), the two congruences give one reduced
class modulo \(4d\), with main term
\((\operatorname{Li}(2X)-\operatorname{Li}(X))/(2\varphi(d))\).
Thus the multiplicative density at an odd prime \(q\) is
\(1/(q-1)\). Its sifting product is \(\asymp1/\log z\) and it
satisfies the dimension-one interval product condition, since
its factors differ from \(1-1/q\) by \(O(1/q^2)\).
Bombieri–Vinogradov bounds the sum of absolute progression
errors for \(d<X^{1/3}\) by \(O(X/(\log X)^3)\).
Replacing the prime-count main term by the displayed
logarithmic integral over all these moduli gives the same
error by the prime number theorem and
\(\sum_{d\le X^{1/3}}\varphi(d)^{-1}\ll\log X\).
With \(z=X^{1/9}\) the lower-sieve parameter is \(3\), so its
positive main term is \(\gg X/(\log X)^2\). For large \(X\),
\(z>(2X)^{1/10}\), which proves the assertion.
[Rosser–Iwaniec lower sieve, Kumchev–Tolev §4.2, equations (4.21) and (4.24)](https://tigerweb.towson.edu/akumchev/a17.pdf).
[Bombieri–Vinogradov in the prime-count form](https://arxiv.org/abs/2510.10853v3).

For each such \(p\) take \(t=1/p\), so \(h(t)=\log p\).
In any abelian field of degree at most \(p^{1/10}\), the
inertia at \(p\) is tame and its order divides both \(p-1\)
and the degree. Since \(v_2(p-1)=1\), its order is at most
two. It follows that
\[
                  r_L(t)\ge\tfrac12 h(t).                    \tag{TR6}
\]
This eliminates a uniform abelian lift with \(o(h)\) residual
radical at any fixed polynomial degree bound in \(h\).
It applies regardless of the allowed discriminant cost or
auxiliary ramification. A wild remedy at \(p\) requires
degree at least \(p\), and by (TR3) its excess in the sum of
log root discriminant and local radical is at least
\((1-1/p)\log p\).

## The two endpoint marks for fixed real covers

Let a finite cover of the projective line over a number field
\(K\) be branched only at \(0,1,\infty\). At one real embedding,
suppose an unbranched fibre above \(0<t<1\) has all of its
geometric points real.
Then every geometric ramification index above \(0\) or \(1\)
is at most two. Indeed the number of real sheets of the
finite étale cover over the real interval \((0,1)\) is
constant, so fibres arbitrarily near either endpoint still
have all sheets real. A nonreal point over the endpoint
would supply nonreal sheets nearby. At a real point of
local index \(e\), use local coordinates with
\(t-t_0=\pm x^e\), absorbing the real analytic unit after
rescaling. At most two of the \(e\) nearby sheets are real,
so \(e\le2\).

After excluding the finitely many bad reduction primes of
the cover, specialization at a prime where \(t\) approaches
\(0\) or \(1\) therefore has tame ramification of order at
most two at every prime of the fibre field. This also holds
in a compositum of such fibre fields outside their union of
bad primes, because tame inertia is cyclic and all indices
divide two. In particular at the primes meeting \(0\) or \(1\)
outside that set at least half of the original normalized
logarithmic radical contribution remains on the lift.
[Darmon–Granville, Proposition 3.2](https://www.math.mcgill.ca/darmon/pub/Articles/Research/12.Granville/pub12.pdf).

The third mark can have arbitrarily large ramification in a
totally real fibre. Let \(K\) be totally real and suppose all
conjugates of \(t\in K\) are in \((0,1)\). Use the monic
Chebyshev polynomial \(P_N\) defined by
\(P_N(z+z^{-1})=z^N+z^{-N}\) and specialize at
\(P_N(x)=2-4t\). Under every real embedding the two roots
\(w,w^{-1}\) of \(w+w^{-1}=2-4t\) lie on the unit circle,
and so do the \(N\)-th roots of \(w\). Their traces are
distinct because \(w\ne w^{-1}\). Each root therefore defines
a totally real field \(L/K\) of relative degree at most \(N\).

At a prime \(\mathfrak p\nmid2N\) where \(m=-v_{\mathfrak p}(t)>0\),
the two local roots \(w,w^{-1}\) have valuations \(\pm m\)
and lie in \(K_{\mathfrak p}\). This follows by taking
the square root of \(1-4/(2-4t)^2\) by Hensel's lemma.
An \(N\)-th root \(z\) has inertia order \(N/\gcd(N,m)\).
The trace \(x=z+z^{-1}\) has negative valuation. Taking a
square root of \(1-4/x^2\) recovers \(z\) from \(x\) in the
same local field. Thus the inertia index is \(N/\gcd(N,m)\)
at every prime of every field factor. At the other primes
away from \(2N\), \(w\) is a unit and any inertia has order at
most two from the quadratic equation for \(w\). The extension
is unramified there unless \(t\) meets \(0\) or \(1\), since
\((2-4t)^2-4=-16t(1-t)\).

Write \(r_i\) for the weighted radical at cusp \(i\) over \(K\),
with \(i\in\{0,1,\infty\}\) and the weights of the Legendre
comparison. For every such \(L\),
\[
 r_\infty-\log(2N)-h(t)/N
       \le\delta(L/K)
       \le r_\infty+\tfrac12(r_0+r_1)+O(\log(2N)).             \tag{TR7}
\]
The denominator depth sum is at most \(h(t)\), so the lower
bound follows from the exact tame inertia formula there.
For the upper bound the tame endpoint contribution is at
most half its radical. The splitting field over \(K\) is
contained in \(K(w,\mu_N,\sqrt[N]{w})\), a Galois extension
of degree at most \(2N\varphi(N)\). The normalized different
bound from (PD2) gives \(O(\log(2N))\) at the omitted
primes, uniformly over the base. The residual radical at
the denominator primes away from \(2N\) is at most \(h(t)/N\).
Taking \(N\) to be a power of two of size \(h(t)+2\) leaves
bounded radical on the entire denominator support, since
only \(2\) was omitted. For \(t=a/c\) this transfers all but
bounded radical weight from \(c\), with degree \(O(\log c)\),
while the two endpoint marks still impose their lower bound.

## Consequence for the modular height route

Pasten's totally real comparison in Theorem 18.10 contains
the Shimura parametrization degree and a coefficient depending
on the base degree. The logarithmic conductor bound in
Theorem 18.13 assumes his modular-degree conjecture and
modularity, and still has a degree-dependent coefficient.
The unconditional estimate in Theorem 18.14 involves a
positive power of the conductor norm.
[Pasten, Theorems 18.10, 18.13 and 18.14](https://arxiv.org/abs/1705.09251).

The cyclic construction supplies totally real fields and
bounded normalized cusp radical with near-sharp discriminant,
at the cost of a potentially exponential field degree. The
abelian degree obstruction and the two-endpoint constraint
locate the loss in two ways of reducing that degree. A
norm-one trace construction now supplies totally real fields
of linear degree and solvable normal closure, with a single
possible auxiliary prime. The pointwise height inequality
still needs degree uniformity and a sharp discriminant
coefficient.
[Linear-degree norm-one lifts](ABC_TOTALLY_REAL_NORM_ONE_LIFTS.md#a-single-trace-gives-an-optimal-linear-degree-bound).
The mixed-signature Legendre isogeny fields already have the
low-degree discriminant property.
[Legendre isogeny discriminant comparison](ABC_PROJECTIVE_DIVISION_FIELD_DISCRIMINANTS.md#the-pointwise-modular-discriminant-obligation).
