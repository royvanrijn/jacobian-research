# Obstructions to direct models for simultaneous norm correlations

The simultaneous Gaussian--Eisenstein weight has an exact character-divisor
expansion. Truncating its divisors at any fixed power of the component size,
however, leaves a full-order error in both the first and second absolute
moments. The truncation is also not a pointwise minorant, even inside the
required residue class. These facts identify a cancellation problem in a
natural route to the pointwise binary correlation; they do not rule out a
refined divisor method.

A separate algebraic obstruction concerns the three quadrics defining the
simultaneous norm equations: no invertible linear changes of variables and
linear recombination of equations identify them with the standard binary
Gram system in a four-dimensional quadratic space. This excludes that direct
identification with a degree-two Siegel theta representation, while leaving
other theta and spectral constructions open.

The unrestricted pointwise lower bound and fixed-prime marginals in the
[pointwise sieve note](SIMULTANEOUS_NORM_POINTWISE_SIEVE.md) remain open.

## The exact expansion

Let \(h(n)=r_G(n)r_E(n)/24\), where \(r_G,r_E\) count all ordered integer
representations by \(x^2+y^2\) and \(u^2-uv+v^2\), respectively. Set

\[
 h_6(n)=h(n)\mathbf1_{(n,6)=1},\qquad k=h_6*\mu.
\]

Here \(*\) is Dirichlet convolution. Thus \(h_6=1*k\). A supported value
of \(h_6\) is \(1\pmod{12}\). The Euler factors established in the
[case study](SIMULTANEOUS_GAUSSIAN_EISENSTEIN_SUMS.md#what-the-global-correlation-actually-is)
give, for \(\Re s>1\),

\[
 \sum_{n\ge1}\frac{k(n)}{n^s}
 =(1-2^{-2s})(1-3^{-2s})
 \frac{L(s,\chi_{-4})L(s,\chi_{-3})L(s,\chi_{12})}
      {L(2s,\chi_{12})}.
 \tag{1}
\]

The characters are the primitive quadratic characters of the indicated
discriminants. For \(p\ge5\), the local coefficients are

\[
 k(p^e)=
 \begin{cases}
 2e+1,&p\equiv1\ (12),\\
 e+1,&p\equiv5,7\ (12),\quad e\text{ even},\\
 -e,&p\equiv5,7\ (12),\quad e\text{ odd},\\
 (-1)^e,&p\equiv11\ (12).
 \end{cases}
 \tag{2}
\]

Indeed \(k(p^e)=h_6(p^e)-h_6(p^{e-1})\). At \(p=2,3\),
\(k(p)=-1\) and \(k(p^e)=0\) for \(e\ge2\).
In particular \(k(p)=3\) on \(p\equiv1\pmod{12}\) and \(-1\) on
the other prime classes. Define

\[
 h_D(n)=\sum_{\substack{d\mid n\\d\le D}}k(d).
 \tag{3}
\]

## Full-order remainder at every fixed power cutoff

**Theorem.** For each fixed \(0<\gamma<1\), put \(D=X^\gamma\).
As \(X\to\infty\),

\[
 \sum_{\substack{X/2<n\le X\\n\equiv1\ (12)}}
       |h_6(n)-h_D(n)|\gg_\gamma X,
 \qquad
 \sum_{\substack{X/2<n\le X\\n\equiv1\ (12)}}
       |h_6(n)-h_D(n)|^2\gg_\gamma X(\log X)^3.
 \tag{4}
\]

These are the natural first- and second-moment scales of \(h_6\).
Projecting \(h_D\) onto the residue class \(1\pmod{12}\) does not change
either lower bound.

To prove this, use the one-variable Frobenian mean theorem,
[Loughran--Matthiesen, arXiv:1904.12845v2, Lemma 2.8](https://arxiv.org/html/1904.12845v2#S2).
The functions \(h_6\) and \(h_6^2\) have prime means \(1\) and \(4\):
their prime values are \(4\) and \(16\), respectively, on a prime set of
density \(1/4\), and zero elsewhere. Both are nonnegative multiplicative
functions with polynomial prime-power growth and divisor bounds. Thus,
for positive constants \(c_1,c_2\), partial summation gives

\[
 \sum_{m\le y}\frac{h_6(m)}m\sim c_1\log y,
 \qquad
 \sum_{m\le y}\frac{h_6(m)^2}m\sim\frac{c_2}{4}(\log y)^4.
 \tag{5}
\]

Choose \(0<\varepsilon<\min(\gamma,1-\gamma,1/2)\). Consider
\(m\le X^\varepsilon\) and primes \(p\equiv1\pmod{12}\) in

\[
 X/(2m)<p\le X/m.
\]

For sufficiently large \(X\), these satisfy \(p>D\), \(p>\sqrt X\),
\(m<D\) and \(p>m\). Every divisor of \(pm\) at most \(D\) is
therefore a divisor of \(m\), and all divisors of \(m\) are included.
Consequently

\[
 h_D(pm)=h_6(m),\qquad h_6(pm)=4h_6(m).
 \tag{6}
\]

When \(h_6(m)>0\), the integer \(pm\) is \(1\pmod{12}\).
Distinct pairs \((m,p)\) give distinct integers \(pm\), since \(p>\sqrt X\)
is their unique prime factor above \(\sqrt X\).
The fixed-progression prime number theorem, equivalently
[Loughran--Matthiesen, Lemma 2.4(1)](https://arxiv.org/html/1904.12845v2#S2),
gives \(\gg X/(m\log X)\) such primes uniformly for \(m\le X^\varepsilon\).
This uniformity uses only that all interval endpoints tend to infinity at
least as fast as \(X^{1-\varepsilon}/2\); the modulus is the fixed integer 12.
Each selected integer contributes \(3h_6(m)\), or \(9h_6(m)^2\), to the
corresponding absolute-error moment. Summing and using (5) proves (4).

The same argument gives, for each fixed \(s\ge1\),

\[
 \sum_{\substack{X/2<n\le X\\n\equiv1\ (12)}}
 |h_6(n)-h_D(n)|^s
 \gg_{\gamma,s}X(\log X)^{4^{s-1}-1},
 \tag{7}
\]

because \(h_6^s\) has prime mean \(4^{s-1}\). Thus the missing large-prime
cofactor terms carry full-order mean mass, not just rare large values.

## The truncation is not a pointwise minorant

A concrete square-root example is

\[
 n=5\cdot7\cdot11\cdot17\cdot29=189805\equiv1\pmod{12},
 \qquad\lfloor\sqrt n\rfloor=435.
\]

All five prime factors are forbidden, so \(h_6(n)=0\). Among its divisors
at most 435 are the unit, all five single primes, nine prime pairs
(only \(17\cdot29=493\) is excluded), and one triple
\(5\cdot7\cdot11=385\). No divisor with more factors occurs. By (2),

\[
 h_{\sqrt n}(n)=1-5+9-1=4>0.
 \tag{8}
\]

In fact, for every fixed \(0<\gamma<1\), infinitely many squarefree
\(n\equiv1\pmod{12}\) satisfy \(h_6(n)=0<h_{n^\gamma}(n)\).
Choose all primes below in a common interval \([Y,2Y]\); the fixed-modulus
prime number theorem supplies the required distinct primes as \(Y\to\infty\).
For \(\gamma<1/2\), take two primes \(5\pmod{12}\); no nonunit divisor
is included, giving sum 1. For \(\gamma=1/2\), take five primes in classes
\(5,5,5,7,11\); exactly the subsets of size at most two are included,
giving \(1-5+10=6\). For \(\gamma>1/2\), choose an even integer \(\ell\)
with \(1<(1-\gamma)\ell<2\), and take \(\ell\) primes \(5\pmod{12}\).
Their product is \(1\pmod{12}\), and precisely subsets of size at most
\(\ell-2\) are included. Their sum is
\(\sum_{j=0}^{\ell-2}(-1)^j\binom\ell j=\ell-1>0\).
Such an even \(\ell\) exists because its permitted interval has length
greater than two. Strict gaps between powers of \(Y\) justify the stated
subset cutoffs despite the fixed factors \(2^\ell\).

## The absolute CRT rounding budget

Expanding a truncated binary convolution gives

\[
 \sum_{a+b=N}h_D(a)h_D(b)=\sum_{d,e\le D}k(d)k(e)
 \#\{a,b>0:a+b=N,\ d\mid a,\ e\mid b\}.
\]

The inner count vanishes unless \((d,e)\mid N\). Otherwise it equals
\(N(d,e)/(de)+O(1)\). Replacing each rounding error by its absolute
constant bound incurs the budget

\[
 \sum_{\substack{d,e\le D\\(d,e)\mid N}}|k(d)k(e)|
 \asymp D^2\log D,
 \tag{9}
\]

uniformly in the positive integer \(N\). This is a statement about that
absolute-error majorant, not a lower bound for the actual signed error.

Here is a justification of both bounds. Put \(G=|k|\). By (2), it is
nonnegative Frobenian with prime mean \(3/2\), \(G(n)\le\tau(n)^2\),
and \(G(pn)\le3G(n)\) at every prime. Lemma 2.8 gives
\(\sum_{d\le D}G(d)\asymp D(\log D)^{1/2}\), proving the upper bound.
For the lower bound, restrict to coprime pairs, which are included for every
\(N\). After excluding all primes at most a fixed \(Y\), the remaining
\(G\)-mass is asymptotic to \(c(Y)D(\log D)^{1/2}\), where
\(c(Y)>0\) and \(c(Y)\asymp(\log Y)^{-3/2}\). This follows from the
same mean theorem and its finite Euler factors; the last estimate uses the
prime mean and the convergent higher Euler-factor errors.
Pairs in this population sharing a prime greater than \(Y\) have total
weight at most

\[
 \sum_{Y<p\le D}\left(\sum_{\substack{d\le D\\p\mid d}}G(d)\right)^2
 \le9\sum_{Y<p\le D}\left(\sum_{u\le D/p}G(u)\right)^2
 \ll D^2\log D\sum_{p>Y}p^{-2}
 \ll\frac{D^2\log D}{Y}.
\]

Choose one sufficiently large fixed \(Y\) so this loss is smaller than
the positive pair mass; this is possible since
\(Yc(Y)^2\asymp Y/(\log Y)^3\to\infty\). Coprime pairs retain a
positive multiple of \(D^2\log D\), proving (9).

At \(D\ge\sqrt N\), this absolute rounding budget is not \(o(N)\).
At a smaller fixed-power cutoff, the budget can be small, but (4) shows that
the omitted divisors still have full-order absolute moments. Small-divisor
truncation followed only by this rounding bound and a global
Cauchy--Schwarz remainder estimate therefore does not give the required
\(o(N)\) binary correlation error. Additional signed cancellation is needed.

## A limited Hecke-recursion diagnostic

At an allowed prime, \(h(p^e)=(e+1)^2\) is not an ordinary unramified
\(\mathrm{GL}_2\) Hecke eigen-coefficient sequence. An order-two recurrence
with eigenvalue \(h(p)=4\) and \(h(1)=1\) would force its determinant
parameter to be \(16-h(p^2)=7\). It would then predict
\(h(p^3)=4\cdot9-7\cdot4=8\), whereas the actual value is 16.
Rescaling by \(p^{te}\) leaves the mismatch. This excludes only that direct
eigen-sequence model; it does not exclude a non-eigenfunction construction,
higher-rank automorphic methods or a different spectral model for the binary
correlation.

## The simultaneous norm net is not a standard binary Gram system

Write
\[
 G(x)=x_1^2+x_2^2,\qquad E(u)=u_1^2-u_1u_2+u_2^2.
\]
For binary vectors \(x,y,u,v\), the simultaneous norm equations are
\[
 G(x)+G(y)=N,\qquad G(x)=E(u),\qquad G(y)=E(v).
\]
Their homogeneous quadratic parts span the net
\[
 F_1=G(x)+G(y),\qquad F_2=G(x)-E(u),\qquad
 F_3=G(y)-E(v)
 \tag{10}
\]
in eight variables. Here a net is the three-dimensional linear space of
quadratic forms, with projective parameter space \(\mathbf P^2\).

**Proposition.** Over \(\mathbb C\), no invertible linear change of the
eight variables, together with an invertible linear recombination of the
three forms in (10), converts this net into the standard binary Gram net
\[
 q(s),\qquad 2\langle s,t\rangle_Q,\qquad q(t),
 \qquad s,t\in\mathbb C^4,
 \tag{11}
\]
for any nondegenerate symmetric matrix \(Q\), where
\(q(s)=s^{\mathsf T}Qs\) and
\(\langle s,t\rangle_Q=s^{\mathsf T}Qt\).
Consequently such a conversion is also impossible over \(\mathbb Q\).

**Proof.** The symmetric matrix of \(E\) is
\[
 E_2=\begin{pmatrix}1&-1/2\\-1/2&1\end{pmatrix},
 \qquad \det E_2=3/4.
\]
In the variable order \((x,y,u,v)\), the symmetric matrix of
\(AF_1+BF_2+CF_3\) is
\[
 M(A,B,C)=\operatorname{diag}
   ((A+B)I_2,(A+C)I_2,-BE_2,-CE_2).
\]
Hence its determinant is exactly
\[
 \det M(A,B,C)=\frac9{16}
       \{BC(A+B)(A+C)\}^2.
 \tag{12}
\]
Its projective degeneracy curve therefore has four distinct line components,
each with multiplicity two.

The corresponding matrix for the Gram net (11) is
\[
 M_Q(A,B,C)=
 \begin{pmatrix}AQ&BQ\\BQ&CQ\end{pmatrix}.
\]
The off-diagonal blocks are \(BQ\) because the middle form in (11) is
\(2s^{\mathsf T}Qt\). The determinant identity for this block matrix gives
\[
 \det M_Q(A,B,C)=(\det Q)^2(AC-B^2)^4.
 \tag{13}
\]
For example, this follows by writing the matrix as the tensor product of
\(\begin{pmatrix}A&B\\B&C\end{pmatrix}\) with \(Q\).
The reduced degeneracy curve is the smooth conic \(AC-B^2=0\): its gradient
\((C,-2B,A)\) has no projective zero. Over \(\mathbb C\), this conic is
irreducible, since a degree-two product of linear forms has a projective
singular point, including when the two linear factors coincide. Choosing the unscaled middle Gram form instead merely
rescales the parameter \(B\), leaving this conclusion unchanged.

A variable change \(z=Tz'\), with \(T\in\mathrm{GL}_8(\mathbb C)\),
replaces every matrix of the net by \(T^{\mathsf T}MT\), multiplying its
determinant by the nonzero constant \((\det T)^2\). An invertible linear
recombination of its three forms acts by a projective linear change of
\((A:B:C)\). Thus these two operations preserve the number and degrees of
the irreducible components of the reduced projective degeneracy curve.
Four lines cannot become one irreducible conic, proving the proposition.

This is an obstruction only to the stated linear identification with a
standard rank-four binary Gram system. It does not exclude nonlinear
substitutions, birational parametrizations, auxiliary variables, trace
integrations, other automorphic realizations, or Siegel and higher-rank theta
methods in general. The count of three quadratic equations in eight variables
alone therefore does not supply that particular Gram representation.

## A bounded norm value prevents strong approximation off infinity

The fixed-target norm variety also rules out a stronger geometric shortcut.
This obstruction concerns density in the finite adeles; it does not imply
failure of integral existence or of rational weak approximation.

**Proposition SimNormBoundedValue.** Let \(X\) be an affine \(\mathbb Q\)-variety with an affine
integral model, and let \(A\) be an integral regular function. Suppose
\(A(X(\mathbb R))\cap\mathbb Z\) is finite, \(X(\mathbb Z_\ell)\ne\varnothing\)
for every finite prime \(\ell\), and \(A(X(\mathbb Z_p))\) contains a
nonempty \(p\)-adic open ball for some prime \(p\). Then \(X(\mathbb Q)\)
is not dense in \(X(\mathbb A_f)\): strong approximation off infinity fails.

**Proof.** The set \(S=A(X(\mathbb Z))\) is finite. Choose a congruence
ball \(B\subset A(X(\mathbb Z_p))\) disjoint from \(S\). The set
\[
 \{x\in X(\mathbb Z_p):A(x)\in B\}
 \times\prod_{\ell\ne p}X(\mathbb Z_\ell)
 \tag{14}
\]
is nonempty and open in the finite adeles. A rational point in it has every
affine coordinate integral at every finite prime, hence is an integral
point. Its \(A\)-value would belong to both \(S\) and \(B\), which is
impossible. \(\square\)

Apply this to every positive \(N\equiv2\pmod{12}\) and the variety
\[
 V_N:\quad x^2+y^2=E(u,v)=a,\qquad
 z^2+w^2=E(r,s)=b,\qquad a+b=N,
 \tag{15}
\]
where \(a,b\) denote the displayed polynomial values in eight coordinates.
Positive definiteness gives \(0\le a\le N\) on real points. In fact
\(V_N(\mathbb R)\) is compact, since
\(E(u,v)\ge(u^2+v^2)/2\).

Here are integral local points at every prime. For \(p\ge5\), both binary
norms represent every \(p\)-adic unit: modulo \(p\), a split norm is a
product and an inert norm is the surjective finite-field norm; a nonzero
value has a nonzero gradient and lifts by Hensel. Choose \(a\) modulo \(p\)
outside \(\{0,N\}\), and represent the units \(a,N-a\) independently by
both norms.
At \(2\), choose \(a=1,b=N-1\), both \(1\pmod4\). Every unit
\(t\equiv1\pmod4\) is a Gaussian norm: if \(t\equiv1\pmod8\) it is a
square; otherwise \(t-4\) is a square and \(t=(\sqrt{t-4})^2+2^2\).
Every odd \(t\) is an Eisenstein norm, because
\(y^2-y+1-t=0\) has a root modulo \(2\) with odd derivative.
At \(3\), the same \(a,b\) are \(1\pmod3\), hence squares and Eisenstein
norms. Every \(3\)-adic unit is a Gaussian norm: the residue solutions
\((1,0)\) and \((1,1)\) for values \(1,2\) have a unit partial derivative.

The required open image is explicit at \(13\). Lift
\(i=5\pmod{13}\), \(\omega=3\pmod{13}\) to
\(i^2=-1\), \(\omega^2+\omega+1=0\) in \(\mathbb Z_{13}\).
For every \(t\in\mathbb Z_{13}\), put
\[
 x=\frac{1+t}{2},\quad y=\frac{1-t}{2i},\qquad
 v=\frac{t-1}{\omega^2-\omega},\quad u=1-\omega v.
 \tag{16}
\]
All denominators are units. The Gaussian norm is \(t\), while
\[
 E(u,v)=(u+\omega v)(u+\omega^2v)=1\cdot t.
\]
Apply (16) to \(t\) and \(N-t\). This gives a section
\(\mathbb Z_{13}\to V_N(\mathbb Z_{13})\) whose \(a\)-value is \(t\).
Thus the proposition applies to every admissible \(N\).

Allowing units in the quartic-element model does not remove this particular
obstruction. For \(K=\mathbb Q(\zeta_{12})\), let
\[
 W_N:\quad N_{K/\mathbb Q}(\alpha)+N_{K/\mathbb Q}(\beta)=N
 \tag{17}
\]
in integral-basis coordinates. At the real place the norm is a product of
two squared complex absolute values, so its \(\alpha\)-value still lies in
\([0,N]\). At \(13\), the integral algebra is \(\mathbb Z_{13}^4\), and
\((1,1,1,t)\), \((1,1,1,N-t)\) give a section of (17). For
\(N=2\cdot13^m\), \(m\ge0\), the global integral points
\(\alpha=\beta=(2-\zeta_{12})^m\) provide local points at all other primes:
\(N_{K/\mathbb Q}(2-\zeta_{12})=\Phi_{12}(2)=13\).
Hence strong approximation off infinity also fails for these arbitrarily
large admissible quartic targets. No all-target integral assertion follows.

## What general norm-form theorems leave unresolved

[Browning--Matthiesen, Theorem 1.3, printed p. 1385](https://numdam.org/item/10.24033/asens.2348.pdf)
proves rational Hasse principle and weak approximation for norm equations
whose homogeneous linear forms are pairwise nonproportional.
Their Theorem 5.2 and Remark 5.3, printed p. 1401, give a full-dimensional
\(s\)-variable count with error \(o(T^s)\), also allowing affine shifts
\(O(T)\) while retaining the condition on linear parts.
The forms \(t,N-t\) have proportional linear parts. Homogenization to
\(u,Nv-u\) restores the hypothesis but leaves \(v\) free; an
\(o(T^2)\) error does not isolate the fixed slice \(v=1\). For fixed \(N\),
positivity bounds \(u\) on that slice; varying \(N\) also varies a
coefficient of the homogeneous system. A fixed-target slice estimate, or
suitable uniformity in the target and a shrinking body or growing modulus,
is an additional obligation.

Rational normalization remains possible. In the quartic model, add
\(v=N_{K/\mathbb Q}(z)\ne0\). The three forms \(u,Nv-u,v\) are pairwise
nonproportional, so the rational theorem applies when its local hypotheses
hold. Division by \(z\) changes the norm sum \(Nv\) to \(N\).
This does not ensure that the resulting coordinates are integral.
The strong-approximation obstruction above does not address this rational
Hasse-principle or weak-approximation route, nor approximation away from a
finite prime where denominators are permitted.

[Shute, *Polynomials represented by norm forms via the beta sieve*,
Theorems 1.1 and 1.3 and Assumption 2.2](https://www.cambridge.org/core/journals/journal-of-the-institute-of-mathematics-of-jussieu/article/polynomials-represented-by-norm-forms-via-the-beta-sieve/967CC74D681A1832C0AA219B4DF1A6CB)
has hypotheses that fail for this field. In the regular \(V_4\)-action of
\(\operatorname{Gal}(K/\mathbb Q)\), the identity has cycle lengths
\(1,1,1,1\), and each other element has lengths \(2,2\).
The proportion \(T(G)\) with noncoprime cycle lengths is therefore \(3/4\),
whereas Theorem 1.1 requires \(T(G)\le0.39/(k+j+1)\le0.39\), with
\(k,j\) nonnegative factor counts; the stated improvements remain below
\(1/2\). For Theorem 1.3, the bad unramified primes are precisely
\(p\equiv5,7,11\pmod{12}\), of density \(3/4\) by the fixed-progression
prime number theorem. Assumption 2.2 covers them, outside a finite
set, by progressions of total density at most
\(\sum_j1/(q_j-1)\), and requires
\(\deg(F)\sum_j1/(q_j-1)\le0.32380\). This is impossible for nonconstant
\(F\). These are precise applicability failures of those statements,
not an obstruction to other norm-form methods or rational existence.

The imported mean and prime-counting statements above were checked in the
actual cited lemmas, with polynomial prime-power growth and divisor bounds
verifying their hypotheses. The norm-form applicability comparisons use the
published statements at the exact locators above; the bounded-value
proposition and its local applications have direct proofs here.
These internal checks establish only the stated analytic, algebraic and
finite-adelic obstructions. No pointwise binary lower bound, fixed-prime
marginal asymptotic, novelty determination, external mathematical review or
formal verification is claimed. No finite coverage search is used.
