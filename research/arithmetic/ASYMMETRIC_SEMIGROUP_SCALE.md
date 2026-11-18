# Every-target counts for the asymmetric prime-factor restrictions

The asymmetric weighted convolution can be turned into an unweighted count
at the expected logarithmic scale. The useful step is to restore small
allowed prime factors after a lower sieve has made the remaining weights
bounded. Common factors of the target compensate for the local choices
that would otherwise be lost there.

The inputs are [AsymPointwise](ASYMMETRIC_GAUSSIAN_EISENSTEIN_SUMS.md) and
its [uniform small-coefficient extension](ASYMMETRIC_UNIFORM_COEFFICIENTS.md).
The argument below gives comparability, rather than an asymptotic with an
identified leading constant. It has no balance restriction or numerical
cutoff. Publication priority and external specialist review remain open.

## Theorem AsymSemigroupScale

Include 1 in the sets
\[
 C=\{n\ge1:p\mid n\Rightarrow p\equiv1\pmod{12}\},\qquad
 A=\{n\ge1:p\mid n\Rightarrow p\equiv1\pmod3\}.
\]
Define the ordered unweighted count and its target factor by
\[
 r_{C,A}(N)=\#\{(a,b)\in C\times A:a+b=N\},\qquad
 J_E(N)=\prod_{\substack{p\mid N,\ p\ge5\\p\equiv5,11\ (12)}}
                 \frac{p-1}{p-2}.
\]
There are constants \(c,C_1>0\) and an integer \(N_0\) such that every
\(N\ge N_0\), \(N\equiv2\pmod6\), satisfies
\[
 \boxed{\displaystyle
 c\frac{N J_E(N)}{(\log N)^{5/4}}
 \le r_{C,A}(N)\le
 C_1\frac{N J_E(N)}{(\log N)^{5/4}}.}
 \tag{1}
\]
Both inequalities are uniform in all prime valuations of the target.
The first summand satisfies both splitting restrictions; the second
satisfies the Eisenstein restriction. No \(C+C\) conclusion is included.

## 1. The uniform analytic input and finite predicates

Put \(U=\tau^2\mathbf1_C\), \(V=\tau\mathbf1_A\). Theorem AsymUniform
provides fixed \(\vartheta,\varepsilon_0>0\) for the normalized families
\[
 F_B(x)=\frac{U(Bx)}{U(B)},\qquad
 G_D(y)=\frac{V(Dy)}{V(D)}.
\]
Uniformly for \(B,u\in C\), \(D,v\in A\), each at most
\(M^\vartheta\), its estimate for
\(\sum_{ux+vy=M}F_B(x)G_D(y)\) has error
\(O(M^{1-\varepsilon_0})\). Its main is the real length \(M/(uv)\)
times the explicit normalized local product. Coefficient gcds, target
valuations and zero local factors are all retained.

Choose small fixed \(\kappa,\xi>0\), after these analytic constants, and
write \(B_0=N^\kappa\), \(L=N^\xi\). Consider squarefree pairwise
coprime triples \((d,u,v)\) with
\[
 d\mid N,\quad d,u\in C,\quad v\in A,\quad (uv,N)=1,
 \qquad duv\le B_0.
 \tag{2}
\]
At \(M=N/d\), count \(x+y=M\) with weight
\(U(x)V(y)/(4^{\omega(u)}2^{\omega(v)})\), imposing
\(v_p(x)=1\) for \(p\mid u\), \(v_p(y)=1\) for \(p\mid v\),
and \(\ell\mid xy\), where \(\ell\le L\) is squarefree and
\((\ell,6uv)=1\).

These predicates are covered by AsymUniform. Exact valuations expand
through divisibility by \(u d_u\), \(v d_v\), with \(d_u\mid u\),
\(d_v\mid v\). At each \(p\mid\ell\), use
\[
 \mathbf1_{p\mid xy}=\mathbf1_{p\mid x}+
    \mathbf1_{p\mid y}-\mathbf1_{p\mid x,\ p\mid y}.
\]
Every resulting coefficient and shift is individually at most \(B_0^2L\).
A forbidden shift makes its term identically zero; otherwise the term
belongs to the analytic family. Thus it suffices that
\[
 2\kappa+\xi<(1-\kappa)\vartheta.
 \tag{3}
\]
The expansion has only \(N^{o(1)}\) terms and divisor-bounded scalar
prefactors. Its total error for each predicate is
\(N^{1-\varepsilon_0+o(1)}\). Changes of variables in the local measures
give exactly the same imposed-predicate main terms. In particular the
both-divisibility term uses the noncoprime coefficient part of AsymUniform.

## 2. Exact rough local masses

For \(p\ge5\), put \(t=p^{-1}\), \(q_p=1-t\),
\[
 A_i=U(p^i),\quad B_i=V(p^i),\quad
 M_{A,p}=q_p\sum_{i\ge0}A_it^i,\quad
 M_{B,p}=q_p\sum_{i\ge0}B_it^i.
\]
Let \(I_{p,e}\) be the local integral in AsymPointwise. The mass where
both components are units on \(x+y=M\) is
\[
 c_e(p)=\begin{cases}1-2/p,&e=0,\\1-1/p,&e\ge1,\end{cases}
 \qquad e=v_p(M).
\]
Consequently the probability of a bad roughness event is
\(g_e(p)=1-c_e(p)/I_{p,e}\). It is zero in classes 5 and 11, and also
in class 7 when \(e\ge1\). Directly from the local geometric series,
\[
 \begin{array}{c|cc}
  &e=0&e\ge1\\\hline
 p\equiv1\ (12)&6/p+O(p^{-2})&8/p+O(p^{-2})\\
 p\equiv7\ (12)&2/p+O(p^{-2})&0
 \end{array}
 \tag{4}
\]
The bounds are uniform in \(e\). Thus \(0\le g_e(p)<1\) and
\(g_e(p)\le K_0/p\) for one fixed \(K_0\). For every finite set of
small primes, \(1-g_e(p)\) is bounded below uniformly in \(e\): the
integrals have uniformly bounded polynomial-geometric tails and
\(c_e(p)\ge3/5\).

Define
\[
 \mathcal V(z)=\prod_{5\le p\le z}
       \frac{1-2/p}{M_{A,p}M_{B,p}},\qquad
 E_p=\frac{p-1}{p-2},\quad
 E_{p,e}=\begin{cases}1,&e=0,\\E_p,&e\ge1,\end{cases}
 \qquad a_p=\frac{p-1}{p(p-2)}.
\]
Fixed-progression prime Mertens gives
\[
 \mathcal V(z)\asymp(\log z)^{-2}.
 \tag{5}
\]
Indeed its factors are \(1-6/p+O(p^{-2})\) in class 1,
\(1-2/p\) in class 7, and \(1+O(p^{-2})\) in classes 5 and 11.
Each class has prime density \(1/4\).

Take \(z=N^\gamma\), \(w=N^\lambda<z\), and restrict all cofactor
primes to \(p\le w\). With the positive base residues \(\rho_U,\rho_V\),
put
\[
 C_N(z)=6\rho_U\rho_V
   \prod_{p>z}\frac{I_{p,v_p(N)}}{M_{A,p}M_{B,p}}.
\]
This is bounded above and below by fixed positive constants uniformly in
\(N\). Away from the target the tail is \(1+O(1/z)\). At most
\(\log N/\log z\) target primes exceed \(z\), and each changes its
factor by \(1+O(1/z)\), uniformly in its valuation. Dividing by a permitted
\(d\) leaves this tail unchanged.

At a unit target prime, a first valuation-one directive has mass
\(q_pU(p)/p\), and a second directive has mass \(q_pV(p)/p\).
Dividing by 4 or 2 leaves \(q_p/p\); relative to \(1-2/p\), this is
exactly \(a_p\). Therefore the exact main after roughness up to \(z\)
for a triple (2) is
\[
 \mathcal M_{d,u,v}=C_N(z)\frac Nd\,\mathcal V(z)
      \prod_{\substack{p\mid N/d\\p\le z}}E_p
      \prod_{p\mid uv}a_p.
 \tag{6}
\]
Directive primes replace the unit-unit condition by the stated exact
valuation. Every other prime up to \(z\), including primes of \(d\),
is required to divide neither residual component.

## 3. An elementary lower sieve at a fixed power level

Choose a fixed \(p_0\) so that \(g_e(p)\le1/2\) above it; treat
\(p\le p_0\) exactly at a fixed divisor cost. Split the other primes
into blocks
\((z^{2^{-j-1}},z^{2^{-j}}]\). In each block, uniformly in targets
and directives, the sum of its local bad densities is at most a fixed
\(M_0\). The true independent-model unit probability is at least
\(e^{-2M_0}\).

For block \(j\), let \(U_j\) be even Bonferroni through degree
\(2(s+j)\), and \(L_j\) odd Bonferroni through degree \(2(s+j)+1\).
Then \(U_j\ge\mathbf1_{\rm no\ bad}\ge L_j\) and \(U_j\ge0\).
A valid lower polynomial for all blocks is
\[
 \Lambda^-=\prod_jU_j-
     \sum_j(U_j-L_j)\prod_{i\ne j}U_i
       \le\mathbf1_{\rm no\ bad\ in\ any\ block}.
 \tag{7}
\]
To check the inequality, telescope the excess of the product of the
upper polynomials over the product of the exact indicators. Each excess
is bounded by the corresponding \(U_j-L_j\).

In the independent local model,
\[
 \mathbb E(U_j-L_j)\le
   \frac{M_0^{2(s+j)+1}}{(2(s+j)+1)!}.
\]
Divide by \(\mathbb E U_j\ge e^{-2M_0}\) and choose one fixed \(s\)
so their sum is at most \(1/4\). Factoring \(\prod_j\mathbb E U_j\)
shows that (7) has mean at least half the full unit probability. This
argument avoids multiplying possibly negative lower polynomials.

Every divisor in (7) is squarefree, with absolute coefficient at most
\(J+1\), where \(J=O(\log\log N)\) is the block count, and size at most
\[
 C_0z^S,\qquad S=\sum_{j\ge0}(2(s+j)+1)2^{-j}=4s+6.
\]
Choose \(\gamma=\xi/(S+1)\); then the support is at most \(L=N^\xi\)
for all large \(N\). Applying the predicate estimates to (7) gives at
least \(\mathcal M_{d,u,v}/2\), apart from their summed errors.
There are at most \(B_0N^{o(1)}\) triples (2), so the total error is
\[
 O\bigl(N^{1-\varepsilon_0+\kappa+\xi+o(1)}\bigr).
 \tag{8}
\]
Choose \(\kappa+\xi<\varepsilon_0/4\), as well as (3). This is
negligible compared with \(N/(\log N)^A\) for every fixed \(A\).

## 4. Restoring the logarithmic mass

Temporarily sum (6) over all valid squarefree cofactors supported on
\(p\le w\), without the size budget. The model becomes a product.
At a class-1 prime outside the target its options are \(1,a_p,a_p\).
At a class-1 prime with \(e=v_p(N)\ge1\), the options from omitting
or including it in \(d\) sum to
\[
 E_p+p^{-1}E_{p,e-1}\ge E_p+1/p=1+2a_p.
 \tag{9}
\]
For \(e=1\) equality holds; for \(e\ge2\) the excess is
\(1/[p(p-2)]\). At a class-7 prime outside the target the options are
\(1,a_p\). If it divides the target there is no cofactor option, but
\[
 E_p-(1+a_p)=\frac1{p(p-2)}>0.
 \tag{10}
\]
At classes 5 and 11 only the target factor remains. Prime Mertens now
gives the lower model mass
\[
 \begin{split}
 &\gg N J_E(N)\mathcal V(z)
  \prod_{\substack{p\le w\\p\equiv1\ (12)}}(1+2a_p)
  \prod_{\substack{p\le w\\p\equiv7\ (12)}}(1+a_p)\\
 &\gg N J_E(N)(\log N)^{-2}(\log w)^{3/4}
  \gg\frac{N J_E(N)}{(\log N)^{5/4}}.
 \end{split}\tag{11}
\]
Target factors above \(z\) have ratio \(1+o(1)\) uniformly. The two
restoration exponents are \(2/4\) and \(1/4\).

Under the normalized actual positive model, a prime is included in
\(d,u\), or \(v\) with probability \(O(1/p)\), uniformly in target
valuations. Hence \(\mathbb E\log(duv)\le C_2\log w\). Choose
\(\lambda<\min(\gamma,\kappa/(4C_2))\). Markov then retains at least
three quarters of the actual total model on \(duv\le B_0\). Applying
Markov only to the comparison expression in (11) would not suffice.
The lower sieve and (8) therefore retain the order in (11).

## 5. Injective recovery and bounded weights

Send a residual pair to the original values \((dx,dy)\). Since
\(d\in C\), these lie in \(C\times A\) and sum to \(N\). At primes
\(p\le w\), the recoverable nonzero valuation patterns are
\((1,1)\) for \(d\), \((1,0)\) for \(u\), and \((0,1)\) for \(v\).
All other valuations there are absent by residual roughness. Thus an
original ordered pair determines its tuple and residual pair uniquely,
including when the target has higher prime valuations.

Write \(x=ux'\), \(y=vy'\). Exact valuation one gives
\[
 \frac{U(x)V(y)}{4^{\omega(u)}2^{\omega(v)}}=U(x')V(y').
\]
All prime factors of \(x',y'\) exceed \(N^\gamma\), so each has at
most \(1/\gamma\) factors with multiplicity. Since
\((e+1)^2\le4^e\) and \(e+1\le2^e\),
\[
 U(x')V(y')\le8^{1/\gamma}.
\]
Injectivity turns (11) into the lower unweighted bound (1). No weight
\(U(d)V(d)\) is inserted: the counted weights are on the residual
values, while the lifting step concerns values only. There is no remaining
large bad-prime tail; the indicators in \(U,V\) enforce all restrictions.

## 6. Matching upper bound

Sieve the integer \(a\in(0,N)\). At class-5 or class-11 primes exclude
\(a\equiv0,N\pmod p\); these coincide exactly when \(p\mid N\).
At class-7 primes exclude only \(a\equiv0\pmod p\). The unit product
has dimension
\(2(2/4)+1/4=5/4\), with target coincidences multiplying it by \(E_p\)
only in classes 5 and 11. Ignoring the fixed conditions at 2 and 3 changes
only an upper constant.

Use the product of the nonnegative even block polynomials from Section 3.
It majorizes the sifted indicator, and its model is at most a fixed
multiple of the true unit product: each block excess is bounded by its
next odd elementary sum, and the relative excesses are summable. Choose
a fixed-power sieve cutoff with divisor level at most \(N^{1/2}\).
CRT counts have the form
\[
 \frac{N\nu(\ell)}\ell+O(\nu(\ell)),\qquad\nu(p)\le2.
\]
The coefficients are bounded, so their total endpoint error is at most
\(O(N^{1/2}\log^C N)\), using
\(\sum_{\ell\le D}2^{\omega(\ell)}\ll D\log(2D)\).
The model product is \(\ll J_E(N)/(\log N)^{5/4}\); large target primes
only increase the full \(J_E(N)\). This proves the upper bound in (1).

## Dependencies and evidence

The only analytic progression input beyond AsymPointwise is AsymUniform,
with its explicitly chosen small family power. Prime Mertens for the
fixed progressions is the standard input recorded in
[Loughran--Matthiesen, Lemma 2.4](https://arxiv.org/abs/1904.12845).
The lower sieve, local compensation, truncation and removal of weights
are proved above. No almost-every-target theorem is used to infer a
pointwise conclusion.

Separate internal mathematical reviews reconstructed the uniform input
and the cofactor argument. The retained finite checker tests selected
exact local compensations, lower-polynomial signs, and injective fibres,
including common factors at targets divisible by \(13^2\). Those finite
controls do not verify the global analytic asymptotic. External review,
formal verification, an identified leading constant, effectivity and a
historical priority claim are outside the result.
